import math
from . import units
from .Element import Element
from .ambient.calculations import seawater_density
from .globals import maxsteps, UserInputError
from .ambient.Ambient import Ambient
from .helpers import num_format
from .params.DiffuserParameters import DiffuserParameters
from .params.ModelParameters import Model
from .vectors import magnitude, MAGNITUDE_GRAVITY


_TEMPLATE_AMBIENT  = Ambient()
_TEMPLATE_ELEMENT  = Element()
_TEMPLATE_DIFFUSER = DiffuserParameters()
_TEMPLATE_MODEL    = [
    'density_diff',
    'net_dilution',
    'mass_pollutant',
    'um3isoplet',
    'iso_diameter',
    'travel',
    'steps',
    'reverals'
    # DEBUG
    # 'entrainment_limspc',
    # 'dist_separation'
]
_TEMPLATE_MODEL_BROOKS = [
    'diffusivity',
    'width',
    'adj_width',
    'farw',
    'step',
    'run_time',
    'total_time',
    'increment'
]


class Output:

    def __init__(self, copy=None, model_template=None):
        if model_template is None:
            model_template = _TEMPLATE_MODEL
        self._model_template = model_template
        self._parameters     = []
        self._is_running     = False
        self._headers        = []
        self._units          = []
        self._output_units   = []
        self._outputs        = []
        self._memos          = []
        if not copy or not isinstance(copy, Output):
            return
        for i, p in enumerate(copy._parameters):
            self.add_parameter(
                p[0],
                p[1],
                copy._headers[i],
                copy._units[i],
                copy._output_units[i]
            )

    def memo(self, message):
        self._memos.append(message)

    @property
    def memos(self):
        return self._memos

    @property
    def length(self):
        return len(self._parameters)

    def headers(self):
        yield from (
            {
                'regime':      p[0],
                'name':        p[1],
                'label':       self._headers[i],
                'units':       self._units[i],
                'in_units':    self._output_units[i],
                'units_label': (
                    "" if self._units[i] == units.Unitless else self._units[i].label(
                        self._output_units[i] if self._output_units[i] is not None else 1
                    )
                )
            }
            for i, p in enumerate(self._parameters)
        )

    def outputs(self):
        yield from self._outputs

    def _get_regime(self, regime):
        regime = regime.lower().strip()
        match regime:
            case 'ambient':
                template = _TEMPLATE_AMBIENT
            case 'element':
                template = _TEMPLATE_ELEMENT
            case 'diffuser':
                template = _TEMPLATE_DIFFUSER
            case 'model':
                template = self._model_template
            case _:
                raise UserInputError(f"Invalid parameter regime ({regime})")
        return regime, template

    def _validate_parameter(self, regime, parameter_name):
        regime, template = self._get_regime(regime)
        if isinstance(template, (list, tuple)):
            if parameter_name not in template:
                raise AttributeError()
        else:
            getattr(template, parameter_name)
        return regime, template

    def add_parameter(self, regime, parameter_name, header=None, unit_type=units.Unitless, in_units=None):
        if self._is_running:
            raise Exception("Cannot adjust parameters once model is running")
        try:
            regime, template = self._validate_parameter(regime, parameter_name)
        except AttributeError:
            raise UserInputError(f"Parameter ({parameter_name}) not valid for {regime}")
        for p_regime, p_name in self._parameters:
            if p_regime == regime and p_name == parameter_name:
                return
        assert issubclass(unit_type, units.Units)
        self._parameters.append([regime, parameter_name])
        self._headers.append(header if header else parameter_name)
        self._units.append(unit_type)
        self._output_units.append(in_units)

    def has_parameter(self, regime, parameter_name):
        for (i_regime, i_param) in self._parameters:
            if i_regime == regime and i_param == parameter_name:
                return True
        return False

    def remove_parameter(self, regime, parameter_name):
        if self._is_running:
            raise Exception("Cannot adjust parameters once model is running")
        try:
            regime, template = self._validate_parameter(regime, parameter_name)
        except:
            return
        remove_i = -1
        for i, pair in enumerate(self._parameters):
            if pair[0] == regime and pair[1] == parameter_name:
                remove_i = i
                break
        if remove_i >= 0:
            del self._parameters[remove_i]
            del self._headers[remove_i]
            del self._units[remove_i]
            del self._output_units[remove_i]

    def _initialize(self):
        if self._is_running:
            return
        self._is_running = True

    def output(self, model):
        raise Exception("Method output() called from abstract class")

    def _output_values(self, model):
        output_values = []
        for regime, pname in self._parameters:
            match regime:
                case 'ambient':
                    from_object = model.ambient
                case 'element':
                    from_object = model.element
                case 'diffuser':
                    from_object = model.diff_params
                case 'model':
                    from_object = model
                case _:
                    output_values.append(None)
                    continue
            try:
                output_values.append(getattr(from_object, pname))
            except:
                output_values.append(None)

        for i, in_units in enumerate(self._output_units):
            if in_units is None or output_values[i] is None:
                continue
            if math.isnan(output_values[i]):
                output_values[i] = None
                continue
            output_values[i] = units.convert(
                output_values[i],
                self._units[i],
                from_units=1,
                to_units=in_units,
                model_params=model.model_params,
                celsius=model.element.temperature,
                psu=model.element.salinity,
                depth=model.element.depth
            )
            if math.isnan(output_values[i]):
                output_values[i] = None
            elif (
                self._units[i] == units.Angle
                and in_units == units.Angle.DEGREES
                and output_values[i] > 180
                and self._parameters[i][1] in ("vertical_angle",)
            ):
                # special case for vertical angle in negative degrees
                output_values[i] -= 360

        return output_values


class OutputUM3(Output):

    def __init__(self, copy=None):
        super().__init__(copy)

    def froude_no(self, model):
        if model.element.density == 0:
            return 0, "indeterminate, density=0"
        density_diff = model.ambient.density - model.element.density
        if density_diff == 0:
            return 9.999e15, "pure jet"
        n_density_diff = density_diff / model.element.density
        sign = 1
        if n_density_diff < 0:
            n_density_diff = abs(n_density_diff)
            sign = -1
        froude_no = (
            sign*model.diff_params.effluent_flow
            / (model.diff_params.num_ports*math.pi*(0.5*model.element.diameter)**2)
            / math.sqrt(MAGNITUDE_GRAVITY*model.element.diameter*n_density_diff)
        )
        return froude_no, ""

    def strat_no(self, model):
        arg = seawater_density(
            ambient_cond=model.ambient_stack[0],
            at_equilibrium=model.model_params.at_equilibrium
        )
        return (
            (model.ambient.density - arg)
            *model.element.diameter
            /model.diff_params.depth
            /(model.ambient.density - model.element.density)
        )

    def spcg_no(self, model):
        if model.diff_params.num_ports <= 1:
            return 0
        return model.diff_params.port_spacing/model.element.diameter

    def k_no(self, model):
        return model.element.speed/model.ambient.current_speed

    def output(self, model, denomproduct=0):
        if not self._is_running:
            self._initialize()

        if model.step == 0:
            # self.memo("Simulation:")
            if model.model_params != Model.PDS:
                froude_no, froude_st = self.froude_no(model)
                self.memo(f"Froude No:    {num_format(froude_no)} {froude_st}")
                self.memo(f"Strat No:     {self.strat_no(model):.2e}")
                if model.diff_params.num_ports > 1:
                    self.memo(f"Spcg No:      {num_format(self.spcg_no(model))}")
                self.memo(f"k:            {num_format(self.k_no(model))}")
                self.memo(f"effluent den: {num_format(model.element.density-1000)} sigmaT")
                self.memo(f"effluent vel: {num_format(model.element.speed)} m/s")
                self.memo(f"bottom depth: {num_format(model.model_params.bottom_depth)} m")
                # TODO: this is usually done under condition of model != justlook and model != pds
                if magnitude(model.v_ambient) < 1e-5:
                    self.memo("Current is very small, flow regime may be transient.")
                if abs(froude_no) < 1:
                    self.memo("Absolute value Froude No. < 1, possible intrusion and/or plume diameter reduction")

        output_values = self._output_values(model)

        statuses = []
        if model.step > 0 and model.model_params != Model.PDS and model.density_diff*model.last_density_diff < 0:
            statuses.append("trap level")
        if model.check_status_changed('has_surfaced'):
            statuses.append("surface")
        if model.check_status_changed('atvmxz'):
            statuses.append("MZ dis")
        if model.check_status_changed('bottomed'):
            statuses.append("bottom hit")
        if model.check_status_changed('reached_cig'):
            statuses.append('ceiling')
        if model.status['reversal']:
            statuses.append('local maximum rise or fall')
        if model.check_status_changed('close_isoplet'):
            statuses.append('isopleth closed')
        if model.check_status_changed('stream_limit'):
            statuses.append('stream limit reached')
        if model.check_status_changed('has_merged'):
            statuses.append('merging')
        if model.status['is_lap'] and not model.status['is_lap_old']:
            statuses.append('begin overlap')
        if model.step > maxsteps:
            statuses.append('max iterations')
        if not model.status['is_lap'] and model.status['is_lap_old']:
            statuses.append('end overlap')
        if model.element.dilution > model.model_params.max_dilution:
            statuses.append('stop dilution reached')
        if denomproduct < 0 < model.step and model.element.diameter > model.ambient.depth:
            val = 0.25*model.element.diameter*magnitude(model.element.v_velocity)/model.ambient.depth
            statuses.append(f"matched energy radial vel = {num_format(val)} m/s")

        self._outputs.append({
            'step':   model.step,
            'values': output_values,
            'status': ", ".join(statuses)
        })


class OutputFarField(Output):

    def __init__(self, copy=None):
        super().__init__(copy, _TEMPLATE_MODEL_BROOKS)

    def output(self, model):
        if not self._is_running:
            self._initialize()

        output_values = self._output_values(model)

        self._outputs.append({
            'step':   model.step,
            'values': output_values
        })
