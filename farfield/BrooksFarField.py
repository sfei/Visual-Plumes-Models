import math
import numpy as np
from .. import units
from ..Element import Element
from ..Output import OutputFarField
from ..ambient.Ambient import Ambient
from ..ambient.AmbientHandler import AmbientHandler
from ..ambient.AmbientStore import AmbientStore
from ..helpers import num_format
from ..params.DiffuserStore import DiffuserStore
from ..params.ModelParameters import ModelParameters, FarfieldDiffusivity
from ..timeseries.TimeseriesHandler import TimeseriesHandler
from ..vectors import DEGREES_TO_RADIAN, ZERO_VECTOR, K_UNIT_VECTOR, GRAVITY_VECTOR, make_vector, rescale_vector, \
    magnitude, project_vector, angle, unit_vector, is_zero_vector


def run_far_field(umunit, brooks_ff_model, timeseries):
    if not umunit.model_params.brooks_far_field:
        return
    if umunit.status['atvmxz']:
        outputit = OutputFarField()
        outputit.memo('Mixing Zone reached in near-field, no far-field calculation attempted')
        return outputit
    elif umunit.ambient.ff_velocity == 0:
        outputit = OutputFarField()
        outputit.memo('No farfield prediction when far vel = 0.')
        return outputit
    elif umunit.um3isoplet <= 0:
        outputit = OutputFarField()
        outputit.memo('Isopleth closed in near-field, no far-field prediction necessary.')
        return outputit
    else:
        # TODO: assumed surfaced, move conditions to, if not
        umunit.element.depth             = 0.0
        umunit.element.v_velocity[2]     = 0.0
        umunit.element.v_surface_tdsp[2] = 0.0
        umunit.element.total_surf_dsp    = magnitude(umunit.element.v_surface_tdsp)

        angel = angle(
            project_vector(umunit.orig_element.v_velocity, GRAVITY_VECTOR),
            project_vector(umunit.v_ambient, GRAVITY_VECTOR)
        )
        cos_angel = abs(math.cos(angel))
        min_cos = math.cos(70*DEGREES_TO_RADIAN)
        if cos_angel < min_cos:
            cos_angel = min_cos
        width = umunit.element.diameter
        width += round(umunit.diff_params.num_ports - 1)*umunit.diff_params.port_spacing*cos_angel

        messages = []
        if umunit.element.diameter < umunit.diff_params.port_spacing and umunit.diff_params.num_ports > 1:
            messages.append('Plumes not merged, Brooks method may be overly conservative.')

        # TODO: what is the point of these?
        # faro = umunit.element.dilution
        # xold = umunit.element.v_surface_tdsp[0]
        # plotting code skipped (523-540)0
        # (oc) in farfield no centerline plot, so plot must go to outline

        # o.c. code just assumed isopleth was in concentration
        match umunit.diffuser_store.isopleth.units:
            case units.Isopleth.CONCENTRATION:
                isopleth_conc = umunit.diff_params.isopleth
            case _:
                isopleth_conc = -1

        return brooks_ff_model.run(
            messages,
            umunit.diffuser_store,
            umunit.v_source,
            umunit.diff_params.acute_mixing_zone,
            umunit.ambient.kt,
            isopleth_conc,
            umunit.element,
            width,
            umunit.mass_pollutant,
            umunit.orig_element.mass,
            umunit.ambient_stack,
            umunit.ambient_store,
            timeseries,
            umunit.graphit
        )


class BrooksFarFieldRun:

    def __init__(self, model_params, element, v_shore, v_source, acute_mixing_zone, mass_pollutant, width, diff_store,
                 ambient_store):
        assert isinstance(model_params, ModelParameters)
        assert isinstance(element, Element)
        assert isinstance(diff_store, DiffuserStore)
        assert isinstance(ambient_store, AmbientStore)

        self.model_params       = model_params
        self.v_shore            = v_shore
        self.v_perp_shore       = rescale_vector(1, K_UNIT_VECTOR*v_source)  # TODO: used anywhere?
        self.ambient            = Ambient()
        self.v_ambient          = None
        self.element            = element
        # formerly `So`  (just use element.dilution)
        # formerly `pm4` (just use element.mass)
        # formerly `dm`  (just use element.d_mass)
        # formerly `Vf`  (just use element.v_velocity)
        # formerly `Rf`  (just use element.v_surface_tdsp)
        self.mass_pollutant     = mass_pollutant
        self.width              = width
        self.adj_width          = self.width
        self.diffusivity        = 0  # formerly `arg` or `v4o3`
        self.farw               = 0
        self.step               = 0
        self.run_time           = 1  # far-field model run time -- (oc) erf time in seconds -- formerly `tim`
        self.total_time         = element.total_time  # total run time (based on case and UM3 simulated time)
        self.increment          = self.model_params.ff_increment  # formerly `ink`
        self.outputit           = OutputFarField()

        # set at stop depth last given in plume model run
        self.ambient.depth = self.element.depth

        while self.element.total_surf_dsp > self.increment and self.model_params.ff_increment > 0:
            self.increment += self.model_params.ff_increment
        if self.increment > acute_mixing_zone:
            self.increment = acute_mixing_zone

        match model_params.farfield_diffusivity:
            case FarfieldDiffusivity.POWER_4_3:
                diff_col_name = "4/3 eddy diffusivity"
            case FarfieldDiffusivity.CONSTANT:
                diff_col_name = "Eddy diffusivity"
            case _:
                diff_col_name = "Eddy diffusivity"
        self.outputit.add_parameter('model',    'run_time',         'Time',                     units.Time,             units.Time.HOURS)
        self.outputit.add_parameter('ambient',  'ff_velocity',      'Current speed',            units.Speed,            ambient_store.ff_velocity.units)
        self.outputit.add_parameter('ambient',  'ff_dir',           'Current direction',        units.Angle,            ambient_store.ff_dir.units)
        self.outputit.add_parameter('model',    'adj_width',        'Width',                    units.Length,           units.Length.METERS)
        self.outputit.add_parameter('element',  'total_surf_dsp',   'Distance',                 units.Length,           units.Length.METERS)
        self.outputit.add_parameter('ambient',  'bg_conc',          'Background concentration', units.Concentration,    ambient_store.bg_conc.units)
        self.outputit.add_parameter('element',  'concentration',    'Concentration',            units.Concentration,    diff_store.concentration.units)
        self.outputit.add_parameter('element',  'dilution',         'Dilution',                 units.Unitless,         units.Unitless.UNITLESS)
        self.outputit.add_parameter('ambient',  'kt',               'Decay rate',               units.DecayRate,        units.DecayRate.PER_SECOND)
        self.outputit.add_parameter('ambient',  'ff_diff_coeff',    'Dispersion',               units.EddyDiffusivity,  units.EddyDiffusivity.DIFFUSIVITY)
        self.outputit.add_parameter('model',    'diffusivity',      diff_col_name,              units.Unitless,         units.Unitless.UNITLESS)

    def output(self, dum_si=None, dum_conc=None, dum_sec=None, dum_displace=None):
        # certain outputs use dummy values, so save old values
        og_dilution           = self.element.dilution
        og_concentration      = self.element.concentration
        og_run_time           = self.run_time
        og_total_displacement = self.element.total_surf_dsp

        if dum_si is not None:
            self.dilution = dum_si
        if dum_conc is not None:
            self.concentration = dum_conc
        if dum_sec is not None:
            self.run_time = dum_sec
        if dum_displace is not None:
            self.element.total_surf_dsp = dum_displace

        # adjusted width
        self.adj_width = 8.0*self.ambient.ff_diff_coeff*self.run_time*self.ambient.ff_velocity
        self.adj_width = self.adj_width/(self.ambient.ff_velocity*math.pow(self.width, 2.0/3.0))
        self.adj_width = self.width*math.pow(1.0+self.adj_width, 1.5)

        self.outputit.output(self)

        # reset params
        if dum_si is not None:
            self.dilution = og_dilution
        if dum_conc is not None:
            self.concentration = og_concentration
        if dum_sec is not None:
            self.run_time = og_run_time
        if dum_displace is not None:
            self.element.total_surf_dsp = og_total_displacement


class BrooksFarField:

    def __init__(self, model_params):
        assert isinstance(model_params, ModelParameters)
        self.decay_memory = [0]*24  # memory of previous decay rates
        self.model_params = model_params
        self.v_shore = make_vector(
            x=(self.model_params.dist_to_shore * math.cos(DEGREES_TO_RADIAN*self.model_params.dir_to_shore)),
            y=(self.model_params.dist_to_shore * math.sin(DEGREES_TO_RADIAN*self.model_params.dir_to_shore))
        )
        # Only used in ambient fill, as replacement for umunit.ambient.temperature in decay-rate conversions. Probably
        # unnecessary as only ff decay rate used but whatever for now.
        self.ambient = None

    def get_ambient(self, ambient_handler, ambient_stack, ambient_store, timeseries, ff):
        for layer in ambient_stack:
            assert isinstance(layer, Ambient)
        assert isinstance(ambient_store, AmbientStore)
        if timeseries:
            assert isinstance(timeseries, TimeseriesHandler)
            timeseries.set_time(ff.run_time)
            ambient_ts_stacks = timeseries.get_ambient()
        else:
            ambient_ts_stacks = {}

        bottom_depth = self.model_params.bottom_depth
        if ambient_store.z.z_is_depth:
            max_z = max((amb.z for amb in ambient_stack))
            if max_z > bottom_depth:
                bottom_depth = max_z
        for vkey in ambient_store._input_vars_:
            ambient_store_at = ambient_store.get(vkey)
            if not ambient_store_at.from_time_series or not ambient_store_at.z_is_depth:
                continue
            max_z = max((amb.z for amb in ambient_ts_stacks[vkey]))
            if max_z > bottom_depth:
                bottom_depth = max_z
        self.model_params.bottom_depth = bottom_depth

        ambient_handler.fill(
            model_params=self.model_params,
            ambient_stack=ambient_stack,
            ambient_ts_stacks=ambient_ts_stacks
        )
        ambient_handler.ambient_level(
            umunit=None,
            ambient_cond=ff.ambient,
            depth=ff.ambient.depth,
            bottom_depth=bottom_depth
        )

        mag_shore = magnitude(self.v_shore)
        v_unit_shore = rescale_vector(1, self.v_shore)
        ff.v_ambient = make_vector(
            x=ff.ambient.ff_velocity*math.cos(DEGREES_TO_RADIAN*ff.ambient.ff_dir),
            y=ff.ambient.ff_velocity*math.sin(DEGREES_TO_RADIAN*ff.ambient.ff_dir),
            z=0
        )
        if self.model_params.farfield_diffusivity == FarfieldDiffusivity.POWER_4_3:
            dot_disp_shore = np.dot(ff.element.v_surface_tdsp, v_unit_shore)
            if dot_disp_shore >= 0 and dot_disp_shore > mag_shore:
                ff.v_ambient = np.dot(ff.v_ambient, ff.v_perp_shore) * ff.v_perp_shore
        else:
            ff.v_ambient = (
                (mag_shore - np.dot(v_unit_shore, ff.element.v_surface_tdsp))/mag_shore
                * (np.dot(ff.v_ambient, v_unit_shore)*v_unit_shore)
            )
            ff.v_ambient += np.dot(ff.v_ambient, ff.v_perp_shore)*ff.v_perp_shore

    def get_arg(self, ff_diff_coeff, width, time):
        w2 = width**2
        a = 8.0*ff_diff_coeff*math.pow(width, 4.0/3.0)*time
        if self.model_params.farfield_diffusivity == FarfieldDiffusivity.POWER_4_3:
            return 1.5/((1.0 + a/w2)**3 - 1.0)
        else:
            return w2/(2.0*a)

    @staticmethod
    def estbackmult(decay_memory, dt):
        def f(i):
            return math.exp(-decay_memory[i-1]*dt)
        return (
            f(1)*(1+f(2)*(1+f(3)*(1+f(4)*(1+f(5)*(1+f(6)*(1+f(7)*(1+f(8)*(1+f(9)
            *(1+f(10)*(1+f(11)*(1+f(12)*(1+f(13)*(1+f(14)*(1+f(15)*(1+f(16)*(1+f(17)
            *(1+f(18)*(1+f(19)*(1+f(20)*(1+f(21)*(1+f(22)*(1+f(23)*(1+f(24))))))))))))))))))))))))
        )

    def run(self, messages, diff_store, v_source, acute_mixing_zone, kt, isopleth_conc, element, width,
            mass_pollutant, orig_mass, ambient_stack, ambient_store, timeseries, graphit):
        assert isinstance(ambient_store, AmbientStore)
        for layer in ambient_stack:
            assert isinstance(layer, Ambient)
        if timeseries:
            assert isinstance(timeseries, TimeseriesHandler)

        ff = BrooksFarFieldRun(self.model_params, element, self.v_shore, v_source, acute_mixing_zone, mass_pollutant,
                               width, diff_store, ambient_store)

        for m in messages:
            ff.outputit.memo(m)

        # why graph twice? already done in UMUnit
        # if self.model_params.use_shore_vector:
        #     graphit.graph_vector(v_source, self.v_shore)

        # decay rate memory, separate array for per-case-run, and individual memory per run
        # note: using decay rate at end of UM3 model run, ff specific decay rate will get calculated when grab ambient
        self.decay_memory = [kt] + self.decay_memory[1:]
        run_decay_memory = self.decay_memory

        ambient_handler = AmbientHandler(self.model_params, ambient_store)
        self.get_ambient(ambient_handler, ambient_stack, ambient_store, timeseries, ff)

        if self.model_params.farfield_diffusivity == FarfieldDiffusivity.POWER_4_3:
            ff.outputit.memo(f"4/3 Power Law. Farfield dispersion based on wastefield width of {num_format(width)} m")
        else:
            ff.outputit.memo(f"Const Eddy Diffusivity. Farfield dispersion based on wastefield width of {num_format(width)} m")

        # determine timestep (separate from umunit.timestep)
        has_timeseries = (
            ambient_store.bg_conc.from_time_series
            or ambient_store.decay_rate.from_time_series
            or ambient_store.ff_velocity.from_time_series
            or ambient_store.ff_dir.from_time_series
            or ambient_store.ff_diff_coeff.from_time_series
        )
        if has_timeseries:
            # TODO: this code means time-series must be supplied for all or none of the above.. and in same time increments
            dt = ambient_store.ff_velocity.ts_increment*3600
        else:
            dt = ff.increment/ff.ambient.ff_velocity if ff.increment > 0 else 360

        odo = K_UNIT_VECTOR
        start_concentration = ff.element.concentration  # formerly cvpolo
        start_dilution = ff.element.dilution  # formerly Sid
        ff.output()

        for i in range(10000):
            if ff.element.total_surf_dsp >= acute_mixing_zone:
                break
            if isopleth_conc >= 0 and ff.element.concentration < isopleth_conc:
                break

            # move conditions over timeseries
            if has_timeseries and i > 0:
                self.get_ambient(ambient_handler, ambient_stack, ambient_store, timeseries, ff)

            # move timestep
            ff.step += 1
            ff.total_time += dt
            ff.run_time += dt

            # adjust dilution, mass, velocity
            ff.diffusivity = self.get_arg(ff.ambient.ff_diff_coeff, width, ff.run_time)
            new_dilution = start_dilution/math.erf(ff.diffusivity**0.5)  # formerly `newS`
            ff.element.d_mass = (new_dilution - ff.element.dilution)*orig_mass  # (oc) pm0 was wrong per Kenwyn experience 21 July 2001
            ff.element.v_velocity = (
                # (oc) _s( _sx(pm4,Vf), _sx(dm,Uam) )  *  1/(pm4+dm)
                (ff.element.mass*ff.element.v_velocity + ff.element.d_mass*ff.v_ambient)
                / (ff.element.mass + ff.element.d_mass)
                + ZERO_VECTOR  # TODO: isn't vector sum with zero vector pointless?
            )
            # (oc) bug Henry 2017 Salas, simplified next line
            v_prev_displacement        = np.copy(ff.element.v_surface_tdsp)
            prev_displacement          = ff.element.total_surf_dsp
            ff.element.v_surface_tdsp += dt*ff.ambient.ff_velocity*(ff.element.v_surface_tdsp/prev_displacement)
            ff.element.total_surf_dsp  = magnitude(ff.element.v_surface_tdsp)

            # next chain of this thing
            run_decay_memory = [ff.ambient.kt] + run_decay_memory[1:]

            backmult = self.estbackmult(run_decay_memory, dt)
            backvpol = ff.ambient.bg_conc
            if not self.model_params.estimate_ff_background:
                backmult = 0

            ff.mass_pollutant       += ff.element.d_mass*backvpol  # (oc) bug 2013 post-SF still needs work
            ff.element.mass         += ff.element.d_mass
            ff.element.concentration = ff.mass_pollutant/ff.element.mass
            ff.element.dilution      = new_dilution

            while ff.element.total_surf_dsp > ff.increment:
                time2ink = dt*(ff.increment - prev_displacement)/(ff.element.total_surf_dsp - prev_displacement)
                diffusivity = self.get_arg(ff.ambient.ff_diff_coeff, width, ff.run_time-dt+time2ink)
                dumconc = (
                    (mass_pollutant + ff.element.d_mass*backvpol/2.2*backmult)
                    *math.exp(-ff.ambient.kt*time2ink)
                    /(ff.element.mass - ff.element.d_mass)
                )
                if diffusivity >= 0:
                    dumSi = ff.element.dilution / math.erf(diffusivity**0.5)
                else:
                    dumSi = ff.element.dilution / math.erf(-(abs(diffusivity))**0.5)
                ff.output(
                    dum_conc=dumconc,
                    dum_si=dumSi,
                    dum_sec=ff.run_time-dt+time2ink,
                    dum_displace=ff.increment
                )
                ff.increment += self.model_params.ff_increment

            ff.farw = ff.ambient.ff_diff_coeff*ff.run_time*ff.ambient.ff_velocity
            ff.farw /= ff.ambient.ff_velocity*width**(2.0/3.0)
            if self.model_params.farfield_diffusivity == FarfieldDiffusivity.POWER_4_3:
                ff.farw = width*math.pow(1+8*ff.farw, 1.5)
            else:
                ff.farw = width*math.pow(1+24*ff.farw, 0.5)

            # BdyPlumeEndPlot() and plotconc()
            if 0 < isopleth_conc < ff.element.concentration:
                odo = (
                    0.5*ff.farw
                    * ((math.log(ff.element.concentration) - math.log(isopleth_conc))/6.0)**0.5
                    * unit_vector(np.cross(ff.element.v_velocity, GRAVITY_VECTOR))
                )
                graphit.plot_ff_concentration(v_prev_displacement, odo)
            elif is_zero_vector(odo):
                odo = ZERO_VECTOR
                adj_displacement = v_prev_displacement + rescale_vector(
                    (start_concentration - isopleth_conc)/(start_concentration - ff.element.concentration),
                    ff.element.v_surface_tdsp - v_prev_displacement
                )
                graphit.plot_ff_concentration(adj_displacement, odo)

            # (oc) Salas bug 2014 omission of dilution plotting in far-field
            graphit.plot_dilution(
                magnitude(project_vector(
                    v_source - ff.element.v_surface_tdsp,
                    GRAVITY_VECTOR
                )),
                ff.element.dilution
            )

            if self.model_params.output_all_ff_increments:
                ff.output()

        ff.output()
        # ff.outputit.memo(f"count: {count}")
        # TODO: fmark = fmarkorg
        # TODO: farreset() -- maybe not needed
        # TODO: lifepens()

        return ff.outputit
