from enum import Enum
from .. import units
from ..helpers import convert_float, is_truthy
from ..globals import UserInputError


class Interpolation(Enum):
    CONSTANT       = 0
    LINEAR_TO_ZERO = 1
    EXTRAPOLATED   = 2


class AmbientStoreSubset:
    """ Store class for desc/metadata about a specific ambient parameter """

    def __init__(self, var_key, default_unit=units.Unitless.UNITLESS):
        # (oc) ambient file order
        self.var_key           = var_key
        self.from_time_series  = False                   # (formerly ffin) time-series file
        self.ts_increment      = 0                       # (formerly fdel) time-series time increment (assumed in hours)
        # self.ts_period         = 0                       # (formerly fper) time-series cycling period (assumed in hours)
        # self.ts_file_units     = default_unit            # (formerly fufi) time-series measurement unit (?) TODO: just use `units`?
        self.ts_depth_units    = units.Length.METERS     # units for depth values in time-series
        self.z_is_depth        = True                    # (formerly fdoh) depth or height
        self.extrapolation_sfc = Interpolation.CONSTANT  # (formerly fsex) surface extrapolation (sfc)
        self.extrapolation_btm = Interpolation.CONSTANT  # (formerly fbex) bottom extrapolation (btm)
        self.units             = default_unit            # (formerly futa) measurement unit
        self._vars_ = list(vars(self).keys())
        self._vars_.remove('var_key')

    @staticmethod
    def from_dict(var_key, input_dict):
        c = AmbientStoreSubset(var_key)
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = AmbientStoreSubset("")
        ret_dict = {}
        for key in template._vars_:
            ret_dict[key] = getattr(template, key)
        return ret_dict

    def to_dict(self):
        ret_dict = {}
        for key in self._vars_:
            ret_dict[key] = getattr(self, key)
        return ret_dict

    def copy(self):
        c = AmbientStoreSubset(self.var_key)
        for key in c._vars_:
            setattr(c, key, getattr(self, key))
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)


# replaces `ambstore` (main.formatted.pas:710) which was a multidimensional array
# could be replaced with pandas dataframe..
# seems to be a dataset of all ambient conditions? Or at least by depth, but the var names are completely impenetrable
class AmbientStore:
    """ Store class for desc/metadata about ambient parameters """
    
    def __init__(self):
        # (oc) ambient file order
        # (formerly `valen`) seems to be z-coordinate (depth)
        self.z              = AmbientStoreSubset('z', units.Length.METERS)
        # (formerly `vaspd`) current speed (flat component only -- no vertical component assumed)
        self.current_speed  = AmbientStoreSubset('current_speed', units.Speed.METERS_PER_SECOND)
        # (formerly `vadir`) current direction in degrees (azimuthal only -- no vertical component assumed)
        self.current_dir    = AmbientStoreSubset('current_dir', units.Angle.DEGREES)
        # (formerly `vasal`) ambient salinity
        self.salinity       = AmbientStoreSubset('salinity', units.Salinity.PRACTICAL_SALINITY_UNITS)
        # (formerly `vatem`) ambient temperature
        self.temperature    = AmbientStoreSubset('temperature', units.Temperature.CELSIUS)
        # (formerly `vapol`) background concentration
        self.bg_conc        = AmbientStoreSubset('bg_conc', units.Concentration.KILOGRAM_PER_KILOGRAM)
        # (formerly `varat`) pollutant decay rate
        self.decay_rate     = AmbientStoreSubset('decay_rate', units.DecayRate.PER_SECOND)
        # (formerly `vafar`) far-field velocity
        self.ff_velocity    = AmbientStoreSubset('ff_velocity', units.Speed.METERS_PER_SECOND)
        # (formerly `vafad`) far-field current direction
        self.ff_dir         = AmbientStoreSubset('ff_dir', units.Angle.DEGREES)
        # (formerly `vadis`) far-field diffusion coefficient
        self.ff_diff_coeff  = AmbientStoreSubset('ff_diff_coeff', units.Unitless.UNITLESS)
        self._vars_ = list(vars(self).keys())
        self._input_vars_ = self._vars_[:]
        self._input_vars_.remove('z')

    @staticmethod
    def from_dict(input_dict):
        c = AmbientStore()
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = AmbientStore()
        ret_dict = {}
        for key in template._vars_:
            ret_dict[key] = getattr(template, key)
        return ret_dict

    def to_dict(self):
        ret_dict = {}
        for key in self._vars_:
            ret_dict[key] = getattr(self, key)
        return ret_dict

    def copy(self):
        c = AmbientStore()
        for key in c._vars_:
            setattr(c, key, getattr(self, key).copy())
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)

    @staticmethod
    def _validation_error(var_name, err_message):
        raise UserInputError(f"Invalid ambient value: {var_name} {err_message}")

    def validate(self):
        for vname in self._vars_:
            vsubset = self.get(vname)
            vname = vname.replace("_", " ")
            if not vsubset:
                self._validation_error(vname, "is empty")
            if not isinstance(vsubset, AmbientStoreSubset):
                self._validation_error(vname, "is invalid")
            if not is_truthy(vsubset.units, zero_is_true=True):
                self._validation_error(vname, "empty or invalid units")
            if not vsubset.extrapolation_sfc or not isinstance(vsubset.extrapolation_sfc, Interpolation):
                self._validation_error(vname, "empty or invalid extrapolation sfc")
            if not vsubset.extrapolation_btm or not isinstance(vsubset.extrapolation_btm, Interpolation):
                self._validation_error(vname, "empty or invalid extrapolation btm")
            if not vsubset.from_time_series:
                continue
            if not is_truthy(vsubset.ts_depth_units, zero_is_true=True):
                self._validation_error(vname, "time series units is empty or invalid")
            self.set(vname, convert_float(
                value=vsubset.ts_increment,
                allow_zero=False,
                allow_negative=False,
                error_handler=lambda msg: self._validation_error(f"{vname} time increment is", msg)
            ))
