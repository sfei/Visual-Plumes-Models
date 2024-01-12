from .. import units
from ..helpers import convert_float, is_truthy
from ..globals import UserInputError


class DiffuserStoreSubset:
    """ Store class for desc/metadata about a specific diffuser parameter """

    def __init__(self, default_unit=units.Unitless.UNITLESS):
        # (oc) ambient file order
        self.from_time_series = False         # (formerly ffin) time-series file
        self.ts_increment     = 0             # (formerly fdel) time-series time increment (assumed in hours)
        # self.ts_period        = 0             # (formerly fper) time-series cycling period (assumed in hours)
        # self.ts_file_units    = default_unit  # (formerly fufi) time-series measurement unit
        # TODO: doesn't need units? putting it in because feels like it should have
        self.units            = default_unit  # (formerly futa) measurement unit
        self._vars_           = list(vars(self).keys())

    @staticmethod
    def from_dict(input_dict):
        c = DiffuserStoreSubset()
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = DiffuserStoreSubset()
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
        c = DiffuserStoreSubset()
        for key in c._vars_:
            setattr(c, key, getattr(self, key))
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)


class DiffuserStore:
    """ Store class for desc/metadata about diffuser parameters """

    def __init__(self):
        self.diameter            = DiffuserStoreSubset(units.Length.METERS)          # (formerly vdia) diameter
        # conduit depth? I think in PDS it reuses vhor (port elevation) var
        # self.port_elevation      = DiffuserStoreSubset(units.Length.METERS)        # (formerly ????) port elevation
        self.offset_x            = DiffuserStoreSubset(units.Length.METERS)          # (formerly vorx) diffuser x-position offset
        self.offset_y            = DiffuserStoreSubset(units.Length.METERS)          # (formerly vory) diffuser y-position offset
        self.vertical_angle      = DiffuserStoreSubset(units.Angle.DEGREES)          # (formerly vang) vertical angle (in degrees)
        self.horizontal_angle    = DiffuserStoreSubset(units.Angle.DEGREES)          # (formerly vhor) horizontal angle
        self.num_ports           = DiffuserStoreSubset(units.Unitless.UNITLESS)      # (formerly vnum) number of ports
        self.port_spacing        = DiffuserStoreSubset(units.Length.METERS)          # (formerly vspc) port spacing
        # these don't seem to be used in UM3
        self.start_time          = DiffuserStoreSubset(units.Time.HOURS)             # (formerly antime) diffuser start time
        self.end_time            = DiffuserStoreSubset(units.Time.HOURS)             # (formerly abtime) diffuser end time
        self.time_increment      = DiffuserStoreSubset(units.Time.HOURS)             # (formerly tinc) time increment
        # continuing..
        self.acute_mixing_zone   = DiffuserStoreSubset(units.Length.METERS)          # (formerly vmxz) acute mixing zone
        self.isopleth            = DiffuserStoreSubset(units.Isopleth.CONCENTRATION) # (formerly viso) isopleth (can be salinity, concentration, temperature, or speed)
        # self.chronic_mixing_zone = DiffuserStoreSubset(units.Length.METERS)        # (formerly viso) chronic mixing zone (region of interest in PDS)
        self.depth               = DiffuserStoreSubset(units.Length.METERS)                        # (formerly vdep) port depth
        self.effluent_flow       = DiffuserStoreSubset(units.FlowRate.CUBIC_METERS_PER_SECOND)     # (formerly vflo) effluent flow
        self.salinity            = DiffuserStoreSubset(units.Salinity.PRACTICAL_SALINITY_UNITS)    # (formerly vsal) salinity
        self.temperature         = DiffuserStoreSubset(units.Temperature.CELSIUS)                  # (formerly vtem) temperature
        self.concentration       = DiffuserStoreSubset(units.Concentration.KILOGRAM_PER_KILOGRAM)  # (formerly vpol) effluent concentration
        self._vars_ = list(vars(self).keys())

    @staticmethod
    def from_dict(input_dict):
        c = DiffuserStore()
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = DiffuserStore()
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
        c = DiffuserStore()
        for key in c._vars_:
            setattr(c, key, getattr(self, key).copy())
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)

    @staticmethod
    def _validation_error(var_name, err_message):
        raise UserInputError(f"Invalid diffuser store value: {var_name} {err_message}")

    def validate(self):
        for vname in self._vars_:
            vsubset = self.get(vname)
            vname = vname.replace("_", " ")
            if not vsubset:
                self._validation_error(vname, "is empty")
            if not isinstance(vsubset, DiffuserStoreSubset):
                self._validation_error(vname, "is invalid")
            if not is_truthy(vsubset.units, zero_is_true=True):
                self._validation_error(vname, "empty or invalid units")
            if not vsubset.from_time_series:
                continue
            self.set(vname, convert_float(
                value=vsubset.ts_increment,
                allow_zero=False,
                allow_negative=False,
                error_handler=lambda msg: self._validation_error(f"{vname} time increment is", msg)
            ))
