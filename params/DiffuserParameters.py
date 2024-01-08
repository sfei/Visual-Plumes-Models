from ..globals import UserInputError
from .AbstractParameters import AbstractParameters


class DiffuserParameters(AbstractParameters):

    def __init__(self):
        super().__init__()

        self.diameter            = 0.0  # (formerly vdia) diameter
        # conduit depth? I think in PDS it reuses vhor (port elevation) var
        # self.port_elevation      = 0.0  # (formerly ????) port elevation
        self.offset_x            = 0.0  # (formerly vorx) diffuser x-position offset
        self.offset_y            = 0.0  # (formerly vory) diffuser y-position offset
        self.vertical_angle      = 0.0  # (formerly vang) vertical angle (in degrees)
        self.horizontal_angle    = 0.0  # (formerly vhor) horizontal angle
        self.num_ports           = 0.0  # (formerly vnum) number of ports
        self.port_spacing        = 0.0  # (formerly vspc) port spacing
        # time parameters (antime/abtime/tinc) moved to timeseries handler
        # these don't seem to be used in UM3
        self.start_time          = 0.0  # (formerly antime) diffuser start time
        self.end_time            = 0.0  # (formerly abtime) diffuser end time
        self.time_increment      = 0.0  # (formerly tinc) time increment
        # continuing..
        self.acute_mixing_zone   = 0.0  # (formerly vmxz) acute mixing zone
        # self.chronic_mixing_zone = 0.0  # (formerly ????) chronic mixing zone (region of interest in PDS)
        self.isopleth            = 0.0  # (formerly viso) isopleth (can be salinity, concentration, temperature, or speed)
        # (oc) diffuser file order
        self.depth               = 0.0  # (formerly vdep) port depth
        self.effluent_flow       = 0.0  # (formerly vflo) effluent flow
        self.salinity            = 0.0  # (formerly vsal) salinity
        self.temperature         = 0.0  # (formerly vtem) temperature
        self.concentration       = 0.0  # (formerly vpol) effluent concentration
        self._vars_ = list(vars(self).keys())

    @staticmethod
    def from_dict(input_dict):
        c = DiffuserParameters()
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = DiffuserParameters()
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
        c = DiffuserParameters()
        for key in c._vars_:
            setattr(c, key, getattr(self, key))
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)

    @staticmethod
    def _check_needs_value(name, store):
        if not store:
            return True
        meta = store.get(name)
        return not meta or not meta.from_time_series

    def validate(self, store=None):
        self.diameter          = self._validate(self.diameter, "diffuser diameter", allow_zero=False, allow_negative=False)
        self.offset_x          = self._validate(self.offset_x, "diffuser offset x", default_if_none=0)
        self.offset_y          = self._validate(self.offset_y, "diffuser offset y", default_if_none=0)
        self.vertical_angle    = self._validate(self.vertical_angle, "diffuser vertical angle", default_if_none=0)
        self.horizontal_angle  = self._validate(self.horizontal_angle, "diffuser horizontal angle", default_if_none=0)
        self.num_ports         = self._validate(self.num_ports, "number of ports", allow_zero=False, allow_negative=False, as_integer=True)
        if self.num_ports > 1:
            self.port_spacing  = self._validate(self.port_spacing, "port spacing", allow_zero=False, allow_negative=False)
        self.start_time        = self._validate(self.start_time, "start time", allow_negative=False, default_if_none=0)
        self.end_time          = self._validate(self.end_time, "end time", allow_negative=False, default_if_none=0)
        self.time_increment    = self._validate(self.time_increment, "time increment", allow_negative=False, default_if_none=0)
        if self.start_time != self.end_time and self.time_increment == 0:
            raise UserInputError("Null or invalid time increment for given start/end times")
        self.acute_mixing_zone = self._validate(self.acute_mixing_zone, "mixing zone distance")
        self.isopleth          = self._validate(self.isopleth, "isopleth value")
        self.depth             = self._validate(self.depth, "port depth")
        # below allow timeseries as inputs
        if self._check_needs_value('effluent_flow', store):
            self.effluent_flow = self._validate(self.effluent_flow, "effluent flow", allow_zero=False)
        if self._check_needs_value('salinity', store):
            self.salinity      = self._validate(self.salinity, "effluent salinity")
        if self._check_needs_value('temperature', store):
            self.temperature   = self._validate(self.temperature, "effluent temperature")
        if self._check_needs_value('concentration', store):
            self.concentration = self._validate(self.concentration, "effluent concentration", allow_negative=False)
