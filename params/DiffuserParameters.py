class DiffuserParameters:

    def __init__(self):
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
