class Ambient:
    """ Struct-like class for holding ambient values at a specific depth level. """
    
    def __init__(self):
        # (oc) ambient file order
        self.z              = None  # (formerly `valen`) z-coordinate
        self.depth          = None  # z converted to depth (because z might be given in height from bottom)
        self.current_speed  = None  # (formerly `vaspd`) current speed (flat component only -- no vertical component assumed)
        self.current_dir    = None  # (formerly `vadir`) current direction in degrees (azimuthal only -- no vertical component assumed)
        self.salinity       = None  # (formerly `vasal`) ambient salinity
        self.temperature    = None  # (formerly `vatem`) ambient temperature
        self.bg_conc        = None  # (formerly `vapol`) background concentration
        self.decay_rate     = None  # (formerly `varat`) pollutant decay rate
        self.ff_velocity    = None  # (formerly `vafar`) far-field velocity
        self.ff_dir         = None  # (formerly `vafad`) far-field current direction
        self.ff_diff_coeff  = None  # (formerly `vadis`) far-field diffusion coefficient
        self.density        = 0.0   # (formerly `vaden`) ambient seawater density (calculated)
        self.kt             = 0.0   # (formerly `kt`) pollutant decay rate after conversion
        self._vars_ = list(vars(self).keys())
        self._input_vars_ = self._vars_[:]
        self._input_vars_.remove('z')
        self._input_vars_.remove('density')
        self._input_vars_.remove('kt')

    @staticmethod
    def from_dict(input_dict):
        c = Ambient()
        for key in c._vars_:
            if key in input_dict:
                setattr(c, key, input_dict[key])
        return c

    @staticmethod
    def as_dict():
        template = Ambient()
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
        c = Ambient()
        for key in c._vars_:
            setattr(c, key, getattr(self, key))
        return c

    def get(self, key):
        return getattr(self, key)

    def set(self, key, value):
        return setattr(self, key, value)