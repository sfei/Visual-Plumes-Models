from ..ambient.Ambient import Ambient
from .AmbientTimeseries import AmbientTimeseries
from .DiffuserTimeseries import DiffuserTimeseries
from .TimeseriesStore import TimeseriesStore


class DiffuserTSVars:
    def __init__(self):
        self.depth         = None  # (formerly vdep) port depth
        self.effluent_flow = None  # (formerly vflo) effluent flow
        self.salinity      = None  # (formerly vsal) salinity
        self.temperature   = None  # (formerly vtem) temperature
        self.concentration = None  # (formerly vpol) effluent concentration
        self._vars_        = list(vars(self).keys())


class AmbientTSVars:
    def __init__(self):
        self.current_speed = None  # (formerly `vaspd`) current speed (flat component only -- no vertical component assumed)
        self.current_dir   = None  # (formerly `vadir`) current direction in degrees (azimuthal only -- no vertical component assumed)
        self.salinity      = None  # (formerly `vasal`) ambient salinity
        self.temperature   = None  # (formerly `vatem`) ambient temperature
        self.bg_conc       = None  # (formerly `vapol`) background concentration
        self.decay_rate    = None  # (formerly `varat`) pollutant decay rate
        self.ff_velocity   = None  # (formerly `vafar`) far-field velocity
        self.ff_dir        = None  # (formerly `vafad`) far-field current direction
        self.ff_diff_coeff = None  # (formerly `vadis`) far-field diffusion coefficient
        self._vars_        = list(vars(self).keys())


class TimeseriesHandler:

    def __init__(self):
        self.start_time     = 0.0  # (formerly antime) diffuser start time
        self.end_time       = 0.0  # (formerly abtime) diffuser end time
        self.time_increment = 0.0  # (formerly tinc) time increment
        self.units          = TimeseriesStore()
        self.diffuser       = DiffuserTSVars()
        self.ambient        = AmbientTSVars()

    def _forward(self, ts_handler, target_time_secs):
        # loop timeseries through values until meeting target time (in seconds)
        if target_time_secs < ts_handler.time():
            ts_handler.reset()
        while target_time_secs >= ts_handler.next_time():
            ts_handler.next()

    def set_time(self, target_time_secs):
        # move all timeseries to target time
        for vkey in self.diffuser._vars_:
            ts_handler = getattr(self.diffuser, vkey)
            if ts_handler:
                assert isinstance(ts_handler, DiffuserTimeseries)
                self._forward(ts_handler, target_time_secs)
        for vkey in self.ambient._vars_:
            ts_handler = getattr(self.ambient, vkey)
            if ts_handler:
                assert isinstance(ts_handler, AmbientTimeseries)
                self._forward(ts_handler, target_time_secs)

    def get_diffuser_info(self):
        diffuser_ts_info = {}
        for vkey in self.diffuser._vars_:
            ts_handler = getattr(self.diffuser, vkey)
            if not ts_handler:
                continue
            assert isinstance(ts_handler, DiffuserTimeseries)
            diffuser_ts_info[vkey] = {
                'seconds': ts_handler._time,
                'length':  ts_handler._length,
                'store':   ts_handler._store
            }
        return diffuser_ts_info

    def get_ambient_info(self):
        ambient_ts_info = {}
        for vkey in self.ambient._vars_:
            ts_handler = getattr(self.ambient, vkey)
            if not ts_handler:
                continue
            assert isinstance(ts_handler, AmbientTimeseries)
            ambient_ts_info[vkey] = {
                'seconds': ts_handler._time,
                'length':  ts_handler._length,
                'depths':  ts_handler._depths,
                'store':   ts_handler._store
            }
        return ambient_ts_info

    def get_diffuser(self):
        diffuser_ts_values = {}
        for vkey in self.diffuser._vars_:
            ts_handler = getattr(self.diffuser, vkey)
            if not ts_handler:
                continue
            assert isinstance(ts_handler, DiffuserTimeseries)
            # only one value per line makes sense to me
            diffuser_ts_values[vkey] = (ts_handler.values()[0], ts_handler._at)
        return diffuser_ts_values

    def get_ambient(self):
        # adapted from procedure fwdambfiles in main.pas
        ambient_ts_stacks = {}
        ambient_ts_indices = {}
        for vkey in self.ambient._vars_:
            ts_handler = getattr(self.ambient, vkey)
            if not ts_handler:
                continue
            assert isinstance(ts_handler, AmbientTimeseries)
            # convert to ambient stacks (all other values as None type, because any variable could have different depth levels)
            # TODO(?): is it a waste to use the Ambient object with bunch of empty variables or benefit to using common class?
            ambient_stack = []
            for val, z in zip(ts_handler.values(), ts_handler.depths()):
                ambient = Ambient()
                ambient.z = z
                ambient.set(vkey, val)
                ambient_stack.append(ambient)
            ambient_ts_stacks[vkey] = ambient_stack
            ambient_ts_indices[vkey] = ts_handler._at
        return ambient_ts_stacks, ambient_ts_indices
