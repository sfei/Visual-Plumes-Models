import numpy as np
from ..helpers import convert_float
from ..globals import UserInputError
from ..ambient.Ambient import Ambient
from .AmbientTimeseries import AmbientTimeseries
from .DiffuserTimeseries import DiffuserTimeseries
from .TimeseriesStore import TimeseriesStore


class DiffuserTSVars:
    """ Basically just a struct for diffuser variables that can come as timeseries. """
    def __init__(self):
        self.depth         = None  # (formerly vdep) port depth
        self.effluent_flow = None  # (formerly vflo) effluent flow
        self.salinity      = None  # (formerly vsal) salinity
        self.temperature   = None  # (formerly vtem) temperature
        self.concentration = None  # (formerly vpol) effluent concentration
        self._vars_        = list(vars(self).keys())


class AmbientTSVars:
    """ Basically just a struct for ambient variables that can come as timeseries. """
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

    @staticmethod
    def _validation_error(var_name, err_message):
        raise UserInputError(f"Error in timeseries inputs: {var_name} {err_message}")

    def _validate(self, value, var_name, allow_nan=False, allow_zero=True, allow_negative=True):
        return convert_float(
            value=value,
            allow_nan=allow_nan,
            allow_zero=allow_zero,
            allow_negative=allow_negative,
            error_handler=lambda msg: self._validation_error(f"{var_name} is", msg)
        )

    def _validate_values(self, ts_handler, var_name, allow_zero=True, allow_negative=True):
        if not ts_handler:
            return
        ts_handler.parse()
        for vals in ts_handler._lines:
            for v in vals:
                if np.isnan(v):
                    self._validation_error(var_name, "is NaN")
                if not allow_zero and v == 0.0:
                    self._validation_error(var_name, "is zero")
                if not allow_negative and v < 0:
                    self._validation_error(var_name, "is negative")

    def validate(self):
        self.start_time     = self._validate(self.start_time, "start time", allow_negative=False)
        self.end_time       = self._validate(self.end_time, "end time", allow_negative=False)
        self.time_increment = self._validate(self.time_increment, "time increment", allow_zero=False, allow_negative=False)
        if self.end_time < self.start_time:
            self._validation_error("end time", "occurs before start time")

        self._validate_values(self.diffuser.depth, 'port depth', allow_negative=False)
        self._validate_values(self.diffuser.effluent_flow, 'effluent flowrate', allow_zero=False, allow_negative=False)
        self._validate_values(self.diffuser.salinity, 'effluent salinity')
        self._validate_values(self.diffuser.temperature, 'effluent temperature')
        self._validate_values(self.diffuser.concentration, 'effluent concentration', allow_negative=False)

        self._validate_values(self.ambient.current_speed, 'current speed')
        self._validate_values(self.ambient.current_dir, 'current direction')
        self._validate_values(self.ambient.salinity, 'ambient temperature')
        self._validate_values(self.ambient.temperature, 'ambient temperature')
        self._validate_values(self.ambient.bg_conc, 'background concentration', allow_negative=False)
        self._validate_values(self.ambient.decay_rate, 'decay rate', allow_negative=False)
        self._validate_values(self.ambient.ff_velocity, 'far-field velocity')
        self._validate_values(self.ambient.ff_dir, 'far-field direction')
        self._validate_values(self.ambient.ff_diff_coeff, 'far-field coefficient')

    def _forward(self, ts_handler, target_time_secs):
        """ Move timeseries datum to target time. Note that timeseries handlers will loop through data, if too short,
        assuming infinitely cyclical forwards in time.
        Args:
            ts_handler: The AmbientTimeseries or DiffuserTimeseries representing the datum.
            target_time_secs: The target time in seconds.
        """
        # loop timeseries through values until meeting target time (in seconds)
        if target_time_secs < 0:
            raise Exception("Invalid target time for timeseries (negative time).")
        if target_time_secs < ts_handler.time():
            ts_handler.reset()
        while target_time_secs >= ts_handler.next_time():
            ts_handler.next()

    def set_time(self, target_time_secs):
        """ Move all timeseries data to the targeted time.
        Args:
            target_time_secs: The target time in seconds.
        """
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
        """ Get dictionary (by var name as key) of timeseries info for available diffuser vars. Each value is subdict of
        'seconds', 'length', and 'store' for current time of datum, length of data/rows (i.e. data cycle length), and
        the relevant DiffuserStoreSubset with metadata. """
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
        """ Get dictionary (by var name as key) of timeseries info for available ambient vars. Each value is subdict of
        'seconds', 'length', and 'store' for current time of datum, length of data/rows (i.e. data cycle length), and
        the relevant AmbientStoreSubset with metadata. """
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
        """ Get dictionary of current diffuser values (by var name as key). """
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
        """ Get current ambient values. Returns pair of dictionaries (both with var name as key). First is the data 
        stack of Ambient objects by depth layer. Only those variables for existing timeseries data will be filled (rest
        will remain None). Second is the timeseries indices (for debugging purposes, if needed). """
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
