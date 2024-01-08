import csv
from ..globals import UserInputError
from .. import units
from ..ambient.AmbientStore import AmbientStoreSubset
from ..params.DiffuserStore import DiffuserStoreSubset


class AbstractTimeseries:
    """ This is abstract parent class for a singular timeseries dataset handler. Defines functions for getting data and
    changing/incrementing the time. """

    def __init__(self, filepath, store_info):
        """
        Args:
            filepath: String filepath to the data. Assumed CSV.
            store_info: Either AmbientStoreSubset or DiffuserStoreSubset with metadata for this datum.
        """
        assert isinstance(store_info, (AmbientStoreSubset, DiffuserStoreSubset))
        self._store     = store_info  # the store/metadata for this type, which can be for ambient or diffuser
        self._raw       = None        # the raw timeseries data by lines/rows
        self._lines     = None        # the timeseries data by lines/rows after converting to numerics
        self._length    = 0           # the amount of timeseries data, AKA number of rows, AKA time increments
        self._at        = 0           # index of current row/time-increment position in this dataset
        self._time      = 0           # the equivalent time (in seconds) for the current index given by self._at
        self._increment = 1           # the time increment (in seconds)
        self._next_time = 1           # the equivalent time (in seconds) for the next increment
        with open(filepath, 'r', encoding='utf-8') as ts_file:
            # just read to raw for now, handling parsing in function that can be called within exception handling
            reader = csv.reader(ts_file)
            self._raw = tuple(row for row in reader)

    def parse(self):
        if self._lines is not None:
            raise Exception("Called parse() on timeseries instance twice.")

        self._lines = []
        for i, row in enumerate(self._raw):
            try:
                self._lines.append(tuple(float(n) for n in row))
            except ValueError:
                raise UserInputError(f"Error in timeseries file: error reading numeric value (at line={i+1})")
        if not len(self._lines):
            raise UserInputError("Error in timeseries file: empty file.")
        last_row = self._lines[-1]
        while not len(last_row):
            self._lines.pop()
            last_row = self._lines[-1]
        self._length = len(self._lines)
        if self._length < 2:
            raise UserInputError("Error in timeseries file: empty file or not enough lines.")

        self._increment = units.Time.convert(self._store.ts_increment, units.Time.HOURS, units.Time.SECONDS)
        self._next_time = self._increment

    def length(self):
        """ Get the length of the data, which is the number of rows and equivalent to the number of time increments. """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._length

    def reset(self):
        """ Reset index/increment to the start (time=0). """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        self._at = 0
        self._time = 0
        self._next_time = self._increment

    def next(self, dont_increment_time=False):
        """ Move current index/time by one increment forwards. If past the length of data available, will loop back to
        the beginning, assuming infinitely cyclical forwards in time.
        Args:
            dont_increment_time: (Optional) If true, increments the index without moving the current time. Not recommended unless fixing some syncing issue.
        """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        self._at += 1
        if self._at >= self._length:
            self._at = 0
        if not dont_increment_time:
            self._time = self._next_time
            self._next_time += self._increment

    def time(self):
        """ Get the current time (in seconds) of the data position. """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._time

    def next_time(self):
        """ Get the next time (in seconds) of the data position after one increment. """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._next_time

    def increment_secs(self):
        """ Get the time increment length (in seconds). """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._increment

    def values(self):
        """ Get the current values of this datum at its current index/time-increment. """
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._lines[self._at]
