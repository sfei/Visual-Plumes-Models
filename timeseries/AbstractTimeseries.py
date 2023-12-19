import csv
from ..globals import UserInputError
from .. import units
from ..ambient.AmbientStore import AmbientStoreSubset
from ..params.DiffuserStore import DiffuserStoreSubset


class AbstractTimeseries:

    def __init__(self, filepath, store_info):
        assert isinstance(store_info, (AmbientStoreSubset, DiffuserStoreSubset))

        self._store     = store_info
        self._lines     = []
        self._length    = 0
        self._at        = 0
        self._time      = 0
        self._increment = 1
        self._next_time = 1

        with open(filepath, 'r', encoding='utf-8') as ts_file:
            reader = csv.reader(ts_file)
            i = 0
            for row in reader:
                i += 1
                try:
                    self._lines.append(tuple(float(n) for n in row))
                except ValueError:
                    raise UserInputError(f"Error in timeseries file: error reading numeric value (at line={i + 1})")
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
        return self._length

    def reset(self):
        self._at = 0
        self._time = 0
        self._next_time = self._increment

    def next(self, dont_increment_time=False):
        self._at += 1
        if self._at >= self._length:
            self._at = 0
        if not dont_increment_time:
            self._time = self._next_time
            self._next_time += self._increment

    def time(self):
        return self._time

    def next_time(self):
        return self._next_time

    def increment_secs(self):
        return self._increment

    def values(self):
        return self._lines[self._at]
