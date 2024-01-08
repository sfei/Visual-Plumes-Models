from ..globals import UserInputError
from .AbstractTimeseries import AbstractTimeseries


class AmbientTimeseries(AbstractTimeseries):

    def __init__(self, filepath, ambient_store_info):
        super().__init__(filepath, ambient_store_info)
        self._depths = None
        self._levels = 0

    def parse(self):
        super(AmbientTimeseries, self).parse()

        # first line is list of depth levels
        self._depths = self._lines.pop(0)
        self._levels = len(self._depths)
        self._length = len(self._lines)

        try:
            self._depths = tuple(float(n) for n in self._depths)
        except ValueError:
            raise UserInputError("Error in timeseries file: invalid depth value.")
        if not self._levels:
            raise UserInputError("Error in timeseries file: no depth levels given in first line.")
        if self._levels < 2:
            raise UserInputError("Error in timeseries file: at least two depth levels must be provided.")
        last_depth = None
        for depth in self._depths:
            if depth < 0:
                raise UserInputError("Error in timeseries file: negative depth or height value encountered.")
            if last_depth is not None:
                if self._store.z_is_depth:
                    if depth <= last_depth:
                        raise UserInputError("Error in timeseries file: depth layers must be given in order of increasing depth (or decreasing height)")
                elif depth >= last_depth:
                    raise UserInputError("Error in timeseries file: depth layers must be given in order of increasing depth (or decreasing height)")
            last_depth = depth

        for i, row in enumerate(self._lines):
            if len(row) != self._levels:
                raise UserInputError(f"Error in timeseries file: length of row does not match number of depth levels (at line={i+1})")

    def depths(self):
        if self._lines is None:
            raise Exception("Timeseries data attempt to use before parse() called")
        return self._depths
