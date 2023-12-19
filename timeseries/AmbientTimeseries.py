from ..globals import UserInputError
from .AbstractTimeseries import AbstractTimeseries


class AmbientTimeseries(AbstractTimeseries):

    def __init__(self, filepath, ambient_store_info):
        super().__init__(filepath, ambient_store_info)

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

        for i, row in enumerate(self._lines):
            if len(row) != self._levels:
                raise UserInputError(f"Error in timeseries file: length of row does not match number of depth levels (at line={i+1})")

    def depths(self):
        return self._depths
