from .. import units


class TimeseriesStore:
    def __init__(self):
        self.start_time     = units.Time.SECONDS
        self.end_time       = units.Time.SECONDS
        self.time_increment = units.Time.SECONDS
