from .AbstractTimeseries import AbstractTimeseries


class DiffuserTimeseries(AbstractTimeseries):

    def __init__(self, filepath, ambient_store_info):
        super().__init__(filepath, ambient_store_info)

