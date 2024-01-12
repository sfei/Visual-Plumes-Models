import math
from .vectors import DEGREES_TO_RADIAN


COS_20_DEGREES = math.cos(20*DEGREES_TO_RADIAN)

# TODO: would NaN be better? (don't think so, seems to need to be big number, but why so specific?)
missing = round((math.pi-1)*1e8-1)

# maximum number of steps for plume model are hard-coded
maxsteps = 2000

class UserInputError(Exception):
    """ Custom exception class for user input error """
    pass
