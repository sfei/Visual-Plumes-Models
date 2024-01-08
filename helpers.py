import numpy as np
from .globals import UserInputError


FALSEY_STRINGS = ("No", "None", "null", "False", "F")
FALSEY_STRINGS_LOWER = tuple(s.lower() for s in FALSEY_STRINGS)


def num_format(num):
    ''' Convert numeric into string with formatting and precision handled via heuristic.
    Args:
        num: Input number value
    Returns: (str)
    '''
    abs_num = abs(num)
    if abs_num < 1e-10:
        return "0.0"
    if abs_num < 0.01:
        return f"{num:.2e}"
    if abs_num < 1.0:
        return f"{num:.3f}"
    if abs_num < 10:
        return f"{num:.2f}"
    if abs_num >= 100000:
        return f"{num:.2e}"
    if abs_num >= 1000:
        return f"{num:,.0f}"
    return f"{num:,.1f}"


def shift(arr, places, fill_value=np.nan):
    ''' Shift array values down (to the right, if places arg is positive), which loses values shifting out of the array
    range (they do not wrap) and fills in values in now 'empty' spaces as desired. Uses what was determined to be the
    fastest method (just creating a new array and copying in values), x5 faster than np.roll(), courtesy of
    https://stackoverflow.com/a/42642326
    Args:
        arr: The array to do the operation on.
        places: The number of spaces to shift the values. Can be negative to shift left.
        fill_value: (Optional, default=np.nan) The fill value for the new spaces after shifting.
    Returns: (np.array) Note that array is generally a copy, but may be the same array if places=0.
    '''
    if places == 0:
        # if no change, returns same array
        # for our current uses this is fine, but note inconsistency if it matters that return is always a copy
        return arr
    result = np.empty_like(arr)
    if places > 0:
        result[:places] = fill_value
        result[places:] = arr[:-places]
    else:
        result[places:] = fill_value
        result[:places] = arr[-places:]
    return result


def is_truthy(value, case_sensitive=False, skip_string_values=False, zero_is_true=False, nan_is_true=False):
    ''' Heueristic truthy checking that allows for different ways to evaluate numbers and strings. By default, numbers
    equal to zero or NaN are false, strings that can be converted to numbers are treated as such, and string values (case
    insensitive) equal to "No", "None", "null", "False", or "F" are false.
    Args:
        value: The value to check.
        case_sensitive: (Optional, default=False) Whether string checking is case sensitive.
        skip_string_values: (Optional, default=False) If true, all strings (unless blank/empty) are treated as true.
        zero_is_true: (Optional, default=False) If true, numbers and strings representing numbers that equal zero are true.
        nan_is_true: (Optional, default=False) If true, numbers and strings representing numbers that are NaN are true.
    Returns: (bool)
    '''
    # numeric handling
    if isinstance(value, (int, float)):
        return _is_truthy_number(value, zero_is_true, nan_is_true)
    # string handling
    if isinstance(value, str):
        # empty or only whitespace is always false
        value = value.strip()
        if not value:
            return False
        # if no special string handling, skip conversion attempts
        if skip_string_values:
            return True
        # convert to numeric if applicable
        try:
            return _is_truthy_number(float(value), zero_is_true, nan_is_true)
        except ValueError:
            pass
        # check false-y string messages
        if not case_sensitive:
            value = value.lower()
            return bool(value not in FALSEY_STRINGS_LOWER)
        return bool(value not in FALSEY_STRINGS)
    # all other types
    return True if value else False


def _is_truthy_number(value, zero_is_true, nan_is_true):
    if zero_is_true and value == 0.0:
        return True
    if nan_is_true and np.isnan(value):
        return True
    return bool(value)


def is_empty(value):
    ''' Return true if value is None or empty string.
    Args:
        value: The value to check
    Returns: (bool)
    '''
    # none type or empty string
    if value is None:
        return True
    if isinstance(value, str) and value.strip() == "":
        return True
    return False


def convert_float(value, allow_nan=False, allow_zero=True, allow_negative=True, error_handler=None):
    ''' Converts given value to a float/numeric type, with options for validation and constraints. By default, raises a
    UserInputError if it fails the validations or cannot be converted to numeric.
    Args:
        value: The value to convert.
        allow_nan: (Optional, default=False) If false, NaN values will fail validation.
        allow_zero: (Optional, default=True) If false, zero values will fail validation.
        allow_negative: (Optional, default=True) If false, negative values will fail validation.
        error_handler: (Optional) If provided, error message will be passed to this function instead of raising a UserInputError.
    Returns: (float)
    '''
    if is_empty(value):
        if error_handler:
            return error_handler("empty")
        raise UserInputError("empty value in numeric input")
    try:
        value = float(value)
    except TypeError:
        if error_handler:
            return error_handler("invalid type")
        raise UserInputError("invalid type in numeric input")
    except ValueError:
        if error_handler:
            return error_handler("invalid")
        raise UserInputError("invalid value in numeric input")
    if not allow_nan and np.isnan(value):
        if error_handler:
            return error_handler("NaN")
        raise UserInputError("invalid NaN value in numeric input")
    if not allow_zero and value == 0:
        if error_handler:
            return error_handler("zero")
        raise UserInputError("zero value in input")
    if not allow_negative and value < 0:
        if error_handler:
            return error_handler("negative")
        raise UserInputError("negative value in input")
    return value
