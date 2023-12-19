import numpy as np


FALSEY_STRINGS = ("None", "False", "F")
FALSEY_STRINGS_LOWER = tuple(s.lower() for s in FALSEY_STRINGS)


def num_format(num):
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
    # fastest shift, x5 faster than np.roll(), courtesy of https://stackoverflow.com/a/42642326
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
            value_num = float(value)
            return _is_truthy_number(value_num, zero_is_true, nan_is_true)
        except ValueError:
            pass
        # check false-y string messages
        if not case_sensitive:
            value = value.lower()
            compare_with = FALSEY_STRINGS_LOWER
        else:
            compare_with = FALSEY_STRINGS
        return bool(value not in compare_with)
    # all other types
    return True if value else False


def _is_truthy_number(value, zero_is_true, nan_is_true):
    if zero_is_true and value == 0.0:
        return True
    if nan_is_true and np.isnan(value):
        return True
    return False
