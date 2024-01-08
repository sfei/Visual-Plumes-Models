from ..helpers import convert_float
from ..globals import UserInputError


class AbstractParameters:
    """ Common functions for parameter classes to inherit. Getters, settings, and validation functions. """

    def __init__(self):
        pass

    def get(self, key):
        """ Quick, class method shortcut for getattr() on this object. """
        return getattr(self, key)

    def set(self, key, value):
        """ Quick, class method shortcut for setattr() on this object. """
        return setattr(self, key, value)

    @staticmethod
    def _validation_error(var_name, err_message):
        """ Simple function for raising UserInputError. Encapsulated here so can reference in lambda functions. """
        raise UserInputError(f"{var_name} {err_message}")

    def _validate(self, value, var_name, allow_nan=False, allow_zero=True, allow_negative=True, as_integer=False, default_if_none=None, min=None, max=None):
        """ Simple validation function for numeric input values.
        Args:
            value: The value to check.
            var_name: The value name, formatted for printing if an error occurs.
            allow_nan: (Optional, default=False) If false, NaN values will fail validation.
            allow_zero: (Optional, default=True) If false, zero values will fail validation.
            allow_negative: (Optional, default=True) If false, negative values will fail validation.
            as_integer: (Optional, default=False) If true, value will be converted to integer type (and checked as int).
            default_if_none: (Optional) If not None, None-type values will be converted to this default value.
            min: (Optional) If set, validates against this minimum value, inclusive.
            max: (Optional) If set, validates against this maximum value, inclusive.
        Returns: The validate value
        """
        if default_if_none is not None and value is None:
            return default_if_none
        value = convert_float(
            value=value,
            allow_nan=allow_nan,
            allow_zero=allow_zero,
            allow_negative=allow_negative,
            error_handler=lambda msg: self._validation_error(f"{var_name} is", msg)
        )
        if min is not None and value < min:
            raise UserInputError(f"{var_name} is below minimum ({min})")
        if max is not None and value > max:
            raise UserInputError(f"{var_name} is above maximum ({max})")
        if as_integer:
            int_value = int(value)
            if int_value != value:
                raise UserInputError(f"{var_name} is not integer")
            value = int_value
        return value

    def _validate_enum(self, value, var_name, enum_class):
        """ Simple validation for enumerate types. Simply makes sure it is not empty and a member of the Enum class.
        Args:
            value: The value to check.
            var_name: The value name, formatted for printing if an error occurs.
            enum_class: An object that inherits class Enum.
        Returns: The validate value
        """
        if not value or not isinstance(value, enum_class):
            raise UserInputError(f"{var_name} is null or invalid")
        return value
