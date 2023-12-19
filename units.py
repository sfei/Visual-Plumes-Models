import math
from .ambient import calculations as ambient


class Units:
    @staticmethod
    def convert(value, ufrom, uto=1):
        if ufrom == uto:
            return value
        return value*(uto/ufrom)
    @staticmethod
    def label(uindex=1):
        return ""
    @staticmethod
    def validate(value, inunits):
        return True


class Unitless(Units):
    UNITLESS = 1
    @staticmethod
    def convert(value, ufrom=1, uto=1):
        return value


class Length(Units):
    METERS      = 1
    CENTIMETERS = 100
    FEET        = 3.2808
    INCHES      = 39.37
    FATHOMS     = 0.5468
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Length.METERS:
                return "m"
            case Length.CENTIMETERS:
                return "cm"
            case Length.FEET:
                return "ft"
            case Length.INCHES:
                return "in"
            case Length.FATHOMS:
                return "ftm"
            case _:
                return ""


class Mass(Units):
    KILOGRAMS = 1
    GRAMS     = 1000
    POUNDS    = 2.20462
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Mass.KILOGRAMS:
                return "kg"
            case Mass.GRAMS:
                return "gr"
            case Mass.POUNDS:
                return "lbs"
            case _:
                return ""
    # TODO: does validation need to ensure non-negative values? doesn't make sense for mass but as weight it could..


class Speed(Units):
    METERS_PER_SECOND      = 1
    CENTIMETERS_PER_SECOND = 100
    FATHOMS_PER_SECOND     = 1/0.514444
    KNOTS                  = 1.94384
    MILES_PER_HOUR         = 2.2369363
    FEET_PER_SECOND        = 3.2808
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Speed.METERS_PER_SECOND:
                return "m/s"
            case Speed.CENTIMETERS_PER_SECOND:
                return "cm/s"
            case Speed.FATHOMS_PER_SECOND:
                return "ftm/s"
            case Speed.KNOTS:
                return "kts"
            case Speed.MILES_PER_HOUR:
                return "mph"
            case Speed.FEET_PER_SECOND:
                return "ft/s"
            case _:
                return ""


class Angle(Units):
    DEGREES   = 1
    RADIANS   = math.pi/180
    N_DEGREES = -1
    N_RADIANS = -math.pi/180
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Angle.DEGREES:
                return "deg"
            case Angle.N_DEGREES:
                return "s-deg"
            case Angle.RADIANS:
                return "rad"
            case Angle.N_RADIANS:
                return "s-rad"
            case _:
                return ""
    @staticmethod
    def fix(value, units):
        # fit within min-max (0-360 for degrees)
        incr = 360 if abs(units) == 1 else 2.0*math.pi
        if value < 0:
            while value < 0:
                value += incr
        else:
            while value >= incr:
                value -= incr
        return value
    @staticmethod
    def convert(value, ufrom, uto=DEGREES):
        if ufrom == uto:
            return Angle.fix(value, ufrom)
        match ufrom:
            case Angle.N_DEGREES:
                value = 90 - value
            case Angle.N_RADIANS:
                value = math.pi*0.5 - value
        value = super(Angle, Angle).convert(Angle.fix(value, ufrom), abs(ufrom), abs(uto))
        match uto:
            case Angle.N_DEGREES:
                value = 90 - value
            case Angle.N_RADIANS:
                value = math.pi*0.5 - value
        return Angle.fix(value, uto)


class Salinity(Units):
    PRACTICAL_SALINITY_UNITS  = 1
    MILLIMHO_PER_CENTIMETER   = 2
    KILOGRAMS_PER_CUBIC_METER = 3
    SIGMA_T                   = 4
    POUNDS_PER_CUBIC_FOOT     = 5
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Salinity.PRACTICAL_SALINITY_UNITS:
                return "PSU"
            case Salinity.MILLIMHO_PER_CENTIMETER:
                return "mmhos/cm"
            case Salinity.KILOGRAMS_PER_CUBIC_METER:
                return "kg/m3"
            case Salinity.SIGMA_T:
                return "sigmaT"
            case Salinity.POUNDS_PER_CUBIC_FOOT:
                return "lbs/ft3"
            case _:
                return ""
    @staticmethod
    def convert(value, ufrom, uto=PRACTICAL_SALINITY_UNITS, celsius=None, at_equilibrium=None, depth=None):
        if value < 0:
            raise Exception("Salinity < 0 (arises when units are converted and salinity is near zero)")
        if ufrom == uto:
            return value
        do_convert = False
        match ufrom:
            case Salinity.MILLIMHO_PER_CENTIMETER:
                if celsius is None:
                    raise Exception("Salinity conversion requires celsius parameter.")
                value, is_good = ambient.mho_salinity(value, celsius)
                if not is_good:
                    raise Exception("Out of range, re-establish units with \"Change label only\"\nIf multiple cases, some values may have been converted.")
            case Salinity.KILOGRAMS_PER_CUBIC_METER:
                value -= 1000
                do_convert = True
            case Salinity.SIGMA_T:
                do_convert = True
            case Salinity.POUNDS_PER_CUBIC_FOOT:
                value = value/0.062427961 - 1000
                do_convert = True
        if do_convert:
            if celsius is None:
                raise Exception("Salinity conversion requires celsius parameter.")
            if at_equilibrium is None:
                raise Exception("Salinity conversion requires at_equilibrium parameter.")
            value = ambient.salinity(celsius, value, at_equilibrium=at_equilibrium, depth=depth)
        match uto:
            case Salinity.MILLIMHO_PER_CENTIMETER:
                if celsius is None:
                    raise Exception("Salinity conversion requires celsius parameter.")
                value, is_good = ambient.mho_salinity(value, celsius)
                if not is_good:
                    raise Exception("Out of range, re-establish units with \"Change label only\"\nIf multiple cases, some values may have been converted.")
                return value
            case Salinity.KILOGRAMS_PER_CUBIC_METER:
                return value + 1000
            case Salinity.SIGMA_T:
                return ambient.salinity(celsius, value, at_equilibrium=at_equilibrium, depth=depth)
            case Salinity.POUNDS_PER_CUBIC_FOOT:
                return (value + 1000)*0.062427961
            case _:
                return value
    @staticmethod
    def validate(value, inunits, celsius=None, at_equilibrium=None, depth=None):
        if inunits == Salinity.SIGMA_T:
            value = Salinity.convert(
                value,
                inunits,
                Salinity.PRACTICAL_SALINITY_UNITS,
                celsius=celsius,
                at_equilibrium=at_equilibrium,
                depth=depth
            )
        return value >= 0


class Temperature(Units):
    CELSIUS    = 1
    KELVIN     = 2
    FAHRENHEIT = 3
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Temperature.CELSIUS:
                return "°C"
            case Temperature.KELVIN:
                return "°K"
            case Temperature.FAHRENHEIT:
                return "°F"
            case _:
                return ""
    @staticmethod
    def convert(value, ufrom, uto=CELSIUS):
        if ufrom == uto:
            return value
        match ufrom:
            case Temperature.FAHRENHEIT:
                value = 5.0/9.0*(value - 32.0)
            case Temperature.KELVIN:
                value -= 273.15
        match uto:
            case Temperature.CELSIUS:
                return value
            case Temperature.KELVIN:
                return value + 273.15
            case Temperature.FAHRENHEIT:
                return 9.0/5.0*value + 32
        return value
    @staticmethod
    def validate(value, inunits):
        value = Temperature.convert(value, inunits, Temperature.KELVIN)
        return value >= 0


class Concentration(Units):
    KILOGRAM_PER_KILOGRAM = 1
    PARTS_PER_MILLION     = 1e6
    PARTS_PER_BILLION     = 1e9
    PERCENT               = 100
    COLONIES_PER_100ML    = -1
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Concentration.KILOGRAM_PER_KILOGRAM:
                return "kg/kg"
            case Concentration.PARTS_PER_MILLION:
                return "ppm"
            case Concentration.PARTS_PER_BILLION:
                return "ppb"
            case Concentration.PERCENT:
                return "%"
            case Concentration.COLONIES_PER_100ML:
                return "col/dl"
            case _:
                return ""
    @staticmethod
    def convert(value, ufrom, uto=KILOGRAM_PER_KILOGRAM):
        return super(Concentration, Concentration).convert(value, abs(ufrom), abs(uto))
    @staticmethod
    def validate(value, inunits):
        return value >= 0


class DecayRate(Units):
    PER_SECOND  = 1
    PER_DAY     = 86400
    T90_HOUR    = -2.30258509/3600
    LY_PER_HOUR = -1
    PER_HOUR    = 3600
    @staticmethod
    def label(uindex=1):
        match uindex:
            case DecayRate.PER_SECOND:
                return "s-1"
            case DecayRate.PER_DAY:
                return "d-1"
            case DecayRate.T90_HOUR:
                return "T90h"
            case DecayRate.LY_PER_HOUR:
                return "ly/hr"
            case DecayRate.PER_HOUR:
                return "hr-1"
            case _:
                return ""
    @staticmethod
    def convert(value, ufrom, uto=PER_SECOND, bacteria_model=None, celsius=None, psu=None, depth=None):
        if ufrom == uto:
            return value
        match ufrom:
            case DecayRate.PER_SECOND:
                pass
            case DecayRate.T90_HOUR:
                value = 1e10 if value == 0 else ufrom / value
            case DecayRate.LY_PER_HOUR:
                if bacteria_model is None:
                    raise Exception("Decay rate conversion from ly/hr requires bacteria_model parameter.")
                if celsius is None:
                    raise Exception("Decay rate conversion from ly/hr requires celsius parameter.")
                if psu is None:
                    raise Exception("Decay rate conversion from ly/hr requires psu parameter.")
                if depth is None:
                    raise Exception("Decay rate conversion from ly/hr requires depth parameter.")
                value = ambient.mancini(
                    model=bacteria_model,
                    arg=value,
                    salinity=psu,
                    temperature=celsius,
                    depth=depth,
                    topersec=True
                )
            case _:
                value /= ufrom
        match uto:
            case DecayRate.T90_HOUR:
                return 1.0 / 1e10 if value == 0 else uto / value
            case DecayRate.LY_PER_HOUR:
                if bacteria_model is None:
                    raise Exception("Decay rate conversion from ly/hr requires bacteria_model parameter.")
                if celsius is None:
                    raise Exception("Decay rate conversion from ly/hr requires celsius parameter.")
                if psu is None:
                    raise Exception("Decay rate conversion from ly/hr requires psu parameter.")
                return ambient.mancini(
                    model=bacteria_model,
                    arg=value,
                    salinity=psu,
                    temperature=celsius,
                    depth=depth,
                    topersec=False
                )
            case _:
                return value*uto
    @staticmethod
    def validate(value, inunits):
        return value >= 0


class Density(Units):
    KILOGRAMS_PER_CUBIC_METER  = 1
    SIGMA_T                    = -1
    GRAMS_PER_CUBIC_CENTIMETER = 0.001
    POUNDS_PER_CUBIC_FOOT      = 0.062428
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Density.KILOGRAMS_PER_CUBIC_METER:
                return "kg/m3"
            case Density.SIGMA_T:
                return "sigmaT"
            case Density.GRAMS_PER_CUBIC_CENTIMETER:
                return "g/cm3"
            case Density.POUNDS_PER_CUBIC_FOOT:
                return "lbs/ft3"
            case _:
                return ""
    @staticmethod
    def convert(value, ufrom, uto=KILOGRAMS_PER_CUBIC_METER):
        if ufrom == uto:
            return value
        if ufrom == Density.SIGMA_T:
            value = 1000.0 + value
        else:
            value = super(Density, Density).convert(value, ufrom, 1)
        return value - 1000.0 if uto == Density.SIGMA_T else value
    @staticmethod
    def validate(value, inunits):
        return value >= 0


class Time(Units):
    SECONDS = 1
    MINUTES = 1/60.0
    HOURS   = 1/3600.0
    DAYS    = 1/86400.0
    @staticmethod
    def label(uindex=1):
        match uindex:
            case Time.SECONDS:
                return "s"
            case Time.MINUTES:
                return "min"
            case Time.HOURS:
                return "hr"
            case Time.DAYS:
                return "d"
            case _:
                return ""


class FlowRate(Units):
    CUBIC_METERS_PER_SECOND = 1
    MEGALITERS_PER_DAY      = 86.4
    MEGAGALLONS_PER_DAY     = 22.824465
    CUBIC_FEET_PER_SECOND   = 35.314667
    BARRELS_PER_DAY         = 86400/0.15898284
    @staticmethod
    def label(uindex=1):
        match uindex:
            case FlowRate.CUBIC_METERS_PER_SECOND:
                return "m3/s"
            case FlowRate.MEGALITERS_PER_DAY:
                return "MLD"
            case FlowRate.MEGAGALLONS_PER_DAY:
                return "MGD"
            case FlowRate.CUBIC_FEET_PER_SECOND:
                return "ft3/s"
            case FlowRate.BARRELS_PER_DAY:
                return "bbl/d"
            case _:
                return ""
    @staticmethod
    def validate(value, inunits):
        return value >= 0


# special case, can't be converted of itself, just defines what class of measurement it is in
class Isopleth(Units):
    CONCENTRATION = -10
    SALINITY      = -20
    TEMPERATURE   = -30
    SPEED         = -40
    @staticmethod
    def convert(value, ufrom, uto=1):
        return ufrom
    @staticmethod
    def label(uindex=1):
        return ""


class EddyDiffusivity(Units):
    DIFFUSIVITY = 1
    @staticmethod
    def convert(value, ufrom, uto=1):
        return value
    @staticmethod
    def label(uindex=1):
        return "m0.67/s2"


def from_var_name(vkey):
    # TODO: update when vbot gets renamed (if ever)
    if vkey in ('z', 'depth', 'diameter', 'radius', 'height', 'x_velocity', 'y_velocity', 'port_spacing',
                'port_elevation', 'acute_mixing_zone', 'chronic_mixing_zone', 'depth', 'bottom_depth', 'x_displacement',
                'y_displacement'):
        return Length
    if vkey in ('speed', 'current_speed'):
        return Speed
    if vkey in ('current_dir', 'horizontal_angle', 'vertical_angle'):
        return Angle
    if vkey == 'salinity':
        return Salinity
    if vkey == 'temperature':
        return Temperature
    if vkey == 'density':
        return Density
    if vkey in ('concentration', 'bg_conc', 'v4o3'):
        return Concentration
    if vkey == 'decay_rate':
        return DecayRate
    if vkey == 'ff_velocity':
        return Speed
    if vkey == 'ff_dir':
        return Angle
    if vkey == 'ff_diff_coeff':
        # technically not unitless, but only one unit option available (m0.67/s^2)
        return Unitless
    if vkey == 'density':
        return Density
    if vkey in ('start_time', 'end_time', 'time_increment', 'total_time'):
        return Time
    if vkey == 'chronic_mixing_zone':
        return Length
    if vkey == 'effluent_flow':
        return FlowRate
    if vkey == 'isopleth':
        return Isopleth
    if vkey in ('eddy_diffusivity', 'diffusivity'):
        return EddyDiffusivity
    return Unitless


def convert(value, units_or_var_name, from_units, to_units=1, model_params=None, celsius=None, psu=None, depth=None):
    if not isinstance(units_or_var_name, str):
        if not issubclass(units_or_var_name, Units):
            raise Exception("Invalid unit type provided (not subclass of Units)")
        unit_handler = units_or_var_name
    else:
        unit_handler = from_var_name(units_or_var_name)
    if unit_handler is Salinity:
        assert model_params is not None
        assert celsius is not None
        return unit_handler.convert(
            value,
            from_units,
            # TODO: it's important that temperature has already been converted to celsius
            celsius=celsius,
            at_equilibrium=model_params.at_equilibrium,
            depth=depth
        )
    elif unit_handler is DecayRate:
        assert model_params is not None
        assert celsius is not None
        assert psu is not None
        return unit_handler.convert(
            value,
            from_units,
            bacteria_model=model_params.bacteria_model,
            celsius=celsius,
            psu=psu,
            depth=depth
        )
    return unit_handler.convert(value, from_units, to_units)
