import math
from .Ambient import Ambient
from ..globals import missing, UserInputError
from ..params.ModelParameters import BacteriaModel, ModelParameters


# Used for seawater density function
Ar = (
    999.842594,   -9.095290e-3,  -1.120083e-6,   8.24493e-1,    7.6438e-5,     5.3875e-9,     1.0227e-4,    4.8314e-4,
    6.793952e-2,   1.001685e-4,   6.536332e-9,  -4.0899e-3,    -8.2467e-7,    -5.72466e-3,   -1.6546e-6
)
Kr = (
    19652.21,      148.4206,      1.360477e-2,   3.239908,      1.16092e-4,    8.50935e-5,    5.2787e-8,    54.6746,
    1.09987e-2,    7.944e-2,     -5.3009e-4,    -1.0981e-5,     1.91075e-4,    2.0816e-8,    -2.327105,    -5.155288e-5,
    1.43713e-3,   -5.77905e-7,   -6.12293e-6,   -0.603459,     -6.167e-5,      1.6483e-2,     2.2838e-3,   -1.6078e-6,
    -9.9348e-7,    9.1697e-10
)
# Used for seawater temperature function

# Used for conductivity function
mho = (
    (0,  9.341, 17.456, 25.238, 32.851),
    (0, 10.816, 20.166, 29.090, 37.778),
    (0, 12.361, 23.010, 33.137, 42.962),
    (0, 13.967, 25.967, 37.351, 48.367),
    (0, 15.628, 29.027, 41.713, 53.963),
    (0, 17.345, 32.188, 46.213, 59.732),
    (0, 19.127, 35.458, 50.856, 65.667)
)


# Used for seawater density function
def rho(salinity, temperature, pressure):
    # TODO: optimized by multiples of temperature powers, check no transcribing errors
    sal_p15 = salinity**(1.5)
    rhoST = (
        # x temperature^0
        Ar[0] + salinity*(Ar[3] + Ar[7]*salinity) + Ar[13]*sal_p15 + temperature*(
            # x temperature^1
            Ar[6]*sal_p15 + Ar[8] + Ar[11]*salinity + temperature*(
                # x temperature^2
                Ar[1] + Ar[4]*salinity + Ar[14]*sal_p15 + temperature*(
                    # x temperature^3
                    Ar[9] + Ar[12]*salinity + temperature*(
                        # x temperature^4        x temperature^5
                        Ar[2] * Ar[5]*salinity + temperature*Ar[10]
                    )
                )
            )
        )
    )
    prs_p2 = pressure**2
    KSTP = (
        # x temperature^0
        Kr[0]
        + Kr[3]*pressure
        + Kr[5]*prs_p2
        + Kr[7]*salinity
        + sal_p15*(Kr[9] + Kr[12]*pressure)
        + pressure*salinity*(Kr[22] + Kr[24]*pressure)
        # x temperature^1
        + temperature*(
            Kr[1]
            + pressure*(
                salinity*(Kr[11] + Kr[13]*pressure)
                + Kr[16]
                + Kr[18]*pressure
            )
            + Kr[19]*salinity
            + Kr[21]*sal_p15
            # x temperature^2
            + temperature*(
                pressure*(Kr[4] + Kr[6]*pressure)
                + Kr[8]*salinity
                + Kr[1]*sal_p15
                + Kr[14]
                + pressure*(Kr[23]*salinity + Kr[25]*pressure)
                # x temperature^3
                + temperature*(
                    Kr[2]
                    + Kr[17]*pressure
                    + Kr[20]*salinity
                    # x temperature^4
                    + temperature*Kr[15]
               )
            )
        )
    )
    return rhoST / (1 - pressure / KSTP)


#-----------------------------------------------------------------------------------
# FUNCTION REMAPS
#-----------------------------------------------------------------------------------
# sigmat(sal, t)                  -> seawater_density(salinity, temperature, at_equilibrium, ambient_cond, in_sigma)
# sigmasal(con, tm)               -> salinity(conductivity, temperature, at_equilibrium, params)
# ambientlevel(z)                 -> ambient_level(params, ambient_store, z)
# interpolateambient(z, temponly) -> interpolate_ambient(params, ambient_store, z, temp_only=False)
# mancini(arg, S, T, z, topersec) -> mancini(arg, salinity, temperature, z, topersec)
# mhocon(sal, tm)                 -> mho_conductivity(salinity, temperature)
# mhosal(con, tm)                 -> mho_salinity(conductivity, temperature)


#-----------------------------------------------------------------------------------
# ambient condition interpolation and value reading
#-----------------------------------------------------------------------------------
def seawater_density(salinity=None, temperature=None, at_equilibrium=False, ambient_cond=None, depth=None,
                     in_sigma=False):
    if temperature is None or salinity is None:
        assert ambient_cond is not None
        assert isinstance(ambient_cond, Ambient)
    if temperature is None:
        temperature = ambient_cond.temperature
    if salinity is None:
        salinity = ambient_cond.salinity
    if at_equilibrium:
        sig0  = (((6.8e-6*salinity) - 4.82e-4)*salinity + 0.8149)*salinity - 0.093
        b     = 1.0e-6*temperature * ((0.01667*temperature - 0.8164)*temperature + 18.03)
        a     = 0.001*temperature * ((0.0010843*temperature - 0.09818)*temperature + 4.7867)
        sumt  = (temperature - 3.98)**2 *(temperature + 283.0) / (503.57*(temperature + 67.26))
        sigma = (sig0 + 0.1324) * (1.0 - a + b * (sig0 - 0.1324)) - sumt
    else:
        if depth is None:
            if ambient_cond is None:
                raise Exception("Density calculation not at equilibirum requires depth parameter or Ambient")
            assert isinstance(ambient_cond, Ambient)
            depth = ambient_cond.depth
        sigma = rho(salinity, temperature, depth*0.1) - 1000
    return sigma if in_sigma else 1000+sigma


def salinity(temperature, density, at_equilibrium=False, depth=None):
    if density < -250:
        raise Exception("Density confusion? Salinity set to 0")
    s1 = 15.0
    s2 = s1 + 1e-8
    salinity = 0
    for i in range(21):
        if i >= 20:
            return missing
        sig1 = seawater_density(s1, temperature, at_equilibrium=at_equilibrium, depth=depth)
        sig2 = seawater_density(s2, temperature, at_equilibrium=at_equilibrium, depth=depth)
        ps = (sig2 - sig1) / 1e-8
        ds = 1e-7 if ps == 0 else (density - sig1)/ps
        s1 += ds
        s2 += ds
        salinity = s1 + ds
    if salinity < 0:
        salinity = 0
        # memo(f"Small negative salinity ({salinity:.6f}) set to zero")  # (oc) sigmasal bug
    return salinity


def mancini(model, arg, salinity, temperature, depth, topersec):
    assert isinstance(model, ModelParameters)
    lim = 1 if depth <= 0 else (1 - math.exp(-model.light_absorb_coeff*depth))/model.light_absorb_coeff/depth
    # these are random constant that were defined globally in main.pas
    sphr = 3600
    secpday = 86400
    match model.bacteria_model:
        case BacteriaModel.COLIFORM_MANCINI:
            const_a = 1.07**(temperature-20)
            if salinity < 0:
                dum = 0.8*const_a
            elif salinity < 33:
                dum = (0.8+0.006*100*salinity/33)*const_a
            else:
                dum = 1.4*const_a
            if topersec:
                return (dum + arg*lim) / secpday
            else:
                return (arg*secpday - dum) / lim
        case BacteriaModel.COLIFORM_301H:
            const_a = 2.303 * math.exp(2.303*(0.0295*temperature - 2.292))
            if topersec:
                return (2.303*const_a + 1.24*arg*lim*0.04184) / sphr
            else:
                return (sphr*arg - const_a)/1.24/0.04184/lim
        case BacteriaModel.ENTEROCCOUS_301H:
            if topersec:
                return (0.5262/24+1.24*arg*lim*0.04184)/sphr
            else:
                return (sphr*arg-0.5262/24)/1.24/0.04184/lim
        case _:
            raise UserInputError("Unknown or unrecognized bacteria model")


def mho_conductivity(salinity, temperature):
    t = int(temperature)
    s = int(salinity/10)  # 10=delsal
    if 0 <= t <= 5 and 0 <= s <= 3:
        a = (salinity - s*10)/10  # 10=delsal
        lo = a*(mho[t][s+1] - mho[t][s]) + mho[t][s]
        hi = a*(mho[t+1][s+1] - mho[t+1][s]) + mho[t+1][s]
        mhocon = (temperature - t*5)/5*(hi - lo) + lo  # 5=deltem
        return mhocon, True
    # not good value
    return 0, False


def mho_salinity(conductivity, temperature):
    t = int(temperature/5)  # 5=delstem
    if t < 0 or t > 5:
        # not good value
        return 0, False
    a = (temperature - t*5)/5
    for s in range(0, 4):
        lo = a*(mho[t+1][s] - mho[t][s]) + mho[t][s]
        hi = a*(mho[t+1][s+1] - mho[t][s+1]) + mho[t][s+1]
        salinity = (conductivity - lo)/(hi - lo)*10 + s*10  # 10=delsal
        if lo <= conductivity or hi > conductivity:
            return salinity, True
    # not good value
    return 0, False
