"""
Contains the following atmospheric functions:

 - density = atm_density(alt, mach)
 - mach = atm_mach(alt, velocity)
 - velocity = atm_velocity(alt, mach)
 - pressure = atm_pressure(alt)
 - temperature = atm_temperature(alt)
 - sos = atm_speed_of_sound(alt)
 - mu = atm_dynamic_viscosity_mu(alt)
 - nu = atm_kinematic_viscosity_nu(alt)
 - eas = atm_equivalent_airspeed(alt, mach)
 - rho, machs, velocity = make_flfacts_alt_sweep_constant_mach(
       mach, alts, eas_limit=1000.,
       alt_units='m', velocity_units='m/s', density_units='kg/m^3',
       eas_units='m/s')
 - rho, machs, velocity = make_flfacts_mach_sweep(
       alt, machs, eas_limit=1000.,
       alt_units='m', velocity_units='m/s', density_units='kg/m^3',
       eas_units='m/s')

All the default units are in English units because the source equations
are in English units.

"""
from __future__ import annotations
import sys
from math import log, exp
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.utils.convert import (
    convert_altitude, convert_density, convert_pressure, convert_velocity,
    _altitude_factor, _pressure_factor, _velocity_factor)
#from pyNastran.utils import deprecated
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNfloat


def get_alt_for_density(density: float, density_units: str='slug/ft^3',
                        alt_units: str='ft', nmax: int=20, tol: float=5.) -> float:
    """
    Gets the altitude associated with a given air density.

    Parameters
    ----------
    density : float
        the air density in slug/ft^3
    density_units : str; default='slug/ft^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    alt_units : str; default='ft'
        sets the units for the output altitude; ft, m, kft
    nmax : int; default=20
        max number of iterations for convergence
    tol : float; default=5.
        tolerance in alt_units

    Returns
    -------
    alt : float
        the altitude in feet

    """
    tol = convert_altitude(tol, alt_units, 'ft')
    dalt = 500. # ft
    alt_old = 0.
    alt_final = 5000.
    n = 0

    #density_scale = _density_factor(density_units, "slug/ft^3")

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < nmax:
        alt_old = alt_final
        alt1 = alt_old
        alt2 = alt_old + dalt
        rho1 = atm_density(alt1, density_units=density_units)
        rho2 = atm_density(alt2, density_units=density_units)
        m = dalt / (rho2 - rho1)
        alt_final = m * (density - rho1) + alt1
        n += 1
    if abs(alt_final - alt_old) > tol:
        raise RuntimeError(f'Did not converge; Check your units; n=nmax={nmax}\n'
                           f'target alt={alt_final} alt_current={alt1}')
    alt_out = convert_altitude(alt_final, 'ft', alt_units)
    return alt_out


def get_alt_for_eas_with_constant_mach(equivalent_airspeed: float, mach: float,
                                       velocity_units: str='ft/s', alt_units: str='ft',
                                       nmax: int=20, tol: float=5.) -> float:
    """
    Gets the altitude associated with a equivalent airspeed.

    Parameters
    ----------
    equivalent_airspeed : float
        the equivalent airspeed in velocity_units
    mach : float
        the mach to hold constant
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    nmax : int; default=20
        max number of iterations for convergence
    tol : float; default=5.
        tolerance in alt_units

    Returns
    -------
    alt : float
        the altitude in alt units

    """
    equivalent_airspeed = convert_velocity(equivalent_airspeed, velocity_units, 'ft/s')
    tol = convert_altitude(tol, alt_units, 'ft')
    dalt = 500.
    alt_old = 0.
    alt_final = 5000.
    n = 0

    R = 1716.
    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)
    k = np.sqrt(T0 / p0)
    #eas = a * mach * sqrt((p * T0) / (T * p0)) = a * mach * sqrt(p / T) * k

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < nmax:
        alt_old = alt_final
        alt1 = alt_old
        alt2 = alt_old + dalt
        T1 = atm_temperature(alt1)
        T2 = atm_temperature(alt2)
        press1 = atm_pressure(alt1)
        press2 = atm_pressure(alt2)
        sos1 = np.sqrt(1.4 * R * T1)
        sos2 = np.sqrt(1.4 * R * T2)
        eas1 = sos1 * mach * np.sqrt(press1 / T1) * k
        eas2 = sos2 * mach * np.sqrt(press2 / T2) * k
        m = dalt / (eas2 - eas1)
        alt_final = m * (equivalent_airspeed - eas1) + alt1
        n += 1

    if n > nmax - 1:
        print(f'n = {n}')
    alt_final = convert_altitude(alt_final, 'ft', alt_units)
    return alt_final

def get_alt_for_q_with_constant_mach(q: float, mach: float,
                                     pressure_units: str='psf', alt_units: str='ft',
                                     nmax: int=20, tol: float=5.) -> float:
    """
    Gets the altitude associated with a dynamic pressure.

    Parameters
    ----------
    q : float
        the dynamic pressure lb/ft^2 (SI=Pa)
    mach : float
        the mach to hold constant
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    nmax : int; default=20
        max number of iterations for convergence
    tol : float; default=5.
        tolerance in alt_units

    Returns
    -------
    alt : float
        the altitude in alt_units

    """
    pressure = 2 * q / (1.4 * mach ** 2) # gamma = 1.4
    alt = get_alt_for_pressure(
        pressure, pressure_units=pressure_units, alt_units=alt_units, nmax=nmax, tol=tol)
    return alt

def get_alt_for_pressure(pressure: float,
                         pressure_units: str='psf', alt_units: str='ft',
                         nmax: int=20, tol: float=5.) -> float:
    """
    Gets the altitude associated with a pressure.

    Parameters
    ----------
    pressure : float
        the pressure lb/ft^2 (SI=Pa)
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    nmax : int; default=20
        max number of iterations for convergence
    tol : float; default=5.
        altitude tolerance in alt_units

    Returns
    -------
    alt : float
        the altitude in alt_units

    """
    pressure = convert_pressure(pressure, pressure_units, 'psf')
    tol = convert_altitude(tol, alt_units, 'ft')
    dalt = 500.
    alt_old = 0.
    alt_final = 5000.
    n = 0

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < nmax:
        alt_old = alt_final
        alt1 = alt_old
        alt2 = alt_old + dalt
        press1 = atm_pressure(alt1)
        press2 = atm_pressure(alt2)
        m = dalt / (press2 - press1)
        alt_final = m * (pressure - press1) + alt1
        n += 1

    if n > nmax - 1:
        print(f'n = {n}')
    #if abs(alt_final - alt_old) > tol:
        #raise RuntimeError('Did not converge; Check your units; n=nmax=%s\n'
                           #'target alt=%s alt_current=%s' % (nmax, alt_final, alt1))

    alt_final = convert_altitude(alt_final, 'ft', alt_units)
    return alt_final

def _feet_to_alt_units(alt_units: str) -> float:
    """helper method"""
    if alt_units == 'm':
        factor = 0.3048
    elif alt_units == 'ft':
        factor = 1.
    else:
        raise RuntimeError(f'alt_units={alt_units!r} is not valid; use [m, ft, kft]')
    return factor

def atm_temperature(alt: float,
                    alt_units: str='ft',
                    temperature_units: str='R') -> float:
    r"""
    Freestream Temperature \f$ T_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    temperature_units : str; default='R'
        the altitude units; R, K

    Returns
    -------
    T : float
        temperature in degrees Rankine or Kelvin (SI)

    .. note ::
        from BAC-7006-3352-001-V1.pdf  # Bell Handbook of Aerodynamic Heating\n
        page ~236 - Table C.1\n
        These equations were used because they are valid to 300k ft.\n
        Extrapolation is performed above that.

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    if z < 36151.725:
        T = 518.0 - 0.003559996 * z
    elif z < 82344.678:
        T = 389.988
    elif z < 155347.756:
        T = 389.988 + .0016273286 * (z - 82344.678)
    elif z < 175346.171:
        T = 508.788
    elif z < 249000.304:
        T = 508.788 - .0020968273 * (z - 175346.171)
    elif z < 299515.564:
        T = 354.348
    else:
        #print("alt=%i kft > 299.5 kft" % (z / 1000.))
        T = 354.348
        #raise AtmosphereError("altitude is too high")

    if temperature_units == 'R':
        factor = 1.
    elif temperature_units == 'K':
        factor = 5. / 9.
    else:
        raise RuntimeError(f'temperature_units={temperature_units!r} is not valid; use [R, K]')

    T2 = T * factor
    return T2

def atm_pressure(alt: float,
                 alt_units: str='ft',
                 pressure_units: str='psf') -> float:
    r"""
    Freestream Pressure \f$ p_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa

    Returns
    -------
    pressure : float
        Returns pressure in pressure_units

    .. note ::
        from BAC-7006-3352-001-V1.pdf  # Bell Handbook of Aerodynamic Heating\n
        page ~236 - Table C.1\n
        These equations were used b/c they are valid to 300k ft.\n
        Extrapolation is performed above that.\n

    """
    z = convert_altitude(alt, alt_units, 'ft')
    if z < 36151.725:
        ln_pressure = 7.657389 + 5.2561258 * log(1 - 6.8634634E-6 * z)
    elif z < 82344.678:
        ln_pressure = 6.158411 - 4.77916918E-5 * (z - 36151.725)
    elif z < 155347.756:
        ln_pressure = 3.950775 - 11.3882724 * log(1.0 + 4.17276598E-6 * (z - 82344.678))
    elif z < 175346.171:
        ln_pressure = 0.922461 - 3.62635373E-5*(z - 155347.756)
    elif z < 249000.304:
        ln_pressure = 0.197235 + 8.7602095 * log(1.0 - 4.12122002E-6 * (z - 175346.171))
    elif z < 299515.564:
        ln_pressure = -2.971785 - 5.1533546650E-5 * (z - 249000.304)
    else:
        #print("alt=%i kft > 299.5 kft" % (z / 1000.))
        ln_pressure = -2.971785 - 5.1533546650E-5 * (z - 249000.304)

    p = exp(ln_pressure)

    factor = _pressure_factor('psf', pressure_units)
    return p * factor

def atm_dynamic_pressure(alt: float, mach: float,
                         alt_units: str='ft',
                         pressure_units: str='psf') -> float:
    r"""
    Freestream Dynamic Pressure  \f$ q_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    pressure_units : str; default='psf'
        the pressure units; psf, psi, Pa, kPa, MPa

    Returns
    -------
    dynamic_pressure : float
        Returns dynamic pressure in pressure_units

    The common method that requires many calculations...
    \f[  \large q = \frac{1}{2} \rho V^2  \f]
    \f[  \large p = \rho R T  \f]
    \f[  \large M = \frac{V}{a}  \f]
    \f[  \large a = \sqrt{\gamma R T}  \f]
    so...
    \f[  \large q = \frac{\gamma}{2} p M^2  \f]

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    p = atm_pressure(z)
    q = 0.7 * p * mach ** 2

    factor = _pressure_factor('psf', pressure_units)
    q2 = q * factor
    return q2

def atm_speed_of_sound(alt: float,
                       alt_units: str='ft',
                       velocity_units: str='ft/s',
                       gamma: float=1.4) -> float:
    r"""
    Freestream Speed of Sound  \f$ a_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    speed_of_sound, a : float
        Returns speed of sound in velocity_units

    \f[  \large a = \sqrt{\gamma R T}  \f]

    """
    # converts everything to English units first
    z = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature(z)
    R = 1716. # 1716.59, dir air, R=287.04 J/kg*K

    a = (gamma * R * T) ** 0.5
    factor = _velocity_factor('ft/s', velocity_units) # ft/s to m/s
    a2 = a * factor
    return a2

def atm_velocity(alt: float, mach: float,
                 alt_units: str='ft',
                 velocity_units: str='ft/s') -> float:
    r"""
    Freestream Velocity  \f$ V_{\infty} \f$

    Parameters
    ----------
    alt : float
        altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    velocity : float
        Returns velocity in velocity_units

    \f[ \large V = M a \f]

    """
    a = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
    V = mach * a # units=ft/s or m/s
    return V

def atm_equivalent_airspeed(alt: float,
                            mach: float,
                            alt_units: str='ft',
                            eas_units: str='ft/s') -> float:
    """
    Freestream equivalent airspeed

    Parameters
    ----------
    alt : float
        altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    eas_units : str; default='ft/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    Returns
    -------
    eas : float
        equivalent airspeed in eas_units

    EAS = TAS * sqrt(rho/rho0)
    p = rho * R * T
    rho = p/(RT)
    rho/rho0 = p/T * T0/p0
    TAS = a * M
    EAS = a * M * sqrt(p/T * T0/p0)
    EAS = a * M * sqrt(p*T0 / (T*p0))

    """
    z = convert_altitude(alt, alt_units, 'ft')
    a = atm_speed_of_sound(z)
    #V = mach * a # units=ft/s or m/s

    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)

    T = atm_temperature(z)
    p = atm_pressure(z)

    eas = a * mach * np.sqrt((p * T0) / (T * p0))
    eas2 = convert_velocity(eas, 'ft/s', eas_units)
    return eas2

def atm_mach(alt: float,
             V: float,
             alt_units: str='ft',
             velocity_units: str='ft/s') -> float:
    r"""
    Freestream Mach Number

    Parameters
    ----------
    alt : float
        altitude in alt_units
    V : float
        Velocity in velocity_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    velocity_units : str; default='ft/s'
        the velocity units; ft/s, m/s, in/s, knots

    Returns
    -------
    mach : float
        Mach Number \f$ M \f$

    \f[ \large M = \frac{V}{a} \f]

    """
    a = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
    mach = V / a
    return mach

def atm_density(alt: float,
                R: float=1716.,
                alt_units: str='ft',
                density_units: str='slug/ft^3') -> float:
    r"""
    Freestream Density   \f$ \rho_{\infty} \f$

    Parameters
    ----------
    alt : float
        altitude in feet or meters
    R : float; default=1716.
        gas constant for air in english units (???)
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    density_units : str; default='slug/ft^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3

    Returns
    -------
    rho : float
        density \f$ \rho \f$ in density_units

    Based on the formula P=pRT
    \f[ \large \rho=\frac{p}{R T} \f]

    """
    z = convert_altitude(alt, alt_units, 'ft')
    #z = alt * _altitude_factor(alt_units, 'ft')
    P = atm_pressure(z)
    T = atm_temperature(z)

    rho = P / (R * T)
    rho2 = convert_density(rho, 'slug/ft^3', density_units)
    return rho2

def atm_kinematic_viscosity_nu(alt: float,
                               alt_units: str='ft',
                               visc_units: str='ft^2/s') -> float:
    r"""
    Freestream Kinematic Viscosity \f$ \nu_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    visc_units : str; default='slug/ft^3'
        the kinematic viscosity units; ft^2/s, m^2/s

    Returns
    -------
    nu : float
        kinematic viscosity \f$ \nu_{\infty} \f$ in visc_units

    \f[ \large \nu = \frac{\mu}{\rho} \f]

    .. seealso::  sutherland_viscoscity
    .. todo:: better debug

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density(z)
    mu = atm_dynamic_viscosity_mu(z)
    nu = mu / rho

    if visc_units == 'ft^2/s':
        factor = 1.
    elif visc_units == 'm^2/s':
        factor = _feet_to_alt_units(alt_units) ** 2
    else:  # pragma: no cover
        raise NotImplementedError(f'visc_units={visc_units!r}; use [m^2/s, ft^2/s]')
    return nu * factor

def atm_dynamic_viscosity_mu(alt: float,
                             alt_units: str='ft',
                             visc_units: str='(lbf*s)/ft^2') -> float:
    r"""
    Freestream Dynamic Viscosity  \f$ \mu_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    visc_units : str; default='(lbf*s)/ft^2'
        the viscosity units; (lbf*s)/ft^2, (N*s)/m^2, Pa*s

    Returns
    -------
    mu : float
        dynamic viscosity  \f$ \mu_{\infty} \f$ in (lbf*s)/ft^2 or (N*s)/m^2 (SI)

    .. seealso::  sutherland_viscoscity

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    T = atm_temperature(z)
    mu = sutherland_viscoscity(T)  # (lbf*s)/ft^2

    # same units as pressure, except multiplied by seconds
    if visc_units == '(lbf*s)/ft^2':
        factor = 1.
    elif visc_units in ['(N*s)/m^2', 'Pa*s']:
        factor = 47.88026
    else:  # pragma: no cover
        raise NotImplementedError(f'visc_units={visc_units!r}; use [(N*s)/m^2, Pa*s, (lbf*s)/ft^2]')
    return mu * factor

def atm_unit_reynolds_number2(alt: float, mach: float,
                              alt_units: str='ft',
                              reynolds_units: str='1/ft') -> float:
    r"""
    Returns the Reynolds Number per unit length.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    reynolds_units : str; default='1/ft'
        the altitude units; 1/ft, 1/m, 1/in

    Returns
    -------
    ReynoldsNumber/L : float
        the Reynolds Number per unit length

    \f[ \large Re_L = \frac{ \rho V}{\mu} = \frac{p M a}{\mu R T} \f]

    .. note ::
        this version of Reynolds number directly calculates the base quantities, so multiple
        calls to atm_press and atm_temp are not made

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    gamma = 1.4
    R = 1716.
    p = atm_pressure(z)
    T = atm_temperature(z)
    #p = rhoRT
    a = (gamma * R * T) ** 0.5
    mu = sutherland_viscoscity(T)
    ReL = p * a * mach / (mu * R * T)

    ReL *= _reynolds_factor('1/ft', reynolds_units)
    return ReL

def atm_unit_reynolds_number(alt: float, mach: float,
                             alt_units: str='ft',
                             reynolds_units: str='1/ft') -> float:
    r"""
    Returns the Reynolds Number per unit length.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    reynolds_units : str; default='1/ft'
        the altitude units; 1/ft, 1/m, 1/in

    Returns
    -------
    ReynoldsNumber/L : float
        Reynolds number per unit length in reynolds_units

    \f[ \large Re   = \frac{ \rho V L}{\mu} \f]
    \f[ \large Re_L = \frac{ \rho V  }{\mu} \f]

    """
    z = alt * _altitude_factor(alt_units, 'ft')
    rho = atm_density(z)
    V = atm_velocity(z, mach)
    mu = atm_dynamic_viscosity_mu(z)

    ReL = (rho * V) / mu

    ReL *= _reynolds_factor('1/ft', reynolds_units)
    return ReL


def _reynolds_factor(reynolds_units_in: str, reynolds_units_out: str) -> float:
    """helper method"""
    factor = 1.0
    # units to 1/feet
    if reynolds_units_in == '1/m':
        factor *= 0.3048
    elif reynolds_units_in == '1/ft':
        pass
    elif reynolds_units_in == '1/in':
        factor *= 12.
    else:
        msg = f'reynolds_units_in={reynolds_units_in!r} is not valid; use [1/ft, 1/m, 1/in]'
        raise RuntimeError(msg)

    # 1/ft to 1/m
    if reynolds_units_out == '1/m':
        factor /= 0.3048
    elif reynolds_units_out == '1/ft':
        pass
    elif reynolds_units_out == '1/in':
        factor /= 12.
    else:
        msg = f'reynolds_units_out={reynolds_units_out!r} is not valid; use [1/m, 1/in, 1/ft]'
        raise RuntimeError(msg)
    return factor

def sutherland_viscoscity(T: float) -> float:
    r"""
    Helper function that calculates the dynamic viscosity \f$ \mu \f$ of air at
    a given temperature.

    Parameters
    ----------
    T : float
        Temperature T is in Rankine

    Returns
    -------
    mu : float
       dynamic viscosity  \f$ \mu \f$ of air in (lbf*s)/ft^2

    .. note ::
        prints a warning if T>5400 deg R

    Sutherland's Equation\n
    From Aerodynamics for Engineers 4th Edition\n
    John J. Bertin 2002\n
    page 6 eq 1.5b\n

    """
    if T < 225.: # Rankine
        viscosity = 8.0382436E-10 * T
    else:
        if T > 5400.:
            msg = f'WARNING:  viscosity - Temperature is too large (T>5400 R) T={T}\n'
            sys.stderr.write(msg)
        viscosity = 2.27E-8 * (T ** 1.5) / (T + 198.6)
    return viscosity

#def make_flfacts_alt_sweep(mach: float, alts: np.ndarray,
                           #eas_limit: float=1000.,
                           #alt_units: str='m',
                           #velocity_units: str='m/s',
                           #density_units: str='kg/m^3',
                           #eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    #deprecated('make_flfacts_alt_sweep', 'make_flfacts_alt_sweep_constant_mach', '1.4',
               #levels=[0, 1, 2])
    #out = make_flfacts_alt_sweep_constant_mach(
        #mach, alts,
        #eas_limit=eas_limit, alt_units=alt_units,
        #velocity_units=velocity_units,
        #density_units=density_units,
        #eas_units=eas_units)
    #return out

def make_flfacts_alt_sweep_constant_mach(mach: float, alts: np.ndarray,
                                         eas_limit: float=1000.,
                                         alt_units: str='m',
                                         velocity_units: str='m/s',
                                         density_units: str='kg/m^3',
                                         eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Makes a sweep across altitude for a constant Mach number.

    Parameters
    ----------
    alts : list[float]
        Altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    """
    rho = np.array([atm_density(alt, R=1716., alt_units=alt_units,
                                density_units=density_units)
                    for alt in alts])
    machs = np.ones(len(alts)) * mach
    sos = np.array([atm_speed_of_sound(alt, alt_units=alt_units,
                                       velocity_units=velocity_units)
                    for alt in alts])
    velocity = sos * machs
    rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      alt_units=alt_units,
                                      density_units=density_units,
                                      velocity_units=velocity_units,
                                      eas_units=eas_units,)
    return rho, machs, velocity


def make_flfacts_tas_sweep_constant_alt(alt: float, tass: np.ndarray,
                                        eas_limit: float=1000.,
                                        alt_units: str='m',
                                        velocity_units: str='m/s',
                                        density_units: str='kg/m^3',
                                        eas_units: str='m/s') -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """TODO: not validated"""
    assert tass[0] <= tass[-1], tass

    rhoi = atm_density(alt, R=1716., alt_units=alt_units,
                      density_units=density_units)
    nvel = len(tass)
    rho = np.ones(nvel, dtype=tass.dtype) * rhoi

    sosi = atm_speed_of_sound(alt, alt_units=alt_units,
                              velocity_units=velocity_units)
    machs = tass / sosi

    velocity = tass # sosi * machs
    rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      alt_units=alt_units,
                                      density_units=density_units,
                                      velocity_units=velocity_units,
                                      eas_units=eas_units)
    return rho, machs, velocity

def make_flfacts_mach_sweep_constant_alt(alt: float, machs: list[float],
                                         eas_limit: float=1000.,
                                         alt_units: str='m',
                                         velocity_units: str='m/s',
                                         density_units: str='kg/m^3',
                                         eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Makes a sweep across Mach number for a constant altitude.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    machs : list[float]
        Mach Number \f$ M \f$
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    """
    assert machs[0] <= machs[-1], machs

    machs = np.asarray(machs)
    rho = np.ones(len(machs)) * atm_density(alt, R=1716., alt_units=alt_units,
                                            density_units=density_units)
    sos = np.ones(len(machs)) * atm_speed_of_sound(alt, alt_units=alt_units,
                                                   velocity_units=velocity_units)
    velocity = sos * machs
    rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      alt_units=alt_units,
                                      density_units=density_units,
                                      velocity_units=velocity_units,
                                      eas_units=eas_units,)
    return rho, machs, velocity

def make_flfacts_alt_sweep_constant_mach(mach: float, alts: list[float],
                                         eas_limit: float=1000.,
                                         alt_units: str='m',
                                         velocity_units: str='m/s',
                                         density_units: str='kg/m^3',
                                         eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Makes a sweep across Mach number for a constant altitude.

    Parameters
    ----------
    mach : float
        Mach Number \f$ M \f$
    alts : list[float]
        Altitude in alt_units
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    """
    assert alts[0] >= alts[-1], alts

    alts = np.asarray(alts)
    rho = [atm_density(alt, R=1716., alt_units=alt_units, density_units=density_units)
           for alt in alts]
    sos = [atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
           for alt in alts]

    rho = np.array(rho, dtype=alts.dtype)
    sos = np.array(sos, dtype=alts.dtype)

    velocity = sos * mach
    machs = np.ones(len(alts)) * mach
    rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      alt_units=alt_units,
                                      density_units=density_units,
                                      velocity_units=velocity_units,
                                      eas_units=eas_units,)
    return rho, machs, velocity

def make_flfacts_eas_sweep_constant_alt(alt: float, eass: list[float],
                                        alt_units: str='m',
                                        velocity_units: str='m/s',
                                        density_units: str='kg/m^3',
                                        eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Makes a sweep across equivalent airspeed for a constant altitude.

    Parameters
    ----------
    alt : float
        Altitude in alt_units
    eass : list[float]
        Equivalent airspeed in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    """
    assert eass[0] <= eass[-1], eass

    # convert eas to output units
    eass = np.atleast_1d(eass) * _velocity_factor(eas_units, velocity_units)
    rho = atm_density(alt, R=1716., alt_units=alt_units,
                      density_units=density_units)
    sos = atm_speed_of_sound(alt, alt_units=alt_units,
                             velocity_units=velocity_units)
    rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)
    velocity = eass * np.sqrt(rho0 / rho)
    machs = velocity / sos

    nvelocity = len(velocity)
    rhos = np.ones(nvelocity, dtype=velocity.dtype) * rho
    assert len(rhos) == len(machs)
    assert len(rhos) == len(velocity)
    return rhos, machs, velocity

#def make_flfacts_eas_sweep_constant_mach(mach: float, eass: list[float],
                                         #alt_units: str='m',
                                         #velocity_units: str='m/s',
                                         #density_units: str='kg/m^3',
                                         #eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    #"""
    #Makes a sweep across equivalent airspeed for a constant altitude.

    #Parameters
    #----------
    #mach : float
        #Constant mach number
    #eass : list[float]
        #Equivalent airspeed in eas_units
    #alt_units : str; default='m'
        #the altitude units; ft, kft, m
    #velocity_units : str; default='m/s'
        #the velocity units; ft/s, m/s, in/s, knots
    #density_units : str; default='kg/m^3'
        #the density units; slug/ft^3, slinch/in^3, kg/m^3
    #eas_units : str; default='m/s'
        #the equivalent airspeed units; ft/s, m/s, in/s, knots

    #"""
    #assert eass[0] <= eass[-1], eass

    ## convert eas to output units
    #eass = np.atleast_1d(eass) * _velocity_factor(eas_units, velocity_units)
    #rho = atm_density(alt, R=1716., alt_units=alt_units,
                      #density_units=density_units)
    #sos = atm_speed_of_sound(alt, alt_units=alt_units,
                             #velocity_units=velocity_units)
    #rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)
    #velocity = eass * np.sqrt(rho0 / rho)
    #machs = velocity / sos

    #nvelocity = len(velocity)
    #rhos = np.ones(nvelocity, dtype=velocity.dtype) * rho
    #assert len(rhos) == len(machs)
    #assert len(rhos) == len(velocity)
    #return rhos, machs, velocity

#def make_flfacts_eas_sweep_constant_mach(
            #gamma: float=1.4,
            #alt_units='m',
            #velocity_units='m/s',
            #density_units='kg/m^3',
            #eas_units='m/s',
            #pressure_units='Pa'):
    #"""
    #eas = tas * sqrt(rho/rho0)
    #ainf*Minf = V
    #eas = ainf*Minf * sqrt(rho_inf/rho0)
    #rho = p/RT
    #eas = ainf*Minf * sqrt(p_inf/(R*T_inf*rho0))
        #= sqrt(gamma*R*Tinf) * Minf * sqrt(p_inf/(R*T_inf*rho0))
        #= Minf * sqrt(gamma*p_inf/rho0)
    #"""
    #return rho, vel, mach


def make_flfacts_eas_sweep_constant_mach(machs: np.ndarray,
                                         eass: np.ndarray,
                                         gamma: float=1.4,
                                         alt_units: str='ft',
                                         velocity_units: str='ft/s',
                                         density_units: str='slug/ft^3',
                                         #pressure_units: str='psf',
                                         eas_units: str='knots') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Veas = Vtrue * sqrt(rho/rho0)
    V = a * mach
    a = sqrt(gamma*R*T)
    p = rho*R*T -> rho = p/(R*T)
    V = Veas * sqrt(rho0 / rho)

    Veas = mach * sqrt(gamma*R*T) * sqrt(p/(R*T*rho0))
    Veas = mach * sqrt(gamma*p / rho0)
    Veas^2 / mach^2 = gamma * p / rho0
    p = Veas^2 / mach^2 * rho0/gamma
    """
    assert len(machs) == len(eass)
    assert len(machs) > 0, machs
    mach = np.asarray(machs)
    eas = np.asarray(eass)
    nmach = len(mach)

    # get eas in ft/s and density in slug/ft^3,
    # so pressure is in psf
    #
    # then pressure in a sane unit (e.g., psi/psf/Pa)
    # without a wacky conversion
    eas_fts = convert_velocity(eas, eas_units, 'ft/s')
    rho0_english = atm_density(0., R=1716., alt_units=alt_units,
                               density_units='slug/ft^3')
    pressure_psf = (eas_fts / machs) ** 2 * rho0_english / gamma  # psf

    rho = np.zeros(nmach, machs.dtype)
    alt_ft = np.zeros(nmach, machs.dtype)
    for i, pressure_psfi in enumerate(pressure_psf):
        alt_fti = get_alt_for_pressure(pressure_psfi, pressure_units='psf',
                                       alt_units='ft', nmax=30, tol=1.)
        rhoi = atm_density(alt_fti, R=1716., alt_units='ft', density_units=density_units)
        rho[i] = rhoi
        alt_ft[i] = alt_fti
    alt = convert_altitude(alt_ft, 'ft', alt_units)

    # eas = Vtas * sqrt(rho/rho0)
    # Vtas = eas * sqrt(rho0/rho)
    rho0 = convert_density(rho0_english, 'slug/ft^3', density_units)
    velocity_fts = eas_fts * np.sqrt(rho0 / rho)
    velocity = convert_velocity(velocity_fts, 'ft/s', velocity_units)

    #rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      #alt_units=alt_units,
                                      #density_units=density_units,
                                      #velocity_units=velocity_units,
                                      #eas_units=eas_units,)
    assert len(rho) == len(machs)
    assert len(rho) == len(velocity)
    return rho, mach, velocity, alt


def _limit_eas(rho: NDArrayNfloat, machs: NDArrayNfloat, velocity: NDArrayNfloat,
               eas_limit: float=1000.,
               alt_units: str='m',
               velocity_units: str='m/s',
               density_units: str='kg/m^3',
               eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """limits the equivalent airspeed"""
    assert alt_units != '', alt_units
    assert velocity_units != '', velocity_units
    assert density_units != '', density_units
    assert eas_units != '', eas_units

    if eas_limit:
        rho0 = atm_density(0., alt_units=alt_units, density_units=density_units)

        # eas in velocity units
        eas = velocity * np.sqrt(rho / rho0)
        kvel = _velocity_factor(eas_units, velocity_units)
        eas_limit_in_velocity_units = eas_limit * kvel

        i = np.where(eas < eas_limit_in_velocity_units)
        rho = rho[i]
        machs = machs[i]
        velocity = velocity[i]

    if len(rho) == 0:
        #print('machs min: %.0f max: %.0f' % (machs.min(), machs.max()))
        #print('vel min: %.0f max: %.0f in/s' % (velocity.min(), velocity.max()))
        #print('EAS min: %.0f max: %.0f in/s' % (eas.min(), eas.max()))
        raise RuntimeError('EAS limit is too struct and has removed all the conditions.\n'
                           'Increase eas_limit or change the mach/altude range\n'
                           '  EAS: min=%.3f max=%.3f limit=%s %s' % (
                               eas.min() / kvel,
                               eas.max() / kvel,
                               eas_limit, eas_units))
    return rho, machs, velocity
