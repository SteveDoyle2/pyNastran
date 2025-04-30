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
 - rho, machs, velocity = make_flfacts_mach_sweep_constant_alt(
       alt, machs, eas_limit=1000.,
       alt_units='m', velocity_units='m/s', density_units='kg/m^3',
       eas_units='m/s')

All the default units are in English units because the source equations
are in English units.

"""
from __future__ import annotations
import sys
from math import log, exp
from itertools import count
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
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3
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
    alt1 = alt_old
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
    velocity_units : str; default=ft/s
        the velocity units; ft/s, m/s, in/s, knots
    nmax : int; default=20
        max number of iterations for convergence
    tol : float; default=5.
        tolerance in alt_units

    Returns
    -------
    alt : float
        the altitude in alt units

    """
    equivalent_airspeed = convert_velocity(
        equivalent_airspeed, velocity_units, 'ft/s')
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
    else:  # pragma: no cover
        raise RuntimeError(f'alt_units={alt_units!r} is not valid; use [m, ft, kft]')
    return factor

def get_alt_for_mach_eas(mach: float,
                         eas: float, eas_units: str='knots',
                         alt_units: str='ft',
                         nmax: int=20, tol: float=5.) -> float:
    """
    Gets the altitude associated with an equivalent airspped.

    Parameters
    ----------
    mach : float
        the Mach number
    eas : float
        the equivalent airspeed in eas_units
    eas_units : str; default='knots'
        the equivalent airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s
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
    eas = convert_velocity(eas, eas_units, 'ft/s')
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
        eas1 = atm_equivalent_airspeed(alt1, mach, alt_units='ft', eas_units='ft/s')
        eas2 = atm_equivalent_airspeed(alt2, mach, alt_units='ft', eas_units='ft/s')

        # y = (y2-y1)/(x2-x1)*(x-x1) + y1
        # y = m * (x-x1) + y1
        m = dalt / (eas2 - eas1)
        alt_final = m * (eas - eas1) + alt1
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
    else:  # pragma: no cover
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
    gamma: float; default=1.4
        the specific heat ratio Cp/Cv; gamma=1.4 for air

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
        the equivalent airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s

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

def atm_calibrated_airspeed(alt: float,
                            mach: float,
                            alt_units: str='ft',
                            cas_units: str='ft/s') -> float:
    """
    Calibrated airspeed

    Parameters
    ----------
    alt : float
        altitude in alt_units
    mach : float
        Mach Number \f$ M \f$
    alt_units : str; default='ft'
        the altitude units; ft, kft, m
    cas_units : str; default='ft/s'
        the calibrated airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s

    Returns
    -------
    cas : float
        Calibrated airspeed in cas_units

    CAS = EAS * sqrt(p/p0)
    https://aerotoolbox.com/airspeed-conversions/
    """
    eas = atm_equivalent_airspeed(alt, mach, alt_units=alt_units,
                                  eas_units=cas_units)

    z = convert_altitude(alt, alt_units, 'ft')
    p0 = atm_pressure(0.)  # psf
    p = atm_pressure(z)  # psf
    a0 = atm_speed_of_sound(z, velocity_units=cas_units)
    qc = p * ((1 + 0.2*mach**2)**3.5 - 1)
    mach_comp = np.sqrt(5 * (qc/p0+1)**(1/3.5) - 1)
    cas = a0 * mach_comp
    return cas


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
        the velocity units; ft/s, in/s, knots, m/s, cm/s, mm/s

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
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3

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
    mu = sutherland_viscoscity(T)

    #a = (gamma * R * T) ** 0.5
    #ReL = p * a * mach / (mu * R * T)
    ReL = p * mach / mu * (gamma / (R * T)) ** 0.5

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

def _make_flfacts_tas_sweep_constant_eas(eas: float, tass: np.ndarray,
                                         alt_units: str='m',
                                         velocity_units: str='m/s',
                                         density_units: str='kg/m^3',
                                         eas_units: str='m/s') -> tuple[np.ndarray, np.ndarray,
                                                                        np.ndarray]:  # pragma: no cover
    """
    Veas = Vtas*sqrt(rho/rho0)
    Veas/Vtas = sqrt(rho/rho0)
    (Veas/Vtas)^2 = rho/rho0
    rho = rho0 * (eas / tass) ** 2
    """
    rho0 = atm_density(0.0, R=1716., alt_units=alt_units,
                       density_units=density_units)
    ntas = len(tass)
    mach = np.full(ntas, np.nan, dtype=tass.dtype)
    rho = rho0 * (eas / tass) ** 2
    for i, (rhoi, tasi) in zip(count(), rho, tass):
        alt = get_alt_for_density(
            rhoi, density_units=density_units,
            alt_units=alt_units, nmax=20, tol=5.)
        sos = atm_speed_of_sound(alt, alt_units=alt_units,
                                 velocity_units=velocity_units)
        machi = tasi / sos
        mach[i] = machi
    velocity = tass
    return rho, mach, velocity

def _make_flfacts_alt_sweep_constant_eas(eas: float, alts: np.ndarray,
                                         alt_units: str='m',
                                         velocity_units: str='m/s',
                                         density_units: str='kg/m^3',
                                         eas_units: str='m/s') -> tuple[np.ndarray, np.ndarray,
                                                                        np.ndarray]:  # pragma: no cover
    """
    Veas = Vtas * sqrt(rho/rho0)
    Vtas = Veas * sqrt(rho0/rho)
    Vtas = sos * Mach
    Mach = Vtas / sos = Veas/sos * sqrt(rho0/rho)
    """
    eas_velocity_units = convert_velocity(eas, eas_units, velocity_units)
    rho, sos = _rho_sos_for_alts(
        alts, alt_units=alt_units,
        density_units=density_units,
        velocity_units=velocity_units)
    rho0 = atm_density(0.0, R=1716., alt_units=alt_units,
                       density_units=density_units)
    velocity = eas_velocity_units * np.sqrt(rho0/rho)
    mach = velocity / sos
    return rho, mach, velocity

def make_flfacts_alt_sweep_constant_tas(tas: float, alts: np.ndarray,
                                        alt_units: str='m',
                                        velocity_units: str='m/s',
                                        density_units: str='kg/m^3',
                                        eas_limit: float=1000.,
                                        eas_units: str='m/s') -> tuple[np.ndarray, np.ndarray,
                                                                       np.ndarray]:    # pragma: no cover
    """
    Veas = Vtas * sqrt(rho/rho0)
    Vtas = Veas * sqrt(rho0/rho)
    Vtas = sos * Mach
    Mach = Vtas / sos = Veas/sos * sqrt(rho0/rho)
    """
    rho, sos = _rho_sos_for_alts(
        alts, alt_units=alt_units,
        density_units=density_units,
        velocity_units=velocity_units)
    mach = tas / sos
    velocity = tas * np.ones(len(sos), dtype=sos.dtype)

    rho, machs, velocity = _limit_eas(rho, mach, velocity, eas_limit,
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
        the velocity units; ft/s, in/s, knots, m/s, cm/s, mm/s
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s

    """
    assert machs[0] <= machs[-1], machs

    machs = np.asarray(machs)
    one = np.ones(len(machs))
    rho = one * atm_density(alt, R=1716., alt_units=alt_units,
                            density_units=density_units)
    sos = one * atm_speed_of_sound(alt, alt_units=alt_units,
                                   velocity_units=velocity_units)
    velocity = sos * machs
    rho, machs, velocity = _limit_eas(rho, machs, velocity, eas_limit,
                                      alt_units=alt_units,
                                      density_units=density_units,
                                      velocity_units=velocity_units,
                                      eas_units=eas_units,)
    return rho, machs, velocity

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
    mach : float
        Mach Number \f$ M \f$
    alts : list[float]
        Altitude in alt_units
    eas_limit : float
        Equivalent airspeed limiter in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, in/s, knots, m/s, cm/s, mm/s
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s

    """
    rho, sos = _rho_sos_for_alts(
        alts, alt_units=alt_units,
        density_units=density_units,
        velocity_units=velocity_units)
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
        the velocity units; ft/s, in/s, knots, m/s, cm/s, mm/s
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, in/s, knots, m/s, cm/s, mm/s

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


def make_flfacts_eas_sweep_constant_mach(mach: float,
                                         eass: np.ndarray,
                                         gamma: float=1.4,
                                         alt_units: str='ft',
                                         velocity_units: str='ft/s',
                                         density_units: str='slug/ft^3',
                                         eas_units: str='knots') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """
    Makes a sweep across equivalent airspeed for a constant altitude.

    Parameters
    ----------
    mach : float
        Constant mach number
    eass : list[float]
        Equivalent airspeed in eas_units
    alt_units : str; default='m'
        the altitude units; ft, kft, m
    velocity_units : str; default='m/s'
        the velocity units; ft/s, m/s, in/s, knots
    density_units : str; default='kg/m^3'
        the density units; slug/ft^3, slinch/in^3, kg/m^3, g/cm^3, Mg/mm^3
    eas_units : str; default='m/s'
        the equivalent airspeed units; ft/s, m/s, in/s, knots

    Veas = Vtas * sqrt(rho/rho0)
    a * mach = Vtas
    a = sqrt(gamma*R*T)
    p = rho*R*T -> rho = p/(R*T)
    Vtas = Veas * sqrt(rho0 / rho)

    Veas = mach * sqrt(gamma*R*T) * sqrt(p/(R*T*rho0))
    Veas = mach * sqrt(gamma*p / rho0)
    Veas^2 / mach^2 = gamma * p / rho0
    p = Veas^2 / mach^2 * rho0/gamma
    """
    assert isinstance(mach, float), type(mach)
    nvel = len(eass)
    assert nvel > 0, eass

    eas = np.asarray(eass)
    machs = np.ones(nvel, dtype=eas.dtype) * mach

    # get eas in ft/s and density in slug/ft^3,
    # so pressure is in psf
    #
    # then pressure in a sane unit (e.g., psi/psf/Pa)
    # without a wacky conversion
    eas_fts = convert_velocity(eas, eas_units, 'ft/s')
    rho0_english = atm_density(0., R=1716., alt_units=alt_units,
                               density_units='slug/ft^3')
    pressure_psf = (eas_fts / mach) ** 2 * rho0_english / gamma  # psf

    # lookup altitude by pressure to get dentisy
    # could be faster if we reused the pressure instead of ignoring it
    rho = np.zeros(nvel, eas.dtype)
    alt_ft = np.zeros(nvel, eas.dtype)
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
    return rho, machs, velocity, alt


def _limit_eas(rho: NDArrayNfloat, machs: NDArrayNfloat, velocity: NDArrayNfloat,
               eas_limit: float=1000.,
               alt_units: str='m',
               velocity_units: str='m/s',
               density_units: str='kg/m^3',
               eas_units: str='m/s') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """limits the equivalent airspeed"""
    assert len(rho) > 0, rho
    assert len(machs) > 0, machs
    assert len(velocity) > 0, velocity
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

def _rho_sos_for_alts(alts: np.ndarray,
                          alt_units: str='m',
                          density_units: str='kg/m^3',
                          velocity_units: str='m/s') -> tuple[np.ndarray, np.ndarray]:
    """gets the density and speed of sound arrays for a set of altitudes"""
    assert alts[0] >= alts[-1], alts
    alts = np.asarray(alts)
    rho_list = [atm_density(alt, R=1716., alt_units=alt_units, density_units=density_units)
                for alt in alts]
    sos_list = [atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
                for alt in alts]
    rho = np.array(rho_list, dtype=alts.dtype)
    sos = np.array(sos_list, dtype=alts.dtype)
    return rho, sos

def create_atmosphere_table(quantities: list[str],
                            alt_min: float, alt_max: float,
                            nalt: int=0,
                            dalt: float=0.0,
                            ) -> np.array:
    """
    Parameters
    ----------
    alt_min: float
    alt_max: float
    nalt: int
    dalt: float

    out = create_atmosphere_table(quantities, alt_min=0, alt_max=10000, dalt=1000.)
    out[:, 0]  # alt
    [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]

    """
    allowed_quantities = [
        'density', 'pressure', 'temperature', 'velocity',
        'speed_of_sound', 'dynamic_viscosity',
    ]
    for quantity in quantities:
        assert quantity in allowed_quantities, f'{quantity} is not a valid quantity; use {allowed_quantities}'
    assert nalt == 0 or dalt == 0.0, f'nalt={nalt} dalt={dalt}'
    assert nalt > 0 or dalt > 0.0, f'nalt={nalt} dalt={dalt}'
    if nalt == 0:
        alts = np.arange(alt_min, alt_max+dalt, dalt)
    else:
        alts = np.linspace(alt_min, alt_max, num=nalt)

    nalti = len(alts)
    out_list = [alts]
    alt_units = 'ft'
    density_units = 'slug/ft^3'
    pressure_units = 'psf'
    temperature_units = 'R'
    velocity_units = 'ft/s'
    dynamic_viscosity_units = '(lbf*s)/ft^2'
    R = 1716.0  # english
    for quantity in quantities:
        out = np.full(nalti, np.nan, dtype='float64')
        if quantity == 'density':
            for ialt, alt in enumerate(alts):
                out[ialt] = atm_density(alt, R, alt_units=alt_units, density_units=density_units)
        elif quantity == 'pressure':
            for ialt, alt in enumerate(alts):
                out[ialt] = atm_pressure(alt, alt_units=alt_units, pressure_units=pressure_units)
        elif quantity == 'temperature':
            for ialt, alt in enumerate(alts):
                out[ialt] = atm_temperature(alt, alt_units=alt_units, temperature_units=temperature_units)
        elif quantity == 'speed_of_sound':
            for ialt, alt in enumerate(alts):
                out[ialt] = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
        elif quantity == 'dynamic_viscosity':
            for ialt, alt in enumerate(alts):
                val = atm_dynamic_viscosity_mu(alt, alt_units=alt_units, visc_units=dynamic_viscosity_units)
        else:
            raise NotImplementedError(quantity)
        out_list.append(out)

    out_array = np.column_stack(out_list)
    return out_array


def main():  # pragma: no cover
    import os
    dirname = os.path.dirname(__file__)
    quantities = ['density', 'pressure', 'speed_of_sound']
    x = create_atmosphere_table(
        quantities,
        0., 50000.,
        dalt=1000.)
    csv_filename = os.path.join(dirname, 'atmosphere.csv')
    quantity_list = ['alt'] + quantities
    header = ','.join(quantity_list)
    np.savetxt(csv_filename, x, delimiter=',', header=header)


if __name__ == '__main__':  # pragma: no cover
    main()
