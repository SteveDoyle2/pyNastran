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
"""
from __future__ import print_function
import sys
from math import log, exp
import numpy as np

def _update_alt(alt, alt_units='ft', debug=False):
    """
    Converts altitude alt_units to feet

    Paramters
    ---------
    alt : float
        altitude in feet or meters
    alt_units : str; default='ft'
        sets the units for altitude
        TODO: remove default

    Returns
    -------
    alt2 : float
        altitude in feet
    """
    if alt_units == 'ft':
        factor = 1.
    elif alt_units == 'm':
        factor = 1. / 0.3048
    elif alt_units == 'kft':
        factor = 1000.
    else:
        raise RuntimeError('alt_units=%r is not valid; use [ft, m, kft]' % alt_units)
    alt2 = alt * factor

    if debug:
        if alt_units == 'ft':
            SI = False
        else:
            SI = True
        if SI:
            print("z = %s [m] = %s [ft]"  % (alt, alt2))
        else:
            print("z = %s [m] = %s [ft]" % (alt * _feet_to_alt_units(alt_units), alt2))
    return alt2

def _update_velocity(velocity, velocity_units='ft/s', debug=False):
    """
    Converts altitude alt_units to feet

    Paramters
    ---------
    velocity : float
        altitude in feet or meters
    velocity_units : str; default='ft/s'
        sets the units for altitude
        TODO: remove default

    Returns
    -------
    velocity2 : float
        velocity in feet/s
    """
    if velocity_units == 'ft/s':
        factor = 1.
    elif velocity_units == 'm/s':
        factor = 1. / 0.3048
    elif velocity_units == 'knots':
        factor = 1.68781
    else:
        msg = 'velocity_units=%r is not valid; use [ft/s, m/s, knots]' % velocity_units
        raise RuntimeError(msg)
    velocity2 = velocity * factor

    #if debug:
        #if SI:
            #print("z = %s [m] = %s [ft]"  % (alt, alt2))
        #else:
            #print("z = %s [m] = %s [ft]" % (alt * _feet_to_meters(True), alt2))
    return velocity2


def get_alt_for_density(density):
    """
    Gets the altitude associated with a given air density.

    Parameters
    ----------
    density : float
        the air density in slug/ft^3

    Returns
    -------
    alt : float
        the altitude in feet
    """
    dalt = 500.
    alt_old = 0.
    alt_final = 5000.
    n = 0
    tol = 5. # ft

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < 20:
        alt_old = alt_final
        alt1 = alt_old
        alt2 = alt_old + dalt
        rho1 = atm_density(alt1)
        rho2 = atm_density(alt2)
        m = dalt / (rho2 - rho1)
        alt_final = m * (density - rho1) + alt1
        n += 1
    if n > 18:
        print('n = %s' % n)
    return alt_final


def get_alt_for_eas_mach(equivalent_airspeed, mach, SI=False):
    """
    Gets the altitude associated with a equivalent airspeed.

    Parameters
    ----------
    equivalent_airspeed : float
        the equivalent airspeed in ft/s (SI=m/s)
    mach : float
        the mach to hold constant
    SI : bool
        should SI units be used; default=False

    Returns
    -------
    alt : float
        the altitude in ft (SI=m)
    """
    if SI:
        equivalent_airspeed /= 0.3048  # m/s to ft/s
    dalt = 500.
    alt_old = 0.
    alt_final = 5000.
    n = 0
    tol = 5. # ft

    R = 1716.
    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)
    k = np.sqrt(T0 / p0)
    #eas = a * mach * sqrt((p * T0) / (T * p0)) = a * mach * sqrt(p / T) * k

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < 20:
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

    if n > 18:
        print('n = %s' % n)
    if SI:
        alt_final *= 0.3048  # feet to meters
    return alt_final

def get_alt_for_q_mach(q, mach, SI=False):
    """
    Gets the altitude associated with a equivalent airspeed.

    Parameters
    ----------
    q : float
        the dynamic pressure lb/ft^2 (SI=Pa)
    mach : float
        the mach to hold constant
    SI : bool
        should SI units be used; default=False

    Returns
    -------
    alt : float
        the altitude in ft (SI=m)
    """
    pressure = 2 * q / (1.4 * mach ** 2) # gamma = 1.4
    alt = get_alt_for_pressure(pressure, SI=SI)
    return alt

def get_alt_for_pressure(pressure, SI=False):
    """
    Gets the altitude associated with a equivalent airspeed.

    Parameters
    ----------
    pressure : float
        the pressure lb/ft^2 (SI=Pa)
    SI : bool
        should SI units be used; default=False

    Returns
    -------
    alt : float
        the altitude in ft (SI=m)
    """
    if SI:
        pressure /= 47.880259  # Pa to psf
    dalt = 500.
    alt_old = 0.
    alt_final = 5000.
    n = 0
    tol = 5. # ft

    # Newton's method
    while abs(alt_final - alt_old) > tol and n < 20:
        alt_old = alt_final
        alt1 = alt_old
        alt2 = alt_old + dalt
        press1 = atm_pressure(alt1)
        press2 = atm_pressure(alt2)
        m = dalt / (press2 - press1)
        alt_final = m * (pressure - press1) + alt1
        n += 1

    if n > 18:
        print('n = %s' % n)
    if SI:
        alt_final *= 0.3048  # feet to meters
    return alt_final

def _feet_to_alt_units(alt_units):
    if alt_units == 'm':
        factor = 0.3048
    elif alt_units == 'ft':
        factor = 1.
    else:
        raise RuntimeError('alt_units=%r is not valid; use [ft, m]' % alt_units)
    return factor

def _convert_alt(alt, alt_units_in, alt_units_out):
    factor = 1.0
    # units to feet
    if alt_units_in == 'm':
        factor /= 0.3048
    elif alt_units_in == 'ft':
        pass
    elif alt_units_in == 'kft':
        factor *= 1000.
    else:
        raise RuntimeError('alt_units_in=%r is not valid; use [ft, m, kft]' % alt_units_in)
    #print('alt=%.1f alt_units_in=%s alt_mid=%.1f ft' % (
        #alt, alt_units_in, alt*factor))

    # ft to m
    if alt_units_out == 'm':
        factor *= 0.3048
    elif alt_units_out == 'ft':
        pass
    elif alt_units_out == 'kft':
        factor /= 1000.
    else:
        raise RuntimeError('alt_units_out=%r is not valid; use [ft, m, kft]' % alt_units_out)
    #print('alt=%.1f alt_units_in=%s alt_units_out=%s alt2=%.1f' % (
        #alt, alt_units_in, alt_units_out, alt*factor))
    return alt * factor

def _convert_velocity(velocity, velocity_units_in, velocity_units_out):
    factor = 1.0
    if velocity_units_in == 'm/s':
        factor /= 0.3048
    elif velocity_units_in == 'ft/s':
        pass
    elif velocity_units_in == 'knots':
        factor *= 1.68781
    else:
        msg = 'velocity_units_in=%r is not valid; use [ft/s, m/s, knots]' % velocity_units_in
        raise RuntimeError(msg)
    #print('velocity=%.1f velocity_units_in=%s velocity_mid=%.1f' % (
        #velocity, velocity_units_in, velocity * factor))

    if velocity_units_out == 'm/s':
        factor *= 0.3048
    elif velocity_units_out == 'ft/s':
        pass
    elif velocity_units_out == 'knots':
        factor /= 1.68781
    else:
        msg = 'velocity_units_in=%r is not valid; use [ft/s, m/s, knots]' % velocity_units_in
        raise RuntimeError(msg)
    #print('velocity=%.1f velocity_units_in=%s velocity_units_out=%s velocity2=%.1f' % (
        #velocity, velocity_units_in, velocity_units_out, velocity * factor))
    return velocity * factor

def _feet_s_to_velocity_units(velocity_units):
    if velocity_units == 'm/s':
        factor = 0.3048
    elif velocity_units == 'ft/s':
        factor = 1.
    elif velocity_units == 'knots':
        factor = 1. / 1.68781
    else:
        raise RuntimeError('alt_units=%r is not valid; use [ft, m]' % velocity_units)
    return factor

def _rankine_to_kelvin(SI):
    if SI:
        factor = 5 / 9.
    else:
        factor = 1.
    return factor

def _psf_to_pressure_units(pressure_units):
    if pressure_units == 'psf':
        factor = 1.
    elif pressure_units == 'Pa':
        factor = 47.880259
    else:
        raise RuntimeError('pressure_units=%r is not valid; use [psf, Pa]' % pressure_units)
    return factor

def atm_temperature(alt, alt_units='ft', temperature_units='R', debug=False):
    r"""
    Freestream Temperature \f$ T_{\infty} \f$

    Paramters
    ---------
    alt : bool
        Altitude in feet or meters (SI)
    SI : bool; default=False
        returns temperature in SI units if True

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
    z = _update_alt(alt, alt_units)
    if z < 36151.725:
        T = 518.0-0.003559996 * z
    elif z < 82344.678:
        T = 389.988
    elif z < 155347.756:
        T = 389.988+.0016273286 * (z - 82344.678)
    elif z < 175346.171:
        T = 508.788
    elif z < 249000.304:
        T = 508.788-.0020968273 * (z - 175346.171)
    elif z < 299515.564:
        T = 354.348
    else:
        print("alt=%i kft > 299.5 kft" % (z / 1000.))
        T = 354.348
        #raise AtmosphereError("altitude is too high")

    if temperature_units == 'R':
        factor = 1.
    elif temperature_units == 'K':
        factor = 5. / 9.
    else:
        raise RuntimeError('temperature_units=%r is not valid; use [ft, m]' % temperature_units)

    T2 = T * factor
    #if debug:
        #if SI:
            #print("z = %s [m] = %s [ft]"  % (alt, z))
            #print("T = %s [K] = %s [R]"  % (T2, T))
        #else:
            #print("z = %s [m] = %s [ft]" % (alt * _feet_to_meters(True), z))
            #print("T = %s [K] = %s [R]" % (T * _rankine_to_kelvin(True), T2))
    return T2

def atm_pressure(alt, alt_units='ft', pressure_units='psf', debug=False):
    r"""
    Freestream Pressure \f$ p_{\infty} \f$
    Paramters
    ---------
    alt : float
        Altitude in feet or meters (SI)
    SI : bool; default=False
        returns pressure in SI units if True

    Returns
    -------
    pressure : float
        Returns pressure in psf or Pa (SI)

    .. note ::
        from BAC-7006-3352-001-V1.pdf  # Bell Handbook of Aerodynamic Heating\n
        page ~236 - Table C.1\n
        These equations were used b/c they are valid to 300k ft.\n
        Extrapolation is performed above that.\n
    """
    z = _update_alt(alt, alt_units=alt_units)
    if z < 36151.725:
        lnP = 7.657389 + 5.2561258 * log(1 - 6.8634634E-6 * z)
    elif z < 82344.678:
        lnP = 6.158411 - 4.77916918E-5 * (z-36151.725)
    elif z < 155347.756:
        lnP = 3.950775 - 11.3882724 * log(1.0 + 4.17276598E-6 * (z - 82344.678))
    elif z < 175346.171:
        lnP = 0.922461 - 3.62635373E-5*(z - 155347.756)
    elif z < 249000.304:
        lnP = 0.197235 + 8.7602095 * log(1.0 - 4.12122002E-6 * (z - 175346.171))
    elif z < 299515.564:
        lnP = -2.971785 - 5.1533546650E-5 * (z - 249000.304)
    else:
        print("alt=%i kft > 299.5 kft" % (z / 1000.))
        lnP = -2.971785 - 5.1533546650E-5 * (z - 249000.304)

    p = exp(lnP)

    if pressure_units == 'psf':
        factor = 1.
    elif pressure_units == 'Pa':
        factor = 47.880259 # psf to Pa
    else:
        raise RuntimeError('pressure_units=%r is not valid; use [psf, Pa]' % pressure_units)

    #if debug:
        #ft_to_m = _feet_to_meters(True)
        #if SI:
            #print("z    = %s [m]  = %s [ft]" % (alt, z))
            #print("Patm = %g [Pa] = %g [psf]" % (p * factor, p))
        #else:
            #print("z    = %s [m]  = %s [ft]" % (alt * ft_to_m, z))
            #print("Patm = %g [Pa] = %g [psf]" % (p * _psf_to_pascals(True), p))
    return p*factor

def atm_dynamic_pressure(alt, mach, alt_units='ft', pressure_units='psf', debug=False):
    r"""
    Freestream Dynamic Pressure  \f$ q_{\infty} \f$

    Parameters
    ----------
    alt : float
        Altitude in feet or meters (SI)
    mach : float
        Mach Number \f$ M \f$
    SI : bool
        returns dynamicPressure in SI units if True (default=False)

    Returns
    -------
    dynamic_Pressure : float
        Returns dynamic pressure in lb/ft^2 or Pa (SI).

    The common method that requires many calculations...
    \f[  \large q = \frac{1}{2} \rho V^2  \f]
    \f[  \large p = \rho R T  \f]
    \f[  \large M = \frac{V}{a}  \f]
    \f[  \large a = \sqrt{\gamma R T}  \f]
    so...
    \f[  \large q = \frac{\gamma}{2} p M^2  \f]
    """
    z = _update_alt(alt, alt_units)
    p = atm_pressure(z)
    q = 0.7 * p * mach ** 2

    factor = _psf_to_pressure_units(pressure_units)
    q2 = q * factor

    #if debug:
        #ft_to_m = _feet_to_meters(True)
        #if SI:
            #print("z = %s [m]   = %s [ft]" % (alt, z))
            #print("p = %s [psf] = %s [Pa]" % (p, p * factor))
            #print("q = %s [psf] = %s [Pa]" % (q, q2))
        #else:
            #print("z = %s [m]   = %s [ft]" % (alt * ft_to_m, z))
            #print("p = %s [psf] = %s [Pa]" % (p, p * _psf_to_pascals(True)))
            #print("q = %s [psf] = %s [Pa]" % (q, q * _psf_to_pascals(True)))
    return q2

def atm_speed_of_sound(alt, alt_units='ft', velocity_units='ft/s', gamma=1.4, debug=False):
    r"""
    Freestream Speed of Sound  \f$ a_{\infty} \f$

    Parameters
    ----------
    alt : bool
        Altitude in feet or meters (SI)
    SI : bool; default=False
        convert to SI units

    Returns
    -------
    speed_of_sound, a : float
        Returns speed of sound in ft/s or m/s (SI).
   \f[  \large a = \sqrt{\gamma R T}  \f]
    """
    # converts everything to English units first
    z = _update_alt(alt, alt_units)
    T = atm_temperature(z)
    R = 1716. # 1716.59, dir air, R=287.04 J/kg*K

    a = (gamma * R * T) ** 0.5
    factor = _feet_s_to_velocity_units(velocity_units) # ft/s to m/s
    a2 = a * factor

    #if debug:
        #ft_to_m = _feet_to_meters(True)
        #if SI:
            #print("z = %s [m]   = %s [ft]" % (alt, z))
            #print("T = %s [K]   = %s [R]" % (T / 1.8, T))
            #print("a = %s [m/s] = %s [ft/s]" % (a2, a))
        #else:
            #print("z = %s [m]   = %s [ft]" % (alt * ft_to_m, z))
            #print("T = %s [K]   = %s [R]" % (T / 1.8, T))
            #print("a = %s [m/s] = %s [ft/s]" % (a * _feet_to_meters(True), a2))
    return a2

def atm_velocity(alt, mach, alt_units='ft', velocity_units='ft/s', debug=False):
    r"""
    Freestream Velocity  \f$ V_{\infty} \f$
    alt : float
        altitude in feet or meters
    SI : bool; default=False
        convert velocity to SI units
    Mach : float
        Mach Number \f$ M \f$

    Returns
    velocity : float
        Returns velocity in ft/s or m/s (SI).

    \f[ \large V = M a \f]
    """
    a = atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
    V = mach * a # units=ft/s or m/s

    #if debug:
        #ft_to_m = _feet_to_alt_units('m')
        #if SI:
            #print("z = %s [m]   = %s [ft]"  % (alt, alt))
            #print("a = %s [m/s] = %s [ft/s]"  % (a, a / ft_to_m))
            #print("V = %s [m/s] = %s [ft/s]"  % (V, V / ft_to_m))
        #else:
            #print("z = %s [m]   = %s [ft]" % (alt * ft_to_m, alt))
            #print("a = %s [m/s] = %s [ft/s]" % (a * ft_to_m, a))
            #print("V = %s [m/s] = %s [ft/s]" % (V * ft_to_m, V))
    return V

def atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='ft/s', debug=False):
    """
    EAS = TAS * sqrt(rho/rho0)
    p = rho * R * T
    rho = p/(RT)
    rho/rho0 = p/T * T0/p0
    TAS = a * M
    EAS = a * M * sqrt(p/T * T0/p0)
    """
    z = _update_alt(alt, alt_units)
    a = atm_speed_of_sound(z)
    #V = mach * a # units=ft/s or m/s

    z0 = 0.
    T0 = atm_temperature(z0)
    p0 = atm_pressure(z0)

    T = atm_temperature(z)
    p = atm_pressure(z)

    eas = a * mach * np.sqrt((p * T0) / (T * p0))
    if eas_units == 'ft/s':
        pass
    elif eas_units == 'knots':
        eas /= 1.68781 # ft/s to knots
    elif eas_units == 'm/s':
        ft_to_m = _feet_to_alt_units('m')
        eas *= ft_to_m
    else:
        raise NotImplementedError(eas_units)


    #if debug:
        #if SI:
            #print("z = %s [m]   = %s [ft]"  % (alt, z))
            #print("a = %s [m/s] = %s [ft/s]"  % (a * ft_to_m, a))
            #print("eas = %s [m/s] = %s [ft/s]"  % (eas, eas / ft_to_m))
        #else:
            #ft_to_m = _feet_to_meters(True)
            #print("z = %s [m]   = %s [ft]" % (alt * ft_to_m, alt))
            #print("a = %s [m/s] = %s [ft/s]" % (a * ft_to_m, a))
            #print("eas = %s [m/s] = %s [ft/s]" % (eas * ft_to_m, eas))
    return eas

def atm_mach(alt, V, alt_units='ft', velocity_units='ft/s', debug=False):
    r"""
    Freestream Mach Number

    Parameters
    ----------
    alt : float
        altitude in feet or meters
    V : float
        Velocity in ft/s or m/s (SI)
    SI : bool; default=False
        convert velocity to SI units

    Returns
    -------
    mach : float
        Mach Number \f$ M \f$

    \f[ \large M = \frac{V}{a} \f]
    """
    z = _update_alt(alt, alt_units)
    a = atm_speed_of_sound(z, alt_units='ft', velocity_units=velocity_units)
    mach = V / a

    if debug:
        print("z = %.1f [m] = %.1f [ft] = %.1f [%s]"  % (
            _convert_alt(alt, alt_units, 'm'),
            z, # ft
            alt, alt_units))
        print("a = %.3f [m/s] = %.3f [ft/s] = %.3f [%s]"  % (
            _convert_velocity(a, velocity_units, 'm/s'),
            _convert_velocity(a, velocity_units, 'ft/s'),
            a, velocity_units))
        print("V = %.3f [m/s] = %.3f [ft/s] = %.3f [%s]"  % (
            _convert_velocity(V, velocity_units, 'm/s'),
            _convert_velocity(V, velocity_units, 'ft/s'),
            V, velocity_units))
        print("M = %.3f"  % (mach))
    return mach

def atm_density(alt, R=1716., alt_units='ft', density_units='slug/ft^3', debug=False):
    r"""
    Freestream Density   \f$ \rho_{\infty} \f$

    Parameters
    ----------
    Parameters
    ----------
    alt : float
        altitude in feet or meters
    SI : bool; default=False
        convert velocity to SI units
    R : float; default=1716.
        gas constant for air in english units (???)

    Returns
    -------
    rho : float
        density \f$ \rho \f$ in slug/ft^3 or kg/m^3 (SI).

    Based on the formula P=pRT
    \f[ \large \rho=\frac{p}{R T} \f]
    """
    z = _update_alt(alt, alt_units)
    P = atm_pressure(z)
    T = atm_temperature(z)

    # going from slug/ft^3 to kg/m^3
    if density_units == 'slug/ft^3':
        factor = 1.
    elif density_units == 'kg/m^3':
        factor = 515.378818
    #elif density_units == 'slug/in^3':
        #factor = None
    else:
        raise NotImplementedError(density_units)

    #if debug:
        #rho = P / (R * T)
        #ft_to_m = _feet_to_alt_units('m')
        #if SI:
            #pressure_units = 'Pa'
            #print("z    = %s [m] = %s [ft]" % (alt, z))
            #print("Patm = %g [Pa] = %g [psf]" % (P * _psf_to_pressure_units(pressure_units), P))
            #print("T    = %s [K] = %s [R]" % (T / 1.8, T))
            #print("rho  = %e [kg/m^3] = %e [slug/ft^3]" % (rho * 515.378818, rho))
        #else:
            #pressure_units = 'Pa'
            #print("z    = %s [m] = %s [ft]" % (alt * ft_to_m, z))
            #print("Patm = %g [Pa] = %g [psf]" % (P * _psf_to_pressure_units(pressure_units), P))
            #print("T    = %s [K] = %s [R]" % (T / 1.8, T))
            #print("rho  = %e [kg/m^3] = %e [slug/ft^3]" % (rho * 515.378818, rho))

    return P / (R * T) * factor

def atm_kinematic_viscosity_nu(alt, alt_units='ft', visc_units='ft^2/s', debug=False):
    r"""
    Freestream Kinematic Viscosity \f$ \nu_{\infty} \f$

    Parameters
    ----------
    alt : bool
        Altitude in feet or meters (SI)
    SI : bool; default=False
        convert to SI units

    Returns
    -------
    nu : float
        kinematic viscosity \f$ \nu_{\infty} \f$ in ft^2/s or m^2/s (SI)

    \f[ \large \nu = \frac{\mu}{\rho} \f]

    .. see ::  SutherlandVisc
    .. todo:: better debug
    """
    z = _update_alt(alt, alt_units)
    rho = atm_density(z)
    mu = atm_dynamic_viscosity_mu(z)
    nu = mu / rho
    if debug:  # doesnt work unless US units
        print("atm_nu - rho=%g [slug/ft^3] mu=%e [lb*s/ft^2] nu=%e [ft^2/s]" % (rho, mu, nu))

    if visc_units == 'ft^2/s':
        factor = 1.
    elif visc_units == 'm^2/s':
        factor = _feet_to_alt_units(alt_units) ** 2
    return nu * factor

def atm_dynamic_viscosity_mu(alt, SI=False):
    r"""
    Freestream Dynamic Viscosity  \f$ \mu_{\infty} \f$

    Parameters
    ----------
    alt : bool
        Altitude in feet or meters (SI)
    SI : bool; default=False
        convert to SI units

    Returns
    -------
    mu : float
        dynamic viscosity  \f$ \mu_{\infty} \f$ in (lbf*s)/ft^2 or (N*s)/m^2 (SI)

    .. see ::  SutherlandVisc
    @ todo units...
    """
    z = _update_alt(alt, SI)
    T = atm_temperature(z)
    mu = sutherland_viscoscity(T)
    if SI:
        return mu * 47.88026
    return mu

def atm_unit_reynolds_number2(alt, mach, alt_units='ft', ReL_units='1/ft', debug=False):
    r"""
    Returns the Reynolds Number per unit length

    Parameters
    ----------
    alt : bool
        Altitude in feet or meters (SI)
    mach : float
        Mach Number \f$ M \f$
    SI : bool; default=False
        convert to SI units

    Returns
    -------
    ReynoldsNumber/L : float
        1/ft or 1/m (SI)

    \f[ \large Re_L = \frac{ \rho V}{\mu} = \frac{p M a}{\mu R T} \f]

    .. note ::
        this version of Reynolds number directly caculates the base quantities, so multiple
        calls to atm_press and atm_temp are not made
    """
    z = _update_alt(alt, alt_units)
    #print "z = ",z
    gamma = 1.4
    R = 1716.
    p = atm_pressure(z)
    T = atm_temperature(z)
    #p = rhoRT
    a = (gamma * R * T) ** 0.5
    mu = sutherland_viscoscity(T)
    ReL = p * a * mach / (mu * R * T)

    if debug:
        print("---atm_UnitReynoldsNumber2---")
        print("z  = %s [m]   = %s [ft]"  % (alt * _feet_to_alt_units('m'), z))
        print("a  = %s [m/s] = %s [ft/s]"  % (a * _feet_to_alt_units('m'), a))
        rho = p / (R * T)
        print("rho = %s [kg/m^3] = %s [slug/ft^3]"  % (rho * 515.378818, rho))
        print("M  = %s"  % mach)
        print("V  = %s [m/s] = %s [ft/s]"  % (a * mach * _feet_to_alt_units('m'), a * mach))
        print("T  = %s [K] = %s [R]" % (T * 5 / 9., T))
        print("mu = %s [(N*s)/m^2] = %s [(lbf*s)/ft^2]" % (mu * 47.88026, mu))
        print("Re = %s [1/m] = %s [1/ft]" % (ReL / 0.3048, ReL))

    # convert ReL in 1/ft to 1/m
    if ReL_units == '1/ft':
        factor = 1.
    elif ReL_units == '1/m':
        factor = 1. / .3048
    else:
        raise NotImplementedError(ReL_units)
    return ReL * factor

def atm_unit_reynolds_number(alt, mach, SI=False, debug=False):
    r"""
    Returns the Reynolds Number per unit length

    Parameters
    ----------
    alt : bool
        Altitude in feet or meters (SI)
    mach : float
        Mach Number \f$ M \f$
    SI : bool; default=False
        convert to SI units

    Returns
    -------
    ReynoldsNumber/L : float
        1/ft or 1/m (SI)

    \f[ \large Re   = \frac{ \rho V L}{\mu} \f]
    \f[ \large Re_L = \frac{ \rho V  }{\mu} \f]
    """
    z = _update_alt(alt, SI)
    rho = atm_density(z)
    V = atm_velocity(z, mach)
    mu = atm_dynamic_viscosity_mu(z)

    ReL = (rho * V) / mu

    if debug:
        print("---atm_UnitReynoldsNumber---")
        print("z  = %s [m]   = %s [ft]"  % (alt * _feet_to_alt_units('m'), z))
        print("rho = %s [kg/m^3] = %s [slug/ft^3]"  % (rho * 515.378818, rho))
        print("V  = %s [m/s] = %s [ft/s]"  % (V * _feet_to_alt_units('m'), V))
        print("mu = %s [(N*s)/m^2] = %s [(lbf*s)/ft^2]" % (mu * 47.88026, mu))
        print("Re = %s [1/m] = %s [1/ft]" % (ReL / 0.3048, ReL))

    if SI:
        return ReL / .3048  # convert ReL in 1/ft to 1/m
    return ReL

def sutherland_viscoscity(T):
    r"""
    Helper function that calculates the dynamic viscosity \f$ \mu \f$ of air at
    a given temperature

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
    .. todo:: Consider raising an error instead of writing to stderr
               and letting the function return an answer.

    Sutherland's Equation\n
    From Aerodynamics for Engineers 4th Edition\n
    John J. Bertin 2002\n
    page 6 eq 1.5b\n
    """
    if T < 225.: # Rankine
        viscosity = 8.0382436E-10 * T
    else:
        if T > 5400.:
            msg = "WARNING:  viscosity - Temperature is too large (T>5400) T=%s\n" % T
            sys.stderr.write(msg)
        viscosity = 2.27E-8 * (T ** 1.5) / (T + 198.6)
    return viscosity
