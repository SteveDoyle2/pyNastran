from __future__ import annotations
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.atmosphere import (
    atm_density, atm_speed_of_sound,
    get_alt_for_pressure)
from pyNastran.utils.convert import (
    convert_altitude, convert_density, 
    convert_velocity, _velocity_factor,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNfloat


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

    velocity = tass  # sosi * machs
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
                                                                        np.ndarray]:
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
                                                                       np.ndarray]:
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
                                         eas_min: float = 0.0,
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
                                        eas_min: float=0.0,
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
                                         minus_eas: Optional[list[float]]=None,
                                         alt_units: str='ft',
                                         velocity_units: str='ft/s',
                                         density_units: str='slug/ft^3',
                                         eas_units: str='knots',
                                         ) -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
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
    gamma : float; default=1.4
        the gas constant
    minus_eas : float; default=0.0
        tag the velocity with a -1

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
    if minus_eas is None:
        minus_eas = []
    assert isinstance(mach, float), type(mach)
    nvel = len(eass)
    assert nvel > 0, eass

    eas = np.asarray(eass)  # knots or other
    ieas_mins = []
    for minus_easi in minus_eas:
        deas = eas - minus_easi
        ieas_min = np.argmin(np.abs(deas))
        ieas_mins.append(ieas_min)

    machs = np.ones(nvel, dtype=eas.dtype) * mach

    # get eas in ft/s and density in slug/ft^3,
    # so pressure is in psf
    #
    # then pressure in a sane unit (e.g., psi/psf/Pa)
    # without a wacky conversion
    assert eas.min() > 0, eas
    eas_fts = convert_velocity(eas, eas_units, 'ft/s')
    rho0_english: float = atm_density(
        0., R=1716., alt_units=alt_units,
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
        assert np.isfinite(alt_ft[i]), alt_fti
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
    for ieas_min in ieas_mins:
        eas[ieas_min] *= -1
        eas_fts[ieas_min] *= -1
        velocity[ieas_min] *= -1
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
