def convert_altitude(alt: float, alt_units_in: str, alt_units_out: str) -> float:
    """nominal unit is ft"""
    if alt_units_in == alt_units_out:
        return alt
    return alt * _altitude_factor(alt_units_in, alt_units_out)


def _altitude_factor(alt_units_in: str, alt_units_out: str) -> float:
    """helper method for convert_altitude"""
    factor = 1.0
    # units to feet
    if alt_units_in == 'm':
        factor /= 0.3048
    elif alt_units_in == 'ft':
        pass
    elif alt_units_in == 'kft':
        factor *= 1000.
    else:  # pragma: no cover
        raise RuntimeError(f'alt_units_in={alt_units_in!r} is not valid; use [m, ft, kft]')

    # ft to m
    if alt_units_out == 'm':
        factor *= 0.3048
    elif alt_units_out == 'ft':
        pass
    elif alt_units_out == 'kft':
        factor /= 1000.
    else:  # pragma: no cover
        raise RuntimeError(f'alt_units_out={alt_units_out!r} is not valid; use [1/m, 1/in, 1/ft]')
    return factor


def convert_density(density: float, density_units_in: str,
                    density_units_out: str) -> float:
    """
    nominal unit is slug/ft^3
    TODO: change this to SI
    """
    if density_units_in == density_units_out:
        return density
    return density * _density_factor(density_units_in, density_units_out)


def _density_factor(density_units_in: str, density_units_out: str) -> float:
    """helper method for convert_density"""
    factor = 1.0
    if density_units_in == 'slug/ft^3':
        pass
    elif density_units_in == 'slinch/in^3':
        factor *= 12**4
    elif density_units_in == 'kg/m^3':
        factor /= 515.378818
    elif density_units_in == 'Mg/mm^3':
        factor /= 5.15378818e-10
    elif density_units_in == 'g/cm^3':
        factor *= (1000. / 515.378818)
    else:
        msg = f'density_units_in={density_units_in!r} is not valid; use ['\
              'kg/m^3, g/cm^3, Mg/mm^3, slinch/in^3, slug/ft^3]'
        raise RuntimeError(msg)

    # data is now in slug/ft^3
    if density_units_out == 'slug/ft^3':
        pass
    elif density_units_out == 'slinch/in^3':
        factor /= 12**4
    elif density_units_out == 'kg/m^3':
        factor *= 515.378818
    elif density_units_out == 'Mg/mm^3':
        factor *= 5.15378818e-10
    elif density_units_out == 'g/cm^3':
        factor /= (1000. / 515.378818)
    else:
        msg = f'density_units_out={density_units_out!r} is not valid; use ['\
              'kg/m^3, Mg/mm^3, g/cm^3, Mg/mm^3, slinch/in^3, slug/ft^3]'
        raise RuntimeError(msg)
    return factor


def convert_temperature(temperature: float, temperature_units_in: str,
                        temperature_units_out: str) -> float:
    """nominal unit is C"""
    if temperature_units_in == temperature_units_out:
        return temperature
    return _temperature_factor(temperature, temperature_units_in, temperature_units_out)


def _temperature_factor(temperature: float, temperature_units_in: str,
                        temperature_units_out: str) -> float:
    """helper method for convert_temperature"""
    temperature2 = temperature
    # 9/5 = 1.8
    if temperature_units_in == 'C':
        pass
    elif temperature_units_in == 'F':
        temperature2 = (temperature2 - 32) / 1.8
    elif temperature_units_in == 'K':
        temperature2 -= 273.15
    elif temperature_units_in == 'R':
        temperature2 = (temperature2 - 491.67) / 1.8
    else:
        msg = f'temperature_units_in={temperature_units_in!r} is not valid; use [C, R, F, K]'
        raise RuntimeError(msg)

    # data is now in C
    if temperature_units_out == 'C':
        pass
    elif temperature_units_out == 'F':
        temperature2 = temperature2 * 1.8 + 32
    elif temperature_units_out == 'K':
        temperature2 += 273.15
    elif temperature_units_out == 'R':
        temperature2 = temperature2 * 1.8 + 491.67
    else:
        msg = f'temperature_units_out={temperature_units_out!r} is not valid; use [C, R, F, K]'
        raise RuntimeError(msg)

    return temperature2


def convert_pressure(pressure: float, pressure_units_in: str,
                     pressure_units_out: str) -> float:
    """nominal unit is psf"""
    if pressure_units_in == pressure_units_out:
        return pressure
    return pressure * _pressure_factor(pressure_units_in, pressure_units_out)


def _pressure_factor(pressure_units_in: str, pressure_units_out: str) -> float:
    """helper method for convert_pressure"""
    factor = 1.0
    if pressure_units_in == 'psf':
        pass
    elif pressure_units_in == 'psi':
        factor *= 144
    elif pressure_units_in == 'bar':
        factor *= 2088.545633
    elif pressure_units_in == 'Pa':
        factor /= 47.880172
    elif pressure_units_in == 'kPa':
        factor *= 20.88543815038
    elif pressure_units_in == 'MPa':
        factor *= 20885.43815038
    else:
        msg = f'pressure_units_in={pressure_units_in!r} is not valid; use [Pa, kPa, MPa, bar, psf, psi]'
        raise RuntimeError(msg)

    if pressure_units_out == 'psf':
        pass
    elif pressure_units_out == 'psi':
        factor /= 144
    elif pressure_units_out == 'bar':
        factor /= 2088.545633
    elif pressure_units_out == 'Pa':
        factor *= 47.880172
    elif pressure_units_out == 'kPa':
        factor /= 20.88543815038
    elif pressure_units_out == 'MPa':
        factor /= 20885.43815038
    else:
        raise RuntimeError(f'pressure_units_out={pressure_units_out} is not valid; use [Pa, kPa, MPa, bar, psf, psi]')
    return factor


def convert_velocity(velocity: float, velocity_units_in: str,
                     velocity_units_out: str) -> float:
    """nominal unit is ft/s"""
    if velocity_units_in == velocity_units_out:
        return velocity
    return velocity * _velocity_factor(velocity_units_in, velocity_units_out)


def _velocity_factor(velocity_units_in: str, velocity_units_out: str) -> float:
    """helper method for convert_velocity"""
    factor = 1.0
    if velocity_units_in == 'm/s':
        factor /= 0.3048
    elif velocity_units_in == 'cm/s':
        factor /= 30.48
    elif velocity_units_in == 'mm/s':
        factor /= 304.8
    elif velocity_units_in == 'ft/s':
        pass
    elif velocity_units_in == 'in/s':
        factor /= 12.
    elif velocity_units_in == 'knots':
        factor *= 1.68781
    else:  # pragma: no cover
        msg = f'velocity_units_in={velocity_units_in!r} is not valid; use [m/s, cm/s, mm/s, in/s, ft/s, knots]'
        raise RuntimeError(msg)

    if velocity_units_out == 'm/s':
        factor *= 0.3048
    elif velocity_units_out == 'cm/s':
        factor *= 30.48
    elif velocity_units_out == 'mm/s':
        factor *= 304.8
    elif velocity_units_out == 'ft/s':
        pass
    elif velocity_units_out == 'in/s':
        factor *= 12.
    elif velocity_units_out == 'knots':
        factor /= 1.68781
    else:  # pragma: no cover
        msg = f'velocity_units_out={velocity_units_out!r} is not valid; use [m/s, cm/s, mm/s, in/s, ft/s, knots]'
        raise RuntimeError(msg)
    return factor


def convert_area(area: float, area_units_in: str, area_units_out: str) -> float:
    """nominal unit is ft^2"""
    if area_units_in == area_units_out:
        return area
    return area * _area_factor(area_units_in, area_units_out)


def _area_factor(area_units_in: str, area_units_out: str) -> float:
    """helper method for convert_area"""
    assert '^' in area_units_in, area_units_in
    assert '^' in area_units_out, area_units_out
    length_units_in, two_in = area_units_in.split('^')
    length_units_out, two_out = area_units_out.split('^')
    assert two_in == '2', two_in
    assert two_out == '2', two_out
    length_factor = _length_factor(length_units_in, length_units_out)
    area_factor = length_factor ** 2
    return area_factor


def convert_length(length: float, length_units_in: str,
                   length_units_out: str) -> float:
    """nominal unit is ft"""
    if length_units_in == length_units_out:
        return length
    return length * _length_factor(length_units_in, length_units_out)


def _length_factor(length_units_in: str, length_units_out: str) -> float:
    """helper method for convert_length"""
    factor = 1.0
    if length_units_in == 'm':
        factor /= 0.3048
    elif length_units_in == 'cm':
        factor /= 30.48
    elif length_units_in == 'mm':
        factor /= 304.8
    elif length_units_in == 'ft':
        pass
    elif length_units_in == 'in':
        factor /= 12.
    else:
        msg = f'length_units_in={length_units_in!r} is not valid; use [m, cm, mm, in, ft]'
        raise RuntimeError(msg)

    if length_units_out == 'm':
        factor *= 0.3048
    elif length_units_out == 'cm':
        factor *= 30.48
    elif length_units_out == 'mm':
        factor *= 304.8
    elif length_units_out == 'ft':
        pass
    elif length_units_out == 'in':
        factor *= 12.
    else:
        msg = f'length_units_out={length_units_out!r} is not valid; use [m, cm, mm, in, ft]'
        raise RuntimeError(msg)
    return factor


def convert_mass(mass: float, mass_units_in: str, mass_units_out: str) -> float:
    """nominal unit is slug"""
    if mass_units_in == mass_units_out:
        return mass
    return mass * _mass_factor(mass_units_in, mass_units_out)


def _mass_factor(mass_units_in: str, mass_units_out: str) -> float:
    """helper method for convert_mass"""
    factor = 1.0
    if mass_units_in == 'kg':
        factor *= 0.0685218
    elif mass_units_in == 'slug':
        pass
    elif mass_units_in == 'slinch':
        factor *= 12.
    else:
        msg = f'mass_units_in={mass_units_in!r} is not valid; use [kg, slug, slinch]'
        raise RuntimeError(msg)

    if mass_units_out == 'kg':
        factor /= 0.0685218
    elif mass_units_out == 'slug':
        pass
    elif mass_units_out == 'slinch':
        factor /= 12.
    else:
        msg = f'mass_units_out={mass_units_out!r} is not valid; use [kg, slug, slinch]'
        raise RuntimeError(msg)
    return factor


def convert_force(force: float, force_units_in: str, force_units_out: str) -> float:
    """nominal unit is lbf"""
    if force_units_in == force_units_out:
        return force
    return force * _force_factor(force_units_in, force_units_out)


def _force_factor(force_units_in: str, force_units_out: str) -> float:
    """helper method for convert_force"""
    factor = 1.0
    if force_units_in == 'MN':
        factor *= 224809.4
    elif force_units_in == 'N':
        factor *= 0.2248094
    elif force_units_in == 'cN':
        factor *= 0.0022480894
    elif force_units_in == 'mN':
        factor *= 0.0002248094
    elif force_units_in == 'lbf':
        pass
    else:
        msg = f'force_units_in={force_units_in!r} is not valid; use [MN, N, cN, mN, lbf]'
        raise RuntimeError(msg)

    if force_units_out == 'MN':
        factor /= 224809.4
    elif force_units_out == 'N':
        factor /= 0.2248094
    elif force_units_out == 'cN':
        factor /= 0.0022480894
    elif force_units_out == 'mN':
        factor /= 0.0002248094
    elif force_units_out == 'lbf':
        pass
    else:
        msg = f'force_units_out={force_units_out!r} is not valid; use [MN, N, cN, mN, lbf]'
        raise RuntimeError(msg)
    return factor
