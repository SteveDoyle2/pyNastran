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
    else:
        raise RuntimeError(f'alt_units_in={alt_units_in!r} is not valid; use [m, ft, kft]')

    # ft to m
    if alt_units_out == 'm':
        factor *= 0.3048
    elif alt_units_out == 'ft':
        pass
    elif alt_units_out == 'kft':
        factor /= 1000.
    else:
        raise RuntimeError(f'alt_units_out={alt_units_out!r} is not valid; use [1/m, 1/in, 1/ft]')
    return factor

def convert_density(density: float, density_units_in: str, density_units_out: str) -> float:
    """nominal unit is slug/ft^3"""
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
    else:
        msg = f'density_units_in={density_units_in!r} is not valid; use [kg/m^3, slinch/in^3, slug/ft^3]'
        raise RuntimeError(msg)

    # data is now in slug/ft^3
    if density_units_out == 'slug/ft^3':
        pass
    elif density_units_out == 'slinch/in^3':
        factor /= 12**4
    elif density_units_out == 'kg/m^3':
        factor *= 515.378818
    else:
        msg = f'density_units_out={density_units_out!r} is not valid; use [kg/m^3, slinch/in^3, slug/ft^3]'
        raise RuntimeError(msg)
    return factor

def convert_pressure(pressure: float, pressure_units_in: str, pressure_units_out: str) -> float:
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
    elif pressure_units_in == 'Pa':
        factor /= 47.880172
    elif pressure_units_in == 'kPa':
        factor *= 20.88543815038
    elif pressure_units_in == 'MPa':
        factor *= 20885.43815038
    else:
        msg = f'pressure_units_in={pressure_units_in!r} is not valid; use [Pa, kPa, MPa, psf, psi]'
        raise RuntimeError(msg)

    if pressure_units_out == 'psf':
        pass
    elif pressure_units_out == 'psi':
        factor /= 144
    elif pressure_units_out == 'Pa':
        factor *= 47.880172
    elif pressure_units_out == 'kPa':
        factor /= 20.88543815038
    elif pressure_units_out == 'MPa':
        factor /= 20885.43815038
    else:
        raise RuntimeError(f'pressure_units_out={pressure_units_out} is not valid; use [Pa, kPa, MPa, psf, psi]')
    return factor

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

def convert_velocity(velocity: float, velocity_units_in: str, velocity_units_out: str) -> float:
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
    elif velocity_units_in == 'ft/s':
        pass
    elif velocity_units_in == 'in/s':
        factor /= 12.
    elif velocity_units_in == 'knots':
        factor *= 1.68781
    else:
        msg = f'velocity_units_in={velocity_units_in!r} is not valid; use [m/s, cm/s, in/s, ft/s, knots]'
        raise RuntimeError(msg)

    if velocity_units_out == 'm/s':
        factor *= 0.3048
    elif velocity_units_out == 'cm/s':
        factor *= 30.48
    elif velocity_units_out == 'ft/s':
        pass
    elif velocity_units_out == 'in/s':
        factor *= 12.
    elif velocity_units_out == 'knots':
        factor /= 1.68781
    else:
        msg = f'velocity_units_out={velocity_units_out!r} is not valid; use [m/s, cm/s, in/s, ft/s, knots]'
        raise RuntimeError(msg)
    return factor

def convert_length(velocity: float, length_units_in: str, length_units_out: str) -> float:
    """nominal unit is ft"""
    if length_units_in == length_units_out:
        return velocity
    return velocity * _length_factor(length_units_in, length_units_out)

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
    if mass_units_in == 'm':
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
    return velocity * _force_factor(force_units_in, force_units_out)

def _force_factor(force_units_in: str, force_units_out: str) -> float:
    """helper method for convert_force"""
    factor = 1.0
    if force_units_in == 'MN':
        factor *= 2248094.
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
        factor /= 2248094.
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

