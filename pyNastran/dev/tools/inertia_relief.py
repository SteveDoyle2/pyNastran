from typing import Optional
import numpy as np
from cpylog import SimpleLogger
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt

def inertia_relief3(xyz: np.ndarray,
                    exterior_force: np.ndarray,
                    mass: np.ndarray,
                    g: Optional[float]=None,
                    xyz_cg: Optional[float]=None,
                    inertia_cg: Optional[np.ndarray]=None,
                    include_linear_inertia: bool=True,
                    include_angular_inertia: bool=True,
                    mass_units: str='slinch',
                    length_units: str='in',
                    force_units: str='lbf',
                    accel_units: Optional[str]=None,
                    inertia_units: Optional[str]=None,
                    case: str='',
                    show: bool=True,
                    locations: dict[str, float]=None,
                    plot_differential: bool=True,
                    plot_integrated: bool=True):
    """
    Parameters
    ----------
    xyz : (n, 3) float array
        the xyz station
    xyz_cg : (3, ) float array
        the center of gravity of the vehicle
    exterior_force : (n, 3) float array
        differential force
    exterior_moment : (n, 3) float array
        differential moment
    mass : (n, ) float array
        differential mass
        mass has units of slinch, slug, kg, lbm
    inertia_cg : (n, ) float array
        differential inertia about the cg
        inertia has units of slinch-in^2, slug-ft^2, kg-m^2
    g : float; default=None -> sum(F)/m
        the gravity to use for inertia relief
        probably should be positive (a=F/m) to apply inertia relief
        in the opposite direction
        slinch: 12*32.174=386.088
    mass_units / length_units / force_units: str
        used for prints
    accel_units / inertia_units : Optional[str]
        automatically calculated; used for prints
    locations : dict[str, float] | None
        locations to identify key locations (e.g., nose, wing LE, tail)
    """
    log = SimpleLogger()
    if locations is None:
        locations = {}
    assert isinstance(locations, dict), locations
    assert isinstance(xyz, np.ndarray), type(xyz)
    assert isinstance(mass, np.ndarray), type(mass)
    assert isinstance(exterior_force, np.ndarray), type(exterior_force)
    nx = len(mass)

    moment_units, accel_units, alpha_units, inertia_units = _get_units(
        mass_units, length_units, force_units,
        accel_units=accel_units,
        inertia_units=inertia_units)

    sum_mass = mass.sum()
    xyz_cgi, dxyz_cg, inertia_cgi, sum_inertia_cg = _get_mass_cg_inertia3(
        xyz, xyz_cg, mass, sum_mass,
        inertia_cg=inertia_cg)
    del inertia_cg, xyz_cg

    exterior_moment = exterior_force * dxyz_cg  # fintegrate(exterior_force, x=x)
    sum_exterior_force = exterior_force.sum(axis=0)
    sum_exterior_moment = exterior_moment.sum(axis=0)
    assert len(sum_exterior_force) == 3, sum_exterior_force
    assert len(sum_exterior_moment) == 3, sum_exterior_moment

    #  F = m * a
    #  a = F / m
    exterior_accel = sum_exterior_force / sum_mass
    assert len(exterior_accel) == 3, exterior_accel

    # Mcg = Icg * alpha_cg
    # rxx, ryy, rzz, rxy, rxz, ryz
    #print(inertia_cgi.shape)
    #--------------------------------------------
    ixx = inertia_cgi[:, 0]
    iyy = inertia_cgi[:, 1]
    izz = inertia_cgi[:, 2]
    ixy = inertia_cgi[:, 3]
    ixz = inertia_cgi[:, 4]
    iyz = inertia_cgi[:, 5]
    #ixx, iyy, izz, ixy, ixz, iyz = inertia_cgi
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = sum_inertia_cg

    # adding a small amount of noise allows us to invert Imat
    tol = 1e-16
    Imat = np.array([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz],
    ], dtype='float64') + tol * np.eye(3)

    # adding a small amount of noise allows us to invert Imati
    Imati = np.zeros((nx, 3, 3), dtype='float64')
    Imati[:, 0, 0] = ixx + tol
    Imati[:, 0, 1] = Imati[:, 1, 0] = ixy
    Imati[:, 0, 2] = Imati[:, 2, 0] = ixz
    Imati[:, 1, 1] = iyy + tol
    Imati[:, 1, 2] = Imati[:, 2, 1] = iyz
    Imati[:, 2, 2] = izz + tol
    #print(Imat)
    assert Imat.shape == (3, 3), Imat.shape
    #--------------------------------------------
    # [Mx]    [Ixx, Ixy, Ixz] [alpha_x]
    # [My]  = [Ixy, Iyy, Iyz] [alpha_y]
    # [Mz]    [Ixz, Iyz, Izz] [alpha_z]
    #
    # [Ixx, Ixy, Ixz]^-1 [Mx]   [alpha_x]
    # [Ixy, Iyy, Iyz]    [My] = [alpha_y]
    # [Ixz, Iyz, Izz]    [Mz]   [alpha_z]
    sum_exterior_moment = sum_exterior_moment.reshape(3, 1)
    exterior_alpha_cg = np.linalg.solve(Imat, sum_exterior_moment) #/ sum_inertia_cg
    sum_exterior_moment = sum_exterior_moment.flatten()
    exterior_alpha_cg = exterior_alpha_cg.flatten()
    assert len(exterior_alpha_cg) == 3, exterior_alpha_cg

    if g is not None:
        accel = g

    if include_linear_inertia:
        # Fi = mi * a
        inertial_linear_force = -mass[:, np.newaxis] * exterior_accel
        #inertial_linear_moment = fintegrate(inertial_linear_force, x=x)
        inertial_linear_moment = inertial_linear_force * dxyz_cg
    else:
        inertial_linear_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
        inertial_linear_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    assert inertial_linear_force.shape == (nx, 3), inertial_linear_force.shape
    assert inertial_linear_moment.shape == (nx, 3), inertial_linear_moment.shape

    #--------
    # Mx_cg            alpha_cg
    # My_cg = [Imat] @ alpha_cg
    # Mz_cg            alpha_cg
    linear_force = exterior_force + inertial_linear_force
    linear_moment = exterior_moment + inertial_linear_moment
    sum_linear_force = linear_force.sum(axis=0)
    sum_linear_moment = linear_moment.sum(axis=0)
    assert len(sum_linear_force) == 3, sum_linear_force
    assert len(sum_linear_moment) == 3, sum_linear_moment

    #linear_alpha_cg = sum_linear_moment / sum_inertia_cg
    linear_alpha_cg = np.linalg.solve(Imat, sum_linear_moment).flatten() #/ sum_inertia_cg

    inertial_angular_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    log.info(f'dxyz_cg={type(dxyz_cg)} dxyz_cg.shape={str(dxyz_cg.shape)}')
    if include_angular_inertia:
        # Mcgi = -Icgi * alpha_cg
        icg_positive = np.abs(dxyz_cg) > 0
        #icg_positive = dxcg.abs() > 0
        assert linear_alpha_cg.shape == (3,), linear_alpha_cg.shape
        #inertial_angular_moment = -inertia_cgi * linear_alpha_cg
        inertial_angular_moment = -Imati @ linear_alpha_cg
        inertial_angular_force[icg_positive] = inertial_angular_moment[icg_positive] / dxyz_cg[icg_positive]
    else:
        inertial_angular_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    #-------------------------------------------------------------------------------------
    total_force = exterior_force + inertial_linear_force + inertial_angular_force
    total_moment = exterior_moment + inertial_linear_moment + inertial_angular_moment
    total_force_sum = total_force.sum(axis=0)
    total_moment_sum = total_moment.sum(axis=0)

    total_accel = total_force_sum / sum_mass
    #total_alpha_cg = total_moment.sum(axis=0) / sum_inertia_cg
    IImat = Imat[np.newaxis, :, :]
    #print(IImat.shape)
    #print(total_moment_sum.shape)
    total_alpha_cg = np.linalg.solve(Imat, total_moment_sum) #/ sum_inertia_cg

    assert len(total_force) == nx
    assert len(total_moment) == nx

    assert len(total_accel) == 3, total_accel
    assert len(total_alpha_cg) == 3, total_alpha_cg
    # print(f'  total_accel = {total_accel} {accel_units}')
    # print(f'  total_alpha_cg = {total_alpha_cg} {accel_units}')

    data_balanced = {
        f'xyz_cg ({length_units})': xyz_cgi,
        f'mass ({mass_units})': sum_mass,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,

        f'total_force ({force_units})': total_force_sum,
        f'total_moment_cg ({moment_units})': total_moment_sum,

        f'exterior_accel ({accel_units})': exterior_accel,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_balanced, indent='  ')
    log.info('Summary Balanced Loads')
    log.info('\n' + msg)

    print('done')


def inertia_relief1(x: np.ndarray,
                    exterior_force: np.ndarray,
                    mass: np.ndarray,
                    g: Optional[float]=None,
                    x_cg: Optional[float]=None,
                    inertia_cg: Optional[np.ndarray]=None,
                    include_linear_inertia: bool=True,
                    include_angular_inertia: bool=True,
                    mass_units: str='slinch',
                    length_units: str='in',
                    force_units: str='lbf',
                    accel_units: Optional[str]=None,
                    inertia_units: Optional[str]=None,
                    case: str='',
                    show: bool=True,
                    locations: dict[str, float]=None,
                    plot_differential: bool=True,
                    plot_integrated: bool=True):
    """
    Parameters
    ----------
    x : (n, ) float array
        the xyz station
    x_cg : float
        the center of gravity of the vehicle
    exterior_force : (n, ) float array
        differential force
    exterior_moment : (n, ) float array
        differential moment
    mass : (n, ) float array
        differential mass
        mass has units of slinch, slug, kg, lbm
    inertia_cg : (n, ) float array
        differential inertia about the cg
        inertia has units of slinch-in^2, slug-ft^2, kg-m^2
    g : float; default=None -> sum(F)/m
        the gravity to use for inertia relief
        probably should be positive (a=F/m) to apply inertia relief
        in the opposite direction
        slinch: 12*32.174=386.088
    mass_units / length_units / force_units: str
        used for prints
    accel_units / inertia_units : Optional[str]
        automatically calculated; used for prints
    locations : dict[str, float] | None
        locations to identify key locations (e.g., nose, wing LE, tail)
    """
    if locations is None:
        locations = {}
    assert isinstance(locations, dict), locations
    assert isinstance(x, np.ndarray), type(x)
    assert isinstance(mass, np.ndarray), type(mass)
    assert isinstance(exterior_force, np.ndarray), type(exterior_force)
    assert x.ndim == 1, x
    assert mass.ndim == 1, mass.shape

    moment_units, accel_units, alpha_units, inertia_units = _get_units(
        mass_units, length_units, force_units,
        accel_units=accel_units,
        inertia_units=inertia_units)

    sum_mass = mass.sum()
    xcgi, dxcg, inertia_cgi, sum_inertia_cg = _get_mass_cg_inertia1(
        x, x_cg, mass, sum_mass,
        inertia_cg=inertia_cg)

    exterior_moment = exterior_force * dxcg # fintegrate(exterior_force, x=x)
    sum_exterior_force = exterior_force.sum()
    sum_exterior_moment = exterior_moment.sum()

    #  F = m * a
    #  a = F / m
    exterior_accel = sum_exterior_force / sum_mass

    # Mcg = Icg * alpha_cg
    exterior_alpha_cg = sum_exterior_moment / sum_inertia_cg

    data = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'exterior_accel ({accel_units})': exterior_accel,
        f'g ({accel_units})': g,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
    }
    msg = _write_summary(data, indent='  ')

    word = ''
    if case:
        word = f' for {case}'
    print(f'Input Summary (about CG){word}')
    print(msg)

    if g is not None:
        accel = g

    if include_linear_inertia:
        # Fi = mi * a
        inertial_linear_force = -mass * exterior_accel
        #inertial_linear_moment = fintegrate(inertial_linear_force, x=x)
        inertial_linear_moment = inertial_linear_force * dxcg
    else:
        inertial_linear_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
        inertial_linear_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    # Mcg = Icg * alpha_cg
    linear_force = exterior_force + inertial_linear_force
    linear_moment = exterior_moment + inertial_linear_moment
    sum_linear_force = linear_force.sum()
    sum_linear_moment = linear_moment.sum()
    linear_alpha_cg = sum_linear_moment / sum_inertia_cg

    #---------------------------------------------------------------
    data = {
        #f'xcg ({length_units})': xcgi,
        #f'mass ({mass_units})': total_mass,
        #f'accel ({accel_units})': accel,
        f'F_linear ({force_units})': sum_linear_force,
        f'M_linear ({moment_units})': sum_linear_moment,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'linear_alpha_cg ({alpha_units})': linear_alpha_cg,
    }
    msg = _write_summary(data, indent='  ')


    print('Summary after Linear Inertia Applied')
    print(msg)
    inertial_angular_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    print('dxcg =', type(dxcg), dxcg)
    if include_angular_inertia:
        # Mcgi = -Icgi * alpha_cg
        icg_positive = np.abs(dxcg) > 0
        #icg_positive = dxcg.abs() > 0
        inertial_angular_moment = -inertia_cgi * linear_alpha_cg
        inertial_angular_force[icg_positive] = inertial_angular_moment[icg_positive] / dxcg[icg_positive]
    else:
        inertial_angular_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    #-------------------------------------------------------------------------------------
    total_force = exterior_force + inertial_linear_force + inertial_angular_force
    total_moment = exterior_moment + inertial_linear_moment + inertial_angular_moment
    total_accel = total_force.sum() / sum_mass
    total_alpha_cg = total_moment.sum() / sum_inertia_cg
    nx = len(x)
    assert len(total_force) == nx
    assert len(total_moment) == nx

    data_linear = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_accel ({accel_units})': exterior_accel,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_linear, indent='  ')
    print('Linear Summary')
    print(msg)

    data_angle = {
        f'xcg ({length_units})': xcgi,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
    }
    msg = _write_summary(data_angle, indent='  ')
    print('Angular Summary')
    print(msg)

    data_balanced = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,

        f'total_force ({force_units})': total_force.sum(),
        f'total_moment_cg ({moment_units})': total_moment.sum(),

        f'exterior_accel ({accel_units})': exterior_accel,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_balanced, indent='  ')
    print('Summary Balanced Loads')
    print(msg)
    #print(f'  xcg  ({length_units}) = {xcgi:.3f}')
    #print(f'  inertia_cg ({inertia_units}) = {total_inertia_cg:.3f}')
    #print(f'  exterior_moment_cg ({force_units}) = {total_exterior_moment:.3f}')
    #print(f'  alpha_cg ({accel_units}) = {total_alpha_cg:.3f}')
    out = (
        exterior_force,  inertial_linear_force, inertial_angular_force, total_force,
        exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
    )
    if plot_differential:
        plot_inertia(x,
                     exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                     exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                     locations,
                     exterior_accel, exterior_alpha_cg,
                     accel_units=accel_units,
                     integrate=False,
                     case=case)
    if plot_integrated:
        plot_inertia(x,
                     exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                     exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                     locations,
                     exterior_accel, exterior_alpha_cg,
                     accel_units=accel_units,
                     integrate=True,
                     case=case)
    if show:
        plt.show()
    return out

def _get_units(mass_units: str, length_units: str,
               force_units: str,
               accel_units: Optional[str]=None,
               inertia_units: Optional[str]=None,
               ) -> tuple[str, str, str, str]:
    if inertia_units is None:
        inertia_units = f'{mass_units}-{length_units}^2'
    if accel_units is None:
        accel_units = f'{length_units}/s^2'
    moment_units = f'{length_units}-{force_units}'
    alpha_units = 'rad/s^2'
    return moment_units, accel_units, alpha_units, inertia_units

def _get_mass_cg_inertia1(x: np.ndarray,
                          x_cg: np.ndarray,
                          mass: np.ndarray,
                          sum_mass: float,
                          inertia_cg=None):
    assert x.ndim == 1, x
    assert mass.ndim == 1, mass.shape
    x_cgi = _get_x_cg(x, x_cg, mass, sum_mass)
    dx_cg, inertia_cgi = _get_inertia1(x, x_cgi, mass, inertia_cg)
    del inertia_cg
    sum_inertia_cg = inertia_cgi.sum()
    return x_cgi, dx_cg, inertia_cgi, sum_inertia_cg

def _get_mass_cg_inertia3(xyz: np.ndarray,
                          xyz_cg: np.ndarray,
                          mass: np.ndarray,
                          sum_mass: float,
                          inertia_cg=None):
    assert xyz.ndim == 2, xyz
    assert mass.ndim == 1, mass.shape
    xyz_cgi = _get_xyz_cg(xyz, xyz_cg, mass, sum_mass)
    assert len(xyz_cgi) == 3, xyz_cgi
    dxyz_cg, inertia_cgi = _get_inertia3(xyz, xyz_cgi, mass, inertia_cg)

    sum_inertia_cg = inertia_cgi.sum(axis=0)
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = sum_inertia_cg
    assert len(sum_inertia_cg) == 6, sum_inertia_cg
    return xyz_cgi, dxyz_cg, inertia_cgi, sum_inertia_cg


def _get_x_cg(x: np.ndarray,
              x_cg: Optional[float],
              mass: np.ndarray,
              sum_mass: float) -> float:
    if x_cg is None:
        x_cgi = (x * mass).sum() / sum_mass
    else:
        assert isinstance(x_cg, float), x_cg
        x_cgi = x_cg
    return x_cgi

def _get_xyz_cg(xyz: np.ndarray,
                xyz_cg: Optional[float],
                mass: np.ndarray,
                sum_mass: float) -> float:
    if xyz_cg is None:
        xyz_cgi = (xyz * mass[:, np.newaxis]).sum(axis=0) / sum_mass
    else:
        assert isinstance(xyz_cg, np.ndarray), xyz_cg
        xyz_cgi = xyz_cg
    assert len(xyz_cgi) == 3, xyz_cgi
    return xyz_cgi

def _get_inertia1(x: np.ndarray, xcg: float,
                  mass: np.ndarray,
                  inertia_cg: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    dxcg = x - xcg
    if inertia_cg is not None:
        inertia_cgi = inertia_cg
    else:
        inertia_cgi = mass * dxcg ** 2
    return dxcg, inertia_cgi

def _get_inertia3(xyz: np.ndarray,
                  xyz_cg: np.ndarray,
                  mass: np.ndarray,
                  inertia_cg: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    if inertia_cg is not None:
        inertia_cgi = inertia_cg
    else:
        dxyz_cg = xyz - xyz_cg[np.newaxis, :]
        dx = dxyz_cg[:, 0]
        dy = dxyz_cg[:, 1]
        dz = dxyz_cg[:, 2]
        dx2 = dx * dx
        dy2 = dy * dy
        dz2 = dz * dz
        rxx = dy2 + dz2
        ryy = dx2 + dz2
        rzz = dx2 + dy2
        rxy = dx * dy
        ryz = dy * dz
        rxz = dx * dz
        radius_gyration2 = np.column_stack([rxx, ryy, rzz, rxy, rxz, ryz])
        inertia_cgi = mass[:, np.newaxis] * radius_gyration2
    return dxyz_cg, inertia_cgi

def _write_summary(data: dict[str, float], indent: str='  ') -> str:
    len_max_name = max((len(key) for key in data))
    msg = ''
    for key, value in data.items():
        if isinstance(value, float):
            fmt = indent + '%%-%ds = %%.3f\n' % len_max_name
        else:
            fmt = indent + '%%-%ds = %%s\n' % len_max_name
        msgi = fmt % (key, value)
        msg += msgi
    return msg

def fintegrate(y, x=None) -> np.ndarray:
    if x[0] > x[-1]:
        sign = -1
    else:
        sign = 1
    xy = np.hstack([0., cumulative_trapezoid(y, x=sign*x)])
    return xy

def plot_inertia(x,
                 exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                 exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                 locations: dict[str, float],
                 accel: float, alpha: float,
                 accel_units:str='in/s^2',
                 alpha_units:str='rad/s^2',
                 length_units:str='in', force_units: str='lbf',
                 integrate: bool=True,
                 case: str='') -> None:
    fig = plt.figure()
    axes = fig.subplots(nrows=2, ncols=2)
    moment_units = f'{length_units}-{force_units}'

    word = ''
    if case:
        word = f': {case}'
    if integrate:
        #dword = ''
        exterior_force = fintegrate(exterior_force, x=x)
        inertial_linear_force = fintegrate(inertial_linear_force, x=x)
        inertial_angular_force = fintegrate(inertial_angular_force, x=x)
        exterior_moment = fintegrate(exterior_moment, x=x)
        inertial_linear_moment = fintegrate(inertial_linear_moment, x=x)
        inertial_angular_moment = fintegrate(inertial_angular_moment, x=x)
        total_moment = fintegrate(total_moment, x=x)
        fig.suptitle(f'Integrated Loads{word}')
    else:
        fig.suptitle(f'Differential Loads{word}')
        #dword = 'd'

    fma = exterior_force + inertial_linear_force
    mia = exterior_moment + inertial_linear_moment

    # force
    axes[0, 0].set_title(f'accel {accel:.6f} ({accel_units})')
    axes[0, 0].set_ylabel(f'Force ({force_units})')
    axes[0, 0].plot(x, exterior_force, label=f'F = exterior [{exterior_force.min():.0f}, {exterior_force.max():.0f}]', linestyle='--', )
    axes[0, 0].plot(x, inertial_linear_force, label=f'-ma = linear [{inertial_linear_force.min():.0f}, {inertial_linear_force.max():.0f}]', linestyle='--', )
    axes[0, 0].plot(x, fma, label=f'F-ma [{fma.min():.0f}, {fma.max():.0f}]')

    axes[0, 1].set_ylabel(f'Force ({force_units})')
    axes[0, 1].plot(x, fma, label=f'F-ma [{fma.min():.0f}, {fma.max():.0f}]', linestyle='--', )
    axes[0, 1].plot(x, inertial_angular_force, label=f'-I*alpha/dx = angular [{inertial_angular_force.min():.0f}, {inertial_angular_force.max():.0f}]', linestyle='--', )
    axes[0, 1].plot(x, total_force, label=f'total [{total_force.min():.0f}, {total_force.max():.0f}]')

    #----------------
    # moment
    axes[1, 0].set_title(f'alpha {alpha:.6f} ({alpha_units})')
    axes[1, 0].set_ylabel(f'Moment ({moment_units})')
    axes[1, 0].plot(x, exterior_moment, label=f'M = exterior [{exterior_moment.min():.3g}, {exterior_moment.max():.3g}]', linestyle='--', )
    axes[1, 0].plot(x, inertial_linear_moment, label=f'-I*alpha = linear [{inertial_linear_moment.min():.3g}, {inertial_linear_moment.max():.3g}]', linestyle='--', )
    axes[1, 0].plot(x, mia, label=f'M-I*alpha [{mia.min():.3g}, {mia.max():.3g}]')

    axes[1, 1].set_ylabel(f'Moment ({moment_units})')
    axes[1, 1].plot(x, mia, label=f'(F-ma)*dx [{mia.min():.3g}, {mia.max():.3g}]', linestyle='--', )
    axes[1, 1].plot(x, inertial_angular_moment, label=f'-I*alpha = angular [{inertial_angular_moment.min():.3g}, {inertial_angular_moment.max():.3g}]', linestyle='--', )
    axes[1, 1].plot(x, total_moment, label=f'total [{total_moment.min():.3g}, {total_moment.max():.3g}]')

    for ax in axes.ravel():
        for key, value in locations.items():
            ax.axvline(value, label=key, linestyle='--', color='k')
        ax.set_xlabel('Station (in)')
        ax.grid(True)
        ax.legend()

def test_inertia1():
    x = np.linspace(0., 10., num=51)
    mass = np.ones(x.shape, x.dtype)
    mass = np.ones(x.shape, x.dtype) + 0.1 * x
    #exterior_force = np.zeros(x.shape, x.dtype)
    #exterior_force[0] = 10.
    #exterior_force[1] = 5.
    #exterior_force = np.sin(x)
    exterior_force = mass * 3.
    #exterior_force = 2 * x

    inertia_relief1(
        x, exterior_force, mass, # inertia_cg=inertia_cg,
        xyz_cg=None,
        include_linear_inertia=True, include_angular_inertia=True,
        mass_units='slinch', length_units='in', force_units='lbf',
        plot_differential=True,
        plot_integrated=True,
        #integrate=False,
    )
    # inertia_relief(
    #     x, exterior_force, mass,
    #     xcg=None,
    #     include_linear_inertia=True, include_angular_inertia=True,
    #     mass_units='slinch', length_units='in', force_units='lbf',
    #     plot_differential=True,
    #     plot_integrated=True,
    #     #integrate=True,
    # )
    plt.show()

def test_inertia3():
    nx = 51
    xyz = np.zeros((nx, 3), dtype='float64')
    x = np.linspace(0., 10., num=51)
    xyz[:, 0] = x
    mass = np.ones(x.shape, dtype=x.dtype)
    mass = np.ones(x.shape, dtype=x.dtype) + 0.1 * x

    exterior_force = np.zeros((nx, 3), dtype='float64')
    exterior_force[:, 0] = mass * 3.
    inertia_relief3(xyz, exterior_force, mass)

if __name__ == '__main__':
    #test_inertia1()
    test_inertia3()
