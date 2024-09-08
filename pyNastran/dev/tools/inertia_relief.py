from typing import Optional
import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt

def inertia_relief(x: np.ndarray,
                   exterior_force: np.ndarray,
                   mass: np.ndarray,
                   g: Optional[float]=None,
                   xcg: Optional[float]=None,
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
        the station
    xcg : float
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

    sum_mass = mass.sum()
    xcgi = _get_xcg(x, xcg, mass, sum_mass)
    del xcg

    dxcg, inertia_cgi = _get_inertia(x, xcgi, mass, inertia_cg)
    del inertia_cg

    sum_inertia_cg = inertia_cgi.sum()

    exterior_moment = exterior_force * dxcg # fintegrate(exterior_force, x=x)
    sum_exterior_force = exterior_force.sum()
    sum_exterior_moment = exterior_moment.sum()

    if inertia_units is None:
        inertia_units = f'{mass_units}-{length_units}^2'
    if accel_units is None:
        accel_units = f'{length_units}/s^2'
    moment_units = f'{length_units}-{force_units}'
    alpha_units = 'rad/s^2'

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

def _get_xcg(x: np.ndarray,
             xcg: Optional[float],
             mass: np.ndarray,
             sum_mass: float) -> float:
    if xcg is None:
        xcgi = (x * mass).sum() / sum_mass
    else:
        assert isinstance(xcg, float), xcg
        xcgi = xcg
    return xcgi

def _get_inertia(x: np.ndarray, xcg: float,
                 mass: np.ndarray,
                 inertia_cg: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    dxcg = x - xcg
    if inertia_cg is not None:
        inertia_cgi = inertia_cg
    else:
        inertia_cgi = mass * dxcg ** 2
    return dxcg, inertia_cgi

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

    inertia_relief(
        x, exterior_force, mass, # inertia_cg=inertia_cg,
        xcg=None,
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

if __name__ == '__main__':
    test_inertia1()
