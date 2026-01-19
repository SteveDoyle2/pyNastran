from typing import Callable
import functools
import numpy as np


def underflow(func: Callable):
    """A decorator that handles underflow by upcasting the data."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except FloatingPointError:
            args2 = []
            is_float32 = False
            for arg in args:
                if isinstance(arg, np.ndarray):
                    if arg.dtype.name == 'float32':
                        is_float32 = True
                        args2.append(arg.astype('float64'))
                    # elif arg.dtype.name == 'complex64':
                    #     is_float32 = True
                    #     args2.append(arg.astype('complex128'))
                    else:
                        args2.append(arg)
                else:
                    args2.append(arg)
            if not is_float32:
                raise
            results = func(*args2, **kwargs)
            if isinstance(results, np.ndarray):
                result = results.astype('float32')
            elif isinstance(results, tuple):
                result_list = []
                for arg in results:
                    if isinstance(arg, np.ndarray):
                        if arg.dtype.name == 'float64':
                            result_list.append(arg.astype('float32'))
                        else:
                            result_list.append(arg)
                    else:
                        result_list.append(arg)
                result = tuple(result_list)
            else:  # pragma: no cover
                raise TypeError(results)
        return result
    return wrapper


@underflow
def max_shear_2d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 txy: np.ndarray,
                 is_stress: bool) -> np.ndarray:
    """
    per nastran elements.pdf

    not verified for stress/strain
    same as Tresca?
    """
    if is_stress:
        max_sheari = np.sqrt((oxx - oyy) ** 2 / 4 + txy ** 2)
    else:
        # max_shear = np.sqrt((oxx - oyy) ** 2 + txy ** 2)
        max_sheari = np.sqrt((oxx - oyy) ** 2 + txy ** 2)
    return max_sheari

@underflow
def max_shear(omax: np.ndarray,
              omin: np.ndarray, is_stress: bool) -> np.ndarray:
    """
    not verified for stress/strain
    same as Tresca?
    """
    if is_stress:
        max_sheari = np.abs(omax - omin) / 2.
    else:
        max_sheari = np.abs(omax - omin)
    return max_sheari


@underflow
def von_mises_2d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 txy: np.ndarray,
                 is_stress: bool) -> np.ndarray:
    """verified for stress/strain"""
    if is_stress:
        ovm = np.sqrt(oxx**2 + oyy**2 - oxx*oyy + 3*(txy**2) )
    else:
        # center = (oxx + oyy) / 2
        # radius = np.sqrt((oxx - oyy) ** 2 / 4 + txy ** 2 / 4)
        # eig_max = center + radius
        # eig_min = center - radius
        # return von_mises_2d(eig_max, eig_min, txy*0, is_stress=True)
        ovm = np.sqrt(4/9 * (oxx**2 + oyy**2 - oxx*oyy) + (txy**2)/3 )
        assert np.all(np.isfinite(ovm))
    return ovm


def ovm_shear_2d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 txy: np.ndarray,
                 o1: np.ndarray,
                 o2: np.ndarray,
                 is_von_mises: bool,
                 is_stress: bool):
    if is_von_mises:
        ovm_shear = von_mises_2d(oxx, oyy, txy, is_stress)
    else:
        # max shear
        ovm_shear = max_shear_2d(oxx, oyy, txy, is_stress)
    return ovm_shear

@underflow
def principal_2d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 txy: np.ndarray,
                 is_stress: bool) -> tuple[np.ndarray, np.ndarray]:
    """
    https://www.simscale.com/docs/simwiki/fea-finite-element-analysis/principal-stress-and-principal-strain/
    """
    if is_stress:
        center = (oxx + oyy) / 2
        radius = np.sqrt((oxx - oyy) ** 2 / 4 + txy ** 2)
    else:
        center = (oxx + oyy) / 2
        radius = 0.5 * np.sqrt((oxx - oyy) ** 2 + txy ** 2)
    eig_max = center + radius
    eig_min = center - radius
    return eig_max, eig_min


@underflow
def von_mises_3d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 ozz: np.ndarray,
                 txy: np.ndarray,
                 tyz: np.ndarray,
                 txz: np.ndarray,
                 is_stress: bool) -> np.ndarray:
    """
    verified for von mises stress/strain
    """
    if is_stress:
        # octehedral = 1 / 3 * np.sqrt(
        #     (oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2 +
        #     6. * (txy ** 2 + tyz ** 2 + txz ** 2))
        # ovm_shear = 3 * octehedral / 2 ** 0.5
        # ovm_shear = 1 / 2 ** 0.5 * np.sqrt(
        #     (oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2 +
        #     6. * (txy ** 2 + tyz ** 2 + txz ** 2))

        vm = np.sqrt(
            0.5 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 +(oxx - ozz) ** 2) +
            3 * (txy**2 + tyz**2 + txz**2))
    else:
        # octaedral = np.sqrt(
        #     1/9 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2) +
        #     1/6 * (txy ** 2 + tyz ** 2 + txz ** 2))
        vm = np.sqrt(
            2/9 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2) +
            1/3 * (txy**2 + tyz**2 + txz**2))
        # vm = np.sqrt(
        #     0.5 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2) +
        #     3 * (txy ** 2 + tyz ** 2 + txz ** 2))
    return vm


def ovm_shear_3d(oxx: np.ndarray, oyy: np.ndarray, ozz: np.ndarray,
                 txy: np.ndarray, tyz: np.ndarray, txz: np.ndarray,
                 o1: np.ndarray, o3: np.ndarray,
                 is_von_mises: bool,
                 is_stress: bool):
    """
    not verified for max shear stress/strain
    octehedral stress/strain not supported
    """
    if is_stress:
        if is_von_mises:
            # von mises
            # ovm_shear = np.sqrt((oxx - oyy)**2 + (oyy - ozz)**2 + (oxx - ozz)**2 +
            #                     3. * (txy**2 + tyz**2 + txz ** 2))
            # octehedral = 1/3 * np.sqrt((oxx - oyy)**2 + (oyy - ozz)**2 + (oxx - ozz)**2 +
            #                     6. * (txy**2 + tyz**2 + txz ** 2))
            # ovm_shear = 3 * octehedral / 2 ** 0.5
            ovm_shear = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz, is_stress)
        else:
            # max shear
            ovm_shear = (o1 - o3) / 2.
    else:
        if is_von_mises:
            ovm_shear = von_mises_3d(oxx, oyy, ozz, txy, tyz, txz, is_stress)
        else:
            # max shear
            ovm_shear = (o1 - o3) / 2.
    return ovm_shear


def principal_components_3d(ntimes: int, nelements_nnodes: int,
                            oxx: np.ndarray, oyy: np.ndarray, ozz: np.ndarray,
                            txy: np.ndarray, tyz: np.ndarray, txz: np.ndarray,
                            is_stress: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    is this validated for stress and strain?

    # seems we should use no negative signs on the cross terms
    https://www.engapplets.vt.edu/Mohr/java/nsfapplets/MohrCircles2-3D/Theory/theory.htm
    """
    # if self.is_stress:
    #     a = oxx + oyy + ozz
    #     b = oxx*oyy + oyy*ozz + oxx*ozz - (txy**2 + tyz**2 + txz**2)
    #     c = oxx*oyy*ozz + 2*txy*tyz*txz - oxx*tyz**2 - oyy*txz**2 - ozz*txy**2

    a_matrix = np.full((ntimes * nelements_nnodes, 3, 3), np.nan, dtype='float64')
    #print(a_matrix.shape, oxx.shape)
    a_matrix[:, 0, 0] = oxx
    a_matrix[:, 1, 1] = oyy
    a_matrix[:, 2, 2] = ozz

    # we're only filling the lower part of the A matrix
    if is_stress:
        a_matrix[:, 1, 0] = txy
        a_matrix[:, 2, 0] = txz
        a_matrix[:, 2, 1] = tyz
    else:
        a_matrix[:, 1, 0] = txy / 2.
        a_matrix[:, 2, 0] = txz / 2.
        a_matrix[:, 2, 1] = tyz / 2.

    eigs = np.linalg.eigvalsh(a_matrix)  # array = (..., M, M) array

    dtype = oxx.dtype
    o1 = eigs[:, 2].astype(dtype)  # max
    o2 = eigs[:, 1].astype(dtype)  # mid
    o3 = eigs[:, 0].astype(dtype)  # min
    return o1, o2, o3
