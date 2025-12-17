import numpy as np


def max_shear(omax: np.ndarray,
              omin: np.ndarray) -> np.ndarray:
    """
    not verified for stress/strain
    same as Tresca?
    """
    max_shear = (omax - omin) / 2.
    return max_shear

def von_mises_2d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 txy: np.ndarray) -> np.ndarray:
    """not verified for stress/strain"""
    ovm = np.sqrt(oxx**2 + oyy**2 - oxx*oyy +3*(txy**2) )
    return ovm

def von_mises_3d(oxx: np.ndarray,
                 oyy: np.ndarray,
                 ozz: np.ndarray,
                 txy: np.ndarray,
                 tyz: np.ndarray,
                 txz: np.ndarray) -> np.ndarray:
    """not verified for stress/strain"""
    vm = np.sqrt(
        0.5 * ((oxx - oyy) ** 2 + (oyy - ozz) ** 2 +(oxx - ozz) ** 2) +
        3 * (txy**2 + tyz**2 + txz**2))
    return vm
