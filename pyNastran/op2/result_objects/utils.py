import numpy as np


def real_modes_to_omega_freq(eigns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    omega is also called radians
    cycles is also called frequency
    """
    omega_radians = np.sqrt(np.abs(eigns))
    abs_freqs = omega_radians / (2 * np.pi)
    return omega_radians, abs_freqs


def complex_damping_frequency(eigr: np.ndarray,
                              eigi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """eigenvalue = eigr + eigi*1j"""
    assert isinstance(eigr, np.ndarray), eigr
    eigr[eigr == -0.] = 0.
    eigi[eigi == -0.] = 0.
    damping = np.zeros(len(eigr), dtype=eigr.dtype)
    if 0:  # pragma: no cover
        denom = np.sqrt(eigr ** 2 + eigi ** 2)
        inonzero = np.where(denom != 0)[0]
        if len(inonzero):
            damping[inonzero] = -eigr[inonzero] / denom[inonzero]

        # not sure
        abs_freqs = np.sqrt(np.abs(eigi)) / (2 * np.pi)
    else:
        # flutter
        # eig = omega*zeta + omega*1j = eigr + eigi*1j
        # freq = eigi/(2*pi)
        # zeta = 2*eigr/eigi
        abs_freqs = abs(eigi) / (2 * np.pi)
        inonzero = np.where(eigi != 0)[0]
        if len(inonzero):
            damping[inonzero] = 2 * eigr[inonzero] / eigi[inonzero]
    return damping, abs_freqs
