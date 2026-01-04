import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase

def get_freq_damping(model: BDF,
                     omegan: np.ndarray,
                     ndof_g: int) -> np.ndarray:
    """
    The different types of damping are:
     - linear/nonlinear: time domain; based on
           MAT GE, CBUSH B; viscous and/or structural
     - modal: freq domain; based on constant/interpolation
              based on modal damping
     - zeta: freq domain; based on constant/interpolation
             based on G, zeta, or Q
    This function is intended to handle the modal/zeta

    Paramaeters
    -----------
    model : BDF
        the model object
    omegan : (nmode,) float ndarray
        the natural frequencies
    """
    # Q, dynamic amplification factor
    # zeta = 1/(2*Q)

    # G, structural damping
    # zeta = G/2
    #
    # 3% structural damping
    G = 0.03
    zeta = G / 2
    zomegan2 = 2 * zeta * omegan

    # TODO: parse inputs for damping
    Cstr_gg = np.eye(ndof_g)
    np.fill_diagonal(Cstr_gg, zomegan2)
    return Cstr_gg, zomegan2
