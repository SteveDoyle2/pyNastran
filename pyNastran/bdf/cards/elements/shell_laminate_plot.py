import numpy as np
import matplotlib.pyplot as plt
from pyNastran.bdf.cards.properties.shell import (
    PCOMP, get_Qbar_matrix, get_Sbar_matrix)
from pyNastran.bdf.cards.materials import MAT8


def plot_equivalent_lamina_vs_theta(pcomp: PCOMP,
                                    mid_ref: MAT8,
                                    thetad: np.ndarray,
                                    plot: bool=False,
                                    show: bool=False,
                                    close: bool=True,
                                    png_filename_base=None):
    """plots a PCOMP mid vs. theta"""
    e22 = mid_ref.e22
    g12 = mid_ref.g12
    theta = np.radians(thetad)

    Ex = []
    Ey = []
    Gxy = []
    Q66 = []
    nu_xy = []
    for thetai in theta:
        Sbar = get_Sbar_matrix(mid_ref, thetai)
        Qbar = get_Qbar_matrix(mid_ref, thetai)
        Sbar = np.linalg.inv(Qbar)
        Exi = 1 / Sbar[0, 0]
        Eyi = 1 / Sbar[1, 1]
        Gxyi = 1 / Sbar[2, 2]
        Q66i = Qbar[2, 2]
        nu_xyi = -Sbar[0, 1] * Exi

        #Gxyi = 1 / Q66i
        Ex.append(Exi)
        Ey.append(Eyi)
        Gxy.append(Gxyi)
        Q66.append(Q66i)
        nu_xy.append(nu_xyi)
    Ex = np.array(Ex)
    Ey = np.array(Ey)
    Gxy = np.array(Gxy)
    Q66 = np.array(Q66)
    nu_xy = np.array(nu_xy)
    out = {
        'Ex' : Ex,
        'Ey' : Ey,
        'Gxy' : Gxy,
        'nu_xy' : nu_xy,
    }

    min_max_theta = [thetad.min(), thetad.max()]

    if plot:
        #from pyNastran.gui.matplotlib_backend import matplotlib_backend
        #import matplotlib
        #matplotlib.use(matplotlib_backend)
        fig1 = plt.figure(1)
        ax = fig1.gca()

        ax.plot(thetad, Q66/Q66.max(), label='Q66=%g' % Q66.max())
        if Ex.max() != 0.:
            ax.plot(thetad, Ex/Ex.max(), label='Ex/|Ex|=%g' % Ex.max())
        if Ey.max() != 0.:
            ax.plot(thetad, Ey/Ey.max(), label='Ey/|Ey|=%g' % Ey.max())
        if Gxy.max() != 0.:
            ax.plot(thetad, Gxy/Gxy.max(), label='Gxy=%g' % Gxy.max())
        if e22 != 0.:
            ax.plot(thetad, Ex/e22, label='Ex/E2=%g' % Ex.max())
        if g12 != 0.:
            ax.plot(thetad, Gxy/g12, label='Gxy/G12=%g' % Gxy.max())
        ax.plot(thetad, Ex, label='Ex=%g' % Ex.max(), linestyle='--')
        ax.plot(thetad, Ey, label='Ey=%g' % Ey.max(), linestyle='--')
        ax.set_xlim(min_max_theta)
        ax.legend()
        ax.grid()
        ax.set_ylabel(r'$Q_{66}$')
        ax.set_xlabel(r'$\theta$')
        #----------------------------
        fig2 = plt.figure(2)
        ax = fig2.gca()
        ax.plot(thetad, nu_xy, label=r'$\nu_{xy}=\mathrm{%g}$' % nu_xy.max())
        ax.set_xlim(min_max_theta)
        ax.legend()
        ax.grid()
        ax.set_ylabel(r'$\nu_{xy}$')
        ax.set_xlabel(r'$\theta$')
        if png_filename_base:
            png_filename1 = png_filename_base + '_stiffness.png'
            png_filename2 = png_filename_base + '_nu.png'
            fig1.savefig(png_filename1)
            fig2.savefig(png_filename2)
        if show:
            plt.show()
    return out
