import os
import numpy as np
import matplotlib.pyplot as plt

def plot_equivalent_lamina_vs_theta(pcomp, mid_ref, thetad,
                                    plot=False, show=False,
                                    close=True, png_filename=None):
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
        Sbar = pcomp.get_Sbar_matrix(mid_ref, thetai)
        Qbar = pcomp.get_Qbar_matrix(mid_ref, thetai)
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
        ax.plot(thetad, Ex/Ex.max(), label='Ex=%g' % Ex.max())
        ax.plot(thetad, Ey/Ey.max(), label='Ey=%g' % Ey.max())
        ax.plot(thetad, Gxy/Gxy.max(), label='Gxy=%g' % Gxy.max())
        ax.plot(thetad, Ex/e22, label='Ex/E2=%g' % Ex.max())
        ax.plot(thetad, Ey, label='Ey=%g' % Ey.max())
        ax.plot(thetad, Gxy/g12, label='Gxy/G12=%g' % Gxy.max())
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
        if png_filename:
            base, ext = os.path.splitext(png_filename)
            png_filename1 = base + '_stiffness' + ext
            png_filename2 = base + '_nu' + ext
            fig1.savefig(png_filename1)
            fig2.savefig(png_filename2)
        if show:
            plt.show()
    return out
