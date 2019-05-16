import numpy as np
import matplotlib.pyplot as plt

def plot_material_properties_vs_theta(pcomp, mid_ref, thetad,
                                      plot=False, show=False, png_filename=None):
    """plots a PCOMP mid vs. theta"""
    e22 = mid_ref.e22
    g12 = mid_ref.g12
    theta = np.radians(thetad)

    Ex = []
    Ey = []
    Gxy = []
    Q66 = []
    nu_xy = []
    for thetai in  theta:
        Qbar = pcomp.get_Q_matrix(mid_ref, thetai)
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

    min_max_theta = [thetad.min(), thetad.max()]

    if plot:
        #from pyNastran.gui.matplotlib_backend import matplotlib_backend
        #import matplotlib
        #matplotlib.use(matplotlib_backend)
        fig = plt.figure(1)
        ax = fig.gca()

        #ax.plot(thetad, Q66/Q66.max(), label='Q66=%g' % Q66.max())
        #ax.plot(thetad, Ex/Ex.max(), label='Ex=%g' % Ex.max())
        #ax.plot(thetad, Ey/Ey.max(), label='Ey=%g' % Ey.max())
        #ax.plot(thetad, Gxy/Gxy.max(), label='Gxy=%g' % Gxy.max())
        ax.plot(thetad, Ex/e22, label='Ex/E2=%g' % Ex.max())
        ax.plot(thetad, Ey, label='Ey=%g' % Ey.max())
        ax.plot(thetad, Gxy/g12, label='Gxy/G12=%g' % Gxy.max())
        #ax.set_xlim(min_max_theta)
        ax.legend()
        ax.grid()
        ax.set_xlabel('Q66')
        ax.set_xlabel('theta')
        #----------------------------
        fig = plt.figure(2)
        ax = fig.gca()
        ax.plot(thetad, nu_xy, label='\nu xy=%g' % nu_xy.max())
        ax.set_xlim(min_max_theta)
        ax.legend()
        ax.grid()
        ax.set_xlabel('\nu xy')
        ax.set_xlabel('theta')
        if png_filename:
            fig.savefig(png_filename)
        if show:
            plt.show()

