import numpy as np
np.set_printoptions(precision=3, threshold=1000, linewidth=1000, suppress=True)

def main():
    """
    https://www.reddit.com/r/fea/comments/k2m75f/whas_is_wrong_with_my_stiffness_matrix/
    https://www.math.purdue.edu/~caiz/math615/matlab_fem.pdf
    http://what-when-how.com/the-finite-element-method/fem-for-3d-solids-finite-element-method-part-2/
    http://web.archive.org/web/20170809053607/http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
    http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
    """
    fdtype = 'float64'
    #clc
    #clear all
    # 8 gauss points in total
    k = (1 / 3) ** 0.5
    xhi = np.array([-k,-k,-k,-k, k, k, k, k], dtype=fdtype)
    eta = np.array([-k,-k, k, k,-k,-k ,k, k], dtype=fdtype)
    my  = np.array([-k, k,-k, k,-k, k,-k, k], dtype=fdtype)
    ngp = 8 # 2 gauss points in each direction

    #xhi = np.array([0.])
    #eta = np.array([0.])
    #my = np.array( [0.])
    #ngp = 1 # gauss points at the centroid

    #coord = np.array([
       #[0., 0.,  0.],
       #[20., 0.,  0.],
       #[20., 20., 0.],
       #[0.,  20., 0.],
       #[0.,  0.,  20.],
       #[20., 0.,  20.],
       #[20., 20., 20.],
       #[0.,  20., 20], ])

    coord = np.array([
       [-1., -1., -1.],
       [1.,  -1., -1.],
       [1.,   1., -1.],
       [-1.,  1., -1.],

       [-1., -1.,  1.],
       [1.,  -1.,  1.],
       [1.,   1.,  1.],
       [-1.,  1.,  1.]], dtype=fdtype)
    ncoord = coord.shape[0]
    ndof = 3 * ncoord
    print('ncoord =', ncoord)
    print('ndof =', ndof)

    # original
    E = 210.e3
    #poisson = 0.3

    #E = 10000.  # was 210.e3
    #poisson = 0.2 # was 0.3

    # test2
    #E = 32.
    poisson = 1 / 3
    #poisson = 0.
    Dcoeff = E/((1 + poisson) * (1 - 2 * poisson))
    d11 = (1 - poisson)
    d12 = poisson
    d44 = ((1 - 2 * poisson)/2)
    D = Dcoeff * np.array([
       [d11, d12, d12,   0, 0, 0],
       [d12, d11, d12,   0, 0, 0],
       [d12, d12, d11,   0, 0, 0],
       [0, 0, 0,       d44,   0,   0],
       [0, 0, 0,         0, d44,   0],
       [0, 0, 0,         0,   0, d44],
    ], dtype=fdtype)
    N = np.zeros(8, dtype=fdtype)
    dNdXhi = np.zeros(8, dtype=fdtype)
    dNdEta = np.zeros(8, dtype=fdtype)
    dNdMy = np.zeros(8, dtype=fdtype)
    dNdMy = np.zeros(8, dtype=fdtype)
    dNdMy = np.zeros(8, dtype=fdtype)
    dNdMy = np.zeros(8, dtype=fdtype)
    dNdXhidEtadMy = np.zeros((3, 8), dtype=fdtype)

    zero = np.zeros(8, dtype=fdtype)
    Me = np.zeros((ndof, ndof), dtype=fdtype)
    Ke = np.zeros((ndof, ndof), dtype=fdtype)
    Ke2 = np.zeros((ndof, ndof), dtype=fdtype)

    # preallocate B and ix
    B = np.zeros((6, 24), dtype=fdtype)
    ix = np.arange(0, 24, 3)
    assert len(ix) == 8, len(ix)

    # if this was a real problem i would have to loop over each element
    # but i only have 1 element, so only need to loop over the Gauss points
    for i in range(ngp):
        shape_functions(xhi[i], eta[i], my[i],
                        N, dNdXhi, dNdEta, dNdMy)
        # find Jacobian for each Gauss point
        dNdXhidEtadMy[1-1,:] = dNdXhi
        dNdXhidEtadMy[2-1,:] = dNdEta
        dNdXhidEtadMy[3-1,:] = dNdMy

        # J = [dNdXhi*coord(:,1),dNdXhi*coord(:,1),dNdXhi*coord(:,1);
        J = dNdXhidEtadMy @ coord
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)

        # find shape functions derivate matrix wrt x, y & z
        dN_dxyz = Jinv @ dNdXhidEtadMy

        # construct B matrix
        Nx = dN_dxyz[1-1,:]
        Ny = dN_dxyz[2-1,:]
        Nz = dN_dxyz[3-1,:]
        #Nx *= 0; Ny *= 0; Nz *= 0
        #Nx += 1; Ny += 2; Nz += 3
        print('dN/dx = ', dN_dxyz)
        print(Nx.shape, Ny.shape, Nz.shape, zero.shape)

        #[array([-0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ]),
         #array([ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               #-0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ]),
         #array([ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               #-0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125]),
         #array([-0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125,
               #-0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ]),
         #array([ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               #-0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125,
               #-0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125]),
         #array([-0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125,
                #0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               #-0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125])]

        # top 3 rows
        B[0, ix] = Nx
        B[1, ix+1] = Ny
        B[2, ix+2] = Nz

        #http://web.archive.org/web/20170809053607/http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
        #http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf

        # [N1x 0   0   N2x 0   0   ... ]
        # [0   N1y 0   0   N2y 0   ... ]
        # [0   0   N1z 0   0   N2z ... ]
        # [N1y N1x 0   N2y N2x 0   ... ]
        # [0   N1z N1y 0   N2z N2y ... ]
        # [N1z 0   N1x N2z 0   N2x ... ]

        # column 1
        B[3, ix] = Ny
        B[5, ix] = Nz

        # column 2
        B[3, ix+1] = Nx
        B[4, ix+1] = Nz

        # column 3
        B[4, ix+2] = Ny
        B[5, ix+2] = Nx

        Bi = [
            np.hstack([Nx,   zero, zero]),
            np.hstack([zero, Ny,   zero]),
            np.hstack([zero, zero, Nz]),

            np.hstack([zero, Nz,   Ny]),
            np.hstack([Nz,   zero, Nx]),
            np.hstack([Ny,   Nx, zero]),
            #np.hstack([Ny,   Nx,   zero]),
            #np.hstack([zero, Nz,   Ny]),
            #np.hstack([Nz,   zero, Nx]),
        ]
        Ni = np.vstack([
            np.hstack([Nx,   zero, zero]),
            np.hstack([zero, Ny,   zero]),
            np.hstack([zero, zero, Nz]),
        ])

        #for bi in Bi:
            #print(bi.shape)
        B2 = np.vstack(Bi)
        assert B2.shape == (6, 24), B2.shape
        assert B.shape == B2.shape, B.shape

        # compute Gauss point contribution to element Ke matrix
        Ki = B.T @ D @ B
        Ki2 = B2.T @ D @ B2

        Mi = Ni.T @ Ni
        Ke2 += Ki2 * detJ
        Ke += Ki * detJ
        Me += Mi * detJ

    # node 1-4 dof 1-3 = 0:
    Kg = Ke
    Kconstrained = Kg[12:, :][:, 12:].copy()

    # apply boundary conditions on nodes 1-4 --> dofs 1-12
    for i in range(12): # 3*4
        Kg[:, i] = np.zeros(ndof, dtype=fdtype)
        Kg[i, :] = np.zeros(ndof, dtype=fdtype)
        Kg[i, i] = 1

    # on the right side, we'll apply uniaxial tension
    # nodes 5-8 (in the x degree of freedom)
    Fg = np.zeros((ndof, 1))
    #F[-1, 0] = 100e3 # old load
    idof = np.arange(12, ndof, 3)
    L = 2.
    Pi = 10000.
    Fg[idof, :] = Pi
    #print('F =', F.T)
    Fconstrained = Fg[12:, :].copy()

    Frotated = Fg.reshape(ncoord, 3)
    print('Frotated =\n', Frotated)

    # area of the cube is 1 in^2
    # length=1.0
    # the displacement is:
    # PL / AE = P / E
    u_expected = 4 * Pi * L / E
    print('u_expected = ', u_expected)
    #f = np.array([zeros(1,23) 100e3]).reshape(3, 1)
    u = np.linalg.solve(Kg, Fg)
    ucon = np.linalg.solve(Kconstrained, Fconstrained)
    ux = ucon[::3]
    uy = ucon[1::3]
    uz = ucon[2::3]

    #u /= 1.0e+14
    print('ucon =', ucon)
    print('u =', u)

    print('ux =', ux)
    print('uy =', uy)
    print('uz =', uz)
    # This yields the solution u:

    # u = 1.0e+14*[0 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1.0293, -1.0531, -1.0511, -1.0503, -1.3835, -1.3835, -0.3304, -0.3304, -1.3791, -1.3895, -0.3304, -0.3304,

    x = 1


def shape_functions(xhi, eta, my,
                    N, dNdXhi, dNdEta, dNdMy):
    # shape functions, 1/8=0.125
    N[1-1] = (1-xhi) * (1-eta) * (1-my)*0.125
    N[2-1] = (1+xhi) * (1-eta) * (1-my)*0.125
    N[3-1] = (1+xhi) * (1+eta) * (1-my)*0.125
    N[4-1] = (1-xhi) * (1+eta) * (1-my)*0.125
    N[5-1] = (1-xhi) * (1-eta) * (1+my)*0.125
    N[6-1] = (1+xhi) * (1-eta) * (1+my)*0.125
    N[7-1] = (1+xhi) * (1+eta) * (1+my)*0.125
    N[8-1] = (1-xhi) * (1+eta) * (1+my)*0.125

    #N1 = (1-xhi[i]) * (1-eta[i]) * (1-my[i])*0.125
    #dN1 / dxhi = a
    #(1 - xhi) -> -1
    #(1 + xhi) -> 1
    # derive shape functions wrt xhi - good
    dNdXhi[1-1] = -(1-eta) * (1-my)*0.125
    dNdXhi[2-1] =  (1-eta) * (1-my)*0.125
    dNdXhi[3-1] =  (1+eta) * (1-my)*0.125
    dNdXhi[4-1] = -(1+eta) * (1-my)*0.125
    dNdXhi[5-1] = -(1-eta) * (1+my)*0.125
    dNdXhi[6-1] =  (1-eta) * (1+my)*0.125
    dNdXhi[7-1] =  (1+eta) * (1+my)*0.125
    dNdXhi[8-1] = -(1+eta) * (1+my)*0.125

    # derive shape functions  wrt eta - good
    dNdEta[1-1] = -(1-xhi) * (1-my)*0.125
    dNdEta[2-1] = -(1+xhi) * (1-my)*0.125
    dNdEta[3-1] =  (1+xhi) * (1-my)*0.125
    dNdEta[4-1] =  (1-xhi) * (1-my)*0.125
    dNdEta[5-1] = -(1-xhi) * (1+my)*0.125
    dNdEta[6-1] = -(1+xhi) * (1+my)*0.125
    dNdEta[7-1] =  (1+xhi) * (1+my)*0.125
    dNdEta[8-1] =  (1-xhi) * (1+my)*0.125

    # derive shape functions wrt my - good
    dNdMy[1-1]  = -(1-xhi) * (1-eta)*0.125
    dNdMy[2-1]  = -(1+xhi) * (1-eta)*0.125
    dNdMy[3-1]  = -(1+xhi) * (1+eta)*0.125
    dNdMy[4-1]  = -(1-xhi) * (1+eta)*0.125
    dNdMy[5-1]  =  (1-xhi) * (1-eta)*0.125
    dNdMy[6-1]  =  (1+xhi) * (1-eta)*0.125
    dNdMy[7-1]  =  (1+xhi) * (1+eta)*0.125
    dNdMy[8-1]  =  (1-xhi) * (1+eta)*0.125

    return

if __name__ == '__main__':  # pragma: no cover
    main()
