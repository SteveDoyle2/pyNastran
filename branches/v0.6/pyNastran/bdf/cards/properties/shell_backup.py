## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 21:32:56 2012

@author: steve
"""

from numpy import transpose, array
import copy

from math import sin, cos, radians

class ShellPropertyBackup(Property):
    def __init__(self):
        pass

    def ABD(self):
        raise NotImplementedError('ABD is not implemented')
        z = self.get_z_locations()
        z0 = z[:-1]
        z1 = z[1:]
        dz = z1-z0
        dz2 = (z1**2 - z0**2)/2
        dz3 = (z1**3 - z0**3)/3
        #Qb = self.Qbar()
        
        A = zeros((3,3), dtype='float64')
        B = A
        D = A
        
        Q = zeros(3,3)

        #T = array([ [cc, ss,  2*s*c],
        #            [ss, cc, -2*s*c],
        #            [-sc, sc, cc-ss]], dtype='float64')
        #Q = array([ [Q11, Q12, 0.],
        #            [Q12, Q22, 0.],
        #            [0., 0., Q66]], dtype='float64')
        #iT = T.inv()
        #Qbar = dot(dot(iT, Q), transpose(iT))
        #Q45 = array([ [Q44, 0.],
        #              [0., Q55]], dtype='float64')
        
        
        #S11 = 1/E1
        #S12 = -nu12/E1
        #S22 = 1/E2
        #S66 = 1/G12
        #S = array([ [S11, S12, 0.],
        #            [S12, S22, 0.],
        #            [0., 0., S66]], dtype='float64')
        #Sbar = dot(dot(iT, S), transpose(iT))
        #Sb11 = S11*cccc + (2*S12*S66)*ss*cc + S22*ssss
        #Sb22 = S11*ssss + (2*S12*S66)*ss*cc + S22*cccc
        #Sb66 = 2*(2*S11+2*S22+4*S12-S66)*ss*cc+S66*(ssss+cccc)
        #Sb12 = S12*(ssss+cccc) + (S11+S22-S66)*ss*cc
        #Sb16 = (2*S11-2*S12-S66)*s*ccc - (2*S22-2*S12-S66)*sss*c
        #Sb26 = (2*S11-2*S12-S66)*sss*c - (2*S22-2*S12-S66)*s*ccc
        
        #Qb16 = 0. # for isotropic
        #Qb26 = 0. # for isotropic

        for iply in range(self.nPlies()):
            theta = radians(self.Theta(iply))
            s = sin(theta)
            c = cos(theta)
            ss = s * s
            sss = ss * s
            ssss = ss * ss
            cc = c * c
            ccc = cc * c
            cccc = cc * cc

            mat = self.Material(iply)
            E = mat.E()
            G = mat.G()
            nu = mat.Nu()
            
            E11 = E
            E22 = E
            G12 = G
            G23 = G
            G31 = 0.01
            nu12 = nu
            nu21 = E2/E1 * nu12
            
            Q11 = E1     /(1-nu12*nu21)
            Q12 = nu21*E1/(1-nu12*nu21)
            Q22 = E2     /(1-nu12*nu21)
            Q44 = G23
            Q55 = G31
            Q66 = G12

            Qb11 = Q11*cccc + 2*(Q12+2Q66)*ss*cc + Q22*ssss
            Qb22 = Q11*ssss + 2*(Q12+2Q66)*ss*cc + Q22*cccc
            Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66) * ss*cc + Q66*(ssss + cccc)

            Qb12 = (Q11 - Q12 - 4*Q66)*ss*cc + Q12*(ssss+cccc)
            Qb16 = (Q11 - Q12 - 2*Q66)*s*ccc + (Q12-Q22-2*Q66)*sss*c
            Qb26 = (Q11 - Q12 - 2*Q66)*sss*c + (Q12-Q22-2*Q66)*s*ccc
            
            Qb44 = Q44*cc + Q55*ss
            Qb55 = Q44*ss + Q55*cc
            Qb45 = (Q55-Q44)*ss*cc
            Qb = array([ [Qb11, Qb12, Qb166],
                         [Qb12, Qb22, Qb26],
                         [Qb16, Qb26, Qb66]], dtype='float64')
            Qb45 = array([ [Q44, 0.],
                          [0., Q55]], dtype='float64')

            A += Qb * dz[iply]
            B += Qb * dz2[iply]
            D += Qb * dz3[iply]

            A45 += Qb45 * dz[iply]
            B45 += Qb45 * dz2[iply]
            D45 += Qb45 * dz3[iply]
            #N += Q * alpha * dz
            #M += 1/2. * Q * alpha * dz2
        
        if self.isSymmetrical():
            B = zeros(Q.shape)
            B45 = zeros(Q.shape)
        return A, B, D, A45, B45, C45

    def FGH(self):
        raise NotImplementedError('FGH is not implemented')
        (A, B, D, A45, B45, C45) = self.ABD()
        A = matrix(A)
        B = matrix(B)
        D = matrix(D)

        iD = D.inv()
        iA = A.inv()
        if max(B)==min(B) and B[0,0]==0.0:
            F = iA
            G = B
            H = iD
        else:
            F = inv(A-B*iD*B)
            DBiAB = D-B*iA*B
            iDBiAB = DBiAB.inv()
            G = -iA*B*iDBiAB
            H = iDBiAB
        return F, G, H

    def S(self):
        r"""
        Calculates the compliance matrix for a lamina
        \f[ [Q] = [S]^{-1}  \f]
        """
        return self.Q().inv()

    def ABDH(self):
        r"""
        tranforms load to strain/bending curvature taken at \f$ z=0 \f$
        \f[ \left[
          \begin{array}{c}
              Nx  \\
              Ny  \\
              Nz  \\
              Mx  \\
              My  \\
              Mz  \\
          \end{array} \right] =
          \left[
          \begin{array}{cc}
              A & B  \\
              B & D  \\
          \end{array} \right]
          \left[
          \begin{array}{c}
              \epsilon_{xx}  \\
              \epsilon_{yy}  \\
                \gamma_{xy}  \\
                \kappa_{xx}  \\
                \kappa_{yy}  \\
                \kappa_{xy}  \\
          \end{array} \right]
        \f]

        ::
          [Nx] = [            ] [ e_xx0    ]
          [Ny] = [  [A]   [B] ] [ e_yy0    ]
          [Nz] = [            ] [ gamma_xy0]
          [Mx] = [            ] [ k_xx0    ]
          [My] = [  [B]   [D] ] [ k_yy0    ]
          [Mz] = [            ] [ k_xy0    ]
        
        \f[ A_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k  -z_{k-1}   \right) = \Sigma_{k=1}^N (Q_{ij})_k t_k                                                                                    \f]
        \f[ B_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k^2-z_{k-1}^2 \right) = \Sigma_{k=1}^N (Q_{ij})_k                           \left( \overline{z} t_k                      \right)         \f]
        \f[ D_{ij} = \Sigma_{k=1}^N (\overline{Q_{ij}})_k \left( z_k^3-z_{k-1}^3 \right) = \Sigma_{k=1}^N (Q_{ij})_k                           \left( \overline{z}^2 t_k + \frac{t_k^3}{12} \right)         \f]
        \f[ H_{ij} =                                                                       \Sigma_{k=1}^N (Q_{ij})_k \left( t_k -\frac{4}{t^2} \left( \overline{z}^2 t_k + \frac{t_k^3}{12} \right) \right) \f]

        p. 138 of "Introduction to Composite Material Design"
        """
        raise NotImplementedError()
        A = zeros(9, 'd')
        B = copy.deepcopy(A)
        D = copy.deepcopy(A)
        H = copy.deepcopy(A)

        for (i, layer) in enumerate(self.layers):
            t = layer.t
            z0 = layer.z
            #z1 = z0+t
            zbar = (2 * z0 + t) / 2.  # (z1+z0)/2.
            Qraw = layer.Q()  # needs E11, E22, G12, G13, nu12, nu21
            qlayer = self.Qall(Qout, layer.thetad)
            A += qlayer * t
            B += qlayer * t * z

            Dfactor = t * zbar * zbar + t ** 3 / 12.
            D += qlayer * Dfactor
            H += qlayer * (t - 4. / t ** 2 * Dfactor)
        B = 0.5 * B
        D = D / 3.
        H = H * 5. / 4.

    def Qall(self, thetad):
        r"""
        Caculates the laminate tranformation  stiffness \f$ [Q]_{all} \f$
        \f[ [Q]_{all} = [T]^{-1} [Q] [R][T][R]^{-1}  \f]
        \f[ [Q]_{all} = [T]^{-1} [Q] [T]^{-T}        \f]

        p. 123 of "Introduction to Composite Material Design"
        """
        raise NotImplementedError()
        theta = radians(thetad)
        ct = cos(theta)
        c2t = ct * ct
        c3t = ct * c2t
        c4t = c2t * c2t

        st = sin(theta)
        s2t = st * st
        s3t = st * s2t
        s4t = s2t * s2t

        s2c2t = s2t * c2t
        #s4tpc4t = s4t+c4t

        Q11a = Q11 * c4t + 2 * (Q12 + 2 * Q66) * s2c2t + Q22 * s4t
        Q12a = (Q11 + Q22 - 4 * Q66) * s2c2t + Q12(s4t + c4t)
        Q22a = Q11 * s4t + 2 * (Q12 + 2 * Q66) * s2c2t + Q22 * c4t
        Q16a = (Q11-Q12-2*Q66)*st*c3t + (Q12-Q22+2*Q66)*s3t*ct
        Q26a = (Q11-Q12-2*Q66)*s3t*ct + (Q12-Q22+2*Q66)*st*c3t
        Q66a = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2c2t + Q66(s4t + c4t)
        Q44a = Q44 * c2t + Q55 * s2t
        Q55a = Q44 * s2t + Q55 * c2t
        Q45a = (Q55 - Q44) * st * ct
        return array([Q11a, Q12a, Q22a, Q16a, Q26a, Q66a, Q44a, Q55a, Q45a])

    def Q(self):
        r"""
        Calculates the stiffness matrix \f$ [Q] \f$ for a lamina
        .. todo:: is this done?
        p. 114 "Introduction to Composite Material Design"
        """
        raise NotImplementedError()
        nu12 = self.nu12
        E1 = self.E1()
        E2 = self.E2()
        delta = 1 - nu12 * nu21
        Q11 = E1 / delta
        Q12 = nu12 * E2 / delta
        Q22 = E2 / delta
        Q66 = G12
        Q44 = G23
        Q55 = G13
        Qout = (Q11, Q22, Q12, Q44, Q55, Q66)
        return Qout

    def T(self, theta):
        r"""
        calculates the Transformation matricies \$ [T] \$ and  \f$ [T]^{-1} \f$
        :param self:            the object pointer
        :param theta:           in radians...
        :returns Tinv:          the inverse transformation matrix
        :returns TinvTranspose: the transposed inverse transformation matrix
        .. todo:: document better

        tranformation matrix  \f$ [T] \f$
        \f[ [T] = \left[
          \begin{array}{ccc}
              m^2 & n^2 &  2mn    \\
              n^2 & m^2 & -2mn    \\
              -mn & mn  & m^2-n^2
          \end{array} \right]
        \f]
        
        ::
                   [ m^2  n^2        2mn]
          [T]    = [ n^2  m^2       -2mn]   # transformation matrix
                   [ -mn   mn    m^2-n^2]

        inverse transformation matrix \f$ [T]^{-1} \f$
        \f[ [T]^{-1} = \left[
          \begin{array}{ccc}
              m^2 & n^2 & -2mn    \\
              n^2 & m^2 &  2mn    \\
              mn  & -mn & m^2-n^2
          \end{array} \right]
        \f]
        
        ::
                   [ m^2  n^2       -2mn]
          [T]^-1 = [ n^2  m^2        2mn]   # inverse transform
                   [ mn   -mn    m^2-n^2]

        \f[ \left[
            \begin{array}{c}
                 \sigma_{xx}   \\
                 \sigma_{yy}   \\
                 \sigma_{xy}
            \end{array} \right]

             = [T]^{-1} [Q] [R][T]

            \left[ \begin{array}{c}
                 \epsilon_{xx} \\
                 \epsilon_{yy} \\
                 \frac{1}{2} \gamma_{xy}
            \end{array} \right]
        \f]

        p.119 "Introduction to Composite Material Design"
        """
        n = cos(theta)
        m = sin(theta)
        Tinv = zeros((6, 6), 'd')
        mm = m * m
        nn = n * n
        mn = m * n
        Tinv[0][0] = Tinv[1][1] = mm
        Tinv[1][0] = Tinv[0][1] = nn
        Tinv[0][2] = -2 * mn
        Tinv[1][2] = 2 * mn

        Tinv[2][0] = mn
        Tinv[2][1] = -mn
        Tinv[2][2] = mm - nn
        TinvT = transpose(Tinv)
        return (Tinv, TinvT)

class PCOMP_backup(self):
    def __init__(self):
        pass

    def D(self):
        D = zeros([3, 3])
        #isSym = self.isSymmetrical()
        for iply in xrange(len(self.plies)):
            theta = self.Theta(iply)
            #t    = self.Thickness(iply)
            mat = self.Material(iply)
            Di = mat.Dplate()
            transform = self.T(theta)
            D += Di * transform
        return D

