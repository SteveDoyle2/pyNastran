"""
defines the GridPointWeight class
"""
from __future__ import print_function
from six import StringIO

from numpy import zeros

from pyNastran.utils import object_attributes, object_methods


class GridPointWeight(object):
    def __init__(self):
        """
        .. seealso:: http://www.6dof.com/index.php?option=com_content&view=article&id=298:output-from-the-grid-point-weight-generator&catid=178:courses-and-trainings&Itemid=61
        """
        # The Grid Point Weight Generator (GPWG) module computes the rigid body
        # mass properties of an entire structure with respect to a user specified point and with
        # respect to the center of mass. Output from the module is requested by a PARAM
        # GRDPNT card in the Bulk Data Deck which specifies from which grid point mass
        # computations are to be referenced. Optionally, the absence of a specific grid point
        # (i.e. PARAM, GRDPNT, 0) automatically causes the origin of the basic
        # coordinate system to be utilized as a reference. The mass properties are initially
        # defined in the basic coordinate system. Subsequently, the mass properties are
        # transformed to principal mass axes and to principal inertia axes. The actual printout
        # is composed of several elements. These are:
        self.reference_point = None

        # M0 RIGID BODY MASS MATRIX IN BASIC COORDINATE SYSTEM
        # This is the rigid body mass matrix of the entire structure in the basic coordinate
        # system with respect to a reference point chosen by the analyst.
        self.MO = None

        # S TRANSFORMATION MATRIX FOR SCALAR MASS PARTITION
        # S is the transformation from the basic coordinate system to the set of principal axes
        # for the 3 x 3 scalar mass partition of the 6 x 6 mass matrix. The principal axes for
        # just the scalar partition are known as the principal mass axes.
        self.S = None

        self.mass = None

        # XC.G. YC.G. ZC.G.
        # It is possible in NASTRAN to assemble a structural model having different values of
        # mass in each coordinate direction at a grid point. This can arise, for example, by
        # assembling scalar mass components or from omitting some components by means of bar
        # element pin flags. Consequently three distinct mass systems are assembled one in each of
        # the three directions of the principal mass axes (the S system). This third tabulation
        # has five columns. The first column lists the axis direction in the S coordinates. The
        # second column lists the mass associated with the appropriate axis direction. The final
        # three columns list the x, y, and z coordinate distances from the reference point to the
        # center of mass for each of the three mass systems.
        self.cg = None

        # I(S) INERTIAS RELATIVE TO C.G.
        # This is the 3 x 3 mass moment of inertia partition with respect to the center of
        # gravity referred to the principal mass axes (the S system). This is not necessarily a
        # diagonal matrix because the determination of the S system does not involve second
        # moments. The values of inertias at the center of gravity are found from the values at
        # the reference point employing the parallel axes rule.
        self.IS = None

        # I(Q) PRINCIPAL INERTIAS
        # The principal moments of inertia at the center of gravity are displayed in matrix
        # form with reference to the Q system of axes. The Q system is obtained from an eigenvalue
        # analysis of the I(S) matrix.
        self.IQ = None

        # Q TRANSFORMATION MATRIX I(Q) = QT*IBAR(S)*Q
        # Q is the coordinate transformation between the S axes and the Q axes. IBAR(S) is the
        # same as I(s) except that the signs of the offdiagonal terms are reversed.
        self.Q = None

    def object_attributes(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def set_grid_point_weight(self, reference_point, MO, S, mass, cg, IS, IQ, Q):
        self.reference_point = reference_point
        self.MO = MO
        self.S = S
        self.mass = mass
        self.cg = cg
        self.IS = IS
        self.IQ = IQ
        self.Q = Q

    def get_stats(self, short=True):
        if self.reference_point is None:
            return ''
        if short:
            msg = ('GridPointWeight: ref_point=%s mass=%g; '
                   '[reference_point, M0, S, mass, cg, IS, IQ, Q]\n' % (
                       self.reference_point, self.mass.max()))
        else:
            msg = (
                'GridPointWeight:'
                '  reference_point=%s\n'
                '  mass=[%10g %10g %10g]\n'
                '  cg  =[%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n\n'

                '  IS  =[%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n\n'

                '  IQ  =[%10g %10s %10s]\n'
                '       [%10s %10g %10s]\n'
                '       [%10s %10s %10g]\n\n'

                '  Q  = [%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n'
                '       [%10g %10g %10g]\n' % (
                    self.reference_point, self.mass[0], self.mass[1], self.mass[2],
                    self.cg[0, 0], self.cg[0, 1], self.cg[0, 2],
                    self.cg[1, 0], self.cg[1, 1], self.cg[1, 2],
                    self.cg[2, 0], self.cg[2, 1], self.cg[2, 2],

                    self.IS[0, 0], self.IS[0, 1], self.IS[0, 2],
                    self.IS[1, 0], self.IS[1, 1], self.IS[1, 2],
                    self.IS[2, 0], self.IS[2, 1], self.IS[2, 2],

                    self.IQ[0], '', '',
                    '', self.IQ[1], '',
                    '', '', self.IQ[2],

                    self.Q[0, 0], self.Q[0, 1], self.Q[0, 2],
                    self.Q[1, 0], self.Q[1, 1], self.Q[1, 2],
                    self.Q[2, 0], self.Q[2, 1], self.Q[2, 2],
                    )
            )
        return msg

    def __repr__(self):
        f = StringIO()
        page_stamp = 'PAGE %i'
        page_num = 1
        self.write_f06(f, page_stamp, page_num)
        msg = f.getvalue()
        return msg

    def read_grid_point_weight(self, lines):
        """
         0-                                 REFERENCE POINT =        0
         1-                                           M O
         2- *  2.338885E+05  2.400601E-13 -7.020470E-15 -1.909968E-11  2.851745E+06 -5.229834E+07 *
         3- *  2.400601E-13  2.338885E+05 -2.520547E-13 -2.851745E+06  2.151812E-10  2.098475E+08 *
         4- * -7.020470E-15 -2.520547E-13  2.338885E+05  5.229834E+07 -2.098475E+08 -1.960403E-10 *
         5- * -1.909968E-11 -2.851745E+06  5.229834E+07  2.574524E+10 -5.566238E+10 -4.054256E+09 *
         6- *  2.851745E+06  2.151812E-10 -2.098475E+08 -5.566238E+10  2.097574E+11 -2.060162E+09 *
         7- * -5.229834E+07  2.098475E+08 -1.960403E-10 -4.054256E+09 -2.060162E+09  2.336812E+11 *
         8-                                           S
         9-                      *  1.000000E+00  0.000000E+00  0.000000E+00 *
        10-                      *  0.000000E+00  1.000000E+00  0.000000E+00 *
        11-                      *  0.000000E+00  0.000000E+00  1.000000E+00 *
        12-        DIRECTION
        13-     MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
        14-             X            2.338885E+05     -8.166148E-17  2.236038E+02  1.219276E+01
        15-             Y            2.338885E+05      8.972118E+02  9.200164E-16  1.219276E+01
        16-             Z            2.338885E+05      8.972118E+02  2.236038E+02 -8.381786E-16
        17-                                         I(S)
        18-                    *  1.401636E+10  8.739690E+09  1.495636E+09 *
                               *  8.739690E+09  2.144496E+10  1.422501E+09 *
                               *  1.495636E+09  1.422501E+09  3.370946E+10 *
                                                    I(Q)
                               *  3.389001E+10                             *
                               *                8.073297E+09               *
                               *                              2.720748E+10 *
                                                     Q
                               * -3.599259E-02 -8.305739E-01  5.557441E-01 *
                               * -8.850329E-02 -5.512702E-01 -8.296194E-01 *
                               *  9.954254E-01 -7.904533E-02 -5.366689E-02 *

        .. note::
            pyNastran's BDF mass_properties method uses the following (not
            totally correct as there technically isn't one xcg):

                               DIRECTION
                          MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
                                  X            mass              0.000000E+00  ycg           zcg
                                  Y            mass              xcg           0.000000E+00  zcg
                                  Z            nass              xcg           ycg           0.000000E+00

            The inertias are close to I(S), but not exact as the method doesn't
            use the mass matrix, but is close for sufficiently complex models.  The terms are:
                *  Ixx  Ixy  Ixz *
                *  Iyx  Iyy  Iyz *
                *  Izz  Izy  Izz *
            or inertia = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        """
        self.reference_point = int(lines[0].split('=')[1])
        assert lines[1] == 'M O', lines[1]

        self.MO = zeros((6, 6), dtype='float64')
        self.S  = zeros((3, 3), dtype='float64')

        self.mass = zeros(3, dtype='float64')
        self.cg = zeros((6, 6), dtype='float64')

        self.IS = zeros((3, 3), dtype='float64')
        self.IQ = zeros(3, dtype='float64')
        self.Q  = zeros((3, 3), dtype='float64')

        #========================================
        # MO
        n = 2
        for i in range(6):
            line = lines[n + i][1:-1]  # get rid of the * characters
            sline = line.split()
            for j in range(6):
                self.MO[i, j] = sline[j]
        #print("MO =", self.MO)
        n += i + 1

        #========================================
        # S
        assert lines[n] == 'S', lines[n]
        n += 1
        for i in range(3):
            line = lines[n + i][1:-1]  # get rid of the * characters
            sline = line.split()
            for j in range(3):
                self.S[i, j] = sline[j]
        #print("S =", self.S)
        n += i + 1

        #========================================
        assert lines[n] == 'DIRECTION', lines[n]
        n += 2
        for i in range(3):
            line = lines[n + i][1:]  # get rid of the * characters
            sline = line.split()

            self.mass[i] = sline[0]
            for j in range(3):
                self.cg[i, j] = sline[j + 1]

        #print("mass =", self.mass)
        #print("mass =", self.cg)
        n += 3

        #========================================
        assert lines[n] == 'I(S)', lines[n]
        n += 1
        for i in range(3):
            line = lines[n + i][1:-1]  # get rid of the * characters
            sline = line.split()
            for j in range(3):
                self.IS[i, j] = sline[j]
        #print("IS =", self.IS)
        n += i + 1

        #========================================
        assert lines[n] == 'I(Q)', lines[n]
        n += 1
        for i in range(3):
            sline = lines[n + i][1:-1].strip().split()  # get rid of the * characters
            self.IQ[i] = sline[0]
        #print("IQ =", self.IQ)
        n += i + 1

        #========================================
        # S
        assert lines[n] == 'Q', lines[n]
        n += 1
        for i in range(3):
            line = lines[n + i][1:-1]  # get rid of the * characters
            sline = line.split()
            for j in range(3):
                self.Q[i, j] = sline[j]
        #print("Q =", self.Q)
        n += i

    def write_f06(self, f06_file, page_stamp, page_num):
        """
        writes the f06

        Parameters
        ----------
        f06_file : file / StringIO
            a file-like object
        page_stamp : str
            the page formatter (e.g., 'PAGE %i')
        page_num : int
            the active page number

        Returns
        -------
        page_num : int
            the new page number
        """
        if self.reference_point is None:
            return page_num
        msg = ['                           O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R']
        msg.append('0                                                     REFERENCE POINT =        %i' % self.reference_point)

        # MO
        msg.append('                                                                M O')
        for i in range(6):
            msg.append('                      * %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E *' % tuple(self.MO[i, :]))

        msg.append('                                                                 S')
        for i in range(3):
            msg.append('                                           * %13.6E %13.6E %13.6E *' % tuple(self.S[i, :]))

        msg.append('                               DIRECTION')
        msg.append('                          MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.')
        msg.append('                                  X            %12.6E     %13.6E %13.6E %13.6E' % (self.mass[0], self.cg[0, 0], self.cg[0, 1], self.cg[0, 2]))
        msg.append('                                  Y            %12.6E     %13.6E %13.6E %13.6E' % (self.mass[1], self.cg[1, 0], self.cg[1, 1], self.cg[1, 2]))
        msg.append('                                  Z            %12.6E     %13.6E %13.6E %13.6E' % (self.mass[2], self.cg[2, 0], self.cg[2, 1], self.cg[2, 2]))

        msg.append('                                                                I(S)')
        for i in range(3):
            msg.append('                                           * %13.6E %13.6E %13.6E *' % tuple(self.IS[i, :]))

        msg.append('                                                                I(Q)')
        msg.append('                                           * %13.6E %13s %13s *' % (self.IQ[0], '', ''))
        msg.append('                                           * %13s %13.6E %13s *' % ('', self.IQ[1], ''))
        msg.append('                                           * %13s %13s %13.6E *' % ('', '', self.IQ[2]))


        msg.append('                                                                 Q')
        for i in range(3):
            msg.append('                                           * %13.6E %13.6E %13.6E *' % tuple(self.Q[i, :]))
        msg.append('\n' + page_stamp % page_num + '\n')
        f06_file.write('\n'.join(msg))
        #print('\n'.join(msg))
        return page_num + 1
