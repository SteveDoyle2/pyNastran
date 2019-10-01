"""defines the GridPointWeight class"""
from io import StringIO

from pyNastran.utils import object_attributes, object_methods


class GridPointWeight:
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
        return page_num + 1
