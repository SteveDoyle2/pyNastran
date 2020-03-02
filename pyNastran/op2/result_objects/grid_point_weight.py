"""defines the GridPointWeight class"""
from io import StringIO
from struct import pack
import numpy as np

from pyNastran.utils import object_attributes, object_methods
from pyNastran.op2.result_objects.op2_objects import _write_table_header
from pyNastran.op2.op2_interface.write_utils import export_to_hdf5

float_types = (float, np.float32)
integer_types = (int, np.int32)

#                                                                       ?           ?       ?                                     ?
#good = (4, 2, 4, 8, 1464878927, 538976327, 8, 4, -1, 4, 4, 7, 4, 28, 101, 0, 0, 0, 0,   0, 1, 28, 4, -2, 4, 4, 1, 4, 4, 0, 4, 4, 2, 4, 8,  1464878927, 538976327, 8, 4,                   -3, 4, 4, 1, 4, 4, 0, 4, 4, 146, 4)
#bad  = (4, 2, 4, 8, 1464878927, 538976327, 8, 4, -1, 4, 4, 7, 4, 28, 102, 0, 0, 0, 512, 0, 0, 28, 4, -2, 4, 4, 1, 4, 4, 0, 4, 4, 7, 4, 28, 1464878927, 538976327, 9, 27, 19, 0, 1, 28, 4, -3, 4, 4, 1, 4, 4)
class GridPointWeight:
    def __init__(self, reference_point, MO, S, mass, cg, IS, IQ, Q,
                 approach_code=1, table_code=13,
                 title='', subtitle='', label='',
                 superelement_adaptivity_index=''):
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
        self.reference_point = reference_point
        assert isinstance(reference_point, int), f'reference_point={reference_point!r}; type={type(reference_point)}'

        # M0 RIGID BODY MASS MATRIX IN BASIC COORDINATE SYSTEM
        # This is the rigid body mass matrix of the entire structure in
        # the basic coordinate system with respect to a reference point
        # chosen by the analyst.
        self.MO = MO
        assert MO.shape == (6, 6), MO.shape

        # S TRANSFORMATION MATRIX FOR SCALAR MASS PARTITION
        # S is the transformation from the basic coordinate system to the
        # set of principal axes for the 3 x 3 scalar mass partition of the
        # 6 x 6 mass matrix. The principal axes for just the scalar
        # partition are known as the principal mass axes.
        self.S = S
        assert S.shape == (3, 3), S.shape

        self.mass = mass
        assert mass.shape == (3,), mass.shape

        # XC.G. YC.G. ZC.G.
        # It is possible in NASTRAN to assemble a structural model having
        # different values of mass in each coordinate direction at a grid
        # point. This can arise, for example, by assembling scalar mass
        # components or from omitting some components by means of bar
        # element pin flags. Consequently three distinct mass systems are
        # assembled one in each of the three directions of the principal
        # mass axes (the S system). This third tabulation has five columns.
        # The first column lists the axis direction in the S coordinates.
        # The second column lists the mass associated with the appropriate
        # axis direction. The final three columns list the x, y, and z
        # coordinate distances from the reference point to the center of
        # mass for each of the three mass systems.
        self.cg = cg
        assert cg.shape == (3, 3), cg.shape

        # I(S) INERTIAS RELATIVE TO C.G.
        # This is the 3 x 3 mass moment of inertia partition with respect
        # to the center of gravity referred to the principal mass axes
        # (the S system).
        #
        # This is not necessarily a diagonal matrix because the
        # determination of the S system does not involve second moments.
        # The values of inertias at the center of gravity are found from
        # the values at the reference point employing the
        # parallel axes rule.
        self.IS = IS
        assert IS.shape == (3, 3), IS.shape

        # I(Q) PRINCIPAL INERTIAS
        # The principal moments of inertia at the center of gravity are displayed
        # in matrix form with reference to the Q system of axes. The Q system is
        # obtained from an eigenvalue analysis of the I(S) matrix.
        self.IQ = IQ
        assert IQ.shape == (3, ), f'IQ.shape={IQ.shape}; IQ=\n{IQ}'

        # Q TRANSFORMATION MATRIX I(Q) = Q^T*IBAR(S)*Q
        # Q is the coordinate transformation between the S axes and the Q axes.
        # IBAR(S) is the same as I(s) except that the signs of the off-diagonal
        # terms are reversed.
        self.Q = Q
        assert Q.shape == (3, 3), Q.shape

        self.title = title
        self.subtitle = subtitle
        self.label = label
        self.superelement_adaptivity_index = subtitle
        self.approach_code = approach_code
        self.table_code = table_code

    def export_to_hdf5(self, group, log) -> None:
        """exports the object to HDF5 format"""
        export_to_hdf5(self, group, log)

    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode,
                                 keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def __eq__(self, weight):
        msg = ''
        if not self.reference_point == weight.reference_point:
            msg += f'reference_point: {self.reference_point} -> {weight.reference_point}\n'
        if not np.array_equal(self.MO, weight.MO):
            msg += f'reference_point: {self.MO} -> {weight.MO}\n'
        if not np.array_equal(self.S, weight.S):
            msg += f'reference_point: {self.S} -> {weight.S}\n'
        if not np.array_equal(self.mass, weight.mass):
            msg += f'reference_point: {self.mass} -> {weight.mass}\n'
        if not np.array_equal(self.cg, weight.cg):
            msg += f'reference_point: {self.cg} -> {weight.cg}\n'
        if not np.array_equal(self.IS, weight.IS):
            msg += f'reference_point: {self.IS} -> {weight.IS}\n'
        if not np.array_equal(self.IQ, weight.IQ):
            msg += f'reference_point: {self.IQ} -> {weight.IQ}\n'
        if not np.array_equal(self.Q, weight.Q):
            msg += f'reference_point: {self.Q} -> {weight.Q}'
        if msg:
            raise ValueError('GridPointWeight:\n' + msg)
        return True

    def get_stats(self, key='', short=True):
        key2 = f'[{key!r}]'
        if short:
            msg = (f'GridPointWeight{key2}: ref_point=%s mass=%g; '
                   '[reference_point, M0, S, mass, cg, IS, IQ, Q]\n' % (
                       self.reference_point, self.mass.max()))
        else:
            msg = (
                f'GridPointWeight{key2}:'
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

    def _write_table_3(self, op2_file, fascii, new_result, table_name, itable=-3):
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2_file.write(pack(b'%ii' % len(header), *header))
        fascii.write('table_3_header = %s\n' % header)
        #op2_file.write(pack('12i', *[4, itable, 4,
                              #4, 1, 4,
                              #4, 0, 4,
                              #4, 146, 4,
                              #]))

        approach_code = self.approach_code
        table_code = self.table_code
        #isubcase = self.isubcase
        #random_code = self.random_code
        #format_code = 1
        isubcase = 0
        num_wide = 79 # self.num_wide
        #acoustic_flag = self.acoustic_flag if hasattr(self, 'acoustic_flag') else 0
        reference_point = self.reference_point
        reference_point = 22
        #thermal = self.thermal
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')  # missing superelement_adaptivity_index
        label = b'%-128s' % self.label.encode('ascii')

        #1, 13, 0, 0, 0, 0, 0, 0, 0, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ftable3 = b'i' * 50 + b'128s 128s 128s'

        #print(self.get_stats())
        table3 = [
            approach_code, table_code, reference_point, isubcase, 0,
            0, 0, 0, 0, num_wide,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]

        n = 0
        from itertools import count
        for i, val, ftable3i in zip(count(), table3, ftable3.decode('ascii')):
            assert val is not None, 'i=%s val=%s ftable3i=%s\n%s' % (i, val, ftable3i, self.get_stats())
            if isinstance(val, integer_types):
                n += 4
                assert ftable3i == 'i', 'i=%s val=%s type=%s' % (i, val, ftable3i)
            elif isinstance(val, float_types):
                n += 4
                assert ftable3i == 'f', 'i=%s val=%s type=%s' % (i, val, ftable3i)
            else:
                n += len(val)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'

        #op2_file.write(pack(fascii, '%s header 3c' % table_name, fmt, data))
        fascii.write('%s header 3c = %s\n' % (table_name, data))

        #j = 7
        #print(ftable3[:j])
        #print(table3[:j])
        #pack(ftable3[:j], *table3[:j])
        op2_file.write(pack(fmt, *data))

    def write_op2(self, op2_file, op2_ascii, date, endian=b'<'):
        itable = -1
        import inspect
        #allowed_tables = [
            #'OUGV1', 'BOUGV1', 'BOPHIG', 'BOPG1',
            #'OUPV1',
            #'OQP1', 'OQMG1', 'OQG1', 'OQGV1', 'OPNL1',
            #'OPG1', 'OPGV1',
            #'OAGATO1', 'OAGCRM1', 'OAGNO1', 'OAGPSD1', 'OAGRMS1',
            #'OQGPSD1',
            #'OCRPG', 'OCRUG', 'OUG1',
            #'OUGV1PAT',
        #]
        #assert self.table_name in allowed_tables, self.table_name
        table_name = 'OGPWG'

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        subtable_name = b'OGPWG'
        if itable == -1:
            _write_table_header(op2_file, op2_ascii, date, table_name, subtable_name)
            itable = -3

        #s = Struct(op2_format)

        #unused_node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]
        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        #nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        #ntotal = nnodes * (2 + 6)

        #print('shape = %s' % str(self.data.shape))
        #assert nnodes > 1, nnodes
        #assert ntotal > 1, ntotal

        #unused_device_code = self.device_code
        #fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #for itime in range(self.ntimes):

        # good = (4, -4, 4, 4, 1, 4, 4, 0, 4, 4, 78, 4, 312,    1057795080, 0, 0, 0, 0, 0, 0, 1057795080, 0, 0, 0, 1111254630, 0, 0, 1057795080, 0, -1036229018, 0, 0, 0, 0, 1143715840, -1451229184, 0, 0, 0, -1036229018, -1451229184, 1169886464, 0, 0, 1111254630, 0, 0, 0, 1171293184, 1065353216)
        #bad   = (4, -4, 4, 4, 1, 4, 4, 0, 4, 4, 78, 4, 312,    1057795080, 0, 0, 0, 0, 0, 0, 1057795080, 0, 0, 0, 1111254630, 0, 0, 1057795080, 0, -1036229018, 0, 0, 0, 0, 1143715840, -1451229184, 0, 0, 0, -1036229018, -1451229184, 1169886464, 0, 0, 1111254630, 0, 0, 0, 1171293184, 1065353216, 0, 0, 0, 1065353216, 0, 0, 0, 1065353216, 1057795080, 0, -2147483648, 0, 1057795080, 1118530999, 0, -2147483648, 1057795080, 1118530999, 0, 0, 1143715840, 696254464, 0, 696254464, 1156812654, 0, 0, 0, 1160033719, 1143715840, 1156812654, 1160033719, 1065353216, 0, 0, 0, 1065353216, 0, 0, 0, 1065353216, 312, -5, 4, 4, 1, 4, 4, 0, 4, 4, 0, 4, 4, 4, 2, 4, 8, 1447515471, 538976305)
        ntotal = 78  # 79 * 4 = 316
        new_result = True
        self._write_table_3(op2_file, op2_ascii, new_result, table_name, itable)

        # record 4
        itable -= 1
        header = [4, itable, 4,
                  4, 1, 4,
                  4, 0, 4,
                  4, ntotal, 4,
                  4*ntotal]
        op2_file.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('r4 [4, 0, 4]\n')
        op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
        op2_ascii.write('r4 [4, %i, 4]\n' % (4*ntotal))

        # -------------------------------------------------------
        fmt = endian + b'78f'
        mcg = np.zeros((3, 4), dtype=self.cg.dtype)
        mcg[:, 0] = self.mass
        mcg[:, 1:] = self.cg
        data = (self.MO.ravel().tolist() + self.S.ravel().tolist() +
                mcg.ravel().tolist() + self.IS.ravel().tolist() + self.IQ.ravel().tolist() +
                self.Q.ravel().tolist())
        assert None not in data, data
        msgi = pack(fmt, *data)
        op2_file.write(msgi)
        # -------------------------------------------------------
        itable -= 1
        header = [4 * ntotal,
                  4, itable, 4,
                  4, 1, 4,
                  4, 0, 4,
                  4, 0, 4, ]
        op2_file.write(pack(endian + b'13i', *header))
        op2_ascii.write('footer = %s\n' % header)
        return itable


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

def make_grid_point_weight(reference_point, MO,
                           approach_code=1, table_code=13,
                           title='', subtitle='', label='',
                           superelement_adaptivity_index='') -> None:
    """creates a grid point weight table"""
    Mtt_ = MO[:3, :3]
    Mrr_ = MO[3:, 3:]
    Mtr_ = MO[:3, 3:]
    Mtd = np.diag(Mtt_)
    delta = np.linalg.norm(Mtd)
    e_ = [Mtt_[0, 1], Mtt_[0, 2], Mtt_[1, 2]]
    epsilon = np.linalg.norm(e_)
    #print(Mtd)
    #print(Mte)
    #print(Mtd)

    #print(Mte)
    if epsilon/delta > 0.001:
        unused_eigvals, S = np.linalg.eigh(Mtt_)
        #print('S1:')
        #print(S)
        #print(f'eigvals = {eigvals}')
        msg = (
            '*** USER WARNING MESSAGE 3042 MODULE = GPWG\n'
            f'INCONSISTENT SCALAR MASSES HAVE BEEN USED. EPSILON/DELTA = {epsilon/delta:.7E}\n')
        print(msg)
        #print('S*:')
        #print(S)
    else:
        S = np.eye(3, dtype=Mtt_.dtype)


    #iswap = []
    #if np.sign(S[0, 0]) != np.sign(Mtt_[0, 0]):
        #mass[0] *= -1.
        #S[:, 0] *= -1
        #iswap.append(0)

    #if np.sign(S[1, 1]) != np.sign(Mtt_[1, 1]):
        #iswap.append(1)
    #if S[1, 1] < 0.:
        #iswap.append(1)
    #if S[2, 2] < 0.:
        #iswap.append(2)

    #if abs(S[0, 0]) < np.abs(S[:, 0]).max():
        #print('swap0')
        #iswap.append(0)
    #if S[0, 0] < 0.:
        #S[:, 0] *= -1
    #if abs(S[1, 1]) < np.abs(S[1, 0]).max():
        #iswap.append(1)
    #if abs(S[0, 0]) < np.abs(S[1, 0]):
        #iswap.append(2)

    #if iswap:
        #i1, i2 = iswap
        #print(f'swap; {iswap}')
        #m1, m2 = mass[iswap]
        #mass[[i2, i1]] = mass[iswap]
        #S[[i2, i1], :] = S[iswap, :]
        #S[:, [i2, i1]] = S[:, iswap]
    #print('S2:')
    #print(S)
    #print('S3:')
    #S = np.array([
        #[.432, 0.902, 0.],
        #[-.902, .432, 0.],
        #[0., 0., 1.],
    #])
    #print(S)

    Mtt = S.T @ Mtt_ @ S # Mt
    Mtr = S.T @ Mtr_ @ S
    Mrr = S.T @ Mrr_ @ S # Mr = I(Q)?
    #print('---------')
    #print(Mtt)

    #cg = Mtr / mass[:, np.newaxis]
    mx, my, mz = Mtt[0, 0], Mtt[1, 1], Mtt[2, 2]
    mass = np.hstack([mx, my, mz])
    cg = np.vstack([
        [Mtr[0, 0], -Mtr[0, 2], Mtr[0, 1]],  # Xx, Yx, Zx
        [Mtr[1, 2], Mtr[1, 1], -Mtr[1, 0]],  # Xy, Yy, Zy
        [-Mtr[2, 1], Mtr[2, 0], -Mtr[2, 2]], # Xz, Yz, Zz
    ])
    if abs(mx) > 0:
        cg[0, :] /= mx
    if abs(my) > 0:
        cg[1, :] /= my
    if abs(mz) > 0:
        cg[2, :] /= mz
    Xx, Yx, Zx = cg[0, :]
    Xy, Yy, Zy = cg[1, :]
    Xz, Yz, Zz = cg[2, :]
    IS = np.zeros((3, 3), dtype=MO.dtype)
    IS[0, 0] = Mrr[0, 0] - my * Zy ** 2 - mz * Yz ** 2
    IS[1, 1] = Mrr[1, 1] - mz * Xz ** 2 - mx * Zx ** 2
    IS[2, 2] = Mrr[2, 2] - mx * Yx ** 2 - my * Xy ** 2

    IS[1, 0] = IS[0, 1] = -Mrr[0, 1] - mz * Xz * Yz
    IS[2, 0] = IS[0, 2] = -Mrr[0, 2] - my * Xy * Zy
    IS[2, 1] = IS[1, 2] = -Mrr[1, 2] - mx * Yx * Zx

    neg_off_diag = np.array([
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]], dtype=MO.dtype)
    IQi = neg_off_diag * IS
    #print(IQi)
    IQ, Q = np.linalg.eigh(IQi)
    #IQ = np.diag(Q.T @ IS @ Q)

    weight = GridPointWeight(
        reference_point, MO, S, mass, cg, IS, IQ, Q,
        approach_code=approach_code, table_code=table_code,
        title=title, subtitle=subtitle, label=label,
        superelement_adaptivity_index=superelement_adaptivity_index)
    return weight
