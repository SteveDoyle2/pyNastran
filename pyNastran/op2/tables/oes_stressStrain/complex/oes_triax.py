import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import OES_Object
#from pyNastran.f06.f06_formatting import write_imag_floats_13e, write_float_13e


class ComplexTriaxArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 4 # CTRIAX6 = centroid + 3

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def build(self) -> None:
        """sizes the vectorized attributes of the ComplexPlateArray"""
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        nnodes = self.nnodes_per_element

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #print('****nelements =', self.nelements, self.ntotal, self.element_name)
        #print('element_type=%r ntimes=%s nelements=%s nnodes=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, nnodes, self.ntotal, self.subtitle))

        nlayers = 1
        self.ntotal = self.nelements * nnodes# * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self.eids = zeros(self.ntotal, dtype=idtype)
        self.element_node = zeros((self.ntotal, 2), idtype)
        #self.element_cid = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes * nlayers == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes, self.ntotal)
            raise RuntimeError(msg)

        self.fiber_curvature = zeros(self.ntotal, 'float32')
        # [e_radial, e_azimuthal, e_axial, e_shear]
        self.data = zeros((self.ntimes, self.ntotal, 4), 'complex64')

    def get_headers(self) -> list[str]:
        return self._get_headers()

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nnodes, self.table_name))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 3] where 3=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid, Nid\n' % str(self.code_information())
            for (eid1, nid1), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s), (%s, %s)\n' % (eid1, nid1, eid2, nid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, (eid, nid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (oxx1, oyy1, txy1) = t1
                        (oxx2, oyy2, txy2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                                [oxx1.real, oxx1.imag, oyy1.real, oyy1.imag, txy1.real, txy1.imag, ], # atol=0.0001
                                [oxx2.real, oxx2.imag, oyy2.real, oyy2.imag, txy2.real, txy2.imag, ], atol=0.075):
                            ni = len(str(eid)) + len(str(nid))
                        #if not np.array_equal(t1, t2):
                            msg += ('(%s %s)  (%s, %sj, %s, %sj, %s, %sj)\n'
                                    '%s     (%s, %sj, %s, %sj, %s, %sj)\n' % (
                                        eid, nid,
                                        oxx1.real, oxx1.imag, oyy1.real, oyy1.imag,
                                        txy1.real, txy1.imag,
                                        ' ' * ni,
                                        oxx2.real, oxx2.imag, oyy2.real, oyy2.imag,
                                        txy2.real, txy2.imag,
                                    ))
                            msg += ('%s     (%s, %sj, %s, %sj, %s, %sj)\n'
                                    % (
                                        ' ' * ni,
                                        oxx1.real - oxx2.real, oxx1.imag - oxx2.imag,
                                        oyy1.real - oyy2.real, oyy1.imag - oyy2.imag,
                                        txy1.real - txy2.real, txy1.imag - txy2.imag,
                                    ))

                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_element_sort1(self, dt, eid):
        self._times[self.itime] = dt
        self.eids[self.ielement] = eid
        self.ielement += 1
        #print('eid =', eid)

    def add_sort1(self, dt, eid, loc, rs, azs, As, ss):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s loc=%s R=%s azimuth=%s axial=%s shear=%s' % (
            #self.itotal, dt, eid, loc,
            #rs, azs, As, ss))

        # dt, eid, loc, rs, azs, As, ss
        assert isinstance(eid, int), eid
        self.data[self.itime, self.itotal] = [rs, azs, As, ss]
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        #self.fiber_curvature[self.itotal] = fdr
        #self.ielement += 1
        self.itotal += 1

class ComplexTriaxStressArray(ComplexTriaxArray):
    def _get_headers(self) -> list[str]:
        return ['o_radial', 'o_azimuthal', 'o_axial', 'o_shear']

class ComplexTriaxStrainArray(ComplexTriaxArray):
    def _get_headers(self) -> list[str]:
        return ['e_radial', 'e_azimuthal', 'e_axial', 'e_shear']


'      COMPLEX EIGENVALUE = -8.256059E-03,  5.505277E-01'
'                              C O M P L E X    S T R A I N S    I N   T R I A X 6   E L E M E N T S'
'                                                            (REAL/IMAGINARY)'
'   ELEMENT'
'       ID   GRID ID         R A D I A L              A Z I M U T H A L               A X I A L                   S H E A R'
'      5301        0  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00'
'               5301  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00'
'               5303  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00'
'               5305  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00'
''
'      5311        0  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00  0.000000E+00/ 0.000000E+00'
