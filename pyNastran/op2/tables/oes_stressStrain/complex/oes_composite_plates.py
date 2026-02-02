# import warnings
from itertools import zip_longest
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header, build_pandas_transient_element_node
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import (
    write_float_12e, write_imag_floats_12e,
    # write_float_13e,
    write_imag_floats_13e,
)

#BASIC_TABLES = {
    #'OES1X', 'OES1',
    #'OES2',
    #'OSTR1X',
#}
#VM_TABLES = {'OESVM1', 'OESVM2',
             #'OSTRVM1', 'OSTRVM2'}

class ComplexLayeredCompositesVMArray(OES_Object):
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
        assert self.table_name is not None, self.data_code

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def nnodes_per_element(self) -> int:
        return 1

    #@property
    #def nnodes(self):
        #return self.nnodes_per_element()

    def build(self) -> None:
        r"""sizes the vectorized attributes of the ComplexPlateArray

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqrdbx111.op2
        name;      nelements  numwide ndata size ntotal     nelements     nnodes nlayers

        """
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        # print(self._ntotals, self.ntotal)
        #print(self.code_information())

        # self.names = []
        # self.nelements //= nnodes
        self.nelements //= self.ntimes
        # print('element_type=%r ntimes=%s nelements=%s nnodes=%s ntotal=%s subtitle=%s' % (
        #     self.element_type, self.ntimes, self.nelements, nnodes, self.ntotal, self.subtitle))

        # self.ntotal = self.nelements * nnodes * 2
        # self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        # print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        idtype, cfdtype = get_complex_times_dtype(self.size)

        if self.is_sort1:
            ntimes = self.ntimes
            nlayers = self.ntotal
            # print(f'  SORT1: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        else:
            raise NotImplementedError(self.code_information())
        # elif self.is_sort2:
        #     nelements = self.ntimes
        #     nlayers = nelements * 2 * nnodes
        #     ntimes = self.ntotal
        #     print(f'  SORT2: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        # print("nelements=%s nlayers=%s ntimes=%s" % (nelements, nlayers, ntimes))

        self._times = np.zeros(ntimes, dtype=self.analysis_fmt)
        # self.ntotal = self.nelements * nnodes

        # the number is messed up because of the offset for the element's properties
        # if not self.nelements * nnodes * 2 == self.ntotal:
        #     msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
        #         self.ntimes, self.nelements, nnodes,
        #         self.nelements * nnodes, self.ntotal)
        #     raise RuntimeError(msg)

        # [o1, o2, t12, o1z, o2z, ovm]
        self.data = np.zeros((ntimes, nlayers, 6), dtype=cfdtype)

        self.element_layer = np.zeros((nlayers, 2), dtype=idtype)
        # print(self.data.shape, self.element_node.shape)


    def build_dataframe(self) -> None:
        """creates a pandas dataframe"""
        self.data_frame = build_dataframe(self)

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        assert table.element_layer.shape == self.element_layer.shape
        assert table.data.shape == self.data.shape

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_layer in enumerate(self.element_layer):
                    (eid, layer) = eid_layer
                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, t121, t131, t231, ovm1) = t1
                    (oxx2, oyy2, t122, t132, t232, ovm2) = t2
                    if not np.allclose(t1, t2):
                        raise RuntimeError((t1, t2))
                        # msg += (
                        #     # (oxx1, oyy1, t121, t131, t231, ovm1) = t1
                        # '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s)\n'
                        #     '%s      (%s, %s, %s, %s, %s, %s, %s)\n' % (
                        #         eid, nid,
                        #         oxx1, oyy1, ozz1, txy1, tyz1, txz1, ovm1,
                        #         ' ' * (len(str(eid)) + len(str(nid)) + 2),
                        #         oxx2, oyy2, ozz2, txy2, tyz2, txz2, ovm2))
                        # i += 1
                        # if i > 10:
                        #     print(msg)
                        #     raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, ply_id,
                  o1, o2, t12, o1z, o2z, ovm):
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(ply_id, int), ply_id
        self.data[self.itime, self.itotal] = [o1, o2, t12, o1z, o2z, ovm]
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nlayers = self.element_layer.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nlayers=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nlayers, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nlayers=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nlayers, self.table_name))
        msg.append('  eType, cid\n')
        headers = self._get_headers()
        nheaders = len(headers)
        headers_str = ', '.join(headers)
        msg.append(f'  data: [ntimes, nlayers, {nheaders}] where {nheaders}=[{headers_str}]\n')
        msg.append(f'  element_layer.shape = {self.element_layer.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True) -> int:
        if header is None:
            header = []
        msg_temp, nnodes = _get_composite_plate_msg(self, is_mag_phase, is_sort1)
        # if not self.is_von_mises:
        #     warnings.warn(f'{self.class_name} doesnt support writing von Mises')
        #     f06_file.write(f'{self.class_name} doesnt support writing von Mises\n')

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

            # '      FREQUENCY =  3.000000E+01'
            # '                C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'
            # '                                                         (MAGNITUDE/PHASE)                          OPTION = BILIN  '
            # ' '
            # '    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -'
            # '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY'
            # '0         6    CEN/4   0.0             0.0          /   0.0            0.0          /   0.0            0.0          /   0.0'
            # '                      -1.000000E+00    0.0          /   0.0            0.0          /   0.0            0.0          /   0.0'
            # ''

            #if self.element_type == 1: # 144-CQUAD4 bilinear
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, is_mag_phase, 'CEN/4')
            #elif self.element_type in [33, 74, 227, 228]:
                ## CQUAD4 linear, CTRIA3, CTRIAR linear, CQUADR linear
                #self._write_f06_tri3_transient(f06_file, itime, is_mag_phase)
            #elif self.element_type == 64:  #CQUAD8
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 5, is_mag_phase, 'CEN/8')
            #elif self.element_type == 82:  # CQUADR
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 5, is_mag_phase, 'CEN/8')
            #elif self.element_type == 70:  # CTRIAR
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 3, is_mag_phase, 'CEN/3')
            #elif self.element_type == 75:  # CTRIA6
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 3, is_mag_phase, 'CEN/6')
            #else:
                #raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            # o11, o22, t12, t1z, t2z, ovm
            o11 = self.data[itime, :, 0]
            o22 = self.data[itime, :, 1]
            t12 = self.data[itime, :, 2]
            t1z = self.data[itime, :, 3]
            t2z = self.data[itime, :, 4]

            eids = self.element_layer[:, 0]
            layers = self.element_layer[:, 1]
            ovm = self.data[itime, :, 5].real
            for (eid, layer,
                 do11, do22, dt12, dt1z, dt2z, dovm) in zip(eids, layers,
                                          o11, o22, t12, t1z, t2z, ovm):
                ovmr = write_float_12e(dovm)
                [o11r, o22r, t12r, t1zr, t2zr,
                 o11i, o22i, t12i, t1zi, t2zi,] = write_imag_floats_12e([
                    do11, do22, dt12, dt1z, dt2z,
                ], is_mag_phase)
                # """
                #    ELEMENT      PLY STRESSES IN FIBER AND MATRIX DIRECTIONS   INTER-LAMINAR  STRESSES
                #      ID          ID   NORMAL-1     NORMAL-2     SHEAR-12    SHEAR XZ-MAT  SHEAR YZ-MAT  VON MISES
                #       1014        1  8.13423E+03  1.75179E+02  1.50284E+03  -1.79942E+00  3.44091E+00
                #                     -1.47375E+02 -5.35187E+00 -2.76251E+01   2.44372E-01 -1.24947E-01  8.45992E+03
                # """
                f06_file.write(f'  {eid:8d} {layer:8d} {o11r:<12s} {o22r:<12s} {t12r:<12s} {t1zr:<12s} {t2zr:<12s}\n'
                                 f'                    {o11i:<12s} {o22i:<12s} {t12i:<12s} {t1zi:<12s} {t2zi:<12s} {ovmr:s}\n')
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

def _get_composite_plate_msg(self, is_mag_phase: bool=True,
                             is_sort1: bool=True) -> tuple[list[str], int]:
    if self.is_von_mises:
        von = 'VON'
        mises = 'MISES'
    else:
        von = 'MAX'
        mises = 'SHEAR'

    if self.is_strain:
        words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % von,
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]
    else:
        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % von,
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

    if self.element_type == 95:  # CQUAD4
        nnodes = 4
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
    #elif self.element_type == 96:  # CQUAD8
        #nnodes_per_element = 1
    elif self.element_type == 97:  # CTRIA3
        nnodes = 3
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
    elif self.element_type == 96:  # QUAD8
        nnodes = 8
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words

    elif self.element_type == 98:  # CTRIA6
        nnodes = 6
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
    elif self.element_type == 232:  # CQUADR
        nnodes = 4
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )\n'] + words
    elif self.element_type == 233:  # CTRIAR
        nnodes = 6
        if self.is_strain:
            msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )\n'] + words
        else:
            msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )\n'] + words
    else:  # pragma: no cover
        msg = 'element_name=%s element_type=%s' % (self.element_name, self.element_type)
        raise NotImplementedError(msg)
    return msg, nnodes


class ComplexLayeredCompositesArray(OES_Object):
    """MSC"""
    def __init__(self, data_code, is_sort1: bool, isubcase: int, dt):
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
        assert self.table_name is not None, self.data_code

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def nnodes_per_element(self) -> int:
        return 1

    #@property
    #def nnodes(self):
        #return self.nnodes_per_element()

    def build(self) -> None:
        r"""sizes the vectorized attributes of the ComplexPlateArray12

        C:\MSC.Software\msc_nastran_runs\1-elm--pcomp-mat2-pdistb.op2

        """
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        # print(self._ntotals, self.ntotal)
        #print(self.code_information())

        # self.names = []
        # self.nelements //= nnodes
        self.nelements //= self.ntimes
        # print('element_type=%r ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
        #     self.element_type, self.ntimes, self.nelements, self.ntotal, self.subtitle))

        # self.ntotal = self.nelements * nnodes * 2
        # self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        # print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        idtype, cfdtype = get_complex_times_dtype(self.size)
        if self.is_sort1:
            ntimes = self.ntimes
            nlayers = self.ntotal
            # print(f'  SORT1: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        else:
            raise NotImplementedError(self.code_information())
        # elif self.is_sort2:
        #     nelements = self.ntimes
        #     nlayers = nelements * 2 * nnodes
        #     ntimes = self.ntotal
        #     print(f'  SORT2: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        # print("nelements=%s nlayers=%s ntimes=%s" % (nelements, nlayers, ntimes))

        self._times = np.zeros(ntimes, dtype=self.analysis_fmt)
        # self.ntotal = self.nelements * nnodes

        # the number is messed up because of the offset for the element's properties
        # if not self.nelements * nnodes * 2 == self.ntotal:
        #     msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
        #         self.ntimes, self.nelements, nnodes,
        #         self.nelements * nnodes, self.ntotal)
        #     raise RuntimeError(msg)

        # ELEMENT   PLY   STRESSES IN FIBER AND MATRIX DIRECTIONS     INTER-LAMINAR  STRESSES
        #      ID    ID       NORMAL-1      NORMAL-2      SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT
        #       1     1   -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05
        #                 -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05
        # [oxx, oyy, t12, t1z, t2z]
        self.data = np.zeros((ntimes, nlayers, 5), dtype=cfdtype)

        self.element_layer = np.zeros((nlayers, 2), dtype=idtype)
        # print(self.data.shape, self.element_node.shape)

    def build_dataframe(self) -> None:
        """creates a pandas dataframe"""
        self.data_frame = build_dataframe(self)

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        assert table.element_layer.shape == self.element_layer.shape
        assert table.data.shape == self.data.shape

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_layer in enumerate(self.element_layer):
                    (eid, layer) = eid_layer
                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, t121, t131, t231, ovm1) = t1
                    (oxx2, oyy2, t122, t132, t232, ovm2) = t2
                    if not np.array_equal(t1, t2):
                        raise RuntimeError((t1, t2))
        return True

    def add_sort1(self, dt, eid, ply_id,
                  oxx, oyy, t12, t1z, t2z):
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        # print(self.element_types2, element_type, self.element_types2.dtype)
        # print('itotal=%s dt=%s eid=%s ply_id=%s oxx=%s' % (self.itotal, dt, eid, ply_id, oxx))
        assert isinstance(ply_id, int), ply_id
        self.data[self.itime, self.itotal] = [oxx, oyy, t12, t1z, t2z]
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]
        nelements = self.nelements
        ntimes = self.ntimes
        nlayers = self.element_layer.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nlayers=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nlayers, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nlayers=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nlayers, self.table_name))
        msg.append('  eType, cid\n')
        headers = self._get_headers()
        nheaders = len(headers)
        headers_str = ', '.join(headers)
        msg.append(f'  data: [ntimes, nlayers, {nheaders}] where {nheaders}=[{headers_str}]\n')
        msg.append(f'  element_layer.shape = {self.element_layer.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True) -> int:
        """
              FREQUENCY =  2.700000E+01
                         C O M P L E X   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
                                                                  (REAL/IMAGINARY)
           ELEMENT              PLY   STRESSES IN FIBER AND MATRIX DIRECTIONS     INTER-LAMINAR  STRESSES
                ID              ID       NORMAL-1      NORMAL-2      SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT
                 1               1   -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05
                                      1.307266E-01  1.698439E+00  3.178526E-03    3.383883E-02 -8.803315E-06
        """
        element_type = self.element_type
        if element_type == 95:
            element_name = 'Q U A D 4 '
        elif element_type == 97:
            element_name = 'T R I A 3 '
        elif element_type == 98:
            element_name = 'T R I A 6 '
        elif element_type == 232:
            element_name = 'Q U A D R '
        elif element_type == 233:
            element_name = 'T R I A R '
        else:  # pragma: no cover
            raise RuntimeError(self.code_information())

        if self.is_stress:
            msg_temp = [
                f'                 C O M P L E X   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ({element_name})\n'
                '                                                          (REAL/IMAGINARY)\n'
                '   ELEMENT              PLY   STRESSES IN FIBER AND MATRIX DIRECTIONS     INTER-LAMINAR  STRESSES\n'
                '        ID              ID       NORMAL-1      NORMAL-2      SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT\n'
                # '         1               1   -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
                # '                             -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
            ]
        else:
            msg_temp = [
                f'                 C O M P L E X   S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ({element_name})\n'
                '                                                          (REAL/IMAGINARY)\n'
                '   ELEMENT              PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS      INTER-LAMINAR  STRAINS\n'
                '        ID              ID       NORMAL-1      NORMAL-2      SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT\n'
                # '         1               1   -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
                # '                             -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
            ]
        if header is None:
            header = []
        # msg_temp, nnodes = _get_composite_plate_msg(self, is_mag_phase, is_sort1)
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

            # '      FREQUENCY =  3.000000E+01'
            # '                C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'
            # '                                                         (MAGNITUDE/PHASE)                          OPTION = BILIN  '
            # ' '
            # '    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -'
            # '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY'
            # '0         6    CEN/4   0.0             0.0          /   0.0            0.0          /   0.0            0.0          /   0.0'
            # '                      -1.000000E+00    0.0          /   0.0            0.0          /   0.0            0.0          /   0.0'
            eids = self.element_layer[:, 0]
            layers = self.element_layer[:, 1]
            oxxs = self.data[itime, :, 0]
            oyys = self.data[itime, :, 1]
            txys = self.data[itime, :, 2]
            t1zs = self.data[itime, :, 3]
            t2zs = self.data[itime, :, 4]
            for (eid, layr, oxx, oyy, txy, t1z, t2z) in zip_longest(eids, layers, oxxs, oyys, txys, t1zs, t2zs):
                [oxxr, oyyr, txyr, t1zr, t2zr,
                 oxxi, oyyi, txyi, t1zi, t2zi] = write_imag_floats_13e(
                     [oxx, oyy, txy, t1z, t2z], is_mag_phase)
                f06_file.write(
                    f'  {eid:8d}        {layr:8d}   {oxxr} {oyyr} {txyr}   {t1zr} {t2zr}\n'
                     f'                             {oxxi} {oyyi} {txyi}   {t1zi} {t2zi}\n')

            # '         1               1   -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
            # '                             -8.713408E-01 -1.132072E+01 -2.119228E-02   -2.255483E-01  5.870259E-05'
            #if self.element_type == 1: # 144-CQUAD4 bilinear
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, is_mag_phase, 'CEN/4')
            #elif self.element_type in [33, 74, 227, 228]:
                ## CQUAD4 linear, CTRIA3, CTRIAR linear, CQUADR linear
                #self._write_f06_tri3_transient(f06_file, itime, is_mag_phase)
            #elif self.element_type == 64:  #CQUAD8
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 5, is_mag_phase, 'CEN/8')
            #elif self.element_type == 82:  # CQUADR
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 5, is_mag_phase, 'CEN/8')
            #elif self.element_type == 70:  # CTRIAR
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 3, is_mag_phase, 'CEN/3')
            #elif self.element_type == 75:  # CTRIA6
                #self._write_f06_quad4_bilinear_transient(f06_file, itime, 3, is_mag_phase, 'CEN/6')
            #else:
                #raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexLayeredCompositeStressVMArray(ComplexLayeredCompositesVMArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        ComplexLayeredCompositesVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        assert self.is_stress, self.stress_bits

    @property
    def headers(self) -> list[str]:
        headers = ['o11', 'o22', 't12', 't1z', 't2z', 'ovm']
        return headers


class ComplexLayeredCompositeStrainVMArray(ComplexLayeredCompositesVMArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        ComplexLayeredCompositesVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        assert self.is_strain, self.stress_bits

    @property
    def headers(self) -> list[str]:
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z', 'evm']
        return headers


class ComplexLayeredCompositeStressArray(ComplexLayeredCompositesArray, StressObject):
    def __init__(self, data_code, is_sort1: bool, isubcase: int, dt):
        ComplexLayeredCompositesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)
        assert self.is_stress, self.stress_bits

    @property
    def headers(self) -> list[str]:
        headers = ['o11', 'o22', 't12', 't1z', 't2z']
        return headers


class ComplexLayeredCompositeStrainArray(ComplexLayeredCompositesArray, StrainObject):
    def __init__(self, data_code, is_sort1: bool, isubcase: int, dt):
        ComplexLayeredCompositesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    @property
    def headers(self) -> list[str]:
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z']
        return headers


def build_dataframe(self: (ComplexLayeredCompositeStressArray | ComplexLayeredCompositeStressVMArray |
                           ComplexLayeredCompositeStrainArray | ComplexLayeredCompositeStrainVMArray)):
    """creates a pandas dataframe"""
    headers = self._get_headers()
    column_names, column_values = build_dataframe_transient_header(self)
    data_frame = build_pandas_transient_element_node(
        self, column_values, column_names,
        headers, self.element_layer, self.data,
        names=['ElementID', 'Layer', 'Item'])
    #print(data_frame)
    return data_frame
