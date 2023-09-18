import warnings

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_imag_floats_13e, write_float_13e

#BASIC_TABLES = {
    #'OES1X', 'OES1',
    #'OES2',
    #'OSTR1X',
#}
#VM_TABLES = {'OESVM1', 'OESVM2',
             #'OSTRVM1', 'OSTRVM2'}

class ComplexLayeredCompositesArray(OES_Object):
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
        #print(self._ntotals, self.ntotal)
        #print(self.code_information())

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #print('element_type=%r ntimes=%s nelements=%s nnodes=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, nnodes, self.ntotal, self.subtitle))

        #self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        idtype, cfdtype = get_complex_times_dtype(self.size)

        if self.is_sort1:
            ntimes = self.ntimes
            nlayers = self.ntotal
            #print(f'  SORT1: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        else:
            raise NotImplementedError(self.code_information())
        #elif self.is_sort2:
            #nelements = self.ntimes
            #nlayers = nelements * 2 * nnodes
            #ntimes = self.ntotal
            #print(f'  SORT2: ntimes={ntimes} nlayers={nlayers} {self.element_name}-{self.element_type}')
        #print("nelements=%s nlayers=%s ntimes=%s" % (nelements, nlayers, ntimes))

        self._times = zeros(ntimes, dtype=self.analysis_fmt)
        #self.ntotal = self.nelements * nnodes

        # the number is messed up because of the offset for the element's properties
        #if not self.nelements * nnodes * 2 == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                #self.ntimes, self.nelements, nnodes,
                #self.nelements * nnodes, self.ntotal)
            #raise RuntimeError(msg)

        # [o1a, o2a, t12a, o1za, o2za,
        # o1b, o2b, t12b, o1zb, e2zb, ovm]
        self.data = zeros((ntimes, nlayers, 11), dtype=cfdtype)

        self.element_layer = zeros((nlayers, 2), dtype=idtype)
        #print(self.data.shape, self.element_node.shape)

    #def build_dataframe(self) -> None:
        #"""creates a pandas dataframe"""
        #headers = self.get_headers()
        #column_names, column_values = self._build_dataframe_transient_header()

        #data_frame = self._build_pandas_transient_element_node(
            #column_values, column_names,
            #headers, self.element_node, self.data)
        ##print(data_frame)
        #self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        return True

    def _get_headers(self):
        headers = ['o1a', 'o2a', 't12a', 'o1za', 'o2za',
                   'o1b', 'o2b', 't12b', 'o1zb', 'e2zb', 'ovm']
        return headers

    def add_sort1(self, dt, eid, ply_id,
                  o1a, o2a, t12a, o1za, o2za,
                  o1b, o2b, t12b, o1zb, e2zb, ovm):
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(ply_id, int), ply_id
        self.data[self.itime, self.itotal] = [o1a, o2a, t12a, o1za, o2za,
                                              o1b, o2b, t12b, o1zb, e2zb, ovm]
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.itotal += 1

    def add_sort1_real(self, dt, eid, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear) -> None:
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(node_id, int), node_id
        self.data[self.itime, self.itotal] = [oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear]
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.itotal += 1
        #self.ielement += 1

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
        if self.is_von_mises:
            warnings.warn(f'{self.class_name} doesnt support writing von Mises')
            f06_file.write(f'{self.class_name} doesnt support writing von Mises\n')

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
            # [o1a, o2a, t12a, o1za, o2za,
            # o1b, o2b, t12b, o1zb, e2zb, ovm]
            #self.data = zeros((ntimes, nlayers, 11), dtype=cfdtype)

            #self.element_layer = zeros((nlayers, 2), dtype=idtype)

            #fds = self.fiber_distance
            o1a = self.data[itime, :, 0]
            o2a = self.data[itime, :, 1]
            t12a = self.data[itime, :, 2]
            o1za = self.data[itime, :, 3]
            o2za = self.data[itime, :, 4]

            o1b = self.data[itime, :, 5]
            o2b = self.data[itime, :, 6]
            t12b = self.data[itime, :, 7]
            o1zb = self.data[itime, :, 8]
            e2zb = self.data[itime, :, 9]

            ovm = self.data[itime, :, 10]
            fds = o1a

            eids = self.element_layer[:, 0]
            layer = self.element_layer[:, 1]

            ilayer0 = True
            for eid, layer, fd, do1a, do2a, dt12a in zip(eids, layer, fds, o1a, o2a, t12a):
                fdr = write_float_13e(fd.real)
                [do1ar, do2ar, dt12ar,
                 do1ai, do2ai, dt12ai,] = write_imag_floats_13e(
                     [do1a, do2a, dt12a], is_mag_phase)

                #print(do1ar, do2ar, dt12ar, do1ai, do2ai, dt12ai)
                ilayer0 = not ilayer0
                continue
                #if node == 0 and ilayer0:
                    #f06_file.write('0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % (
                        #eid, cen, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
                #elif ilayer0:    # TODO: assuming 2 layers?
                    #f06_file.write('   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % (
                        #'', node, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
                #else:
                    #f06_file.write('   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n\n' % (
                        #'', '', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
                ilayer0 = not ilayer0

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexLayeredCompositeStressArray(ComplexLayeredCompositesArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexLayeredCompositesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)
        assert self.is_stress, self.stress_bits

class ComplexLayeredCompositeStrainArray(ComplexLayeredCompositesArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexLayeredCompositesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

def _get_composite_plate_msg(self, is_mag_phase=True, is_sort1=True) -> tuple[list[str], int]:
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

