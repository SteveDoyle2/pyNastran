from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types
from six.moves import range

import numpy as np
from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, write_imag_floats_13e, get_key0, write_float_13e
try:
    import pandas as pd
except ImportError:
    pass


class ComplexPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        return get_nnodes(self)

    def build(self):
        #print('data_code = %s' % self.data_code)
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = self.get_nnodes()

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        #self.element_cid = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes * 2 == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        self.fiber_curvature = zeros(self.ntotal, 'float32')
        # [oxx, oyy, txy]
        self.data = zeros((self.ntimes, self.ntotal, 3), 'complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element_node, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
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
            if self.is_sort1():
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
                                oxx1.real, oxx1.imag, oyy1.real, oyy1.imag, txy1.real, txy1.imag,
                                ' ' * ni,
                                oxx2.real, oxx2.imag, oyy2.real, oyy2.imag, txy2.real, txy2.imag,
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
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_new_eid_sort1(self, dt, eid, node_id, fdr, oxx, oyy, txy):
        self.add_eid_sort1(dt, eid, node_id, fdr, oxx, oyy, txy)

    def add_sort1(self, dt, eid, gridC, fdr, oxx, oyy, txy):
        self.add_eid_sort1(dt, eid, gridC, fdr, oxx, oyy, txy)

    def add_new_node_sort1(self, dt, eid, gridc, fdr, oxx, oyy, txy):
        self.add_eid_sort1(dt, eid, gridc, fdr, oxx, oyy, txy)

    def add_eid_sort1(self, dt, eid, node_id, fdr, oxx, oyy, txy):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(node_id, int), node_id
        self.data[self.itime, self.itotal] = [oxx, oyy, txy]
        self.element_node[self.itotal, :] = [eid, node_id]  # 0 is center
        self.fiber_curvature[self.itotal] = fdr
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self._write_f06_quad4_bilinear_transient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.element_type == 33:  # CQUAD4 linear
                self._write_f06_tri3_transient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.element_type == 74: # CTRIA3
                self._write_f06_tri3_transient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 64:  #CQUAD8
                self._write_f06_quad4_bilinear_transient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type == 82:  # CQUADR
                self._write_f06_quad4_bilinear_transient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type == 70:  # CTRIAR
                self._write_f06_quad4_bilinear_transient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 75:  # CTRIA6
                self._write_f06_quad4_bilinear_transient(f, itime, 3, is_mag_phase, 'CEN/6')
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_f06_tri3_transient(self, f, itime, n, is_magnitude_phase, cen):
        """
        CQUAD4 linear
        CTRIA3
        """
        fds = self.fiber_curvature
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        ilayer0 = True
        for eid, node, fdr, doxx, doyy, dtxy in zip(eids, nodes, fds, oxx, oyy, txy):
            vals = write_float_13e(fdr)
            fdr = vals[0]
            [oxxr, oyyr, txyr,
             oxxi, oyyi, txyi,] = write_imag_floats_13e([doxx, doyy, dtxy], is_magnitude_phase)

            if ilayer0:    # TODO: assuming 2 layers?
                f.write('0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0

    def _write_f06_quad4_bilinear_transient(self, f, itime, n, is_magnitude_phase, cen):
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        fds = self.fiber_curvature
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        ilayer0 = True
        for eid, node, fd, doxx, doyy, dtxy in zip(eids, nodes, fds, oxx, oyy, txy):
            fdr = write_float_13e(fd)
            [oxxr, oyyr, txyr,
             oxxi, oyyi, txyi,] = write_imag_floats_13e([doxx, doyy, dtxy], is_magnitude_phase)

            if node == 0 and ilayer0:
                f.write('0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % (eid, cen, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            elif ilayer0:    # TODO: assuming 2 layers?
                f.write('   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % ('', node, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n\n' % ('', '', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0

def _get_plate_msg(self, is_mag_phase=True, is_sort1=True):
    if self.is_von_mises():
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress():
        if self.is_fiber_distance():
            grid_msg_temp = ['    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            grid_msg_temp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
    else:
        if self.is_fiber_distance():
            grid_msg_temp = ['    ELEMENT              FIBER                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            grid_msg_temp = ['    ELEMENT              FIBRE                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']


    if is_mag_phase:
        mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
    else:
        mag_real = ['                                                          (REAL/IMAGINARY)\n', ' \n']

    #if self.is_fiber_distance():
        #quadMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -',
                       #'      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY']
    #else:
        #pass

#'0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
#'                       2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'

    #quadMsg  = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'] + formWord + quadMsgTemp
    #quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )'] + formWord + quadMsgTemp
    #quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'] + formWord + quadMsgTemp

    ## TODO: validation on header formatting...
    if self.is_stress():
        cquad4_bilinear = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']
        cquad4_linear = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad8 = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria3 = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        ctria6 = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']
    else:
        cquad4_bilinear = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']
        cquad4_linear = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']
        cquad8 = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria3 = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']
        ctria6 = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

    msg = []
    is_bilinear = False
    if self.element_type == 144: # CQUAD4
        is_bilinear = True
        msg += cquad4_linear + mag_real + grid_msg_temp
    elif self.element_type == 33: # CQUAD4
        is_bilinear = False
        msg += cquad4_bilinear + mag_real + fiber_msg_temp
    elif self.element_type == 64:  #CQUAD8
        msg += cquad8 + mag_real + grid_msg_temp
        is_bilinear = True
    elif self.element_type == 82:  # CQUADR
        msg += cquadr + mag_real + grid_msg_temp
        is_bilinear = True

    elif self.element_type == 74: # CTRIA3
        msg += ctria3 + fiber_msg_temp
    elif self.element_type == 75:  # CTRIA6
        msg += ctria6 + grid_msg_temp
        is_bilinear = True
    elif self.element_type == 70:  # CTRIAR
        msg += ctriar + grid_msg_temp
        is_bilinear = True
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

    nnodes = get_nnodes(self)
    return msg, nnodes, is_bilinear

def get_nnodes(self):
    if self.element_type in [64, 82, 144]:  # ???, CQUADR, CQUAD4 bilinear
        nnodes = 4 # + 1 centroid
    elif self.element_type in [70, 75]:   #???, CTRIA6
        nnodes = 3 # + 1 centroid
    elif self.element_type in [74, 33]:  # CTRIA3, CQUAD4 linear
        nnodes = 1
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
    return nnodes + 1  # + 1 centroid

class ComplexPlateStressArray(ComplexPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['oxx', 'oyy', 'txy']

    def get_headers(self):
        return self._get_headers()

class ComplexPlateStrainArray(ComplexPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain(), self.stress_bits

    def _get_headers(self):
        return ['exx', 'eyy', 'exy']

    def get_headers(self):
        return self._get_headers()
