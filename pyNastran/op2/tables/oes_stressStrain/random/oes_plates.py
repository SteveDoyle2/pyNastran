from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_float_13e

BASIC_RANDOM_TABLES = [
    'OESATO1', 'OESCRM1', 'OESPSD1', 'OESRMS1', 'OESNO1',
    'OESATO2', 'OESCRM2', 'OESPSD2', 'OESRMS2', 'OESNO2',
    'OSTRATO1', 'OSTRCRM1', 'OSTRPSD1', 'OSTRRMS1', 'OSTRNO1',
    'OSTRATO2', 'OSTRCRM2', 'OSTRPSD2', 'OSTRRMS2', 'OSTRNO2',
]

class RandomPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
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
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        return get_nnodes(self)

    def build(self):
        """sizes the vectorized attributes of the ComplexPlateArray"""
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        nnodes = self.get_nnodes()

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #print('element_type=%r ntimes=%s nelements=%s nnodes=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, nnodes, self.ntotal, self.subtitle))

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

        self.element = zeros((self.ntotal, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes * 2 == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        self.fiber_curvature = zeros(self.ntotal, 'float32')

        # [oxx, oyy, txy]
        nresults = 3
        if self._is_nx_random():
            # ovm
            nresults += 1

        self.data = zeros((self.ntimes, self.ntotal, nresults), 'float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values,
                                   major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid, Nid\n' % str(self.code_information())
            for eid1, eid2 in zip(self.element, table.element):
                msg += '(%s, %s), (%s, %s)\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (oxx1, oyy1, txy1) = t1
                        (oxx2, oyy2, txy2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                                [oxx1, oyy1, txy1], # atol=0.0001
                                [oxx2, oyy2, txy2], atol=0.075):
                            ni = len(str(eid)) + len(str(eid))
                        #if not np.array_equal(t1, t2):
                            msg += ('%s  (%s, %s, %s)\n'
                                    '%s  (%s, %s, %s)\n' % (
                                        eid,
                                        oxx1, oyy1, txy1,
                                        ' ' * ni,
                                        oxx2, oyy2, txy2,
                                    ))
                            msg += ('%s  (%s, %s, %s)\n'
                                    % (
                                        ' ' * ni,
                                        oxx1 - oxx2,
                                        oyy1 - oyy2,
                                        txy1 - txy2,
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

    def add_new_eid_sort1(self, dt, eid, fd, oxx, oyy, txy):
        self.add_eid_sort1(dt, eid, fd, oxx, oyy, txy)

    def add_sort1(self, dt, eid, fd, oxx, oyy, txy):
        """unvectorized method for adding SORT1 transient data"""
        self.add_eid_sort1(dt, eid, fd, oxx, oyy, txy)

    def add_eid_sort1(self, dt, eid, fd, oxx, oyy, txy):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.data[self.itime, self.itotal] = [oxx, oyy, txy]
        self.element[self.itotal, :] = eid  # 0 is center
        self.fiber_curvature[self.itotal] = fd
        #self.ielement += 1
        self.itotal += 1
    #---------------------------------------------------------------------------

    def add_new_eid_ovm_sort1(self, dt, eid, fd, oxx, oyy, txy, ovm):
        self.add_eid_ovm_sort1(dt, eid, fd, oxx, oyy, txy, ovm)

    def add_ovm_sort1(self, dt, eid, fd, oxx, oyy, txy, ovm):
        """unvectorized method for adding SORT1 transient data"""
        self.add_eid_ovm_sort1(dt, eid, fd, oxx, oyy, txy, ovm)

    def add_eid_ovm_sort1(self, dt, eid, fd, oxx, oyy, txy, ovm):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, dt, eid, node_id, oxx))

        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.data[self.itime, self.itotal] = [oxx, oyy, txy, ovm]
        self.element[self.itotal, :] = eid  # 0 is center
        self.fiber_curvature[self.itotal] = fd
        #self.ielement += 1
        self.itotal += 1

    #---------------------------------------------------------------------------

    def add_new_node_sort1(self, dt, eid, fd, oxx, oyy, txy):
        self.add_eid_sort1(dt, eid, fd, oxx, oyy, txy)

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self._write_f06_quad4_bilinear_transient(f06_file, itime)
            elif self.element_type == 33:  # CQUAD4 linear
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 74: # CTRIA3
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 64:  #CQUAD8
                self._write_f06_quad4_bilinear_transient(f06_file, itime)
            elif self.element_type == 82:  # CQUADR
                self._write_f06_quad4_bilinear_transient(f06_file, itime)
            elif self.element_type == 70:  # CTRIAR
                self._write_f06_quad4_bilinear_transient(f06_file, itime)
            elif self.element_type == 75:  # CTRIA6
                self._write_f06_quad4_bilinear_transient(f06_file, itime)
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_f06_tri3_transient(self, f06_file, itime):
        """
        CQUAD4 linear
        CTRIA3
        """
        fds = self.fiber_curvature
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element

        ilayer0 = True
        for eid, fd, oxx, oyy, txy in zip(eids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(oxx)
            soyy = write_float_13e(oyy)
            stxy = write_float_13e(txy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

    def _write_f06_quad4_bilinear_transient(self, f06_file, itime):
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

        eids = self.element

        ilayer0 = True
        for eid, fd, oxx, oyy, dxy in zip(eids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(oxx)
            soyy = write_float_13e(oyy)
            stxy = write_float_13e(txy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

    def _is_nx_random(self):
        if self.table_name in BASIC_RANDOM_TABLES:
            is_nx_random = False
        elif self.table_name in ['OESXRMS1', ]:
            is_nx_random = True
        else:
            msg = 'self.table_name=%s self.table_name_str=%s' % (self.table_name, self.table_name_str)
            raise NotImplementedError(msg)
        return is_nx_random

def _get_plate_msg(self, is_mag_phase=True, is_sort1=True):
    #if self.is_von_mises:
        #von_mises = 'VON MISES'
    #else:
        #von_mises = 'MAX SHEAR'

    if self.is_stress:
        if self.is_fiber_distance:
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
        if self.is_fiber_distance:
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

    ## TODO: validation on header formatting...

    if self.is_stress:
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
        msg += ctria3 + mag_real + fiber_msg_temp
    elif self.element_type == 75:  # CTRIA6
        msg += ctria6 + mag_real + grid_msg_temp
        is_bilinear = True
    elif self.element_type == 70:  # CTRIAR
        msg += ctriar + mag_real + grid_msg_temp
        is_bilinear = True
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

    nnodes = get_nnodes(self)
    msg = [
        '                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
        '                                             ( POWER SPECTRAL DENSITY FUNCTION )\n'
        ' \n'
        '                    FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n'
        '    FREQUENCY      DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n'
        #'0  2.000000E+01 -5.000000E-02              1.925767E-05                      1.404795E-04                     1.097896E-03'
        #'                 5.000000E-02              1.925766E-05                      1.404794E-04                     1.097896E-03'
    ]
    return msg, nnodes, is_bilinear

def get_nnodes(self):
    if self.element_type in [64, 82, 144]:  # ???, CQUADR, CQUAD4 bilinear
        nnodes = 4 + 1 # centroid
    elif self.element_type in [70, 75]:   #???, CTRIA6
        nnodes = 3 + 1 # centroid
    elif self.element_type in [144, 74, 33]:  # CTRIA3, CQUAD4 linear
        nnodes = 1
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
    return nnodes

class RandomPlateStressArray(RandomPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        headers = ['oxx', 'oyy', 'txy']
        if self._is_nx_random():
            headers.append('ovm')
        return headers

    def get_headers(self):
        return self._get_headers()

class RandomPlateStrainArray(RandomPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        headers = ['exx', 'eyy', 'exy']
        if self._is_nx_random():
            headers.append('ovm')
        return headers

    def get_headers(self):
        return self._get_headers()
