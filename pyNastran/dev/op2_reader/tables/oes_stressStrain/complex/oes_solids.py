#pylint: disable=C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import numpy as np
from numpy import zeros, concatenate
try:
    import pandas as pd  # type: ignore
except ImportError:
    pass

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexSolidArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def combine(self, results):
        #print('ComplexSolid combine')
        #print('data.shape1 =', self.data.shape)
        #self.data = vstack(data)

        data = [self.nelements] + [result.nelements for result in results]
        self.nelements = sum(data)

        data = [self.ntotal] + [result.ntotal for result in results]
        self.ntotal = sum(data)

        data = [self.element_types3] + [result.element_types3 for result in results]
        self.element_types3 = concatenate(data, axis=0)

        data = [self.element_node] + [result.element_node for result in results]
        self.element_node = concatenate(data, axis=0)

        data = [self.element_cid] + [result.element_cid for result in results]
        self.element_cid = concatenate(data, axis=0)

        data = [self.data] + [result.data for result in results]
        self.data = concatenate(data, axis=1)

        self.data = concatenate(data, axis=1)

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        if self.element_type == 39: # CTETRA
            nnodes = 5
        elif self.element_type == 68: # CPENTA
            nnodes = 7
        elif self.element_type == 67: # CHEXA
            nnodes = 9
        else:
            raise NotImplementedError(self.element_name)
        return nnodes

    def build(self):
        """sizes the vectorized attributes of the ComplexSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = self.get_nnodes()

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.element_types2 = array(self.nelements, dtype='|S8')
        #self.element_types3 = zeros((self.nelements, 2), dtype='int32')


        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        self.element_cid = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        if self.result_flag == 0:
            # [oxx, oyy, ozz, txy, tyz, txz]
            self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')
        else:
            # oxx
            self.data = zeros((self.ntimes, self.ntotal, 1), 'complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid_nid in enumerate(self.element_node):
                        eid, nid = eid_nid
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        d = t1 - t2
                        if not np.allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
                                           [tx2.real, tx2.imag, ty2.real, ty2.imag], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                tx1.real, tx1.imag, ty1.real, ty1.imag,
                                tx2.real, tx2.imag, ty2.real, ty2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
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

    def add_eid_sort1(self, element_num, element_type, dt, eid, cid, ctype, nodef):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #self.element_types2[self.ielement] = string_(element_type)   # TODO: save this...
        #self.element_types2[self.ielement] = element_type

        #try:
        if self.ielement < self.nelements:
            self.element_cid[self.ielement] = [eid, cid]
            #self.element_types3[self.ielement, :] = [element_num, nodef]
        #except IndexError:
            #pass
            #print('element_types3', self.element_types3)

        #self.node_element_cid[self.itotal] = []
        #self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        #print("etype=%s ctype=%s nodef=%s" % (element_type, ctype, nodef))
        self.ielement += 1
        #self.itotal += 1

    def add_node_sort1(self, dt, eid, grid, inode, ex, ey, ez, etxy, etyz, etzx):
        if self.result_flag == 0:
            self.data[self.itime, self.itotal, :] = [ex, ey, ez, etxy, etyz, etzx]
        else:
            self.data[self.itime, self.itotal, 0] = ex
        self.element_node[self.itotal, :] = [eid, grid]
        self.itotal += 1

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nnodes = self.element_node.shape[0]
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  element_cid.shape = %s\n' % str(self.element_cid.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06
        ntimes = self.data.shape[0]

        cid = 0
        for itime in range(ntimes):
            dt = self._times[itime]

            #print('eids=', eids)

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]

            eids2 = self.element_node[:, 0]
            nodes = self.element_node[:, 1]

            for deid, node, doxx, doyy, dozz, dtxy, dtyz, dtxz in zip(eids2, nodes, oxx, oyy, ozz, txy, tyz, txz):
                # TODO: cid not supported
                [oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                 oxxi, oyyi, ozzi, txyi, tyzi, txzi,] = write_imag_floats_13e([doxx, doyy, dozz,
                                                                               dtxy, dtyz, dtxz], is_mag_phase)
                if node == 0:  # CENTER
                    f06_file.write(
                        '0 %12i %11sGRID CS %2i GP\n'
                        '0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n'
                        '    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (
                            deid, cid, nnodes,
                            'CENTER', oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                            '', oxxi, oyyi, ozzi, txyi, tyzi, txzi,
                    ))
                else:
                    f06_file.write(
                        '0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n'
                        '    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (
                            node, oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                            '', oxxi, oyyi, ozzi, txyi, tyzi, txzi,
                    ))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSolidStressArray(ComplexSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        return headers

def _get_msgs(self, is_mag_phase, is_sort1):
    if is_mag_phase:
        mag_phase = '                                                          (MAGNITUDE/PHASE)'
    else:
        mag_phase = '                                                          (REAL/IMAGINARY)'

    if self.is_stress:
        base_msg = [
            mag_phase,
            '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )', ]
        hexa_msg = ['                 C O M P L E X   S T R E S S E S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )', ]
        penta_msg = ['                 C O M P L E X   S T R E S S E S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )', ]
    else:
        base_msg = [
            mag_phase,
            '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',]
        hexa_msg = ['                 C O M P L E X     S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )',]
        penta_msg = ['                 C O M P L E X     S T R A I N S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )',]

    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg


def get_f06_header(self, is_mag_phase=True, is_sort1=True):
    tetra_msg, penta_msg, hexa_msg = _get_msgs(self, is_mag_phase, is_sort1)

    if self.element_type == 39:  # CTETRA
        return tetra_msg, 4
    elif self.element_type == 67:  # CHEXA
        return hexa_msg, 8
    elif self.element_type == 68:  # CPENTA
        return penta_msg, 6
    else:
        raise NotImplementedError('complex solid stress/strain name=%r Type=%s' % (self.element_name, self.element_type))


class ComplexSolidStrainArray(ComplexSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers
