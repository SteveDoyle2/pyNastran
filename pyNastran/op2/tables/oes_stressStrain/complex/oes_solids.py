#pylint: disable=C0301,C0111
from typing import List

import numpy as np
from numpy import zeros, concatenate

from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexSolidArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
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

    @property
    def nnodes_per_element(self):
        if self.element_type == 39: # CTETRA
            nnodes = 5
        elif self.element_type == 68: # CPENTA
            nnodes = 7
        elif self.element_type == 67: # CHEXA
            nnodes = 9
        elif self.element_type == 255: # CPYRAM
            nnodes = 6
        else:
            raise NotImplementedError(self.element_name)
        return nnodes

    def build(self):
        """sizes the vectorized attributes of the ComplexSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        nnodes = self.nnodes_per_element

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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        #self.element_types2 = array(self.nelements, dtype='|S8')
        #self.element_types3 = zeros((self.nelements, 2), dtype='int32')


        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), dtype=idtype)
        self.element_cid = zeros((self.nelements, 2), dtype=idtype)

        # the number is messed up because of the offset for the element's properties

        #if not self.nelements * nnodes == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           #self.nelements, nnodes,
                                                                           #self.nelements * nnodes,
                                                                           #self.ntotal)
            #raise RuntimeError(msg)

        # [oxx, oyy, ozz, txy, tyz, txz]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        # Freq                  0.00001  10.00000 20.00000 30.00000                 40.00000 50.00000 60.00000
        # ElementID NodeID Item
        # 1         0      oxx        0j       0j       0j       0j    (3200.0806+6017.714j)       0j       0j
        #                  oyy        0j       0j       0j       0j    (410.68146+772.2816j)       0j       0j
        #                  ozz        0j       0j       0j       0j    (0.306115+0.5756457j)       0j       0j
        #                  txy        0j       0j       0j       0j  (-120.69606-226.96753j)       0j       0j
        #                  tyz        0j       0j       0j       0j  (0.70554054+1.3267606j)       0j       0j
        #                  txz        0j       0j       0j       0j     (5193.834+9766.943j)       0j       0j
        # 2                oxx        0j       0j       0j       0j    (8423.371+15840.051j)       0j       0j
        #                  oyy        0j       0j       0j       0j    (-3364.359-6326.637j)       0j       0j
        #                  ozz        0j       0j       0j       0j  (-74931.664-140908.11j)       0j       0j
        #                  txy        0j       0j       0j       0j  (-261.20972-491.20178j)       0j       0j
        #                  tyz        0j       0j       0j       0j   (121.57285+228.61633j)       0j       0j
        #                  txz        0j       0j       0j       0j     (5072.678+9539.112j)       0j       0j
        #import pandas as pd
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_element_node(column_values, column_names,
                                                               headers, self.element_node, self.data)

        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        #data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        #print(data_frame)
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
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
        self.data[self.itime, self.itotal, :] = [ex, ey, ez, etxy, etyz, etzx]
        self.element_node[self.itotal, :] = [eid, grid]
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
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

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (
                self.__class__.__name__, nelements, nnodes, self.table_name))
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

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.element_cid.shape[0]
        #print(self.element_cid)
        #print(self.element_node)

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids = self.element_node[:, 0]
        eids_device = eids * 10 + self.device_code

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 4s 2i 12f')
            struct2 = Struct(endian + b'i 12f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)

        #cid = 0
        unused_msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1=True)
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
            op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]

            eids2 = self.element_node[:, 0]
            nodes = self.element_node[:, 1]
            assert len(eids2) == len(oxx)
            #print(len(eids2), self.num_wide)
            #(eid_device, cid, ctype, nodef) = out
            #(grid,
             #exr, eyr, ezr, etxyr, etyzr, etzxr,
            #exi, eyi, ezi, etxyi, etyzi, etzxi) = out

            nwide = 0
            ielement = -1
            for eid_device, deid, node, doxx, doyy, dozz, dtxy, dtyz, dtxz in zip(eids_device, eids2, nodes, oxx, oyy, ozz, txy, tyz, txz):

                if node == 0:  # CENTER
                    ielement += 1
                    cid = self.element_cid[ielement, 1]
                    data = [eid_device, cid, b'GRID', nnodes, node,
                            doxx.real, doyy.real, dozz.real, dtxy.real, dtyz.real, dtxz.real,
                            doxx.imag, doyy.imag, dozz.imag, dtxy.imag, dtyz.imag, dtxz.imag]
                    #op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, tuple(data)))
                    op2.write(struct1.pack(*data))

                    op2_ascii.write(
                        '0 %12i %11sGRID CS %2i GP\n'
                        '0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n'
                        '    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (
                            deid, cid, nnodes,
                            'CENTER', doxx.real, doyy.real, dozz.real, dtxy.real, dtyz.real, dtxz.real,
                            '', doxx.imag, doyy.imag, dozz.imag, dtxy.imag, dtyz.imag, dtxz.imag,
                    ))
                else:
                    data = [node,
                            doxx.real, doyy.real, dozz.real, dtxy.real, dtyz.real, dtxz.real,
                            doxx.imag, doyy.imag, dozz.imag, dtxy.imag, dtyz.imag, dtxz.imag]
                    op2.write(struct2.pack(*data))
                    op2_ascii.write(
                        '0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n'
                        '    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (
                            node, doxx.real, doyy.real, dozz.real, dtxy.real, dtyz.real, dtxz.real,
                            '', doxx.imag, doyy.imag, dozz.imag, dtxy.imag, dtyz.imag, dtxz.imag,
                    ))
                nwide += len(data)
                #print(data)
            assert nwide == ntotal, 'nwide=%s ntotal=%s' % (nwide, ntotal)

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexSolidStressArray(ComplexSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
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
        pyram_msg = ['                                      C O M P L E X   S T R E S S E S   I N   C P Y R A M      E L E M E N T S   ', ]

    else:
        base_msg = [
            mag_phase,
            '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',]
        hexa_msg = ['                 C O M P L E X     S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )',]
        penta_msg = ['                 C O M P L E X     S T R A I N S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )',]
        pyram_msg = ['                                       C O M P L E X   S T R A I N S   I N   C P Y R A M      E L E M E N T S   ', ]
    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    pyram_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg, pyram_msg


def get_f06_header(self, is_mag_phase=True, is_sort1=True):
    tetra_msg, penta_msg, hexa_msg, pyram_msg = _get_msgs(self, is_mag_phase, is_sort1)

    if self.element_type == 39:  # CTETRA
        return tetra_msg, 4
    elif self.element_type == 67:  # CHEXA
        return hexa_msg, 8
    elif self.element_type == 68:  # CPENTA
        return penta_msg, 6
    elif self.element_type == 255:  # CPYRAM
        return pyram_msg, 6
    else:
        raise NotImplementedError('complex solid stress/strain name=%r Type=%s' % (self.element_name, self.element_type))


class ComplexSolidStrainArray(ComplexSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers
