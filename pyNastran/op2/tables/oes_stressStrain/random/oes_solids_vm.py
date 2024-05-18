# pylint: disable=C0301,R0913,R0914,R0904,C0111,R0201,R0902
from itertools import count

import numpy as np
from numpy import zeros, where, searchsorted

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
from pyNastran.op2.result_objects.op2_objects import get_times_dtype


class RandomSolidVMArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def get_headers(self):
        raise NotImplementedError()

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def build(self):
        """sizes the vectorized attributes of the RealSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        if self.is_sort1:
            ntimes = self.ntimes
            nelements = self.nelements
            ntotal = self.ntotal
        else:
            nelements = self.ntimes
            ntimes = self.nelements
            ntotal = self.ntotal
            dtype = self._get_analysis_code_dtype()

        self._times = zeros(ntimes, dtype=self.analysis_fmt)

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((ntotal, 2), dtype=idtype)
        self.element_cid = zeros((nelements, 2), dtype=idtype)

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = zeros((self.ntotal, nnodes, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz, vm]
        self.data = zeros((self.ntimes, self.ntotal, 7), 'float32')
        self.nnodes = self.element_node.shape[0] // self.nelements
        #self.data = zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        # TODO: cid?
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']

    def add_eid_sort1(self, unused_etype, cid, dt, eid, unused_node_id,
                      oxx, oyy, ozz, txy, tyz, txz, vm):
        assert cid >= -2, cid
        assert eid >= 0, eid

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center

        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz]
        #self.data[self.itime, self.ielement, 0, :] = [oxx, oyy, ozz, txy, tyz, txz]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        if self.ielement == self.nelements:
            self.ielement = 0
        self.element_cid[self.ielement, :] = [eid, cid]
        self.itotal += 1
        self.ielement += 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_nid in enumerate(self.element_node):
                    (eid, nid) = eid_nid
                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, ozz1, txy1, tyz1, txz1, vm1) = t1
                    (oxx2, oyy2, ozz2, txy2, tyz2, txz2, vm2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s)\n'
                            '%s      (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid,
                                oxx1, oyy1, ozz1, txy1, tyz1, txz1, vm1,
                                ' ' * (len(str(eid)) + len(str(nid)) + 2),
                                oxx2, oyy2, ozz2, txy2, tyz2, txz2, vm2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_node_sort1(self, dt, eid, unused_inode, node_id, oxx, oyy, ozz, txy, tyz, txz, vm):
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, vm]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        #self.data[self.itime, self.ielement-1, self.inode, :] = [oxx, oyy, ozz, txy, tyz, txz]

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        self.element_node[self.itotal, :] = [eid, node_id]
        #self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
        self.itotal += 1

    @property
    def nnodes_per_element(self) -> int:
        if self.element_type == 39: # CTETRA
            nnodes = 4
        elif self.element_type == 67: # CHEXA
            nnodes = 8
        elif self.element_type == 68: # CPENTA
            nnodes = 6
        elif self.element_type == 255:  # CPYRAM
            nnodes = 5
        else:  # pragma: no cover
            raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))
        return nnodes

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        try:
            nnodes_per_element = self.element_node.shape[0] // nelements
        except ZeroDivisionError:
            nnodes_per_element = '???'
        nnodes = self.element_node.shape[0]

        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n  nnodes_per_element=%s (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n  nodes_per_element=%i (including centroid)\n'
                       % (self.__class__.__name__, nelements, nnodes, nnodes_per_element))
            ntimes_word = '1'
        msg.append('  eType, cid\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  element_cid.shape = {self.element_cid.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        #print(''.join(msg))
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element_node[:, 0])
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp: str='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        nnodes, msg_temp = _get_f06_header_nnodes(self, is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        unused_cids3 = self.element_cid[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]
            vm = self.data[itime, :, 6]

            cnnodes = nnodes + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, dvm in zip(
                count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, vm):

                unused_j = where(eids3 == deid)[0]
                #cid = cids3[j]
                cid = 0
                [oxxi, oyyi, ozzi, txyi, tyzi, txzi, vmi] = write_floats_13e(
                    [doxx, doyy, dozz, dtxy, dtyz, dtxz, dvm])

                if i % cnnodes == 0:
                    #print('deid, cid, nnodes = %s, %s, %s' % (deid, cid, nnodes))
                    f06_file.write('0  %8s    %8iGRID CS  %i GP\n' % (deid, cid, nnodes))
                    f06_file.write(
                        #               center   oxx    oyy     ozz     txy     tyz     txz
                        '0              %8s   %-13s   %-13s   %-13s   %-13s   %-13s   %-13s   %-13s\n'
                        % ('CENTER', oxxi, oyyi, ozzi, txyi, tyzi, txzi, vmi))
                else:
                    f06_file.write(
                        '0              %8s   %-13s   %-13s   %-13s   %-13s   %-13s   %-13s   %-13s\n'
                        % (node_id, oxxi, oyyi, ozzi, txyi, tyzi, txzi, vmi))
                i += 1
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RandomSolidVMStressArray(RandomSolidVMArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomSolidVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'ovm']
        return headers


class RandomSolidVMStrainArray(RandomSolidVMArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomSolidVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', 'evm']
        return headers

def _get_solid_msgs(self):
    # '                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )'
    # '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  1.00E+00 )'
    # '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------'
    # '      ELEMENT-ID   GRID-ID      NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX   VON MISES'
    # '0            6701        0GRID CS  8 GP'
    # '0                   CENTER     6.413E+01   2.752E-01   6.413E+01     6.496E+00   6.496E+00   3.842E+01   9.359E+01'
    # '0                    30000     5.367E+00   1.069E+01   5.367E+00     1.134E-15   1.134E-15   4.473E-15   5.321E+00'

    assert self.table_name in ['OESXNO1', 'OESXRMS1'], self.table_name
    if self.is_stress:
        base_msg = [
            '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  1.00E+00 )',
            '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------\n',
            '      ELEMENT-ID   GRID-ID      NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX   VON MISES\n',
        ]
        tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        pyram_msg = ['                    S T R E S S E S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
    else:
        base_msg = [
            '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  1.00E+00 )',
            '0                   CORNER      --------------------------CENTER AND CORNER POINT STRAINs---------------------------\n',
            '      ELEMENT-ID   GRID-ID      NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX   VON MISES\n',
        ]
        tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        pyram_msg = ['                      S T R A I N S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]
        hexa_msg  = ['                        S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
    tetra_msg += base_msg
    penta_msg += base_msg
    pyram_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, pyram_msg, penta_msg, hexa_msg

def _get_f06_header_nnodes(self, is_mag_phase=True) -> tuple[int, str]:
    tetra_msg, pyram_msg, penta_msg, hexa_msg = _get_solid_msgs(self)
    if self.element_type == 39: # CTETRA
        msg = tetra_msg
        nnodes = 4
    elif self.element_type == 67: # CHEXA
        msg = hexa_msg
        nnodes = 8
    elif self.element_type == 68: # CPENTA
        msg = penta_msg
        nnodes = 6
    elif self.element_type == 255:  # CPYRAM
        msg = pyram_msg
        nnodes = 5
    else:  # pragma: no cover
        msg = 'element_name=%s self.element_type=%s' % (self.element_name, self.element_type)
        raise NotImplementedError(msg)
    return nnodes, msg
