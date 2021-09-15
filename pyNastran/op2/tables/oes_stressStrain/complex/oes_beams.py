from typing import List

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexBeamArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        #if is_sort1:
            ##sort1
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def build(self):
        """sizes the vectorized attributes of the ComplexCBeamArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        #nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))
        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        #self.element = array(self.nelements, dtype='|S8')

        #self.ntotal = self.nelements * nnodes
        if self.is_sort1 == 1:
            ntimes = self.ntimes
            ntotal = self.ntotal
        else:
            ntimes = self.ntotal
            ntotal = self.ntimes
            #print(f'complex beam: ntimes={ntimes} ntotal={ntotal}')

        # ntotal is nelements
        self._times = zeros(ntimes, dtype)
        self.element_node = zeros((ntotal, 2), idtype)
        self.sd = zeros(ntotal, 'float32')

        # the number is messed up because of the offset for the element's properties
        #if self.nelements * nnodes != self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           #self.nelements, nnodes,
                                                                           #self.nelements * nnodes,
                                                                           #self.ntotal)
            #raise RuntimeError(msg)

        #[sxc, sxd, sxe, sxf]
        self.data = zeros((ntimes, ntotal, 4), cfdtype)

    def finalize(self):
        #enode_sum = self.element_node.sum(axis=1)
        eids = self.element_node[:, 0]
        ueids = np.unique(eids)
        #if len(enode_sum) != self.nelements:
            #msg = 'self.element_node.shape=%s len(enode_sum)=%s' % (str(self.element_node.shape), len(enode_sum))
            #raise RuntimeError(msg)

        #ielem_nonzero = np.where(enode_sum != 0)[0]
        ielem_sdzero = np.searchsorted(eids, ueids)
        isd_nonzero = np.where(self.sd > 0.0)[0]

        # returns unique sorted set
        inonzero = np.union1d(ielem_sdzero, isd_nonzero)
        self.nelements = len(inonzero)
        self.element_node = self.element_node[inonzero, :]
        self.sd = self.sd[inonzero]
        self.data = self.data[:, inonzero, :]

    def _build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values,
                                   major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

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
                    for ieid, (eid, grid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (sxc1, sxd1, sxe1, sxf1) = t1
                        (sxc2, sxd2, sxe2, sxf2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                            [sxc1.real, sxd1.real, sxe1.real, sxf1.real],
                            [sxc2.real, sxd2.real, sxe2.real, sxf2.real], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s %-4s (%s, %s, %s, %s)\n      (%s, %s, %s, %s)\n' % (
                                eid, grid,
                                sxc1.real, sxd1.real, sxe1.real, sxf1.real,
                                sxc2.real, sxd2.real, sxe2.real, sxf2.real,
                                )
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

    def add_sort1(self, dt, eid, grid, sd,
                  exc, exd, exe, exf):
        """adds the non-vectorized data

        201 60000 0.0 (2.1205626410392142e-07-7.067835650076404e-09j)
        201 0 0.5 (6.486846615416653e-08-2.162066081723424e-09j)
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 0 0.0 0j
        201 60001 1.0 (-8.231933179558837e-08+2.7437032645849513e-09j)
        202 60001 0.0 (-1.8220836750515446e-07+6.073004765738688e-09j)
        """
        assert isinstance(eid, integer_types) and eid > 0, f'dt={dt} eid={eid} grid={grid}'
        assert isinstance(eid, integer_types) and grid >= 0, f'dt={dt} eid={eid} grid={grid}'
        #print(eid, grid, sd, exc)
        self.element_node[self.itotal] = [eid, grid]
        self.sd[self.itotal] = sd
        self.data[self.itime, self.itotal, :] = [exc, exd, exe, exf]

        self._times[self.itime] = dt
        self.itotal += 1


    def add_sort2(self, dt, eid, grid, sd,
                  exc, exd, exe, exf):
        """adds the non-vectorized data"""
        assert isinstance(eid, integer_types) and eid > 0, f'dt={dt} eid={eid} grid={grid}'
        assert isinstance(eid, integer_types) and grid >= 0, f'dt={dt} eid={eid} grid={grid}'
        itotal = self.itime
        itime = self.itotal
        #print(f'dt={dt} eid={eid} itime={itime} itotal={itotal}')
        self.element_node[itotal] = [eid, grid]
        self.sd[itotal] = sd
        self.data[itime, itotal, :] = [exc, exd, exe, exf]

        self._times[itime] = dt
        self.itotal += 1

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nelements2 = self.element_node.shape[0]
        assert nelements, nelements2
        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, self.table_name))
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, self.table_name))
        msg.append(
            '  eType\n'
            '  data: [ntimes, nnodes, 4] where 4=[%s]\n'
            '  element_node.shape = %s\n'
            '  sd.shape = %s\n'
            '  data.shape = %s\n' % (
                str(', '.join(self.get_headers())),
                str(self.element_node.shape).replace('L', ''),
                str(self.sd.shape).replace('L', ''),
                str(self.data.shape).replace('L', ''),
        ))

        msg.append('  CBEAM\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1,
                  is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)'

        name = self.data_code['name']
        if self.is_sort1:
            line1 = '                    STAT DIST/     LOCATION         LOCATION         LOCATION         LOCATION'
            line2 = '   ELEMENT-ID  GRID   LENGTH          C                D                E                F\n'
        else:
            raise NotImplementedError('sort2 f06 writing')
            #line1 = '                                       LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n'
            #if name == 'freq':
                #line2 = '           FREQUENCY                       1              2              3              4             AXIAL STRESS\n'
            #else:
                #msg = 'name=%r\n\n%s' % (name, self.code_information())
                #raise RuntimeError(msg)

        if self.is_stress:
            stress_strain = '                         C O M P L E X   S T R E S S E S   I N   B E A M   E L E M E N T S   ( C B E A M )'
        else:
            stress_strain = '                         C O M P L E X    S T R A I N S    I N   B E A M   E L E M E N T S   ( C B E A M )'

        msg_temp = [
            stress_strain,
            mag_phase,
            line1,
            line2,
        ]
        if self.is_sort1:
            self._write_sort1_as_sort1(f06_file, name, header, page_stamp, msg_temp, page_num,
                                       is_mag_phase=is_mag_phase)
        else:
            raise NotImplementedError()
        return page_num - 1

    def _write_sort1_as_sort1(self, f06_file, name, header, page_stamp, msg_temp, page_num,
                              is_mag_phase=False):
        ntimes = self.data.shape[0]

        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (name, dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write('\n'.join(msg))

            sxc = self.data[itime, :, 0]
            sxd = self.data[itime, :, 1]
            sxe = self.data[itime, :, 2]
            sxf = self.data[itime, :, 3]
            #[sxc, sxd, sxe, sxf]

            eids = self.element_node[:, 0]
            nids = self.element_node[:, 1]
            ueids = np.unique(eids)
            i_sd_zero = np.searchsorted(eids, ueids, side='left')
            i_sd_one = np.searchsorted(eids, ueids, side='right')
            i_wrong_index = np.where(i_sd_one == len(eids))
            i_sd_one[i_wrong_index] = i_sd_one[i_wrong_index] - 1
            i_sd_one -= 1
            inid = np.union1d(i_sd_zero, i_sd_one)
            #print(inid)
            #print(eids)
            nid_type = np.zeros(len(eids), dtype='bool')
            i_sd_zero_all = np.zeros(len(eids), dtype='bool')
            nid_type[inid] = 1
            i_sd_zero_all[i_sd_zero] = 1
            #print(nid_type)
            for eid, nid, i_sd_zeroi, nid_typei, sd, sxc, sxd, sxe, sxf in zip(eids, nids, i_sd_zero_all, nid_type, self.sd, sxc, sxd, sxe, sxf):
                vals = (sxc, sxd, sxe, sxf)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (sxcr, sxdr, sxer, sxfr,
                 sxci, sxdi, sxei, sxfi) = vals2

                if nid_typei or 1:
                    if i_sd_zeroi:
                        f06_file.write('0  %8s\n' % eid)
                        f06_file.write(
                            '0%8s  %8i   %.3f     %-13s    %-13s    %-13s    %s\n'
                            ' %8s   %8s  %5s     %-13s    %-13s    %-13s    %s\n' % (
                                '', nid, sd,
                                sxcr, sxdr, sxer, sxfr,
                                '', '', '',
                                sxci, sxdi, sxei, sxfi))
                    else:
                        f06_file.write(
                            '0%8s  %8i   %.3f     %-13s    %-13s    %-13s    %s\n'
                            ' %8s   %8s  %5s     %-13s    %-13s    %-13s    %s\n' % (
                                '', nid, sd,
                                sxcr, sxdr, sxer, sxfr,
                                '', '', '',
                                sxci, sxdi, sxei, sxfi))
                else:
                    f06_file.write(
                        ' %8s  %8s   %.3f     %-13s    %-13s    %-13s    %s\n'
                        ' %8s   %8s  %5s     %-13s    %-13s    %-13s    %s\n' % (
                            '', '', sd,
                            sxcr, sxdr, sxer, sxfr,
                            '', '', '',
                            sxci, sxdi, sxei, sxfi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        #long_form = False
        #if nids.min() == 0:
            #long_form = True

        eids_device = eids * 10 + self.device_code
        ueids = np.unique(eids)
        #ieid = np.searchsorted(eids, ueids)
        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = len(ueids)

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 9f')
            struct2 = Struct(endian + b'i 9f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write(f'nelements={nelements:d}\n')

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        ueids = np.unique(eids)
        neids = len(eids)

        #print('eids', eids)
        #print('nids', nids)
        #print('sd', self.sd)
        # C:\MSC.Software\msc_nastran_runs\cc145.op2
        #inid  0        2  3    4
        #eids [201 201 201 202 202]
        #nids [60000     0 60001 60001 60002]
        #sd [0.  0.5 1.  0.  1. ]
        i_sd_zero = np.searchsorted(eids, ueids, side='left')  # first location of x/xb=0.0 (sd)
        i_sd_one = np.searchsorted(eids, ueids, side='right')  # location of x/xb=1.0 (sd)

        #print('i_sd_zero', i_sd_zero)
        #print('i_sd_one', i_sd_one)
        #i_sd_zero [0 3]
        #i_sd_one [3 5]

        # we want to make the union [0, 3] between i_sd_zero and i_sd_one
        # but that [5] index is bad (it's the length of the array); so let's get rid of it

        i_sd_one -= 1
        #print('i_sd_one', i_sd_one)
        inid = np.union1d(i_sd_zero, i_sd_one)  # the indices of nid 1/2

        #print('inid', inid)
        #inid [0 2 3]

        is_nid = np.zeros(neids, dtype='bool')
        i_sd_zero_all = np.zeros(neids, dtype='bool')
        is_nid[inid] = 1
        i_sd_zero_all[i_sd_zero] = 1

        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            sxc = self.data[itime, :, 0]
            sxd = self.data[itime, :, 1]
            sxe = self.data[itime, :, 2]
            sxf = self.data[itime, :, 3]
            #[sxc, sxd, sxe, sxf]

            nwide = 0
            icount = 0
            ielement = 0
            #print('eids_device =', eids_device)
            for eid, nid, i_sd_zeroi, is_nidi, sd, sxc, sxd, sxe, sxf in zip(eids, nids, i_sd_zero_all, is_nid, self.sd, sxc, sxd, sxe, sxf):
                if icount == 0:
                    # write eid and node 1
                    # xxb = 0.0
                    eid_device = eids_device[ielement]
                    nid = nids[ielement]
                    data = [eid_device, nid, sd.real,
                            sxc.real, sxd.real, sxe.real, sxf.real,
                            sxc.imag, sxd.imag, sxe.imag, sxf.imag,] # 10
                    #print(data)
                    op2_file.write(struct1.pack(*data))
                    ielement += 1
                    icount = 1
                elif nid > 0 and icount > 0:
                    # final line (sd=1.0)
                    #
                    # 11 total nodes, with 1, 11 getting an nid; the other 9 being
                    # xxb sections
                    data = [0, 0.,
                            0., 0., 0., 0.,
                            0., 0., 0., 0.,]
                    #print('***adding %s\n' % (10-icount))
                    for unused_i in range(10 - icount):
                        op2_file.write(struct2.pack(*data))
                        nwide += len(data)

                    eid_device2 = eids_device[ielement]
                    assert eid_device == eid_device2
                    nid = nids[ielement]
                    data = [nid, sd.real,
                            sxc.real, sxd.real, sxe.real, sxf.real,
                            sxc.imag, sxd.imag, sxe.imag, sxf.imag,] # 10
                    #print(data)
                    op2_file.write(struct2.pack(*data))
                    ielement += 1
                    icount = 0
                else:
                    # intermediate stations
                    data = [0, sd.real,
                            sxc.real, sxd.real, sxe.real, sxf.real,
                            sxc.imag, sxd.imag, sxe.imag, sxf.imag,]  # 10
                    #print(data)
                    op2_file.write(struct2.pack(*data))
                    #raise RuntimeError(f'CBEAM OES op2 writer; nid={nid} icount={icount}')
                    #op2_file.write(struct2.pack(*data))
                    ielement += 1
                    icount += 1

                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                nwide += len(data)

            assert ntotal == nwide, 'ntotal=%s nwide=%s' % (ntotal, nwide)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

class ComplexBeamStressArray(ComplexBeamArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['sxc', 'sxd', 'sxe', 'sxf']
        return headers

class ComplexBeamStrainArray(ComplexBeamArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['exc', 'exd', 'exe', 'exf']
        return headers
