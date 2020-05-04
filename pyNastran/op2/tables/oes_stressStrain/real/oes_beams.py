from itertools import count
from typing import List

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class RealBeamArray(OES_Object):
    """
    common class used by:
     - RealBeamStressArray
     - RealBeamStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, unused_dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        self.data = None
        self.element_node = None
        self.xxb = None

        #if not is_sort1:
            #raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealBeamArray"""
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 2:
            nnodes_per_element = 10
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        self.ntotal = self.nelements  #* 2  # for A/B
        #self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes,
            #self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        _times = zeros(self.ntimes, dtype=dtype)
        element_node = zeros((self.ntotal, 2), dtype=idtype)

        # sxc, sxd, sxe, sxf
        # smax, smin, MSt, MSc
        xxb = zeros(self.ntotal, dtype=fdtype)
        data = zeros((self.ntimes, self.ntotal, 8), dtype=fdtype)

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_node = group.create_dataset('element_node', data=element_node)
            self.xxb = group.create_dataset('xxb', data=xxb)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_node = element_node
            self.xxb = xxb
            self.data = data

    def finalize(self):
        sd = self.data[0, :, 0].real
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)
        #self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]
        self.xxb = self.xxb[i]

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, self.element_node, self.data)
            #data_frame = pd.Panel(self.data, items=column_values,
                                  #major_axis=element_node, minor_axis=headers).to_frame()
            #data_frame.columns.names = column_names
            #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            # Static            sxc  sxd  sxe  sxf  smax  smin    MS_tension  MS_compression
            # ElementID NodeID
            # 12        22      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            #           26      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            index = pd.MultiIndex.from_arrays(self.element_node.T, names=['ElementID', 'NodeID'])
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
            data_frame.columns.names = ['Static']
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
                    for ieid, (eid, unused_nid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (axial_stress1, equiv_stress1, total_strain1, eff_plastic_creep_strain1,
                         eff_creep_strain1, linear_torsional_stress1) = t1
                        (axial_stress2, equiv_stress2, total_strain2, eff_plastic_creep_strain2,
                         eff_creep_strain2, linear_torsional_stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                axial_stress1, equiv_stress1, total_strain1,
                                eff_plastic_creep_strain1, eff_creep_strain1,
                                linear_torsional_stress1,

                                axial_stress2, equiv_stress2, total_strain2,
                                eff_plastic_creep_strain2, eff_creep_strain2,
                                linear_torsional_stress2)
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

    #def add_new_eid(self, dt, eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc):
        #self.add_new_eid_sort1(dt, eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc)

    def add_new_eid_sort1(self, dt, eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc):
        assert isinstance(eid, integer_types), eid
        assert eid >= 0, eid
        self._times[self.itime] = dt
        self.element_node[self.itotal] = [eid, grid]
        self.xxb[self.itotal] = sd
        self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                                                 smax, smin, mst, msc]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, unused_dt, eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc):
        """unvectorized method for adding SORT1 transient data"""
        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out

        self.element_node[self.itotal, :] = [eid, grid]
        self.xxb[self.itotal] = sd
        self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                                                 smax, smin, mst, msc]
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
        nnodes = self.nnodes
        ntotal = self.ntotal
        #nlayers = 2
        nelements = self.ntotal // self.nnodes  # // 2

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, nnodes, ntotal))
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  xxb.shape = %s\n' % str(self.xxb.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1,
                  is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        xxbs = self.xxb
        assert len(eids) == len(nids)
        assert len(eids) == len(xxbs)
        #print('CBEAM ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            sxcs = self.data[itime, :, 0]
            sxds = self.data[itime, :, 1]
            sxes = self.data[itime, :, 2]
            sxfs = self.data[itime, :, 3]
            smaxs = self.data[itime, :, 4]
            smins = self.data[itime, :, 5]
            smts = self.data[itime, :, 6]
            smcs = self.data[itime, :, 7]
            assert len(eids) == len(sxcs)

            eid_old = None
            xxb_old = None
            for (eid, nid, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc) in zip(
                    eids, nids, xxbs, sxcs, sxds, sxes, sxfs, smaxs, smins, smts, smcs):
                if eid != eid_old:
                    f06_file.write('0  %8i\n' % eid)
                if xxb == xxb_old:
                    continue
                # #if eid != eid_old and xxb != xxb_old:
                    #continue
                vals = [sxc, sxd, sxe, sxf, smax, smin, smt, smc]
                vals2 = write_floats_13e(vals)
                [sxc, sxd, sxe, sxf, smax, smin, smt, smc] = vals2
                f06_file.write('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (
                    nid, xxb, sxc, sxd, sxe, sxf,
                    smax, smin, smt, smc.strip()))
                eid_old = eid
                xxb_old = xxb

            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        xxbs = self.xxb
        #print(xxbs)

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

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 9f')
            struct2 = Struct(endian + b'i 9f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
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

            sxcs = self.data[itime, :, 0]
            sxds = self.data[itime, :, 1]
            sxes = self.data[itime, :, 2]
            sxfs = self.data[itime, :, 3]
            smaxs = self.data[itime, :, 4]
            smins = self.data[itime, :, 5]
            smts = self.data[itime, :, 6]
            smcs = self.data[itime, :, 7]

            #eid_old = None
            #xxb_old = None
            icount = 0
            nwide = 0
            ielement = 0
            #print('------------')
            #print(self.element_node.shape, self.data.shape)
            for (unused_i, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc) in zip(
                    count(), xxbs, sxcs, sxds, sxes, sxfs, smaxs, smins, smts, smcs):

                if icount == 0:
                    eid_device = eids_device[ielement]
                    nid = nids[ielement]
                    data = [eid_device, nid, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc] # 11
                    op2.write(struct1.pack(*data))
                    ielement += 1
                    icount = 1
                elif xxb == 1.0:
                    # 11 total nodes, with 1, 11 getting an nid; the other 9 being
                    # xxb sections
                    data = [0, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
                    #print('***adding %s\n' % (10-icount))
                    for unused_j in range(10 - icount):
                        op2.write(struct2.pack(*data))
                        nwide += len(data)

                    eid_device2 = eids_device[ielement]
                    assert eid_device == eid_device2
                    nid = nids[ielement]
                    data = [nid, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc] # 11
                    op2.write(struct2.pack(*data))
                    ielement += 1
                    icount = 0
                else:
                    data = [0, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc]  # 10
                    op2.write(struct2.pack(*data))
                    icount += 1

                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                nwide += len(data)

            assert ntotal == nwide, 'ntotal=%s nwide=%s' % (ntotal, nwide)

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealNonlinearBeamArray(OES_Object):
    """tested by elements/loadstep_elements.op2"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealNonlinearBeamArray"""
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 94:
            nnodes_per_element = 10
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        self.ntotal = self.nelements  #* 2  # for A/B
        #self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes,
            #self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 3), dtype=idtype)

        #gridA, CA, long_CA, eqS_CA, tE_CA, eps_CA, ecs_CA,
        #       DA, long_DA, eqS_DA, tE_DA, eps_DA, ecs_DA,
        #       EA, long_EA, eqS_EA, tE_EA, eps_EA, ecs_EA,
        #       FA, long_FA, eqS_FA, tE_FA, eps_FA, ecs_FA,
        #gridB, CB, long_CB, eqS_CB, tE_CB, eps_CB, ecs_CB,
        #       DB, long_DB, eqS_DB, tE_DB, eps_DB, ecs_DB,
        #       EB, long_EB, eqS_EB, tE_EB, eps_EB, ecs_EB,
        #       FB, long_FB, eqS_FB, tE_FB, eps_FB, ecs_FB,
        #self.xxb = zeros(self.ntotal, dtype='float32')
        self.data = zeros((self.ntimes, self.ntotal, 5), dtype=fdtype)

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        #nlayers = 2
        nelements = self.ntotal // self.nnodes  # // 2

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, nnodes, ntotal))
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def add_new_eid_sort1(self, dt, eid, grid_a,
                          unused_ca, long_ca, eqs_ca, te_ca, eps_ca, ecs_ca,
                          unused_da, long_da, eqs_da, te_da, eps_da, ecs_da,
                          unused_ea, long_ea, eqs_ea, te_ea, eps_ea, ecs_ea,
                          unused_fa, long_fa, eqs_fa, te_fa, eps_fa, ecs_fa,
                          grid_b,
                          unused_cb, long_cb, eqs_cb, te_cb, eps_cb, ecs_cb,
                          unused_db, long_db, eqs_db, te_db, eps_db, ecs_db,
                          unused_eb, long_eb, eqs_eb, te_eb, eps_eb, ecs_eb,
                          unused_fb, long_fb, eqs_fb, te_fb, eps_fb, ecs_fb):
        #assert eid >= 0, eid
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        #(grid_a,
         #unused_ca, long_ca, eqs_ca, te_ca, eps_ca, ecs_ca,
         #unused_da, long_da, eqs_da, te_da, eps_da, ecs_da,
         #unused_ea, long_ea, eqs_ea, te_ea, eps_ea, ecs_ea,
         #unused_fa, long_fa, eqs_fa, te_fa, eps_fa, ecs_fa,
         #grid_b,
         #unused_cb, long_cb, eqs_cb, te_cb, eps_cb, ecs_cb,
         #unused_db, long_db, eqs_db, te_db, eps_db, ecs_db,
         #unused_eb, long_eb, eqs_eb, te_eb, eps_eb, ecs_eb,
         #unused_fb, long_fb, eqs_fb, te_fb, eps_fb, ecs_fb,) = out[1:]

        self.element_node[self.itotal] = [eid, grid_a, 0]
        self.element_node[self.itotal + 1] = [eid, grid_a, 1]
        self.element_node[self.itotal + 2] = [eid, grid_a, 2]
        self.element_node[self.itotal + 3] = [eid, grid_a, 3]
        self.element_node[self.itotal + 4] = [eid, grid_b, 4]
        self.element_node[self.itotal + 5] = [eid, grid_b, 5]
        self.element_node[self.itotal + 6] = [eid, grid_b, 6]
        self.element_node[self.itotal + 7] = [eid, grid_b, 7]

        self.data[self.itime, self.itotal, :] = [long_ca, eqs_ca, te_ca, eps_ca, ecs_ca]
        self.data[self.itime, self.itotal + 1, :] = [long_da, eqs_da, te_da, eps_da, ecs_da]
        self.data[self.itime, self.itotal + 2, :] = [long_ea, eqs_ea, te_ea, eps_ea, ecs_ea]
        self.data[self.itime, self.itotal + 3, :] = [long_fa, eqs_fa, te_fa, eps_fa, ecs_fa]
        self.data[self.itime, self.itotal + 4, :] = [long_cb, eqs_cb, te_cb, eps_cb, ecs_cb]
        self.data[self.itime, self.itotal + 5, :] = [long_db, eqs_db, te_db, eps_db, ecs_db]
        self.data[self.itime, self.itotal + 6, :] = [long_eb, eqs_eb, te_eb, eps_eb, ecs_eb]
        self.data[self.itime, self.itotal + 7, :] = [long_fb, eqs_fb, te_fb, eps_fb, ecs_fb]
        self.itotal += 8
        #print('CBEAM-94:  out=%s' % str(out))
        self.ielement += 1

    #def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        #return ind

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
                    for ieid, (eid, nid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (force1, stress1) = t1
                        (force2, stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s)\n  (%s, %s)\n' % (
                                eid,
                                force1, stress1,
                                force2, stress2)
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

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1,
                  is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        locs = self.element_node[:, 2]
        #xxbs = self.xxb
        #print('CBEAM ntimes=%s ntotal=%s' % (ntimes, ntotal))
        loc_map = ['C', 'D', 'E', 'F',
                   'C', 'D', 'E', 'F',]
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            longs = self.data[itime, :, 0]
            eqss = self.data[itime, :, 1]
            tes = self.data[itime, :, 2]
            epss = self.data[itime, :, 3]
            ecss = self.data[itime, :, 4]

            #msg = ['                        N O N L I N E A R   S T R E S S E S   I N   B E A M   E L E M E N T S     ( C B E A M )\n',
            #' \n',
            #'          ELEMENT    GRID     POINT        STRESS          EQUIVALENT        TOTAL STRAIN      EFF. STRAIN       EFF. CREEP\n',
            #'             ID       ID                                     STRESS                          PLASTIC/NLELAST       STRAIN\n',]
            #'0               1         1     C        1.738817E+03      1.738817E+03      5.796055E-05      0.0               0.0\n',
            #'                                D        1.229523E+03      1.229523E+03      4.098411E-05      0.0               0.0\n',
            for (eid, nid, loc, longi, eqs, te, eps, ecs) in zip(
                eids, nids, locs, longs, eqss, tes, epss, ecss):

                vals = [longi, eqs, te, eps, ecs]
                vals2 = write_floats_13e(vals)
                [longi, eqs, te, eps, ecs] = vals2
                if loc == 0:
                    f06_file.write('0  %14i  %8i  %4s       %13s     %13s     %13s %13s %s\n' % (
                        eid, nid, 'C', longi, eqs, te, eps, ecs.rstrip()))
                elif loc == 4:
                    f06_file.write('   %14s  %8i  %4s       %13s     %13s     %13s %13s %s\n' % (
                        '', nid, 'C', longi, eqs, te, eps, ecs.rstrip()))
                else:
                    loci = loc_map[loc]
                    f06_file.write('   %14s  %8s  %4s       %13s     %13s     %13s %13s %s\n' % (
                        '', '', loci, longi, eqs, te, eps, ecs.rstrip()))
            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num


class RealBeamStressArray(RealBeamArray, StressObject):
    """Real CBEAM Stress"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
            'smax', 'smin', 'MS_tension', 'MS_compression'
        ]
        return headers

    def _get_msgs(self):
        if self.element_type == 2:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
            '                    STAT DIST/\n',
            '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        return msg


class RealBeamStrainArray(RealBeamArray, StrainObject):
    """Real CBEAM Strain"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
            'smax', 'smin', 'MS_tension', 'MS_compression'
        ]
        return headers

    def _get_msgs(self):
        if self.element_type == 2:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
            '                    STAT DIST/\n',
            '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        return msg


class RealNonlinearBeamStressArray(RealNonlinearBeamArray, StressObject):
    """tested by elements/loadstep_elements.op2"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealNonlinearBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = [
            'longitudinal_stress', 'equivalent_stress',
            'total_strain', 'equivalent_plastic_strain', 'equivalent_creep_strain'
        ]
        return headers

    def _get_msgs(self):
        if self.element_type == 94:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = ['                        N O N L I N E A R   S T R E S S E S   I N   B E A M   E L E M E N T S     ( C B E A M )\n',
               ' \n',
               '          ELEMENT    GRID     POINT        STRESS          EQUIVALENT        TOTAL STRAIN      EFF. STRAIN       EFF. CREEP\n',
               '             ID       ID                                     STRESS                          PLASTIC/NLELAST       STRAIN\n',]
        #'0               1         1     C        1.738817E+03      1.738817E+03      5.796055E-05      0.0               0.0\n',
        #'                                D        1.229523E+03      1.229523E+03      4.098411E-05      0.0               0.0\n',

        #msg = ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
        #                '                    STAT DIST/\n',
        #                '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        return msg
