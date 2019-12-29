from typing import List

import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class RandomBendArray(OES_Object):
    """
    common class used by:
     - RandomBendStressArray
     - RandomBendStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, unused_dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

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
        nnodes = 1

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 69:
            nnodes_per_element = 2
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
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
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        # sxc, sxd, sxe, sxf
        self.angle = zeros(self.ntotal, dtype='float32')
        self.data = zeros((self.ntimes, self.ntotal, 4), dtype='float32')

    def finalize(self):
        sd = self.data[0, :, 0].real
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)
        #self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node,
                                       minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']

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

    #def add_new_eid_sort1(self, dt, eid, grid, angle, sxc, sxd, sxe, sxf):
        #assert isinstance(eid, integer_types), eid
        #assert eid >= 0, eid
        #self._times[self.itime] = dt
        #self.element_node[self.itotal] = [eid, grid]
        #self.angle[self.itotal] = angle
        #self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf]
        #self.itotal += 1
        #self.ielement += 1

    def add_sort1(self, dt, eid, grid, angle, sxc, sxd, sxe, sxf):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.element_node[self.itotal, :] = [eid, grid]
        self.angle[self.itotal] = angle
        self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf]
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
        msg.append('  angle.shape = %s\n' % str(self.angle.shape).replace('L', ''))
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
        angles = self.angle
        #print('CBEAM ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            sxcs = self.data[itime, :, 0]
            sxds = self.data[itime, :, 1]
            sxes = self.data[itime, :, 2]
            sxfs = self.data[itime, :, 3]

            eid_old = None
            xxb_old = None
            for (eid, nid, angle, sxc, sxd, sxe, sxf) in zip(
                eids, nids, angles, sxcs, sxds, sxes, sxfs):
                if eid != eid_old:
                    f06_file.write('0  %8i\n' % eid)
                #if xxb == xxb_old:
                    #continue
                # #if eid != eid_old and xxb != xxb_old:
                    #continue
                vals = [sxc, sxd, sxe, sxf]
                vals2 = write_floats_13e(vals)
                [sxc, sxd, sxe, sxf] = vals2
                f06_file.write('%19s   %4.3f   %12s %12s %12s %s\n' % (
                    nid, angle, sxc, sxd, sxe, sxf.strip()))
                eid_old = eid
                #xxb_old = xxb

            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num

    def _get_msgs(self):
        if self.element_type == 69:
            pass
        else:
            raise NotImplementedError(self.element_type)

        #is_stress = True
        is_strain = 'OSTR' in self.table_name
        is_stress = not is_strain
        if self.table_name in ['OESATO2', 'OSTRATO2']:
            table_header = '                                                 ( AUTO-CORRELATION FUNCTION )'
        elif self.table_name in ['OESCRM2', 'OSTRCRM2']:
            table_header = '                                               ( CUMULATIVE ROOT MEAN SQUARE )'
        elif self.table_name in ['OESPSD2', 'OSTRPSD2']:
            table_header = '                                             ( POWER SPECTRAL DENSITY FUNCTION )'
        #elif self.table_name in ['OESNO2', 'OSTRNO2']:
            #pass
        #elif self.table_name in ['OESNO2', 'OSTRNO2']:
            #pass
        else:
            raise NotImplementedError(self.table_name)

        if is_stress:
            stress_strain = '                                  S T R E S S E S   I N   B E N D   E L E M E N T S        ( C B E N D )'
        else:
            stress_strain = '                                   S T R A I N S    I N   B E N D   E L E M E N T S        ( C B E N D )'

        #'                                               ( CUMULATIVE ROOT MEAN SQUARE )'
        '                        CIRC.      LOCATION         LOCATION         LOCATION         LOCATION'
        '   FREQUENCY   GRID END  ANG.         C                D                E                F'
        #'0 2.500E+00    6901   A    0     2.253653E+02     1.017869E+01     1.981943E+02     1.699232E+01'
        #'0              6904   B    0     2.785208E+01     4.538458E+01     1.257464E+01     2.858128E+01'

        msg = [
            table_header,
            stress_strain,
            '                        CIRC.      LOCATION         LOCATION         LOCATION         LOCATION'
            '   FREQUENCY   GRID END  ANG.         C                D                E                F\n'
        ]
        return msg


class RandomBendStressArray(RandomBendArray, StressObject):
    """Random CBEND Stress"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
        ]
        return headers


class RandomBendStrainArray(RandomBendArray, StrainObject):
    """Random CBEND Strain"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
        ]
        return headers
