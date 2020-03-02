from itertools import count
from typing import List

import numpy as np
from numpy import zeros, searchsorted
from pyNastran.utils.numpy_utils import integer_types

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class RealBar10NodesArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealBar10NodesArray"""
        #print("self.ielement =", self.ielement)
         #print('RealBar10NodesArray isubcase=%s ntimes=%s nelements=%s ntotal=%s' % (
            # self.isubcase, self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 100:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        #self.ntotal = self.nelements  #* 2  # for A/B
        #self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'

        _times = zeros(self.ntimes, dtype=dtype)
        element = zeros(self.ntotal, dtype='int32')

        #[sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
        data = zeros((self.ntimes, self.ntotal, 9), dtype='float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element = element
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Time                      0.0           0.5           1.0
            #ElementID Item
            #11        sd     0.000000e+00  0.000000e+00  0.000000e+00
            #          sxc    0.000000e+00  1.000876e-01  2.200609e-01
            #          sxd    0.000000e+00 -1.000876e-01 -2.200609e-01
            #          sxe    0.000000e+00  1.000876e-01  2.200609e-01
            #          sxf    0.000000e+00 -1.000876e-01 -2.200609e-01
            #          axial  1.000000e-03 -5.566661e-05  4.833280e-05
            #          smax   1.000000e-03  1.000319e-01  2.201092e-01
            #          smin   1.000000e-03 -1.001432e-01 -2.200126e-01
            #          MS     1.401298e-45  1.401298e-45  1.401298e-45
            #          sd     5.000000e-01  5.000000e-01  5.000000e-01
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            # >=25.0
            #Static      sd  sxc  sxd  sxe  sxf     axial      smax      smin            MS
            #ElementID
            #10         0.0  0.0  0.0  0.0  0.0  0.003300  0.003300  0.003300  1.401298e-45
            #10         1.0  0.0  0.0  0.0  0.0 -0.000033 -0.000033 -0.000033  1.401298e-45
            #
            # <=24.2
            #ElementID Item
            #10        sd     0.000000e+00
                #sxc    0.000000e+00
                #sxd    0.000000e+00
                #sxe    0.000000e+00
                #sxf    0.000000e+00
                #axial  3.300000e-03
                #smax   3.300000e-03
                #smin   3.300000e-03
                #MS     1.401298e-45
                #sd     1.000000e+00
                #sxc    0.000000e+00
                #sxd    0.000000e+00
                #sxe    0.000000e+00
                #sxf    0.000000e+00
                #axial -3.333333e-05
                #smax  -3.333333e-05
                #smin  -3.333333e-05
                #MS     1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
            #data_frame = pd.Panel(self.data, major_axis=self.element,
                                  #minor_axis=headers).to_frame()
            #data_frame.columns.names = ['Static']
            #data_frame.index.names = ['ElementID', 'Item']
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
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1) = t1
                        (axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1,
                                axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2)
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

    def add_new_eid_sort1(self, etype, dt, eid,
                          sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS):
        self._times[self.itime] = dt
        #print('isubcase=%s itotal=%s ieid=%s eid=%s' % (self.isubcase, self.itotal, self.ielement, eid))
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
        self.itotal += 1
        self.ielement += 1

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
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ##ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        #print('CBAR ntimes=%s ntotal=%s' % (ntimes, ntotal))
        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(f06_file, header, page_stamp, msg, page_num)
        else:
            raise RuntimeError()
        return page_num

    def _write_sort1_as_sort1(self, f06_file, header, page_stamp, msg, page_num):
        ntimes = self.data.shape[0]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            sd = self.data[itime, :, 0]
            sxc = self.data[itime, :, 1]
            sxd = self.data[itime, :, 2]
            sxe = self.data[itime, :, 3]
            sxf = self.data[itime, :, 4]
            axial = self.data[itime, :, 5]
            smax = self.data[itime, :, 6]
            smin = self.data[itime, :, 7]
            MS = self.data[itime, :, 8]

            for (i, eid, sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi) in zip(
                    count(), eids, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS):

                vals = [sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi]
                vals2 = write_floats_13e(vals)
                [sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi] = vals2
                f06_file.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s %s\n'
                               % (eid, sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi))

            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num

class RealBar10NodesStressArray(RealBar10NodesArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBar10NodesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        #if self.is_fiber_distance:
            #fiber_dist = 'fiber_distance'
        #else:
            #fiber_dist = 'fiber_curvature'

        #if self.is_von_mises:
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        headers = ['sd', 'sxc', 'sxd', 'sxe', 'sxf', 'axial', 'smax', 'smin', 'MS']
        return headers

    def _get_msgs(self):
        msg = [
            '                         S T R E S S   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )\n'
            '0    ELEMENT  STATION    SXC           SXD           SXE           SXF            AXIAL          S-MAX         S-MIN         M.S.-T\n'
            '       ID.     (PCT)                                                                                                         M.S.-C\n'
            #'            1   0.000   4.919032E+05 -4.348710E+05 -4.348710E+05  4.919032E+05   0.0            4.919032E+05 -4.348710E+05 \n'
        ]
        return msg


class RealBar10NodesStrainArray(RealBar10NodesArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBar10NodesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        #if self.is_fiber_distance:
            #fiber_dist = 'fiber_distance'
        #else:
            #fiber_dist = 'fiber_curvature'

        #if self.is_von_mises:
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        headers = ['sd', 'sxc', 'sxd', 'sxe', 'sxf', 'axial', 'smax', 'smin', 'MS']
        return headers

    def _get_msgs(self):
        msg = [
            '                         S T R A I N   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )\n'
            '0    ELEMENT  STATION    SXC           SXD           SXE           SXF            AXIAL          S-MAX         S-MIN         M.S.-T\n'
            '       ID.     (PCT)                                                                                                         M.S.-C\n'
            #'            1   0.000   4.919032E+05 -4.348710E+05 -4.348710E+05  4.919032E+05   0.0            4.919032E+05 -4.348710E+05 \n'
        ]
        return msg
