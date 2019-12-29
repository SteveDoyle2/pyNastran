from typing import List
import numpy as np
from numpy import zeros, searchsorted, ravel

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header #, get_key0


class RealTriaxArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

        self.itime = 0
        self.itotal = 0
        self.element_node = None

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> int:
        if self.element_type == 53:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)
        return nnodes_per_element

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)
        #return headers

    def build(self):
        """sizes the vectorized attributes of the RealTriaxArray"""
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        unused_nnodes_per_element = self.nnodes_per_element
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
        element_node = zeros((self.ntotal, 2), dtype='int32')

        # [radial, azimuthal, axial, shear, omax, oms, ovm]
        data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_node = group.create_dataset('element_node', data=element_node)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_node = element_node
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            # LoadStep                             1.0
            # ElementID NodeID Item
            #30011     0      radial     2.000005e+02
            #          30011  azimuthal  2.000005e+02
            #          30012  axial      1.960005e+02
            #          30013  shear     -1.441057e-09
            #30012     0      omax       2.000005e+02
            #...                                  ...
            #30021     30023  axial      1.934256e+02
            #30022     0      shear      0.000000e+00
            #          30021  omax       1.973730e+02
            #          30023  oms        1.973730e+00
            #          30024  ovm        3.947461e+00
            column_names, column_values = self._build_dataframe_transient_header()
            names = ['ElementID', 'NodeID', 'Item']
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, element_node, self.data, from_tuples=False, from_array=True,
                names=names,
            )

            #column_names, column_values = self._build_dataframe_transient_header()
            #data_frame = self._build_pandas_transient_element_node(
                #column_values, column_names,
                #headers, self.element_node, self.data)
        else:
            #                    radial  azimuthal     axial     shear      omax       oms       ovm
            #ElementID NodeID
            #5301      0      -0.018626  -0.090677 -0.007052  0.010027 -0.090677  0.044707  0.080379
            #          5301   -0.000184  -0.050802 -0.015103  0.025020 -0.050802  0.034634  0.062511
            #          5303   -0.040201  -0.117615  0.015897  0.010630 -0.117615  0.067729  0.117565
            #          5305   -0.021945  -0.166931 -0.028821 -0.005979 -0.166931  0.074223  0.142052
            #5311      0      -0.016238  -0.072509 -0.006677  0.010041 -0.072509  0.036086  0.064018
            #          5311   -0.006616  -0.089126 -0.007863  0.012473 -0.089126  0.047187  0.084695
            #          5313   -0.032047  -0.038660  0.001105  0.011260 -0.038660  0.021614  0.041742
            #          5315   -0.017285  -0.129091 -0.017932  0.002921 -0.129091  0.057211  0.111599
            df1 = pd.DataFrame(element_node).T
            df1.columns = ['ElementID', 'NodeID']
            df2 = pd.DataFrame(self.data[0], columns=headers)
            data_frame = df1.join(df2).set_index(['ElementID', 'NodeID'])
            #self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
            #self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid1, Eid2\n' % str(self.code_information())
            for (eid1, nid1), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s) (%s, %s)\n' % (eid1, nid1, eid2, nid2)
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
                        (radial1, azimuthal1, axial1, shear1, omax1, oms1, ovm1) = t1
                        (radial2, azimuthal2, axial2, shear2, omax2, oms2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                radial1, azimuthal1, axial1, shear1, omax1, oms1, ovm1,
                                radial2, azimuthal2, axial2, shear2, omax2, oms2, ovm2)
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

    def add_sort1(self, dt, eid, nid, radial, azimuthal, axial, shear, omax, oms, ovm):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [radial, azimuthal, axial, shear, omax, oms, ovm]
        self.itotal += 1
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.ntotal
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        (ntimes, unused_ntotal) = self.data.shape[:2]
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #[radial, azimuthal, axial, shear, omax, oms, ovm]
            radial = self.data[itime, :, 0]
            azimuthal = self.data[itime, :, 1]
            axial = self.data[itime, :, 2]
            shear = self.data[itime, :, 3]
            omax = self.data[itime, :, 4]
            oms = self.data[itime, :, 5]
            ovm = self.data[itime, :, 6]

            for (eid, nid, radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi) in zip(
                    eids, nids, radial, azimuthal, axial, shear, omax, oms, ovm):

                vals = [radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi]
                vals2 = write_floats_13e(vals)
                [radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi] = vals2
                f06_file.write(
                    '0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    % (eid, nid, radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num


class RealTriaxStressArray(RealTriaxArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTriaxArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['radial', 'azimuthal', 'axial', 'shear', 'omax', 'oms', 'ovm']
        return headers

    def _get_msgs(self) -> List[str]:
        if self.element_type == 53:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
            '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
            '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',
            #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
            #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
        ]
        return msg

class RealTriaxStrainArray(RealTriaxArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTriaxArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['radial', 'azimuthal', 'axial', 'shear', 'omax', 'oms', 'ovm']
        return headers

    def _get_msgs(self) -> List[str]:
        if self.element_type == 53:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
            '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
            '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',
            #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
            #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
        ]
        return msg
