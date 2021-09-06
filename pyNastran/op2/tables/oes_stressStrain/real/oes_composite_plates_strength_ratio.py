from typing import List
import numpy as np
from numpy import zeros, searchsorted, unique, ravel

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, OES_Object) # StrainObject,
from pyNastran.f06.f06_formatting import write_float_12e, write_float_13e, _eigenvalue_header
from pyNastran.femutils.utils import pivot_table


class RealCompositePlateStrengthRatioArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        self.element_layer = None

        #if is_sort1:
            #if dt is not None:
                #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealCompositePlateArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 95:  # CQUAD4
            nnodes_per_element = 1
        elif self.element_type == 96:  # CQUAD8
            nnodes_per_element = 1
        elif self.element_type == 97:  # CTRIA3
            nnodes_per_element = 1
        elif self.element_type == 98:  # CTRIA6
            nnodes_per_element = 1
        #elif self.element_type == 232:  # CQUADR
            #nnodes_per_element = 1
        #elif self.element_type == 233:  # CTRIAR
            #nnodes_per_element = 1
        else:  # pragma: no cover
            msg = 'element_name=%s element_type=%s' %(self.element_name, self.element_type)
            raise NotImplementedError(msg)

        self.nnodes = nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        if self.is_sort1:
            ntimes = self.ntimes
            ntotal = self.ntotal
        else:
            ntimes = self.ntotal
            ntotal = self.ntimes

        _times = zeros(ntimes, dtype=dtype)
        element_layer = zeros((ntotal, 2), dtype=idtype)

        #[strength_ratio_ply, failure_index_bonding, strength_ratio_bonding]
        data = zeros((ntimes, ntotal, 3), dtype=fdtype)
        # [failure_theory, flag]
        failure_theory_flag = np.empty((ntotal, 2), dtype='U8')

        #if self.load_as_h5:
        #    #for key, value in sorted(self.data_code.items()):
        #        #print(key, value)
        #    group = self._get_result_group()
        #    self._times = group.create_dataset('_times', data=_times)
        #    self.element_layer = group.create_dataset('element_layer', data=element_layer)
        #    self.data = group.create_dataset('data', data=data)
        #else:
        self._times = _times
        self.element_layer = element_layer
        self.failure_theory_flag = failure_theory_flag
        self.data = data

    #def build_dataframe(self):
    #    """
    #    major-axis - the axis
    #
    #    mode   1     2   3
    #    freq  1.0   2.0 3.0
    #    T1
    #    T2
    #    T3
    #    R1
    #    R2
    #    R3
    #
    #    major_axis / top = [
    #        [1, 2, 3],
    #        [1.0, 2.0, 3.0]
    #    ]
    #    minor_axis / headers = [T1, T2, T3, R1, R2, R3]
    #    name = mode
    #    """
    #    import pandas as pd
    #
    #    headers = self.get_headers()
    #    if self.nonlinear_factor not in (None, np.nan):
    #        #Mode                                  1             2             3
    #        #Freq                       1.482246e-10  3.353940e-09  1.482246e-10
    #        #Eigenvalue                -8.673617e-19  4.440892e-16  8.673617e-19
    #        #Radians                    9.313226e-10  2.107342e-08  9.313226e-10
    #        #ElementID Layer Item
    #        #16        1     o11       -1.052490e-13  3.106268e-08  1.121784e-13
    #        #                o22        4.804592e-13  1.855033e-07 -9.785236e-13
    #        #                t12        4.436908e-14  4.873383e-09  4.387037e-15
    #        #                t1z        8.207617e-14  2.501582e-08 -1.056211e-13
    #        #                t2z       -5.918040e-14 -1.112469e-08  1.255247e-13
    #        #                angle      8.569244e+01  8.819442e+01  2.304509e-01
    #        #                major      4.838012e-13  1.856569e-07  1.121961e-13
    #        #                minor     -1.085910e-13  3.090905e-08 -9.785411e-13
    #        #                max_shear  2.961961e-13  7.737391e-08  5.453687e-13
    #        #          2     o11       -6.490381e-14  2.856533e-08  4.105937e-14
    #        # columns
    #        #[(1, 1.4822459136312394e-10, -8.673617379884035e-19, 9.313225746154785e-10)
    #         #(2, 3.353939638127037e-09, 4.440892098500626e-16, 2.1073424255447017e-08)
    #         #(3, 1.4822459136312394e-10, 8.673617379884035e-19, 9.313225746154785e-10)]
    #        column_names, column_values = self._build_dataframe_transient_header()
    #        data_frame = self._build_pandas_transient_element_node(
    #            column_values, column_names,
    #            headers, self.element_layer, self.data)
    #    else:
    #        element_layer = [self.element_layer[:, 0], self.element_layer[:, 1]]
    #        # Static                 o11        o22        t12        t1z  ...     angle       major       minor   max_shear
    #        # ElementID Layer                                              ...
    #        # 16        1     -2193.9639   1773.909  -2325.400  5.477e+02  ... -65.32178   284.30176  -326.28027   56.329102
    #        #           2     -1843.9912   1465.191  -2445.139  1.277e+03  ... -62.41302   276.41992  -314.80713   52.761230
    #        #           3     -1260.6953    952.560  -2646.621  1.451e+03  ... -56.48576   271.34375  -302.74707   68.154541
    #        #           4      -444.0792    235.137  -2926.092 -0.000e+00  ... -48.08685   284.21777  -305.24219   46.322998
    #        # 17        1     -1546.0195   4338.887  -2750.557  3.610e+02  ... -68.65797   542.58496  -263.74561   28.316406
    #        #           2     -1597.4194   4303.379  -2707.898  9.309e+02  ... -68.34154   535.37598  -265.98535   04.518066
    #        #           3     -1683.7607   4245.215  -2634.891  1.393e+03  ... -69.88499   524.96875  -268.98779   65.647705
    #        #           4     -1802.0312   4163.371  -2531.777  1.295e+03  ... -69.39493   509.74609  -273.14307   12.944336
    #        #           5     -1956.2432   4058.359  -2400.559  2.975e-13  ... -70.02080   489.06738  -279.80811   48.243652
    #        #
    #        #element_layer = self.element_layer #???
    #        index = pd.MultiIndex.from_arrays(element_layer, names=['ElementID', 'Layer'])
    #        data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
    #        data_frame.columns.names = ['Static']
    #        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.element_layer, table.element_layer):
            assert self.element_layer.shape == table.element_layer.shape, 'element_layer shape=%s table.shape=%s' % (
                self.element_layer.shape, table.element_layer.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Layer)\n'
            for (eid, layer1), (eid2, layer2) in zip(self.element_layer, table.element_layer):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, layer1, eid2, layer2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_layer):
                    (eid, layer) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    nt1 = filter_nan(t1)
                    nt2 = filter_nan(t2)
                    assert len(nt1) == len(nt2), 't1=%s t2=%s' % (t1, t2)
                    #print(nt1, nt2)

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(nt1, nt2):
                        (strength_ratio_ply1, failure_index_bonding1, strength_ratio_bonding1) = t1
                        (strength_ratio_ply2, failure_index_bonding2, strength_ratio_bonding2) = t2
                        msg += (
                            '(%s, %s)    (%s, %s, %s)'
                            '  (%s, %s, %s)\n' % (
                                eid, layer,
                                strength_ratio_ply1, failure_index_bonding1, strength_ratio_bonding1,
                                strength_ratio_ply2, failure_index_bonding2, strength_ratio_bonding2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_eid_sort1(self, etype, dt, eid: int,
                      failure_theory: str, ply_id: int,
                      strength_ratio_ply: float, failure_index_bonding: float,
                      strength_ratio_bonding: float, min_sr_bonding_fi_bonding: float,
                      flag: str):
        self._times[self.itime] = dt
        #layer = 0
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.failure_theory_flag[self.itotal, :] = [failure_theory, flag]
        #(eid_device, failure_theory, ply_id, strength_ratio_ply, failure_index_bonding,
        # strength_ratio_bonding, min_sr_bonding_fi_bonding, flag) = out
        #print(strength_ratio_ply, failure_index_bonding, strength_ratio_bonding)
        #asf
        self.data[self.itime, self.itotal, :] = [strength_ratio_ply, failure_index_bonding, strength_ratio_bonding]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, dt, eid: int,
                      failure_theory: str, ply_id: int, strength_ratio_ply: float, failure_index_bonding: str,
                      strength_ratio_bonding: float, min_sr_bonding_fi_bonding: float,
                      flag: str):
        """unvectorized method for adding SORT1 transient data"""
        assert eid is not None
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.element_layer[self.itotal, :] = [eid, ply_id]
        self.failure_theory_flag[self.itotal, :] = [failure_theory, flag]
        #(eid_device, failure_theory, ply_id, strength_ratio_ply, failure_index_bonding,
        # strength_ratio_bonding, min_sr_bonding_fi_bonding, flag) = out
        self.data[self.itime, self.itotal, :] = [strength_ratio_ply, failure_index_bonding, strength_ratio_bonding]
        self.itotal += 1

    #def add_eid_sort2(self, etype, dt, eid, layer, o11, o22, t12, t1z, t2z,
                          #angle, major, minor, ovm):
        #itime = self.itotal
        #itotal = self.itime
        #self._times[itime] = dt
        #self.element_layer[itotal, :] = [eid, layer]
        #self.data[itime, itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        #self.itotal += 1
        #self.ielement += 1

    #def add_sort2(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                  #major, minor, ovm):
        #"""unvectorized method for adding SORT2 transient data"""
        #assert eid is not None
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #itime = self.itotal
        #itotal = self.itime
        #self._times[itime] = dt
        ##self.element_layer[itotal, :] = [eid, layer]
        #self.data[self.itime, itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        #self.itotal += 1

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            msg = [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]
            return msg

        nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = len(unique(self.element_layer[:, 0]))

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, ntotal))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element_layer.shape = {self.element_layer.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        ctria3_msg, ctria6_msg, cquad4_msg, cquad8_msg = self._get_msgs()

        if self.element_type == 95:  # CQUAD4
            msg = cquad4_msg
            nnodes = 4
        #elif self.element_type == 96:  # CQUAD8
            #msg = cquad8_msg
            #nnodes = 8
        elif self.element_type == 97:  # CTRIA3
            msg = ctria3_msg
            #nnodes = 3
        #elif self.element_type == 98:  # CTRIA6
            #msg = ctria6_msg
            #nnodes = 6
        else:
            raise NotImplementedError(self.element_name)

        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_layer[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element_layer[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        """
        '          S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'
        '   ELEMENT  FAILURE    PLY  SRP-STRENGTH RATIO FOR PLY  SRB-STRENGTH RATIO FOR BONDING  STRENGTH RATIO FOR ELEMENT     FLAG'
        '     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MIN OF SRP,SRB FOR ALL PLIES'
        '        57   TSAI-WU     1      1.420811E+00      '
        '                                                                 8.824254E+05                                               '
        '                         2      1.130841E+01      '
        '                                                                 8.670931E+05                                               '
        '                         3      4.838508E+00      '
        '                                                                 8.183240E+05                                               '
        '                         4      2.312551E+00      '
        '                                                                 8.036290E+05                                             '
        '                         5      2.316453E+00      '
        '                                                                 8.183240E+05                                               '
        '                         6      4.871836E+00      '
        '                                                                 8.670931E+05                                               '
        '                         7      1.148123E+01      '
        '                                                                 8.824254E+05                                               '
        '                         8      1.443441E+00      '
        '                                                                                               1.420811E+00                 '
        '        59   TSAI-WU     1      1.069008E+00      '
        """
        #raise NotImplementedError('write_f06')
        if header is None:
            header = []
        #msg, nnodes, is_bilinear = self._get_msgs()
        #if self.is_von_mises:
            #von = 'VON'
            #mises = 'MISES'
        #else:
            #von = 'MAX'
            #mises = 'SHEAR'
        if self.element_type == 95:  # CQUAD4
            etype = '          S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'
        elif self.element_type == 96:  # CQUAD8
            etype = '          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'
        elif self.element_type == 97:  # CTRIA3
            etype = '          S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'
        elif self.element_type == 98:  # CTRIA6
            etype = '          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'
        else:
            print(''.join(self.get_stats()))
            msg = 'element_name=%s element_type=%s' % (self.element_name, self.element_type)
            raise NotImplementedError(msg)

        msg = [
            etype,
            '   ELEMENT  FAILURE    PLY  SRP-STRENGTH RATIO FOR PLY  SRB-STRENGTH RATIO FOR BONDING  STRENGTH RATIO FOR ELEMENT     FLAG\n'
            '     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MIN OF SRP,SRB FOR ALL PLIES\n'
        ]

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_layer[:, 0]
        layers = self.element_layer[:, 1]
        element_by_layer, ueids = pivot_table(layers, eids, layers)
        nlayers = np.nanmax(element_by_layer, axis=1)
        assert len(nlayers) == len(ueids)

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            # min_sr_bonding_fi_bonding
            #[strength_ratio_ply, failure_index_bonding, strength_ratio_bonding]
            strength_ratio_ply = self.data[itime, :, 0]
            failure_index_bonding = self.data[itime, :, 1]
            strength_ratio_bonding = self.data[itime, :, 2]

            i = 0
            is_nans_sr_bonding = np.isnan(strength_ratio_bonding)
            #is_nans_fi_bonding = np.isnan(failure_index_bonding)
            #print(is_nans.tolist())
            #print("nlayers =", nlayers.tolist())
            if 1:
                for eid, nlayer in zip(ueids, nlayers):
                    i0 = i
                    i1 = i + nlayer
                    #ft = failure_theory[i]
                    sr_ply = strength_ratio_ply[i]
                    sr_plys = write_float_13e(sr_ply)
                    #sr_bonding = strength_ratio_bonding[i]
                    #sr_bondings = write_float_13e(sr_bonding)

                    #fi_bondings = '             ' if is_nans_fi_bonding[i] else write_float_13e(failure_index_bonding[i])
                    sr_bondings = '             ' if is_nans_sr_bonding[i] else write_float_13e(strength_ratio_bonding[i])

                    ft = self.failure_theory_flag[i, 0]
                    flag = self.failure_theory_flag[i, 1]

                    # layer1
                    f06_file.write('  %8d  %8s %5d     %s      \n'
                                   '                                                                %s                                               \n' % (
                                       eid, ft, 1, sr_plys, sr_bondings))
                    i += 1

                    for unused_jlayer in range(1, nlayer-1):
                        ilayer = layers[i]
                        sr_ply = strength_ratio_ply[i]
                        sr_plys = write_float_13e(sr_ply)
                        #sr_bonding = strength_ratio_bonding[i]
                        #sr_bondings = write_float_13e(sr_bonding)
                        #fi_bondings = '             ' if is_nans_fi_bonding[i] else write_float_13e(failure_index_bonding[i])
                        sr_bondings = '             ' if is_nans_sr_bonding[i] else write_float_13e(strength_ratio_bonding[i])

                        f06_file.write('                     %5d     %s      \n'
                                       '                                                                %s                                               \n' % (
                                           ilayer, sr_plys, sr_bondings))
                        i += 1

                    # final
                    sr_ply = strength_ratio_ply[i]
                    sr_plys = write_float_13e(sr_ply)

                    sr_ply = strength_ratio_ply[i]
                    sr_plys = write_float_13e(sr_ply)
                    all_data = strength_ratio_ply[i0:i1]
                    min_sr_plys = nanmin13s(all_data)
                    #min_sr_ply = np.nanmin(strength_ratio_ply[i0:i1])
                    #min_sr_plys = write_float_13e(min_sr_ply)
                    #min_sr_plys = ' ?.??????E+00'
                    f06_file.write('                     %5d     %s      \n'
                                   '                                                                                              %13s                 \n' % (
                                       nlayer, sr_plys, min_sr_plys))
                    i += 1
            else:
                # loop over nlayers max and ueids
                failure_theory = self.failure_theory_flag[:, 0]
                flags = self.failure_theory_flag[:, 1]
                for eid, layer, ft, sr_ply, fi_bonding, sr_bonding, flag in zip(
                        eids, layers, failure_theory, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flags):
                    # TODO: super slow way to get the max layer id
                    irows = np.where(eids == eid)[0]
                    nlayers = len(irows)
                    #ft = 'TSAI-WU'

                    sr_plys = write_float_12e(sr_ply)
                    sr_bondings = write_float_12e(sr_bonding)
                    if layer == 1:
                        #eid=3;
                        #sr_ply=1.34007; fi_bonding='' sr_bonding=34482.3; min(sr_bonding, fi_bonding)=nan flag=''
                        #print(eid, layer, sr_plys, fi_bonding, sr_bondings)
                        f06_file.write('  %8d  %8s %5d      %s      \n'
                                       '                                                                 %s                                               \n' % (
                                           eid, ft, layer, sr_plys, sr_bondings))
                    elif layer == nlayers:
                        f06_file.write('                     %5d      %s      \n'
                                       '                                                                                               ?.??????E+00                 \n' % (
                                           layer, sr_plys))
                    else:
                        f06_file.write('                     %5d      %s      \n'
                                       '                                                                 %s                                               \n' % (
                                           layer, sr_plys, sr_bondings))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    #def write_op2(self, op2_file, op2_ascii, itable, new_result,
    #              date, is_mag_phase=False, endian='>'):
    #    """writes an OP2"""
    #    import inspect
    #    from struct import Struct, pack
    #    frame = inspect.currentframe()
    #    call_frame = inspect.getouterframes(frame, 2)
    #    op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')
    #
    #    if itable == -1:
    #        self._write_table_header(op2_file, op2_ascii, date)
    #        itable = -3
    #
    #    #print("nnodes_all =", nnodes_all)
    #    #msg.append(f'  element_node.shape = {self.element_node.shape}\n')
    #    #msg.append(f'  data.shape={self.data.shape}\n')
    #
    #    eids = self.element_layer[:, 0]
    #    layers = self.element_layer[:, 1]
    #    eids_device = eids * 10 + self.device_code
    #
    #    nelements = len(np.unique(eids))
    #    #print('nelements =', nelements)
    #    # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
    #    #ntotal = ((nnodes * 21) + 1) + (nelements * 4)
    #
    #    ntotali = self.num_wide
    #    nlayers = self.data.shape[1]
    #    ntotal = ntotali * nlayers
    #
    #    #print('shape = %s' % str(self.data.shape))
    #    #assert self.ntimes == 1, self.ntimes
    #
    #    #device_code = self.device_code
    #    op2_ascii.write(f'  ntimes = {self.ntimes}\n')
    #
    #    #fmt = '%2i %6f'
    #    #print('ntotal=%s' % (ntotal))
    #    #assert ntotal == 193, ntotal
    #
    #    #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
    #    op2_ascii.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
    #    op2_ascii.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
    #
    #    if not self.is_sort1:
    #        raise NotImplementedError('SORT2')
    #
    #    fdtype = self.data.dtype
    #    if self.size == 4:
    #        pass
    #    else:
    #        print(f'downcasting {self.class_name}...')
    #        #idtype = np.int32(1)
    #        fdtype = np.float32(1.0)
    #
    #    data_out = np.empty((nlayers, 11), dtype=fdtype)
    #    data_out[:, 0] = eids_device.view(fdtype)
    #    data_out[:, 1] = layers.view(fdtype)
    #
    #    op2_ascii.write(f'nelements={nelements:d}\n')
    #    ntimes = self.data.shape[0]
    #
    #    for itime in range(ntimes):
    #        self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)
    #
    #        # record 4
    #        #print('stress itable = %s' % itable)
    #        itable -= 1
    #        header = [4, itable, 4,
    #                  4, 1, 4,
    #                  4, 0, 4,
    #                  4, ntotal, 4,
    #                  4 * ntotal]
    #        op2_file.write(pack('%ii' % len(header), *header))
    #        op2_ascii.write('r4 [4, 0, 4]\n')
    #        op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
    #        op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')
    #
    #        #dt = self._times[itime]
    #        #header = _eigenvalue_header(self, header, itime, ntimes, dt)
    #        #f06_file.write(''.join(header + msg))
    #
    #        # [eid_device, layer, o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
    #        # [                   o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
    #        data_out[:, 2:] = self.data[itime, :, :]
    #        assert data_out.size == ntotal, f'data_out.shape={data_out.shape} size={data_out.size}; ntotal={ntotal}'
    #        op2_file.write(data_out)
    #
    #        itable -= 1
    #        header = [4 * ntotal,]
    #        op2_file.write(pack('i', *header))
    #        op2_ascii.write('footer = %s\n' % header)
    #        new_result = False
    #    return itable


class RealCompositePlateStressStrengthRatioArray(RealCompositePlateStrengthRatioArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCompositePlateStrengthRatioArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    @property
    def nnodes_per_element(self):
        return 1

    @property
    def is_stress(self):
        return True

    @property
    def is_strain(self):
        return False

    def get_headers(self) -> List[str]:
        #if self.is_von_mises:
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        headers = ['strength_ratio_ply', 'strength_ratio_bonding']
        #headers = ['o11', 'o22', 't12', 't1z', 't2z', 'angle', 'major', 'minor', ovm]
        return headers


#class RealCompositePlateStrainStrengthRatioArray(RealCompositePlateStrengthRatioArray, StrainObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealCompositePlateStrengthRatioArray.__init__(self, data_code, is_sort1, isubcase, dt)
        #StrainObject.__init__(self, data_code, isubcase)

    #@property
    #def is_stress(self) -> bool:
        #return False

    #@property
    #def is_strain(self) -> bool:
        #return True

    #def get_headers(self) -> List[str]:
        #if self.is_von_mises:
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        #headers = ['strength_ratio_ply_b', 'strength_ratio_bonding_b']
        #return headers

def nanmin13s(nparray):
    if np.any(np.isfinite(nparray)):
        return write_float_13e(np.nanmin(nparray))
    return '             '

def filter_nan(my_1d_array):
    is_not_nan = ~np.isnan(my_1d_array)
    return my_1d_array[is_not_nan]
