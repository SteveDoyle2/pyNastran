# coding: utf-8
#pylint disable=C0103
from itertools import count
import warnings
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.op2_interface.write_utils import to_column_bytes, view_dtype, view_idtype_as_fdtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object,
    oes_real_data_code, get_scode,
    set_static_case, set_modal_case, set_transient_case)
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
from pyNastran.op2.errors import SixtyFourBitError

NUM_WIDE_CENTROID = 17
NUM_WIDE_CORNER = 17
ELEMENT_NAME_TO_NUM_WIDE = {
    'CTRIA3': NUM_WIDE_CENTROID,
    'CQUAD4': NUM_WIDE_CENTROID,

    #2 + 17 * nnodes_all
    'CQUAD4-144': 2 + NUM_WIDE_CORNER*4,
    'CQUADR': 2 + NUM_WIDE_CORNER*5,
    'CQUAD8': 2 + NUM_WIDE_CORNER*5,

    'CTRIA6': 2 + NUM_WIDE_CORNER*4,
    'CTRIAR': 2 + NUM_WIDE_CORNER*4,

    'CTRIAR_LINEAR': NUM_WIDE_CENTROID,
    'CQUADR_LINEAR': NUM_WIDE_CENTROID,
}

ELEMENT_NAME_TO_ELEMENT_TYPE = {
    'CTRIA3': 74,
    'CQUAD4': 33,
    'CQUAD4-144': 144,

    'CQUADR': 82,
    'CQUAD8': 64,
    'CTRIA6': 75,
    'CTRIAR': 70,

    'CTRIAR_LINEAR': 227,
    'CQUADR_LINEAR': 228,
}
class RealPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        #if is_sort1:
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
        if self.element_type in [33, 74, 83, 227, 228]:
            nnodes_per_element = 1
        elif self.element_type == 144:
            nnodes_per_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_per_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_per_element = 5
        elif self.element_type == 70:  # CTRIAR
            nnodes_per_element = 4
        elif self.element_type == 75:  # CTRIA6
            nnodes_per_element = 4
        else:
            raise NotImplementedError(f'name={self.element_name!r} type={self.element_type}')
        return nnodes_per_element

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            return False
        elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            return True
        else:
            raise NotImplementedError(f'name={self.element_name} type={self.element_type}')

    def build(self):
        """sizes the vectorized attributes of the RealPlateArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        #nnodes = 2
        #ntotal = 99164
        # 99164 / 2 = 49582
        # nelements = 49582
        # nnodes = 49582 * 2 = 99164

        #self.names = []

        #factor = self.size // 4
        nnodes_per_element = self.nnodes_per_element
        #print(self.code_information())
        #print('nnodes_per_element =', nnodes_per_element)
        nlayers_per_element = 2 * nnodes_per_element

        #print('nnodes_per_element[%s, %s] = %s' % (
            #self.isubcase, self.element_type, nnodes_per_element))
        self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.nelements //= self.ntimes
        #self.ntotal //= factor
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes,
            #self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        if self.is_sort1:
            ntimes = self.ntimes
            nlayers = self.nelements
        else:
            # NUMBER OF CQUAD4   ELEMENTS =      956
            # NUMBER OF CTRIA3   ELEMENTS =       27
            #***nelements=956 nlayers_per_element=2 ntimes=201 nlayers=1912
            #***nelements=27  nlayers_per_element=2 ntimes=201

            #print(self.ntimes, self.nelements, self.ntotal, self._ntotals)
            nelements = self.ntimes # good
            ntimes = self._ntotals[0] // nlayers_per_element
            nlayers = nelements * nlayers_per_element
            #print(f'***nelements={nelements} nlayers_per_element={nlayers_per_element} ntimes={ntimes} -> nlayers={nlayers}')
            #nelements = self._ntotals[0]  # good
            #nlayers += 1
            #ntimes = nlayers // nelements
            assert nlayers % nelements == 0
            #print('***', self.element_name, nlayers)
            #assert nelements == 4, self.ntimes
            #nelements = 4
            #nelements = = self.ntimes // 2
            #print(f'ntimes={ntimes} nelements={nelements} nlayers={nlayers}; '
                  #f'nlayers_per_element={nlayers_per_element}')
            #bbb
            #assert ntimes == 1, ntimes
            #print(self.code_information())

        if self.analysis_code == 1:
            #ntimes = 1
            if ntimes != 1:
                # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\acc002.op2
                warnings.warn(f'ntimes != 1; {self.element_name}-{self.element_type}\n'
                              f'ntimes={ntimes} _ntotals={self._ntotals} '
                              f'sort_method={self.sort_method} nlayers_per_element={nlayers_per_element} nlayers={nlayers}')

        assert nlayers >= 2, self.code_information()
        _times = np.zeros(ntimes, dtype=dtype)
        element_node = np.zeros((nlayers, 2), dtype=idtype)

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        data = np.zeros((ntimes, nlayers, 8), dtype=fdtype)
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
        #print(self.element_node.shape, self.data.shape)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        nelements = self.element_node.shape[0] // 2
        if self.is_fiber_distance:
            fiber_distance = ['Top', 'Bottom'] * nelements
        else:
            fiber_distance = ['Mean', 'Curvature'] * nelements
        fd = np.array(fiber_distance, dtype='unicode')

        node = pd.Series(data=self.element_node[:, 1])
        node.replace(to_replace=0, value='CEN', inplace=True)
        element_node = [
            self.element_node[:, 0],
            node,
            fd,
        ]

        if self.nonlinear_factor not in (None, np.nan):
            # Mode                                                 1             2             3
            # Freq                                      1.482246e-10  3.353940e-09  1.482246e-10
            # Eigenvalue                               -8.673617e-19  4.440892e-16  8.673617e-19
            # Radians                                   9.313226e-10  2.107342e-08  9.313226e-10
            # ElementID NodeID Location Item
            # 8         0      Top      fiber_distance -1.250000e-01 -1.250000e-01 -1.250000e-01
            #                           oxx             7.092928e-12 -3.259632e-06 -9.558293e-12
            #                           oyy             3.716007e-12 -2.195630e-06 -5.435632e-12
            #                           txy            -7.749725e-14  1.438695e-07 -6.269848e-13
            #                           angle          -1.313964e+00  8.243371e+01 -8.154103e+01
            #                           omax            7.094705e-12 -2.176520e-06 -5.342388e-12
            #                           omin            3.714229e-12 -3.278742e-06 -9.651537e-12
            #                           von_mises       6.146461e-12  2.889834e-06  8.374427e-12
            #                  Bottom   fiber_distance  1.250000e-01  1.250000e-01  1.250000e-01
            #                           oxx            -7.530338e-12  2.134777e-06  1.063986e-11
            #                           oyy            -4.434658e-12 -9.347183e-07  6.212209e-12
            #                           txy             2.291380e-12 -5.399188e-07 -4.161393e-12
            #                           angle           6.201962e+01 -9.690845e+00 -3.099370e+01
            #                           omax           -3.217317e-12  2.226978e-06  1.313966e-11
            #                           omin           -8.747680e-12 -1.026920e-06  3.712415e-12
            #                           von_mises       7.663484e-12  2.881133e-06  1.173255e-11
            # 9         0      Top      fiber_distance -1.250000e-01 -1.250000e-01 -1.250000e-01
            #
            #LoadStep                                         1.0
            #ElementID NodeID Location Item
            #2001      CEN    Top      fiber_distance   -0.635000
            #                 Bottom   oxx              26.197712
            #2007      CEN    Top      oyy              65.378319
            #                 Bottom   txy             -28.221191
            #2008      CEN    Top      angle           -62.383610
            #...                                              ...
            #2024      CEN    Bottom   txy             -28.961452
            #2025      CEN    Top      angle           -21.011902
            #                 Bottom   omax            -23.810177
            #2033      CEN    Top      omin           -110.334686
            #                 Bottom   von_mises       100.566292
            #
            column_names, column_values = self._build_dataframe_transient_header()
            names = ['ElementID', 'NodeID', 'Location', 'Item']
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, element_node, self.data, from_tuples=False, from_array=True,
                names=names,
            )
        else:
            # option B - nice!
            df1 = pd.DataFrame(element_node).T
            df1.columns = ['ElementID', 'NodeID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            data_frame = df1.join(df2)
            data_frame = data_frame.reset_index().set_index(['ElementID', 'NodeID', 'Location'])
        self.data_frame = data_frame

    @classmethod
    def _add_case(cls,
                  table_name, element_name, isubcase,
                  is_sort1, is_random, is_msc,
                  random_code, title, subtitle, label):
        is_strain = 'Strain' in cls.__name__

        num_wide = ELEMENT_NAME_TO_NUM_WIDE[element_name]
        data_code = oes_real_data_code(table_name,
                                       element_name, num_wide,
                                       is_sort1=is_sort1, is_random=is_random,
                                       random_code=random_code,
                                       title=title, subtitle=subtitle, label=label,
                                       is_msc=is_msc)

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #stress / strain -> 1, 3
        # stress_bits[2] = 0 # curvature
        # stress_bits[4] = 1 # von mises (vs. max shear)
        fiber = 1
        von_mises = 0
        if is_strain:
            strain = 1
        else:
            strain = 0
        stress_bits = [0, strain, fiber, strain, von_mises]
        s_code = get_scode(stress_bits)

        # stress
        assert stress_bits[1] == stress_bits[3], 'stress_bits=%s' % (stress_bits)
        data_code['stress_bits'] = stress_bits
        data_code['s_code'] = s_code

        element_type = ELEMENT_NAME_TO_ELEMENT_TYPE[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        return data_code

    @classmethod
    def add_static_case(cls, table_name, element_name, nnodes, element_node, fiber, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_static_case(cls, is_sort1, isubcase, data_code,
                              set_element_node_fiber_case, (nnodes, element_node, fiber, data))
        return obj

    @classmethod
    def add_modal_case(cls, table_name, element_name: str, nnodes, element_node, fiber, data, isubcase,
                       modes, eigns, cycles,
                       is_sort1=True, is_random=False, is_msc=True,
                       random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_modal_case(cls, is_sort1, isubcase, data_code,
                             set_element_node_fiber_case, (nnodes, element_node, fiber, data),
                             modes, eigns, cycles)
        return obj

    @classmethod
    def add_transient_case(cls, table_name, element_name, nnodes, element_node, fiber, data, isubcase,
                           times,
                           is_sort1=True, is_random=False, is_msc=True,
                           random_code=0, title='', subtitle='', label=''):
        data_code = cls._add_case(
            table_name, element_name,
            isubcase, is_sort1, is_random, is_msc,
            random_code, title, subtitle, label)
        obj = set_transient_case(cls, is_sort1, isubcase, data_code,
                                 set_element_node_fiber_case, (nnodes, element_node, fiber, data),
                                 times)
        return obj

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, element_nodei in enumerate(self.element_node):
                    (eid, nid) = element_nodei
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_eid_sort1(self, dt, eid, node_id,
                          fiber_dist1, oxx1, oyy1, txy1, angle1, major_principal1, minor_principal1, ovm1,
                          fiber_dist2, oxx2, oyy2, txy2, angle2, major_principal2, minor_principal2, ovm2):
        assert isinstance(eid, integer_types), eid
        assert isinstance(node_id, integer_types), node_id
        self._times[self.itime] = dt
        #assert self.itotal == 0, oxx
        self.element_node[self.itotal, :] = [eid, node_id]
        self.element_node[self.itotal+1, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist1, oxx1, oyy1, txy1, angle1,
                                                 major_principal1, minor_principal1, ovm1]
        self.data[self.itime, self.itotal+1, :] = [fiber_dist2, oxx2, oyy2, txy2, angle2,
                                                   major_principal2, minor_principal2, ovm2]
        self.itotal += 2
        self.ielement += 2

    def add_sort1(self, dt, eid, node_id,
                  fiber_dist1, oxx1, oyy1, txy1, angle1, major_principal1, minor_principal1, ovm1,
                  fiber_dist2, oxx2, oyy2, txy2, angle2, major_principal2, minor_principal2, ovm2):
        assert eid is not None, eid
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, integer_types), node_id
        self.element_node[self.itotal, :] = [eid, node_id]
        self.element_node[self.itotal+1, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist1, oxx1, oyy1, txy1, angle1,
                                                 major_principal1, minor_principal1, ovm1]
        self.data[self.itime, self.itotal+1, :] = [fiber_dist2, oxx2, oyy2, txy2, angle2,
                                                   major_principal2, minor_principal2, ovm2]
        self.itotal += 2
        #self.ielement += 2

    def _get_sort2_itime_ilower_iupper_from_itotal(self, dt, eid: int, nid: int,
                                                   debug=False) -> tuple[int, int, int]:
        ntimes = self.data.shape[0]

        # the monotonic element index (no duplicates)
        ielement = self.itime

        # itime = self.itime
        itime = self.ielement
        #ie_upper = self.ielement
        #ie_lower = self.ielement + 1
        #itotal = self.itotal
        #inid = 0
        nnodes = self.nnodes_per_element
        #itime = self.ielement // nnodes

        #ilayer = self.itotal % 2 == 0 # 0/1
        #ielement_inid = self.itotal // ntimes
        inid = self.itotal // (2 * ntimes)
        if self.element_type in {33, 74, 227, 228}:
            #  CQUAD4-33, CTRIA3-74, CTRIAR-227, CQUADR-228
            assert inid == 0, (self.element_name, self.element_type, inid)
            #print('inid', inid)
        elif self.element_type in [64, 144]:  # CQUAD8, CQUAD4-144
            assert inid in (0, 1, 2, 3, 4), (self.element_name, self.element_type, inid)
        elif self.element_type == 75:  # CQUAD8
            assert inid in (0, 1, 2, 3), (self.element_name, self.element_type, inid)
        else:
            raise NotImplementedError((self.element_name, self.element_type, inid))


        #inid = self.ielement % nnodes
        #itotal = self.itotal
        #if itime >= self.data.shape[0]:# or itotal >= self.element_node.shape[0]:
        ielement = self.itime
        #if self.element_name == 'CQUAD8':
            #print(f'*SORT2 {self.element_name}: itime={itime} ielement={ielement} ilayer={ilayer}  inid={inid} itotal={itotal} dt={dt} eid={eid} nid={nid}')
            #print(f'*SORT2 {self.element_name}: itime={itime} ielement={ielement} ilayer=False inid={inid} itotal={itotal+1} dt={dt} eid={eid} nid={nid}')
            #print(self.data.shape)
            #print(self.element_node.shape)
        #else:
        #aaa
        #print(itime, inid, ielement)

        #ibase = 2 * ielement # ctria3/cquad4-33
        if debug:
            print(f'ielement={ielement} nnodes={nnodes} inid={inid}')
        ibase = 2 * (ielement * nnodes + inid)
        #ibase = ielement_inid
        ie_upper = ibase
        ie_lower = ibase + 1

        #if self.element_name == 'CTRIAR': # and self.table_name == 'OESATO2':
        #debug = False
        #if self.element_name == 'CTRIAR': # and self.table_name in ['OSTRRMS1', 'OSTRRMS2']:
            #debug = True
        #if debug:
            #print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_upper={ie_upper} ielement={self.itime} inid={inid} nid={nid} itotal={itotal} dt={dt} eid={eid} nid={nid}')
            #print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_lower={ie_lower} ielement={self.itime} inid={inid} nid={nid} itotal={itotal+1} dt={dt} eid={eid} nid={nid}')
        return itime, ie_upper, ie_lower

    def add_new_eid_sort2(self, dt, eid, node_id,
                          fiber_dist1, oxx1, oyy1, txy1, angle1, major_principal1, minor_principal1, ovm1,
                          fiber_dist2, oxx2, oyy2, txy2, angle2, major_principal2, minor_principal2, ovm2):
        assert isinstance(eid, integer_types), eid
        assert isinstance(node_id, integer_types), node_id
        #itime, itotal = self._get_sort2_itime_ielement_from_itotal()


        itime, ie_upper, ie_lower = self._get_sort2_itime_ilower_iupper_from_itotal(dt, eid, node_id)
        try:
            #print(f'SORT2: itime={itime} -> dt={dt};   ie_upper={ie_upper} -> eid={eid} ({self.element_name})')
            self._times[itime] = dt
            #assert self.itotal == 0, oxx
            #if itime == 0:
            self.element_node[ie_upper, :] = [eid, node_id]  # 0 is center
            self.element_node[ie_lower, :] = [eid, node_id]  # 0 is center
        except Exception:
            itime, ie_upper, ie_lower = self._get_sort2_itime_ilower_iupper_from_itotal(
                dt, eid, node_id, debug=True)
            print(f'SORT2: itime={itime} -> dt={dt};   ie_upper={ie_upper} -> eid={eid} ({self.element_name})')
            raise
        #print(self.element_node)
        #self.data[self.itime, ie_upper, :] = [fiber_dist1, oxx1, oyy1, txy1, angle1,
                                              #major_principal1, minor_principal1, ovm1]
        #self.data[self.itime, ie_lower, :] = [fiber_dist2, oxx2, oyy2, txy2, angle2,
                                              #major_principal2, minor_principal2, ovm2]
        self.itotal += 2
        #self.ielement += 1

    def add_sort2(self, dt, eid, node_id,
                  fiber_dist1, oxx1, oyy1, txy1, angle1, major_principal1, minor_principal1, ovm1,
                  fiber_dist2, oxx2, oyy2, txy2, angle2, major_principal2, minor_principal2, ovm2):
        assert eid is not None, eid
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, integer_types), node_id
        itime, ie_upper, ie_lower = self._get_sort2_itime_ilower_iupper_from_itotal(dt, eid, node_id)
        #print(f'SORT2b: itime={itime} -> dt={dt};   ie_upper={ie_upper} -> eid={eid} nid={node_id}')
        #print(self.element_node.shape)
        #if itime == 0:
        self.element_node[ie_upper, :] = [eid, node_id]
        self.element_node[ie_lower, :] = [eid, node_id]
        #print(self.element_node.tolist())
        self.data[itime, ie_upper, :] = [fiber_dist1, oxx1, oyy1, txy1, angle1,
                                       major_principal1, minor_principal1, ovm1]
        self.data[itime, ie_lower, :] = [fiber_dist2, oxx2, oyy2, txy2, angle2,
                                         major_principal2, minor_principal2, ovm2]
        self.itotal += 2
        #self.ielement += 2

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  data.shape={self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg.append(f'  s_code: {self.s_code}\n')
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = np.searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg, nnodes, cen = _get_plate_msg(self)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #cen_word = 'CEN/%i' % nnodes
        cen_word = cen
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            major_principal = self.data[itime, :, 5]
            minor_principal = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            is_linear = self.element_type in {33, 74, 227, 228, 83}
            is_bilinear = self.element_type in {64, 70, 75, 82, 144}
            for (i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 count(), eids, nids, fiber_dist, oxx, oyy, txy, angle, major_principal, minor_principal, ovm):
                [fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                    [fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                ilayer = i % 2
                # tria3
                if is_linear:  # CQUAD4, CTRIA3, CTRIAR linear, CQUADR linear
                    if ilayer == 0:
                        f06_file.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:
                        f06_file.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                elif is_bilinear:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    # bilinear
                    if nid == 0 and ilayer == 0:  # CEN
                        f06_file.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 0:
                        f06_file.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            '', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 1:
                        f06_file.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % (
                            '', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                else:  # pragma: no cover
                    msg = 'element_name=%s self.element_type=%s' % (
                        self.element_name, self.element_type)
                    raise NotImplementedError(msg)

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def get_nnodes_bilinear(self):
        """gets the number of nodes and whether or not the element has bilinear results"""
        is_bilinear = False
        if self.element_type == 74:
            nnodes = 3
        elif self.element_type == 33:
            nnodes = 4
        elif self.element_type == 144:
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 82:  # CQUADR
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 64:  # CQUAD8
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 75:  # CTRIA6
            nnodes = 3
            is_bilinear = True
        elif self.element_type == 70:  # CTRIAR
            nnodes = 3
            is_bilinear = True
        elif self.element_type == 227:  # CTRIAR-linear
            nnodes = 3
            is_bilinear = False
        elif self.element_type == 228:  # CQUADR-linear
            nnodes = 4
            is_bilinear = False
        else:
            raise NotImplementedError(f'name={self.element_name} type={self.element_type}')
        return nnodes, is_bilinear

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        nnodes, is_bilinear = self.get_nnodes_bilinear()
        if is_bilinear:
            nnodes_all = nnodes + 1
            ntotal = 2 + 17 * nnodes_all
        else:
            nnodes_all = nnodes
        #print("nnodes_all =", nnodes_all)
        #cen_word_ascii = f'CEN/{nnodes:d}'
        cen_word_bytes = b'CEN/'
        idtype = self.element_node.dtype
        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            pass
        else:
            print(f'downcasting {self.class_name}...')
            #cen_word_bytes = b'CEN/    '
            idtype = np.int32(1)
            fdtype = np.float32(1.0)

        #msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        #msg.append(f'  data.shape={self.data.shape}\n')

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        max_id = self.element_node.max()
        if max_id > 99999999:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max id={max_id}')

        eids_device = eids * 10 + self.device_code

        nelements = len(np.unique(eids))
        nlayers = len(eids)

        #print('nelements =', nelements)
        #print('nlayers =', nlayers)
        nnodes_per_element = nlayers // nelements // 2
        assert nnodes_per_element in [1, 4, 5], nnodes_per_element
        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements
        assert nnodes > 1, nnodes

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        op2_ascii.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')  # 1+16

        op2_ascii.write('  #elementi = [eid_device, node1, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n') # 1 + 17*5
        op2_ascii.write('  #elementi = [            node2, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n') # 17
        op2_ascii.write('  #elementi = [            node3, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n') # 17
        op2_ascii.write('  #elementi = [            node4, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n') # 17
        op2_ascii.write('  #elementi = [            node5, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                               fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n') # 17

        if not self.is_sort1:
            raise NotImplementedError('SORT2')

        #struct_isi8f = Struct('i 4s i 8f')
        #struct_i8f = Struct(endian + b'i8f')
        #struct_8f = Struct(endian + b'8f')

        nelements_nnodes = len(nids) // 2
        is_centroid = self.element_type in [33, 74, 227, 228]
        is_nodes = self.element_type in [64, 70, 75, 82, 144]
        if is_centroid:
            eids_device2 = to_column_bytes([eids_device[::2]], idtype).view(fdtype)
            assert len(eids_device2) == nelements
        elif is_nodes:
            cen_word_array_temp = np.full((nelements, 1), cen_word_bytes)
            cen_word_array = cen_word_array_temp.view(fdtype)
            eids_device2 = view_idtype_as_fdtype(eids_device[::2*nnodes_per_element].reshape(nelements, 1),
                                                 fdtype)
            nids2 = view_idtype_as_fdtype(nids[::2].reshape(nelements_nnodes, 1),
                                          fdtype)

        #nheader = 15
        struct_i = Struct('i')
        struct_13i = Struct('13i')
        op2_ascii.write(f'nelements={nelements:d}\n')
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(struct_13i.pack(*header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            datai_ = self.data[itime, :, :]
            # 16 values
            # [eid_device, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
            # [            fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
            if is_centroid:
                datai2_ = datai_.reshape(nelements, 16)
                datai = view_dtype(datai2_, fdtype)
                data_out = np.hstack([eids_device2, datai])
            elif is_nodes:
                # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                # bilinear
                #16
                datai2_ = datai_.reshape(nelements*nnodes_per_element, 16)
                datai = view_dtype(datai2_, fdtype)
                nids_data = np.hstack([nids2, datai]).reshape(nelements, nnodes_per_element*17)
                data_out = np.hstack([eids_device2, cen_word_array, nids_data])
            else:  # pragma: no cover
                msg = f'element_name={self.element_name} element_type={self.element_type}'
                raise NotImplementedError(msg)
            assert data_out.size == ntotal, f'data_out.shape={data_out.shape} size={data_out.size}; ntotal={ntotal}'
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]

            op2_file.write(struct_i.pack(*header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealPlateStressArray(RealPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        fiber_dist = 'fiber_distance' if self.is_fiber_distance else 'fiber_curvature'
        ovm = 'von_mises' if self.is_von_mises else 'max_shear'
        headers = [fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        return headers


class RealPlateStrainArray(RealPlateArray, StrainObject):
    """
    used for:
     - RealPlateStressArray
     - RealPlateStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        fiber_dist = 'fiber_distance' if self.is_fiber_distance else 'fiber_curvature'
        ovm = 'von_mises' if self.is_von_mises else 'max_shear'
        headers = [fiber_dist, 'exx', 'eyy', 'exy', 'angle', 'emax', 'emin', ovm]
        return headers


def _get_plate_msg(self):
    von_mises = 'VON MISES' if self.is_von_mises else 'MAX SHEAR'

    if self.is_stress:
        if self.is_fiber_distance:
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        cquad4_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                        S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp
        #cquadr_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        ctria3_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp
    else:
        if self.is_fiber_distance:
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        cquad4_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp

        #cquadr_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp
        ctria3_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp

    if self.element_type in [74, 83]:
        msg = ctria3_msg
        nnodes = 3
        cen = 'CEN/3'
    elif self.element_type == 33:
        msg = cquad4_msg
        nnodes = 4
        cen = 'CEN/4'
    #elif self.element_type == 228:
        #msg = cquadr_msg
        #nnodes = 4
        #cen = None # 'CEN/4'

    elif self.element_type == 144:
        msg = cquad4_bilinear_msg
        nnodes = 4
        cen = 'CEN/4'
    elif self.element_type in [82, 228]:  # CQUADR bilinear, CQUADR linear
        msg = cquadr_msg
        nnodes = 4
        cen = 'CEN/4'
    elif self.element_type == 64:  # CQUAD8
        msg = cquad8_msg
        nnodes = 4
        cen = 'CEN/8'
    elif self.element_type == 75:  # CTRIA6
        msg = ctria6_msg
        nnodes = 3
        cen = 'CEN/6'
    elif self.element_type in [70, 227]:
        # 70: CTRIAR bilinear
        # 227: CTRIAR linear
        msg = ctriar_msg
        nnodes = 3
        cen = 'CEN/3'
    else:  # pragma: no cover
        raise NotImplementedError(f'name={self.element_name} type={self.element_type}')
    return msg, nnodes, cen

def set_element_node_fiber_case(cls, data_code, is_sort1, isubcase,
                                nnodes, element_node, fiber_distance, data, times):
    assert element_node.ndim == 2, element_node.shape
    assert element_node.shape[1] == 2, element_node.shape
    assert fiber_distance.ndim == 1, fiber_distance.shape
    assert element_node[:, 0].min() > 0, element_node
    assert element_node[:, 1].min() == 0, element_node
    ntimes = data.shape[0]
    nlayers = data.shape[1]
    dt = times[0]
    obj = cls(data_code, is_sort1, isubcase, dt)
    obj.element_node = element_node
    obj.fiber_distance = fiber_distance
    obj.data = data

    #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm
    #assert data.shape[2] == 8, data.shape

    # nlayers: includes all plies
    # assumed 2 layers per node
    # a CTRIA3 has 1 node; a CQUAD4-144 has 4+1=5 nodes (centroid is included)
    # nelements should be the *actual* number of elements, not the size of element_node
    #
    # TODO: element_node should not contain duplicates?
    nelements = nlayers // (nnodes * 2)
    assert nelements >= 1, nelements

    obj.ntimes = ntimes
    obj.ntotal = nnodes
    obj.nelements = nelements
    obj._times = times
    obj.nnodes = nnodes
    #obj.update_data_components()
    return obj
