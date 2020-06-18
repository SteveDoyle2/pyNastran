#pylint disable=C0301
from struct import Struct, pack
from abc import abstractmethod
import inspect
from typing import List

import numpy as np
from numpy import zeros, searchsorted, allclose

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.op2.result_objects.op2_objects import BaseElement, get_times_dtype
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_12e,
    write_float_13e, # write_float_12e,
    _eigenvalue_header,
)
from pyNastran.op2.op2_interface.write_utils import set_table3_field


SORT2_TABLE_NAME_MAP = {
    'OEF2' : 'OEF1',
    'OEFATO2' : 'OEFATO1',
    'OEFCRM2' : 'OEFCRM1',
    'OEFPSD2' : 'OEFPSD1',
    'OEFRMS2' : 'OEFRMS1',
    'OEFNO2' : 'OEFNO1',
}
TABLE_NAME_TO_TABLE_CODE = {
    'OEF1' : 4,
}

class ForceObject(BaseElement):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = np.nan
        self.element = None
        self._times = None
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    def finalize(self):
        """it's required that the object be in SORT1"""
        self.set_as_sort1()

    def set_as_sort1(self):
        """the data is in SORT1, but the flags are wrong"""
        if self.is_sort1:
            return
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_bits[1] = 0 # sort1
        self.sort_method = 1
        assert self.is_sort1 is True, self.is_sort1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError()

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def _write_table_3(self, op2, op2_ascii, new_result, itable, itime): #itable=-3, itime=0):
        import inspect
        from struct import pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        #print('new_result=%s itable=%s' % (new_result, itable))
        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('table_3_header = %s\n' % header)

        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        element_type = self.element_type
        assert isinstance(self.element_type, int), self.element_type
        #[
            #'aCode', 'tCode', 'element_type', 'isubcase',
            #'???', '???', '???', 'load_set'
            #'format_code', 'num_wide', 's_code', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', 'Title', 'subtitle', 'label']
        #random_code = self.random_code
        format_code = self.format_code
        s_code = 0 # self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        #oCode = 0
        load_set = 0
        #print(self.code_information())

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0
        if self.analysis_code == 1:
            field5 = self.loadIDs[itime]
        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            field6 = self.eigns[itime]
            field7 = self.cycles[itime]
            assert isinstance(field6, float), type(field6)
            assert isinstance(field7, float), type(field7)
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 5:
            try:
                field5 = self.freqs[itime]
            except AttributeError:  # pragma: no cover
                print(self)
                raise
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            if hasattr(self, 'times'):
                field5 = self.times[itime]
            #elif hasattr(self, 'dts'):
                #field5 = self.times[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find times or dts on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 7:  # pre-buckling
            field5 = self.loadIDs[itime] # load set number
        elif self.analysis_code == 8:  # post-buckling
            if hasattr(self, 'lsdvmns'):
                field5 = self.lsdvmns[itime] # load set number
            elif hasattr(self, 'loadIDs'):
                field5 = self.loadIDs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find lsdvmns or loadIDs on analysis_code=8')

            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            assert isinstance(field6, float_types), type(field6)
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        elif self.analysis_code == 9:  # complex eigenvalues
            field5 = self.modes[itime]
            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            field7 = self.eigis[itime]
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.load_steps[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            field5 = self.loadIDs[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, element_type, isubcase, field5,
            field6, field7, load_set, format_code, num_wide,
            s_code, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        n = 0
        for v in table3:
            if isinstance(v, (int, float, np.float32)):
                n += 4
            elif isinstance(v, str):
                #print(len(v), v)
                n += len(v)
            else:
                #print('write_table_3', v)
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2.write(pack(fmt, *data))


class RealForceObject(ForceObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        ForceObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False


"""
          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )
   ELEMENT  FAILURE    PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG
     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES
         3   STRAIN      1        20345.4805   -2
                                                                       7.1402
                         2            0.9025   -12
                                                                       7.1402
                         3        20342.2461   -2
                                                                                                 20345.4805             ***
         4   STRAIN      1        16806.9277   -2
                                                                      38.8327
                         2            0.9865   -2
                                                                      38.8327
                         3        16804.4199   -2

          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )
   ELEMENT  FAILURE    PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG
     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES
         5   STRAIN      1        21850.3184   -2
                                                                  166984.4219
                         2            0.7301   -2
                                                                  166984.4219
                         3        21847.9902   -2
                                                                                                166984.4219             ***
         6   STRAIN      1        18939.8340   -2
                                                                  130371.3828
                         2            0.7599   -1
                                                                  130371.3828
                         3        18937.7734   -2
          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
   ELEMENT  FAILURE    PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG
     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES
         1   STRAIN      1        18869.6621   -2
                                                                      16.2471
                         2            1.0418   -2
                                                                      16.2471
                         3        18866.6074   -2
                                                                                                 18869.6621             ***
1   CC227: CANTILEVERED COMPOSITE PLATE 3 LAYER SYMM PLY    CC227          DECEMBER   5, 2011  MSC.NASTRAN  6/17/05   PAGE    15
      FAILURE CRITERION IS STRAIN, STRESS ALLOWABLES, LIST STRESSES
0

          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )
   ELEMENT  FAILURE    PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG
     ID      THEORY     ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES
         2   STRAIN      1        14123.7451   -2
                                                                      31.4861
                         2            1.0430   -2
                                                                      31.4861
                         3        14122.1221   -2
"""
class FailureIndicesArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific

    def build(self):
        """sizes the vectorized attributes of the FailureIndices"""
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)

        self._times = zeros(self.ntimes, dtype=dtype)
        self.failure_theory = np.full(self.nelements, '', dtype='U8')
        self.element_layer = zeros((self.nelements, 2), dtype=idtype)

        #[failure_stress_for_ply, interlaminar_stress, max_value]
        self.data = zeros((self.ntimes, self.nelements, 3), dtype=fdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        element_layer = [self.element_layer[:, 0], self.element_layer[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            # Time                                                              0.00          0.05
            # ElementID NodeID Item
            # 2         1      failure_index_for_ply (direct stress/strain)      0.0  5.431871e-14
            #           2      failure_index_for_bonding (interlaminar stresss)  0.0  3.271738e-16
            #           3      max_value                                         NaN           NaN
            #           1      failure_index_for_ply (direct stress/strain)      0.0  5.484873e-30
            #           2      failure_index_for_bonding (interlaminar stresss)  0.0  3.271738e-16
            #           3      max_value                                         NaN           NaN
            #           1      failure_index_for_ply (direct stress/strain)      0.0  5.431871e-14
            #           2      failure_index_for_bonding (interlaminar stresss)  NaN           NaN
            #           3      max_value                                         0.0  5.431871e-14
            column_names, column_values = self._build_dataframe_transient_header()
            names = ['ElementID', 'Layer', 'Item']
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names, headers,
                element_layer, self.data, names=names,
                from_tuples=False, from_array=True)

            #column_names, column_values = self._build_dataframe_transient_header()
            #data_frame = pd.Panel(self.data, items=column_values,
                                  #major_axis=element_layer, minor_axis=headers).to_frame()
            #data_frame.columns.names = column_names
            #data_frame.index.names = ['ElementID', 'Layer', 'Item']
        else:
            #Static           failure_index_for_ply (direct stress/strain)  failure_index_for_bonding (interlaminar stresss)     max_value
            #ElementID Layer
            #101       1                                      7.153059e-07                                               0.0           NaN
            #          2                                      1.276696e-07                                               0.0           NaN
            #          3                                      7.153059e-07                                               NaN  7.153059e-07
            element_layer = [self.element_layer[:, 0], self.element_layer[:, 1]]
            index = pd.MultiIndex.from_arrays(element_layer, names=['ElementID', 'Layer'])

            data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def get_headers(self) -> List[str]:
        #headers = ['eid', 'failure_theory', 'ply', 'failure_index_for_ply (direct stress/strain)',
                   #'failure_index_for_bonding (interlaminar stresss)', 'failure_index_for_element', 'flag']
        headers = ['failure_index_for_ply (direct stress/strain)',
                   'failure_index_for_bonding (interlaminar stresss)', 'max_value']
        return headers

    def __eq__(self, table):  # pragma: no cover
        return True

    def add_sort1(self, dt, eid, failure_theory, ply_id, failure_stress_for_ply, flag,
                  interlaminar_stress, max_value, nine):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element_layer[self.ielement] = [eid, ply_id]
        self.failure_theory[self.ielement] = failure_theory
        self.data[self.itime, self.ielement, :] = [failure_stress_for_ply, interlaminar_stress, max_value]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes = self.data.shape[0]
        nelements = self.data.shape[1]
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        return [] # raise NotImplementedError('this should be overwritten by %s' % (self.__class__.__name__))

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        f06_file.write('skipping FailureIndices f06\n')

        return page_num
        #NotImplementedError(self.code_information())
        #asd
        #if self.is_sort1:
            #page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp)
        #else:
            #raise NotImplementedError(self.code_information())
            #page_num = self._write_sort2_as_sort2(header, page_stamp, page_num, f06_file, msg_temp)

        #'          F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'
        #'   ELEMENT  FAILURE        PLY   FP=FAILURE INDEX FOR PLY    FB=FAILURE INDEX FOR BONDING   FAILURE INDEX FOR ELEMENT      FLAG\n'
        #'     ID      THEORY         ID  (DIRECT STRESSES/STRAINS)     (INTER-LAMINAR STRESSES)      MAX OF FP,FB FOR ALL PLIES\n'
        #'         1   HOFFMAN       101      6.987186E-02      \n'
        #'                                                                     1.687182E-02                                              \n'
        #'                           102      9.048269E-02      \n'
        #'                                                                     1.721401E-02                                               \n'
        #return page_num


class RealSpringDamperForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific
        #if not is_sort1:
            #raise NotImplementedError('SORT2')

    @classmethod
    def add_static_case(cls, table_name, element_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):

        analysis_code = 1 # static
        data_code = oef_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['loadIDs'] = [0] # TODO: ???
        data_code['data_names'] = []

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #if is_stress:
            #data_code['stress_bits'] = [0, 0, 0, 0]
            #data_code['s_code'] = 0
        #else:
            #data_code['stress_bits'] = [0, 1, 0, 1]
            #data_code['s_code'] = 1 # strain?
        element_name_to_element_type = {
            'CELAS1' : 11,
            'CELAS2' : 12,
            'CELAS3' : 13,
            'CELAS4' : 14,
        }

        element_type = element_name_to_element_type[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        #data_code['load_set'] = 1

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.element = element
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    def build(self):
        """sizes the vectorized attributes of the RealSpringDamperForceArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self.build_data(self.ntimes, self.nelements, dtype, idtype, fdtype)

    def build_data(self, ntimes, nelements, dtype, idtype, fdtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements

        self._times = zeros(ntimes, dtype=dtype)
        self.element = zeros(nelements, dtype=idtype)

        #[force]
        self.data = zeros((ntimes, nelements, 1), dtype=fdtype)


    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                               1             2             3
            #Freq                    1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue             -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians                 9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #30        spring_force  2.388744e-19 -1.268392e-10 -3.341473e-19
            #31        spring_force  2.781767e-19 -3.034770e-11 -4.433221e-19
            #32        spring_force  0.000000e+00  0.000000e+00  0.000000e+00
            #33        spring_force  0.000000e+00  0.000000e+00  0.000000e+00
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static      spring_force
            #ElementID
            #30                   0.0
            #31                   0.0
            #32                   0.0
            #33                   0.0
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame


    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor not in (None, np.nan):
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid1, Eid2\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
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
                        (force1, stress1) = t1
                        (force2, stress2) = t2
                        if not allclose(t1, t2):
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

    def add_sort1(self, dt, eid, force):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #print('dt=%s eid=%s' % (dt, eid))
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes = self.data.shape[0]
        nelements = self.data.shape[1]
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        raise NotImplementedError('this should be overwritten by %s' % (self.__class__.__name__))

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp)
        else:
            raise NotImplementedError(self.code_information())
            #page_num = self._write_sort2_as_sort2(header, page_stamp, page_num, f06_file, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        nwrite = len(eids)
        nrows = nwrite // 4
        nleftover = nwrite - nrows * 4

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))
            stress = self.data[itime, :, 0]

            out = []
            for eid, stressi in zip(eids, stress):
                out.append([eid, write_float_13e(stressi)])

            for i in range(0, nrows * 4, 4):
                f06_file.write('    %10i  %13s    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2] + out[i + 3])))

            i = nrows * 4
            if nleftover == 3:
                f06_file.write('    %10i  %13s    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1] + out[i + 2])))
            elif nleftover == 2:
                f06_file.write('    %10i  %13s    %10i  %13s\n' % (
                    tuple(out[i] + out[i + 1])))
            elif nleftover == 1:
                f06_file.write('    %10i  %13s\n' % tuple(out[i]))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

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

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'if')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            force = self.data[itime, :, 0]

            for eid, forcei in zip(eids_device, force):
                data = [eid, forcei]
                op2_ascii.write('  eid=%s force=%s\n' % tuple(data))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealSpringForceArray(RealSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def get_headers(self) -> List[str]:
        headers = ['spring_force']
        return headers

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 11:  # CELAS1
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )\n']
        elif self.element_type == 12:  # CELAS2
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n']
        elif self.element_type == 13:  # CELAS3
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )\n']
        elif self.element_type == 14:  # CELAS4
            msg = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )\n']
        else:  # pragma: no cover
            msg = 'element_name=%s element_type=%s' % (self.element_name, self.element_type)
            raise NotImplementedError(msg)

        msg += [
            '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n'
            '        ID.                              ID.                              ID.                              ID.\n'
        ]
        return msg


class RealDamperForceArray(RealSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def get_headers(self) -> List[str]:
        headers = ['damper_force']
        return headers

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 20:  # CDAMP1
            msg = ['                              F O R C E S   I N   S C A L A R   D A M P E R S        ( C D A M P 1 )\n']
        elif self.element_type == 21:  # CDAMP2
            msg = ['                              F O R C E S   I N   S C A L A R   D A M P E R S        ( C D A M P 2 )\n']
        elif self.element_type == 22:  # CDAMP3
            msg = ['                              F O R C E S   I N   S C A L A R   D A M P E R S        ( C D A M P 3 )\n']
        elif self.element_type == 23:  # CDAMP4
            msg = ['                              F O R C E S   I N   S C A L A R   D A M P E R S        ( C D A M P 4 )\n']
        else:  # pragma: no cover
            msg = 'element_name=%s element_type=%s' % (self.element_name, self.element_type)
            raise NotImplementedError(msg)

        if is_sort1:
            msg += [
                '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n'
                '        ID.                              ID.                              ID.                              ID.\n'
            ]
        else:
            msg += [
                '                         AXIAL                                                       AXIAL\n'
                '       TIME              FORCE         TORQUE                      TIME              FORCE         TORQUE\n'
            ]
        return msg


class RealRodForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific

    @classmethod
    def add_static_case(cls, table_name, element_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):

        analysis_code = 1 # static
        data_code = oef_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['loadIDs'] = [0] # TODO: ???
        data_code['data_names'] = []

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #if is_stress:
            #data_code['stress_bits'] = [0, 0, 0, 0]
            #data_code['s_code'] = 0
        #else:
            #data_code['stress_bits'] = [0, 1, 0, 1]
            #data_code['s_code'] = 1 # strain?
        element_name_to_element_type = {
            'CROD' : 1,
            'CTUBE' : 3,
            'CONROD' : 10,
        }

        element_type = element_name_to_element_type[element_name]
        data_code['element_name'] = element_name
        data_code['element_type'] = element_type
        #data_code['load_set'] = 1

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.element = element
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def get_headers(self) -> List[str]:
        headers = ['axial', 'torsion']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                    '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        crod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

    def build(self):
        """sizes the vectorized attributes of the RealRodForceArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self.build_data(self.ntimes, self.nelements, float_fmt='float32')

    def build_data(self, ntimes, nelements, float_fmt='float32'):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        #self.ntotal = ntimes * nelements
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(ntimes, dtype=dtype)
        self.element = zeros(nelements, dtype='int32')

        #[axial_force, torque]
        self.data = zeros((ntimes, nelements, 2), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def add_sort1(self, dt, eid, axial, torque):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1
        if self.ielement == self.nelements:
            self.ielement = 0

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        if 'CROD' in self.element_name:
            msg = crod_msg
        elif 'CONROD' in self.element_name:
            msg = conrod_msg
        elif 'CTUBE' in self.element_name:
            msg = ctube_msg
        return self.element_name, msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [axiali, torsioni] = write_floats_13e([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %s\n' % tuple(out[i] + out[i + 1])
                f06_file.write(out_line)
            if is_odd:
                out_line = '      %8i   %-13s  %s\n' % tuple(out[-1])
                f06_file.write(out_line)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor not in (None, np.nan):
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
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

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i2f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            #print('eids3', eids3)
            for eid, axiali, torsioni in zip(eids_device, axial, torsion):
                data = [eid, axiali, torsioni]
                op2_ascii.write('  eid=%s axial=%s torsion=%s\n' % tuple(data))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealCBeamForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        #ForceObject.__init__(self, data_code, isubcase)
        RealForceObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        self.itime = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            ##sort1
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def build(self):
        """sizes the vectorized attributes of the RealCBeamForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 11

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
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype)
        self.element = zeros(self.ntotal, idtype)
        self.element_node = zeros((self.ntotal, 2), idtype)

        # the number is messed up because of the offset for the element's properties
        if not (self.nelements * nnodes) == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.data = zeros((self.ntimes, self.ntotal, 8), fdtype)

    def finalize(self):
        sd = self.data[0, :, 0]
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)

        #self.nelements = len(self.element) // 11
        self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        element_location = [
            self.element_node[:, 0],
            self.data[0, :, 0],
        ]
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                                           1             2             3
            #Freq                                1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue                         -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians                             9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Location Item
            #12        0.0      bending_moment1  1.505494e-13 -2.554764e-07 -5.272747e-13
            #                   bending_moment2 -2.215085e-13 -2.532377e-07  3.462328e-13
            #                   shear1           1.505494e-13 -2.554763e-07 -5.272747e-13
            #                   shear2          -2.215085e-13 -2.532379e-07  3.462328e-13
            #                   axial_force      1.294136e-15 -1.670896e-09  4.759476e-16
            #                   total_torque    -4.240346e-16  2.742446e-09  1.522254e-15
            #                   warping_torque   0.000000e+00  0.000000e+00  0.000000e+00
            #          1.0      bending_moment1  0.000000e+00 -1.076669e-13  1.009742e-28
            #                   bending_moment2 -5.048710e-29  1.704975e-13  0.000000e+00
            #                   shear1           1.505494e-13 -2.554763e-07 -5.272747e-13
            #                   shear2          -2.215085e-13 -2.532379e-07  3.462328e-13
            #                   axial_force      1.294136e-15 -1.670896e-09  4.759476e-16
            #                   total_torque    -4.240346e-16  2.742446e-09  1.522254e-15
            #                   warping_torque   0.000000e+00  0.000000e+00  0.000000e+00
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers[1:], element_location, self.data[:, :, 1:], from_tuples=False, from_array=True)
            data_frame.index.names = ['ElementID', 'Location', 'Item']
        else:
            df1 = pd.DataFrame(element_location).T
            df1.columns = ['ElementID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            data_frame = df1.join([df2])
        #self.data_frame = data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID'])
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1

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

    def add_new_element_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        return self.add_sort1(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.element[self.itotal] = eid
        self.element_node[self.itotal, :] = [eid, nid]
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
        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, self.table_name))
        #msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1, self.is_sort2))
        msg.append('  CBEAM\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #name = self.data_code['name']
        #if name == 'freq':
            #name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)
        #if is_sort1:
        msg_temp = [
            '                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
            '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
            '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']
        #else:
            #raise NotImplementedError('CBEAM-SORT2')

        if self.is_sort1:
            #assert self.is_sort1 is True, str(self)
            #if is_sort1:
            page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, msg_temp)
            #else:
                #self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, msg_temp)
        else:
            print(f'skipping {self.__class__.__name__} because its sort2')
            #assert self.is_sort1 is True, str(self)
        return page_num - 1

    def get_headers(self) -> List[str]:
        headers = [
            'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
            'axial_force', 'total_torque', 'warping_torque', ]
        return headers

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, msg_temp):
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        long_form = False
        if nids.min() == 0:
            msg = header + [
                '                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']

            long_form = True

        #times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            if self.nonlinear_factor not in (None, np.nan):
                dt = self._times[itime]
                dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                header[1] = dt_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            #sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq
            assert self.is_sort1 is True, str(self)
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            for eid, nid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, nids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                vals = (bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi)
                vals2 = write_floats_13e(vals)
                (sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq) = vals2

                if long_form:
                    f06_file.write('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (
                        eid, sdi, sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq))
                else:
                    if sdi == 0.:
                        f06_file.write('0  %8i\n' % eid)
                    f06_file.write('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (
                        nid, sdi, sbm1i, sbm2i, sts1i, sts2i, safi, sttrqi, swtrq))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        #long_form = False
        #if nids.min() == 0:
            #long_form = True
        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        #xxbs = self.xxb
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

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 8f')
            struct2 = Struct(endian + b'i 8f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)
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

            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            icount = 0
            nwide = 0
            ielement = 0
            assert len(eids) == len(sd)
            for eid, nid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, nids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                if icount == 0:
                    eid_device = eids_device[ielement]
                    nid = nids[ielement]
                    data = [eid_device, nid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi] # 10
                    op2.write(struct1.pack(*data))
                    ielement += 1
                    icount = 1
                elif nid > 0 and icount > 0:
                    # 11 total nodes, with 1, 11 getting an nid; the other 9 being
                    # xxb sections
                    data = [0, 0., 0., 0., 0., 0., 0., 0., 0.]
                    #print('***adding %s\n' % (10-icount))
                    for unused_i in range(10 - icount):
                        op2.write(struct2.pack(*data))
                        nwide += len(data)

                    eid_device2 = eids_device[ielement]
                    assert eid_device == eid_device2
                    nid = nids[ielement]
                    data = [nid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi] # 9
                    op2.write(struct2.pack(*data))
                    ielement += 1
                    icount = 0
                else:
                    raise RuntimeError('CBEAM op2 writer')
                    #data = [0, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc]  # 10
                    #op2.write(struct2.pack(*data))
                    #icount += 1

                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                nwide += len(data)

            assert ntotal == nwide, 'ntotal=%s nwide=%s' % (ntotal, nwide)

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealCShearForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if not is_sort1:
            #raise NotImplementedError('SORT2')

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'force41', 'force21', 'force12', 'force32', 'force23', 'force43',
            'force34', 'force14',
            'kick_force1', 'shear12', 'kick_force2', 'shear23',
            'kick_force3', 'shear34', 'kick_force4', 'shear41',
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the RealCShearForceArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[force41, force21, force12, force32, force23, force43,
        # force34, force14,
        # kick_force1, shear12, kick_force2, shear23,
        # kick_force3, shear34, kick_force4, shear41]
        self.data = zeros((self.ntimes, self.ntotal, 16), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                              1             2             3
            #Freq                   1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue            -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians                9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #22        force41     -4.025374e-14  2.935730e-08  1.017620e-13
            #          force21     -4.025374e-14  2.935730e-08  1.017620e-13
            #          force12      4.025374e-14 -2.935730e-08 -1.017620e-13
            #          force32      4.025374e-14 -2.935730e-08 -1.017620e-13
            #          force23     -4.025374e-14  2.935730e-08  1.017620e-13
            #          force43     -4.025374e-14  2.935730e-08  1.017620e-13
            #          force34      4.025374e-14 -2.935730e-08 -1.017620e-13
            #          force14      4.025374e-14 -2.935730e-08 -1.017620e-13
            #          kick_force1 -0.000000e+00  0.000000e+00  0.000000e+00
            #          shear12     -8.050749e-14  5.871460e-08  2.035239e-13
            #          kick_force2 -0.000000e+00  0.000000e+00  0.000000e+00
            #          shear23     -8.050749e-14  5.871460e-08  2.035239e-13
            #          kick_force3 -0.000000e+00  0.000000e+00  0.000000e+00
            #          shear34     -8.050749e-14  5.871460e-08  2.035239e-13
            #          kick_force4 -0.000000e+00  0.000000e+00  0.000000e+00
            #          shear41     -8.050749e-14  5.871460e-08  2.035239e-13
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor not in (None, np.nan):
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid1, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (force41a, force14a, force21a, force12a, force32a, force23a, force43a, force34a, kick_force1a, kick_force2a, kick_force3a, kick_force4a, shear12a, shear23a, shear34a, shear41a) = t1
                    (force41b, force14b, force21b, force12b, force32b, force23b, force43b, force34b, kick_force1b, kick_force2b, kick_force3b, kick_force4b, shear12b, shear23b, shear34b, shear41b) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                force41a, force14a, force21a, force12a, force32a, force23a, force43a, force34a, kick_force1a, kick_force2a, kick_force3a, kick_force4a, shear12a, shear23a, shear34a, shear41a,
                                force41b, force14b, force21b, force12b, force32b, force23b, force43b, force34b, kick_force1b, kick_force2b, kick_force3b, kick_force4b, shear12b, shear23b, shear34b, shear41b
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid,
                  force41, force14, force21, force12, force32, force23, force43, force34,
                  kick_force1, kick_force2, kick_force3, kick_force4,
                  shear12, shear23, shear34, shear41):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            force41, force14, force21, force12, force32, force23, force43, force34,
            kick_force1, kick_force2, kick_force3, kick_force4,
            shear12, shear23, shear34, shear41]
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
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                           F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n'
            ' \n'
            '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n'
            '   ELEMENT        F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n'
            '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            f14 = self.data[itime, :, 0]
            f12 = self.data[itime, :, 1]
            f21 = self.data[itime, :, 2]
            f23 = self.data[itime, :, 3]
            f32 = self.data[itime, :, 4]
            f34 = self.data[itime, :, 5]
            f43 = self.data[itime, :, 6]
            f41 = self.data[itime, :, 7]

            kick1 = self.data[itime, :, 8]
            tau12 = self.data[itime, :, 9]
            kick2 = self.data[itime, :, 10]
            tau23 = self.data[itime, :, 11]
            kick3 = self.data[itime, :, 12]
            tau34 = self.data[itime, :, 13]
            kick4 = self.data[itime, :, 14]
            tau41 = self.data[itime, :, 15]

            #zip_in = [
                #f14, f12, f21, f23, f32, f34, f43, f41,
                #kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41,
            #]
            for (eid, f14i, f12i, f21i, f23i, f32i, f34i, f43i, f41i,
                 kick1i, tau12i, kick2i, tau23i, kick3i, tau34i, kick4i, tau41i) in zip(
                    eids, f14, f12, f21, f23, f32, f34, f43, f41,
                    kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41):
                vals2 = write_floats_12e([
                    f14i, f12i, f21i, f23i, f32i, f34i, f43i, f41i,
                    kick1i, tau12i, kick2i, tau23i, kick3i, tau34i, kick4i, tau41i])
                [
                    f14i, f12i,
                    f21i, f23i,
                    f32i, f34i,
                    f43i, f41i,
                    kick1i, tau12i, kick2i, tau23i,
                    kick3i, tau34i, kick4i, tau41i
                ] = vals2
                f06_file.write(
                    '0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n'
                    '                     %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, f14i, f12i, f21i, f23i, f32i, f34i, f43i, f41i,
                            kick1i, tau12i, kick2i, tau23i, kick3i, tau34i, kick4i, tau41i))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealViscForceArray(RealForceObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if not is_sort1:
            #raise NotImplementedError('SORT2')

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def get_headers(self) -> List[str]:
        headers = ['axial', 'torsion']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the RealViscForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                          1             2             3
            #Freq               1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue        -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians            9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #50        axial            -0.0          -0.0           0.0
            #          torsion           0.0           0.0          -0.0
            #51        axial             0.0          -0.0          -0.0
            #          torsion          -0.0           0.0           0.0
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial  torsion
            #ElementID
            #14           0.0      0.0
            #15           0.0      0.0
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def add_sort1(self, dt, eid, axial, torque):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
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
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if is_sort1:
            msg = [
                '                                   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )\n'
                ' \n'
                '       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n'
                '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n'
            ]
        else:
            msg = [
                '                                   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )\n'
                ' \n'
                '                         AXIAL                                                       AXIAL\n'
                '       TIME              FORCE         TORQUE                      TIME              FORCE         TORQUE\n'
                #'   0.0                0.0            0.0                       1.000000E+00      -5.642718E-04   0.0\n'
                #'   2.000000E+00      -1.905584E-06   0.0                       3.000000E+00       9.472010E-07   0.0\n'
            ]
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [axiali, torsioni] = write_floats_13e([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %s\n' % tuple(out[i] + out[i + 1])
                f06_file.write(out_line)
            if is_odd:
                out_line = '      %8i   %-13s  %s\n' % tuple(out[-1])
                f06_file.write(out_line)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid1, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealPlateForceArray(RealForceObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0
        assert self.element_name != 'RBAR', self.data_code

        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self) -> List[str]:
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        """sizes the vectorized attributes of the RealPlateForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype=idtype)

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype=fdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                       1             2             3
            #Freq            1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue     -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians         9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #8         mx   -5.467631e-14 -1.406068e-07  1.351960e-13
            #          my   -8.983144e-14 -3.912936e-07  9.707208e-14
            #          mxy   2.767353e-13 -4.950616e-08 -5.985472e-13
            #          bmx   7.616284e-14 -2.809588e-08 -1.051987e-13
            #          bmy   4.245138e-14 -6.567249e-09 -6.066584e-14
            #          bmxy -1.233790e-14  3.561397e-09  1.840837e-14
            #          tx    2.601638e-13 -9.601510e-08 -3.611116e-13
            #          ty   -5.825233e-14 -7.382687e-09  9.038553e-14
            #9         mx    5.444685e-15 -1.014145e-07 -4.500100e-14
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1

        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    #def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        #self._times[self.itime] = dt
        #self.element[self.ielement] = eid
        #self.data[self.itime, self.ielement, :] = [axial, SMa, torsion, SMt]
        #self.ielement += 1

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        #self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

    @property
    def nnodes_per_element(self):
        return 1

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
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        if 'CTRIA3' in self.element_name:
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                ' \n'
                '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            ]
            nnodes = 3
        elif 'CQUAD4' in self.element_name:
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        elif 'CTRIAR' in self.element_name:
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'
                ' \n'
                '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            ]
            nnodes = 3
        elif 'CQUADR' in self.element_name:
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            msg = f'element_name={self.element_name} self.element_type={self.element_type}'
            raise NotImplementedError(msg)
        return self.element_name, nnodes, msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element
        cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            if self.element_type in [74, 83, 227, 228]:
                # 74, 83 CTRIA3
                # 227 CTRIAR linear
                # 228 CQUADR linear
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                        [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # ctria3
                    #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                    f06_file.write('   %8i %18s %13s %13s   %13s %13s %13s   %13s %s\n' % (
                        eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))

            elif self.element_type == 33:
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                        [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # cquad4
                    #0         6    CEN/4  1.072685E+01  2.504399E+03 -2.455727E+01 -5.017930E+00 -2.081427E+01 -5.902618E-01 -9.126162E+00  4.194400E+01#
                    #Fmt = '% 8i   ' + '%27.20E   ' * 8 + '\n'
                    #f06_file.write(Fmt % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                    #
                    f06_file.write('0 %8i %8s %13s %13s %13s %13s %13s %13s %13s %s\n' % (
                        eid, cen_word, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            else:
                raise NotImplementedError(f'element_name={self.element_name} element_type={self.element_type}')
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        if 'CTRIA3' in self.element_name:
            nnodes = 3
        elif 'CQUAD4' in self.element_name:
            nnodes = 4
        elif 'CTRIAR' in self.element_name:
            nnodes = 4 # ???
        elif 'CQUADR' in self.element_name:
            nnodes = 5 # ???
        else:  # pragma: no cover
            raise NotImplementedError(self.code_information())

        #print("nnodes_all =", nnodes_all)
        #cen_word_ascii = 'CEN/%i' % nnodes
        #cen_word = b'CEN/%i' % nnodes

        eids = self.element
        #cen_word = 'CEN/%i' % nnodes

        #msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        #msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))

        eids = self.element
        eids_device = eids * 10 + self.device_code

        nelements = len(eids)
        assert nelements > 0, eids
        #print('nelements =', nelements)
        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]

        if self.is_sort1:
            structi = Struct(endian + b'i 8f')
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

            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            nwide = 0
            for eid_device, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                data = [eid_device, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi]
                #[mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                #    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                op2.write(structi.pack(*data))
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data[1:])))
                nwide += len(data)

            assert nwide == ntotal, "nwide=%s ntotal=%s" % (nwide, ntotal)
            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealPlateBilinearForceArray(RealForceObject):  # 144-CQUAD4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self) -> List[str]:
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        """sizes the vectorized attributes of the RealPlateBilinearForceArray"""
         #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        # -MEMBRANE FORCES-   -BENDING MOMENTS- -TRANSVERSE SHEAR FORCES -
        #     FX FY FXY           MX MY MXY            QX QY
        #[fx, fy, fxy,  mx,  my,  mxy, qx, qy]
        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        node = pd.Series(self.element_node[:, 1])
        node.replace({'NodeID': {0:'CEN'}}, inplace=True)
        element_node = [self.element_node[:, 0], node]

        if self.nonlinear_factor not in (None, np.nan):
            # Mode                              1             2             3
            # Freq                   1.482246e-10  3.353940e-09  1.482246e-10
            # Eigenvalue            -8.673617e-19  4.440892e-16  8.673617e-19
            # Radians                9.313226e-10  2.107342e-08  9.313226e-10
            # ElementID NodeID Item
            # 6         0      mx    2.515537e-13 -2.294306e-07 -3.626725e-13
            #                  my    2.916815e-13 -7.220319e-08 -5.030049e-13
            #                  mxy  -2.356622e-14  4.391171e-07 -5.960345e-14
            #                  bmx   4.138377e-14 -1.861012e-08 -5.586283e-14
            #                  bmy   5.991298e-15 -2.471926e-09 -5.400710e-15
            #                  bmxy  4.511364e-15 -1.190845e-09 -5.546569e-15
            #                  tx    1.122732e-13 -5.563460e-08 -1.523176e-13
            #                  ty   -1.164320e-14  4.813929e-09  1.023404e-14
            # 4                mx    3.839208e-13 -4.580973e-07 -4.949736e-13
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, element_node, self.data, from_tuples=False, from_array=True)
        else:
            df1 = pd.DataFrame(element_node).T
            df1.columns = ['ElementID', 'NodeID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            data_frame = df1.join(df2)
            data_frame = data_frame.reset_index().set_index(['ElementID', 'NodeID'])
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                            mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element_node[self.itotal] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        #self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

    @property
    def nnodes_per_element(self):
        if self.element_type == 144:  # CQUAD4
            nnodes_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_element = 5
        elif self.element_type == 75:  # CTRIA6
            nnodes_element = 4
        elif self.element_type == 70:  # CTRIAR
            nnodes_element = 4
        else:
            raise NotImplementedError('element_type=%s element_name=%s' % (self.element_type, self.element_name))
        return nnodes_element

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i ntotal=%i nnodes/element=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, ntotal, self.nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i ntotal=%i nnodes/element=%i\n'
                       % (self.__class__.__name__, nelements, ntotal, self.nnodes_per_element))
            ntimes_word = '1'

        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        # if 'CTRIA3' in self.element_name:
            # msg = [
                # '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                # ' \n'
                # '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                # '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            # ]
            # nnodes = 3

        if self.element_type == 70:
            # CQUAD4
            element_name = 'CTRIAR'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 75:
            # CQUAD4
            element_name = 'CTRIA6'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 64:
            # CQUAD4
            element_name = 'CQUAD8'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 8
        elif self.element_type == 82:
            # CQUAD4
            element_name = 'CQUADR'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        elif self.element_type == 144:
            # CQUAD4
            element_name = 'CQUAD4'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        return element_name, nnodes, msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        cen_word = 'CEN/%i' % nnodes
        if self.element_type  in [64, 82, 144]: # CQUAD8, CQUADR, CQUAD4
            cyci = [0, 1, 2, 3, 4]
            #cyc = cycle([0, 1, 2, 3, 4])  # TODO: this is totally broken...
            nnodes_per_eid = 5
        elif self.element_type  in [70, 75]: # CTRIAR, CTRIA6
            cyci = [0, 1, 2, 3]
            #cyc = cycle([0, 1, 2, 3])  # TODO: this is totally broken...
            nnodes_per_eid = 4
        else:
            raise NotImplementedError(self.element_type)

        # TODO: this shouldn't be neccessary
        cyc = cyci * (len(eids) // nnodes_per_eid)
        assert len(eids) % nnodes_per_eid == 0

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for i, eid, nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(cyc, eids, nids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                # ctria3
                #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                if i == 0:
                    f06_file.write(
                        '0  %8i    %s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, cen_word, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                else:
                    f06_file.write(
                        '            %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            # else:
                # raise NotImplementedError(self.element_type)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        (unused_elem_name, nnodes, unused_msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        #ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        if self.element_type  in [64, 82, 144]: # CQUAD8, CQUADR, CQUAD4
            cyci = [0, 1, 2, 3, 4]
            #cyc = cycle([0, 1, 2, 3, 4])  # TODO: this is totally broken...
            nnodes_per_eid = 5
        elif self.element_type  in [70, 75]: # CTRIAR, CTRIA6
            cyci = [0, 1, 2, 3]
            #cyc = cycle([0, 1, 2, 3])  # TODO: this is totally broken...
            nnodes_per_eid = 4
        else:
            raise NotImplementedError(self.element_type)

        # TODO: this shouldn't be neccessary
        cyc = cyci * (len(eids) // nnodes_per_eid)
        assert len(eids) % nnodes_per_eid == 0

        #print("nnodes_all =", nnodes_all)
        #cen_word_ascii = 'CEN/%i' % nnodes
        cen_word = b'CEN/'

        #msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        #msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        eids_device = eids * 10 + self.device_code

        nelements = len(np.unique(eids))
        #print('nelements =', nelements)
        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]

        if self.is_sort1:
            struct1 = Struct(endian + b'i4s i 8f')
            struct2 = Struct(endian + b'i 8f')
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

            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            nwide = 0
            for i, eid, eid_device, nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(cyc, eids, eids_device, nids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                #[mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_floats_13e(
                #    [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])

                if i == 0:
                    data = [eid_device, cen_word, nnodes, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi]
                    op2.write(struct1.pack(*data))
                    op2_ascii.write(
                        '0  %8i    %s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            eid, cen_word, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                else:
                    data = [nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi]
                    op2.write(struct2.pack(*data))
                    op2_ascii.write(
                        '            %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                            nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))

                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                nwide += len(data)

            assert nwide == ntotal, "nwide=%s ntotal=%s" % (nwide, ntotal)
            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealCBarFastForceArray(RealForceObject):
    """
    34-CBAR
    119-CFAST

    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def get_headers(self) -> List[str]:
        headers = [
            'bending_moment_a1', 'bending_moment_a2',
            'bending_moment_b1', 'bending_moment_b2',
            'shear1', 'shear2',
            'axial', 'torque']
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealCBarForceArray"""
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
        self.is_built = True
        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        if self.is_sort1:
            ntimes = self.ntimes
            nelements = self.nelements
            ntotal = self.ntotal
            self._build(ntimes, nelements, ntotal, self._times_dtype)
        else:
            ntimes = self.nelements
            nelements = self.ntimes
            ntotal = nelements * 2
            name = self.analysis_method + 's'
            self._build(ntimes, nelements, ntotal, self._times_dtype)
            setattr(self, name, self._times)
            self.data_code['name'] = self.analysis_method
            self.data_names[0] = self.analysis_method
            #print(f'data_names -> {self.data_names}')

    def _build(self, ntimes, nelements, ntotal, dtype):
        self.ntimes = ntimes
        self.nelements = nelements
        self.ntotal = ntotal
        #print(f"*ntimes={ntimes} nelements={nelements} ntotal={ntotal} data_names={self.data_names}")
        unused_dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(ntimes, dtype=dtype)
        self.element = zeros(nelements, dtype=idtype)

        #[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
        self.data = zeros((ntimes, ntotal, 8), dtype=fdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def add_sort1(self, dt, eid, bending_moment_a1, bending_moment_a2,
                  bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #[eid, bending_moment_a1, bending_moment_a2,
         #bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque] = data
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            bending_moment_a1, bending_moment_a2,
            bending_moment_b1, bending_moment_b2,
            shear1, shear2, axial, torque]
        self.ielement += 1

    def add_sort2(self, dt, eid, bending_moment_a1, bending_moment_a2,
                  bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque):
        """unvectorized method for adding SORT2 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #[eid, bending_moment_a1, bending_moment_a2,
         #bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque] = data
        itime = self.ielement
        ielement = self.itime
        #print(f'{self.table_name} itime={itime} ielement={ielement} time={dt} eid={eid} axial={axial}')
        self._times[itime] = dt
        self.element[ielement] = eid
        self.data[itime, ielement, :] = [
            bending_moment_a1, bending_moment_a2,
            bending_moment_b1, bending_moment_b2,
            shear1, shear2, axial, torque]
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
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = self._words()

        #msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f06_file.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + words))

            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(
                    eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                [bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi] = write_floats_13e([
                 bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi])
            f06_file.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (
                eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi))
            f06_file.write(page_stamp % page_num)
        return page_num

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1) = t1
                    (bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1,
                            bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
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

        eids = self.element
        eids_device = eids * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

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
            struct1 = Struct(endian + b'i 8f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            for eid_device, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(
                    eids_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):

                data = [eid_device, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    @abstractmethod
    def _words(self) -> List[str]:
        return []

class RealCBarForceArray(RealCBarFastForceArray):  # 34-CBAR
    """34-CBAR"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCBarFastForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @classmethod
    def add_static_case(cls, table_name, element, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):

        analysis_code = 1 # static
        data_code = oef_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['loadIDs'] = [0] # TODO: ???
        data_code['data_names'] = []

        # I'm only sure about the 1s in the strains and the
        # corresponding 0s in the stresses.
        #if is_stress:
            #data_code['stress_bits'] = [0, 0, 0, 0]
            #data_code['s_code'] = 0
        #else:
            #data_code['stress_bits'] = [0, 1, 0, 1]
            #data_code['s_code'] = 1 # strain?

        data_code['element_name'] = 'CBAR'
        data_code['element_type'] = 34
        #data_code['load_set'] = 1

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.element = element
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    def _words(self) -> List[str]:
        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        return words

class RealCWeldForceArray(RealCBarFastForceArray):  # 34-CBAR
    """117-CWELD"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCBarFastForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

class RealCFastForceArrayNX(RealCBarFastForceArray):  # 34-CBAR
    """119-CFAST"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCBarFastForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

class RealConeAxForceArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if not is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'hopa', 'bmu', 'bmv', 'tm', 'su', 'sv'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the RealConeAxForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[hopa, bmu, bmv, tm, su, sv]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element)
            df1.columns = ['ElementID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #print(self.data_frame)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (hopa1, bmu1, bmv1, tm1, su1, sv1) = t1
                    (hopa2, bmu2, bmv2, tm2, su2, sv2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                hopa1, bmu1, bmv1, tm1, su1, sv1,
                                hopa2, bmu2, bmv2, tm2, su2, sv2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, hopa, bmu, bmv, tm, su, sv):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [hopa, bmu, bmv, tm, su, sv]
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
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '             F O R C E S   I N   A X I S - S Y M M E T R I C   C O N I C A L   S H E L L   E L E M E N T S   (CCONEAX)\n'
            ' \n'
            '  ELEMENT     HARMONIC    POINT           BEND-MOMENT       BEND-MOMENT      TWIST-MOMENT           SHEAR            SHEAR\n'
            '   ID.         NUMBER     ANGLE               V                 U                                     V                U\n'
            #'      101        0                       5.864739E-09      1.759422E-09      0.0               0.0               0.0'
            #'      101                0.0000          5.864739E-09      1.759422E-09      0.0               0.0               0.0'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            hopa = self.data[itime, :, 0]
            bmu = self.data[itime, :, 1]
            bmv = self.data[itime, :, 2]
            tm = self.data[itime, :, 3]
            su = self.data[itime, :, 4]
            sv = self.data[itime, :, 5]

            for (eid, hopai, bmui, bmvi, tmi, sui, svi) in zip(
                    eids, hopa, bmu, bmv, tm, su, sv):
                if hopai > 0.1:
                    raise NotImplementedError(hopai)

                vals2 = write_floats_13e([hopai, bmui, bmvi, tmi, sui, svi])
                [hopai, bmui, bmvi, tmi, sui, svi] = vals2

                # TODO: hopa is probably the wrong type
                              # hopa        # Mu       Mv      twist       Vy       Vu
                f06_file.write(' %8i  %-13s  %-13s %-13s     %-13s     %-13s     %-13s     %s\n' % (
                    eid, 0.0, hopai, bmui, bmvi, tmi, sui, svi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBar100ForceArray(RealForceObject):  # 100-CBAR
    """
    CBAR-34s are converted to CBAR-100s when you have PLOAD1s
    (distributed bar loads).  The number of stations by default is 2,
    but with a CBARAO, you can change this (max of 8 points; 6 internal
    points).

    If you use a CBARO without PLOAD1s, you wil turn CBAR-34s into
    CBAR-100s as well.
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if not is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    def get_headers(self) -> List[str]:
        headers = [
            'station', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2', 'axial', 'torque'
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealCBar100ForceArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    #def finalize(self):
        #sd = self.data[0, :, 0]
        #i_sd_zero = np.where(sd != 0.0)[0]
        #i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        #assert i_node_zero.max() > 0, 'CBAR element_node hasnt been filled'
        #i = np.union1d(i_sd_zero, i_node_zero)
        #self.element = self.element[i]
        #self.element_node = self.element_node[i, :]
        #self.data = self.data[:, i, :]

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        element_location = [
            self.element,
            self.data[0, :, 0],
        ]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data[:, :, 1:], items=column_values,
                                       major_axis=element_location, minor_axis=headers[1:]).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Location', 'Item']
        else:
            df1 = pd.DataFrame(element_location).T
            df1.columns = ['ElementID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID'])
        #print(self.data_frame)

    def add_sort1(self, dt, eid, sd, bm1, bm2, ts1, ts2, af, trq):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid

        # station, bending_moment1, bending_moment2, shear1, shear2, axial, torque
        self.data[self.itime, self.ielement, :] = [sd, bm1, bm2, ts1, ts2, af, trq]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            msg = [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
            return msg

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #'                         F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )'
        #'0    ELEMENT  STATION         BEND-MOMENT                      SHEAR FORCE                     AXIAL'
        #'       ID.     (PCT)     PLANE 1        PLANE 2           PLANE 1        PLANE 2               FORCE              TORQUE'
        #'           10   0.000   0.0           -5.982597E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           10   1.000   0.0            1.868857E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.000   0.0            1.868857E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.050   0.0            2.261429E+06      0.0           -7.851454E+03          0.0                0.0'
        #'           11   0.100   0.0            2.654002E+06      0.0           -7.851454E+03          0.0                0.0'
        words = [
            '                         F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )\n'
            '0    ELEMENT  STATION         BEND-MOMENT                      SHEAR FORCE                     AXIAL\n'
            '       ID.     (PCT)     PLANE 1        PLANE 2           PLANE 1        PLANE 2               FORCE              TORQUE\n']
            # '        15893   0.000   1.998833E+02   9.004551E+01      2.316835E+00   1.461960E+00         -2.662207E+03       9.795244E-02'
        #msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + words))

            # sd, bm1, bm2, ts1, ts2, af, trq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            trq = self.data[itime, :, 6]
            for eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, trq):
                [bm1i, bm2i, ts1i, ts2i, afi, trqi] = write_floats_13e([
                    bm1i, bm2i, ts1i, ts2i, afi, trqi])
                f06_file.write(
                    '     %8i   %4.3f  %-13s  %-13s     %-13s  %-13s         %-13s      %s\n' % (
                        eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi))
            f06_file.write(page_stamp % page_num)
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
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

        eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i2f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
        eids = self.element
        #f.write(''.join(words))

        #ntimes = self.data.shape[0]
        struct1 = Struct(endian + b'i7f')
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

            # sd, bm1, bm2, ts1, ts2, af, trq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            trq = self.data[itime, :, 6]
            for eid, eid_device, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi in zip(
                    eids, eids_device, sd, bm1, bm2, ts1, ts2, af, trq):
                [sbm1i, sbm2i, sts1i, sts2i, safi, strqi] = write_floats_13e([
                    bm1i, bm2i, ts1i, ts2i, afi, trqi])
                op2_ascii.write(
                    '     %8i   %4.3f  %-13s  %-13s     %-13s  %-13s         %-13s      %s\n' % (
                        eid, sdi, sbm1i, sbm2i, sts1i, sts2i, safi, strqi))
                data = [eid_device, sdi, bm1i, bm2i, ts1i, ts2i, afi, trqi]
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (sd1, bm11, bm21, ts11, ts21, af1, trq1) = t1
                    (sd2, bm12, bm22, ts12, ts22, af2, trq2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            sd1, bm11, bm21, ts11, ts21, af1, trq1,
                            sd2, bm12, bm22, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealCGapForceArray(RealForceObject):  # 38-CGAP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if not is_sort1:
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'fx', 'sfy', 'sfz', 'u', 'v', 'w', 'sv', 'sw'
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealCGapForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        # [fx, sfy, sfz, u, v, w, sv, sw]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        # LoadStep                 1.0
        # ElementID Item
        # 101       fx    33333.332031
        #           sfy       0.000000
        #           sfz       0.000000
        #           u         0.000115
        #           v         0.000000
        #           w         0.000000
        #           sv        0.000000
        #           sw        0.000000
        # 102       fx       -0.000002
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            # Static               fx  sfy  sfz         u         v    w   sv   sw
            # ElementID
            # 1          1.253610e-10 -0.0  0.0  0.250722 -0.852163  0.0  0.0  0.0
            # 21         1.253610e-10 -0.0  0.0  0.250722 -0.852163  0.0  0.0  0.0
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fx1, sfy1, sfz1, u1, v1, w1, sv1, sw1) = t1
                    (fx2, sfy2, sfz2, u2, v2, w2, sv2, sw2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fx1, sfy1, sfz1, u1, v1, w1, sv1, sw1,
                                fx2, sfy2, sfz2, u2, v2, w2, sv2, sw2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, fx, sfy, sfz, u, v, w, sv, sw):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [fx, sfy, sfz, u, v, w, sv, sw]
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
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                                    F O R C E S   I N   G A P   E L E M E N T S   ( C G A P )\n'
            '  ELEMENT   - F O R C E S  I N  E L E M  S Y S T E M -         - D I S P L A C E M E N T S  I N  E L E M  S Y S T E M -\n'
            '    ID       COMP-X         SHEAR-Y        SHEAR-Z        AXIAL-U        TOTAL-V        TOTAL-W        SLIP-V         SLIP-W\n'
            #'     101   3.333333E+04   0.0            0.0            1.149425E-04   0.0            0.0            0.0            0.0\n'
        ]

        ##(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #ntimes, ntotal = self.data.shape[:1]
        ntimes = self.data.shape[0]

        eids = self.element

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            # [fx, sfy, sfz, u, v, w, sv, sw]
            fx = self.data[itime, :, 0]
            sfy = self.data[itime, :, 1]
            sfz = self.data[itime, :, 2]
            u = self.data[itime, :, 3]
            v = self.data[itime, :, 4]
            w = self.data[itime, :, 5]
            sv = self.data[itime, :, 6]
            sw = self.data[itime, :, 7]

            for (eid, fxi, sfyi, sfzi, ui, vi, wi, svi, swi) in zip(eids, fx, sfy, sfz, u, v, w, sv, sw):
                vals2 = write_floats_12e([fxi, sfyi, sfzi, ui, vi, wi, svi, swi])
                [fxi, sfyi, sfzi, ui, vi, wi, svi, swi] = vals2
                f06_file.write('0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (
                    eid, fxi, sfyi, sfzi, ui, vi, wi, svi, swi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealBendForceArray(RealForceObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)
        self.nelements = 0  # result specific

    def get_headers(self) -> List[str]:
        headers = [
            'bending_moment_1a', 'bending_moment_2a', 'shear_1a', 'shear_2a', 'axial_a', 'torque_a',
            'bending_moment_1b', 'bending_moment_2b', 'shear_1b', 'shear_2b', 'axial_b', 'torque_b',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealBendForceArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.nelements, 3), dtype='int32')

        #[bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a
        # bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b]
        self.data = zeros((self.ntimes, self.nelements, 12), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            # TODO: add NodeA, NodeB
            element = self.element_node[:, 0]
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, element, self.data)
            #data_frame = pd.Panel(self.data, items=column_values,
                                  #major_axis=element, minor_axis=headers).to_frame()
            #data_frame.columns.names = column_names
            #data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element_node)
            df1.columns = ['ElementID', 'NodeA', 'NodeB']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            data_frame = df1.join(df2)
        self.data_frame = data_frame

    def add_sort1(self, dt, eid,
                  nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                  nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)

        self._times[self.itime] = dt
        self.element_node[self.ielement] = [eid, nid_a, nid_b]
        self.data[self.itime, self.ielement, :] = [
            bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
            bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b,
        ]
        self.ielement += 1
        if self.ielement == self.nelements:
            self.ielement = 0

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            '                                 F O R C E S   I N   B E N D   E L E M E N T S        ( C B E N D )\n'
            '                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL\n'
            '   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
            #'0      6901    6901      A     0.0           0.0            0.0           0.0            0.0           -6.305720E-16'
            #'               6902      B    -5.000000E-01  5.000000E-01   1.000000E+00 -1.000000E+00  -5.000000E-07  -1.666537E-07'
        ]
        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nid_a = self.element_node[:, 1]
        nid_b = self.element_node[:, 2]

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            bending_moment_1a = self.data[itime, :, 0]
            bending_moment_2a = self.data[itime, :, 1]
            shear_1a = self.data[itime, :, 2]
            shear_2a = self.data[itime, :, 3]
            axial_a = self.data[itime, :, 4]
            torque_a = self.data[itime, :, 5]
            bending_moment_1b = self.data[itime, :, 6]
            bending_moment_2b = self.data[itime, :, 7]
            shear_1b = self.data[itime, :, 8]
            shear_2b = self.data[itime, :, 9]
            axial_b = self.data[itime, :, 10]
            torque_b = self.data[itime, :, 11]


            for (eid,
                 nid_ai, bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                 nid_bi, bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi) in zip(
                     eids,
                     nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                     nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
                [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                 bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi] = write_floats_13e(
                     [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                      bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi])
                f06_file.write(
                    '0  %8i%8i      A    %13s %13s  %13s %13s  %13s  %13s\n'
                    '           %8i      B    %13s %13s  %13s %13s  %13s  %13s\n' % (
                        eid, nid_ai, bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                        nid_bi, bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi
                    ))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'element_node shape=%s table.shape=%s' % (self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Nid_A, Nid_B\n'
            for (eid1, nida1, nidb1), (eid2, nida2, nidb2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s, %s), (%s, %s, %s)\n' % (eid1, nida1, nidb1, eid2, nida2, nidb2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            eids = self.element_node[:, 0]
            for itime in range(self.ntimes):
                for ie, eid in enumerate(eids):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]

                    (bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                     bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1) = t1
                    (bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                     bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2) = t2

                    if not np.array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                            bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1,

                            bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                            bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealSolidPressureForceArray(RealForceObject):  # 77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if not is_sort1:
            #raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'ax', 'ay', 'az', 'vx', 'vy', 'vz', 'pressure'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the RealSolidPressureForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[ax, ay, az, vx, vy, vz, pressure]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype=fdtype)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (ax1, ay1, az1, vx1, vy1, vz1, pressure1) = t1
                    (ax2, ay2, az2, vx2, vy2, vz2, pressure2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                ax1, ay1, az1, vx1, vy1, vz1, pressure1,
                                ax2, ay2, az2, vx2, vy2, vz2, pressure2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, etype, ax, ay, az, vx, vy, vz, pressure):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [ax, ay, az, vx, vy, vz, pressure]
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
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape

        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file)
        else:
            raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())
        return page_num

    def _write_sort2_as_sort1(self, header, page_stamp, page_num=1, f=None):
        msg_temp = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                    ' \n',
                    '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        ntimes = self.data.shape[0]

        eids = self.element
        etype = self.element_name
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 1]
            vz = self.data[itime, :, 2]
            ax = self.data[itime, :, 3]
            ay = self.data[itime, :, 4]
            az = self.data[itime, :, 5]
            pressure = self.data[itime, :, 5]

            for (eid, vxi, vyi, vzi, axi, ayi, azi, pressurei) in zip(
                    eids, vx, vy, vz, ax, ay, az, pressure):

                vals2 = write_floats_13e([axi, ayi, azi, pressurei])
                [sax, say, saz, spressure] = vals2

                #etype = 'PENPR'
                f.write('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, etype, sax, say, saz, spressure))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_sort1_as_sort1(self, header, page_stamp, page_num=1, f=None):
        msg_temp = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                    ' \n',
                    '    ELEMENT-ID   EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']  # TODO: bad line...
        ntimes = self.data.shape[0]

        eids = self.element
        etype = self.element_name
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 1]
            vz = self.data[itime, :, 2]
            ax = self.data[itime, :, 3]
            ay = self.data[itime, :, 4]
            az = self.data[itime, :, 5]
            pressure = self.data[itime, :, 5]

            for (eid, vxi, vyi, vzi, axi, ayi, azi, pressurei) in zip(
                    eids, vx, vy, vz, ax, ay, az, pressure):

                vals2 = write_floats_13e([axi, ayi, azi, pressurei])
                [sax, say, saz, spressure] = vals2

                #etype = 'PENPR'
                f.write('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, etype, sax, say, saz, spressure))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


# F:\work\pyNastran\examples\Dropbox\move_tpl\beamp11.op2
class RealCBeamForceVUArray(RealForceObject):  # 191-VUBEAM
    """
    **ELTYPE = 191 Beam view element (VUBEAM)**

    2 PARENT I     Parent p-element identification number
    3 COORD  I     Coordinate system identification number
    4 ICORD  CHAR4 Flat/curved and so on

    TCODE,7 = 0 Real
    5 VUGRID   I  VU grid ID for output grid
    6 POSIT    RS x/L position of VU grid identification number
    7 FORCEX   RS Force x
    8 SHEARY   RS Shear force y
    9 SHEARZ   RS Shear force z
    10 TORSION RS Torsional moment x
    11 BENDY   RS Bending moment y
    12 BENDZ   RS Bending moment z


          DIRECT TRANSIENT RESPONSE                                           ADAPTIVITY INDEX=      1
    0                                                                         PVAL ID=       1                       SUBCASE=       1
       VU-ELEMENT ID =  100001002

                               F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )
                        TIME =   0.000000E+00,  P-ELEMENT ID =       1,  OUTPUT COORD. ID =       0,  P OF EDGES =  1

          VUGRID VUGRID DIST/     - BENDING MOMENTS -              -WEB  SHEARS -               AXIAL           TOTAL
            ID.     LENGTH       PLANE 1       PLANE 2          PLANE 1       PLANE 2            FORCE           TORQUE
       111001002     0.333     0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     0.000000E+00     0.000000E+00
       111001003     0.667     0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     0.000000E+00     0.000000E+00

                               F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )
                        TIME =   1.000000E+00,  P-ELEMENT ID =       1,  OUTPUT COORD. ID =       0,  P OF EDGES =  1

          VUGRID VUGRID DIST/     - BENDING MOMENTS -              -WEB  SHEARS -               AXIAL           TOTAL
            ID.     LENGTH       PLANE 1       PLANE 2          PLANE 1       PLANE 2            FORCE           TORQUE
       111001002     0.333     0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     9.982032E-01     0.000000E+00
       111001003     0.667     0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     9.982032E-01     0.000000E+00
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = ['xxb', 'fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealCBeamForceVUArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')
        self.parent_coord = zeros((self.ntotal, 2), dtype='int32')

        #[xxb, fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                                1             2
            #Freq                      1.214849e-07  1.169559e-07
            #Eigenvalue               -5.826450e-13 -5.400125e-13
            #Radians                   7.633119e-07  7.348554e-07
            #ElementID NodeID    Item
            #100001001 111001001 xxb   0.000000e+00  0.000000e+00
            #                    fx   -2.363981e-14  1.091556e-15
            #                    fy    1.041715e-13  5.642625e-14
            #                    fz    2.040026e-14  1.133024e-12
            #                    mx   -1.903338e-15  5.201282e-16
            #                    my    8.236364e-14  2.141288e-13
            #                    mz   -2.247944e-14  9.947491e-14
            #          111001002 xxb   3.333333e-01  3.333333e-01
            #                    fx   -2.278856e-14  5.621586e-16
            #                    fy   -2.092285e-14 -7.903003e-14
            #                    fz   -6.230423e-14  2.318665e-12
            #                    mx    7.229356e-16 -1.005680e-16
            #                    my    3.688462e-14  7.122293e-14
            #                    mz   -7.785338e-15  3.304847e-14
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, self.element_node, self.data)
        else:
            data_frame = pd.Panel(self.data,
                                  major_axis=self.element_node[:, 0], minor_axis=headers).to_frame()
            data_frame.columns.names = ['Static']
            data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor not in (None, np.nan):
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for (eid, nid), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s) (%s, %s)\n' % (eid, nid, eid2, nid2)
            print(msg)
            raise ValueError(msg)
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
                        (fx1, fy1, fz1, mx1, my1, mz1) = t1
                        (fx2, fy2, fz2, mx2, my2, mz2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fx1, fy1, fz1, mx1, my1, mz1,
                                fx2, fy2, fz2, mx2, my2, mz2)
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

    def _add_sort1(self, dt, eid, parent, coord, icord,
                   nid, xxb, fx, fy, fz, mx, my, mz):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt

        #eids = self.element_node[:, 0]
        #nids = self.element_node[:, 1]
        self.element_node[self.ielement] = [eid, nid]
        self.parent_coord[self.ielement] = [parent, coord]
        self.data[self.itime, self.ielement, :] = [xxb, fx, fy, fz, mx, my, mz]
        self.ielement += 1

    #def add_sort1(self, nnodes, dt, data):
        #[eid, parent, coord, icord, forces] = data
        #self.parent[eid] = parent
        #self.coord[eid] = coord
        #self.icord[eid] = icord

        #self.forceX[eid] = {}
        #self.shearY[eid] = {}
        #self.shearZ[eid] = {}
        #self.torsion[eid] = {}
        #self.bendingY[eid] = {}
        #self.bendingZ[eid] = {}


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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  parent_coord.shape = %s\n' % str(self.parent_coord.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def get_f06_header(self):
        msg = [
            '                                 F O R C E S   I N   VU ELEMENTS\n'
            ' \n'
            '                  ELEMENT-ID        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n']
           #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        NotImplementedError(self.code_information())

        msg_temp = self.get_f06_header()

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            for eid, nid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, nids, fx, fy, fz, mx, my, mz):
                [fxi, fyi, fzi, mxi, myi, mzi] = write_floats_13e([fxi, fyi, fzi, mxi, myi, mzi])
                f06_file.write('                    %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (
                    eid, fxi, fyi, fzi, mxi, myi, mzi))
                #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealForceMomentArray(RealForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if not is_sort1:
            #raise NotImplementedError('SORT2; code_info=\n%s' % self.code_information())

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealCBushForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
            #data_frame = pd.Panel(self.data, items=column_values,
                                  #major_axis=self.element, minor_axis=headers).to_frame()
            #data_frame.columns.names = column_names
            #data_frame.index.names = ['ElementID', 'Item']
        else:
            # >=25.0
            #Static         fx   fy   fz   mx   my   mz
            #ElementID
            #1          1000.0  0.0  0.0  0.0  0.0  0.0
            #
            # <=24.2
            #Static               0
            #ElementID Item
            #1         fx    1000.0
            #          fy       0.0
            #          fz       0.0
            #          mx       0.0
            #          my       0.0
            #          mz       0.0
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
            #data_frame = pd.Panel(self.data,
                                  #major_axis=self.element, minor_axis=headers).to_frame()
            #data_frame.columns.names = ['Static']
            #data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (
            self.nonlinear_factor is not None and
            np.isnan(self.nonlinear_factor) and
            np.isnan(table.nonlinear_factor)
        )
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor not in (None, np.nan):
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
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
                        (fx1, fy1, fz1, mx1, my1, mz1) = t1
                        (fx2, fy2, fz2, mx2, my2, mz2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fx1, fy1, fz1, mx1, my1, mz1,
                                fx2, fy2, fz2, mx2, my2, mz2)
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

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [fx, fy, fz, mx, my, mz]
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
        #ntotal = self.ntotal

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
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header()

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                [fxi, fyi, fzi, mxi, myi, mzi] = write_floats_13e([fxi, fyi, fzi, mxi, myi, mzi])
                f06_file.write('                    %8i     %-13s %-13s %-13s %-13s %-13s %s\n' % (
                    eid, fxi, fyi, fzi, mxi, myi, mzi))
                #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
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

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i6f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids_device, fx, fy, fz, mx, my, mz):
                data = [eid, fxi, fyi, fzi, mxi, myi, mzi]
                op2_ascii.write('  eid=%s data=%s\n' % (eid, str(data[1:])))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable



class RealCBushForceArray(RealForceMomentArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceMomentArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def get_f06_header(self):
        msg = [
            '                                 F O R C E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n'
            ' \n'
            '                  ELEMENT-ID        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n']
           #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
        return msg

class RealCBearForceArray(RealForceMomentArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceMomentArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def get_f06_header(self):
        if self.element_type == 0: # 'CBUSH':

            msg = [
                '                                 F O R C E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n'
                ' \n'
                '                  ELEMENT-ID        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n']
               #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
        elif self.element_type == 280: # 'CBEAR':
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\rotbr60f.op2
            msg = [
                '                              F O R C E S   I N   B E A R I N G   E L E M E N T S        ( C B E A R )\n'
                ' \n'
                '                  ELEMENT-ID        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n']
               #'0                        599      0.0           2.000000E+00  3.421458E-14  1.367133E-13 -3.752247E-15  1.000000E+00\n']
        else:
            msg = f'element_name={self.element_name} self.element_type={self.element_type}'
            raise NotImplementedError(msg)
        return msg

class RealCFastForceArrayMSC(RealForceMomentArray):
    """126-CFAST-MSC"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceMomentArray.__init__(self, data_code, is_sort1, isubcase, dt)

class RealForceVU2DArray(RealForceObject):  # 189-VUQUAD, 190-VUTRIA
    """
                    F O R C E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
                 TIME =   2.500000E-03,  P-ELEMENT ID =      11,  OUTPUT COORD. ID =       0,  P OF EDGES =  3  3  3
                       LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000

     VUGRID                - MEMBRANE  FORCES -      - BENDING  MOMENTS -             - TRANSVERSE SHEAR FORCES -
       ID.            FX            FY            FXY MX            MY            MXY             QX            QY
    111001001     1.497761E+02  4.493284E+01  1.755556E+01    4.601057E+02  1.380317E+02  3.821717E+01  1.924838E+01  1.298280E-13
    111001002    -1.284408E+01 -9.296992E-01 -6.540499E+00   -3.195809E+01  9.447659E+01  2.340929E+00 -4.620755E-01  1.074264E+01
    111001003    -1.336704E+02 -4.010114E+01  1.236144E+01    4.553836E+02  1.366151E+02 -2.018654E+01  9.028288E+01 -1.651896E-13
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealForceObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

    def get_headers(self) -> List[str]:
        return ['mfx', 'mfy', 'mfxy', 'bmx', 'bmy', 'bmxy', 'syz', 'szx']

    def build(self):
        """sizes the vectorized attributes of the RealForceVU2DArray"""
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
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype=idtype)

        #[mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype=fdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        return
        #import pandas as pd
        #headers = self.get_headers()
        #assert 0 not in self.element_node
        #if self.nonlinear_factor not in (None, np.nan):
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values,
                                       #major_axis=self.element, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #self.data_frame = pd.Panel(self.data,
                                       #major_axis=self.element, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
        #self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1

        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx) = t1
                    (mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx,
                            mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, parent, coord, icord, theta,
                  vugrid, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #print('adding %s, %s' % (eid, vugrid))
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, vugrid]
        self.data[self.itime, self.itotal, :] = [mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx]
        self.itotal += 1

    #def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        #raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        #self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

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
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        if 'VUTRIA' in self.element_name:
            msg = [
                '                    F O R C E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )'
                '                 TIME =   2.500000E-03,  P-ELEMENT ID =      11,  OUTPUT COORD. ID =       0,  P OF EDGES =  3  3  3'
                '                       LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000'
                ''
                '     VUGRID                - MEMBRANE  FORCES -      - BENDING  MOMENTS -             - TRANSVERSE SHEAR FORCES -                    '
                '       ID.            FX            FY            FXY MX            MY            MXY             QX            QY                   '
                #'    111001001     1.497761E+02  4.493284E+01  1.755556E+01    4.601057E+02  1.380317E+02  3.821717E+01  1.924838E+01  1.298280E-13'
                #'    111001002    -1.284408E+01 -9.296992E-01 -6.540499E+00   -3.195809E+01  9.447659E+01  2.340929E+00 -4.620755E-01  1.074264E+01'
                #'    111001003    -1.336704E+02 -4.010114E+01  1.236144E+01    4.553836E+02  1.366151E+02 -2.018654E+01  9.028288E+01 -1.651896E-13'
            ]
            nnodes = 3
        elif 'VUQUAD' in self.element_name:
            msg = [
                '                 F O R C E S   I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'
                '                 TIME =   2.500000E-03,  P-ELEMENT ID =      11,  OUTPUT COORD. ID =       0,  P OF EDGES =  3  3  3'
                '                       LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000'
                ''
                '     VUGRID                - MEMBRANE  FORCES -      - BENDING  MOMENTS -             - TRANSVERSE SHEAR FORCES -                    '
                '       ID.            FX            FY            FXY MX            MY            MXY             QX            QY                   '
                #'    111001001     1.497761E+02  4.493284E+01  1.755556E+01    4.601057E+02  1.380317E+02  3.821717E+01  1.924838E+01  1.298280E-13'
                #'    111001002    -1.284408E+01 -9.296992E-01 -6.540499E+00   -3.195809E+01  9.447659E+01  2.340929E+00 -4.620755E-01  1.074264E+01'
                #'    111001003    -1.336704E+02 -4.010114E+01  1.236144E+01    4.553836E+02  1.366151E+02 -2.018654E+01  9.028288E+01 -1.651896E-13'
            ]
            #msg = [
                #'                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
                #' \n'
                #'    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                #'      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            #]
            nnodes = 4
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, nnodes, msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        """
                        F O R C E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
                     TIME =   2.500000E-03,  P-ELEMENT ID =      11,  OUTPUT COORD. ID =       0,  P OF EDGES =  3  3  3
                           LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000

         VUGRID                - MEMBRANE  FORCES -      - BENDING  MOMENTS -             - TRANSVERSE SHEAR FORCES -
           ID.            FX            FY            FXY MX            MY            MXY             QX            QY
        111001001     1.497761E+02  4.493284E+01  1.755556E+01    4.601057E+02  1.380317E+02  3.821717E+01  1.924838E+01  1.298280E-13
        111001002    -1.284408E+01 -9.296992E-01 -6.540499E+00   -3.195809E+01  9.447659E+01  2.340929E+00 -4.620755E-01  1.074264E+01
        111001003    -1.336704E+02 -4.010114E+01  1.236144E+01    4.553836E+02  1.366151E+02 -2.018654E+01  9.028288E+01 -1.651896E-13
        """
        if header is None:
            header = []
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        #eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        #cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            ##print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            ##[mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx]
            mfx = self.data[itime, :, 0]
            mfy = self.data[itime, :, 1]
            mfxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            syz = self.data[itime, :, 6]
            szx = self.data[itime, :, 7]

            if self.element_type in [189, 190]: # VUQUAD, VUTRIA
                # TODO: format the data properly
                for nid, mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi in zip(
                    nids, mfx, mfy, mfxy, bmx, bmy, bmxy, syz, szx):
                    [mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi] = write_floats_13e(
                        [mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi])
                    f06_file.write('   %8i  %13s %13s %13s   %13s %13s %13s  %13s %s\n' % (
                        nid, mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi))

            #elif self.element_type == 33:
                #for nid, mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    #[mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi] = write_floats_13e(
                        #[mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi])
                    ##Fmt = '% 8i   ' + '%27.20E   ' * 8 + '\n'
                    ##f06_file.write(Fmt % (eid, mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi))
                    ##
                    #f06_file.write('0 %8i %8s %13s %13s %13s %13s %13s %13s %13s %s\n' % (
                        #eid, mfxi, mfyi, mfxyi, bmxi, bmyi, bmxyi, syzi, szxi))
            else:
                raise NotImplementedError(self.element_type)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


def oef_data_code(table_name, analysis_code,
                  is_sort1=True, is_random=False,
                  random_code=0, title='', subtitle='', label='', is_msc=True):
    sort1_sort_bit = 0 if is_sort1 else 1
    random_sort_bit = 1 if is_random else 0
    sort_method = 1 if is_sort1 else 2
    assert analysis_code != 0, analysis_code
    #if format_code == 1:
        #format_word = "Real"
    #elif format_code == 2:
        #format_word = "Real/Imaginary"
    #elif format_code == 3:
        #format_word = "Magnitude/Phase"
    #DEVICE_CODE_MAP = {
        #1 : "Print",
        #2 : "Plot",
        #3 : "Print and Plot",
        #4 : "Punch",
        #5 : "Print and Punch",
        #6 : "Plot and Punch",
        #7 : "Print, Plot, and Punch",
    #}

    table_code = TABLE_NAME_TO_TABLE_CODE[table_name]
    sort_code = 1 # TODO: what should this be???

    #table_code = tCode % 1000
    #sort_code = tCode // 1000
    tCode = table_code * 1000 + sort_code

    device_code = 2  # Plot
    approach_code = analysis_code * 10 + device_code
    #print(f'approach_code={approach_code} analysis_code={analysis_code} device_code={device_code}')
    data_code = {
        'nonlinear_factor': None,
        'approach_code' : approach_code,
        'analysis_code' : analysis_code,
        'sort_bits': [0, sort1_sort_bit, random_sort_bit], # real, sort1, random
        'sort_method' : sort_method,
        'is_msc': is_msc,
        #'is_nasa95': is_nasa95,
        'format_code': 1, # real
        'table_code': table_code,
        'tCode': tCode,
        'table_name': table_name, ## TODO: should this be a string?
        'device_code' : device_code,
        'random_code' : random_code,
        'thermal': 0,
        'title' : title,
        'subtitle': subtitle,
        'label': label,
        'num_wide' : 8, # displacement-style table
    }
    return data_code
