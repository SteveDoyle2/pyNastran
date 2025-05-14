from typing import Any
import numpy as np
from pyNastran.op2.result_objects.op2_objects import BaseElement, set_as_sort1
from pyNastran.op2.op2_interface.write_utils import set_table3_field
from pyNastran.op2.writer.utils import fix_table3_types

TABLE_NAME_TO_TABLE_CODE = {
    # stress
    'OES1': 5,
    'OES1C': 5,
    'OES1X': 5,
    'OES1X1': 5,
    'OES2': 5,

    # strain
    'OSTR1': 5,
    'OSTR1X': 5,
    'OSTR2': 5,

    'OSTR1C': 5,
}

class OES_Object(BaseElement):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        self.element_type = None
        self.element_name = None
        self.nonlinear_factor = np.nan
        self._times = None
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)
        #self.log.debug("starting OES...element_name=%-6s isubcase=%s" % (self.element_name, self.isubcase))
        #print self.data_code

    def finalize(self):
        """it's required that the object be in SORT1"""
        self.set_as_sort1()

    def _get_sort2_itime_ielement_from_itotal(self) -> tuple[int, int]:
        #print(f'itime={self.itime} ielement={self.ielement} nelements={self.nelements}') # self.itotal
        #itime = self.itotal % self.nelements
        #ielement = self.itime
        #print(f'itime={itime} ielement={ielement} nelements={self.nelements}') # self.itotal
        #itime = self.itotal % self.nelements
        #return itime, ielement
        # ------------------
        itime = self.ielement
        #itotal = self.itotal
        ielement = self.itime
        #print(f'itime={self.itime} ielement={self.ielement} nelements={self.nelements}') # self.itotal
        return itime, ielement

    def set_as_sort1(self):
        """the data is in SORT1, but the flags are wrong"""
        set_as_sort1(self)
        #update_stress_force_time_word(self)

    def _get_analysis_code_dtype(self) -> str:
        if self.analysis_code in [5, 6]:
            if self.size == 4:
                dtype = 'float32'
            else:
                dtype = 'float64'
        else:
            raise NotImplementedError(self.data_code)
        return dtype

    @property
    def is_curvature(self) -> bool:
        if self.is_stress:
            curvature_flag = False
        else:
            # strain only
            curvature_flag = True if self.stress_bits[2] == 0 else False
        if self.s_code in [10, 11, 20, 27]:
            assert curvature_flag, curvature_flag
            return True
        assert not curvature_flag, curvature_flag
        return False

    @property
    def is_fiber_distance(self) -> bool:
        return not self.is_curvature

    @property
    def is_von_mises(self) -> bool:
        return not self.is_max_shear

    @property
    def is_max_shear(self):
        return self.stress_bits[4] == 0

    def _get_headers(self):
        raise NotImplementedError(f'overwrite this {self.class_name}')

    @property
    def is_stress(self):
        raise NotImplementedError(f'overwrite this in {self.class_name}')

    def _write_table_3(self, op2_file, op2_ascii, new_result, itable, itime): #, itable=-3, itime=0):
        import inspect
        from struct import pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        #if itable == -3:
        #print('*writing itable=%s' % itable)
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
        op2_file.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('table_3_header = %s\n' % header)
        #op2_file.write(pack('12i', *header))
        #else:
            #print('***writing itable=%s' % itable)
            #op2_file.write(pack('3i', *[
                ##4, itable, 4,
                ##4, 1, 4,
                ##4, 0, 4,
                #4, 146, 4,
            #]))
        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        element_type = self.element_type
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
        s_code = self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = self.thermal
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        unused_oCode = 0

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0
        if self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
            if np.isnan(field5):  # poor sort2 -> sort1
                raise RuntimeError('field5 in a static case is nan...; do you have SORT2?')
                #field5 = 1

        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            field6 = self.eigns[itime]
            field7 = self.cycles[itime]
            assert isinstance(field6, float), type(field6)
            assert isinstance(field7, float), type(field7)
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7

        #elif self.analysis_code == 3:
            #field5 = self.freqs[itime]
        elif self.analysis_code == 5:
            field5 = self.freqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            field5 = self.dts[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 7:  # pre-buckling
            field5 = self.lsdvmns[itime] # load set number
        elif self.analysis_code == 8:  # post-buckling
            field5 = self.lsdvmns[itime] # load set number
            #if hasattr(self, 'eigns'):
            if hasattr(self, 'eigens'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        elif self.analysis_code == 9:  # complex eigenvalues
            field5 = self.modes[itime]
            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=9')

            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            field7 = self.eigis[itime]
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.lftsfqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            field5 = self.lsdvmns[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, element_type, isubcase, field5,
            field6, field7, self.load_set, format_code, num_wide,
            s_code, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        table3 = fix_table3_types(table3, size=4)
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2_file.write(pack(fmt, *data))


class StressObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)
        assert self.is_stress, self.code_information()
        assert not self.is_strain, self.code_information()

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        #assert dt >= 0.
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating stress...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    @property
    def is_strain(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        if self.table_name_str != 'OES1A':
            assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        return False

    @property
    def is_stress(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        if self.table_name_str != 'OES1A':
            assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        return True


class StrainObject(OES_Object):
    def __init__(self, data_code, isubcase):
        OES_Object.__init__(self, data_code, isubcase)
        assert self.is_strain, self.code_information()
        assert not self.is_stress, self.code_information()

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        self.element_name = self.data_code['element_name']
        if dt is not None:
            #print("updating strain...%s=%s element_name=%s" %
            #     (self.data_code['name'], dt, self.element_name))
            self.dt = dt
            self.add_new_transient(dt)

    @property
    def is_strain(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == 1, 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        return True

    @property
    def is_stress(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == 1, 'class_name=%s is_stress=False scode=%s stress_bits=%s; element_type=%s element_name=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.element_type, self.element_name, self.table_name)
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        return False

def get_scode(stress_bits: list[int]) -> int:
    scode = 0
    for i, bit in enumerate(stress_bits):
        scode += 2 ** i * bit
    return scode


def oes_complex_data_code(table_name: str,
                          element_name: str, num_wide: int,
                          is_sort1: bool=True, is_random: bool=False,
                          random_code=0, title='', subtitle='', label='', is_msc=True):
    dtype_code = 1 # complex
    data_code = _oes_data_code(table_name,
                               element_name, num_wide, dtype_code,
                               is_sort1=is_sort1,
                               is_random=is_random, random_code=random_code,
                               title=title, subtitle=subtitle, label=label, is_msc=is_msc)
    return data_code

def oes_real_data_code(table_name: str,
                       element_name: str, num_wide: int,
                       is_sort1: bool=True, is_random: bool=False,
                       random_code=0, title='', subtitle='', label='', is_msc=True):
    dtype_code = 0 # real
    assert isinstance(element_name, str), element_name
    data_code = _oes_data_code(table_name, element_name, num_wide, dtype_code,
                               is_sort1=is_sort1,
                               is_random=is_random, random_code=random_code,
                               title=title, subtitle=subtitle, label=label, is_msc=is_msc)
    return data_code

def _oes_data_code(table_name: str,
                   element_name: str, num_wide: int, dtype_code: int,
                   is_sort1: bool=True, is_random: bool=False,
                   random_code=0, title='', subtitle='', label='', is_msc=True):
    """
    Parameters
    ----------
    dtype_code : int
      0 : real
      1 : complex
      2 : random
    """
    sort1_sort_bit = 0 if is_sort1 else 1
    random_sort_bit = 1 if is_random else 0
    sort_method = 1 if is_sort1 else 2
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
    #approach_code = analysis_code * 10 + device_code
    #print(f'approach_code={approach_code} analysis_code={analysis_code} device_code={device_code}')
    data_code = {
        'nonlinear_factor': None,
        #'approach_code' : approach_code,
        #'analysis_code' : analysis_code,
        'sort_bits': [dtype_code, sort1_sort_bit, random_sort_bit], # real, sort1, random
        'sort_method' : sort_method,
        'is_msc': is_msc,
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
        'num_wide': num_wide,
        'element_name': element_name,
        #'num_wide' : 8, # displacement-style table
    }
    return data_code

def set_approach_code(analysis_code: int, device_code: int):
    approach_code = analysis_code * 10 + device_code
    return approach_code

def _check_num_wide(data_code: dict[str, Any]):
    if 'num_wide' not in data_code:
        if 'element_name' in data_code:
            element_name = data_code['element_name']
            raise RuntimeError(f'missing num_wide for {element_name}')
        else:
            raise RuntimeError(f'missing num_wide')

def set_static_case(cls, is_sort1, isubcase,
                    data_code, func, args):
    data_code['lsdvmns'] = [0] # TODO: ???
    data_code['data_names'] = []
    data_code['load_set'] = 1
    data_code['analysis_code'] = 1
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])
    times = [None]
    _check_num_wide(data_code)
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, times)
    obj.is_built = True
    obj.get_stats()
    return obj

def set_modal_case(cls, is_sort1, isubcase, data_code,
                   func, args, modes, eigns, cycles):

    data_code['data_names'] = ['modes', 'eigns', 'mode_cycles']
    data_code['load_set'] = 1
    #data_code['lsdvmns'] = [0] # TODO: ???
    data_code['analysis_code'] = 2 # modal
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])

    _check_num_wide(data_code)
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, modes)
    obj.modes = modes
    obj.eigns = eigns

    obj.mode_cycles = cycles
    obj.cycles = cycles
    obj.is_built = True
    obj.get_stats()
    return obj

def set_transient_case(cls, is_sort1: bool, isubcase: int,
                       data_code, func, args, times: np.ndarray):
    data_code['lsdvmns'] = [0] # TODO: ???
    data_code['data_names'] = ['dt']
    data_code['load_set'] = 1
    data_code['analysis_code'] = 6
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])

    _check_num_wide(data_code)
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, times)
    obj.times = times
    obj.dts = times
    obj.is_built = True
    obj.get_stats()
    return obj

def set_post_buckling_case(cls, is_sort1, isubcase,
                           data_code, func, args,
                           modes, eigrs, eigis):
    data_code['lsdvmns'] = [0] # TODO: ???
    data_code['data_names'] = ['modes', 'eigrs', 'eigis']
    data_code['load_set'] = 1
    data_code['analysis_code'] = 8
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])

    _check_num_wide(data_code)
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, modes)
    obj.dts = modes
    obj.modes = modes
    obj.eigrs = eigrs
    obj.eigis = eigis
    obj.is_built = True
    obj.get_stats()
    return obj

def set_freq_case(cls, is_sort1, isubcase,
                  data_code, func, args, freqs):
    #data_code['lsdvmns'] = [0] # TODO: ???
    #data_code['data_names'] = ['dt']
    #data_code['load_set'] = 1
    data_code['load_set'] = 1
    data_code['analysis_code'] = 5
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])
    _check_num_wide(data_code)

    data_code['data_names'] = ['freq']
    data_code['name'] = 'FREQ'
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, freqs)
    obj.freqs = freqs
    obj.is_built = True
    obj.get_stats()
    return obj

def set_complex_modes_case(cls, is_sort1, isubcase,
                           data_code, func, args, modes, eigrs, eigis):
    #data_code['lsdvmns'] = [0] # TODO: ???
    #data_code['data_names'] = ['dt']
    #data_code['load_set'] = 1

    data_code['load_set'] = 1
    data_code['analysis_code'] = 9  # complex eigenvalues
    data_code['approach_code'] = set_approach_code(data_code['analysis_code'],
                                                   data_code['device_code'])
    _check_num_wide(data_code)

    data_code['data_names'] = ['modes', 'eigrs', 'eigis']
    data_code['name'] = 'mode'
    obj = func(cls, data_code, is_sort1, isubcase,
               *args, modes)
    obj.modes = modes
    obj.eigrs = eigrs
    obj.eigis = eigis
    obj.is_built = True
    obj.get_stats()
    return obj

def set_element_case(cls, data_code, is_sort1, isubcase,
                     element, data, times):
    """
    CBAR
    CROD, CONROD, CTUBE
    CELAS1-4, CDAMP1-4
    CSHEAR
    """
    assert element.ndim == 1, element.shape
    ntimes = data.shape[0]
    nnodes = data.shape[1]
    dt = times[0]
    obj = cls(data_code, is_sort1, isubcase, dt)
    obj.element = element
    obj.data = data

    obj.ntimes = ntimes
    obj.ntotal = nnodes
    obj.nelements = nnodes
    obj._times = times
    return obj

def set_element_node_xxb_case(cls, data_code, is_sort1, isubcase,
                              element_node, xxb, data, times):
    """CBEAM"""
    assert element_node.ndim == 2, element_node.shape
    ntimes = data.shape[0]
    nelement_nnodes = data.shape[1]
    dt = times[0]
    obj = cls(data_code, is_sort1, isubcase, dt)
    obj.nnodes = 11
    obj.element = element_node[:, 0]
    obj.xxb = xxb
    obj.element_node = element_node
    obj.data = data

    assert element_node.shape == (nelement_nnodes, 2)

    obj.ntimes = ntimes
    obj.ntotal = nelement_nnodes
    obj.nelements = nelement_nnodes
    obj._times = times
    return obj

def update_stress_force_time_word(obj) -> None:
    if obj.analysis_code == 1:
        name = 'lsdvmn'
        obj._times = np.full(1, np.nan, dtype=obj.data.dtype)
        obj.nonlinear_factor = obj._times[0]
    else:
        try:
            name = obj.analysis_method
        except AttributeError:
            msg = str(obj.get_stats())
            msg += obj.object_stats()
            raise RuntimeError(msg) from e

        #print(f'\n{self.class_name}: name = {name}')
        #_analysis_code_fmt
    #if self.analysis_code == 5:
        #name = 'freq'
    #elif self.analysis_code == 6:
        #name = 'dt'
    #elif self.analysis_code == 10:
        #name = 'dt'
    #else:
        #raise NotImplementedError(self.data_code)

    names = name + 's'
    obj.data_code['name'] = name
    obj.data_code['data_names'][0] = name
    #print(self.element_ids)
    del obj.element_id, obj.element_ids
    #print('timesA =', self._times)
    old_times = None
    if 1:
        old_times = obj._times
        setattr(obj, names, obj._times)
        #print('timesB =', self._times)

    if obj.analysis_code == 1:
        assert len(obj._times) == 1, obj.data.shape
    elif len(obj._times) > 1:
        class_name = obj.__class__.__name__
        #assert obj._times.min() != obj._times.max(), f'{class_name}: old_times={old_times} -> times={obj._times}; data.shape={obj.data.shape}\n{obj.code_information()}'
    #print(obj.object_attributes())

