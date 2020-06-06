import numpy as np
from pyNastran.op2.result_objects.op2_objects import BaseElement
from pyNastran.op2.op2_interface.write_utils import set_table3_field

SORT2_TABLE_NAME_MAP = {
    'OES2' : 'OES1',
    'OES2C' : 'OES1C',
    'OESATO2' : 'OESATO1',
    'OESCRM2' : 'OESCRM1',
    'OESNO2' : 'OESNO1',
    'OESPSD2' : 'OESPSD1',
    'OESRMS2' : 'OESRMS1',
    'OESNLXR2' : 'OESNLXR',

    'OSTR2' : 'OSTR1',
    'OSTR2C' : 'OSTR1C',
    'OSTRATO2' : 'OSTRATO1',
    'OSTRCRM2' : 'OSTRCRM1',
    'OSTRNO2' : 'OSTRNO1',
    'OSTRPSD2' : 'OSTRPSD1',
    'OSTRRMS2' : 'OSTRRMS1',
    'OESVM2' : 'OESVM1',
    'OSTRVM2' : 'OSTRVM1',
    'OESNL2' : 'OESNL1',
    'OSTPSD2C' : 'OSTPSD1C',
    'OESPSD2C' : 'OESPSD1C',
}

TABLE_NAME_TO_TABLE_CODE = {
    # stress
    'OES1': 5,
    'OES1X1': 5,
    'OES2': 5,

    # strain
    'OSTR1': 5,
    'OSTR2': 5,
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

    def set_as_sort1(self):
        """the data is in SORT1, but the flags are wrong"""
        if self.is_sort1:
            return
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_bits[1] = 0 # sort1
        self.sort_method = 1
        assert self.is_sort1 is True, self.is_sort1

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
        return True if self.stress_bits[4] == 0 else False

    def _get_headers(self):
        raise NotImplementedError('overwrite this')

    @property
    def is_stress(self):
        raise NotImplementedError('overwrite this')

    def _write_table_3(self, op2, op2_ascii, new_result, itable, itime): #, itable=-3, itime=0):
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
        op2.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('table_3_header = %s\n' % header)
        #op2.write(pack('12i', *header))
        #else:
            #print('***writing itable=%s' % itable)
            #op2.write(pack('3i', *[
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
        thermal = 0
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

        n = 0
        for v in table3:
            if isinstance(v, (int, float)):
                n += 4
            elif isinstance(v, str):
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
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        return False

    @property
    def is_stress(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
        assert self.stress_bits[1] == 0, 'class_name=%s scode=%s stress_bits=%s' % (class_name, self.s_code, self.stress_bits)
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
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        assert self.stress_bits[1] == 1, 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        return True

    @property
    def is_stress(self) -> bool:
        class_name = self.__class__.__name__
        assert self.stress_bits[1] == self.stress_bits[3], 'class_name=%s scode=%s stress_bits=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.table_name)
        assert self.stress_bits[1] == 1, 'class_name=%s is_stress=False scode=%s stress_bits=%s; element_type=%s element_name=%s; table_name=%r' % (class_name, self.s_code, self.stress_bits, self.element_type, self.element_name, self.table_name)
        return False


def oes_data_code(table_name, analysis_code,
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
