#pylint: disable=C0301,C0111
from __future__ import annotations
import copy
import warnings
from itertools import count
from struct import pack
from typing import TYPE_CHECKING
import numpy as np

from cpylog import SimpleLogger
from pyNastran import stop_on_op2_missed_table
from pyNastran.utils import object_attributes, object_methods, object_stats, simplify_object_keys
from pyNastran.utils.numpy_utils import integer_types

from pyNastran.op2.errors import OverwriteTableError
from pyNastran.op2.op2_interface.op2_codes import Op2Codes, get_sort_method_from_table_name
from pyNastran.op2.op2_interface.write_utils import write_table_header, export_to_hdf5
if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd
Date = tuple[int, int, int]

NULL_GRIDTYPE = {538976288, 1065353216}
GRID_TYPE_INT_TO_STR = {
    1: 'G',  # GRID
    2: 'S',  # SPOINT
    3: 'E',  # EXTRA POINT
    4: 'M',  # MODAL POINT
    7: 'L',  # RIGID POINT (e.g. RBE3)
    0: 'H',  # SECTOR/HARMONIC/RING POINT
    -1: '',
}
GRID_TYPE_TO_STR_MAP = {
    'G': 1,  # GRID
    'S': 2,  # SPOINT
    'L': 7,  # RIGID POINT (e.g. RBE3)
    'H': 0,  # SECTOR/HARMONIC/RING POINT
}

SORT2_TABLE_NAME_MAP = {
    # sort2_name : sort1_name
    # displacement
    'OUGATO2': 'OUGATO1',
    'OUGCRM2': 'OUGCRM1',
    'OUGNO2': 'OUGNO1',
    'OUGPSD2': 'OUGPSD1',
    'OUGRMS2': 'OUGRMS1',

    # velocity
    'OVGATO2': 'OVGATO1',
    'OVGCRM2': 'OVGCRM1',
    'OVGNO2': 'OVGNO1',
    'OVGPSD2': 'OVGPSD1',
    'OVGRMS2': 'OVGRMS1',

    # acceleration
    'OAGATO2': 'OAGATO1',
    'OAGCRM2': 'OAGCRM1',
    'OAGNO2': 'OAGNO1',
    'OAGPSD2': 'OAGPSD1',
    'OAGRMS2': 'OAGRMS1',

    # spc forces
    'OQGATO2': 'OQGATO1',
    'OQGCRM2': 'OQGCRM1',
    'OQGNO2': 'OQGNO1',
    'OQGPSD2': 'OQGPSD1',
    'OQGRMS2': 'OQGRMS1',

    # mpc forces
    'OQMATO2': 'OQMATO1',
    'OQMCRM2': 'OQMCRM1',
    'OQMNO2': 'OQMNO1',
    'OQMPSD2': 'OQMPSD1',
    'OQMRMS2': 'OQMRMS1',

    # load vectors
    'OPGATO2': 'OPGATO1',
    'OPGCRM2': 'OPGCRM1',
    'OPGNO2': 'OPGNO1',
    'OPGPSD2': 'OPGPSD1',
    'OPGRMS2': 'OPGRMS1',

    # pressure
    'OPRATO2': 'OPRATO1',
    'OPRCRM2': 'OPRCRM1',
    'OPRNO2': 'OPRNO1',
    'OPRPSD2': 'OPRPSD1',
    'OPRRMS2': 'OPRRMS1',

    #'OUG2' : 'OUG1',
    'OUGV2': 'OUGV1',
    'OQG2': 'OQG1',
    'OQMG2': 'OQMG1',
    'OPG2': 'OPG1',
    'OPNL2': 'OPNL1',
    'OUXY2': 'OUXY1',
    'OQGGF2': 'OQGGF1',
    'OQGCF2': 'OQGCF1',
    'OUGF2': 'OUGF1',
    # --------------------
    # OES
    'OES2': 'OES1',
    'OES2C': 'OES1C',
    'OESATO2': 'OESATO1',
    'OESCRM2': 'OESCRM1',
    'OESNO2': 'OESNO1',
    'OESPSD2': 'OESPSD1',
    'OESRMS2': 'OESRMS1',
    'OESNLXR2': 'OESNLXR',
    'OESVM2': 'OESVM1',
    'OESNL2': 'OESNL1',
    'OESPSD2C': 'OESPSD1C',

    # OSTR
    'OSTR2': 'OSTR1',
    'OSTR2C': 'OSTR1C',
    'OSTRATO2': 'OSTRATO1',
    'OSTRCRM2': 'OSTRCRM1',
    'OSTRNO2': 'OSTRNO1',
    'OSTRPSD2': 'OSTRPSD1',
    'OSTRRMS2': 'OSTRRMS1',
    'OSTRVM2': 'OSTRVM1',
    'OSTPSD2C': 'OSTPSD1C',

    # OEF
    'OEF2': 'OEF1',
    'OEFATO2': 'OEFATO1',
    'OEFCRM2': 'OEFCRM1',
    'OEFPSD2': 'OEFPSD1',
    'OEFRMS2': 'OEFRMS1',
    'OEFNO2': 'OEFNO1',

    # ONR / OEE
    'ONRGY2': 'ONRGY1',
}

SORT1_TABLES = list(SORT2_TABLE_NAME_MAP.values())
SORT1_TABLES.extend([
    'BOUGV1',
    'OUG1F',
    'BOUGF1',
    'OUG1',
    'OVG1',
    'OAG1',
])
SORT2_TABLES = list(SORT2_TABLE_NAME_MAP.keys())


class BaseScalarObject(Op2Codes):
    """
    The base scalar class is used by:
     - RealEigenvalues
     - BucklingEigenvalues
     - ComplexEigenvalues
     - ScalarObject

    """
    def __init__(self):
        Op2Codes.__init__(self)
        self.is_built = False

        # init value
        #  None - static
        #  int/float - transient/freq/load step/eigenvalue
        self.nonlinear_factor = np.nan
        self.format_code = None
        self.sort_code = None
        self.table_code = None
        self.title = None
        self.subtitle = None
        self.label = None
        self.num_wide = None
        self.device_code = None
        self.table_name = None
        #self.ntimes = 0
        #self.ntotal = 0
        #assert isinstance(self.name, (str, bytes)), 'name=%s type=%s' % (self.name, type(self.name))

    def object_attributes(self, mode: str='public', keys_to_skip=None,
                          filter_properties: bool=False) -> list[str]:
        keys_to_skip = simplify_object_keys(keys_to_skip)

        my_keys_to_skip = [
            'object_methods', 'object_attributes', ', object_stats',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public', keys_to_skip=None) -> list[str]:
        keys_to_skip = simplify_object_keys(keys_to_skip)

        my_keys_to_skip = [
            'object_methods', 'object_attributes', ', object_stats',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_stats(self, mode: str='public', keys_to_skip=None,
                     filter_properties: bool=False) -> str:
        keys_to_skip = simplify_object_keys(keys_to_skip)

        my_keys_to_skip = [
            'object_methods', 'object_attributes', ', object_stats',
        ]
        return object_stats(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                            filter_properties=filter_properties)

    def __eq__(self, table) -> bool:  # pragma: no cover
        #raise NotImplementedError(str(self.get_stats()))
        return False

    def __ne__(self, table) -> bool:
        return not self == table

    @property
    def class_name(self) -> str:
        return self.__class__.__name__

    def get_headers(self):  # pragma: no cover
        raise RuntimeError()

    def _get_stats_short(self):  # pragma: no cover
        raise NotImplementedError('_get_stats_short')

    def build_dataframe(self) -> None:  # pragma: no cover
        """creates a pandas dataframe"""
        print('build_dataframe is not implemented in %s' % self.__class__.__name__)

    def export_to_hdf5(self, group, log: SimpleLogger) -> None:
        """exports the object to HDF5 format"""
        export_to_hdf5(self, group, log)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True) -> int:
        if header is None:
            header = []
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient(header, page_stamp, page_num=page_num, f06_file=f06_file,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = 'write_f06 is not implemented in %s\n' % self.__class__.__name__
        f06_file.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f06_file=None,
                             is_mag_phase=False, is_sort1=True) -> int:
        msg = '_write_f06_transient is not implemented in %s\n' % self.__class__.__name__
        f06_file.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def __repr__(self) -> str:
        return ''.join(self.get_stats())

    def get_stats(self, short: bool=False):
        msg = 'get_stats is not implemented in %s\n' % self.__class__.__name__
        if stop_on_op2_missed_table:
            raise NotImplementedError(msg)
        return msg


class ScalarObject(BaseScalarObject):
    """
    Used by all vectorized objects including:
     - DisplacementArray
     - RodStressArray
    """
    def __init__(self, data_code, isubcase, apply_data_code=True):
        assert 'nonlinear_factor' in data_code, data_code

        # new terms...I think they're all valid...
        self.result_name = None
        self.approach_code = None
        self.analysis_code = None
        self.data = None
        self._times = None

        self.isubcase = None
        self.ogs = None
        self.pval_step = None
        self.name = None
        self.superelement_adaptivity_index = None
        self._count = None
        #--------------------------------

        BaseScalarObject.__init__(self)
        self.isubcase = isubcase

        #--------------------------------
        self.data_frame = None

        # word size:
        #  32-bit: 4
        #  64-bit: 8
        self.size = 4

        # the nonlinear factor; None=static; float=transient
        self.dt = None
        # the number of time steps
        self.ntimes = 0
        # the number of entries in a table for a single time step
        self.ntotal = 0
        # there are a few tables (e.g. GridPointForcesArray) that change
        # length from time=1 to time=2.  As such, we will track the
        # length of each time step
        self._ntotals = []

        self.load_as_h5 = False
        self.h5_file = None
        if 'load_as_h5' in data_code:
            self.load_as_h5 = data_code['load_as_h5']
            del data_code['load_as_h5']
        if 'h5_file' in data_code:
            self.h5_file = data_code['h5_file']
            del data_code['h5_file']

        self.data_code = copy.deepcopy(data_code)
        #self.table_name = self.data_code['table_name']
        # if data code isn't being applied and you don't have
        # parameters that were in data_code (e.g.
        # self.element_name/self.element_type), you need to define
        # your vectorized class as setting data_code
        #
        # see pyNastran.op2.oes_stressStrain.real.oes_bars for an example
        # it's really subtle...
        #
        # since you won't get it, the inherited RealBarStressArray class inits
        # data code, but RealBarArray does not.
        if apply_data_code:
            self.apply_data_code()
            self._set_data_members()
        #print(self.code_information())

    def _get_stats_short(self) -> list[str]:
        msg = []
        class_name = self.__class__.__name__
        if hasattr(self, 'data'):
            unused_data = self.data
            shape = [int(i) for i in self.data.shape]
            headers = self.get_headers()
            headers_str = str(', '.join(headers))
            msg.append(f'{class_name}[{self.isubcase}]; {shape}; [{headers_str}]\n')
        return msg

    def __eq__(self, table) -> bool:  # pragma: no cover
        self._eq_header(table)
        #raise NotImplementedError(self.class_name)
        #raise NotImplementedError(str(self.get_stats()))
        return False

    def _eq_header(self, table) -> None:
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor, 'nonlinear_factor=%s table.nonlinear_factor=%s' % (self.nonlinear_factor, table.nonlinear_factor)
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (
            self.table_name, table.table_name)
        assert self.approach_code == table.approach_code

    def __getstate__(self):
        """we need to remove the saved functions"""
        state = self.__dict__.copy()
        if 'add' in state:
            del state['add']
        if 'add_new_eid' in state:
            del state['add_new_eid']
        if 'class_name' in state:
            del state['class_name']
        if 'nnodes_per_element' in state:
            del state['nnodes_per_element']
        if 'add_node' in state:
            del state['add_node']
        if 'add_eid' in state:
            del state['add_eid']
        if '_add' in state:
            del state['_add']
        #if '_add_new_eid_sort1' in state:
            #del state['_add_new_eid_sort1']
        #if '_add_new_node_sort1' in state:
            #del state['_add_new_node_sort1']
        if '_add_new_eid' in state:
            del state['_add_new_eid']
        if '_add_new_node' in state:
            del state['_add_new_node']
        if 'dataframe' in state:
            del state['dataframe']

        # for key, value in state.items():
        #     if isinstance(value, (int, float, str, np.ndarray, list)) or value is None:
        #         continue
        #     print(' %s = %s' % (key, value))
        # print(state)
        return state

    def _get_result_group(self):
        """gets the h5 result group"""
        code = self._get_code()
        case_name = 'Subcase=%s' % str(code)  # (self.isubcase)
        if case_name in self.h5_file:
            subcase_group = self.h5_file[case_name]
        else:
            subcase_group = self.h5_file.create_group(case_name)
        group = subcase_group.create_group(self.result_name)
        return group

    def _get_code(self) -> tuple[int, int, int, int, int, str, str]:
        code = self.isubcase
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        # if self.binary_debug:
        #     self.binary_debug.write(self.code_information(include_time=True))

        code = (self.isubcase, self.analysis_code, self._sort_method(), self._count, ogs,
                self.superelement_adaptivity_index, self.pval_step)
        # code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
        #         self.superelement_adaptivity_index, self.table_name_str)
        # print('%r' % self.subtitle)
        # self.code = code
        # self.log.debug('code = %s' % str(self.code))
        return code

    #@property
    def _sort_method(self) -> int:
        try:
            sort_method, unused_is_real, unused_is_random = self._table_specs()
        except Exception:
            sort_method = get_sort_method_from_table_name(self.table_name)
        #is_sort1 = self.table_name.endswith('1')
        #is_sort1 = self.is_sort1  # uses the sort_bits
        assert sort_method in [1, 2], 'sort_method=%r\n%s' % (sort_method, self.code_information())
        return sort_method

    @property
    def dataframe(self) -> pd.DataFrame:
        """alternate way to get the dataframe"""
        return self.data_frame

    def apply_data_code(self) -> None:
        #print(self.__class__.__name__)
        if self.table_name is not None and self.table_name != self.data_code['table_name']:
            #print(self.data_code)
            msg = 'old_table_name=%r new_table_name=%r' % (
                self.table_name, self.data_code['table_name'])
            raise OverwriteTableError(msg)
        for key, value in sorted(self.data_code.items()):
            if isinstance(value, bytes):
                print("  key=%s value=%s; value is bytes" % (key, value))
            self.__setattr__(key, value)
            #print("  key=%s value=%s" % (key, value))

    def get_data_code(self, prefix: str='  ') -> list[str]:
        msg = ''
        if 'data_names' not in self.data_code:
            return ['']

        msg += '%ssort1\n' % prefix if self.is_sort1 else '%ssort2\n' % prefix
        for name in self.data_code['data_names']:
            if hasattr(self, name + 's'):
                vals = getattr(self, name + 's')
                name += 's'
                vals_array = np.array(vals)
            elif hasattr(self, name):
                vals = getattr(self, name)
                vals_array = np.array(vals)
            else:
                vals_array = '???'
            #msg.append('%s = [%s]\n' % (name, ', '.join(['%r' % val for val in vals])))
            if isinstance(vals_array, np.ndarray):
                dtypei = vals_array.dtype.name
            else:
                if isinstance(vals_array, str):
                    dtypei = 'str'
                else:
                    dtypei = type(vals_array)
            msg += f'{prefix}{name} = {vals_array}; dtype={dtypei}\n'
        #print("***data_names =", self.data_names)
        return [msg]

    def get_unsteady_value(self):
        name = self.data_code['name']
        return self._get_var(name)

    def _get_var(self, name):
        return getattr(self, name)

    def _set_var(self, name, value):
        return self.__setattr__(name, value)

    def _start_data_member(self, var_name, value_name):
        if hasattr(self, var_name):
            return True
        elif hasattr(self, value_name):
            self._set_var(var_name, [])
            return True
        return False

    def _append_data_member(self, var_name, value_name):
        """
        this appends a data member to a variable that may or may not exist
        """
        has_list = self._start_data_member(var_name, value_name)
        if has_list:
            list_a = self._get_var(var_name)
            if list_a is not None:
                #print("has %s" % var_name)
                value = self._get_var(value_name)
                try:
                    n = len(list_a)
                except Exception:
                    print("listA = ", list_a)
                    raise
                list_a.append(value)
                assert len(list_a) == n + 1

    def _set_data_members(self):
        if 'data_names' not in self.data_code:
            msg = ('No "transient" variable was set for %s ("data_names" '
                   'was not defined in self.data_code).\n' % self.table_name)
            raise NotImplementedError(msg + self.code_information())

        for name in self.data_code['data_names']:
            #print(name)
            self._append_data_member(name + 's', name)

    def update_data_code(self, data_code):
        if not self.data_code or (data_code['nonlinear_factor'] != self.data_code['nonlinear_factor']):
            self.data_code = data_code
            self.apply_data_code()  # take all the parameters in data_code and make them attributes of the class
            self._set_data_members()  # set the transient variables
        #else:
            #print('NF_new=%r NF_old=%r' % (data_code['nonlinear_factor'], self.data_code['nonlinear_factor']))

    def print_data_members(self):
        """
        Prints out the "unique" vals of the case.

        Uses a provided list of data_code['data_names'] to set the values for
        each subcase.  Then populates a list of self.name+'s' (by using
        setattr) with the current value.  For example, if the variable name is
        'mode', we make self.modes.  Then to extract the values, we build a
        list of of the variables that were set like this and then loop over
        them to print their values.

        This way there is no dependency on one result type having ['mode'] and
        another result type having ['mode','eigr','eigi'].
        """
        key_vals = []
        for name in self.data_code['data_names']:
            vals = getattr(self, name + 's')
            key_vals.append(vals)
            #print("%ss = %s" % (name, vals))

        msg = ''
        for name in self.data_code['data_names']:
            msg += '%-10s ' % name
        msg += '\n'

        nmodes = len(key_vals[0])
        for i in range(nmodes):
            for vals in key_vals:
                msg += '%-10g ' % vals[i]
            msg += '\n'
        return msg + '\n'

    def recast_gridtype_as_string(self, grid_type):
        """
        converts a grid_type integer to a string

        Point type (per NX 10; OUG table; p.5-663):
        -1, Harmonic Point (my choice) -> '    ' = 538976288 (as an integer)
        =1, GRID Point
        =2, Scalar Point
        =3, Extra Point
        =4, Modal
        =5, p-elements, 0-DOF
        -6, p-elements, number of DOF
        """
        try:
            grid_type_str = GRID_TYPE_INT_TO_STR[grid_type]
        except KeyError:
            if grid_type in NULL_GRIDTYPE:  # 32/64 bit error...
                warnings.warn(''.join(self.get_stats()))
            raise RuntimeError(f'grid_type={grid_type!r}')
        return grid_type_str

    def cast_grid_type(self, grid_type_str):
        """converts a grid_type string to an integer"""
        try:
            grid_type = GRID_TYPE_TO_STR_MAP[grid_type_str]
        except KeyError:
            raise RuntimeError(f'grid_type={grid_type_str!r}')
        return grid_type

    def update_dt(self, data_code, unused_dt):
        """
        This method is called if the object already exits and a new
        time/freq/load step is found
        """
        self.data_code = data_code
        self.apply_data_code()
        msg = 'update_dt not implemented in the %s class' % self.__class__.__name__
        raise RuntimeError(msg)
        # assert dt>=0.
        # print("updating dt...dt=%s" % dt)
        # if dt is not None:
        #     self.dt = dt
        #     self.add_new_transient()

    def _build_dataframe_transient_header(self):
        """builds the header for the Pandas DataFrame/table"""
        assert isinstance(self.name, (str, bytes)), 'name=%s type=%s' % (self.name, type(self.name))
        # print('self.name = %r' % self.name)
        # name = self.name #data_code['name']
        times = self._times
        utimes = np.unique(times)
        if not len(times) == len(utimes):
            msg = 'WARNING : %s - times=%s unique_times=%s...assuming new values...new_times=%s' % (
                self.__class__.__name__, times, utimes, np.arange(len(times)))
            print(msg)
            #raise RuntimeError(msg)
            times = np.arange(len(times))
        column_names = []
        column_values = []

        data_names = self.data_code['data_names']
        #print(f'data_names = {data_names}')
        for name in data_names:
            #if name == primary_name:
            #times = self.da
            times = np.array(getattr(self, name + 's'))
            if name == 'mode':
                column_names.append('Mode')
                column_values.append(times)

                # if freq not in data_names:
                # if name == 'freq':
                # #if hasattr(self, 'freqs'):
                #     column_names.append('Freq')
                #     column_values.append(self.freqs)
                # elif name == 'eigr':
                #     column_names.append('eigenvalue_real')
                #     column_values.append(self.eigrs)
                # elif hasattr(self, 'eigrs') and 0:
                #     try:
                #         abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
                #     except FloatingPointError:
                #         msg = 'Cant analyze freq = sqrt(eig)/(2*pi)\neigr=%s\n' % (self.eigrs)
                #         abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
                #         msg += 'freq = sqrt(abs(self.eigrs)) / (2 * np.pi)=%s' % abs_freqs
                #         raise FloatingPointError(msg)
                #     column_names.append('Freq')
                #     column_values.append(abs_freqs)
                # else:
                #     pass

                # Convert eigenvalues to frequencies
                # TODO: add damping header
            elif name in ['eign']:
                omega_radians, abs_freqs = real_modes_to_omega_freq(self.eigns)
                column_names.append('Freq')
                column_values.append(abs_freqs)
                column_names.append('Eigenvalue')
                column_values.append(times)
                column_names.append('Radians')
                column_values.append(omega_radians)

            elif name in ['eigr']:
                column_names.append('EigenvalueReal')
                column_values.append(times)

            elif name in ['eigi']:
                column_names.append('EigenvalueImag')
                column_values.append(times)
                eigr = np.array(self.eigrs)
                eigi = np.array(self.eigis)
                damping, frequncy = complex_damping_frequency(eigr, eigi)
                column_names.append('Damping')
                column_values.append(damping)
            elif name in ['mode_cycle']:
                continue
                #column_names.append('mode_cycle(Freq?)')
                #column_values.append(times)
            elif name in ['mode2']:
                continue
                #column_names.append('mode2(Freq?)')
                #column_values.append(times)
            elif name in ['cycle']:
                continue
                #column_names.append('Freq (Cycles/s)')
                #column_values.append(times)

            elif name in ['freq']:
                column_names.append('Freq')
                column_values.append(times)
            elif name in ['dt', 'time']:
                column_names.append('Time')
                column_values.append(times)
            elif name in ['lftsfq', 'lsdvmn', 'load_step', 'loadID', 'load_factor', 'loadIDs']:
                column_names.append('LoadStep')
                column_values.append(times)
            elif name == 'node_id':
                column_names.append('NodeID')
                column_values.append(times)
            elif name == 'element_id':
                column_names.append('ElementID')
                column_values.append(times)
            else:  # pragma: no cover
                msg = 'build_dataframe; name=%r' % name
                print(msg)
                raise NotImplementedError(msg)
        assert len(column_names) > 0, column_names
        assert len(column_names) == len(column_values), 'names=%s values=%s' % (column_names, column_values)
        assert len(self.get_headers()) == self.data.shape[-1], 'headers=%s; n=%s\ndata.headers=%s' % (self.get_headers(), len(self.get_headers()), self.data.shape[-1])
        return column_names, column_values

    def _write_table_header(self, op2_file, fascii,
                            date: tuple[int, int, int],
                            include_date: bool=True,
                            subtable_name_default: bytes=b'OUG1    ') -> None:
        try:
            subtable_name = self.subtable_name
        except AttributeError:
            subtable_name = subtable_name_default
            #print('attrs =', self.object_attributes())
            #raise
            pass
        _write_table_header(op2_file, fascii, date, self.table_name, subtable_name,
                            include_date=include_date)


def real_modes_to_omega_freq(eigns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    omega is also called radians
    cycles is also called frequency
    """
    omega_radians = np.sqrt(np.abs(eigns))
    abs_freqs = omega_radians / (2 * np.pi)
    return omega_radians, abs_freqs


def complex_damping_frequency(eigr: np.ndarray,
                              eigi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """eigenvalue = eigr + eigi*1j"""
    assert isinstance(eigr, np.ndarray), eigr
    eigr[eigr == -0.] = 0.
    eigi[eigi == -0.] = 0.
    damping = np.zeros(len(eigr), dtype=eigr.dtype)
    if 0:  # pragma: no cover
        denom = np.sqrt(eigr ** 2 + eigi ** 2)
        inonzero = np.where(denom != 0)[0]
        if len(inonzero):
            damping[inonzero] = -eigr[inonzero] / denom[inonzero]

        # not sure
        abs_freqs = np.sqrt(np.abs(eigi)) / (2 * np.pi)
    else:
        # flutter
        # eig = omega*zeta + omega*1j = eigr + eigi*1j
        # freq = eigi/(2*pi)
        # zeta = 2*eigr/eigi
        abs_freqs = abs(eigi) / (2 * np.pi)
        inonzero = np.where(eigi != 0)[0]
        if len(inonzero):
            damping[inonzero] = 2 * eigr[inonzero] / eigi[inonzero]
    return damping, abs_freqs


def _write_table_header(op2_file, fascii,
                        date: Date,
                        table_name: str,
                        subtable_name: bytes,
                        include_date: bool=True) -> None:
    endian = b'<'
    table_name = '%-8s' % table_name  # 'BOUGV1  '
    fascii.write(f'{table_name}._write_table_header\n')
    # get_nmarkers- [4, 0, 4]
    # marker = [4, 2, 4]
    # table_header = [8, 'BOUGV1  ', 8]
    write_table_header(op2_file, fascii, table_name)

    # read_markers -> [4, -1, 4]
    # get_nmarkers- [4, 0, 4]
    # read_record - marker = [4, 7, 4]
    # read_record - record = [28, recordi, 28]

    # write_markers(op2_file, fascii, '  %s header1a' % table_name, [-1, 0, 7])
    data_a = [4, -1, 4,]
    # data_a = []
    # data_b = [4, -1, 4,]
    data_c = [4, 7, 4,]
    fmt_header = endian + b'6i'
    data = data_a + data_c
    op2_file.write(pack(fmt_header, *data))

    #-----------------
    table1_fmt = endian + b'9i'
    table1 = [
        28,
        102, 0, 0, 0, 512, 0, 0,
        28,
    ]

    blank = ' ' * len(table_name)
    fascii.write(f'{table_name} header1a_i = {data_a}\n')
    #fascii.write(f'{blank}            = {data_b}\n')
    fascii.write(f'{blank}            = {data_c}\n')

    fascii.write(f'{table_name} header1b = {table1}\n')
    op2_file.write(pack(table1_fmt, *table1))

    #recordi = [subtable_name, month, day, year, 0, 1]

    data = [
        4, -2, 4,
        4, 1, 4,
        4, 0, 4,
    ]
    fascii.write('%s header2a = %s\n' % (table_name, data))
    op2_file.write(pack(endian + b'9i', *data))

    month, day, year = date
    dyear = year - 2000
    if subtable_name and include_date:
        table2 = [
            4, 7, 4,
            28,  # 4i -> 13i
            # subtable,todays date 3/6/2014, 0, 1  ( year=year-2000)
            b'%-8s' % subtable_name, month, day, dyear, 0, 1,
            28,
        ]
        table2_format = '4i 8s 6i'
    elif subtable_name:
        table2 = [
            4, 2, 4,
            8,
            b'%-8s' % subtable_name,
            8,
        ]
        table2_format = '4i 8s i'
    else:
        assert include_date is True, include_date
        table2 = [
            4, 7, 4,
            28,  # 4i -> 13i
            # todays date 3/6/2014, 0, 1  ( year=year-2000)
            month, day, dyear, 0, 1,
            28,
        ]
        table2_format = '4i 6i'

    fascii.write('%s header2b = %s\n' % (table_name, table2))
    op2_file.write(pack(table2_format, *table2))


def get_sort_element_sizes(self, debug: bool=False) -> tuple[int, int, int]:
    if self.is_sort1:
        ntimes = self.ntimes
        nelements = self.nelements
        ntotal = self.ntotal
        #nx = ntimes
        #ny = nnodes
        if debug:
            print(f'SORT1 {self.__class__.__name__} ntimes={ntimes} nelements={nelements}')
    elif self.is_sort2:
        # flip this to sort1
        ntimes = self.ntotal
        nelements = self.ntimes
        ntotal = self.nelements
        #nx = ntimes
        #ny = nnodes
        if debug:
            print(f'***SORT2 {self.__class__.__name__} ntimes={ntimes} nelements={nelements} ntotal={ntotal}')
    else:
        raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())
    #assert nelements == ntotal or ntimes == ntotal, (ntimes, nelements, ntotal)
    return ntimes, nelements, ntotal


def get_sort_node_sizes(self, debug: bool=False) -> tuple[int, int, int]:
    if self.is_sort1:
        #self._nnodes //= self.ntimes
        ntimes = self.ntimes
        nnodes = self._nnodes // self.ntimes
        ntotal = self.ntotal
        #nx = ntimes
        #ny = nnodes
        #print("SORT1 ntimes=%s nelements=%s" % (ntimes, nelements))
    elif self.is_sort2:
        #dt=0.0 nid=3306 itime=0/5=5 inode=0/5=20
        # flip this to sort1
        if debug:
            print("***SORT2 ntotal=%s _nnodes=%s ntimes=%s" % (self.ntotal, self._nnodes, self.ntimes))
        ntimes = self.ntotal
        nnodes = self._nnodes // self.ntotal
        ntotal = nnodes
        if debug:
            print(f'***SORT2 {self.__class__.__name__} ntimes={ntimes} nnodes={nnodes} ntotal={ntotal}')
        #nx = ntimes
        #ny = nnodes
        #print("***SORT2 ntotal=%s nnodes=%s ntimes=%s" % (ntotal, nnodes, ntimes))
    else:
        raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())
    #assert nnodes == ntotal or ntimes == ntotal, (ntimes, nnodes, ntotal)
    return ntimes, nnodes, ntotal


class BaseElement(ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        #--------------------------------
        # TODO: remove ???
        #self.element_name = None
        #self.element = None
        #self.element_node = None
        #self.element_type = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    def _eq_header(self, table):
        ScalarObject._eq_header(self, table)
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if hasattr(self, 'element_name'):
            if self.nonlinear_factor not in (None, np.nan) and not is_nan:
                assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                    self.element_name, self.element_type, self._times, table._times)

        log = SimpleLogger()
        _check_element(self, table, log)
        _check_element_node(self, table, log)

    def _build_pandas_transient_elements(self, column_values: np.ndarray, column_names: list[str],
                                         headers: list[str],
                                         element: np.ndarray, data: np.ndarray):
        """common method to build a transient dataframe"""
        import pandas as pd
        columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

        eid_item = []
        for eid in element:
            for header in headers:
                eid_item.append([eid, header])
        ntimes, nelements = data.shape[:2]

        nheaders = len(headers)
        if self.table_name in ['OEFPSD1']:
            ifilter = ~np.isnan(np.max(data[-1, :, :], axis=1))
            if ifilter.sum():
                warnings.warn(f'filtering NaNs from {self.class_name}')
                data2 = data[:, ifilter, :]
                nelements = data2.shape[1]
                A = data2.reshape(ntimes, nelements*nheaders).T
            else:
                A = data.reshape(ntimes, nelements*nheaders).T
        else:
            A = data.reshape(ntimes, nelements*nheaders).T

        names = ['ElementID', 'Item']
        index = pd.MultiIndex.from_tuples(eid_item, names=names)
        try:
            data_frame = pd.DataFrame(A, columns=columns, index=index)
        except ValueError:
            print('A.shape =', A.shape)
            print('len(element) =', len(element))
            print('columns =', columns)
            raise
        # old
        #data_frame = pd.Panel(data, items=column_values, major_axis=element, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'Item']
        #print(data_frame)
        return data_frame

    def _build_pandas_transient_element_node(self, column_values, column_names, headers,
                                             element_node, data, names=None,
                                             from_tuples=True, from_array=False):
        """common method to build a transient dataframe"""
        # Freq                  0.00001  10.00000 20.00000 30.00000                 40.00000 50.00000 60.00000
        # ElementID NodeID Item
        # 1         0      oxx        0j       0j       0j       0j    (3200.0806+6017.714j)       0j       0j
        #                  oyy        0j       0j       0j       0j    (410.68146+772.2816j)       0j       0j
        #                  ozz        0j       0j       0j       0j    (0.306115+0.5756457j)       0j       0j
        #                  txy        0j       0j       0j       0j  (-120.69606-226.96753j)       0j       0j
        #                  tyz        0j       0j       0j       0j  (0.70554054+1.3267606j)       0j       0j
        #                  txz        0j       0j       0j       0j     (5193.834+9766.943j)       0j       0j
        # 2                oxx        0j       0j       0j       0j    (8423.371+15840.051j)       0j       0j
        #                  oyy        0j       0j       0j       0j    (-3364.359-6326.637j)       0j       0j
        #                  ozz        0j       0j       0j       0j  (-74931.664-140908.11j)       0j       0j
        #                  txy        0j       0j       0j       0j  (-261.20972-491.20178j)       0j       0j
        #                  tyz        0j       0j       0j       0j   (121.57285+228.61633j)       0j       0j
        #                  txz        0j       0j       0j       0j     (5072.678+9539.112j)       0j       0j
        import pandas as pd
        columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

        #print(data.shape)
        ntimes, nelements = data.shape[:2]
        nheaders = len(headers)
        try:
            A = data.reshape(ntimes, nelements*nheaders).T
        except ValueError:  # pragma: no cover
            ntotal = ntimes * nelements * nheaders
            print(f'data.shape={data.shape}; ntimes={ntimes} nelements={nelements} nheaders={nheaders}; ntotal={ntotal}')
            raise

        if names is None:
            names = ['ElementID', 'NodeID', 'Item']

        assert not(from_tuples and from_array)
        if from_tuples:
            nvars = element_node.shape[1]
            assert len(names) == nvars + 1, f'names={names} element_node={element_node} {element_node.shape}'
            eid_nid_item = []
            for eid, nid in element_node:
                for header in headers:
                    eid_nid_item.append([eid, nid, header])
            index = pd.MultiIndex.from_tuples(eid_nid_item, names=names)
        elif from_array:
            nvars = len(element_node)
            assert len(names) == nvars + 1, f'names={names} element_node={element_node} (n={len(element_node)})'
            eid_nid_item = []
            for eid in element_node:
                eidi = np.vstack([eid]*nheaders)
                eid_nid_item.append(eidi.ravel())
            eid_nid_item.append(headers * nelements)
            index = pd.MultiIndex.from_arrays(eid_nid_item, names=names)
        else:  # pragma: no cover
            raise RuntimeError('from_tuple, from_array')
        data_frame = pd.DataFrame(A, columns=columns, index=index)

        #element_node = [element_node[:, 0], element_node[:, 1]]
        #data_frame = pd.Panel(data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        #print(data_frame)
        return data_frame


def get_times_dtype(nonlinear_factor: int | float, size: int,
                    analysis_code_fmt=None) -> tuple[str, str, str]:
    dtype = 'float'
    if isinstance(nonlinear_factor, integer_types):
        dtype = 'int'

    if size == 4:
        dtype += '32'
        fdtype = 'float32'
        idtype = 'int32'
    else:
        dtype += '64'
        fdtype = 'float64'
        idtype = 'int64'

    if analysis_code_fmt:
        dtype = analysis_code_fmt
        return dtype, idtype, fdtype
    return dtype, idtype, fdtype


def get_complex_times_dtype(size: int) -> tuple[str, str]:
    # assert isinstance()
    # dtype = 'float'
    # if isinstance(nonlinear_factor, integer_types):
    #     dtype = 'int'

    if size == 4:
        # dtype += '32'
        cfdtype = 'complex64'
        idtype = 'int32'
    else:
        # dtype += '64'
        cfdtype = 'complex128'
        idtype = 'int64'
    return idtype, cfdtype


def _check_element(table1: BaseElement, table2: BaseElement, log: SimpleLogger) -> None:
    """checks the ``element_node`` variable"""
    if not hasattr(table1, 'element'):
        return
    if table1.element is None:
        return

    if not np.array_equal(table1.element, table2.element):
        assert table1.element.shape == table2.element.shape, 'shape=%s element.shape=%s' % (
            table1.element.shape, table2.element.shape)
        msg = f'table_name={table1.table_name!r} class_name={table1.__class__.__name__}\n'
        msg += '%s\nEid\n' % str(table1.code_information())
        for eid1, eid2 in zip(table1.element, table2.element):
            msg += '%s, %s\n' % (eid1, eid2)
        print(msg)
        raise ValueError(msg)

    element = table1.element
    try:
        eid_min = element.min()
    except TypeError:
        # strain energy element names can be U8 strings
        if element.dtype.type != np.str_:
            print(element, element.dtype)
            raise
        return
    nshape = len(element.shape)
    if eid_min <= 0:
        if nshape == 1:
            msg = (f'{table1}\ntable_name={table1.table_name}\n'
                   f'eids={element}.min = {eid_min}')
            log.error(msg)
        else:
            if table1.table_name not in ['ONRGY1', 'ONRGY2', 'OEKE1']:
                msg = f'table_name = {table1.table_name}\n'
                for i, eidsi in enumerate(element):
                    eid_min = eidsi.min()
                    if eid_min <= 0:
                        msg += f'{table1}\neids[{i}]={eidsi}.min = {eid_min}\n'
                log.error(msg)


def _check_element_node(table1: BaseElement, table2: BaseElement, log: SimpleLogger) -> None:
    """checks the ``element_node`` variable"""
    if not hasattr(table1, 'element_node'):
        return
    if table1.element_node is None:
        return

    if not np.array_equal(table1.element_node, table2.element_node):
        if table1.element_node.shape != table2.element_node.shape:
            msg = (
                f'{table1}\n'
                f'table1.element_node.shape={table1.element_node.shape} '
                f'table2.element_node.shape={table2.element_node.shape}')

            print(table1.element_node.tolist())
            print(table2.element_node.tolist())
            raise ValueError(msg)
        msg = f'table_name={table1.table_name!r} class_name={table1.__class__.__name__}\n'
        msg += '%s\n' % str(table1.code_information())
        for i, (eid1, nid1), (eid2, nid2) in zip(count(), table1.element_node, table2.element_node):
            msg += '%s : (%s, %s), (%s, %s)\n' % (i, eid1, nid1, eid2, nid2)
            if i > 20:
                msg += '...'
                break
        print(msg)
        raise ValueError(msg)

    eids = table1.element_node[:, 0]
    if eids.min() <= 0:
        print(table1.element_node)
        print(table2.element_node)
        log.error(f'{table1}\neids={eids}.min={eids.min()}; n={len(eids)}')
        raise ValueError(f'{table1}\neids={eids}.min={eids.min()}; n={len(eids)}')


def set_as_sort1(obj):
    #print('set_as_sort1: table_name=%r' % obj.table_name)
    if obj.is_sort1:
        if obj.analysis_code == 1:
            pass
        else:
            name = obj.name
            setattr(obj, name + 's', obj._times)
            # if name == 'mode':
            #     print('obj._times', 'modes', obj._times)
        return

    # sort2
    if obj.analysis_code == 1:
        # static...because reasons
        analysis_method = 'N/A'
    else:
        try:
            analysis_method = obj.analysis_method
        except AttributeError:
            print(obj.code_information())
            print(obj.object_stats())
            raise

    # print(obj.get_stats())
    # print(obj.data.shape)
    obj.sort_method = 1
    obj.sort_bits[1] = 0
    bit0, bit1, bit2 = obj.sort_bits
    obj.table_name = SORT2_TABLE_NAME_MAP[obj.table_name]
    obj.sort_code = bit0 + 2*bit1 + 4*bit2
    print(obj.code_information())
    # assert obj.is_sort1

    if analysis_method != 'N/A':
        obj.data_names[0] = analysis_method
        #print(obj.table_name_str, analysis_method, obj._times)
        setattr(obj, obj.analysis_method + 's', obj._times)

        #print('set_as_sort1: table_name=%r' % self.table_name)
        # dt
        obj.data_code['name'] = obj.analysis_method
        #print(self.get_stats())
        obj.name = obj.analysis_method
        del obj.analysis_method
