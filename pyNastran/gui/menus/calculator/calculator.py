from __future__ import annotations
"""
Calculator

+-----+-------------+
|     |             |
| Res |       Funcs |
|     |             |
+-----+-------------+
"""
"""
defines:
 - CalculatorWindow

"""
"""
defines:
 - AnimationWindow

"""
import os
from typing import TYPE_CHECKING

import numpy
#np_funcs2 = []
#np_linalg_funcs2 = []
import scipy, scipy.interpolate

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QPushButton, QLineEdit,
    #QWidget, QButtonGroup,
    QGridLayout, QHBoxLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
    QCheckBox, QGroupBox, QComboBox, QFileDialog)
from qtpy.compat import getexistingdirectory

from pyNastran.gui.menus.python_console import PythonConsoleWidget, PythonConsoleLayout
from pyNastran.gui.utils.qt.pydialog import PyDialog, set_combo_box_text
from pyNastran.gui.utils.qt.checks.qlineedit import (
    check_int, check_float, check_name_str, check_path)
from pyNastran.gui.utils.qt.dialogs import open_file_dialog
from pyNastran.gui.menus.results_sidebar import ResultsWindow
from pyNastran.gui.menus.results_sidebar_utils import (
    get_cases_from_tree, #build_pruned_tree
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow

IS_RESULTS_SELECTOR = True
class CalculatorWindow(PyDialog):
    """
    +-------------------+
    | Animation         |
    +-------------------------+
    | icase   ______          |
    | scale   ______  Default |
    | time    ______  Default |
    |                         |
    | nframes ______  Default |
    | resolu. ______  Default |
    | Dir     ______  Browse  |
    | iFrame  ______          |
    |                         |
    | Animations:             |
    | o Scale, Phase, Time    |
    |                         |
    | x delete images         |
    | x repeat                |  # TODO: change to an integer
    | x make gif              |
    |                         |
    |      Step, RunAll       |
    |         Close           |
    +-------------------------+

    TODO: add key-frame support
    """
    def __init__(self, data, win_parent=None, modules_form=None, fringe_cases=None, is_gui_parent=False):
        PyDialog.__init__(self, data, win_parent)

        # is the parent the gui?
        self.is_gui_parent = False
        self.module_map = modules_form[0]
        self.modules_form = modules_form[1]
        assert isinstance(self.module_map, dict), self.module_map
        assert isinstance(self.modules_form, list), self.modules_form
        self.fringe_cases = fringe_cases
        self.set_font_size(data['font_size'])
        self.istep = 0
        self._animate_type = 'time'

        self._updated_animation = False
        self._active_deformation = 0.
        self._icase_fringe = data['icase_fringe']
        self._icase_disp = data['icase_disp']
        self._icase_vector = data['icase_vector']

        self._default_title = data['title']
        self._default_time = data['time']
        self._default_fps = data['frames/sec']
        self._default_resolution = data['resolution']

        self._scale = data['scale']
        self._default_scale = data['default_scale']
        self._default_is_scale = data['is_scale']

        self._arrow_scale = data['arrow_scale']
        self._default_arrow_scale = data['default_arrow_scale']

        self._phase = data['phase']
        self._default_phase = data['default_phase']

        self._default_dirname = data['dirname']
        self._default_gif_name = os.path.join(self._default_dirname, data['title'] + '.gif')

        self.animation_types = [
            'Animate Scale',
        ]
            #'Animate Phase',
            #'Animate Time',
            #'Animate Frequency Sweep'
        #]

        self.setWindowTitle('Calculator')
        self.create_widgets()
        self.create_layout()
        self.set_connections()

        self.is_gui = False
        self.gui = None  # type: MainWindow
        if hasattr(self.win_parent, '_updated_legend'):
            self.win_parent.is_animate_open = True
            self.is_gui = True
            self.gui = self.win_parent.win_parent

        icase_max = 1000
        if is_gui_parent:
            self.is_gui = True
            self.gui = self.win_parent
            icase_max = max(self.gui.result_cases)  # TODO: update 1000

        #self.icase_fringe_edit.setRange(0, icase_max)
        #self.icase_disp_edit.setRange(1, icase_max)
        #self.icase_vector_edit.setRange(1, icase_max)


    def create_widgets(self):
        """creates the menu objects"""
        self.code_window = PythonConsoleLayout(self)
        #widget = QWidget(self)
        #horizontal_vertical_group = QButtonGroup(widget)
        #horizontal_vertical_group.addButton(self.animate_scale_radio)
        #horizontal_vertical_group.addButton(self.animate_phase_radio)
        #horizontal_vertical_group.addButton(self.animate_time_radio)
        #horizontal_vertical_group.addButton(self.animate_freq_sweeep_radio)

        # bottom buttons
        self.cancel_button = QPushButton("Close")

    def set_connections(self):
        """creates the actions for the menu"""
        self.cancel_button.clicked.connect(self.on_cancel)

    def on_default_title(self):
        """sets the default gif name"""
        self.gif_edit.setText(self._default_title + '.gif')

    def create_layout(self):
        """displays the menu objects"""
        #----------
        # bottom buttons
        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.cancel_button)

        label = QLabel('Python Console:')

        vbox = QVBoxLayout()
        #vbox.addWidget(self.code_window)
        vbox.addWidget(label)
        vbox.addLayout(self.code_window)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        if IS_RESULTS_SELECTOR and self.fringe_cases:
            cases = get_cases_from_tree(self.modules_form)
            parent = self
            name = 'modules'
            data = self.modules_form
            choices = cases
            def on_print(icase):
                module, func_str = self.module_map[icase]
                func = getattr(module, func_str)
                print(func.__doc__)
                print('-------------------------------')
                print(f'on_print(icase={icase}) -> {func_str}')
                try:

                    args = func.__code__.co_varnames
                    #('x', 'y')
                    arg_word = [', '.join(args)]
                    word = f'{module.__name__}.{func_str}({arg_word})'
                    self.code_window.enter_data.add(word)
                    #print(f'on_print(icase={icase}) -> {outi}')
                except AttributeError:
                    word = f'{module.__name__}.{func_str}()'
                    self.code_window.enter_data.add(word)


            module_actions = [
                ## (right_click_msg, callback, validate?)
                ##('Clear Results...', self.on_clear_results, False),
                ##('Apply Results to Fringe...', 'fringe', self.on_fringe, True),
                ##('Apply Results to Displacement...', self.on_disp, True),
                ##('Apply Results to Vector...', self.on_vector, True),
                ('Print...', on_print, True),
                #('Modify...', on_modify, True),
            ]
            def on_add(icase):
                pass
            results_actions = [
                ## (right_click_msg, callback, validate?)
                ##('Clear Results...', self.on_clear_results, False),
                ##('Apply Results to Fringe...', 'fringe', self.on_fringe, True),
                ##('Apply Results to Displacement...', self.on_disp, True),
                ##('Apply Results to Vector...', self.on_vector, True),
                ('Add...', on_add, True),
            ]

            self.modules_widget = ResultsWindow(
                parent, name, data, choices,
                is_single_select=True,
                left_click_callback=None,
                right_click_actions=module_actions,
                include_export_case=False,
                include_clear=False,
                include_delete=False)
            #---------------
            cases = get_cases_from_tree(self.fringe_cases)
            parent = self
            name = 'main'
            data = self.fringe_cases
            choices = cases
            self.results_widget = ResultsWindow(
                parent, name, data, choices,
                is_single_select=True,
                left_click_callback=None,
                right_click_actions=results_actions,
                #include_export_case=True,
                #include_delete=True,
                include_clear=False, include_delete=False)
            #-----------------------------------------------------------------------------
            vbox_modules = QVBoxLayout()
            results_widget_label = QLabel('Modules:')
            vbox_modules.addWidget(results_widget_label)
            vbox_modules.addWidget(self.modules_widget)

            vbox_results = QVBoxLayout()
            results_widget_label = QLabel('Results:')
            vbox_results.addWidget(results_widget_label)
            vbox_results.addWidget(self.results_widget)
            hbox_main = QHBoxLayout()
            hbox_main.addLayout(vbox_modules)
            hbox_main.addLayout(vbox)
            hbox_main.addLayout(vbox_results)
            self.setLayout(hbox_main)
        else:
            self.setLayout(vbox)

    #def on_fringe(self, icase):
        #"""sets the icase fringe"""
        #pass

    #def on_disp(self, icase):
        #"""sets the icase disp"""
        #pass

    #def on_vector(self, icase):
        #"""sets the icase vector"""
        #pass

    #def on_clear_results(self):
        #"""sink for the right click menu"""
        #pass

    def _on_execute_python_button(self, clear=False):
        pass

    def on_validate(self, wipe=False):
        """checks to see if the input is valid"""
        # requires no special validation
        icase_fringe, flag0 = check_int(self.icase_fringe_edit)
        icase_disp, unused_flaga = check_int(self.icase_disp_edit)
        icase_vector, unused_flagb = check_int(self.icase_vector_edit)
        #icase_disp = self._icase_disp
        #icase_vector = self._icase_vector

        scale, flag1 = check_float(self.scale_edit)
        time, flag2 = check_float(self.time_edit)
        fps, flag3 = check_float(self.fps_edit)

        min_value = max_value = None
        flag4 = flag5 = True
        if self.min_value_edit.isEnabled():
            min_value, flag4 = check_float(self.min_value_edit)
        if self.max_value_edit.isEnabled():
            max_value, flag5 = check_float(self.max_value_edit)

        if wipe:
            animate_in_gui = False
            scale = 0.
            flag1 = True
        else:
            animate_in_gui = self.animate_in_gui_checkbox.isChecked()
            if scale == 0.0:
                self.scale_edit.setStyleSheet("QLineEdit{background: red;}")
                flag1 = False

        if animate_in_gui or wipe:
            passed = all([flag0, flag1, flag2, flag3, flag4, flag5])
            magnify, output_dir, gifbase = None, None, None
        else:
            magnify, flag6 = check_int(self.resolution_edit)
            output_dir, flag7 = check_path(self.browse_folder_edit)
            gifbase, flag8 = check_name_str(self.gif_edit)
            passed = all([flag0, flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8])
        return passed, (icase_fringe, icase_disp, icase_vector, scale, time, fps, animate_in_gui,
                        magnify, output_dir, gifbase, min_value, max_value)

    #def on_ok(self):
        #"""click the OK button"""
        #passed = self.on_apply()
        #if passed:
            #self.win_parent._animation_window_shown = False
            #self.close()
            ##self.destroy()

    def on_cancel(self):
        """click the Cancel button"""
        self.out_data['close'] = True
        self.close()

def enable_disable_objects(qt_objects, enable=True):
    for obj in qt_objects:
        obj.setEnabled(enable)

def load_module_functions(module, func_list, exclusion_set):
    for func in dir(module):
        skip = (
            func.startswith('_') or
            func.endswith('_') or
            func.upper() == func or
            func.endswith(('Warning', 'Error')) or
            func in exclusion_set)
        if skip:
            continue
        func_list.append(func)

def main(): # pragma: no cover
    """test example for CalculatorWindow"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    #from pyNastran.gui.menus.legend.animation import AnimationWindow
    data2 = {
        'font_size' : 8,
        'icase_fringe' : 1,
        'icase_disp' : 2,
        'icase_vector' : 3,

        'title' : 'cat',
        'time' : 2,
        'frames/sec' : 30,
        'resolution' : 1,
        'iframe' : 0,
        'is_scale' : False,
        'dirname' : os.getcwd(),
        'scale' : 2.0,
        'default_scale' : 10,

        'arrow_scale' : 3.0,
        'default_arrow_scale' : 30,

        #'phase' : 0.,
        'phase' : None,
        'default_phase' : 120.,
        #'default_phase' : None,

        #'start_time' : 0.,
        #'end_time' : 0.,
        'default_time' : 0.,
        'icase_start' : 10,
        'icase_delta' : 3,
        'stress_min' : 0.,
        'stress_max' : 1000.,
    }
    data2['phase'] = 0.  # uncomment for phase

    form = [
        ['Geometry', None, [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('PropertyID', 2, []),
            ('MaterialID', 3, []),
            ('E', 4, []),
            ('Element Checks', None, [
                ('ElementDim', 5, []),
                ('Min Edge Length', 6, []),
                ('Min Interior Angle', 7, []),
                ('Max Interior Angle', 8, [])],
            ),],
        ],
    ]
    np_funcs = []
    #np_linalg_funcs = []

    modules = [
        ('numpy', numpy, [], []),
        ('numpy.linalg', numpy.linalg, [], []),
        ('scipy', scipy, [], []),
        ('scipy.interpolate', scipy.interpolate, [], []),
        ('scipy.special', scipy.special, [], []),
        #('scipy', scipy, [], []),
        #('scipy', scipy, [], []),
        #('numpy', numpy),
        #('numpy', numpy),
        #('numpy', numpy),
        #('numpy', numpy),
    ]

    #modules_form = [
        #['numpy', None, np_funcs2],
        #['numpy.linalg', None, np_linalg_funcs2],
    #]

    #from pyNastran.utils import object_methods
    ['ALLOW_THREADS', 'AxisError', 'BUFSIZE', 'Bytes0', 'CLIP', 'ComplexWarning', 'DataSource', 'Datetime64', 'ERR_CALL',
     'ERR_DEFAULT', 'ERR_IGNORE', 'ERR_LOG', 'ERR_PRINT', 'ERR_RAISE', 'ERR_WARN', 'FLOATING_POINT_SUPPORT', 'FPE_DIVIDEBYZERO', 'FPE_INVALID', 'FPE_OVERFLOW', 'FPE_UNDERFLOW', 'False_', 'Inf', 'Infinity', 'MAXDIMS', 'MAY_SHARE_BOUNDS', 'MAY_SHARE_EXACT', 'MachAr', 'ModuleDeprecationWarning', 'NAN', 'NINF', 'NZERO', 'NaN', 'PINF', 'PZERO', 'RAISE', 'RankWarning', 'SHIFT_DIVIDEBYZERO', 'SHIFT_INVALID', 'SHIFT_OVERFLOW', 'SHIFT_UNDERFLOW', 'ScalarType', 'Str0', 'Tester', 'TooHardError', 'True_', 'UFUNC_BUFSIZE_DEFAULT', 'UFUNC_PYVALS_NAME', 'Uint64', 'VisibleDeprecationWarning', 'WRAP', '_NoValue', '_UFUNC_API', '__NUMPY_SETUP__', '__all__', '__builtins__', '__cached__', '__config__', '__deprecated_attrs__', '__dir__', '__doc__', '__expired_functions__', '__file__', '__getattr__', '__git_revision__', '__loader__', '__name__', '__package__', '__path__', '__spec__', '__version__', '_add_newdoc_ufunc', '_builtins', '_distributor_init', '_financial_names', '_globals', '_mat', '_pytesttester', 'abs', 'absolute', 'add', 'add_docstring', 'add_newdoc', 'add_newdoc_ufunc', 'alen', 'all', 'allclose', 'alltrue', 'amax', 'amin', 'angle', 'any', 'append', 'apply_along_axis', 'apply_over_axes', 'arange', 'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan', 'arctan2', 'arctanh', 'argmax', 'argmin', 'argpartition', 'argsort', 'argwhere', 'around', 'array', 'array2string', 'array_equal', 'array_equiv', 'array_repr', 'array_split', 'array_str', 'asanyarray', 'asarray', 'asarray_chkfinite', 'ascontiguousarray', 'asfarray', 'asfortranarray', 'asmatrix', 'asscalar', 'atleast_1d', 'atleast_2d', 'atleast_3d', 'average', 'bartlett', 'base_repr', 'binary_repr', 'bincount', 'bitwise_and', 'bitwise_not', 'bitwise_or', 'bitwise_xor', 'blackman', 'block', 'bmat', 'bool8', 'bool_', 'broadcast', 'broadcast_arrays', 'broadcast_shapes', 'broadcast_to', 'busday_count', 'busday_offset', 'busdaycalendar', 'byte', 'byte_bounds', 'bytes0', 'bytes_', 'c_', 'can_cast', 'cast', 'cbrt', 'cdouble', 'ceil', 'cfloat', 'char', 'character', 'chararray', 'choose', 'clip', 'clongdouble', 'clongfloat', 'column_stack', 'common_type', 'compare_chararrays', 'compat', 'complex128', 'complex64', 'complex_', 'complexfloating', 'compress', 'concatenate', 'conj', 'conjugate', 'convolve', 'copy', 'copysign', 'copyto', 'core', 'corrcoef', 'correlate', 'cos', 'cosh', 'count_nonzero', 'cov', 'cross', 'csingle', 'ctypeslib', 'cumprod', 'cumproduct', 'cumsum', 'datetime64', 'datetime_as_string', 'datetime_data', 'deg2rad', 'degrees', 'delete', 'deprecate', 'deprecate_with_doc', 'diag', 'diag_indices', 'diag_indices_from', 'diagflat', 'diagonal', 'diff', 'digitize', 'disp', 'divide', 'divmod', 'dot', 'double', 'dsplit', 'dstack', 'dtype', 'e', 'ediff1d', 'einsum', 'einsum_path', 'emath', 'empty', 'empty_like', 'equal', 'errstate', 'euler_gamma', 'exp', 'exp2', 'expand_dims', 'expm1', 'extract', 'eye', 'fabs', 'fastCopyAndTranspose', 'fft', 'fill_diagonal', 'find_common_type', 'finfo', 'fix', 'flatiter', 'flatnonzero', 'flexible', 'flip', 'fliplr', 'flipud', 'float16', 'float32', 'float64', 'float_', 'float_power', 'floating', 'floor', 'floor_divide', 'fmax', 'fmin', 'fmod', 'format_float_positional', 'format_float_scientific', 'format_parser', 'frexp', 'frombuffer', 'fromfile', 'fromfunction', 'fromiter', 'frompyfunc', 'fromregex', 'fromstring', 'full', 'full_like', 'gcd', 'generic', 'genfromtxt', 'geomspace', 'get_array_wrap', 'get_include', 'get_printoptions', 'getbufsize', 'geterr', 'geterrcall', 'geterrobj', 'gradient', 'greater', 'greater_equal', 'half', 'hamming', 'hanning', 'heaviside', 'histogram', 'histogram2d', 'histogram_bin_edges', 'histogramdd', 'hsplit', 'hstack', 'hypot', 'i0', 'identity', 'iinfo', 'imag', 'in1d', 'index_exp', 'indices', 'inexact', 'inf', 'info', 'infty', 'inner', 'insert', 'int0', 'int16', 'int32', 'int64', 'int8', 'int_', 'intc', 'integer', 'interp', 'intersect1d', 'intp', 'invert', 'is_busday', 'isclose', 'iscomplex', 'iscomplexobj', 'isfinite', 'isfortran', 'isin', 'isinf', 'isnan', 'isnat', 'isneginf', 'isposinf', 'isreal', 'isrealobj', 'isscalar', 'issctype', 'issubclass_', 'issubdtype', 'issubsctype', 'iterable', 'ix_', 'kaiser', 'kron', 'lcm', 'ldexp', 'left_shift', 'less', 'less_equal', 'lexsort', 'lib', 'linalg', 'linspace', 'little_endian', 'load', 'loads', 'loadtxt', 'log', 'log10', 'log1p', 'log2', 'logaddexp', 'logaddexp2', 'logical_and', 'logical_not', 'logical_or', 'logical_xor', 'logspace', 'longcomplex', 'longdouble', 'longfloat', 'longlong', 'lookfor', 'ma', 'mafromtxt', 'mask_indices', 'mat', 'math', 'matmul', 'matrix', 'matrixlib', 'max', 'maximum', 'maximum_sctype', 'may_share_memory', 'mean', 'median', 'memmap', 'meshgrid', 'mgrid', 'min', 'min_scalar_type', 'minimum', 'mintypecode', 'mod', 'modf', 'moveaxis', 'msort', 'multiply', 'nan', 'nan_to_num', 'nanargmax', 'nanargmin', 'nancumprod', 'nancumsum', 'nanmax', 'nanmean', 'nanmedian', 'nanmin', 'nanpercentile', 'nanprod', 'nanquantile', 'nanstd', 'nansum', 'nanvar', 'nbytes', 'ndarray', 'ndenumerate', 'ndfromtxt', 'ndim', 'ndindex', 'nditer', 'negative', 'nested_iters', 'newaxis', 'nextafter', 'nonzero', 'not_equal', 'numarray', 'number', 'obj2sctype', 'object0', 'object_', 'ogrid', 'oldnumeric', 'ones', 'ones_like', 'os', 'outer', 'packbits', 'pad', 'partition', 'percentile', 'pi', 'piecewise', 'place', 'poly', 'poly1d', 'polyadd', 'polyder', 'polydiv', 'polyfit', 'polyint', 'polymul', 'polynomial', 'polysub', 'polyval', 'positive', 'power', 'printoptions', 'prod', 'product', 'promote_types', 'ptp', 'put', 'put_along_axis', 'putmask', 'quantile', 'r_', 'rad2deg', 'radians', 'random', 'ravel', 'ravel_multi_index', 'real', 'real_if_close', 'rec', 'recarray', 'recfromcsv', 'recfromtxt', 'reciprocal', 'record', 'remainder', 'repeat', 'require', 'reshape', 'resize', 'result_type', 'right_shift', 'rint', 'roll', 'rollaxis', 'roots', 'rot90', 'round', 'round_', 'row_stack', 's_', 'safe_eval', 'save', 'savetxt', 'savez', 'savez_compressed', 'sctype2char', 'sctypeDict', 'sctypes', 'searchsorted', 'select', 'set_numeric_ops', 'set_printoptions', 'set_string_function', 'setbufsize', 'setdiff1d', 'seterr', 'seterrcall', 'seterrobj', 'setxor1d', 'shape', 'shares_memory', 'short', 'show_config', 'sign', 'signbit', 'signedinteger', 'sin', 'sinc', 'single', 'singlecomplex', 'sinh', 'size', 'sometrue', 'sort', 'sort_complex', 'source', 'spacing', 'split', 'sqrt', 'square', 'squeeze', 'stack', 'std', 'str0', 'str_', 'string_', 'subtract', 'sum', 'swapaxes', 'sys', 'take', 'take_along_axis', 'tan', 'tanh', 'tensordot', 'test', 'testing', 'tile', 'timedelta64', 'trace', 'tracemalloc_domain', 'transpose', 'trapz', 'tri', 'tril', 'tril_indices', 'tril_indices_from', 'trim_zeros', 'triu', 'triu_indices', 'triu_indices_from', 'true_divide', 'trunc', 'typeDict', 'typecodes', 'typename', 'ubyte', 'ufunc', 'uint', 'uint0', 'uint16', 'uint32', 'uint64', 'uint8', 'uintc', 'uintp', 'ulonglong', 'unicode_', 'union1d', 'unique', 'unpackbits', 'unravel_index', 'unsignedinteger', 'unwrap', 'use_hugepage', 'ushort', 'vander', 'var', 'vdot', 'vectorize', 'version', 'void', 'void0', 'vsplit', 'vstack', 'warnings', 'where', 'who', 'zeros', 'zeros_like']

    exclusion_set = set([
        'os', 'sys', 'True_', 'False_', 'Bytes0', 'Tester', 'add_docstring', 'add_newdoc', 'add_newdoc_ufunc',
        'Str0', 'Tester', 'True_', 'Uint64',
        'bytes0', 'bytes_', 'c_', 'can_cast', 'cast', 'cbrt', 'cdouble', 'cfloat', 'char', 'character', 'chararray',
        'clongdouble', 'clongfloat',
        #'choose', 'clip', 'column_stack', 'common_type', 'compare_chararrays', 'compat', 'complex128', 'complex64', 'complex_', 'complexfloating', 'compress', 'concatenate', 'conj', 'conjugate', 'convolve', 'copy', 'copysign', 'copyto', 'core', 'corrcoef', 'correlate', 'cos', 'cosh', 'count_nonzero', 'cov', 'cross', 'csingle', 'ctypeslib', 'cumprod', 'cumproduct', 'cumsum', 'datetime64', 'datetime_as_string', 'datetime_data', 'deg2rad', 'degrees', 'delete', 'deprecate', 'deprecate_with_doc', 'diag', 'diag_indices', 'diag_indices_from', 'diagflat', 'diagonal', 'diff', 'digitize', 'disp', 'divide', 'divmod', 'dot', 'double', 'dsplit', 'dstack', 'dtype', 'e', 'ediff1d', 'einsum', 'einsum_path', 'emath', 'empty', 'empty_like', 'equal', 'errstate', 'euler_gamma', 'exp', 'exp2', 'expand_dims', 'expm1', 'extract', 'eye', 'fabs', 'fastCopyAndTranspose', 'fft', 'fill_diagonal', 'find_common_type', 'finfo', 'fix', 'flatiter', 'flatnonzero', 'flexible', 'flip', 'fliplr', 'flipud', 'float16', 'float32', 'float64', 'float_', 'float_power', 'floating', 'floor', 'floor_divide', 'fmax', 'fmin', 'fmod', 'format_float_positional', 'format_float_scientific', 'format_parser', 'frexp', 'frombuffer', 'fromfile', 'fromfunction', 'fromiter', 'frompyfunc', 'fromregex', 'fromstring', 'full', 'full_like', 'gcd', 'generic', 'genfromtxt', 'geomspace', 'get_array_wrap', 'get_include', 'get_printoptions', 'getbufsize', 'geterr', 'geterrcall', 'geterrobj', 'gradient', 'greater', 'greater_equal', 'half', 'hamming', 'hanning', 'heaviside', 'histogram', 'histogram2d', 'histogram_bin_edges', 'histogramdd', 'hsplit', 'hstack', 'hypot', 'i0', 'identity', 'iinfo', 'imag', 'in1d', 'index_exp', 'indices', 'inexact', 'inf', 'info', 'infty', 'inner', 'insert', 'int0', 'int16', 'int32', 'int64', 'int8', 'int_', 'intc', 'integer', 'interp', 'intersect1d', 'intp', 'invert', 'is_busday', 'isclose', 'iscomplex', 'iscomplexobj', 'isfinite', 'isfortran', 'isin', 'isinf', 'isnan', 'isnat', 'isneginf', 'isposinf', 'isreal', 'isrealobj', 'isscalar', 'issctype', 'issubclass_', 'issubdtype', 'issubsctype', 'iterable', 'ix_', 'kaiser', 'kron', 'lcm', 'ldexp', 'left_shift', 'less', 'less_equal', 'lexsort', 'lib', 'linalg', 'linspace', 'little_endian', 'load', 'loads', 'loadtxt', 'log', 'log10', 'log1p', 'log2', 'logaddexp', 'logaddexp2', 'logical_and', 'logical_not', 'logical_or', 'logical_xor', 'logspace', 'longcomplex', 'longdouble', 'longfloat', 'longlong', 'lookfor', 'ma', 'mafromtxt', 'mask_indices', 'mat', 'math', 'matmul', 'matrix', 'matrixlib', 'max', 'maximum', 'maximum_sctype', 'may_share_memory', 'mean', 'median', 'memmap', 'meshgrid', 'mgrid', 'min', 'min_scalar_type', 'minimum', 'mintypecode', 'mod', 'modf', 'moveaxis', 'msort', 'multiply', 'nan', 'nan_to_num', 'nanargmax', 'nanargmin', 'nancumprod', 'nancumsum', 'nanmax', 'nanmean', 'nanmedian', 'nanmin', 'nanpercentile', 'nanprod', 'nanquantile', 'nanstd', 'nansum', 'nanvar', 'nbytes', 'ndarray', 'ndenumerate', 'ndfromtxt', 'ndim', 'ndindex', 'nditer', 'negative', 'nested_iters', 'newaxis', 'nextafter', 'nonzero', 'not_equal', 'numarray', 'number', 'obj2sctype', 'object0', 'object_', 'ogrid', 'oldnumeric', 'ones', 'ones_like', 'os', 'outer', 'packbits', 'pad', 'partition', 'percentile', 'pi', 'piecewise', 'place', 'poly', 'poly1d', 'polyadd', 'polyder', 'polydiv', 'polyfit', 'polyint', 'polymul', 'polynomial', 'polysub', 'polyval', 'positive', 'power', 'printoptions', 'prod', 'product', 'promote_types', 'ptp', 'put', 'put_along_axis', 'putmask', 'quantile', 'r_', 'rad2deg', 'radians', 'random', 'ravel', 'ravel_multi_index', 'real', 'real_if_close', 'rec', 'recarray', 'recfromcsv', 'recfromtxt', 'reciprocal', 'record', 'remainder', 'repeat', 'require', 'reshape', 'resize', 'result_type', 'right_shift', 'rint', 'roll', 'rollaxis', 'roots', 'rot90', 'round', 'round_', 'row_stack', 's_', 'safe_eval', 'save', 'savetxt', 'savez', 'savez_compressed', 'sctype2char', 'sctypeDict', 'sctypes', 'searchsorted', 'select', 'set_numeric_ops', 'set_printoptions', 'set_string_function', 'setbufsize', 'setdiff1d', 'seterr', 'seterrcall', 'seterrobj', 'setxor1d', 'shape', 'shares_memory', 'short', 'show_config', 'sign', 'signbit', 'signedinteger', 'sin', 'sinc', 'single', 'singlecomplex', 'sinh', 'size', 'sometrue', 'sort', 'sort_complex', 'source', 'spacing', 'split', 'sqrt', 'square', 'squeeze', 'stack', 'std', 'str0', 'str_', 'string_', 'subtract', 'sum', 'swapaxes', 'sys', 'take', 'take_along_axis', 'tan', 'tanh', 'tensordot', 'test', 'testing', 'tile', 'timedelta64', 'trace', 'tracemalloc_domain', 'transpose', 'trapezoid', 'tri', 'tril', 'tril_indices', 'tril_indices_from', 'trim_zeros', 'triu', 'triu_indices', 'triu_indices_from', 'true_divide', 'trunc',
        'typeDict', 'typecodes', 'typename', 'ubyte', 'ufunc', 'uint', 'uint0', 'uint16', 'uint32', 'uint64', 'uint8', 'uintc', 'uintp'
    ])

    #load_module_functions(numpy, np_funcs, exclusion_set)
    #exclusion_set.update(set(np_funcs))
    #load_module_functions(numpy.linalg, np_linalg_funcs, exclusion_set)

    """
    Operations on field output:
    +	Perform addition.
    −	Perform subtraction or unary negation.
    *	Perform multiplication.
    /	Perform division.
    abs(A)	Take the absolute value.
    acos(A)	Take the arccosine.
    asin(A)	Take the arcsine.
    atan(A)	Take the arctangent.
    cos(A)	Take the cosine.
    degreeToRadian(A)	Convert degrees to radians.
    exp(A)	Take the natural exponential.
    exp10(A)	Take the base 10 exponential.
    log(A)	Take the natural logarithm.
    log10(A)	Take the base 10 logarithm.
    power(FO,F)	Raise a field output object to a power.
    radianToDegree(A)	Convert radians to degrees.
    sin(A)	Take the sine.
    sqrt(A)	Take the square root.
    tan(A)	Take the tangent.
    """
    funcs_dict = {
        #'numpy' : ('Inf', 'Infinity', 'NaN', 'abs', 'add', 'all', 'allclose', 'angle', 'any', 'append', 'arange',
        #           'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan', 'arctan2', 'arctanh', 'argmax', 'argmin',
        #           'argsort', 'argwhere', 'around', 'array', 'array_equal', 'array_equiv', 'asarray', 'atleast_1d', 'atleast_2d',
        #           'atleast_3d', 'average', 'bartlett', 'base_repr', 'binary_repr', 'bincount', 'bitwise_and', 'bitwise_not', 'bitwise_or', 'bitwise_xor', 'blackman', 'block', 'bmat', 'bool8', 'broadcast', 'broadcast_arrays', 'broadcast_shapes', 'broadcast_to', 'busday_count', 'busday_offset', 'busdaycalendar', 'byte', 'byte_bounds', 'ceil', 'choose', 'clip', 'column_stack', 'common_type', 'compare_chararrays', 'compat', 'complex128', 'complex64', 'complexfloating', 'compress', 'concatenate', 'conj', 'conjugate', 'convolve', 'copy', 'copysign', 'copyto', 'core', 'corrcoef', 'correlate', 'cos', 'cosh', 'count_nonzero', 'cov', 'cross', 'csingle', 'ctypeslib', 'cumprod', 'cumproduct', 'cumsum', 'datetime64', 'datetime_as_string', 'datetime_data', 'deg2rad', 'degrees', 'delete', 'deprecate', 'deprecate_with_doc', 'diag', 'diag_indices', 'diag_indices_from', 'diagflat', 'diagonal', 'diff', 'digitize', 'disp', 'divide', 'divmod', 'dot', 'double', 'dsplit', 'dstack', 'dtype', 'e', 'ediff1d', 'einsum', 'einsum_path', 'emath', 'empty', 'empty_like', 'equal', 'errstate', 'euler_gamma',
        #           'exp', 'exp2', 'expand_dims', 'expm1', 'extract', 'eye', 'fabs', 'fastCopyAndTranspose', 'fft', 'fill_diagonal',
        # 'find_common_type', 'finfo', 'fix', 'flatiter', 'flatnonzero', 'flexible', 'flip', 'fliplr', 'flipud', 'float16', 'float32', 'float64', 'float_power', 'floating', 'floor', 'floor_divide', 'fmax', 'fmin', 'fmod', 'format_float_positional', 'format_float_scientific', 'format_parser', 'frexp', 'frombuffer', 'fromfile', 'fromfunction', 'fromiter', 'frompyfunc', 'fromregex', 'fromstring', 'full', 'full_like', 'gcd', 'generic', 'genfromtxt', 'geomspace', 'get_array_wrap', 'get_include', 'get_printoptions', 'getbufsize', 'geterr', 'geterrcall', 'geterrobj', 'gradient', 'greater', 'greater_equal', 'half', 'hamming', 'hanning', 'heaviside', 'histogram', 'histogram2d', 'histogram_bin_edges', 'histogramdd', 'hsplit', 'hstack', 'hypot', 'i0', 'identity', 'iinfo', 'imag', 'in1d', 'index_exp', 'indices', 'inexact', 'inf', 'info', 'infty', 'inner', 'insert', 'int0', 'int16', 'int32', 'int64', 'int8', 'intc', 'integer', 'interp', 'intersect1d', 'intp', 'invert', 'is_busday', 'isclose', 'iscomplex', 'iscomplexobj', 'isfinite', 'isfortran', 'isin', 'isinf', 'isnan', 'isnat', 'isneginf', 'isposinf', 'isreal', 'isrealobj', 'isscalar', 'issctype', 'issubdtype', 'issubsctype', 'iterable', 'kaiser', 'kron', 'lcm', 'ldexp', 'left_shift', 'less', 'less_equal', 'lexsort', 'lib', 'linalg', 'linspace', 'little_endian', 'load', 'loads', 'loadtxt', 'log', 'log10', 'log1p', 'log2', 'logaddexp', 'logaddexp2', 'logical_and', 'logical_not', 'logical_or', 'logical_xor', 'logspace', 'longcomplex', 'longdouble', 'longfloat', 'longlong', 'lookfor', 'ma', 'mafromtxt', 'mask_indices', 'mat', 'math', 'matmul', 'matrix', 'matrixlib', 'max', 'maximum', 'maximum_sctype', 'may_share_memory', 'mean', 'median', 'memmap', 'meshgrid', 'mgrid', 'min', 'min_scalar_type', 'minimum', 'mintypecode', 'mod', 'modf', 'moveaxis', 'msort', 'multiply', 'nan', 'nan_to_num', 'nanargmax', 'nanargmin', 'nancumprod', 'nancumsum', 'nanmax', 'nanmean', 'nanmedian', 'nanmin', 'nanpercentile', 'nanprod', 'nanquantile', 'nanstd', 'nansum', 'nanvar', 'nbytes', 'ndarray', 'ndenumerate', 'ndfromtxt', 'ndim', 'ndindex', 'nditer', 'negative', 'nested_iters', 'newaxis', 'nextafter', 'nonzero', 'not_equal', 'numarray', 'number', 'obj2sctype', 'object0', 'ogrid', 'oldnumeric', 'ones', 'ones_like', 'outer', 'packbits', 'pad', 'partition', 'percentile', 'pi', 'piecewise', 'place', 'poly', 'poly1d', 'polyadd', 'polyder', 'polydiv', 'polyfit', 'polyint', 'polymul', 'polynomial', 'polysub', 'polyval', 'positive', 'power', 'printoptions', 'prod', 'product', 'promote_types', 'ptp', 'put', 'put_along_axis', 'putmask', 'quantile', 'rad2deg', 'radians', 'random', 'ravel', 'ravel_multi_index', 'real', 'real_if_close', 'rec', 'recarray', 'recfromcsv', 'recfromtxt', 'reciprocal', 'record', 'remainder', 'repeat', 'require', 'reshape', 'resize', 'result_type', 'right_shift', 'rint', 'roll', 'rollaxis', 'roots', 'rot90', 'round', 'row_stack', 'savetxt', 'searchsorted', 'setdiff1d', 'setxor1d', 'shape', 'sign', 'sin', 'sinc', 'sinh', 'sometrue', 'sort', 'sqrt', 'square', 'squeeze', 'stack', 'std', 'sum', 'swapaxes', 'tan', 'tanh', 'tensordot', 'tile', 'transpose', 'trapezoid', 'union1d', 'unique', 'vdot', 'vectorize', 'vsplit', 'vstack', 'where', 'zeros'),
        'numpy' : ('abs', 'arccos', 'arcsin', 'arctan', 'cos', 'degrees', 'radians', 'exp', 'log', 'log10', 'power', 'sqrt', 'tan'),
        'scipy.special': ('exp10', ),
    }
    i = 0
    modules_form = []
    module_map = {}
    for mod in modules:
        mod_name, module, func_list, func_list2 = mod
        if mod_name in funcs_dict:
            func_list.extend(funcs_dict[mod_name])
        else:
            # load the functions into func_list
            continue
            load_module_functions(module, func_list, exclusion_set)
        exclusion_set.update(set(func_list))

        # build the lookup table
        modules_form.append([mod_name, None, func_list2])
        for func in func_list:
            module_map[i] = (module, func)
            func_list2.append((func, i, []))
            i += 1

    #for func in np_funcs:
        #func_list2.append((func, i, []))
        #if i > 3:
            #break
    #for func in np_linalg_funcs:
        #np_linalg_funcs2.append((func, i, []))
        #i += 1
    print(np_funcs)

    #modules_form = [
        #['Geometry', None, [
            #('NodeID', 0, []),
            #('ElementID', 1, []),
            #('PropertyID', 2, []),
            #('MaterialID', 3, []),
        #]]
    #]
    #modules_form = ['Modules', None, [
        #['numpy', None, [
            #[
                #('DataSource', 0, []),
                #('Datetime64', 1, []),
                #('Inf', 2, []),
                #('Infinity', 3, []),
            #]]
            #]]
    #]


    #[0, 1, 2, 3, 4, 5, 6, 7, 8]
    main_window = CalculatorWindow(data2,
                                   modules_form=(module_map, modules_form),
                                   fringe_cases=form)
    main_window.show()
    # Enter the main loop
    app.exec_()


if __name__ == "__main__": # pragma: no cover
    main()
