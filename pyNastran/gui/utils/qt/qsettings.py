# coding: utf-8
# pylint: disable=W0201,C0301
import os
import json
import warnings
import traceback

#from pyNastran.gui.qt_version import qt_int, qt_version
import numpy as np
from qtpy import QtCore
from pyNastran.utils.numpy_utils import integer_types, float_types


class QSettingsLike2:
    _tuples = {
        'background_color', 'background_color2',
        'highlight_color', 'shear_moment_torque_color',
        'corner_text_color', 'annotation_color',
        'screen_shape', 'screen_position',
        'nastran_caero_color', 'nastran_rbe_line_color',
        'nastran_plotel_color',
        'cart3d_fluent_include', 'cart3d_fluent_remove',
        'xyz_ref',
        'units_model_in',
    }
    def __init__(self):
        """
        json gui loading location is stored in:
        1) local directory
        2) exe directory
        3) home directory
        The priority goes in that order.

        The last directory would be stored as an additional preference
        of launch directory as:
        a) previous working directory
        b) use working directory

        """
        self.data = {}

        home_dirname = os.path.expanduser('~')
        # local_dirname = os.getcwd()
        # exe_dirname = os.path.dirname(__file__)
        # local_filename = os.path.join(local_dirname, 'pyNastranGUI.json')
        # exe_filename   = os.path.join(exe_dirname,   'pyNastranGUI.json')
        home_filename  = os.path.join(home_dirname,  'pyNastranGUI.json')
        # is_pynastrangui_exe = False
        # if os.path.exists(local_filename):
        #     self._filename = local_filename
        # elif os.path.exists(exe_filename) and is_pynastrangui_exe:
        #     self._filename = exe_filename
        # if os.path.exists(home_filename):
        self._filename = home_filename
        # else:
        #     self._filename = ''
    def clear(self) -> None:
        self.data = {}
    def childKeys(self) -> list[str]:
        return list(self.data.keys())

    def value(self, key: str, default=None):
        assert isinstance(key, str), key
        if key in self.data:
            value = self.data[key]
        elif default is None:
            raise RuntimeError('default=None?')
        else:
            value = default

        if key in {'main_window_geometry', 'main_window_state'}:
            value_bytes = value.encode('ascii')
            value = QtCore.QByteArray(value_bytes)

        if key in self._tuples:
            value = tuple(value)
        return value

    def setValue(self, key: str, value) -> None:
        assert isinstance(key, str), key
        self.data[key] = value

    def load_json(self) -> None:
        if os.path.exists(self._filename):
            with open(self._filename, 'r') as json_file:
                try:
                    self.data = json.load(json_file)
                except Exception as e:
                    warnings.warn(f'failed to load {self._filename}\n{str(e)}')
                    print(traceback.format_exc())
                    # print(traceback.format_exception_only(e))
                    # raise
                    return
        #x = 1

    def save_json(self) -> None:
        data = {}
        for key, value in self.data.items():
            if isinstance(value, str):
                value_out = value
            elif isinstance(value, float_types):
                # json struggles with numpy floats
                value_out = float(value)
            elif isinstance(value, integer_types):
                # json probably struggles with numpy ints
                value_out = int(value)
            elif key in self._tuples:
                value_out = totuple(key, value)
                assert isinstance(value_out, tuple), (key, value_out)
            elif key == 'recent_files':
                value_out = value
            elif value.__class__.__name__ == 'QByteArray':
                value2 = bytes(value.toBase64())
                value_out = value2.decode('ascii')
            else:  # pragma: no cover
                print(f'error...{key!r}={value} type={str(type(value))}')
                raise NotImplementedError(f'key={key!r} value={value!r}')
            data[key] = value_out

        with open(self._filename, 'w') as json_file:
            json.dump(data, json_file, indent=True)


def totuple(key: str, value) -> tuple:
    """json requires lists and arrays be tuples"""
    if isinstance(value, tuple):
        return value

    if isinstance(value, np.ndarray):
        value = tuple(value.tolist())
    elif isinstance(value, list):
        value = tuple(value)
    else:  # pragma: no cover
        print(f'error...{key!r}={value} type={str(type(value))}')
        raise NotImplementedError(f'key={key!r} value={value!r}')
    return value

#x = 1
#QSettingsLike = QtCore.QSettings
QSettingsLike = QSettingsLike2
