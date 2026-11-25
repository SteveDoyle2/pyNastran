"""tests pyNastran/gui/user_formats.py (no Qt required)"""
import os
import sys
import tempfile
import unittest
from unittest import mock

from cpylog import SimpleLogger
import pyNastran
from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
from pyNastran.gui import user_formats as user_formats_module
from pyNastran.gui.user_formats import (
    strip_json_comments, read_user_formats, write_user_formats,
    load_user_format_classes)
from pyNastran.gui.arg_handling import (
    determine_format, get_input_format_to_extension, _validate_format)

PKG_PATH = pyNastran.__path__[0]
CART3D_IO_PATH = os.path.join(PKG_PATH, 'converters', 'cart3d', 'cart3d_io.py')


class TestUserFormats(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.json_filename = os.path.join(self.tmp_dir.name, 'pynastrangui_formats.json')
        self.log = SimpleLogger(level='error')

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

    def _write(self, text: str) -> None:
        with open(self.json_filename, 'w', encoding='utf-8') as json_file:
            json_file.write(text)

    def test_strip_json_comments(self) -> None:
        """
        Tests '#' comment removal (exact string match):
         - full-line and trailing comments are removed; newlines are kept
         - '#' and escaped quotes inside strings are not comments
        """
        text = (
            '# header\n'
            '{"a#b": "c\\"#d", # trailing\n'
            ' "e": 1}\n')
        expected = (
            '\n'
            '{"a#b": "c\\"#d", \n'
            ' "e": 1}\n')
        self.assertEqual(strip_json_comments(text), expected)

    def test_read_write(self) -> None:
        """
        Tests the json reader/writer:
         - missing file -> {}
         - names lowercased, extensions get a leading '.' and are lowercased
         - write/read round trip is exact
         - bad json / missing path / missing class are skipped (not raised)
        """
        self.assertEqual(read_user_formats(self.json_filename, log=self.log), {})
        self._write('# only comments\n# {"cart3d": {}}\n')
        self.assertEqual(read_user_formats(self.json_filename, log=self.log), {})

        self._write(
            '# my formats\n'
            '{\n'
            '    "Cart3D": {"path": "C:/a/#1/cart3d_io.py", "class": "Cart3dIO",  # comment\n'
            '               "extensions": ["TRI", ".triq"]},\n'
            '    # "stl": {"path": "stl_io.py", "class": "STL_IO"},\n'
            '    "no_class": {"path": "foo.py"},\n'
            '    "no_path": {"class": "Foo"},\n'
            '    "not_a_dict": 1\n'
            '}\n')
        user_formats = read_user_formats(self.json_filename, log=self.log)
        expected = {
            'cart3d': {'path': 'C:/a/#1/cart3d_io.py', 'class': 'Cart3dIO',
                       'extensions': ['.tri', '.triq']},
        }
        self.assertEqual(user_formats, expected)

        write_user_formats(user_formats, self.json_filename)
        self.assertEqual(read_user_formats(self.json_filename, log=self.log), expected)

        # single backslashes are invalid json
        self._write('{"cart3d": {"path": "C:\\work\\cart3d_io.py", "class": "Cart3dIO"}}')
        self.assertEqual(read_user_formats(self.json_filename, log=self.log), {})
        self._write('[1, 2]')
        self.assertEqual(read_user_formats(self.json_filename, log=self.log), {})

    def test_load_classes(self) -> None:
        """
        Tests importing the format classes:
         - *.py path and dotted module both load Cart3dIO
         - bad path/module/class and a name without
           get_<name>_wildcard_geometry_results_functions are skipped
        """
        user_formats = {
            'cart3d': {'path': CART3D_IO_PATH, 'class': 'Cart3dIO', 'extensions': []},
            'bad_path': {'path': 'C:/missing/missing_io.py', 'class': 'Foo', 'extensions': []},
            'bad_module': {'path': 'pyNastran.missing_module', 'class': 'Foo', 'extensions': []},
            'bad_class': {'path': CART3D_IO_PATH, 'class': 'Missing', 'extensions': []},
            # Cart3dIO doesn't have get_tri_wildcard_geometry_results_functions
            'tri': {'path': CART3D_IO_PATH, 'class': 'Cart3dIO', 'extensions': []},
        }
        class_map = load_user_format_classes(user_formats, self.log)
        self.assertEqual(list(class_map), ['cart3d'])
        self.assertEqual(class_map['cart3d'].__name__, 'Cart3dIO')
        self.assertIn('pynastrangui_user_format_cart3d', sys.modules)

        user_formats = {
            'cart3d': {'path': 'pyNastran.converters.cart3d.cart3d_io',
                       'class': 'Cart3dIO', 'extensions': []},
        }
        class_map = load_user_format_classes(user_formats, self.log)
        self.assertIs(class_map['cart3d'], Cart3dIO)

    def test_arg_handling(self) -> None:
        """
        Tests the command line hooks (USER_FORMATS_FILENAME is patched):
         - a new format name is allowed for -f
         - user extensions map to the user format and are removed from the built-in
         - an unknown format still fails
        """
        user_formats = {
            'my_cart3d': {'path': CART3D_IO_PATH, 'class': 'Cart3dIO',
                          'extensions': ['.tri', '.mytri']},
        }
        write_user_formats(user_formats, self.json_filename)
        with mock.patch.object(user_formats_module, 'USER_FORMATS_FILENAME', self.json_filename):
            _validate_format(['my_cart3d', 'cart3d', 'nastran'])
            with self.assertRaises(AssertionError):
                _validate_format(['not_a_format'])

            format_to_extension = get_input_format_to_extension()
            self.assertEqual(format_to_extension['my_cart3d'], ['.tri', '.mytri'])
            self.assertEqual(format_to_extension['cart3d'], ['.triq'])
            self.assertEqual(determine_format('model.mytri'), 'my_cart3d')
            self.assertEqual(determine_format('model.tri'), 'my_cart3d')
            self.assertEqual(determine_format('model.triq'), 'cart3d')

        # no user file -> built-ins
        missing = os.path.join(self.tmp_dir.name, 'missing.json')
        with mock.patch.object(user_formats_module, 'USER_FORMATS_FILENAME', missing):
            self.assertEqual(determine_format('model.tri'), 'cart3d')


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
