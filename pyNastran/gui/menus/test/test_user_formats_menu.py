"""tests user formats in the MainWindow and the User Formats menu (requires qtpy)"""
import os
import json
import tempfile
import unittest
from unittest import mock

import pyNastran
from pyNastran import DEV
from pyNastran.gui.menus.test.test_gui_menu import UsesQApplication
from pyNastran.gui.main_window import MainWindow
from pyNastran.gui.menus.user_formats.user_formats import UserFormatsWindow
from pyNastran.gui.user_formats import read_user_formats, write_user_formats
from pyNastran.gui.formats import CLASS_MAP, DEV_FORMATS
from pyNastran.gui.menus.preferences.preferences import PreferencesWindow
from pyNastran.gui.utils.qt.qsettings import QSettingsLike2

PKG_PATH = pyNastran.__path__[0]
CART3D_IO_PATH = os.path.join(PKG_PATH, 'converters', 'cart3d', 'cart3d_io.py')
TRI_FILENAME = os.path.join(PKG_PATH, 'converters', 'cart3d', 'models', 'threePlugs.a.tri')


class TestUserFormatsMenu(UsesQApplication):
    def test_user_formats_gui(self) -> None:
        """
        Tests (exact checks, no tolerances):
         - a user 'cart3d' (loaded from cart3d_io.py) overwrites CLASS_MAP['cart3d']
           in the gui without modifying the module-level CLASS_MAP
         - the user format loads a model
         - the menu adds a row (on_browse), saves the json, and reloads the formats
         - removing the row and applying restores the built-in class
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            json_filename = os.path.join(tmp_dir, 'pynastrangui_formats.json')
            write_user_formats({
                'cart3d': {'path': CART3D_IO_PATH, 'class': 'Cart3dIO', 'extensions': []},
            }, json_filename)
            inputs = {
                'debug': True, 'is_groups': True, 'is_gui': False, 'log': None,
                'geomscript': None, 'postscript': None, 'format': None,
                'user_points': None, 'user_geom': None,
                'user_formats_filename': json_filename,
            }
            gui = MainWindow(inputs)
            gui.is_gui = False
            built_in_cart3d = CLASS_MAP['cart3d']
            user_cart3d = gui.format_class_map['cart3d']
            self.assertIsNot(user_cart3d, built_in_cart3d)
            self.assertEqual(user_cart3d.__module__, 'pynastrangui_user_format_cart3d')
            self.assertIs(CLASS_MAP['cart3d'], built_in_cart3d)
            self.assertIn('cart3d', gui.supported_formats)
            self.assertIn('user_formats', gui.actions)

            gui.on_load_geometry(infile_name=TRI_FILENAME, geometry_format='cart3d',
                                 stop_on_failure=True)

            # menu: add a new format with a dotted module; name comes from the file
            data = {
                'font_size': 8,
                'filename': json_filename,
                'user_formats': read_user_formats(json_filename),
            }
            window = UserFormatsWindow(data, win_parent=gui)
            window.on_add()
            window.on_browse(path=os.path.join(tmp_dir, 'stl_io.py'))
            irow = window.table.rowCount() - 1
            self.assertEqual(window._text(irow, 0), 'stl')
            window.table.item(irow, 1).setText('pyNastran.converters.stl.stl_io')
            window.table.item(irow, 2).setText('STL_IO')
            window.table.item(irow, 3).setText('stl, .STLB')
            self.assertTrue(window.on_apply())
            self.assertEqual(read_user_formats(json_filename)['stl']['extensions'],
                             ['.stl', '.stlb'])
            self.assertIn('stl', gui.supported_formats)

            # partial row -> error, nothing written
            window.on_add()
            window.table.item(window.table.rowCount() - 1, 0).setText('partial')
            self.assertFalse(window.on_apply())
            self.assertNotIn('partial', read_user_formats(json_filename))

            # remove everything -> back to the built-in cart3d
            for irow in reversed(range(window.table.rowCount())):
                window.table.setCurrentCell(irow, 0)
                window.on_remove()
            self.assertEqual(window.table.rowCount(), 0)
            window.on_ok()
            self.assertEqual(read_user_formats(json_filename), {})
            self.assertIs(gui.format_class_map['cart3d'], built_in_cart3d)
            gui.on_set_user_formats()


    def test_format_preferences(self) -> None:
        """
        Tests the Preferences "Formats" tab (exact checks, no tolerances):
         - every format in all_fmt_order has a checkbox; nastran can't be disabled
         - DEV formats default to the DEV flag; others default to enabled
         - disabling stl removes it from supported_formats; only overrides are stored
         - Reset Defaults restores the defaults
         - format_enabled round trips through the json settings file
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            json_filename = os.path.join(tmp_dir, 'pynastrangui_formats.json')
            inputs = {
                'debug': True, 'is_groups': True, 'is_gui': False, 'log': None,
                'geomscript': None, 'postscript': None, 'format': None,
                'user_points': None, 'user_geom': None,
                'user_formats_filename': json_filename,
            }
            gui = MainWindow(inputs)
            gui.is_gui = False
            # don't overwrite ~/pyNastranGUI.json
            gui._save_settings = mock.Mock()
            gui.settings.format_enabled = {}
            gui.update_formats()

            for fmt in DEV_FORMATS:
                self.assertIn(fmt, gui.all_fmt_order)
                self.assertEqual(gui.is_format_enabled(fmt), DEV)
                self.assertEqual(fmt in gui.fmt_order, DEV)
            self.assertIn('stl', gui.supported_formats)

            data = gui.preferences_obj.get_preferences_data()
            window = PreferencesWindow(data, win_parent=gui)
            self.assertEqual(list(window.format_checkboxes), gui.all_fmt_order)
            self.assertFalse(window.format_checkboxes['nastran'].isEnabled())

            # disable stl; flip a dev format
            window.format_checkboxes['stl'].setChecked(False)
            window.format_checkboxes['obj'].setChecked(not DEV)
            self.assertTrue(window.on_apply())
            self.assertEqual(gui.settings.format_enabled, {'stl': False, 'obj': not DEV})
            self.assertNotIn('stl', gui.supported_formats)
            self.assertNotIn('stl', gui.fmt_order)
            self.assertEqual('obj' in gui.fmt_order, not DEV)
            gui._save_settings.assert_called()

            # json round trip
            qsettings = QSettingsLike2()
            qsettings._filename = os.path.join(tmp_dir, 'pyNastranGUI.json')
            qsettings.setValue('format_enabled', gui.settings.format_enabled)
            qsettings.save_json()
            with open(qsettings._filename, 'r') as json_file:
                self.assertEqual(json.load(json_file)['format_enabled'],
                                 {'stl': False, 'obj': not DEV})

            window.on_reset_defaults()
            self.assertTrue(window.on_apply())
            self.assertEqual(gui.settings.format_enabled, {})
            self.assertIn('stl', gui.supported_formats)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
