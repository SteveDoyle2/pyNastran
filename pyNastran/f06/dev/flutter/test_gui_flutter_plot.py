"""
Headless tests for the flutter plot GUI (gui_flutter_plot.py / plot_layout.py)

Every GUI interaction runs inside a real Qt event loop (like a user click),
so regressions such as a nested plt.show() event loop are caught instead of
hanging.  The settings json and output pngs go to a temporary directory,
so ~/plot_flutter.json is never touched.

requires qtpy, a Qt binding, and vtk (pyNastran.gui imports vtkmodules)
"""
import os
import sys
import shutil
import inspect
import tempfile
import unittest
from pathlib import Path

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import pyNastran
from pyNastran.f06.dev.flutter.utils import get_plot_file, PLOT_TYPES, X_PLOT_TYPES

# no windows pop up; must be set before the QApplication is created
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
try:
    from qtpy.QtCore import QTimer, qInstallMessageHandler
    from qtpy.QtWidgets import QApplication
    import vtkmodules  # pylint: disable=unused-import
    from pyNastran.f06.dev.flutter.gui.gui_flutter_plot import FlutterGui
    IS_GUI = True
except (ImportError, RuntimeError):  # qtpy raises a RuntimeError subclass
    IS_GUI = False

PKG_PATH = Path(pyNastran.__path__[0])
AERO_PATH = (PKG_PATH / '..' / 'models' / 'aero').resolve()
F06_0012 = AERO_PATH / '2_mode_flutter' / '0012_flutter.f06'  # PKNL; 1 subcase; 2 modes
F06_BAH = AERO_PATH / 'bah_plane' / 'bah_plane.f06'  # PK; 2 subcases; 10 modes

# the scalar settings that the plot code requires to be floats
FLOAT_SETTINGS = [
    'freq_tol', 'freq_tol_remove', 'mag_tol', 'vl', 'vf',
    'damping_required', 'damping_required_tol', 'damping',
]


@unittest.skipIf(not IS_GUI, 'qtpy, a Qt binding, and vtk are required')
class TestFlutterPlotGui(unittest.TestCase):
    """drives FlutterGui like a user: Load F06, select modes, Run"""
    @classmethod
    def setUpClass(cls) -> None:
        app = QApplication.instance()
        if app is None:
            app = QApplication([])
        cls.app = app

    def setUp(self) -> None:
        self.cwd = os.getcwd()
        self.backend = matplotlib.get_backend()
        # the GUI's real backend; required for the nested event loop test
        plt.switch_backend('QtAgg')

        self._tempdir = tempfile.TemporaryDirectory()
        self.dirname = Path(self._tempdir.name)
        self.json_filename = self.dirname / 'plot_flutter.json'
        self.f06_0012 = self._copy(F06_0012)
        self.f06_bah = self._copy(F06_BAH)

        self.qt_messages = []
        self.guis = []

    def tearDown(self) -> None:
        os.chdir(self.cwd)  # plot() chdirs to the f06 directory
        for gui in self.guis:
            gui.close()
            gui.deleteLater()
        plt.close('all')
        plt.switch_backend(self.backend)
        self._tempdir.cleanup()

    #---------------------------------------------------------------------------
    # helpers
    def _copy(self, f06_filename: Path) -> Path:
        """plot() writes the png next to the f06, so work on a copy"""
        out_filename = self.dirname / f06_filename.name
        shutil.copyfile(f06_filename, out_filename)
        return out_filename

    def _make_gui(self, save_filename=None) -> 'FlutterGui':
        if save_filename is None:
            save_filename = self.json_filename
        gui = FlutterGui(save_filename=str(save_filename))
        self.guis.append(gui)
        return gui

    def _run_in_event_loop(self, func):
        """
        Runs func inside a running Qt event loop (as a button click would be)
        and re-raises any exception, including ones raised inside Qt slots.
        Qt warnings are collected in self.qt_messages.
        """
        app = self.app
        out = {}
        slot_errors = []

        def excepthook(exc_type, exc_value, exc_tb):
            # without a custom hook, PyQt aborts on an exception in a slot
            slot_errors.append(exc_value)

        def message_handler(mode, context, message):
            self.qt_messages.append(str(message))

        def step():
            try:
                out['value'] = func()
            except BaseException as error:  # re-raised below
                out['error'] = error
            finally:
                app.quit()

        watchdog = QTimer()
        watchdog.setSingleShot(True)
        watchdog.timeout.connect(app.quit)

        old_excepthook = sys.excepthook
        sys.excepthook = excepthook
        old_handler = qInstallMessageHandler(message_handler)
        try:
            watchdog.start(120_000)
            QTimer.singleShot(0, step)
            app.exec_()
        finally:
            watchdog.stop()
            qInstallMessageHandler(old_handler)
            sys.excepthook = old_excepthook

        if slot_errors:
            raise slot_errors[0]
        if 'error' in out:
            raise out['error']
        if 'value' not in out:
            self.fail('the event loop timed out')
        return out['value']

    @staticmethod
    def _load(plot_layout, f06_filename: Path, plot_type: str) -> None:
        plot_layout.plot_type_pulldown.setCurrentText(plot_type)
        assert plot_layout.plot_type_pulldown.currentText() == plot_type, plot_type
        plot_layout.f06_filename_edit.setText(str(f06_filename))
        plot_layout.f06_load_button.click()

    @staticmethod
    def _select_rows(plot_layout, rows: list[int]) -> None:
        modes_widget = plot_layout.modes_widget
        for irow in range(modes_widget.count()):
            modes_widget.item(irow).setSelected(irow in rows)

    @staticmethod
    def _legend_modes() -> list[str]:
        """['Mode 5', 'Mode 7'] from the figure 1 legend"""
        axes = plt.figure(1).axes[0]
        labels = [text.get_text() for text in axes.get_legend().get_texts()]
        return [label.split(';')[0] for label in labels if label.startswith('Mode ')]

    #---------------------------------------------------------------------------
    def test_save_filename_default(self) -> None:
        """the default settings file is still ~/plot_flutter.json"""
        default = inspect.signature(FlutterGui.__init__).parameters['save_filename'].default
        assert default == get_plot_file(), f'default={default!r} expected={get_plot_file()!r}'

    def test_fresh_settings_load_run(self) -> None:
        """no saved json: Load F06 + Run works for each plot type"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]
        assert not self.json_filename.exists(), self.json_filename

        plot_types = ['x-damp-freq', 'x-damp-kfreq', 'root-locus', 'zimmerman']
        for plot_type in plot_types:
            with self.subTest(plot_type=plot_type):
                pngs_before = set(self.dirname.glob('*.png'))

                def load_run():
                    self._load(plot_layout, self.f06_0012, plot_type)
                    plot_layout.run_button.click()
                self._run_in_event_loop(load_run)

                fig = plt.figure(1)
                nlines = sum(len(axes.lines) for axes in fig.axes)
                assert nlines > 0, f'{plot_type}: nothing was plotted'
                assert plt.get_fignums() == [1], f'{plot_type}: extra figures={plt.get_fignums()}'
                window = fig.canvas.manager.window
                assert window.isVisible(), f'{plot_type}: the figure was not shown'
                # plot_zimmerman has no png_filename, so it doesn't export
                npngs_expected = 0 if plot_type == 'zimmerman' else 1
                new_pngs = set(self.dirname.glob('*.png')) - pngs_before
                assert len(new_pngs) == npngs_expected, f'{plot_type}: new pngs={new_pngs}'
        assert self.json_filename.exists(), 'settings were not saved to save_filename'

    def test_get_xlim_float_settings(self) -> None:
        """hidden/empty boxes still give float settings for every plot type"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]
        for plot_type in PLOT_TYPES:
            for x_plot_type in X_PLOT_TYPES:
                with self.subTest(plot_type=plot_type, x_plot_type=x_plot_type):
                    def set_validate():
                        plot_layout.plot_type_pulldown.setCurrentText(plot_type)
                        plot_layout.x_plot_type_pulldown.setCurrentText(x_plot_type)
                        plot_layout.validate()
                    self._run_in_event_loop(set_validate)
                    for name in FLOAT_SETTINGS:
                        value = getattr(gui, name)
                        assert isinstance(value, float), f'{name}={value!r}'

    def test_root_locus_selected_modes(self) -> None:
        """root locus only plots the highlighted modes"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]

        def load_select_run():
            self._load(plot_layout, self.f06_bah, 'root-locus')
            self._select_rows(plot_layout, [4, 6])
            plot_layout.run_button.click()
        self._run_in_event_loop(load_select_run)

        modes = self._legend_modes()
        assert modes == ['Mode 5', 'Mode 7'], modes
        ncollections = len(plt.figure(1).axes[0].collections)
        assert ncollections == 2, f'expected 1 scatter per mode; ncollections={ncollections}'

    def test_no_nested_event_loop(self) -> None:
        """Run shows the figure without starting a 2nd Qt event loop"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]
        for plot_type in ['x-damp-freq', 'root-locus']:
            with self.subTest(plot_type=plot_type):
                def load_run():
                    self._load(plot_layout, self.f06_0012, plot_type)
                    plot_layout.run_button.click()
                self._run_in_event_loop(load_run)

                nested = [msg for msg in self.qt_messages if 'already running' in msg]
                assert len(nested) == 0, nested
                window = plt.figure(1).canvas.manager.window
                assert window.isVisible(), 'the figure was not shown'

    def test_bad_f06_filename(self) -> None:
        """a missing f06 turns the box red and loads nothing"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]
        missing_filename = self.dirname / 'missing.f06'

        self._run_in_event_loop(
            lambda: self._load(plot_layout, missing_filename, 'x-damp-freq'))
        assert 'red' in plot_layout.f06_filename_edit.styleSheet()
        assert plot_layout.subcase_edit.count() == 0, plot_layout.subcase_edit.count()
        assert len(plot_layout.responses) == 0, plot_layout.responses

    def test_subcase_switch(self) -> None:
        """changing the subcase updates the modes table and the plot"""
        gui = self._make_gui()
        plot_layout = gui.plot_layout[0]

        def load_switch_run():
            self._load(plot_layout, self.f06_bah, 'root-locus')
            plot_layout.subcase_edit.setCurrentIndex(1)
            plot_layout.run_button.click()
        self._run_in_event_loop(load_switch_run)

        assert plot_layout.subcase_edit.count() == 2, plot_layout.subcase_edit.count()
        assert plot_layout.subcase_edit.currentText() == 'Subcase 2'
        nmodes = len(plot_layout.responses[2].modes)
        assert plot_layout.modes_widget.count() == nmodes, (plot_layout.modes_widget.count(), nmodes)
        title = plt.figure(1)._suptitle.get_text()
        assert title.startswith('Subcase 2'), title

    def test_settings_round_trip(self) -> None:
        """settings saved by one GUI are restored by the next one"""
        gui1 = self._make_gui()
        plot_layout1 = gui1.plot_layout[0]

        def set_load_run():
            plot_layout1.freq_tol_edit.setText('0.5')
            self._load(plot_layout1, self.f06_bah, 'root-locus')
            plot_layout1.run_button.click()
        self._run_in_event_loop(set_load_run)
        assert self.json_filename.exists(), self.json_filename

        gui2 = self._make_gui()
        plot_layout2 = gui2.plot_layout[0]
        assert plot_layout2.plot_type_pulldown.currentText() == 'root-locus'
        assert float(plot_layout2.freq_tol_edit.text()) == 0.5, plot_layout2.freq_tol_edit.text()
        assert gui2.recent_files[0] == os.path.abspath(self.f06_bah), gui2.recent_files
        assert Path(plot_layout2.f06_filename_edit.text()) == self.f06_bah, plot_layout2.f06_filename_edit.text()

    def test_no_save_filename(self) -> None:
        """save_filename='' runs without loading/saving any settings"""
        gui = self._make_gui(save_filename='')
        plot_layout = gui.plot_layout[0]

        def load_run():
            self._load(plot_layout, self.f06_0012, 'x-damp-freq')
            plot_layout.run_button.click()
        self._run_in_event_loop(load_run)

        assert gui.save_filename == '', gui.save_filename
        jsons = list(self.dirname.glob('*.json'))
        assert len(jsons) == 0, jsons
        assert np.any([len(axes.lines) for axes in plt.figure(1).axes])


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
