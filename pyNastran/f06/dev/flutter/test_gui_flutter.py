from pathlib import Path
import unittest

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.f06.dev.flutter.action import Action
from pyNastran.f06.dev.flutter.utils import (
    validate_json, get_plot_flags,
    load_f06_op2, get_png_filename, get_plot_file,
    point_removal_str_to_point_removal,
    get_point_removal_str,
    PLOT_TYPES, X_PLOT_TYPES,)

IS_DEV = pyNastran.DEV
if IS_DEV:
    from pyNastran.f06.dev.flutter.nastran_utils import (
        get_element_table, get_property_table, get_material_table,
        get_table_trees, read_obj, write_obj)
    from pyNastran.dev.bdf_vectorized3.bdf import read_bdf

try:
    import pandas as pd
    IS_PANDAS = True
except ModuleNotFoundError:
    IS_PANDAS = False
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').absolute()
AERO_PATH = MODEL_PATH / 'aero'


class TestGuiFlutter(unittest.TestCase):
    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_split_by_pattern(self):
        from pyNastran.f06.dev.flutter.utils_report import split_by_pattern
        base_filenames = [
            'model_plane_mach_0.5_mgtow_kactuator_100',
            'model_plane_mach_0.2_mgtow_kactuator_50',
            'model_plane_mach_0.2_bdfw_kactuator_100',
            'model_plane_mach_0.2_bdfw_kactuator_50',
        ]
        expected = [
            ['model_plane_mach', '0.5', 'mgtow', 'kactuator', '100'],
            ['model_plane_mach', '0.2', 'mgtow', 'kactuator', '50'],
            ['model_plane_mach', '0.2', 'bdfw', 'kactuator', '100'],
            ['model_plane_mach', '0.2', 'bdfw', 'kactuator', '50']]
        split_paths = split_by_pattern(
            base_filenames, delimiter='_', group_common=True)
        assert split_paths == expected

    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_write_docx_path(self):
        from pyNastran.f06.dev.flutter.utils_report import write_docx_path
        mypath = 'cat/frog/pig/dog.txt'
        out = write_docx_path(mypath, ndir_levels=0)
        assert out == 'dog.txt', out

        out = write_docx_path(mypath, ndir_levels=1)
        assert out == 'pig/dog.txt', out

        out = write_docx_path(mypath, ndir_levels=2)
        assert out == 'frog/pig/dog.txt', out

    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_filenames_to_data_table_good(self):
        from pyNastran.f06.dev.flutter.utils_report import filenames_to_data_table, data_to_dataframe
        base_filenames = [
            'model_plane_mach_0.5_mgtow_kactuator_100.f06',
            'model_plane_mach_0.2_mgtow_kactuator_50.f06',
            'model_plane_mach_0.2_bdfw_kactuator_100.f06',
            'model_plane_mach_0.2_bdfw_kactuator_50.f06',
        ]
        headers, data_table = filenames_to_data_table(base_filenames)
        # print(f'headers = {headers}')
        # print(f'data_table = {data_table}')
        assert headers == ['1', '2', '3', '4', '5', 'File'], headers
        expected = [
            ['model_plane_mach', '0.5', 'mgtow', 'kactuator', '100', 'model_plane_mach_0.5_mgtow_kactuator_100.f06'],
            ['model_plane_mach', '0.2', 'mgtow', 'kactuator', '50', 'model_plane_mach_0.2_mgtow_kactuator_50.f06'],
            ['model_plane_mach', '0.2', 'bdfw', 'kactuator', '100', 'model_plane_mach_0.2_bdfw_kactuator_100.f06'],
            ['model_plane_mach', '0.2', 'bdfw', 'kactuator', '50', 'model_plane_mach_0.2_bdfw_kactuator_50.f06']]
        assert data_table == expected, data_table
        df = data_to_dataframe(
            headers, data_table,
            convert_numeric=True)
        print(df)

    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_filenames_to_data_table_good_pformat(self):
        from pyNastran.f06.dev.flutter.utils_report import (
            filenames_to_data_table, data_to_dataframe,
            get_icases)
        base_filenames = [
            'model_plane_mach_m0p5_mgtow_kactuator_100.f06',  # 3
            'model_plane_mach_0p5_mgtow_kactuator_50.f06',   # 0
            'model_plane_mach_0p2_bdfw_kactuator_100.f06',   # 2
            'model_plane_mach_0p2_bdfw_kactuator_50.f06',    # 1
        ]
        headers, data_table = filenames_to_data_table(base_filenames)
        # print(f'headers = {headers}')
        # print(f'data_table = {data_table}')
        assert headers == ['1', '2', '3', '4', '5', 'File'], headers
        expected = [
            ['model_plane_mach', 'm0p5', 'mgtow', 'kactuator', '100', 'model_plane_mach_m0p5_mgtow_kactuator_100.f06'],
            ['model_plane_mach', '0p5', 'mgtow', 'kactuator', '50', 'model_plane_mach_0p5_mgtow_kactuator_50.f06'],
            ['model_plane_mach', '0p2', 'bdfw', 'kactuator', '100', 'model_plane_mach_0p2_bdfw_kactuator_100.f06'],
            ['model_plane_mach', '0p2', 'bdfw', 'kactuator', '50', 'model_plane_mach_0p2_bdfw_kactuator_50.f06']]
        assert data_table == expected, data_table
        df = data_to_dataframe(headers, data_table, convert_numeric=True)
        mach = df['2'].tolist()
        kact = df['5'].tolist()
        # print(f'mach = {mach}')
        assert np.allclose(mach, [0.5, 0.5, 0.2, 0.2]), mach
        # print(f'fuel = {kact}')
        assert np.allclose(kact, [100, 50, 100, 50]), kact

        #                     mach   fuel            kact
        #                   1    2      3          4    5
        # 3  model_plane_mach  0.2   bdfw  kactuator   50
        # 2  model_plane_mach  0.2   bdfw  kactuator  100
        # 1  model_plane_mach  0.2  mgtow  kactuator   50
        # 0  model_plane_mach  0.5  mgtow  kactuator  100
        imach = '2'
        ifuel = '3'
        ikact = '5'
        # xaxis = imach
        cols = [ifuel, ikact]
        df = df.drop(columns=['File'])
        #    mach   fuel            kact
        #       2      3          4    5
        # 3   0.2   bdfw  kactuator   50
        # 2   0.2   bdfw  kactuator  100
        # 1   0.2  mgtow  kactuator   50
        # 0   0.5  mgtow  kactuator  100
        # [[3 'bdfw' 50]
        #  [2 'bdfw' 100]
        #  [1 'mgtow' 50]
        #  [0 'mgtow' 100]] 4
        icases, icase_dict = get_icases(df, cols)
        assert icases == [[3, 2], [1, 0]], icases

        cols = [ifuel, ikact, imach]
        icases, icase_dict = get_icases(df, cols)
        assert icases == [[1], [3], [0], [2]], icases
        # print(f'icases = {icases}')


    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_filenames_to_data_table_icases(self):
        from pyNastran.f06.dev.flutter.utils_report import (
            filenames_to_data_table, data_to_dataframe,
            get_icases, get_trades)
        base_filenames = [
            'model_plane_mach_0.2_mgtow_kactuator_100.f06',  # 1
            'model_plane_mach_0.2_mgtow_kactuator_50.f06',   # 0
            'model_plane_mach_0.2_bdfw_kactuator_100.f06',   # 3
            'model_plane_mach_0.2_bdfw_kactuator_50.f06',    # 2

            'model_plane_mach_0.5_mgtow_kactuator_100.f06',  # 5
            'model_plane_mach_0.5_mgtow_kactuator_50.f06',   # 4
            'model_plane_mach_0.5_bdfw_kactuator_100.f06',   # 7
            'model_plane_mach_0.5_bdfw_kactuator_50.f06',    # 6
        ]
        headers, data_table = filenames_to_data_table(base_filenames)
        df = data_to_dataframe(headers, data_table, convert_numeric=True)
        imach = '2'
        ifuel = '3'
        ikact = '5'
        cols = [ifuel, ikact, imach]
        # mach = df['2'].tolist()
        # fuel = df['5'].tolist()
        icases, icases_dict = get_icases(df, cols)
        assert icases == [[3, 7], [1, 5], [2, 6], [0, 4]], icases

        log = SimpleLogger(level='debug')
        trade_str = (
            # f'{ikact}, {ifuel}, {imach}; '
            f'{ifuel}, {ikact}, {imach}; '
        )
        is_passed, configs, trades = get_trades(df, trade_str, log)
        assert is_passed, is_passed

        configs_expected = [
            '3=mgtow, 5=100', '3=mgtow, 5=50',
            '3=bdfw, 5=100', '3=bdfw, 5=50',
            '3=mgtow, 5=100', '3=mgtow, 5=50',
            '3=bdfw, 5=100', '3=bdfw, 5=50']
        trades_expected = [
            (
                ['3', '5', '2'],
                {('bdfw', 50): [3, 7],
                 ('mgtow', 50): [1, 5],
                 ('bdfw', 100): [2, 6],
                 ('mgtow', 100): [0, 4]},
            )
        ]
        assert configs == configs_expected
        assert trades == trades_expected
        # print(f'configs = {configs}')
        # print(f'trades = {trades}')

    @unittest.skipIf(not IS_PANDAS, 'pandas is needed')
    def test_filenames_to_data_table_bad(self):
        """bad split, but good enough"""
        from pyNastran.f06.dev.flutter.utils_report import filenames_to_data_table
        base_filenames = [
            'model_plane_mach_0.5_mgtow_kactuator_100.f06',
            'model_plane_mach_0.2_mgtow_kactuator_50.f06',
            'model_plane_mach_0.2_bdfw_kactuator_100.f06',
            'model_plane_mach_0.2_bdfw_kactuator_50.f06',
            'cat.f06',
        ]
        expected = [
            ['model', 'plane', 'mach', '0.5', 'mgtow', 'kactuator', '100', 'model_plane_mach_0.5_mgtow_kactuator_100.f06'],
            ['model', 'plane', 'mach', '0.2', 'mgtow', 'kactuator', '50', 'model_plane_mach_0.2_mgtow_kactuator_50.f06'],
            ['model', 'plane', 'mach', '0.2', 'bdfw', 'kactuator', '100', 'model_plane_mach_0.2_bdfw_kactuator_100.f06'],
            ['model', 'plane', 'mach', '0.2', 'bdfw', 'kactuator', '50', 'model_plane_mach_0.2_bdfw_kactuator_50.f06'],
            ['cat', '', '', '', '', '', '', 'cat.f06'],]
        headers, data_table = filenames_to_data_table(base_filenames)
        assert headers == ['1', '2', '3', '4', '5', '6', '7', 'File'], headers
        assert data_table == expected, data_table

    @unittest.skipIf(not IS_DEV, 'no flutter-dev')
    def test_flutter_nastran_utils(self) -> None:
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        model = read_bdf(bdf_filename, debug='warning')
        get_element_table(model)
        get_property_table(model)
        get_material_table(model)
        get_table_trees(model, model)

    @unittest.skipIf(not IS_DEV, 'no flutter-dev')
    def test_flutter_static_solid_shell_bar(self) -> None:
        bdf_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.bdf'
        obj_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.obj'
        model = read_bdf(bdf_filename, debug='warning')
        get_element_table(model)
        get_property_table(model)
        get_material_table(model)
        get_table_trees(model, model)
        write_obj(model, obj_filename)
        model2 = read_obj(obj_filename)
        assert model == model2
        obj_filename.unlink()

    def test_flutter_action(self) -> None:
        act = Action('cat', 'dog', show=True)
        str(act)

    def test_flutter_load_f06_op2(self) -> None:
        f06_filename = MODEL_PATH / 'aero' / '2_mode_flutter' / '0012_flutter.op2'
        log = SimpleLogger(level='warning')
        in_units = 'si'
        out_units = 'si'
        load_f06_op2(f06_filename, log, in_units, out_units, use_rhoref=False)

        op2_filename = MODEL_PATH / 'bwb' / 'bwb_saero.op2'
        load_f06_op2(op2_filename, log, in_units, out_units, use_rhoref=False)

    @unittest.skipIf(not IS_DEV, 'no flutter-dev')
    def test_plot_flutter_0012_dev(self):
        """tests load_f06_op2 and get_png_filename"""
        log = SimpleLogger(level='warning')
        dirname = AERO_PATH / '2_mode_flutter'
        f06_filename = dirname / '0012_flutter.f06'
        model, responses = load_f06_op2(f06_filename, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 1, responses

        f06_filename2 = 'cat.f06'
        model, responses = load_f06_op2(
            f06_filename2, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 0, responses

        f06_filename2 = 'cat.op2'
        model, responses = load_f06_op2(
            f06_filename2, log,
            in_units='si', out_units='si', use_rhoref=False)
        assert model is None, model
        assert len(responses) == 0, responses

    def test_plot_file(self):
        out = get_plot_file()

        base = 'base'
        x_plot_type = 'eas'
        plot_type = 'Vg'
        export_to_png = True
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, plot_type,
            export_to_png)
        assert png_filename0 == 'base_Vg.png', png_filename0
        assert png_filename == 'base_Vg.png', png_filename
        png_filename0, png_filename = get_png_filename(
            base, x_plot_type, 'x-'+plot_type,
            export_to_png)
        assert png_filename0 == 'base_eas-Vg.png', png_filename0
        assert png_filename == 'base_eas-Vg.png', png_filename

        png_filename0, png_filename = get_png_filename(
            base, 'x-'+x_plot_type, plot_type,
            export_to_png=False)
        assert png_filename0 == 'base_Vg.png', png_filename0
        assert png_filename is None, png_filename

    def test_get_plot_flags(self):
        # X_PLOT_TYPES = ['eas', 'tas', 'rho', 'q', 'mach', 'alt', 'kfreq', 'ikfreq', 'index']
        # PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus', 'modal-participation', 'zimmerman']
        flags = get_plot_flags('x-damp-freq', 'eas')
        for plot_type in PLOT_TYPES:
            for x_plot_type in X_PLOT_TYPES:
                flags = get_plot_flags(plot_type, x_plot_type)
        del flags
        with self.assertRaises(AssertionError):
            get_plot_flags('x-damp-freq???', 'eas')

    def test_flutter_point_removal(self):
        log = SimpleLogger(level='debug')
        msg = '400:410,450:500'
        point_removal = [[400.0, 410.0], [450.0, 500.0]]
        point_removal_list = point_removal_str_to_point_removal(msg, log)
        assert np.allclose(point_removal, point_removal_list), point_removal_list
        out = get_point_removal_str(point_removal_list)
        assert out == msg, f'out={out!r} expected={msg!r}'
        #-----------
        point_removal = [[450.0, -1.0]]
        msg = '450:'
        point_removal_list = point_removal_str_to_point_removal(msg, log)
        assert np.allclose(point_removal, point_removal_list), point_removal_list
        out = get_point_removal_str(point_removal_list)
        assert out == msg, f'out={out!r} expected={msg!r}'
        #---
        point_removal = [[-1.0, 500.0]]
        msg = ':500'
        point_removal_list = point_removal_str_to_point_removal(msg, log)
        assert np.allclose(point_removal, point_removal_list), point_removal_list
        out = get_point_removal_str(point_removal_list)
        assert out == msg, f'out={out!r} expected={msg!r}'
        #---
        # point_removal = []
        msg = ''
        point_removal_list = point_removal_str_to_point_removal(msg, log)
        assert point_removal_list == [], point_removal_list
        out = get_point_removal_str(point_removal_list)
        assert out == msg, f'out={out!r} expected={msg!r}'
        #---
        point_removal = []
        msg = ''
        point_removal_list = point_removal_str_to_point_removal(msg, log)
        assert point_removal_list == point_removal, point_removal_list
        out = get_point_removal_str(point_removal_list)
        assert out == msg, f'out={out!r} expected={msg!r}'

    def test_flutter_validate_json(self):
        log = SimpleLogger(level='debug')
        mydict = {}
        is_valid = validate_json(mydict, log)
        assert not is_valid, is_valid

        mydict = {
            'units_in': 'fakeenglish_in',
            'units_out': 'fakesi',
            'plot_type': 'fakex-damp-freq',
        }
        is_valid = validate_json(mydict, log)
        assert not is_valid, is_valid

        mydict = {
            'units_in': 'english_in',
            'units_out': 'si',
            'plot_type': 'x-damp-freq',
        }
        is_valid = validate_json(mydict, log)
        assert is_valid, is_valid


if __name__ == '__main__':
    unittest.main()
