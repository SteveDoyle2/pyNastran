from pathlib import Path
import unittest

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.op2.op2 import read_op2
from pyNastran.op2.tools.solution_combination import (
    run_load_case_combinations, run_load_case_multi_combinations)

np.seterr(all='raise')
DIRNAME = Path(__file__).parent
PKG_PATH = Path(pyNastran.__path__[0]).parent
MODEL_PATH = PKG_PATH / 'models'


class TestCombination(unittest.TestCase):
    def test_op2_combination_solid_bending(self):
        log = SimpleLogger(level='warning')
        lines = """
        # this is line 0. line 1 has the input subcase ids
                 ,     ,       1,
        # Subcase, Name,      Scale1,
        10,        case10,    1.0,
        20,        case20,    1.2
        30,        "case 30", 1.2,
        """
        op2_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.op2'
        combination_filename = DIRNAME / 'combination_file.txt'

        assert op2_filename.exists(), print_bad_path(op2_filename)
        with open(combination_filename, 'w') as combination_file:
            combination_file.write(lines)

        # combinations = load_combinations([combination_filename], delimiter=',')

        run_load_case_combinations(
            op2_filename, combination_filename, log=log)
        combination_filename.unlink()

    def test_op2_combination_solid_shell_bar(self):
        log = SimpleLogger(level='warning')
        lines = """
        # this is line 0. line 1 has the input subcase ids
        Null     , Subcases,  1,    # subcases
        # Subcase, Name,      Scale1,
        10,        case10,    1.0,  # some comment
        20,        case20,    1.2
        30,        "case 30", 1.2,
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        combination_filename = DIRNAME / 'combination_file.txt'

        assert op2_filename.exists(), print_bad_path(op2_filename)
        with open(combination_filename, 'w') as combination_file:
            combination_file.write(lines)

        # combinations = load_combinations([combination_filename], delimiter=',')

        run_load_case_combinations(
            op2_filename, combination_filename, log=log)
        combination_filename.unlink()

    def test_op2_multi_combination_solid_shell_bar(self):
        log = SimpleLogger(level='warning')
        lines = """
        # this is line 0. line 1 has the input subcase ids
                 , iOP2s,     0,  0 # op2
                 , Subcases,  1,  1 # subcases
        # Subcase, Name,      Scales,
        10,        case10,    1.0, 2.0,
        20,        case20,    1.2, 1.,
        30,        "case 30", 1.2, 0.
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        combination_filename = DIRNAME / 'combination_file.txt'

        assert op2_filename.exists(), print_bad_path(op2_filename)
        with open(combination_filename, 'w') as combination_file:
            combination_file.write(lines)

        model = read_op2(op2_filename, log=log)
        run_load_case_multi_combinations(
            {0: model}, combination_filename, log=log,
            mode='nx')

        run_load_case_multi_combinations(
            [op2_filename], combination_filename, log=log,
            mode='nx', include_results=['displacements'],
        )
        run_load_case_multi_combinations(
            [op2_filename, model], combination_filename, log=log,
            mode='nx', include_results=['displacements'],
        )

        run_load_case_multi_combinations(
            {0: op2_filename}, combination_filename, log=log,
            mode='nx', include_results=['displacements'])
        combination_filename.unlink()

    def test_op2_multi_combination_solid_shell_bar_string(self):
        log = SimpleLogger(level='warning')
        lines = """
        # this is line 0. line 1 has the input subcase ids
                 , iOP2s,     bending,  bending # op2
                 , Subcases,  1,  1 # subcases
        # Subcase, Name,      Scales,
        10,        case10,    1.0, 2.0,
        20,        case20,    1.2, 1.,
        30,        "case 30", 1.2, 0.
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        combination_filename = DIRNAME / 'combination_file.txt'

        assert op2_filename.exists(), print_bad_path(op2_filename)
        with open(combination_filename, 'w') as combination_file:
            combination_file.write(lines)

        model = read_op2(op2_filename, log=log)
        run_load_case_multi_combinations(
            {'bending': model}, combination_filename, log=log,
            mode='nx')

        run_load_case_multi_combinations(
            {'bending': op2_filename}, combination_filename, log=log,
            mode='nx', include_results=['displacements'])
        combination_filename.unlink()


if __name__ == '__main__':
    unittest.main()
