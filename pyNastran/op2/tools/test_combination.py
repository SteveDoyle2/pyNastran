from pathlib import Path
import unittest

from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.op2.tools.solution_combination import run_load_case_combinations
DIRNAME = Path(__file__).parent
PKG_PATH = Path(pyNastran.__path__[0]).parent
MODEL_PATH = PKG_PATH / 'models'

class TestCombintion(unittest.TestCase):
    def test_op2_combination_solid_bending(self):
        log = SimpleLogger(level='warning')
        lines = """
        # this is line 0. line 1 has the input subcase ids
                            1,
        # Subcase Name,      Scale1,
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
                            1,
        # Subcase Name,      Scale1,
        10,        case10,    1.0,
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


if __name__ == '__main__':
    unittest.main()
