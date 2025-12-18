import os
import copy
from typing import Optional
from cpylog import SimpleLogger
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.op2.op2 import read_op2, OP2


def load_combinations(combination_filenames: list[PathLike],
                      delimiter=',') -> list:
    """

    Combination File
    ----------------
    # line 1 is the subcase ids  for the scale factors
    # line 2 is a comment

    # this is line 0. line 1 has the input subcase ids
                        1      2      3
    # Subcase Name      Scale1 Scale2 Scale3
    10        case10    1.0    2.0    3.0
    20        case20    1.2    2.2    3.2
    30        "case 30" 1.2    2.2    3.2
    """
    assert len(combination_filenames), combination_filenames
    for combination_filename in combination_filenames:
        assert os.path.exists(combination_filename), print_bad_path(combination_filename)

    all_combinations = []
    for combination_filename in combination_filenames:
        subcases_input, combinations = _load_combination(
            combination_filename, delimiter=delimiter)
        all_combinations.append((subcases_input, combinations))
    assert len(all_combinations)
    return all_combinations


def _load_combination(combination_filename: PathLike,
                      delimiter=',') -> tuple[list[int], list]:
    """parses a single load case combination file"""
    with open(combination_filename, 'r') as combination_file:
        lines = combination_file.readlines()

    # remove comments and blank lines
    strip_delimiter = delimiter + ' '
    lines2 = [line.strip().split('#')[0].strip(strip_delimiter) for line in lines
              if line.strip().split('#')[0].strip(strip_delimiter)]

    print(f'lines2[0] = {lines2[0]!r}')
    subcases_line = lines2[0].split(delimiter)
    print(f'subcases_line = {subcases_line!r}')
    subcases_input = []
    for i, subcase in enumerate(subcases_line):
        subcases_input.append(int(subcase))
    nsubcases = len(subcases_input)
    print(f'subcases_input = {subcases_input}')

    combinations = []
    for line in lines2[1:]:
        sline = [val.strip() for val in line.split(delimiter)]
        print(f'split line: {sline}')
        subcase_out = int(sline[0])
        label = sline[1].strip('"\'')
        assert len(label) < 20, label
        factors = [float(factor) for factor in sline[2:]]
        nfactors = len(factors)
        assert len(factors) == nsubcases, f'{label} values={factors} nfactors={nfactors} != nsubcases={nsubcases}'
        combination = (subcase_out, label, factors)
        print(combination)
        combinations.append(combination)
    assert len(combinations)
    return subcases_input, combinations


def run_load_case_combinations(op2_filename: PathLike,
                               combination_filenames: list[PathLike] | PathLike,
                               op2_filenames_new: list[PathLike] | PathLike=None,
                               mode: Optional[str]=None,
                               log: Optional[SimpleLogger]=None) -> None:
    """
    Does a linear combination of a single op2 file.

    Supports multiple combination files to:
     - reduce memory usage
     - reduce file size
     - focus on critical cases

    Parameters
    ----------
    op2_filename : PathLike
        path to input file
    combination_filenames : list[PathLike]
        list of paths to input file
    op2_filenames_new : list[PathLike]; default=None
        list of paths to output files
        must be same length as combination filenames
        default -> 'combination_filename.op2'
    mode : Optional[str]
        the nastran format (msc, nx)
    log : Optional[SimpleLogger]
        a logger object

    """
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)

    if not isinstance(combination_filenames, (list, tuple)):
        combination_filenames = [combination_filenames]

    if op2_filenames_new is None:
        op2_filenames_new = [os.path.splitext(filename)[0] + '.op2'
                             for filename in combination_filenames]
    elif isinstance(op2_filenames_new, PathLike):
        op2_filenames_new = [op2_filenames_new]
    else:
        assert isinstance(op2_filenames_new, (list, tuple)), op2_filenames_new

    all_combinations = load_combinations(combination_filenames)
    print(f'op2_filenames_new = {op2_filenames_new}')
    assert len(op2_filenames_new) == len(combination_filenames)

    # load_geometry: bool = False,
    # combine: bool = True,
    # subcases: Optional[list[int]] = None,
    # exclude_results: Optional[list[str]] = None,
    # include_results: Optional[list[str]] = None,
    # debug: Optional[bool] = True,
    # build_dataframe: Optional[bool] = False,
    # skip_undefined_matrices: bool = True,
    model = read_op2(op2_filename, log=log, mode=mode, combine=True)
    mode_out = model._nastran_format

    # check displacements exist
    assert len(all_combinations) == len(op2_filenames_new)
    disp_keys = list(model.displacements)
    for op2_filename_new, (subcases_in, combinations) in zip(op2_filenames_new, all_combinations):
        for subcase in subcases_in:
            try:
                model.displacements[subcase]
            except KeyError:
                raise KeyError(f'cant find subcase={subcase}; allowed={disp_keys}')

    # do the combinations
    for op2_filename_new, (subcases_in, combinations) in zip(op2_filenames_new, all_combinations):
        combine(model, op2_filename_new,
                subcases_in, combinations, model.log, mode_out)


def get_local_factors(label: str,
                      subcases_in: list[int],
                      factors_in: list[float]) -> tuple[list[int], list[float]]:
    """
    shrink down the number of operations
    [result] = 1.0*[case1] + 2.0*[case2] + 0.0*[case3]
    becomes:
    [result] = 1.0*[case1] + 2.0*[case2]
    """
    factors_out = []
    subcases_out = []
    nsubcases_in = len(subcases_in)
    nfactors_in = len(factors_in)
    assert len(subcases_in) == len(factors_in), f'label={label!r} has a different number of subcases ({nsubcases_in}) than factors ({nfactors_in})'
    for subcase, factor in zip(subcases_in, factors_in):
        if factor == 0.0:
            continue
        factors_out.append(factor)
        subcases_out.append(subcase)
    assert len(factors_out) > 0, f'label={label!r} has no factors != 0'
    return subcases_in, factors_out


def combine(model: OP2,
            op2_filename_new: PathLike,
            subcases_in: list[int],
            combinations: list[tuple],
            log: SimpleLogger,
            mode: str) -> None:
    """Writes a single op2 filename"""
    # assert len(subcases_in) == len(combinations), f'subcases_in={subcases_in} combinations={combinations}'

    unallowed_results = [
        'eigenvectors', 'eigenvalues', 'params', 'gpdt', 'bgpdt', 'eqexin',
        'grid_point_weight', 'psds', 'monitor1', 'monitor3',
        'cstm',
    ]
    table_res_types = list((table_type, model.get_result(table_type)) for table_type in sorted(model.get_table_types())
                            if table_type not in unallowed_results and not table_type.startswith('responses.'))

    model2 = OP2(log=log, mode=mode)
    for combination in combinations:
        (subcase_out, label, factors_in) = combination
        assert len(subcases_in) == len(factors_in), f'subcases_in={subcases_in} factors_in={factors_in}'
        subcases, factors = get_local_factors(label, subcases_in, factors_in)
        subcase0 = subcases[0]

        for table_type, res_type in table_res_types:
            slot = model2.get_result(table_type)
            if len(res_type) == 0:
                continue
            # print(res_type)

            # TODO: we could be smarter here and look for
            #       a case that has a 1.0 factor
            # TODO: is there a way to eliminate the deepcopy?
            case0 = res_type[subcase0]
            case_new = copy.deepcopy(case0)
            # case_new.isubcase = subcase_out
            case_new.label = label

            case_new.linear_combination(0.0, update=False)
            # case_new.scale_data(casei.data, factor)
            for subcase, factor in zip(subcases, factors):
                casei = res_type[subcase]
                case_new.linear_combination(factor, casei.data, update=False)

            if hasattr(case_new, 'update_data_components'):
                case_new.update_data_components()
            else:
                print(type(case_new))

            from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
            slot[subcase_out] = case_new
        break # end of first combination

    print(f'op2_filename_new = {op2_filename_new}')
    assert mode is not None, mode
    model2._nastran_format = mode
    model2.write_op2(op2_filename_new)
    return


def main():  # pragma: no cover
    from pathlib import Path
    dirname = Path(__file__).parent

    op2_filename = dirname / 'test load combo-000.op2'
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    combination_filename = dirname / 'combination_file_real.txt'
    run_load_case_combinations(op2_filename, combination_filename)


if __name__ == '__main__':  # pragma: no cover
    main()
