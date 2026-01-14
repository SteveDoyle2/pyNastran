import os
import copy
from typing import Optional
import numpy as np
from cpylog import SimpleLogger
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.op2.op2 import read_op2, OP2

# (imodel, subcases, factors)
# (0, [1, 1], [1.2, 0.0])
CombinationPair = tuple[int, list[int], list[float]]
# [(0, [1, 1], [1.2, 0.0])]
CombinationPairs = list[CombinationPair]

# (subcase_out, label, pairs)
MultiCombination = tuple[int, str, CombinationPairs]

def load_combinations(combination_filenames: list[PathLike],
                      delimiter: str=',') -> list:
    """

    Combination File
    ----------------
    # line 1 is the subcase ids  for the scale factors
    # line 2 is a comment

    # this is line 0. line 1 has the input subcase ids
      Null    Subcases  1      2      3
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
                      delimiter: str=',') -> tuple[list[int], list]:
    lines = _get_combination_lines(combination_filename, delimiter)

    log = SimpleLogger(level='debug')
    log.debug(f'lines[0] = {lines[0]!r}')
    subcases_sline = lines[0].split(delimiter)
    log.debug(f'subcases_sline = {subcases_sline!r}')

    # chop the output_subcase, label columns
    subcases_input = []
    for i, subcase in enumerate(subcases_sline[2:]):
        subcases_input.append(int(subcase))
    nsubcases = len(subcases_input)
    log.debug(f'subcases_input = {subcases_input}')

    combinations = []
    for line in lines[1:]:
        sline = [val.strip() for val in line.split(delimiter)]
        log.debug(f'split line: {sline}')
        subcase_out = int(sline[0])
        label = sline[1].strip('"\'')
        # assert len(label) < 20, label
        factors = [float(factor) for factor in sline[2:]]
        nfactors = len(factors)
        assert len(factors) == nsubcases, f'{label} values={factors} nfactors={nfactors} != nsubcases={nsubcases}'
        combination = (subcase_out, label, factors)
        log.debug(str(combination))
        combinations.append(combination)
    assert len(combinations)
    return subcases_input, combinations


def load_multi_combinations(
        combination_filenames: list[PathLike],
        log: SimpleLogger,
        delimiter: str=',') -> list[tuple[int, list[MultiCombination]]]:
    """

    Combination File
    ----------------
    # line 1 is the subcase ids  for the scale factors
    # line 2 is a comment
      Null,    iOP2,      1,      2,      3,   # op2 ids
      Null,    SubcaseID, 1,      2,      3,   # subcase ids
    # Subcase, Name,      Scale1, Scale2, Scale3
    10,        case10,    1.0,    2.0,    3.0  # comment
    20,        case20,    1.2,    2.2,    3.2
    30,        "case 30", 1.2,    2.2,    3.2
    """
    assert len(combination_filenames), combination_filenames
    for combination_filename in combination_filenames:
        assert os.path.exists(combination_filename), print_bad_path(combination_filename)

    all_combinations = []
    for combination_filename in combination_filenames:
        subcases_input, combinations = _load_multi_combination(
            combination_filename, log, delimiter=delimiter)
        all_combinations.append((subcases_input, combinations))
    assert len(all_combinations)
    return all_combinations


def _get_combination_lines(combination_filename: PathLike,
                           delimiter: str) -> list[str]:
    """remove comments and blank lines"""
    with open(combination_filename, 'r') as combination_file:
        lines = combination_file.readlines()

    lines2 = [line.split('#')[0].strip().rstrip(delimiter) for line in lines
              if line.split('#')[0].strip().rstrip(delimiter)]
    return lines2


def string_where(op2s_input: list[str], op2_input: str) -> list[int]:
    """
    A string capable version for:
    iop2 = np.where(op2s_input == op2_input)[0]
    """
    iop2 = []
    for i, op2_inputi in enumerate(op2s_input):
        if op2_inputi == op2_input:
            iop2.append(i)
    return iop2

def _load_multi_combination(
        combination_filename: PathLike,
        log: SimpleLogger,
        delimiter: str=',') -> tuple[np.ndarray,
                                     tuple[list[int], list[MultiCombination]]]:
    lines = _get_combination_lines(combination_filename, delimiter)
    log.debug(f'lines[0] = {lines[0]!r}')

    sline0 = lines[0].split(delimiter)
    sline1 = lines[1].split(delimiter)
    op2s_sline = sline0[2:]
    subcases_sline = sline1[2:]
    log.debug(f'op2s_sline     = {op2s_sline!r}')
    log.debug(f'subcases_sline = {subcases_sline!r}')
    assert len(op2s_sline) == len(subcases_sline)
    assert len(subcases_sline) > 0, sline1
    op2s_input = []
    subcases_input_list = []
    for iop2, subcase in zip(op2s_sline, subcases_sline):
        # iop2 = int(iop2)
        op2s_input.append(iop2.strip())
        subcases_input_list.append(int(subcase))

    subcases_input = np.array(subcases_input_list, dtype=int)
    iop2s = []
    uop2_input = np.unique(op2s_input)
    for op2_input in uop2_input:
        iop2 = string_where(op2s_input, op2_input)
        # iop2 = np.where(op2s_input == op2_input)[0]
        iop2s.append(iop2)
    nsubcases = len(subcases_input)
    log.debug(f'subcases_input = {subcases_input}')

    combinations = []
    for line in lines[2:]:
        sline = [val.strip() for val in line.split(delimiter)]
        log.debug(f'split line: {sline}')
        subcase_out = int(sline[0])
        label = sline[1].strip('"\'')
        assert len(label) < 20, label
        factors = np.array([float(factor) for factor in sline[2:]])
        nfactors = len(factors)
        assert len(factors) == nsubcases, f'{label} values={factors} nfactors={nfactors} != nsubcases={nsubcases}'
        combination_pairs = []
        for op2_input, iop2 in zip(uop2_input.tolist(), iop2s):
            factorsi = factors[iop2].tolist()
            subcasesi = subcases_input[iop2].tolist()
            combination_pairs.append((op2_input, subcasesi, factorsi))
        assert len(combination_pairs), combination_pairs
        combination = (subcase_out, label, combination_pairs)
        log.debug(str(combination))
        combinations.append(combination)
    assert len(combinations)

    return subcases_input, combinations


def _check_op2_file(op2_filename: PathLike | OP2) -> None:
    if isinstance(op2_filename, PathLike):
        assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    elif isinstance(op2_filename, OP2):
        pass
    else:
        raise TypeError(type(op2_filename))


def run_load_case_combinations(op2_filename: PathLike | OP2,
                               combination_filenames: list[PathLike] | PathLike,
                               exclude_results: Optional[list[str]]=None,
                               include_results: Optional[list[str]]=None,
                               mode: Optional[str]=None,
                               revision: Optional[str]=None,
                               log: Optional[SimpleLogger]=None) -> None:
    """
    Create a series of load case combinations based on a single op2

    Parameters
    ----------
    op2_filename : PathLike
        path to input file
    combination_filenames : list[PathLike]
        list of paths to input file
    mode : Optional[str]
        the nastran format (msc, nx)
    log : Optional[SimpleLogger]
        a logger object

    """
    _check_op2_file(op2_filename)
    if not isinstance(combination_filenames, (list, tuple)):
        combination_filenames = [combination_filenames]

    all_combinations = load_combinations(combination_filenames)
    op2_filenames_new = [os.path.splitext(filename)[0] + '.op2'
                         for filename in combination_filenames]
    print(f'op2_filenames_new = {op2_filenames_new}')

    # load_geometry: bool = False,
    # combine: bool = True,
    # subcases: Optional[list[int]] = None,
    # debug: Optional[bool] = True,
    # build_dataframe: Optional[bool] = False,
    # skip_undefined_matrices: bool = True,
    model = _read_op2(
        op2_filename, mode=mode, log=log,
        exclude_results=exclude_results,
        include_results=include_results)
    mode_out = model._nastran_format
    model._nastran_revision = revision

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
                subcases_in, combinations, model.log, mode_out, revision)


def _read_op2(op2_filename: PathLike | OP2,
              mode: str,
              log: SimpleLogger,
              exclude_results: Optional[list[str]] = None,
              include_results: Optional[list[str]] = None,
              ) -> OP2:
    if isinstance(op2_filename, OP2):
        return op2_filename
    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    model = read_op2(
        op2_filename, log=log,
        exclude_results=exclude_results,
        include_results=include_results,
        mode=mode, combine=True)
    return model

def run_load_case_multi_combinations(
        op2_filenames: list[PathLike | OP2] | dict[int, PathLike | OP2],
        combination_filenames: list[PathLike] | PathLike,
        exclude_results: Optional[list[str]]=None,
        include_results: Optional[list[str]]=None,
        mode: Optional[str]=None,
        revision: Optional[str]=None,
        log: Optional[SimpleLogger]=None) -> None:
    """

    Parameters
    ----------
    op2_filenames : list[PathLike | OP2]; dict[int, PathLike | OP2]
        path to input file or a loaded op2
    combination_filenames : list[PathLike]
        list of paths to input file
    mode : Optional[str]
        the nastran format (msc, nx)
    log : Optional[SimpleLogger]
        a logger object

    """
    log = log if log is not None else SimpleLogger(level='info')
    if isinstance(op2_filenames, (list, tuple)):
        op2_filenames = {i: op2_filename for i, op2_filename in enumerate(op2_filenames)}
    elif isinstance(op2_filenames, dict):
        pass
    else:
        raise TypeError(op2_filenames)
    # log.info(f'op2_filenames (dict) = {op2_filenames}')

    if not isinstance(combination_filenames, (list, tuple)):
        combination_filenames = [combination_filenames]

    all_combinations = load_multi_combinations(combination_filenames, log)
    ncombinations = sum([len(subcasesin_pairs[1])
                         for subcasesin_pairs in all_combinations])
    op2_filenames_new = [os.path.splitext(filename)[0] + '.op2'
                         for filename in combination_filenames]
    log.info(f'op2_filenames_new = {op2_filenames_new}')
    assert len(all_combinations) == len(op2_filenames_new)
    #-------------------------------------------------------

    models, imodel = _load_models(
        op2_filenames, all_combinations, mode=mode, log=log,
        exclude_results=exclude_results,
        include_results=include_results,
    )
    log = models[imodel].log

    # do the combinations
    icombination = 0
    for op2_filename_new, (subcases_in, combination_pairs) in zip(op2_filenames_new, all_combinations):
        icombination = multi_combine(
            models, op2_filename_new,
            combination_pairs, log, mode, revision,
            icombination, ncombinations=ncombinations)


def _load_models(op2_filenames: dict[int, PathLike | OP2],
                 all_combinations: list[tuple[int, list[MultiCombination]]],
                 exclude_results: Optional[list[str]]=None,
                 include_results: Optional[list[str]]=None,
                 log: Optional[SimpleLogger]=None,
                 mode: Optional[str]=None) -> tuple[dict[int, OP2], int]:
    models = {}
    assert len(op2_filenames) > 0, op2_filenames
    for imodel, op2_filename in op2_filenames.items():
        model = _read_op2(
            op2_filename, mode=mode, log=log,
            exclude_results=exclude_results,
            include_results=include_results)
        # model._nastran_revision = revision
        models[imodel] = model

        # check displacements exist
        disp_keys = list(model.displacements)
        for (subcases_in, combinations) in all_combinations:
            for subcase in subcases_in:
                try:
                    model.displacements[subcase]
                except KeyError:
                    log.level = 'info'
                    log.info(str(model.get_op2_stats(short=True)))
                    raise KeyError(f'cant find subcase={subcase} in {subcases_in}; allowed={disp_keys}')
    return models, imodel


def get_local_factors(label: str,
                      subcases_in: list[int],
                      factors_in: list[float],
                      require_cases: bool=True) -> tuple[list[int], list[float]]:
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
    if require_cases:
        assert len(factors_out) > 0, f'label={label!r} has no factors != 0'
    return subcases_in, factors_out


def multi_combine(models: dict[int, OP2],
                  op2_filename_new: PathLike,
                  combinations: list[MultiCombination],
                  log_op2: SimpleLogger,
                  mode: str, revision: Optional[str],
                  icombination: int=0, ncombinations: int=0) -> int:
    assert isinstance(models, dict), models
    models = {str(key).strip(): model for key, model in models.items()}

    icombination0 = icombination
    if ncombinations == 0:
        ncombinations = len(combinations)

    log = SimpleLogger(level='debug')
    out_model = OP2(log=log_op2, mode=mode)
    out_model._nastran_revision = revision
    # table_res_types = _results(out_model)

    imodel0 = list(models)[0]
    model0 = models[imodel0]
    table_res_types0 = _results(model0)

    # find all subcase associated with model0 (imodel0)
    # (imodel, subcases_in, factors_in) in enumerate(imodel_subcases_factors)
    # (10, 'case10', [(0, [1, 1], [1.0, 2.0])])
    combination0 = combinations[0]
    (subcase_out0, label0, imodel_subcases_factors0) = combination0
    subcases0 = [isf[1] for isf in imodel_subcases_factors0
                 if isf[0] == imodel0]
    assert len(subcases0), f'cant find imodel={imodel0}'
    # now grab a single subcase
    subcase0 = subcases0[0][0]
    # print(f'imodel0 = {imodel0}')
    # print(f'subcases0 = {subcases0}')
    # print(f'subcase0 = {subcase0}')

    for table_type, res_type in table_res_types0:
        if len(res_type) == 0:
            continue
        out_slot = out_model.get_result(table_type)

        for icombinationi, combination in enumerate(combinations):
            icombination = icombination0 + icombinationi
            log.info(f'{table_type} {icombination}/{ncombinations} = {combination}')
            (subcase_out, label, imodel_subcases_factors) = combination

            ires0 = 0

            # imodel = imodel_subcases_factors[0]
            # assert isinstance(imodel, int), imodel
            assert isinstance(models, dict), models
            slot = model0.get_result(table_type)

            # pick a displacement case and 0 it out
            # we only 0 out linear combinable data
            # so not fiber_distance
            case0 = slot[subcase0]
            case_new = copy.deepcopy(case0)
            case_new.label = label

            # case.data *= 0
            case_new.linear_combination(0.0, update=False)

            for ires, (imodel, subcases_in, factors_in) in enumerate(imodel_subcases_factors):
                assert len(subcases_in) == len(factors_in), f'subcases_in={subcases_in} factors_in={factors_in}'
                subcases, factors = get_local_factors(
                    label, subcases_in, factors_in, require_cases=False)
                if len(subcases) == 0:
                    continue
                subcase0 = subcases[0]

                model = models[imodel]
                slot = model.get_result(table_type)
                if len(slot) == 0:
                    raise RuntimeError(f'missing {table_type!r} for op2={model.op2_filename} (subcases={subcases_in})')
                    # continue

                # add data
                for subcase, factor in zip(subcases, factors):
                    casei = res_type[subcase]
                    case_new.linear_combination(factor, casei.data, update=False)
                    #casei = caseii * factor
                    #case.data = casei.data
                ires0 += 1

            if hasattr(case_new, 'update_data_components'):
                case_new.update_data_components()
            else:
                log.warning(f'skipping {str(type(case_new))}')
            out_slot[subcase_out] = case_new
            # break # end of first combination
        # break # end of all displacement combinations

    print(f'op2_filename_new = {op2_filename_new}')
    assert mode is not None, mode
    out_model._nastran_format = mode
    out_model.write_op2(op2_filename_new)
    return icombination


def combine(model: OP2,
            op2_filename_new: PathLike,
            subcases_in: list[int],
            combinations: list[tuple],
            log: SimpleLogger,
            mode: str, revision: Optional[str]):
    # assert len(subcases_in) == len(combinations), f'subcases_in={subcases_in} combinations={combinations}'

    table_res_types = _results(model)

    model2 = OP2(log=log, mode=mode)
    model2._nastran_revision = revision
    for combination in combinations:
        (subcase_out, label, factors_in) = combination
        assert len(subcases_in) == len(factors_in), f'subcases_in={subcases_in} factors_in={factors_in}'
        subcases, factors = get_local_factors(
            label, subcases_in, factors_in, require_cases=True)
        subcase0 = subcases[0]

        for table_type, res_type in table_res_types:
            slot = model2.get_result(table_type)
            if len(res_type) == 0:
                continue
            # print(res_type)

            # displacement
            case0 = res_type[subcase0]
            case_new = copy.deepcopy(case0)
            case_new.label = label

            # case.data *= 0
            case_new.linear_combination(0.0, update=False)
            for subcase, factor in zip(subcases, factors):
                casei = res_type[subcase]
                case_new.linear_combination(factor, casei.data, update=False)
                #casei = caseii * factor
                #case.data = casei.data

            if hasattr(case_new, 'update_data_components'):
                case_new.update_data_components()
            else:
                print(type(case_new))

            slot[subcase_out] = case_new
        # break # end of first combination

    print(f'op2_filename_new = {op2_filename_new}')
    assert mode is not None, mode
    model2._nastran_format = mode
    model2.write_op2(op2_filename_new)
    return


def _results(model: OP2) -> list[dict]:
    unallowed_results = [
        'eigenvectors', 'eigenvalues', 'params', 'gpdt', 'bgpdt', 'eqexin',
        'grid_point_weight', 'psds', 'monitor1', 'monitor3',
        'cstm',
    ]
    table_res_types = list((table_type, model.get_result(table_type)) for table_type in sorted(model.get_table_types())
                            if table_type not in unallowed_results and not table_type.startswith('responses.'))
    return table_res_types


def main():  # pragma: no cover
    from pathlib import Path
    dirname = Path(__file__).parent
    op2_filename = dirname / 'test load combo-000.op2'
    combination_filename = dirname / 'combination_file_real.txt'

    assert os.path.exists(op2_filename), print_bad_path(op2_filename)
    combinations = load_combinations([combination_filename], delimiter=',')

    run_load_case_combinations(op2_filename, combination_filename, mode='nx', revision='2306')


if __name__ == '__main__':  # pragma: no cover
    main()
