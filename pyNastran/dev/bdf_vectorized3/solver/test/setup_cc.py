from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf
from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase


def setup_static_case_control(model: BDF,
                              extra_case_lines=None,
                              thermal_load_id: int=0,
                              load_id: int=2) -> None:
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'ESE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
    ]
    if thermal_load_id:
        lines.append(f'  TEMPERATURE(LOAD) = {thermal_load_id}')
    if load_id:
        lines.append(f'  LOAD = {load_id}')
    lines.append('  SPC = 3')

    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 101
    model.case_control_deck = cc


def setup_modes_case_control(model: BDF, extra_case_lines=None,
                             nmodes: int=10):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'ESE(PLOT,PRINT) = ALL',
        'SUBCASE 1',
        #'  LOAD = 2',
        '  SPC = 3',
        '  METHOD = 42',
    ]
    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 103
    model.case_control_deck = cc
    model.add_eigrl(42, nd=nmodes)


def setup_buckling_case_control(model: BDF, extra_case_lines=None,
                                nmodes: int=10, nsubcases=1):
    lines = [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'OLOAD(PLOT,PRINT) = ALL',
        'ESE(PLOT,PRINT) = ALL',
    ]
    if nsubcases == 0:
        lines += [
            'LOAD = 2',
            'SPC = 3',
            'METHOD = 42',]
    elif nsubcases == 1:
        lines += [
            'SUBCASE 1',
            '  LOAD = 2',
            '  SPC = 3',
            '  METHOD = 42',]
    else:
        assert nsubcases == 2, nsubcases
        lines += [
            'SUBCASE 1',
            '  LOAD = 2',
            '  SPC = 3',
            'SUBCASE 2',
            '  STATSUB(PRELOAD) = 1',
            '  METHOD = 42',]

    if extra_case_lines is not None:
        lines += extra_case_lines
    cc = CaseControlDeck(lines, log=model.log)
    model.sol = 105
    model.case_control_deck = cc
    model.add_eigrl(42, nd=nmodes)

