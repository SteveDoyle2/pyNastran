from typing import Any
import numpy as np
from cpylog import get_logger, SimpleLogger

class OptimizationResult():
    def __init__(self):
        self.design_objective = {
            'label': [],
            'in': [],
            'out': [],
            'subcase_id': [],
        }
        #all_desvars = []
        #objective_functions = {}
        self.design_vars = {
            'internal_id' : [],
            'desvar_id' : [],
            'label' : [],
            'xl' : [],
            'xi' : [],
            'x' : [],
            'xu' : [],
        }
    def __repr__(self):
        design_objective = self.design_objective
        internal_id = np.array(self.design_vars['internal_id'], dtype='int32')
        inputs = np.array(design_objective['in'], dtype='float64')
        outputs = np.array(design_objective['out'], dtype='float64')
        obj = (
            'obj:\n'
            f"  label = {design_objective['label']}; n={len(design_objective['label'])}\n"
            f"  id    = {internal_id}\n"
            f"  in    = {inputs}\n"
            f"  out   = {outputs}\n")

        print(obj)
        return f'OptimizatioResult(obj={obj}, ndesvars={len(internal_id)})'

def _read_line_block(i: int, lines: list[str],
                     stop_marker: str='',
                     rstrip: bool=False, strip: bool=False, debug: bool=False,
                     imax=None) -> tuple[int, list[str]]:
    i0 = i
    lines2 = []
    line = lines[i].rstrip()
    #print('cat', rstrip, strip, line)
    if rstrip:
        #print('****rstrip')
        line = lines[i].rstrip()
        while line != stop_marker:
            if debug:
                print(i, line)
            if imax and i > imax:
                raise RuntimeError('i=%i' % i)
            i += 1
            line = lines[i].rstrip()
            lines2.append(line)
    elif strip:
        line = lines[i].strip()
        #print('****strip')
        i0 = 0
        #print('catt', line.strip())
        if line.startswith('1)'):
            lines2.append(line.strip())

        while line != stop_marker or i0 == 0:
            #if debug:
                #print(i, line)
            if imax and i > imax:
                raise RuntimeError('i=%i' % i)
            i += 1
            if 'MSC.NASTRAN' in line:
                break
            if 'F I N A L   A N A L Y S I S' in line:
                if debug:
                    stopping
                _read_line_block(i0, lines, stop_marker=stop,
                                 rstrip=rstrip, strip=strip, debug=False)
            #print(i, line)
            line = lines[i].strip()
            lines2.append(line)
            i0 += 1
    else:
        line = lines[i].rstrip('\n')
        while line != stop_marker:
            if debug:
                print(i, line)
            if imax and i > imax:
                raise RuntimeError('i=%i' % i)
            i += 1
            line = lines[i].rstrip('\n')
            lines2.append(line)
    try:
        lines2.pop()
    except Exception:
        print('***', line)
        print(lines2)
        raise

    return i, lines2

def _read_startswith_line_block(i: int, lines: list[str], stop_marker='',
                                rstrip=False, strip=False, debug=False, imax=None):
    assert strip or rstrip
    lines2 = []
    if rstrip:
        line = lines[i].rstrip()
        while not line.startswith(stop_marker):
            if debug:
                ab
                print(i, line)
            if imax and i > imax:
                raise RuntimeError('i=%i' % i)
            i += 1
            line = lines[i].rstrip()
            lines2.append(lines2)
    elif strip:
        line = lines[i].strip()
        while not line.startswith(stop_marker):
            if debug:
                abc
                print(i, line)
            if imax and i > imax:
                raise RuntimeError('i=%i' % i)
            i += 1
            line = lines[i].strip()
            lines2.append(lines2)
    else:
        raise NotImplementedError(f'rstrip={rstrip} strip={strip}')
    print(f'!!! {line}')
    return i, lines2

def _goto_page(i: int, lines: list[str], debug: bool=False):
    line = lines[i]
    i0 = i
    if debug:
        print('i0', i, line)
    try:
        while 'PAGE' not in line:
            if debug:
                print(i, line)
            i += 1
            line = lines[i]
    except IndexError:
        try:
            _goto_page(i, lines, debug=True)
        except IndexError:
            pass
        raise
    return i

def _read_int_gradient(i: int, line: str, lines: list[str], nlines: int, log: SimpleLogger, debug: bool=True):
    line = lines[i].strip()
    #if debug:
    log.debug(f'{i} read_int_gradient: {line}')
    i, lines2 = _read_line_block(i, lines, stop_marker='', rstrip=False, strip=True,
                                 debug=True, imax=None)
    line = lines[i].strip()
    #for line in lines2:
        #print(line)
    ids = _parse_ints(lines2)
    #print('end read_int_gradient', i, ids)
    return i, ids

def _read_gradient(i: int, line: str, lines: list[str], nlines, log: SimpleLogger, debug: bool=True):
    i += 1
    #if debug:
    if 'GRADIENT OF CONSTRAINT NUMBER' not in line:
        log.debug(f'{i} read_gradient: {line}')
    #print(i, line)
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, stop_marker='', rstrip=False, strip=True,
                                 debug=False, imax=None)
    line = lines[i].strip()
    if debug:
        for line in lines2[:5]:
            log.debug(line)
    grad = _parse_gradient(lines2, debug=debug)
    #print('end gradient', i, grad)
    return i, grad

def _read_gradient_block(i, line, lines, nlines, log, debug=True):
    i += 1
    if debug:
        log.debug(f'{i} read_gradient: {line}')
    #print(i, line)
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, stop_marker='', rstrip=False, strip=True,
                                 debug=False, imax=None)
    line = lines[i].strip()
    return i, lines2

def _read_design_optimization(i, line, lines, nlines,
                              idesign: int,
                              all_results: list[OptimizationResult],
                              log: SimpleLogger):
    optimization_result = all_results[-1]
    assert isinstance(optimization_result, OptimizationResult), optimization_result
    end_of_job = False
    i += 1
    log.debug(f'{i} read_design_optimization')
    #log.debug(f'{i} {line}')
    line = lines[i].strip()
    gradient_words = [
        'CONSTRAINT VALUES (G-VECTOR)',
        #'DESIGN VARIABLE SCAL
        # E FACTORS',
        'DG(I)/DALPHA VECTOR',
        'SEARCH DIRECTION (S-VECTOR)',
        #'GRADIENT OF THE OBJECTIVE FUNCTION (DF-VECTOR)',
        'DECISION VARIABLES (X-VECTOR)',
        #'LOWER BOUNDS ON THE DECISION VARIABLES (XL-VECTOR)',
        #'UPPER BOUNDS ON THE DECISION VARIABLES (XU-VECTOR)',
    ]
    blank_gradients = [
        'LAGRANGE MULTIPLIERS',
        'CONSTRAINTS, G(X)',
    ]
    int_tables = [
        'CONSTRAINT NUMBERS',
        'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)',
        'RETAINED ACTIVE/VIOLATED CONSTRAINT NUMBERS',
    ]
    from collections import defaultdict
    gradients = defaultdict(list)
    stop_gradients = False
    #print('*'*80)
    iteration = 0
    constraint_ids = None
    gvector = None

    design_cycle = None
    while i < nlines:
        #if line:
            #log.debug(line)
        print_line = True
        line0 = line
        if line == '*                D E S I G N    O P T I M I Z A T I O N            *':
            i -= 1
            #log.debug('breaking %s' % line)
            break
        #elif 'F I N A L   A N A L Y S I S' in line:
            #i -= 1
            #log.debug('breaking %s' % line)
            #break

        elif line == '-----   COMPARISON BETWEEN INPUT PROPERTY VALUES FROM ANALYSIS AND DESIGN MODELS -----':
            i, line, lines2 = _read_property_comparison_table(i, lines, log)
            i += 7
        elif line == '-----   COMPARISON BETWEEN INPUT CONNECTIVITY PROPERTY VALUES FROM ANALYSIS AND DESIGN MODELS -----':
            i, line, lines2 = _read_property_table(i, lines)
            #print('*', i, line)
        elif line == '*******   ANALYSIS RESULTS BASED ON THE INITIAL DESIGN   *******':
            i += 1
            line = lines[i].strip()

            i, is_broken = _find_next_table(i, line, lines)
            if is_broken:
                i -= 1
                continue
            i, lines2 = _read_line_block(i, lines, 'DESIGN OPTIMIZATION TOOLS', strip=True, debug=False)
            line = lines[i].strip()

            #i, lines2 = read_line_block(i, lines, 'CONTROL PARAMETERS', strip=True, debug=False, imax=1000)
            #line = lines[i].strip()
            #print('*', i, line)

            #i += 1
            #line = lines[i].strip()
            #print(i, line)
            i, lines2 = _read_startswith_line_block(i, lines, 'GRADIENT CALLS =', strip=True, debug=False)
            line = lines[i].strip()
            log.debug(f'* {i} {line}')
            i = _goto_page(i, lines)
            line = lines[i].strip()
            log.debug(f'* {i} {line}')
        elif '* * * * * BEGIN ONE-DIMENSIONAL SEARCH * * * * *' in line:
            i = _read_1d_search(i, line, lines, nlines, log)
            line = lines[i].strip()
            #log.debug(f'* {i} {line}')
        #elif '-- BEGIN CONSTRAINED OPTIMIZATION: MFD METHOD' in line:
            #i = read_constrained_optimization(i, line, lines, nlines)
            #line = lines[i].strip()
            #print('*', i, line)
        elif '-- ITERATION NUMBER' in line:
            iteration = int(line.split()[-1])
            i = _read_iteration_number(i, line, lines, design_vars, constraint_ids, gvector, nlines, log)
            line = lines[i].strip()
            #print('*', i, line)
        elif line in gradient_words:
            if line in gradients:
                ngrad = len(gradients[line])
                skip_ = [
                    'SEARCH DIRECTION (S-VECTOR)',
                    'DG(I)/DALPHA VECTOR',
                    'LOWER BOUNDS ON THE DECISION VARIABLES (XL-VECTOR)',
                    'UPPER BOUNDS ON THE DECISION VARIABLES (XU-VECTOR)',
                    'DECISION VARIABLES (X-VECTOR)',
                    'CONSTRAINT VALUES (G-VECTOR)',
                ]
                if line in skip_:
                    pass
                    #if ngrad + 1 != iteration:
                        #print(f"'{line}' {ngrad} {iteration}")
                        #stop_gradients = True
                elif ngrad != iteration:
                    log.info(f'{line} {ngrad} {iteration}')
                    stop_gradients = True
                    adfsdf
            i -= 1
            i, grad = read_gradient(i, line, lines, nlines, log, debug=False)
            gradients[line].append(grad)
        elif line in  blank_gradients:
            # there is a blank line
            i += 1
            log.info(line)
            i, lines2 = _read_gradient(i, line, lines, nlines, log, debug=False)
        elif 'DECISION VARIABLES, X' in line:
            i += 1
            i, lines2 = _read_gradient_block(i, line, lines, nlines, log)
            header, ids, data = _get_decision_variables(lines2)
        elif 'THE GRADIENT OF CONSTRAINT NUMBERS' in line:
            # THE GRADIENT OF CONSTRAINT NUMBERS      4779 AND      4778 ARE DEPENDENT
            # CONSTRAINT NUMBER      4779 IS REMOVED FROM THE ACTIVE SET
            log.debug(line)
            i += 1
        elif 'GRADIENT OF CONSTRAINT NUMBER' in line:
            igradient = int(line.split()[-1])
            i, lines2 = _read_gradient(i, line, lines, nlines, log)
        elif line in int_tables:
            i, ids = _read_int_gradient(i, line, lines, nlines, log)
            gradients[line] = ids
        elif 'K-T PARAMETERS, BETA =' in line:
            # K-T PARAMETERS, BETA =  1.00000E+00  MAX. RESIDUAL =  3.36450E-04
            *trash, beta, maxstr, resstr, eq2, residual = line.split()
            del trash, maxstr, resstr, eq2
            beta = float(beta)
            residual = float(residual)
        elif 'ALPHA =' in line:
            # CALCULATED ALPHA =  1.66390E-02
            # PROPOSED ALPHA =  2.06505E-04
            # ALPHA =  2.06505E-04
            alpha = line.split()[-1]
            alpha = float(alpha)

        #elif ')  ' not in line:
            #print('b', i, line)
        elif _is_page_skip(line):
            print_line = False
            try:
                i = _goto_page(i, lines)
            except Exception:
                log.debug(f'{i} {line}')
                raise
            line = lines[i].strip()
        elif '-- SCALAR PROGRAM PARAMETERS' in line:
            i += 22
        elif 'OBJ =' in line:
            # OBJ =  2.63781E+01
            #assert obj is None, line
            obj = line.split()[-1]
            obj = float(obj)
        elif 'OBJECTIVE =' in line:
            # OBJECTIVE =  2.63781E+01
            #assert obj is None, line
            obj = line.split()[-1]
            obj = float(obj)
        elif '-- OPTIMIZATION RESULTS' in line:
            i += 3
            line = lines[i]
            assert 'OBJECTIVE, F(X) =' in line, line
            obj = line.split()[-1]
            obj = float(obj)
            #print(i, line)
            #aaa
        elif '*****   OPTIMIZATION RESULTS BASED ON THE APPROXIMATE MODEL   *****' in line:
            #log.debug(f'{i} {line}')
            i += 3
            line = lines[i].strip()
            #log.debug(f'{i} {line}')
            assert '-----   DESIGN OBJECTIVE   -----' in line, line
            i, line = _read_design_objective(i, lines, optimization_result, log)

            i += 2
            line = lines[i].strip()
            #log.debug(f'{i} {line}')
            i, line = _read_design_variables(
                i, lines,
                design_cycle, idesign, optimization_result, log)
            line = lines[i].strip()
            #print_design_vars(design_vars)
            #print(design_vars.keys())
            #labels = optimization_result.design_vars['label']
            #dfasf


        elif line == '':
            print_line = False
        elif '31)  ' in line:
            asdf
        #elif 'INSPECTION OF CONVERGENCE DATA FOR THE OPTIMAL DESIGN WITH RESPECT TO APPROXIMATE MODELS' in line:
            #print(i, line)
            #adsf
        elif _is_skip_one_line(line):
            print_line = False
            i += 1
        elif _is_info_msg(line):
            print_line = False
            #log.debug(line)
            pass
        #elif 'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (READ MODULE)' in line:
            #i = goto_page(i, lines)
            #line = lines[i].strip()
        elif _is_end_of_job(line):
            end_of_job = True
            break
        elif '-----   DESIGN OBJECTIVE   -----' in line:
            #log.debug(line)
            if idesign != -1:
                log.warning('reset...')
                all_results.append(optimization_result)
                optimization_result = OptimizationResult()
            i, line = _read_design_objective(i, lines, optimization_result, log)
            assert isinstance(i, int), line
        elif '-----   DESIGNED PROPERTIES   -----' in line:
            #log.debug(line)
            i, line, lines2 = _read_designed_properties(i, lines)
            assert isinstance(i, int), line
        elif '-----   DESIGN VARIABLES   -----' in line:
            print_line = False
            #log.debug(line)
            #print('idesign', idesign)
            i, line = _read_design_variables(
                i, lines,
                design_cycle, idesign, optimization_result, log)
            if design_cycle == 3:
                print([res.design_vars['internal_id'][0] for res in all_results])
                print([res.design_vars['xu'][0] for res in all_results])
                print([res.design_vars['xl'][0] for res in all_results])
                print([res.design_vars['x'][0] for res in all_results])
                asdf
            idesign += 1
            line = lines[i].strip()
            #show_design_vars(design_vars, log)
            assert isinstance(i, int), line
            #aaa
        elif '-----   DESIGN CONSTRAINTS ON RESPONSES   -----' in line:
            #log.debug(line)
            i, line, lines2 = _read_design_constraints_of_responses(i, lines)
            assert isinstance(i, int), line
        elif '-----   CONSTRAINTS ON DESIGNED PROPERTIES   -----' in line:
            #log.debug(line)
            i, line, lines2 = _read_constraints_on_designed_properties(i, lines)
            assert isinstance(i, int), line
        elif '-----   WEIGHT AS A FUNCTION OF MATERIAL ID ----' in line:
            #log.debug(line)
            i, line, lines2 = _read_weight_as_a_function_of_material_id(i, lines)
            assert isinstance(i, int), line
        elif '-----    STRESS RESPONSES    -----' in line:
            #log.debug(f'stress response {i} {line}')
            i, line, lines2 = _read_stress_responses(i, lines)
            #log.debug(f'stress response2 {i} {line}')
        elif '-----    COMPOSITE LAMINATE STRESS RESPONSES    -----' in line:
            log.debug(line)
            i, line, lines2 = _read_composite_stress_responses(i, lines)
        elif '-----    WEIGHT RESPONSE    -----' in line:
            #log.debug(line)
            i, line, lines2 = _read_weight_responses(i, lines)
        elif '-----    FLUTTER RESPONSES    -----' in line:
            #log.debug(line)
            i, line, lines2 = _read_flutter_responses(i, lines)
        elif '---- RETAINED DRESP2 RESPONSES ----' in line:
            #log.debug(line)
            i, line, lines2 = _read_retained_dresp2_responses(i, lines)
        elif 'S U M M A R Y   O F   D E S I G N    C Y C L E    H I S T O R Y' in line:
            i, end_of_job = _read_summary_of_design_cycle_history(i, line, lines, nlines, log)
            if end_of_job:
                break
        elif  line == 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R':
            i += 31
            line = lines[i].strip()  # page xx
        elif _is_page_skip(line):
            print_line = False
            i = _goto_page(i, lines)
            line = lines[i].strip()
        elif 'BLOCK SIZE USED ......................' in line:
            i += 6
        elif '-- OPTIMIZATION IS COMPLETE' in line:
            i = _complete_optimization(i, lines, nlines, log)
        elif 'DESIGN OPTIMIZATION TOOLS' in line:
            i = _read_design_optimization_tools(i, line, lines, nlines,
                                                design_vars, log)
        elif 'D E S I G N   C Y C L E' in line:
            # '*      D E S I G N   C Y C L E      11    *'
            if idesign != -1:
                all_results.append(optimization_result)
                optimization_result = OptimizationResult()
            print_line = False
            sline = line.split('*')
            assert len(sline) == 3, sline
            design_cycle_n = sline[1].strip()
            design_cycle = int(design_cycle_n.split('D E S I G N   C Y C L E')[1])
            log.warning(f'design cycle={design_cycle:d} idesign={idesign}')
            idesign = 0
            i += 2
            line = lines[i].strip()
            del design_cycle_n
        elif line.startswith('*** ') or line.startswith('^^^ ') or 'SIMCENTER NASTRAN' in line or line.startswith('****************'):
            print_line = False
            pass
        else:
            print_line = False
            #log.debug(f'b {i} {line!r}')
        if print_line:
            log.info(f'c {i} {line0!r}')

        assert isinstance(i, int), f'i={i} line={line}'
        i += 1
        if i > nlines - 1:
            break
        try:
            line = lines[i].strip()
        except Exception:
            print(i, nlines, line)
            raise
    log.debug('*'*80)
    if stop_gradients:
        raise RuntimeError(f'stop_gradients = {stop_gradients}')
    #if len(all_results):
        #log.warning(f'ncycles={len(all_results)}; all_results[0]={all_results[0]}')
    #else:
        #log.warning(f'ncycles={len(all_results)}; all_results={all_results}')

    assert isinstance(i, int), line
    #print('idesign', idesign)
    #asdf
    return i, end_of_job, idesign

def _read_design_optimization_tools(i, line, lines, nlines, design_vars, log):
    #                     DESIGN OPTIMIZATION TOOLS
    #
    #
    #
    #                       (C) COPYRIGHT, 2000
    #
    #                         VANDERPLAATS R&D
    #
    #                  ALL RIGHTS RESERVED, WORLDWIDE
    #
    #                           VERSION 5.3
    #
    #
    #
    #
    #
    #  CONTROL PARAMETERS
    #
    #  OPTIMIZATION METHOD,                METHOD =       1
    #  NUMBER OF DECISION VARIABLES,          NDV =      17
    #  NUMBER OF CONSTRAINTS,                NCON =     314
    #  PRINT CONTROL PARAMETER,            IPRINT =       7
    #  GRADIENT PARAMETER,                  IGRAD =       1
    #    GRADIENTS ARE SUPPLIED BY THE USER
    #  THE OBJECTIVE FUNCTION WILL BE MINIMIZED
    #
    #
    #  -- SCALAR PROGRAM PARAMETERS
    #
    #  REAL PARAMETERS
    #    1) CT     = -3.00000E-02            8) DX2    =  1.47270E-01
    #    2) CTMIN  =  3.00000E-03            9) FDCH   =  1.00000E-03
    #    3) DABOBJ =  2.93242E+01           10) FDCHM  =  1.00000E-04
    #    4) DELOBJ =  1.00000E-03           11) RMVLMZ =  4.00000E-01
    #    5) DOBJ1  =  1.00000E-01           12) DABSTR =  2.93242E+01
    #    6) DOBJ2  =  5.86485E+04           13) DELSTR =  1.00000E-03
    #    7) DX1    =  1.00000E-02
    #
    #  INTEGER PARAMETERS
    #    1) IGRAD  =      1    6) NGMAX  =     34   11) IPRNT1 =      1
    #    2) ISCAL  =     17    7) IGMAX  =      0   12) IPRNT2 =      2
    #    3) ITMAX  =    100    8) JTMAX  =     20   13) JWRITE =      0
    #    4) ITRMOP =      2    9) ITRMST =      2   14) MAXINT =  2147418112
    #    5) IWRITE =      6   10) JPRINT =      0   15) NSTORE = 116226
    #
    #
    #     STORAGE REQUIREMENTS
    #  ARRAY  DIMENSION    MINIMUM    DESIRED    MAXIMUM       USED
    #    WK  106587806       2517       5804     118644     118644
    #   IWK        717        716                              716
    #
    #
    i += 51
    line = lines[i].strip()
    end_of_job = False

    gradient_words = [
        #'CONSTRAINT VALUES (G-VECTOR)',
        #'DESIGN VARIABLE SCALE FACTORS',
        #'DG(I)/DALPHA VECTOR',
        #'SEARCH DIRECTION (S-VECTOR)',
        #'GRADIENT OF THE OBJECTIVE FUNCTION (DF-VECTOR)',
        #'DECISION VARIABLES (X-VECTOR)',
        #'LOWER BOUNDS ON THE DECISION VARIABLES (XL-VECTOR)',
        #'UPPER BOUNDS ON THE DECISION VARIABLES (XU-VECTOR)',
    ]

    while i < nlines:
        if line == '':
            pass
        elif '-- INITIAL VARIABLES AND BOUNDS' in line:
            i += 2
            line = lines[i].strip()

            if line == 'LOWER BOUNDS ON THE DECISION VARIABLES (XL-VECTOR)':
                i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
                i += 1
                line = lines[i].strip()

            if line == 'DECISION VARIABLES (X-VECTOR)':
                i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
                i += 1
                line = lines[i].strip()

            if line == 'UPPER BOUNDS ON THE DECISION VARIABLES (XU-VECTOR)':
                i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
                i += 1
                line = lines[i].strip()

            log.debug(line)
            #sss
        #elif line in gradient_words:
            ##i -= 1
            #i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
            #gradients[line].append(grad)
        elif '-- INITIAL FUNCTION VALUES' in line:
            i += 2
            line = lines[i].strip()
            obj = _get_last_float(line) #  OBJ

            i += 2
            line = lines[i].strip()

            if line == 'CONSTRAINT VALUES (G-VECTOR)':
                i, gvector = _read_gradient(i, line, lines, nlines, log, debug=False)
                i += 1
                line = lines[i].strip()
            log.debug(line)
        elif '***** WARNING - ONE OR MORE CONSTRAINT VALUES EXCEEDS 10.0 IN MAGNITUDE.' in line:
            log.warning(line)
            #print('keys=', list(gradients.keys()))
            #gvector = gradients['CONSTRAINT VALUES (G-VECTOR)'][-1]
            #constraint_ids = gradients['CONSTRAINT NUMBERS']
            isort = np.argsort(gvector)
            #print('constraint_ids =', constraint_ids)
            log.info(f'gvector = {gvector[isort].max()} {gvector[isort].min()}')
            labels = design_vars['label']
            #print('isort =', isort)
            #print(design_vars)
            #max_constraint = constraint_ids[isort]
            #print('max_constraint =', max_constraint)

            i += 4
            line = lines[i].strip()

        elif '-- BEGIN CONSTRAINED OPTIMIZATION: MFD METHOD' in line:
            i -= 1
            break
        else:
            #print(line)
            log.warning(f'{i} {line}')
            #design_optimization_tools
        i += 1
        line = lines[i].strip()
    return i

def _read_summary_of_design_cycle_history(i, line, lines, nlines, log):
    i += 1
    line = lines[i].strip()
    unused_convergence_type = None
    unused_nanalyses = None
    end_of_job = False
    while i < nlines:
        if line == '':
            pass
        elif line in ['***************************************************************']:
            pass
        elif line == '(HARD CONVERGENCE ACHIEVED)':
            unused_convergence_type = 'HARD'
        elif line.startswith('NUMBER OF FINITE ELEMENT ANALYSES COMPLETED'):
            unused_nanalyses = _get_last_int(line)
        elif line.startswith('NUMBER OF OPTIMIZATIONS W.R.T. APPROXIMATE MODELS'):
            unused_napproximate_optimizations = _get_last_int(line)
        elif 'OBJECTIVE AND MAXIMUM CONSTRAINT HISTORY' in line:
            i, line, lines2 = _read_objective_and_maximum_constraint_history(i, lines)
        elif 'DESIGN VARIABLE HISTORY' in line:
            i, line, lines2 = _read_design_variable_history(i, lines)
            #for line in lines2:
                #log.debug(line.rstrip())
        elif _is_end_of_job(line):
            end_of_job = True
            break
        elif _is_skip_one_line(line):
            i += 1
        elif _is_info_msg(line):
            pass
        #elif is_page_skip(line):
            #i = goto_page(i, lines)
            #line = lines[i].strip()
        #else:
            #log.debug(f'summary {i} {line}')
            #aaa
        i += 1
        line = lines[i].strip()
    return i, end_of_job

def _get_last_int(line):
    return int(line.split()[-1])

def _get_last_float(line):
    return float(line.split()[-1])

def show_design_vars(design_vars, log: SimpleLogger):
    for key, values in sorted(design_vars.items()):
        log.info(f'{key}: {values}')

def _complete_optimization(i, lines, nlines, log):
    i += 1
    line = lines[i].strip()
    while i < nlines:
        if line == '':
            pass
        elif '-- OPTIMIZATION RESULTS':
            i -= 1
            break
        elif 'NUMBER OF ITERATIONS =' in line:
            pass
        elif 'CONSTRAINT TOLERANCE, CT =' in line:
            pass
        elif 'TERMINATION CRITERIA' in line:
            pass
        elif 'RELATIVE CONVERGENCE CRITERION WAS MET FOR' in line and  'CONSECUTIVE ITERATIONS' in line:
            #RELATIVE CONVERGENCE CRITERION WAS MET FOR  2 CONSECUTIVE ITERATIONS
            pass
        elif 'THERE ARE' in line and 'ACTIVE CONSTRAINTS AND' in line and 'VIOLATED CONSTRAINTS' in line:
            #THERE ARE      11 ACTIVE CONSTRAINTS AND       0 VIOLATED CONSTRAINTS
            pass
        elif 'CONSTRAINT NUMBERS' in line:
            log.debug(line)
            i -= 1
            i, lines2 = _read_gradient(i, line, lines, nlines, log, debug=False)
            constraint_ids = _parse_ints(lines2)
            log.info(f'  constraint_ids={constraint_ids}; n={len(constraint_ids)}')
            i += 1
            line = lines[i].strip()
            log.debug(f'{i} {line}')
        elif 'VARIABLE NUMBERS (MINUS INDICATES LOWER BOUND)' in line:
            i -= 1
            i, lines2 = _read_gradient(i, line, lines, nlines, log, debug=False)
            ids = _parse_ints(lines2)
            log.info(f'ids={ids}')
        else:
            log.debug(f'{i} {line}')
            aaa
        i += 1
        line = lines[i].strip()
    return i

def _is_page_skip(line):
    is_page_skip = (
        'News file' in line or
        'MAXIMUM  DISPLACEMENTS' in line or
        'OLOAD    RESULTANT' in line or
        'MAXIMUM  APPLIED LOADS' in line or

        'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y' in line or
        'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O' in line or
        'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O' in line or
        'C A S E    C O N T R O L    E C H O' in line or
        'G R I D   P O I N T   S I N G U L A R I T Y   T A B L E' in line or

        #'     TOLERANCE LIMITS ARE:  ' in line or
        'TOLERANCE LIMITS ARE:  SA =  30.00, IA(MIN) =  30.00, IA(MAX) = 150.00, WF =   0.05, TR =   0.50 (FLAG = LIMIT VIOLATED)' in line or
        'TOLERANCE LIMITS ARE:     SKEW  =  10.00, IA(MAX) = 160.00 (++++ = LIMIT VIOLATED)' in line or

        #'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (READ MODULE)' in line or
        'R E A L   E I G E N V A L U E S' in line or

        'D I S P L A C E M E N T   V E C T O R' in line or
        'E L E M E N T   S T R A I N   E N E R G I E S' in line or

        'FLUTTER  SUMMARY' in line or
        'A ZERO FREQUENCY ROOT HAS EMERGED.  WHEN THE MACH NO., DENSITY AND VELOCITY ARE COMPATIBLE' in line or
        'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' in line or
        'N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S' in line or
        'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S' in line or
        'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' in line or
        'A E R O D Y N A M I C   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' in line or

        'THIS PROGRAM IS CONFIDENTIAL AND A TRADE SECRET OF MSC.SOFTWARE CORPORATION.  THE RECEIPT OR' in line or
        'FINITE ELEMENT GEOMETRY CHECK RESULTS EXCEED TOLERANCE LEVELS FOR THE FOLLOWING ELEMENTS.' in line or
        'INTERMEDIATE MATRIX ... HP' in line
        #'* * * *  D B D I C T   P R I N T  * * * *' in line
    )
    if line.strip() == '':
        assert is_page_skip is False
    return is_page_skip

def _is_skip_one_line(line):
    return 'DATA BLOCK ' in line and ' WRITTEN ON FORTRAN UNIT ' in line and ', TRL =' in line

def _is_info_msg(line):
    is_info_msg = (
        'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS' in line or
        'PLACED AFTER THE VALUE METRIC FOR THE TEST.' in line or
        'QL HOUSEHOLDER METHOD IS AUTOMATICALLY SELECTED .' in line or
        '(NUMBER OF FORTRAN RECORDS WRITTEN =' in line or
        '(TOTAL DATA WRITTEN FOR DATA BLOCK =' in line or
        '(MAXIMUM SIZE OF FORTRAN RECORDS WRITTEN =' in line or
        '(MAXIMUM POSSIBLE FORTRAN RECORD SIZE =' in  line or
        '^^^ NORMAL MODES ANALYSIS COMPLETED FOR ANALYSIS=' in line or
        '*** USER WARNING MESSAGE 324 (XSORSO)' in line or
        '*** USER WARNING MESSAGE 2251 (IFS1P)' in line or
        '*** USER INFORMATION MESSAGE 4109 (OUTPX2)' in line or
        '*** USER INFORMATION MESSAGE 4114 (OUTPX2)' in line or
        '*** USER INFORMATION MESSAGE 4507 (GPFEV1)'  in line or
        '*** USER INFORMATION MESSAGE 5458 (REIG)' in line or
        '*** USER INFORMATION MESSAGE 6464 (DOM12E)' in line or
        '*** USER WARNING MESSAGE 6654 (DOM9PRD)' in line or
        '*** USER INFORMATION MESSAGE 7555 (GMTSTD)' in line or
        '^^^ USER INFORMATION MESSAGE 9024 (FEA)' in line or
        '^^^ USER INFORMATION MESSAGE 9051 (FEA)' in line or
        '^^^ USER INFORMATION MESSAGE 9052 (FEA)' in line or
        '^^^ STIFFNESS, MASS, DAMPING, AND LOAD GENERATION   INITIATED. DESIGN CYCLE NUMBER=' in line or
        '^^^ STATIC AEROELASTIC ANALYSIS INITIATED.  DESIGN CYCLE NUMBER=' in  line or
        '^^^ STATIC AEROELASTIC ANALYSIS COMPLETED.  DESIGN CYCLE NUMBER=' in line or
        '^^^ STATIC AEROELASTIC ANALYSIS - DATA RECOVERY INITIATED.  DESIGN CYCLE NUMBER=' in line or
        '^^^ FLUTTER DATA RECOVERY COMPLETED.    DESIGN CYCLE NUMBER=' in line or

        '^^^ FLUTTER ANALYSIS INITIATED. DESIGN CYCLE NUMBER=' in line or
        '^^^ FLUTTER ANALYSIS COMPLETED. DESIGN CYCLE NUMBER=' in line or
        '^^^ NORMAL MODES ANALYSIS - DATA RECOVERY COMPLETED FOR ANALYSIS=  FLUTTER  AND DESIGN CYCLE NUMBER=' in line or

        #THERE ARE    829 TRIA3    ELEMENTS HAVING STRAIN ENERGY WHICH IS LESS THAN   0.0010 PERCENT
        'ELEMENTS HAVING STRAIN ENERGY WHICH IS LESS THAN' in line or
        'OF THE TOTAL STRAIN ENERGY OF ALL ELEMENTS.' in line or
        'ONE OR MORE MAT1 ENTRIES HAVE UNREASONABLE OR INCONSISTENT VALUES OF E,G OR NU.  ID OF FIRST ONE =' in line or
        'INTERMEDIATE MATRIX ... ' in line or
        'MSC.NASTRAN' in line or
        # 'DATA BLOCK R1TABRG  WRITTEN ON FORTRAN UNIT 12, TRL ='
        'DATA BLOCK ' in line and ' WRITTEN ON FORTRAN UNIT ' in line and ', TRL =' in line or
        # (TOTAL DATA WRITTEN FOR EOF MARKER =      1 WORDS.)
        '(TOTAL DATA WRITTEN FOR EOF MARKER =' in line or
        'ELEMENTS. FIRST EID =' in line or
        'RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN SUPERELEMENT BASIC SYSTEM COORDINATES.' in line or
        '* * * *  D B D I C T   P R I N T  * * * *' in line or
        line == '*'
    )
    if is_info_msg:
        if 'REFERENCE POINT =' in line:
            raise RuntimeError(line)
        elif '-----   DESIGN VARIABLES   -----' in line:
            raise RuntimeError(line)
    return is_info_msg

def _find_next_table(i, line, lines):
    ii = 0
    is_broken = False
    while 'DESIGN OPTIMIZATION TOOLS' not in line:
        line = lines[i].strip()
        #print(i, line)
        if '-----   DESIGN OBJECTIVE   -----' in line:
            is_broken = True
            break
        elif '-----   DESIGN VARIABLES   -----' in line:
            is_broken = True
            break
        if '-----   DESIGNED PROPERTIES   -----' in line:
            is_broken = True
            break
        elif '-----   DESIGN CONSTRAINTS ON RESPONSES   -----' in line:
            is_broken = True
            break

        #if ii > 1000:
            #aaab
        i += 1
        ii += 1
    return i, is_broken

def _read_property_comparison_table(i, lines, log):
    #                   -----   COMPARISON BETWEEN INPUT PROPERTY VALUES FROM ANALYSIS AND DESIGN MODELS -----
    #
    #      -----------------------------------------------------------------------------------------------------------------------------
    #       PROPERTY   PROPERTY   PROPERTY     ANALYSIS         DESIGN          LOWER           UPPER       DIFFERENCE     SPAWNING
    #         TYPE        ID        NAME         VALUE          VALUE           BOUND           BOUND          FLAG          FLAG
    #      -----------------------------------------------------------------------------------------------------------------------------
    #       PBEAM          1001       A(A)    1.000000E+00    1.000000E+00    1.000000E-15    1.000000E+20      NONE           *
    #       PBEAM          1001      I1(A)    3.333333E-01    3.333333E-01    1.000000E-15    1.000000E+20      NONE           *
    #       PBEAM          1001      I2(A)    2.083333E-02    2.083333E-02    1.000000E-15    1.000000E+20      NONE           *
    #       PBEAML         1474    DIM1(A)    5.000000E-01    5.000000E-01    1.000000E-02    2.000000E+00      NONE
    i += 2
    line = lines[i].strip()
    log.debug(f' {i} {line}')
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    #print('*', i, line)
    return i, line, lines2

def _read_property_table(i, lines):
    # ---------------------------------------------------------------------------------------------------------------
    #     ELEMENT    ELEMENT   CONNECTIVITY  ANALYSIS         DESIGN          LOWER           UPPER       DIFFERENCE
    #         TYPE        ID        NAME         VALUE          VALUE           BOUND           BOUND          FLAG
    # ---------------------------------------------------------------------------------------------------------------
    #     CONM2          1982          M    1.294000E-03    1.294000E-03    2.000000E-04    2.588000E-03      NONE
    #     CONM2          1983          M    1.294000E-03    1.294000E-03    2.000000E-04    2.588000E-03      NONE
    #     CONM2          1984          M    1.294000E-03    1.294000E-03    2.000000E-04    2.588000E-03      NONE
    i += 2
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_design_constraints_of_responses(i, lines):
    #                                          -----   DESIGN CONSTRAINTS ON RESPONSES   -----
    #
    #                                           (MAXIMUM RESPONSE CONSTRAINTS MARKED WITH **)
    #
    #     ---------------------------------------------------------------------------------------------------------
    #                                INTERNAL                            INTERNAL
    #      INTERNAL      DCONSTR     RESPONSE     RESPONSE      L/U       REGION       SUBCASE
    #         ID           ID           ID          TYPE        FLAG        ID            ID           VALUE
    #     ---------------------------------------------------------------------------------------------------------
    #             1         1101            2     STRESS       UPPER            91           1      -4.6692E-01
    i += 9
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_constraints_on_designed_properties(i, lines):
    #                                         -----   CONSTRAINTS ON DESIGNED PROPERTIES   -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #      INTERNAL     PROPERTY        PROPERTY           L/U              CYCLE                 INPUT                 OUTPUT
    #         ID           ID             NAME             FLAG             LIMIT                 VALUE                  VALUE
    #     ---------------------------------------------------------------------------------------------------------------------------
    #           281            2               T          LOWER            5.4736E-01           -2.5000E-01           -1.4618E-01
    i += 6
    line = lines[i].strip()
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_weight_as_a_function_of_material_id(i, lines):
    #                                        -----   WEIGHT AS A FUNCTION OF MATERIAL ID ----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                                    MATERIAL       SEID          INPUT            OUTPUT
    #                                       ID                       WEIGHT            WEIGHT
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                                           2           0       4.2986E+04            N/A
    #                                           3           0       0.0000E+00            N/A
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                                      TOTAL                    4.2986E+04            N/A
    i += 6
    line = lines[i].strip()
    #print(line)
    lines2 = []
    while 'TOTAL' not in line:
        lines2.append(line)
        i += 1
        line = lines[i].rstrip()
    lines2.pop()
    lines2.append(line)

    #i, lines2 = _read_line_block(i, lines, strip=True)
    #line = lines[i].strip()
    #for linei in lines2:
        #print(linei)
    return i, line, lines2

def _read_designed_properties(i, lines):
    #                                                -----   DESIGNED PROPERTIES   -----
    #
    #     ----------------------------------------------------------------------------------------------------------
    #         PROPERTY     PROPERTY     PROPERTY     TYPE OF          LOWER                             UPPER
    #           TYPE          ID          NAME       PROPERTY         BOUND            VALUE            BOUND
    #     ----------------------------------------------------------------------------------------------------------
    #         PSHELL           1464            T     DVPREL1           N/A           5.0000E-02          N/A
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_stress_responses(i, lines):
    #                                                 -----    STRESS RESPONSES    -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #        INTERNAL   DRESP1   RESPONSE   ELEMENT    VIEW    COMPONENT      LOWER         INPUT        OUTPUT         UPPER
    #           ID        ID      LABEL        ID     ELM ID      NO.         BOUND         VALUE         VALUE         BOUND
    #     ---------------------------------------------------------------------------------------------------------------------------
    #               2     10036  U4446         5644                   9       N/A        2.6654E+04    2.5291E+04    5.0000E+04
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_composite_stress_responses(i, lines):
    #                                        -----    COMPOSITE LAMINATE STRESS RESPONSES    -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #        INTERNAL   DRESP1   RESPONSE   ELEMENT  COMPONENT   LAMINA       LOWER         INPUT        OUTPUT         UPPER
    #           ID        ID      LABEL        ID       NO.        NO.        BOUND         VALUE         VALUE         BOUND
    #     ---------------------------------------------------------------------------------------------------------------------------
    #             487     10198  RK1200       13410         11         2       N/A        1.5603E+04    1.5609E+04    2.2500E+04
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_flutter_responses(i, lines):
    #                                                 -----    FLUTTER RESPONSES    -----
    #
    #   -----------------------------------------------------------------------------------------------------------------------------
    #  INTERNAL  DRESP1 RESPONSE    MODE                    MACH                      LOWER        INPUT       OUTPUT        UPPER
    #     ID        ID    LABEL     NO.       DENSITY        NO.       VELOCITY       BOUND        VALUE        VALUE        BOUND
    #   -----------------------------------------------------------------------------------------------------------------------------
    #     1575       44 D350            4   6.5350E-08   7.5000E-01   9.3957E+03      N/A       2.5826E-01   2.6043E-01      N/A
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_weight_responses(i, lines):
    #                                                  -----    WEIGHT RESPONSE    -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #          INTERNAL    DRESP1    RESPONSE     ROW       COLUMN         LOWER          INPUT         OUTPUT          UPPER
    #             ID         ID       LABEL        ID         ID           BOUND          VALUE          VALUE          BOUND
    #     ---------------------------------------------------------------------------------------------------------------------------
    #               1      1000      W             3          3              N/A        2.6865E+05    2.6862E+05       N/A
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_retained_dresp2_responses(i, lines):
    #                                                 ---- RETAINED DRESP2 RESPONSES ----
    #
    #     ----------------------------------------------------------------------------------------------------------
    #         INTERNAL      DRESP2      RESPONSE     EQUATION         LOWER                             UPPER
    #            ID           ID         LABEL          ID            BOUND            VALUE            BOUND
    #     ----------------------------------------------------------------------------------------------------------
    #                1        10200     GRAT400            42          N/A           1.9672E-01       2.0000E-01
    i += 6
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_design_variable_history(i, lines):
    #                                                       DESIGN VARIABLE HISTORY
    # ----------------------------------------------------------------------------------------------------------------------------------
    #  INTERNAL |   EXTERNAL   |             |
    #   DV. ID. |    DV. ID.   |    LABEL    |   INITIAL    :      1       :      2       :      3       :      4       :      5       :
    # ----------------------------------------------------------------------------------------------------------------------------------
    #         1 |        900   |  T900       |   6.3111E-02 :   6.1672E-02 :   6.1672E-02 :   6.1256E-02 :   6.1177E-02 :
    i += 5
    line = lines[i].strip()
    #print(i, line)
    i, lines2 = _read_line_block(i, lines, strip=True)
    line = lines[i].strip()
    return i, line, lines2

def _read_objective_and_maximum_constraint_history(i, lines):
    #                                              OBJECTIVE AND MAXIMUM CONSTRAINT HISTORY
    #           ---------------------------------------------------------------------------------------------------------------
    #                              OBJECTIVE FROM           OBJECTIVE FROM          FRACTIONAL ERROR          MAXIMUM VALUE
    #             CYCLE              APPROXIMATE                 EXACT                    OF                       OF
    #             NUMBER            OPTIMIZATION               ANALYSIS              APPROXIMATION             CONSTRAINT
    #           ---------------------------------------------------------------------------------------------------------------
    #
    #             INITIAL                                      2.686645E+05                                      4.002925E-02
    #
    #                   1             2.686527E+05             2.686527E+05            -2.326423E-07             1.836494E-03
    #
    #                   2             2.686161E+05             2.686161E+05             0.000000E+00             2.400868E-03
    #
    #                   3             2.686068E+05             2.686068E+05             0.000000E+00             2.770226E-03
    #
    #                   4             2.686055E+05             2.686055E+05             0.000000E+00             2.677865E-03
    #           ---------------------------------------------------------------------------------------------------------------
    i += 7
    lines2 = []
    i = _goto_page(i, lines)
    line = lines[i].strip()
    return i, line, lines2

def _read_design_objective(i: int, lines: list[str],
                           optimization_result: OptimizationResult, log: SimpleLogger) -> tuple[int, str]:
    #                                                  -----   DESIGN OBJECTIVE   -----
    #
    #     ----------------------------------------------------------------------------------------------------
    #                INTERNAL      TYPE                  MINIMIZE
    #                RESPONSE       OF                      OR       SUPERELEMENT    SUBCASE
    #                   ID       RESPONSE      LABEL     MAXIMIZE         ID            ID          VALUE
    #     ----------------------------------------------------------------------------------------------------
    #                       1     DRESP1     WEIGHT      MINIMIZE             0            0      2.6861E+05
    #
    #                                                  -----   DESIGN OBJECTIVE   -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                INTERNAL      TYPE                  MINIMIZE
    #                RESPONSE       OF                      OR       SUPERELEMENT    SUBCASE        INPUT          OUTPUT
    #                   ID       RESPONSE      LABEL     MAXIMIZE         ID            ID          VALUE           VALUE
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                       4     DRESP2        N/A      MINIMIZE         N/A          N/A        1.8536E+01      2.3529E+01
    log.warning('design objective')
    objective = optimization_result.design_objective
    i += 7
    line = lines[i].rstrip()
    lines2 = []
    while line != '':
        lines2.append(line)
        i += 1
        line = lines[i].rstrip()
    assert len(lines2), lines2
    #print(lines2)

    for linei in lines2:
        sline = linei.split()
        if len(sline) == 7:
            internal_respose_id, response_type, label, min_max, superelement_id, subcase_id, input_value = sline
            input_value = float(input_value)
            output_value = input_value
        elif len(sline) == 8:
            internal_respose_id, response_type, label, min_max, superelement_id, subcase_id, input_value, output_value = sline
            input_value = float(input_value)
            output_value = float(output_value)
        else:
            print(sline, len(sline))
            aaa
        if superelement_id == 'N/A':
            superelement_id = 0
    objective['subcase_id'].append(subcase_id)
    objective['label'].append(label)
    objective['in'].append(input_value)
    objective['out'].append(output_value)
    return i, line

def _read_design_variables(i: int, lines: list[str],
                           design_cycle: int, idesign: int, optimization_result: OptimizationResult,
                           log: SimpleLogger) -> tuple[int, str]:
    """reads the design variables"""
    #log.info(f'_read_design_variables; idesign={idesign:d}')
    #design_vars = optimization_result.design_vars
    #                                                  -----   DESIGN VARIABLES   -----
    #
    #     ---------------------------------------------------------------------------------------------------------
    #            INTERNAL       DESVAR                         LOWER                               UPPER
    #               ID            ID          LABEL            BOUND             VALUE             BOUND
    #     ---------------------------------------------------------------------------------------------------------
    #                   1           900      T900            5.0000E-02        6.3111E-02        6.0000E+00
    #
    #                                                  -----   DESIGN VARIABLES   -----
    #
    #     ---------------------------------------------------------------------------------------------------------------------------
    #            INTERNAL       DESVAR                         LOWER             INPUT            OUTPUT             UPPER
    #               ID            ID          LABEL            BOUND             VALUE             VALUE             BOUND
    #     ---------------------------------------------------------------------------------------------------------------------------
    #                   1          1001      DV_1001         1.0000E-02        5.0000E-01        5.0000E-01        2.0000E+00
    i += 6
    line = lines[i].rstrip()
    lines2 = []
    while line != '':
        if line.startswith((' ^^^', ' ***')):
            i += 1
            line = lines[i].rstrip()
            if 'PAGE' in line:
                break
            continue
        lines2.append(line)
        i += 1
        line = lines[i].rstrip()
        if 'PAGE' in line:
            break

    assert len(lines2), lines2
    #print(lines2)
    nsline = len(lines2[0].split())

    if nsline == 6:
        for linei in lines2:
            if 'PAGE' not in line:
                log.debug(line)
            sline = linei.split()
            try:
                internal_id, desvar_id, label, xl, x, xu = sline
            except ValueError:
                print(line, sline)
                continue
                #raise

            #  INTERNAL       DESVAR                         LOWER                               UPPER
            #     ID            ID          LABEL            BOUND             VALUE             BOUND
            #                                                  xlb               x                xub
            xi = x
            internal_id = int(internal_id)
            desvar_id = int(desvar_id)
            xl = float(xl)
            xi = float(xi)
            x = float(x)
            xu = float(xu)
            optimization_result.design_vars['internal_id'].append(internal_id)
            optimization_result.design_vars['desvar_id'].append(desvar_id)
            optimization_result.design_vars['label'].append(label)
            optimization_result.design_vars['xl'].append(xl)
            optimization_result.design_vars['xi'].append(xi)
            optimization_result.design_vars['x'].append(x)
            optimization_result.design_vars['xu'].append(xu)
    elif nsline == 7:
        for linei in lines2:
            sline = linei.split()
            #  INTERNAL       DESVAR                         LOWER             INPUT            OUTPUT             UPPER
            #     ID            ID          LABEL            BOUND             VALUE             VALUE             BOUND
            #                                                 xlb               xi                x                 xub
            internal_id, desvar_id, label, xl, xi, x, xu = sline
            #print(f'sline7: internal_id={internal_id}')
            internal_id = int(internal_id)
            desvar_id = int(desvar_id)
            xl = float(xl)
            xi = float(xi)
            x = float(x)
            xu = float(xu)
            optimization_result.design_vars['internal_id'].append(internal_id)
            optimization_result.design_vars['desvar_id'].append(desvar_id)
            optimization_result.design_vars['label'].append(label)
            optimization_result.design_vars['xl'].append(xl)
            optimization_result.design_vars['xi'].append(xi)
            optimization_result.design_vars['x'].append(x)
            optimization_result.design_vars['xu'].append(xu)
            if internal_id == 1:
                print(f'design_cycle={design_cycle} idesign={idesign} x={x} xi={xi} xl={xl} xu={xu}')
                print(f'design_vars={optimization_result.design_vars}')
            #print()
            #if internal_id > 100:
                #print(design_vars)
                #asfd
    else:
        raise RuntimeError('line=%s' % lines2[0])
    assert isinstance(i, int), line
    line = lines[i].strip()
    #if len(design_vars['xu']) > 40:
        #asdf
    #design_vars['internal_id'] = np.array(design_vars['internal_id'], dtype='int32')
    #design_vars['desvar_id'] = np.array(design_vars['desvar_id'], dtype='int32')

    #design_vars['label'] = np.array(design_vars['label'], dtype='|U8')
    #design_vars['x'] = np.array(design_vars['x'], dtype='float32')
    #design_vars['xi'] = np.array(design_vars['xi'], dtype='float32')
    #design_vars['xu'] = np.array(design_vars['xu'], dtype='float32')
    #design_vars['xl'] = np.array(design_vars['xi'], dtype='float32')
    #print('desvars', i)
    return i, line

def _read_iteration_number(i: int, line: str, lines: list[str],
                           design_vars, constraint_ids_failing, gvector,
                           nlines: int, log: SimpleLogger):
    #i += 1
    line = lines[i].strip()
    iteration = int(line.split()[-1])
    log.debug('iteration=%s' % iteration)

    i += 1
    line = lines[i].strip()
    # -----------------------
    scaling_info_words = [
        'DESIGN VARIABLE SCALE FACTORS',
    ]
    scaling_info = {}

    while i < nlines:
        line = lines[i].strip()
        if 'BIOMSG: ERROR' in line:
            labels = design_vars['label']
            internal_id = design_vars['internal_id']
            print('constraint_ids_failing', constraint_ids_failing, len(constraint_ids_failing))
            print('constraint_ids', constraint_ids, len(constraint_ids))
            print('labels', len(labels))
            print('internal_id', internal_id, len(internal_id))
            print('gvector', gvector, len(gvector))
            asdf
            break
        if line == '':
            pass
        elif 'SCALING INFORMATION' in line:
            pass
        elif '***** THE NUMBER OF ACTIVE/VIOLATED' in line:
            i += 6
            continue
        elif 'AT THIS POINT, DOT NEEDS AT LEAST' in line:
            i += 2
            continue
        elif '* * * BEFORE MODIFICATION * * *' in line:
            i += 2
            line = lines[i].strip()
            assert 'ACTIVE CONSTRAINTS' in line and 'VIOLATED CONSTRAINTS' in line
            # THERE ARE      13 ACTIVE CONSTRAINTS AND      14 VIOLATED CONSTRAINTS
            i += 1
            line = lines[i].strip()
            assert 'CONSTRAINT NUMBERS' in line, line
            i, constraint_ids = _read_int_gradient(i, line, lines, nlines, log, debug=False)
            log.info(f'  constraint_ids={constraint_ids}; n={len(constraint_ids)}')
            i += 1
            line = lines[i].strip()
            log.debug(f'{i} {line}')
            #if 'BIOMSG: ERROR' in line:
                #break
            assert '* * * AFTER MODIFICATION * * *' in line, line

        elif line in scaling_info_words:
            assert line not in scaling_info, line
            i -= 1
            i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
            scaling_info[line] = grad
        elif 'ACTIVE CONSTRAINTS' in line and 'VIOLATED CONSTRAINTS' in line:
            # THERE ARE      13 ACTIVE CONSTRAINTS AND      14 VIOLATED CONSTRAINTS
            i += 1
            line = lines[i].strip()
            #print('*******', line)
            assert 'CONSTRAINT NUMBERS' in line, line
            #i -= 1
            #i += 1
            #line = lines[i].strip()
            #assert 'CONSTRAINT NUMBERS' in line, line
            i, constraint_ids = _read_int_gradient(i, line, lines, nlines, log, debug=False)
            log.debug(f'  constraint_ids={constraint_ids}; n={len(constraint_ids)}')
            scaling_info['CONSTRAINT NUMBERS'] = constraint_ids
            i += 1
            line = lines[i].strip()
            #print(i, line)
            if 'BIOMSG: ERROR' in line:
                break
            assert 'GRADIENT OF THE OBJECTIVE FUNCTION (DF-VECTOR)' == line, line
            i -= 1
            i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
            i += 1
            line = lines[i].strip()
            #print(i, line)

            for idi in constraint_ids:
                #print(i, line, idi)
                assert 'GRADIENT OF CONSTRAINT NUMBER' in line
                igradient = int(line.split()[-1])
                assert idi == igradient, idi
                i -= 1
                i, grad = _read_gradient(i, line, lines, nlines, log, debug=False)
                i += 1
                line = lines[i].strip()
            i -= 1
            line = lines[i].strip()
            #print('finished gradients', i, line)

            if 'THE GRADIENT OF CONSTRAINT NUMBERS' in lines[i+1]:
                while 'THE GRADIENT OF CONSTRAINT NUMBERS' in lines[i+1]:
                    #THE GRADIENT OF CONSTRAINT NUMBERS      4779 AND      4778 ARE DEPENDENT
                    #CONSTRAINT NUMBER      4779 IS REMOVED FROM THE ACTIVE SET
                    i += 3
            line = lines[i].strip()
            break

        #elif line in int_tables:
            #i, lines2 = _read_gradient(i, line, lines, nlines)
            #ids = _parse_ints(lines2)
            #scaling_info[line] = ids
        else:
            print(line)
            asdf
        i += 1
    #log.debug('end iteration=%s' % iteration)

    scaling_keys = list(scaling_info.keys())
    if scaling_keys:
        log.info(f'scaling_info {scaling_keys}')
    return i

def _read_1d_search(i, line, lines, nlines, log):
    debug = False
    is_constrained_function = False
    print('-'*80)
    log.debug(f'read_1d_search {i} {line}')
    i += 1
    line = lines[i].strip()
    alpha = None
    alpha_max = None
    obj = None
    constraint_variables = [
        'DECISION VARIABLES (X-VECTOR)',
        'CONSTRAINT VALUES (G-VECTOR)',
        'ACTIVE CONSTRAINT VALUES',
        'PERTURBATION OF X-VECTOR',
        'ADJUSTED DECISION VARIABLES (X-VECTOR)',
    ]
    constraints = {}

    while i < nlines:
        line = lines[i].strip()
        #print(i, line)
        if '* * * * * END OF ONE-DIMENSIONAL SEARCH * * * * *' in line:
            break
        elif line == '':
            pass
        #elif 'ACTIVE CONSTRAINT VALUES' in line:

        elif line == 'CONSTRAINED FUNCTION':
            assert is_constrained_function is False
            is_constrained_function = True
        elif 'ALPHA =' in line and 'ALPMAX =' in line:
            # ALPHA =  0.00000E+00    ALPMAX =  2.44288E+01
            assert alpha_max is None, line
            *trash, alpha, alpha_maxstr, eq2, alpha_max = line.split()
            del trash, alpha_maxstr, eq2
            alpha = float(alpha)
            alpha_max = float(alpha_max)
        elif 'OBJ =' in line:
            # OBJ =  2.63781E+01
            #assert obj is None, line
            obj = line.split()[-1]
            obj = float(obj)
        elif line in constraint_variables:
            i -= 1
            i, grad = _read_gradient(i, line, lines, nlines, log, debug=debug)
            constraints[line] = grad
            del grad

        elif 'FIND BRACKETS ON MINIMUM' in line:
            assert is_constrained_function is True
        elif 'ALPHA =' in line and 'PROPOSED' not in line:
            #assert alpha is None, line
            # ALPHA =  2.06505E-04
            alpha = line.split()[-1]
            alpha = float(alpha)
        elif 'BRACKETS, XL =' in line:
            # BRACKETS, XL =  0.00000E+00  XU =  2.06505E-04
            *trash, xl, xustr, eq2, xu = line.split()
            del trash, xustr, eq2
            xl = float(xl)
            xu = float(xu)
        elif '--- RESULTS OF ONE-DIMENSIONAL SEARCH' in line:
            assert is_constrained_function is True
        elif 'OBJECTIVE =' in line:
            # OBJECTIVE =  2.63781E+01
            #assert obj is None, line
            obj = line.split()[-1]
            obj = float(obj)
        else:
            print(line)
            asdf
        i += 1
    print('-'*80)
    return i

def _parse_ints(lines):
    """
      CONSTRAINT NUMBERS
           1         3         6         7         8
    """
    data = []
    for line in lines:
        data.extend(line.split())
    ids = np.array(data, dtype='int32')
    return ids

def _parse_gradient(lines, debug=True):
    # 486)   2.00000E-04   2.00000E-04   2.00000E-04   2.00000E-04   2.00000E-04
    data = []
    iold = 1 - 5
    for line in lines:
        sline = line.split()
        #if debug:
            #print(sline)
        i_paren, *sline = line.split()
        i = int(i_paren[:-1])
        assert i == iold + 5, f'i={i} iold+5={iold+5}'
        data.extend(sline)
        #print('%r' % line)
        iold = i
    #print('***', lines[0])
    #print(data)
    gradient = np.array(data, dtype='float32')
    return gradient

def _get_decision_variables(lines):
    header = lines[0].split()
    #print(header)
    ids = []
    data = []
    for line in lines:
        idi, *datai = line.split()
        ids.append(idi)
        data.append(datai)
    ids = np.array(ids, dtype='int32')
    data = np.array(data, dtype='float32')
    return header, ids, data

def _is_end_of_job(line):
    return '* * * END OF JOB * * *' in line or '* * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *' in line

def read_sol_200(f06_filename: str):
    log = get_logger(log=None, level='debug', encoding='utf-8')

    with open(f06_filename, 'r') as f:
        lines = f.readlines()

    nlines = len(lines)
    i = 0
    line = lines[i].strip()

    all_results = [OptimizationResult()]
    idesign = -1

    while i < nlines:
        #log.debug(f'{i} {line}')
        if line == '*                D E S I G N    O P T I M I Z A T I O N            *':
            i, end_of_job, idesign = _read_design_optimization(i, line, lines, nlines,
                                                               idesign, all_results, log)
            if end_of_job:
                log.debug('end of design optimization')
                break
        #elif line == '*       I N I T I A L   A N A L Y S I S       *':
            #i = read_optimization_cycle(i, line, lines, nlines,log)
        elif  line == 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R':
            i += 31
            line = lines[i].strip()  # page xx
        elif _is_page_skip(line):
            i = _goto_page(i, lines)
            line = lines[i].strip()
        elif _is_skip_one_line(line):
            i += 1
        elif _is_info_msg(line):
            pass
        elif _is_end_of_job(line):
            end_of_job = True
            break

        elif line == '-----   COMPARISON BETWEEN INPUT PROPERTY VALUES FROM ANALYSIS AND DESIGN MODELS -----':
            i, line, lines2 = _read_property_comparison_table(i, lines, log)
            i += 7
        elif line == '-----   COMPARISON BETWEEN INPUT CONNECTIVITY PROPERTY VALUES FROM ANALYSIS AND DESIGN MODELS -----':
            i, line, lines2 = _read_property_table(i, lines)
            log.debug('end of comparison table2')
        elif '*****   OPTIMIZATION RESULTS BASED ON THE APPROXIMATE MODEL   *****' in line:
            asdf
        elif 'LOWER BOUNDS ON THE DECISION VARIABLES' in line:
            asdf
        elif 'GRADIENT OF' in line:
            asdf
        elif 'PBEAML         1477    DIM1(B)' in line:
            aaa
        elif '196)' in line:
            adsf
        elif 'BLOCK SIZE USED ......................' in line:
            i += 6

        elif _is_skip_one_line(line):
            i += 1
        #elif 'F I N A L   A N A L Y S I S' in line:
            #i -= 1
            #log.debug('breaking %s' % line)
            #break

        elif line == '0':
            pass
        elif line:
            pass
            #log.debug(f'a {i} {line}')
        #print(i, line)

        i += 1
        if i > nlines - 1:
            break
        line = lines[i].strip()

    log.warning('finished reading f06 file...')
    return all_results

def plot_sol_200(f06_filename: str, show: bool=True):
    all_results = read_sol_200(f06_filename)

    #--------------------------------------------------------------
    ncycles = len(all_results)
    objs_in = [result.design_objective['in'] for result in all_results]
    objs_out = [result.design_objective['out'] for result in all_results]

    #(48, 1309)
    #(ncycle, ndvar)
    desvarsi = np.vstack([result.design_vars['xi'] for result in all_results])  # rows are design variable
    desvars = np.vstack([result.design_vars['x'] for result in all_results])  # rows are design variable
    #desvar_id = all_results[0].design_vars['desvar_id']
    #print(desvar_id)
    print(desvars.shape)
    labels = all_results[0].design_vars['label']
    #print(labels[0])
    dvi_max = desvarsi.max(axis=0) # len=1309
    dvi_min = desvarsi.min(axis=0) # len=1309
    dmin_maxi = dvi_max - dvi_min
    print(desvars[:, :8])

    dv_max = desvars.max(axis=0) # len=1309
    dv_min = desvars.min(axis=0) # len=1309
    dmin_max = dv_max - dv_min
    #for label, dmin_maxi in zip(labels, dmin_max):
        #print(label, dmin_maxi)
    print(desvars[:, 0])
    #print(desvarsi[:, 0])
    print(np.abs(dmin_max).max())
    #print(dv_max.shape)
    #print(labels)

    if show:
        design_cycle = range(1, ncycles+1)
        import matplotlib.pyplot as plt
        plt.plot(design_cycle, objs_in, 'o-')
        plt.plot(design_cycle, objs_out, 'o-')
        plt.xlabel('Design Cycle')
        plt.ylabel('Objective Function')
        plt.grid()
        plt.show()
    print(all_results[-1])
    return all_results

if __name__ == '__main__':   # pragma: no cover
    bdf_filename = 'optimize_formsc_2.f06'
    plot_sol_200(bdf_filename)

    #bdf_filename = 'model_200.f06'
    #read_sol_200(bdf_filename)

    #bdf_filename = 'crm_200_9.f06'
    #read_sol_200(bdf_filename)
