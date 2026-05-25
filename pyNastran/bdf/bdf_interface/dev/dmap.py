"""
DMAP ALTER utilities for pyNastran.

Provides functions to generate DMAP ALTER cards for exporting Nastran
datablocks (matrices) like KGG, MGG, KDICT, etc.

"""
from __future__ import annotations
from typing import Optional, TYPE_CHECKING
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


# Known matrix datablocks available at various solution stages
STIFFNESS_MATRICES = {'KGG', 'KAA', 'KNN', 'KFF', 'KOO', 'KELM'}
MASS_MATRICES = {'MGG', 'MAA', 'MNN', 'MFF', 'MOO', 'MELM'}
DAMPING_MATRICES = {'BGG', 'BAA', 'BNN', 'BFF', 'BOO', 'BELM'}
LOAD_MATRICES = {'PG', 'PA', 'PL'}
GEOMETRY_TABLES = {'GPDT', 'EQEXIN', 'CSTM', 'GPL'}
DICT_TABLES = {'KDICT', 'MDICT', 'BDICT', 'KDICTP', 'KELMP'}

# DMAP compile targets and alter locations for matrix export.
# Used only for COMPILE/ALTER approach (MSC Nastran style).
_COMPILE_MODULE = {
    101: 'SESTATIC',   # SOL 101 static
    103: 'SEMODES',    # SOL 103 modal
    108: 'SEFREQR',    # SOL 108 frequency response
    111: 'SEMFREQR',   # SOL 111 modal frequency response
    200: 'DESOPT',     # SOL 200 optimization
}


# Predefined MALTER labels in NX Nastran DMAP source.
# Each label is a known insertion point — you cannot create custom ones.
MALTER_LABELS = {
    'KGG': "MALTER:(KGG, BGG, MGG, K4GG, PG)",
    'KAA': "MALTER:AFTER SUPERELEMENT MATRIX AND LOAD REDUCTION TO A-SET,",
    'OUGV1': "MALTER:(OUGV1, OES1, OEF1, ETC.)",
    'OUGV2': "MALTER:(OUGV2, OES2, OEF2, ETC.)",
    'PHASE1_TOP': "MALTER:TOP OF PHASE 1 SUPERELEMENT LOOP AFTER PARAMETERS AND",
    'PHASE1_BOTTOM': "MALTER:BOTTOM OF PHASE 1 SUPERELEMENT LOOP",
    'PREFACE': "MALTER:AFTER PREFACE MODULES",
}


def make_matrix_export_malter(
    matrices: list[str],
    malter_label: str = 'KGG',
    output_unit: int = -12,
) -> list[str]:
    """Generate NX Nastran MALTER lines for matrix export to the .op2 file.

    This is the simplest and most reliable approach for NX Nastran.
    Uses MALTER which inserts OUTPUT2 at a predefined DMAP location.
    Requires PARAM,POST,-1 in bulk data for the .op2 file to be created.

    Parameters
    ----------
    matrices : list[str]
        Matrix names to export (e.g., ['KGG', 'MGG', 'KELM', 'KDICT'])
    malter_label : str
        Predefined MALTER insertion point. Options:
          'KGG' - after global stiffness/mass assembly (KGG, MGG available)
          'KAA' - after A-set reduction (KAA, MAA available)
          'OUGV1' - after data recovery sort1
        Or pass a full custom label string.
    output_unit : int
        Fortran unit for output. -12 means the standard .op2 file (default).
        Negative unit means binary, positive means formatted.

    Returns
    -------
    exec_lines : list[str]
        Executive control lines (MALTER + OUTPUT2 statements)

    """
    if not matrices:
        raise ValueError("matrices list cannot be empty")

    matrices_upper = [m.upper() for m in matrices]

    # resolve label
    if malter_label in MALTER_LABELS:
        label = MALTER_LABELS[malter_label]
    else:
        label = malter_label

    exec_lines = []
    exec_lines.append(f"MALTER '{label}' $")

    for name in matrices_upper:
        exec_lines.append(f'OUTPUT2 {name}//-3/{output_unit}/{name}// $')

    return exec_lines


def add_matrix_export_malter(
    model: BDF,
    matrices: list[str],
    malter_label: str = 'KGG',
    output_unit: int = -12,
) -> None:
    """Add NX Nastran MALTER to a BDF model for matrix export.

    Modifies model.executive_control_lines in-place. Matrices are
    exported to the standard .op2 output file.

    Requires PARAM,POST,-1 in bulk data for the .op2 file to be created.

    Parameters
    ----------
    model : BDF
        The BDF model to modify
    matrices : list[str]
        Matrix names to export (e.g., ['KGG', 'MGG', 'KELM', 'KDICT'])
    malter_label : str
        Predefined MALTER insertion point (see make_matrix_export_malter)
    output_unit : int
        Fortran unit for output. -12 means the standard .op2 file.

    """
    exec_lines = make_matrix_export_malter(matrices, malter_label, output_unit)
    model.executive_control_lines.extend(exec_lines)


def make_matrix_export_alter(
    matrices: list[str],
    output_filename: str = 'matrices.op4',
    output_format: str = 'op4',
    sol: int = 101,
    unit: int = 51,
    form: str = 'FORMATTED',
    compile_module: Optional[str] = None,
    alter_location: Optional[str] = None,
) -> tuple[list[str], list[str]]:
    """Generate DMAP ALTER lines to export matrix datablocks.

    Parameters
    ----------
    matrices : list[str]
        Matrix names to export (e.g., ['KGG', 'MGG', 'KDICT'])
    output_filename : str
        Output filename (e.g., 'kgg.op4' or 'kgg.op2')
    output_format : str
        'op4' for OUTPUT4 or 'op2' for OUTPUT2
    sol : int
        Solution sequence (101, 103, 108)
    unit : int
        Fortran unit number for the output file (default 51)
    form : str
        For op4: 'FORMATTED' (text) or 'UNFORMATTED' (binary)
        For op2: ignored
    compile_module : str, optional
        Override the COMPILE subdmap name (e.g., 'PHASE1A', 'SESTATIC').
        Default is SOL-dependent.
    alter_location : str, optional
        Override the ALTER location string (e.g., 'END', 'BEGIN', '240').
        Default is 'END'.

    Returns
    -------
    system_lines : list[str]
        ASSIGN statements (go before SOL line)
    exec_lines : list[str]
        COMPILE/ALTER/OUTPUT lines (go after SOL line, before CEND)

    """
    if output_format not in ('op4', 'op2'):
        raise ValueError(f"output_format must be 'op4' or 'op2', got {output_format!r}")
    if sol not in _COMPILE_MODULE and compile_module is None:
        raise ValueError(
            f"SOL {sol} not supported for automatic ALTER. "
            f"Supported: {sorted(_COMPILE_MODULE.keys())}")
    if not matrices:
        raise ValueError("matrices list cannot be empty")

    matrices_upper = [m.upper() for m in matrices]
    if compile_module is None:
        compile_module = _COMPILE_MODULE[sol]
    if alter_location is None:
        alter_location = 'END'

    system_lines = []
    exec_lines = []

    if output_format == 'op4':
        system_lines.append(
            f"ASSIGN OUTPUT4='{output_filename}', UNIT={unit}, "
            f"FORM={form}, STATUS=NEW")
        exec_lines.append(f'COMPILE {compile_module}')
        exec_lines.append(f"ALTER '{alter_location}'")

        # OUTPUT4 exports up to 5 matrices per statement
        for i in range(0, len(matrices_upper), 5):
            chunk = matrices_upper[i:i+5]
            # pad to 5 slots with empty
            padded = ','.join(chunk) + ',' * (5 - len(chunk))
            exec_lines.append(f'OUTPUT4 {padded}//0/{unit} $')

    else:  # op2
        system_lines.append(
            f"ASSIGN OUTPUT2='{output_filename}', UNIT={unit}, STATUS=NEW")
        exec_lines.append(f'COMPILE {compile_module}')
        exec_lines.append(f"ALTER '{alter_location}'")

        # OUTPUT2 exports up to 5 matrices per statement
        for i in range(0, len(matrices_upper), 5):
            chunk = matrices_upper[i:i+5]
            padded = ','.join(chunk) + ',' * (5 - len(chunk))
            exec_lines.append(f'OUTPUT2 {padded}//{unit} $')

    return system_lines, exec_lines


def add_matrix_export_alter(
    model: BDF,
    matrices: list[str],
    output_filename: str = 'matrices.op4',
    output_format: str = 'op4',
    unit: int = 51,
    form: str = 'FORMATTED',
    compile_module: Optional[str] = None,
    alter_location: Optional[str] = None,
) -> None:
    """Add DMAP ALTER to a BDF model for exporting matrix datablocks.

    Modifies model.system_command_lines and model.executive_control_lines
    in-place.

    Parameters
    ----------
    model : BDF
        The BDF model to modify
    matrices : list[str]
        Matrix names to export (e.g., ['KGG', 'MGG', 'KDICT'])
    output_filename : str
        Output filename (e.g., 'kgg.op4' or 'kgg.op2')
    output_format : str
        'op4' for OUTPUT4 or 'op2' for OUTPUT2
    unit : int
        Fortran unit number for the output file (default 51)
    form : str
        For op4: 'FORMATTED' (text) or 'UNFORMATTED' (binary)
    compile_module : str, optional
        Override the COMPILE subdmap name
    alter_location : str, optional
        Override the ALTER location string

    """
    sol = model.sol
    if sol is None:
        raise RuntimeError("model.sol must be set before adding ALTER")

    system_lines, exec_lines = make_matrix_export_alter(
        matrices, output_filename, output_format, sol, unit, form,
        compile_module=compile_module, alter_location=alter_location)

    # Add ASSIGN lines to system_command_lines
    for line in system_lines:
        if line not in model.system_command_lines:
            model.system_command_lines.append(line)

    # Insert ALTER lines before CEND in executive_control_lines
    # Find CEND and insert before it
    cend_idx = None
    for i, line in enumerate(model.executive_control_lines):
        if line.strip().upper() == 'CEND':
            cend_idx = i
            break

    if cend_idx is not None:
        for j, exec_line in enumerate(exec_lines):
            model.executive_control_lines.insert(cend_idx + j, exec_line)
    else:
        # No CEND found, just append
        model.executive_control_lines.extend(exec_lines)
        model.executive_control_lines.append('CEND')


def read_dmap(lines):
    lines2 = []
    append_flag = False
    for line in lines:
        line = line.strip()
        if line == '$':
            continue
        elif len(line) == 0:
            continue


        if line.endswith(','):
            lines2.append(line)
            append_flag = True
        elif append_flag:
            lines2[-1] += line
            append_flag = False
        else:
            lines2.append(line)
            append_flag = False

    lines = [line.upper() for line in lines2]
    del lines2

    #for line in lines:
        #print(line)

    sol = None
    files = {}
    lines2 = []
    code = ""
    echo = True
    level = 0
    spaces = ''
    for line in lines:
        line = ' '.join(line.split())

        comment = ''
        if '$' in line:
            line, comment = line.split('$', 1)
            comment = comment.strip()
        line = line.strip()
        if len(line) == 0:
            continue

        if echo:
            print(comment.rstrip())
            print('line=%r' % line.rstrip())

        if line.startswith('ASSIGN'):
            line = line.replace(' =', '=').replace('= ', '=')
            sline = line.split()
            if len(sline) == 3:
                unused_assign, output, unit = sline
                word, unit_num = unit.split('=')
                output, fname = output.split('=')
                assert word == 'UNIT', word
                assert output == 'OUTPUT2', output
                unit_num = int(unit_num)
                files[unit_num] = [output, fname, unit]
                code = 'from pyNastran.op2.op2 import OP2\n'
                code += 'op2 = OP2()\n'
                #code += ''
            else:
                raise SyntaxError(sline)
        elif line.startswith('IF'):
            assert line.endswith(') THEN'), line
            line = line.replace(') THEN', '):')
            if '=' in line:
                code += line.replace('=', '==').lower() + '\n'
            else:
                raise RuntimeError(line)
            level += 1

        elif line.startswith('TYPE PARM,'):
            line = line[5:]
            #print('line = %r' % line)
            sline = line.split(',')

            # broken for:
            #   TYPE PARAM,,CS,Y,ALPHA=(1.,1.)
            #
            # works for
            #    TYPE PARAM,,CS,Y,ALPHA=1.
            unused_base, unused_NDDL, unused_Type, unused_y, word_val = sline
            param_name, default_value = word_val.split('=')
            code += 'model.params[%r].set_value(%s)\n' % (param_name, default_value)

        elif line.startswith('LAMX '):
            sline = line[5:].split('/')
            print('sline =', sline)
            if len(sline) == 4:
                pre, unused_b, unused_c, resflag = sline
                unused_resflg = int(resflag)
                if pre.startswith('FREQMASS,'):
                    raise NotImplementedError('FREQMASS')
                elif pre.startswith(',,LAMA'):
                    code += spaces + 'LAMX = LAMA\n'
                    if resflag == 1:
                        code += spaces + 'LAMX.msg = "Before Augmentation of Residual Vectors"\n'
                    elif resflag == 2:
                        code += spaces + 'LAMX.msg = "After Augmentation of Residual Vectors"\n'
                    else:
                        code += spaces + 'LAMX.msg = ""\n'
                else:
                    raise RuntimeError(line)
        elif line.startswith('ECHOOFF'):
            #echo = False
            pass
        elif line.startswith('MATMOD'):  # super hard
            pass
        elif line.startswith('MATGEN'):
            line = line[7:].lstrip(',/')
            sline = line.split('/')
            print(sline)
            name = sline[0]
            option = int(sline[1])
            if option == 6:
                sline = sline[2:]
                code += spaces + '%s = \n'
            else:
                raise NotImplementedError(option)

        elif line.startswith('SOL '):
            word, sol = line.split()
            sol = int(sol)
        elif line.startswith('DIAG '):
            # request DIAGnostic output
            pass
        elif line.startswith('COMPILE USERDMAP'):
            pass
        elif line.startswith('COMPILE SEDRCVR'):
            pass
        elif line.startswith('ALTER'):
            pass
        elif line.startswith('DMIIN DMI,DMINDX/'):
            sline = line.split('/')
            pre, data, unused_post = sline
            matrices = data.replace(' ', '').rstrip(',')
            matrices_to_load = matrices.split(',')
            print('matrices_to_load =', matrices_to_load)

            code += 'from pyNastran.bdf.bdf import BDF\n'
            code += 'model = BDF()\n'
            code += 'model.read_bdf(%s, xref=True)\n\n' % 'bdf_filename'
            for matrix_name in matrices_to_load:
                code += "%s = model.dmis[%r].get_matrix()\n" % (matrix_name, matrix_name)
            del sline
            code += '\n'
        elif line.startswith('MPYAD '):
            sline = line[6:].replace(' ', '').split(',')
            a_name, b_name, other = sline
            sline2 = other.strip('/').split('/')
            if len(sline2) == 2:
                out_name, transpose_method = sline2
                transpose_method = int(transpose_method)
                if transpose_method == 1:
                    code += '%s = %s.T @ %s\n' % (out_name, a_name, b_name)
                elif transpose_method == 0:
                    code += '%s = %s @ %s\n' % (out_name, a_name, b_name)
                else:
                    raise RuntimeError(transpose_method)
            else:
                raise RuntimeError(sline2)

            print('line=%r' % line)
        elif line.startswith('OUTPUT2 '):
            sline = line[8:].replace(' ', '').split('/')
            if len(sline) == 4:
                output_matrices, unused_dunno, unused_method, unit_num = sline
                matrix_names = output_matrices.split(',')
                code += '\n'
                for name in matrix_names:
                    code += 'op2.write(%s)\n' % name
            else:
                raise RuntimeError(sline)

            #print('line=%r' % line2)
            #adsf
        elif line.startswith('CEND'):
            pass
        else:
            msg = 'line=%r\n\n' % line
            msg += code
            raise SyntaxError(msg)

    print('-----------')
    print(code)

def main():
    lines = [
        "ASSIGN OUTPUT2='myMatrix.op2',",
        '    UNIT=15',
        '$',
        'SOL 100',
        'DIAG 8,44',
        'COMPILE USERDMAP',
        'ALTER 2',
        'DMIIN DMI,DMINDX/A,B,MYDOF,,,,,,,/ $',
        'MPYAD A,B,/ATB/1///$',
        'MPYAD B,A,/BTA/1///$',
        'OUTPUT2 A,B,ATB,BTA,MYDOF//0/15$',
        'CEND',
    ]
    #with open('isat.bdf', 'r') as f:
        #lines = f.readlines()
    read_dmap(lines)

if __name__ == '__main__':  # pragma: no cover
    main()
