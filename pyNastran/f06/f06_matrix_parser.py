from typing import Optional
import numpy as np
import scipy.sparse
from cpylog import SimpleLogger, get_logger

TABLES_2D = {'MKLIST'}
MATRICES_DENSE = {'QHHA'}
SKIP_FLAGS = [
    'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O',
    'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O',
    'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N',
    'M O D E L   S U M M A R Y',
    'C A S E    C O N T R O L    E C H O',
    'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y',
    'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R',
    #'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (READ MODULE)',
    #'R E A L   E I G E N V A L U E S',
    'MAXIMUM  DISPLACEMENTS',
    'FLUTTER  SUMMARY',
    #'* * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *',
]
def read_f06_matrices(f06_filename: str,
                      log: Optional[SimpleLogger]=None,
                      nlines_max: int=1_000_000) -> tuple[dict[str, np.ndarray],
                                                          dict[str, np.ndarray]]:
    log = get_logger(log=log, level='debug', encoding='utf-8')
    with open(f06_filename, 'r') as f06_file:
        tables, matrices = _read_f06_matrices(f06_file, log, nlines_max)
    if len(tables):
        log.info('found the following tables in the f06: %s' % (list(tables)))
    if len(matrices):
        log.info('found the following matrices in the f06: %s' % (list(matrices)))
    return tables, matrices

def _read_f06_matrices(f06_file, log: SimpleLogger, nlines_max: int) -> dict[str, np.ndarray]:
    i = 0
    debug = False
    tables = {}
    matrices = {}
    while True:
        line = f06_file.readline()
        i += 1
        if debug:
            log.debug('i={i} {line.strip()}')
        if '* * * END OF JOB * * *' in line:
            #print("****done****")
            break
        iflags = [datai in line for datai in SKIP_FLAGS]
        if any(iflags):
            flag = SKIP_FLAGS[iflags.index(True)]
            #print('skip', line)
            line = f06_file.readline()
            i += 1
            # JANUARY  26, 2022  SIMCENTER NASTRAN  3/12/20   PAGE     2
            while 'SIMCENTER NASTRAN' not in line and 'PAGE' not in line:
                line = f06_file.readline()
                i += 1
                if i > nlines_max:
                    raise RuntimeError(f'{nlines_max:d} lines in file is max?...\n'
                                       'this will be removed once the parser is better tested')
            #log.debug(f'line = {line.rstrip()}')

        if 'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y' in line:
            while 'R E A L   E I G E N V A L U E S' not in line:
                line = f06_file.readline()
                i += 1
                if i > nlines_max:
                    raise RuntimeError(f'{nlines_max:d} lines in file is max?...\n'
                                       'this will be removed once the parser is better tested')
            line = f06_file.readline()
            line = f06_file.readline()
            line = f06_file.readline()
            line_strip = line.strip()
            i += 3

            #line_strip = '1         1        4.637141E+01        6.809655E+00        1.083790E+00        1.000000E+00        4.637141E+01'
            #R E A L   E I G E N V A L U E S
            #MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
            #    NO.       ORDER                                                                       MASS              STIFFNESS
            #        1         1        4.637141E+01        6.809655E+00        1.083790E+00        1.000000E+00        4.637141E+01
            #        2         2        2.369379E+02        1.539279E+01        2.449838E+00        1.000000E+00        2.369379E+02
            Mhh_list = []
            Khh_list = []
            Bhh_list = []
            while len(line_strip):
                sline = line_strip.split()
                if 'NASTRAN' in sline:
                    break
                assert len(sline) == 7, sline
                mode_num, extraction_order, eigenvalue_str, radians, cycles, gen_mass_str, gen_stiffness_str = sline
                eigenvalue = float(eigenvalue_str)
                gen_mass = float(gen_mass_str)
                Mhh_list.append(gen_mass)
                Bhh_list.append(0.)
                Khh_list.append(eigenvalue)
                line = f06_file.readline()
                line_strip = line.strip()
                i += 1

            Mhh = np.array(Mhh_list, dtype='float64')
            Bhh = np.array(Bhh_list, dtype='float64')
            Khh = np.array(Khh_list, dtype='float64')
            isort = np.argsort(Khh)
            matrices['MHH'] = np.diag(Mhh[isort])
            matrices['BHH'] = np.diag(Bhh[isort])
            matrices['KHH'] = np.diag(Khh[isort])
            del line_strip, gen_mass, Mhh_list, Khh_list, Bhh_list


        if line.startswith('0    TABLE'):
            table_name, table, line, i = _read_table(f06_file, line, i, log)
            tables[table_name] = table
            del table_name, table
            debug = False
        elif line.startswith('0      MATRIX '):
            matrix_name, matrix, line, i = _read_matrix(f06_file, line, i, log)
            matrices[matrix_name] = matrix
            del matrix_name, matrix
        #else:
            #print(line)
        #print('----')
        if i > nlines_max:
            raise RuntimeError(f'{nlines_max:d} lines in file is max?...\n'
                               'this will be removed once the parser is better tested')
        if i % 1000 == 0:
            log.debug(f'i={i}')
    return tables, matrices

def _read_matrix(f06_file, line: str, i: int, log: SimpleLogger) -> tuple[str, np.ndarray, str, int]:
    """
    0      MATRIX QHHA     (GINO NAME 101 ) IS A COMPLEX          100 COLUMN X          10 ROW RECTANG  MATRIX.
    0COLUMN        1      ROWS        1 THRU       10 --------------------------------------------------
    ROW
        1) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01
        6)  1.5360E-04,-1.1042E-04  -4.7606E-04, 2.3069E-04   1.0359E-03,-1.5668E-04  -1.3075E-03, 7.8472E-04   2.3471E-04,-4.8359E-04
    """
    header = line[1:].strip()
    #header_sline = header.split()
    header_left, header_midright = header.split('(')
    table_name = header_left.strip().split()[1]

    # (GINO NAME 101 )
    # IS A COMPLEX          100 COLUMN X          10 ROW RECTANG  MATRIX.
    header_mid, header_right = header_midright.split(')')
    del header_midright

    is_complex = ('COMPLEX' in header_right)
    assert 'COLUMN' in header_right, header_right
    assert 'ROW' in header_right, header_right
    header_sline_right = header_right.split()

    icolumn = header_sline_right.index('COLUMN')
    irow = header_sline_right.index('ROW')

    ncolumns_str = header_sline_right[icolumn-1]
    nrow_str = header_sline_right[irow-1]
    ncolumn = int(ncolumns_str)
    nrow = int(nrow_str)
    log.info(f' - {table_name}: ({nrow},{ncolumn}) is_complex={is_complex}')

    line = f06_file.readline()
    i += 1
    col_list = []
    row_list = []
    data_list = []
    while '0COLUMN' in line:
        # column header
        #    0           1        2         3  4          5
        # 0COLUMN        1      ROWS        1 THRU       10 --------------------------------------------------
        # 0COLUMN      100      ROWS        1 THRU       10
        sline = line[1:].strip().split()
        column = int(sline[1])
        row0 = int(sline[3])
        row1 = int(sline[5])
        #log.debug(line)

        line = f06_file.readline().rstrip()
        i += 1
        assert line.strip() == 'ROW', line

        line = f06_file.readline().rstrip()
        i += 1
        row_lines = [line]
        while ')' in line:
            line = f06_file.readline().rstrip()
            i += 1
            if ')' not in line:
                #print("****breaking on", line.strip())
                break
            row_lines.append(line)

        if is_complex:
            row_indexi, datai = _parse_complex_row_lines(row_lines)
        else:
            raise RuntimeError(is_complex)
        #print(column)
        ndatai = len(datai)
        col_indexi = np.ones(ndatai, dtype='int32') * column
        col_list.append(col_indexi)
        row_list.append(row_indexi)
        data_list.append(datai)
        # ROW
        #    1) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01
        #    6) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01

    cols = np.hstack(col_list) # - 1
    rows = np.hstack(row_list) # - 1
    data = np.hstack(data_list)
    #urows = np.unique(rows)
    #ucols = np.unique(cols)
    scipy.sparse.coo_matrix

    dtype = 'float64'
    if is_complex:
        dtype = 'complex128'
    sparse_matrix = scipy.sparse.coo_matrix(
        (data, (rows-1, cols-1)),
        shape=(nrow, ncolumn), dtype=dtype)

    matrix = sparse_matrix
    if table_name in MATRICES_DENSE:
        matrix = sparse_matrix.toarray()
    return table_name, matrix, line, i

def _parse_complex_row_lines(lines: list[str]) -> tuple[int, int]:
    """
    1) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01
    6) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01
    """
    #lines = [
        #'        1) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01',
        #'        6)  1.5360E-04,-1.1042E-04  -4.7606E-04, 2.3069E-04   1.0359E-03,-1.5668E-04  -1.3075E-03, 7.8472E-04   2.3471E-04,-4.8359E-04',
    #]
    #lines = [
        #'        1) -3.5846E+01,-1.3275E+02',
        #'        6)  1.5360E-04,-1.1042E-04',
    #]

    row_index_list = []
    data_list = []
    for line in lines:
        row_id_str, data_str = line.split(')')
        pairs = [line[10:35], line[35:61], line[61:84], line[84:110], line[110:135]]
        row_i = int(row_id_str)
        for pair in pairs:
            if len(pair) == 0:
                break
            real_str, complex_str = pair.split(',')
            reali = float(real_str)
            complexi = float(complex_str)
            datai = complex(reali, complexi)
            data_list.append(datai)
            row_index_list.append(row_i)
            row_i += 1
    row_index = np.array(row_index_list, dtype='int32')
    data = np.array(data_list, dtype='complex128')
    return row_index, data

def _read_table(f06_file, line: str, i: int, log: SimpleLogger) -> tuple[str, np.ndarray, str, int]:
    """
    '0    TABLE   MKLIST                      LINES CONTAINING BINARY  ZERO HAVE BEEN DELETED.\n'
    ''
    ''
    '0RECORD NO.    0'
    '      1)         MKLI        ST  '
    '                                END OF       2 WORD RECORD.'
    '0RECORD NO.    1'
    '      1)  8.50000E-01 1.20000E+00 8.50000E-01 1.60000E+00 8.50000E-01 1.00000E-03 8.50000E-01 3.00000E-03 8.50000E-01 1.00000E-02'
    '     11)  8.50000E-01 3.00000E-02 8.50000E-01 1.00000E-01 8.50000E-01 2.00000E-01 8.50000E-01 4.00000E-01 8.50000E-01 8.00000E-01'
    '                                END OF      20 WORD RECORD.'
    '0RECORD NO.    2'
    '0END OF FILE'
    """
    table_name = line[1:].split()[1]
    records = []
    irecord_expected = 0
    while '0END OF FILE' not in line:
        while '0RECORD' not in line:
            line = f06_file.readline()
            i += 1
        irecord_str = line.split()[-1]
        irecord = int(irecord_str)
        assert irecord == irecord_expected, f'irecord={irecord} irecord_expected={irecord_expected}'
        irecord_expected += 1

        line = f06_file.readline()
        if '0END OF FILE' in line:
            break
        assert ') ' in line, line
        i += 1
        record_lines = [line]
        while ') ' in line:
            line = f06_file.readline()
            i += 1
            if ') ' in line:
                record_lines.append(line)
        # 'END OF       2 WORD RECORD.'
        records.append(record_lines)

    matrix = _parse_table_records(table_name, records, log)
    # MATRIX QHHA   3280 COLUMN X  40 ROW RECTANG  MATRIX.   (40, 3280)
     #MATRIX QHHA     63 COLUMN X   7 ROW RECTANG  MATRIX.   (7, 63)
    log.info(f' - {table_name}: {matrix.shape}')
    return table_name, matrix, line, i

def _parse_table_records(table_name: str, records: list[str], log: SimpleLogger) -> np.ndarray:
    #records = [
        #['      1)         MKLI        ST  \n'],
        #['      1)  8.50000E-01 1.20000E+00 8.50000E-01 1.60000E+00 8.50000E-01 1.00000E-03 8.50000E-01 3.00000E-03 8.50000E-01 1.00000E-02\n',
         #'     11)  8.50000E-01 3.00000E-02 8.50000E-01 1.00000E-01 8.50000E-01 2.00000E-01 8.50000E-01 4.00000E-01 8.50000E-01 8.00000E-01\n']]
    record1 = records[1]

    # breakup the data into index and data
    matrix_index = []
    matrix_data = []
    for line in record1:
        istart_str, data_str = line.split(') ')
        sline = data_str.strip().split()
        floats = [float(val) for val in sline]
        nfloats = len(floats)
        i0 = int(istart_str) - 1
        matrix_indexi = np.arange(i0, i0+nfloats)
        matrix_index.append(matrix_indexi)
        matrix_data.append(floats)

    matrix_indexi = np.hstack(matrix_index)
    data = np.hstack(matrix_data)
    ndata = len(data)
    if table_name in TABLES_2D:
        assert ndata % 2 == 0, ndata
        nrows = ndata // 2
        ncols = 2
        matrix = data.reshape((nrows, ncols))
    else:
        cols = np.zeros(len(matrix_data), dtype='int32')
        raise RuntimeError(table_name)
        #matrix = scipy.sparse.coo_matrix(
            #(data, matrix_indexi),
            #shape=(mrows, ncols), dtype='f')
    #log.debug(str(matrix))
    return matrix

def main():
    f06_filename = r'C:\NASA\uuaero_flutterexplorer\flutter.alter.f06'
    #f06_filename = r'C:\NASA\uuaero_flutterexplorer\models\jetzero_flutter\flutter.alter.f06'
    #f06_filename = r'C:\NASA\uuaero_flutterexplorer\models\demo\examples\crm\crm_flutter_alter.f06'
    tables, matrices = read_f06_matrices(f06_filename)
    MKLIST = tables['MKLIST']
    QHHA = matrices['QHHA']
    nmodes, nmk_nmodes = QHHA.shape
    assert nmk_nmodes % nmodes == 0
    nmk = nmk_nmodes // nmodes
    qhh_shape = (nmodes, nmodes, nmk)
    QHHr = np.zeros(qhh_shape, dtype='float64')
    QHHi = np.zeros(qhh_shape, dtype='float64')

    # todo: verify
    # todo: we could do this with a axis swap?
    QHHAr = QHHA.real.reshape((nmodes, nmodes, nmk)) # v1
    QHHAi = QHHA.imag.reshape((nmodes, nmodes, nmk))

    #QHHAr = QHHA.real.reshape((nmodes, nmk, nmodes)) # v2
    for jj in range(nmk):
        mk = MKLIST[jj, :]
        qhhar = QHHAr[:, :, jj]  # v1; best guess
        #qhhar = QHHAr[:, jj, :]  # v2
        print(mk)
        print(qhhar)
        print('---------')

        #QHHAr = QHHA.real.reshape((nmodes, nmk, nmodes))
        #qhhar = QHHAr[:, jj, :]
        #[[-6.5440e+01 -5.4153e+01 -1.1214e+02  3.6589e+02  9.3012e+01]
        # [ 1.2640e+01 -2.2370e+00  6.1121e-02 -9.4691e+01  3.2896e+01]
        # [ 2.7129e+01 -7.8625e-01  7.0466e+00 -1.9454e+02  5.4725e+01]
        # [-4.7326e+01 -2.2778e+01 -5.4106e+01  2.9052e+02  1.3453e+01]
        # [-1.6278e+01 -3.4003e+01 -5.9597e+01  4.6454e+01  4.0195e+01]]

        #QHHAr = QHHA.real.reshape((nmodes, nmodes, nmk))
        #qhhar = QHHAr[:, :, jj]
        #[[ 3.6634e+02 -1.1222e+02 -5.4153e+01 -6.5440e+01  9.3012e+01]
        # [-9.4877e+01  9.3530e-02 -2.2370e+00  1.2640e+01  3.2896e+01]
        # [-1.9492e+02  7.1120e+00 -7.8625e-01  2.7129e+01  5.4725e+01]
        # [ 2.9107e+02 -5.4204e+01 -2.2778e+01 -4.7326e+01  1.3453e+01]
        # [ 4.6313e+01 -5.9567e+01 -3.4003e+01 -1.6278e+01  4.0195e+01]]
        qhhai = QHHAi[:, :, jj]
        QHHr[:, :, jj] = qhhar
        #if MKLIST[jj, 1] != 0.0:
        QHHi[:, :, jj] = qhhai / MKLIST[jj, 1]
    x = 1

if __name__ == '__main__':
    main()
    x = 1
