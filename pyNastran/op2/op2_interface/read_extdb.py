"""
read_extdb(op2_reader)
"""
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np

from pyNastran.op2.errors import FortranMarkerError

from pyNastran.op2.op2_interface.utils import (
    mapfmt, reshape_bytes_block,)

from .extdb import read_extdb_extdb, read_extdb_geomx, read_extdb_phip

from pyNastran.op2.op2_interface.read_matrix import (
    read_matrix_mat)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2

def read_tug1(op2_reader: OP2Reader) -> None:
    """TUG1: table UG1"""
    op2: OP2 = op2_reader.op2
    table_name = op2_reader._read_table_name(rewind=False)
    op2.table_name = table_name
    #op2_reader.log.debug('table_name = %r' % op2.table_name)
    if op2.is_debug_file:
        op2.binary_debug.write('_read_superelement - %s\n' % table_name)

    #data = op2_reader._read_record()
    #op2_reader.show_data(data, types='ifsqd', endian=None, force=False)
    #op2_reader.show(data, types='ifsqd', endian=None, force=False)

    op2_reader.read_markers([-1])
    if op2.read_mode == 1:
        data0 = op2_reader._read_record()
        out0 = Struct(b'7i').unpack(data0)
        # (101, 240, 34, 32767, 0, 7, 0)
        # (101, 184, 34, 32767, 0, 1004, 0)
        #(101, 132, 18, 32767, 0, 7, 0)
        aa, bb, ngrid, cc, dd, seven, ee = out0
        #assert ngrid == 34, out
        data_dict = {
            -1: {
                'data': out0,
                '3_ngrid': ngrid,
            },
        }
        analysis_code_expected = 0
        if table_name == b'TUG1':
            analysis_code_expected = 7
            table_subname_expected = 'PHIP'
        elif table_name == b'TES1':
            analysis_code_expected = 5
            table_subname_expected = 'TES'
        elif table_name == b'TEF1':
            analysis_code_expected = 4
            table_subname_expected = 'TEF'
        elif table_name == b'TEE1':
            analysis_code_expected = 5
            table_subname_expected = 'TEE'
        elif table_name == b'TQMG1':
            analysis_code_expected = 39
            table_subname_expected = 'TQMP'
        else:
            raise NotImplementedError(op2.table_name)
        assert seven in {7, 1004, 1005, 39}, f'seven={seven:d} out0={out0}'
    else:
        data = op2_reader._skip_record()

    op2_reader.read_markers([-2, 1, 0])
    if op2.read_mode == 1:
        data = op2_reader._read_record()
        out = Struct('8s i 12s i').unpack(data)
        #print(out)
        table_subname, analysis_code, table_type, precision = out # TODO: verify
        data_dict[-2] = {
            'data': out,
            '1_table_subname': table_subname,
            '2_analysis_code': analysis_code,
            '3_table_type': table_type,
            '4_precision?': precision,
        }
        table_subname = table_subname.strip().decode(op2_reader._encoding)
        table_type = table_type.strip().decode(op2_reader._encoding)
        assert table_subname == table_subname_expected, table_subname # {'PHIP', 'TEF', 'TES', 'TEE', 'TQMP'}, table_subname
        assert analysis_code == analysis_code_expected, out # in {7, 4, 5, 39}, out
        assert table_type == 'REAL', table_type
        assert precision == 1, precision
        op2.log.info(f'table_name={op2.table_name} out0={out0} out={out}')
    else:
        data = op2_reader._skip_record()
    #op2_reader.show_data(data, types='ifsq', endian=None, force=False)

    itable = -3
    op2_reader.read_markers([itable, 1, 0])
    itable -= 1
    marker = op2_reader.get_marker1(rewind=True)
    while marker != 0:
        if op2.read_mode == 1:
            data = op2_reader._read_record()
            ints = np.frombuffer(data, dtype='int32')
            if ints.max() > 500_000_000:
                assert table_name in {b'TES1', b'TEF1', b'TEE1'}, table_name
                strings = np.frombuffer(data, dtype='|S8')

                #op2_reader.show_data(data, types='ifsq', endian=None, force=False)
                #data2 = data.split()
                try:
                    ints = np.array(strings, dtype='int32')
                except ValueError:
                    ints = strings
            else:
                assert table_name in {b'TUG1', b'TEF1', b'TES1', b'TEE1', b'TQMG1'}, table_name
            if ints.dtype.name in {'int32'}:
                ints = ints.reshape(len(ints) // 2, 2)
            #print(ints, ints.shape)
            data_dict[itable+1] = {'data': ints}
        else:
            data = op2_reader._skip_record()

        op2_reader.read_markers([itable, 1, 0])
        itable -= 1
        marker = op2_reader.get_marker1(rewind=True)
        if marker == 0:
            break

        if op2.read_mode == 1:
            data = op2_reader._read_record()

            #header = b'TYPE  IDCOMP ROW    TYPE  IDCOMP ROW    '
            header = data[:40]
            #print(header)
            ints2 = np.frombuffer(data[40:], dtype='int32')
            nints = len(ints2)
            ints2 = ints2.reshape(nints//5, 5)
            #print(ints2)
            #print(ints2.shape)
            data_dict[itable+1] = {'header': header, 'data': ints2}
        else:
            data = op2_reader._skip_record()

        op2_reader.read_markers([itable, 1, 0])
        itable -= 1
        #op2_reader.show_data(data, types='ifsqd', endian=None, force=False)
        marker = op2_reader.get_marker1(rewind=True)

    op2_reader.read_markers([0])
    if op2.read_mode == 1:
        table_name_str = table_name.decode(op2_reader._encoding)
        op2.op2_results.superelement_tables[table_name_str] = data_dict
    return


def read_mef1(op2_reader: OP2Reader) -> None:
    read_matrix_mat(op2_reader)
    return
    op2: OP2 = op2_reader.op2
    endian = op2_reader._endian
    #size = op2_reader.size

    op2.table_name = op2_reader._read_table_name(rewind=False)
    #op2_reader.log.debug('table_name = %r' % op2.table_name)
    if op2.is_debug_file:
        op2.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
    op2_reader.read_markers([-1])
    data = op2_reader._read_record()
    op2_reader.show_data(data, types='i', endian=None, force=False)

    op2_reader.read_markers([-2, 1, 0])
    data = op2_reader._read_record()
    op2_reader.show_data(data, types='is', endian=None, force=False)

    #op2_reader.show_data(data, types='ifs', endian=None, force=False)

    #op2_reader.show(1000, types='ifsqd')

    itable = -3
    itable_next = -4

    s12 = Struct(endian + b'i d')
    s20 = Struct(endian + b'i dd')
    while 1:
        op2_reader.read_markers([itable, 1, ])
        stop_marker = op2_reader.get_marker1(rewind=True)
        if stop_marker == 0:
            break
        assert stop_marker == 1, stop_marker
        stop_marker = op2_reader.get_marker1(rewind=False)

        while 1:
            n = op2.n
            data = op2_reader.read_block4()
            block_value, = op2.struct_i.unpack(data)
            if block_value == itable_next:
                op2.f.seek(n)
                op2.n = n
                break
            data = op2_reader.read_block4()
            ndata = len(data)
            if ndata == 12:
                out = s12.unpack(data)
                print(block_value, out)
            elif ndata == 20:
                out = s20.unpack(data)
                print(block_value, out)
            else:
                out = np.frombuffer(data[4:], dtype='float64')
                print(block_value, out.tolist())
                #op2_reader.show_data(data, types='ifsqd', endian=None, force=False)
                #raise RuntimeError()
            #op2_reader.show(100, types='ifsqd')
        itable -= 1
        itable_next -= 1
    stop_marker = op2_reader.get_marker1(rewind=False)
    #op2_reader.show(1000, types='ifsqd', endian=None, force=False)

    #for table=-3
    ##ints = (
    ##  4, 1, 4,
    ##  4, 2, 4, 12, 4, 858993456, 1039588147, 12,
    ##  4, 2, 4, 12, 6, 0, -1115684864, 12,
    ##  4, 2, 4, 12, 10, 0, 1039138816, 12,
    ##  4, 2, 4, 12, 12, 0, 1039138816, 12,
    #   4, 2, 4, 12, 15, 0, -1109393408, 12,
    #   4, 2, 4, 12, 18, 0, 1039663104, 12,
    #   4, 2, 4, 12, 20, 0, 1039663104, 12,
    #   4, 2, 4, 12, 23, 0, -1110441984, 12,
    #   4, 2, 4, 12, 26, 0, 1039663104, 12,
    #   4, 2, 4, 12, 28, -1717986904, -1109093991, 12,
    #   4, 4, 4, 20, 30, 0, 1032323072, 0, -1109393408, 20,
    #   4, 2, 4, 12, 34, 0, -1107820544, 12,
    #   4, 2, 4, 12, 36, 1717986916, -1106434970, 12,
    #   4, 2, 4, 12, 38, 0, 1032323072, 12,
    #   4, 2, 4, 12, 42, 0, -1107034112, 12,
    #   4, 2, 4, 12, 44, 858993460, -1106079181, 12,
    #   4, 4, 4, 20, 46, 0, 1032323072, 0, 1038090240, 20,
    #   4, 2, 4, 12, 50, 0, -1105920000, 12,
    #   4, 2, 4, 12, 52, 1717986919, -1105264666, 12,
    #   4, 2, 4, 12, 54, 0, 1032585216, 12,
    #   4, 2, 4, 12, 58, 0, -1105297408, 12,
    #   4, 2, 4, 12, 60, 0, -1105014208, 12,
    #   4, 2, 4, 12, 62, 0, 1032060928, 12,
    #   4, 2, 4, 12, 66, 0, -1104949248, 12,
    #   4, 2, 4, 12, 68, -1717986918, -1104605303, 12,
    #   4, 4, 4, 20, 70, 0, 1031995392, 0, 1038352384, 20,
    #   4, 2, 4, 12, 74, -208725392, -1060776204, 12, 4, 2, 4, 12, 76, -1785409415, -1060783274, 12, 4, 4, 4, 20, 78, 1273831424, -1075489352, -82081896, -1057262893, 20, 4, 2, 4, 12, 82, -1785409368, -1060783274)

    #ints = (
    # 4, -4, 4,
    # 4, 1, 4,
    # 4, 1, 4,
    # 4, 2, 4, 12, 1, 0, -1101004800, 12,
    # 4, 2, 4, 12, 3, -858993472, 1045130444, 12,
    # 4, 2, 4, 12, 5, 0, -1108344832, 12, 4, 2, 4, 12, 19, 858993456, 1047976755, 12, 4, 2, 4, 12, 21, 0, -1107296256, 12, 4, 4, 4, 20, 24, 0, -1118830592, 0, 1048051712, 20, 4, 2, 4, 12, 27, -1717986916, 1049062809, 12, 4, 2, 4, 12, 29, 0, -1107296256, 12, 4, 4, 4, 20, 32, 0, -1117782016, 0, 1049362432, 20, 4, 2, 4, 12, 35, 858993458, 1050055219, 12, 4, 2, 4, 12, 37, 0, -1106771968, 12, 4, 4, 4, 20, 40, 0, -1117519872, 0, 1050148864, 20, 4, 2, 4, 12, 43, -1717986918, 1050691865, 12, 4, 2, 4, 12, 45, 0, -1106771968, 12, 4, 4, 4, 20, 48, 0, -1116602368, 0, 1050705920, 20, 4, 2, 4, 12, 51, 0, 1051174016, 12, 4, 2, 4, 12, 53, 0, -1105985536, 12, 4, 4, 4, 20, 56, 0, -1116733440, 0, 1051099136, 20, 4, 2, 4, 12, 59, -858993460, 1051473612, 12, 4, 2, 4, 12, 61, 0, -1106247680, 12, 4, 4, 4, 20, 64, 0, 1034108928, 0, 1051574272, 20, 4, 2, 4, 12, 67, 0, 1052010048, 12, 4, 2, 4, 12, 69, 0, -1105985536, 12, 4, 4, 4, 20, 72, 0, 1036075008, -247389120, -1050902653, 20, 4, 2, 4, 12, 75, 125277002, -1051393547, 12, 4, 2, 4, 12, 77, -993396224, -1059642086, 12, 4, 4, 4, 20, 80, -1898174180, -1065569672, 125276976, -1051393547, 20, 4, 2, 4, 12, 83, -1556319668, -1053694050, 12, 4)

    #data = op2_reader._read_record()
    #op2_reader.show_data(data, types='ifsqd', endian=None, force=False)

    #op2_reader.show(100, types='ifsqd')
    #op2_reader.show_data(data, types='ifsqd', endian=None, force=False)


def read_extdb(op2_reader: OP2Reader) -> None:
    r"""
    fails if a streaming block:
     - nx_spike\extse04c_0.op2

    DB NAME   op2 Name  Description
    =======   ========  ===========
    XSOP2DIR  XSOP2DIR  Table of contents of the op2 file.
    GEOM1X    GEOM1X    GRID point geometry.
    GEOM2X    GEOM2X    SPOINTs.
    GEOM4X    GEOM4X    ASET/ASET1 and QSET/QSET1 entries.
    EXTDB     MATK      Boundary stiffness.
    EXTDB     MATM      Boundary mass.
    EXTDB     MATP      Boundary loads.
    EXTDB     MATV      Boundary fluid-structure partitioning vector.
    EXTDB     TUG1      Displacement OTM table.
    EXTDB     MUG1      Displacement OTM matrix.
    EXTDB     TES1      Stress OTM table.
    EXTDB     MES1      Stress OTM matrix.
    EXTDB     TEF1      Force OTM table.
    EXTDB     MEF1      Force OTM matrix.
    EXTDB     MUG1B     Displacement OTM matrix in basic coordinates system.

    https://help.hexagonmi.com/bundle/MSC_Nastran_2022.2/page/Nastran_Combined_Book/Superelements_User_s_Guide/superOTM/TOC.OUTPUT2.Files.xhtml

    DB_name   op2_name  Description
    XSOP2DIR  XSOP2DIR  Table of contents of the op2 file (always the first datablock to appear).
    GEOM1X    GEOM1X    GRID point geometry.
    GEOM2X    GEOM2X    SPOINTs.
    GEOM4X    GEOM4X    ASET/ASET1 and QSET/QSET1 entries.

    For AVL EXB output
    EXTDB  LAMAAVP  Eigenvalue table for AVL POST.
    EXTDB  MATAPH   Eigenvectors for AVL POST.
    EXTDB  MATAEK   Diagonal matrix of eigenvalues for AVL POST.
    EXTDB  MATAM0   Generalized mass for AVL POST.

    For Adams MNF output
    EXTDB  MATAKA   Boundary stiffness matrix in basic coordinates (Adams POST).
    EXTDB  MATPH2   Eigenvectors in basic coordinates (Adams POST).
    EXTDB  MATAMA   Boundary mass matrix in basic coordinates with WTMASS removed (Adams POST).
    EXTDB  MATK     Boundary stiffness.
    EXTDB  MATM     Boundary mass.
    EXTDB  MATB     Boundary viscous damping.
    EXTDB  MATK4    Boundary structural damping.
    EXTDB  MATP     Boundary loads.
    EXTDB  MATV     Boundary fluid-structure partitioning vector.
    EXTDB  MATGP    Aerodynamic transformation matrix for loads.
    EXTDB  MATGD    Aerodynamic transformation matrix for displacements.
    EXTDB  MATRGA   Unit transformation from boundary to interior DOF.
    EXTDB  MATVAFS  Fluid-structure partitioning vector.
    EXTDB  MATA     Partitioned acoustic coupling.
    EXTDB  MATRV    Residual vector partitioning vector.
    EXTDB  MATPC    Access points.
    EXTDB  MATKSM   Aerodynamic generalized stiffness.
    EXTDB  MATMSM   Aerodynamic generalized mass.

    For rotors
    NAMELIST  NAMELIST  List of rotors.
    MTRXNAME  MTRXNAME  Rotor matrices defined in NAMELIST.

    #-----------------
    DB_name  op2_name  Description
    EXTDB  TUG1     Displacement OTM table.
    EXTDB  MUG1     Displacement OTM matrix.
    EXTDB  MUG1O    Displacement OTM for loaded interior DOF.
    EXTDB  TQG1     SPCFORCE OTM table.
    EXTDB  MKQG1    SPCFORCE stiffness contribution OTM.
    EXTDB  MBQG1    SPCFORCE viscous damping contribution OTM.
    EXTDB  MMQG1    SPCFORCE mass contribution OTM.
    EXTDB  MK4QG1   SPCFORCE structural damping contribution OTM.
    EXTDB  MKQG1O   SPCFORCE stiffness contribution OTM for loaded interior DOF.
    EXTDB  TELAF1   Elastic element force OTM table.
    EXTDB  MELAF1   Elastic element force OTM matrix.
    EXTDB  TES1     Stress OTM table.
    EXTDB  MES1     Stress OTM matrix.
    EXTDB  MES1O    Stress OTM matrix for loaded interior DOF.
    EXTDB  TEF1     Force OTM table.
    EXTDB  MEF1     Force OTM matrix.
    EXTDB  MEF1O    Force OTM matrix for loaded interior DOF.
    EXTDB  MUG1B    Displacement OTM matrix in basic coordinates system.
    EXTDB  MUG1OB   Displacement OTM matrix in basic coordinates system for loaded interior DOF.
    EXTDB  TEE1     Strain OTM table.
    EXTDB  MEE1     Strain OTM matrix.
    EXTDB  MEE1O    Strain OTM matrix for loaded interior DOF.
    EXTDB  TQMG1    MPCFORCE OTM table.
    EXTDB  MKQMG1   MPCFORCE stiffness contribution OTM.
    EXTDB  MBQMG1   MPCFORCE viscous damping contribution OTM.
    EXTDB  MMQMG1   MPCFORCE mass contribution OTM.
    EXTDB  MK4QMG1  MPCFORCE structural damping contribution OTM.
    EXTDB  MKQMG1O  MPCFORCE stiffness contribution OTM for loaded interior DOF.

    Some matrices are also named with pseudo-degree-of-freedom set names.
    W – The set omitted after auto-omit (a-set combines x-set and w-set)
    X – The set retained after auto-omit (complement of w-set)
    J – Superelement interior degrees-of-freedom; for example, KJJ and PJ
    H – Modal degrees-of-freedom; for example, PHDH, MHH, PHF and UHF
    """
    endian = op2_reader._endian
    size = op2_reader.size

    # gotta read the tables at the beginning because we need to flag the tables?
    read_mode = 1

    #otm_tables   = ['TUG1',          'TEF1', 'TES1']
    #otm_matrices = ['MUG1', 'MUG1B', 'MEF1', 'MES1']
    # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\extse04c_cnv1_0.op2
    op2: OP2 = op2_reader.op2
    log: SimpleLogger = op2.log
    op2.table_name = op2_reader._read_table_name(rewind=False)
    #op2_reader.log.debug('table_name = %r' % op2.table_name)
    if op2.is_debug_file:
        op2.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
    op2_reader.read_markers([-1])
    if op2.is_debug_file:
        op2.binary_debug.write('---markers = [-1]---\n')

    #(101, 1, 0, 8, 0, 0, 0)
    data = op2_reader._read_record()

    markers = op2_reader.get_nmarkers(1, rewind=True)
    if op2.is_debug_file:
        op2.binary_debug.write('---marker0 = %s---\n' % markers)

    op2_reader.read_3_markers([-2, 1, 0])
    marker = -3
    if op2.read_mode == read_mode or op2.make_geom is False:
        if op2.read_mode == read_mode and op2.make_geom is False:
            log.warning('reading the EXTRN tables requires the read_op2_geom')
        data = op2_reader._skip_record()
        while 1:
            #print('====================')
            #print(f'***reading {marker}')
            try:
                op2_reader.read_markers([marker, 1])
            except FortranMarkerError:
                #op2.show_ndata(100)
                raise
            nfields = op2_reader.get_marker1(rewind=True)
            if nfields == 0:
                break
            data, ndata = op2_reader._read_record_ndata()
        return

    struct_3i = Struct(mapfmt(endian + b'3i', size))

    # drop XSOP2DIR and PVT0
    iextdb = op2.table_count[b'EXTDB'] + 0

    log.debug('-'*80)
    xsop2dir_name = op2_reader.xsop2dir_names[iextdb]
    if op2.read_mode == 2:
        data, ndata = op2_reader._read_record_ndata()
        #op2.show_data(data, types='ifsqd', endian=None, force=False)

        name = ''
        name1 = ''
        dtype_str = ''
        nfields = ndata // size
        assert ndata % size == 0

        assert size in {4, 8}, size
        if nfields == 2:
            if size == 4:
                name_bytes, = Struct(endian + b'8s').unpack(data)
                name = name_bytes.decode('latin1').rstrip()
            else:
                name_bytes, = Struct(endian + b'16s').unpack(data)
                name = reshape_bytes_block(name_bytes).decode('latin1').strip()
            log.info(f'{marker} A: name={name!r} -> {xsop2dir_name!r}')
            assert name in {'GEOM1', 'GEOM2', 'IGEOM2X', 'GEOM4', 'EXTDB'}, name
        elif nfields == 4:
            if size == 4:
                name_bytes, int1, int2 = Struct(endian + b'8s 2i').unpack(data)
                name = name_bytes.decode('latin1').rstrip()
            else:
                name_bytes, int1, int2 = Struct(endian + b'16s 2q').unpack(data)
                name = reshape_bytes_block(name_bytes).decode('latin1').strip()
            log.info(f'{marker}: B: name={name!r} -> {xsop2dir_name!r} int1={int1} int2={int2}')
        elif nfields == 7:
            # (PHIP, 7, REAL, 1)
            if size == 4:
                name1_bytes, int1, dtype_bytes, int2 = Struct(endian + b'8s i 12s i').unpack(data)
                name1 = name1_bytes.decode('latin1').rstrip()
                dtype_str = dtype_bytes.decode('latin1').strip()
            else:
                #op2_reader.show_data(data, types='ifsdq')
                name1_bytes, int1, dtype_bytes, int2 = Struct(endian + b'16s q 24s q').unpack(data)
                name1 = reshape_bytes_block(name1_bytes).decode('latin1').strip()
                dtype_str = reshape_bytes_block(dtype_bytes).decode('latin1').strip()

            #DTI         TUG1       0
            #            PHIP       ?       7    ?          ?    REAL       1  ENDREC
            #DTI         TUG1       0     696     606   32767       0       7       0
            #            PHIP   32767       7   32767   32767    REAL       1  ENDREC
            dti = op2.add_dti(
                xsop2dir_name,
                {0: ['?', '?', 32767, 0, int1, 0,
                     name1, 32767, int1, 32767, 32767, dtype_str, int2, 'ENDREC',]})
            #print(dti)
            # m.add_dmi(name, form, tin, tout, nrows, ncols, GCj, GCi, Real, Complex=None, comment='')
            #print(name1, int1, name2, int2, 28)

            log.info(f'{marker}: C DTI: name1={name1!r} -> {xsop2dir_name!r} int1={int1} dtype_str={dtype_str!r} int2={int2}')

    #op2_reader.show(200)

    ints_ = []
    doubles_ = []
    while 1:
        #print('====================')
        #print(f'***reading {marker}')
        try:
            op2_reader.read_markers([marker, 1])
        except FortranMarkerError:
            #op2.show_ndata(100)
            raise

        nfields1 = op2_reader.get_marker1(rewind=True)
        if nfields1 == 0:
            nfields1 = op2_reader.get_marker1(rewind=False)
        elif nfields1 == 1:
            #data, ndata = op2_reader._read_record_ndata()
            nfields1 = op2_reader.read_markers([1])

            nfields_test = op2_reader.get_marker1(rewind=True)

            #out=(1, 1.0) nblock=16
            #out=(2, 1.0) nblock=16
            #...
            #out=(606, 1.0) nblock=16
            itest = 0
            while nfields_test > 0:  # nfields_test == 1
                nfields = op2_reader.get_marker1(rewind=False)
                block = op2_reader._read_block()
                nblock = len(block)
                ndouble = (nblock - 4) // 8
                fmt = mapfmt(endian + b'i%dd' % ndouble, size)
                out = Struct(fmt).unpack(block)
                #print(f'out={out}')
                nfields_test = op2_reader.get_marker1(rewind=True)
                ints_.append(out[0])
                doubles_.append(out[1])
                itest += 1
            #log.debug(f'itest={itest}')

            #print('-------')
            #print(f'end of marker={marker}')
            marker -= 1
            #marker = op2_reader.get_marker1(rewind=True)
            continue
        else:
            raise RuntimeError('EXTDB error')

        #if op2.read_mode == 2:
        if len(ints_):
            ints = np.array(ints_, dtype='int32')
            doubles = np.array(doubles_, dtype='float64')
            print(len(ints), ints, doubles)
            #del ints_, doubles_


        #op2.show_ndata(100)
        nfields = op2_reader.get_marker1(rewind=True)
        #print('nfields =', nfields)
        if nfields == 0:
            #if op2.read_mode == 2:
                #op2_reader.log.warning('breaking...')
            #op2_reader.show(200)
            #log.debug(f'ints={ints_} doubles={doubles_}')
            break
        #log.debug(f'ints={ints_} doubles={doubles_}')
        # ----------------------------------------------------------------------

        data, ndata = op2_reader._read_record_ndata()
        if op2.read_mode == 2:
            #if name not in ['GEOM1', 'GEOM2', 'GEOM2X', 'IGEOM2X', 'GEOM4', 'EXTDB']:
                #if ndata != 12:
                    #op2_reader.log.warning(f'--B; ndata={ndata}--')

            if len(name):
                log.info(f'name = {name!r} size={size}')
            nfields = ndata // size
            assert ndata % size == 0

            if size == 4:
                numpy_idtype = 'int32'
                struct_idtype = 'i'
            else:
                numpy_idtype = 'int64'
                struct_idtype = 'q'

            if nfields == 3:
                intsi = struct_3i.unpack(data)
                if intsi == (65535, 65535, 65535):
                    name = ''
                    int1 = 65535
                    int2 = 65535
                elif intsi == (65535, 65535, 25535):
                    name = ''
                    int1 = 65535
                    int2 = 25535
                else:
                    if size == 4:
                        op2_reader.show_data(data, types='ifs', endian=None, force=False)
                        name, int1, int2 = Struct(endian + b'4s 2i').unpack(data)
                        raise RuntimeError(name, int1, int2)
                    else:
                        ints = Struct(endian + b'3q').unpack(data)
                        assert ints == (65535, 65535, 65535), ints
            elif name == 'GEOM1':
                read_extdb_geomx(op2_reader, data, endian, op2.reader_geom1.geom1_map)
            elif name in ['GEOM2', 'GEOM2X', 'IGEOM2X']:
                read_extdb_geomx(op2_reader, data, endian, op2.reader_geom2.geom2_map)
            elif name == 'GEOM4':
                read_extdb_geomx(op2_reader, data, endian, op2.reader_geom4.geom4_map)
            elif name == 'EXTDB':
                read_extdb_extdb(op2_reader, xsop2dir_name, data, endian, numpy_idtype)
            elif name1 == 'PHIP':
                read_extdb_phip(op2_reader, xsop2dir_name, name1, marker, data, endian, struct_idtype, op2.idtype8)
            elif name1 == 'TES':
                read_extdb_phip(op2_reader, xsop2dir_name, name1, marker, data, endian, struct_idtype, op2.idtype8)
            elif name1 == 'TEF':
                read_extdb_phip(op2_reader, xsop2dir_name, name1, marker, data, endian, struct_idtype, op2.idtype8)
            elif name1 == 'TQP':
                read_extdb_phip(op2_reader, xsop2dir_name, name1, marker, data, endian, struct_idtype, op2.idtype8)
            else:
                #op2_reader.show_data(data, types='sqd')
                log.warning(f'EXTDB; name={name!r} name1={name1!r} ndata={ndata}')
                raise RuntimeError(f'EXTDB; name={name!r} name1={name1!r} ndata={ndata}')
        marker -= 1
        #print('--------------------')
    unused_marker_end = op2_reader.get_marker1(rewind=False)
    #aa
    #if op2.read_mode == 2:
        #op2_reader.log.warning('returning...')
    log.debug('-'*80)
    return
