from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.op2.result_objects.gpdt import GPDT, BGPDT

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader


def read_gpdt(op2_reader: OP2Reader) -> None:
    """
    reads the GPDT table

    tested by ???

    """
    #if self.read_mode == 1:
        #read_record = self._skip_record
    #else:
    op2: OP2 = op2_reader.op2
    log: SimpleLogger = op2.log
    read_record = op2_reader._read_record
    skip_record = op2_reader._skip_record
    read_mode = op2_reader.read_mode
    size = op2_reader.size
    endian = op2_reader._endian

    table_name = op2_reader._read_table_name(rewind=False)
    op2.table_name = table_name
    #op2_reader.log.debug('table_name = %r' % table_name)
    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('read_gpdt - %s\n' % table_name)

    op2_reader.read_markers([-1])
    header_data = op2_reader._read_record()  # (103, 117, 0, 0, 0, 0, 0)
    header_ints = np.frombuffer(header_data, op2.idtype8)
    header_floats = np.frombuffer(header_data, op2.fdtype8)

    # what a mess...
    # there are 3 ways to define the header...
    # for each of the 3, there are size 7/10 and single/double precision...
    #[ 102 3       0    0    0    0    0]  # 3 nodes....7 words   C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\extse11s.op2
    #[ 102 28      0    0    0    0    0]  # 16 nodes...28 bytes  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\boltld06.op2
    #[ 102 5523    7    0    1    0    0]  # 5523 nodes...7 words C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\acssnbena1.op2
    #print('gpdt ints =', header_ints)
    #print('gpdt floats =', header_floats)
    #seid = ints[0] # ??? is this a table number>
    #nnodes = ints[1]

    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('---markers = [-1]---\n')
    #print('--------------------')

    op2_reader.read_3_markers([-2, 1, 0])
    op2_reader.read_table_name(['GPDT', 'GPDTS', 'SAGPDT'])

    #print('--------------------')

    op2_reader.read_3_markers([-3, 1, 0])


    ## TODO: no idea how this works...
    result_name = 'gpdt'
    #if op2._results.is_not_saved(result_name):
        #read_mode = 2 #  fake the read mode to skip the table
    #else:
        #op2._results._found_result(result_name)

    if read_mode == 1 and op2._results.is_saved(result_name):
        op2._results._found_result(result_name)
        i = 0
        data = read_record() # nid,cp,x,y,z,cd,ps
        ndata = len(data)



        # assume nodes
        #_get_gpdt_nnodes(op2_reader.size, ndata, header_ints, op2_reader.log)
        nnodes, numwide = _get_gpdt_nnodes2(ndata, header_ints, size)
        if (size, numwide) == (8, 14):
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\z402cdamp1_04.op2
            nid_cp_cd_ps, xyz = _read_gpdt_8_14(op2, data, nnodes)
        elif (size, numwide) == (4, 7):
            nid_cp_cd_ps, xyz = _read_gpdt_4_7(op2, data, nnodes)
        elif (size, numwide) == (4, 10):
            nid_cp_cd_ps, xyz = _read_gpdt_4_10(op2, data, nnodes)
        else:
            op2_reader.show_data(data, types='if')
            #ndata = 84:
            #[102   3   0   0   0   0   0]
            #ints    = (1,   0, 150.0,   0,   0,   0,   0,
                       #2, 0, 1125515264, 1103626240, 0, 0, 0,
                       #3, 0, 1125515264, -1043857408, 0, 0, 0)
            #floats  = (1, 0.0, 150.0, 0.0, 0.0, 0.0, 0.0,
                       #2, 0.0, 150.0, 25.0, 0.0, 0.0, 0.0,
                       #3, 0.0, 150.0, -25.0, 0.0, 0.0, 0.0)
            raise NotImplementedError((size, numwide))

        if 0:  # pragma: no cover
            if nvalues % 7 == 0:
                # mixed ints, floats
                #  0   1   2   3   4   5   6
                # id, cp, x1, x2, x3, cd, ps
                nrows = get_table_size_from_ncolumns('GPDT', nvalues, 7)
                #print('size =', size)
                ntotal = 28
                structi = Struct(endian + b'2i 4f i')
                for j in range(10):
                    edata = data[i:i+ntotal]
                    op2_reader.show_data(edata, types='ifqd')
                    out = structi.unpack(edata)
                    i += ntotal
                    #print(out)
                asdf
                ints = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()
                floats = np.frombuffer(data, op2.fdtype8).reshape(nnodes, 7).copy()
                iints = [0, 1, 5, 6] # [1, 2, 6, 7] - 1
                nid_cp_cd_ps = ints[:, iints]
                xyz = floats[:, 2:5]
            elif nvalues % 10 == 0:
                # mixed ints, doubles
                nrows = get_table_size_from_ncolumns('GPDT', nvalues, 10)
                iints = [0, 1, 8, 9]
                #ifloats = [2, 3, 4, 5, 6, 7]
                idoubles = [1, 2, 3]
                unused_izero = [1, 8, 9]

                #print('ints:')
                #print(ints)
                # nid cp, x,    y,    z,    cd, ps
                # [0, 1,  2, 3, 4, 5, 6, 7, 8,   9]
                # [ ,  ,  1, 1, 2, 2, 3, 3,  ,    ]
                if read_mode == 1:
                    ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 10).copy()
                    #floats = np.frombuffer(data, op2.fdtype).reshape(nrows, 10).copy()
                    doubles = np.frombuffer(data, 'float64').reshape(nnodes, 5).copy()

                    nid_cp_cd_ps = ints[:, iints]
                    xyz = doubles[:, idoubles]
            else:
                raise NotImplementedError(nvalues)

        op2.op2_results.gpdt = GPDT(nid_cp_cd_ps, xyz)
    else:
        skip_record() # nid,cp,x,y,z,cd,ps


    # 1. Scalar points are identified by CP=-1 and words X1 through
    #    PS are zero.
    # 3. or fluid grid points, CD=-1.
    #print(nid_cp_cd_ps)
    #print(xyz)


    isubtable = -4
    markers = op2_reader.get_nmarkers(1, rewind=True)

    if markers[0] != isubtable:
        op2_reader.read_markers([markers[0], 1, 0, 0])
        op2_reader.show(200)
        log.error('unexpected GPDT marker marker=%s; expected=%s' % (
            markers[0], isubtable))
        #markers = op2_reader.get_nmarkers(1, rewind=False)
        return

    op2_reader.read_3_markers([isubtable, 1, 0])
    markers = op2_reader.get_nmarkers(1, rewind=True)
    while markers[0] != 0:
        #markers = op2_reader.get_nmarkers(1, rewind=True)
        #op2_reader.log.debug('GPDT record; markers=%s' % str(markers))
        if read_mode == 1:
            self._skip_record()
        else:
            #log.debug('unexpected GPDT record; markers=%s' % str(markers))
            data = op2_reader._read_record()
            #print('read_mode=%s freqs=%s' % (read_mode, freqs.tolist()))

        markers = op2_reader.get_nmarkers(1, rewind=True)
        op2_reader.read_3_markers([isubtable, 1, 0])
        markers = op2_reader.get_nmarkers(1, rewind=True)
        isubtable -= 1
    del isubtable
    op2_reader.read_markers([0])

def read_bgpdt(op2_reader: OP2Reader) -> None:
    """
    reads the BGPDT, BGPDTS, BGPDTOLD tables

    tested by TestOP2Matrix.test_gpspc

    """
    #if read_mode == 1:
        #read_record = op2_reader._skip_record
    #else:
    read_record = op2_reader._read_record

    op2: OP2 = op2_reader.op2
    read_mode = op2_reader.read_mode
    size = op2_reader.size
    endian = op2_reader._endian

    table_name = op2_reader._read_table_name(rewind=False)
    op2.table_name = table_name
    #op2_reader.log.debug('table_name = %r' % table_name)
    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('read_bgpdt - %s\n' % table_name)

    op2_reader.read_markers([-1])
    header_data = op2_reader._read_record()  # (105, 51, 0, 0, 0, 0, 0)
    header_ints = np.frombuffer(header_data, op2.idtype8)
    header_floats = np.frombuffer(header_data, op2.fdtype8)
    #print('bgpdt ints =', header_ints)
    #print('bgpdt floats =', header_floats)

    #seid = ints[0] # ??? is this a table number?
    nnodes = header_ints[1]  # validated
    #print('nnodes =', nnodes)
    #print(ints)
    if op2_reader.is_debug_file:
        op2_reader.binary_debug.write('---markers = [-1]---\n')
    #print('--------------------')

    op2_reader.read_3_markers([-2, 1, 0])
    op2_reader.read_table_name(['BGPDT', 'BGPDTS', 'BGPDTOLD', 'BGPDTOUT'])

    #print('--------------------')

    op2_reader.read_3_markers([-3, 1, 0])

    #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\s402_sphere_03.op2
    #GRID 1 0 0.0    0.0 0.0      0
    #GRID 2 0 0.0    0.0 0.0      0
    #GRID 3 0 0.0871 0.0 -0.99619 0
    #D = (0.0, 5e-324,
    #     5e-324, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5e-323,
    #     1e-323, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.4e-323,
    #     1.5e-323, 3e-322, 0.0, 0.0, 0.0871, 0.0, -0.99619)
    #L = (0, 1,
    #     1, 61, 0, 0, 0, 0, 0, 0, 7,
    #     2, 61, 0, 0, 0, 0, 0, 0, 13,
    #     3, 61, 0, 0, 0.0871, 0, -0.99619)

    #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\s402_flxslddriver_05.op2
    #doubles (float64) = (0.0, 5e-324, 5e-324, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5e-323, 1e-323, 3e-322, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 6.4e-323, 1.5e-323, 3e-322, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 9.4e-323, 2e-323, 3e-322, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 1.24e-322, 2.5e-323, 3e-322, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 1.53e-322, 3e-323, 3e-322, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 1.83e-322, 3.5e-323, 3e-322, 0.0, 0.0, 60.0, 0.0, 0.0, 0.0, 2.1e-322, 4e-323, 3e-322, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 2.4e-322, 4.4e-323, 3e-322, 0.0, 0.0, 80.0, 0.0, 0.0, 0.0, 2.7e-322, 5e-323, 3e-322, 0.0, 0.0, 90.0, 0.0, 0.0, 0.0, 3e-322, 5.4e-323, 3e-322, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 3.3e-322, 6e-323, 3e-322, 0.0, 0.0, 5.0, 0.0, 20.0, 0.0, 3.6e-322, 6.4e-323, 3e-322, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 3.9e-322, 7e-323, 3e-322, 0.0, 0.0, 10.0, 0.0, 0.0)
    #long long (int64) = (
        #0, 1,
        #1, 61, 0, 0, 0, 0, 0, 0, 7,
        #2, 61, 0, 0, 4621819117588971520, 0, 0, 0, 13,
        #3, 61, 0, 0, 4626322717216342016, 0, 0, 0, 19,
        #4, 61, 0, 0, 4629137466983448576, 0, 0, 0, 25,
        #5, 61, 0, 0, 4630826316843712512, 0, 0, 0, 31,
        #6, 61, 0, 0, 4632233691727265792, 0, 0, 0, 37,
        #7, 61, 0, 0, 4633641066610819072, 0, 0, 0, 43,
        #8, 61, 0, 0, 4634626229029306368, 0, 0, 0, 49,
        #9, 61, 0, 0, 4635329916471083008, 0, 0, 0, 55,
        #10, 61, 0, 0, 4636033603912859648, 0, 0, 0, 61,
        #11, 61, 0, 0, 4636737291354636288, 0, 0, 0, 67,
        #12, 61, 0, 0, 4617315517961601024, 0, 4626322717216342016, 0, 73,
        #13, 61, 0, 0, 4617315517961601024, 0, 0, 0, 79,
        #14, 61, 0, 0, 4621819117588971520, 0, 0)

    if read_mode == 1:
        op2_reader._skip_record()

    elif read_mode == 2:
        i = 0
        data = read_record() # cd,x,y,z
        # xword = 4 * op2_reader.factor
        nvalues4 = len(data) // 4
        assert len(data) % 4 == 0, len(data) % 4
        # assert len(data) % xword == 0

        numwide = nvalues4 // nnodes
        result_name = 'bgpdt'
        if op2._results.is_saved(result_name):
            op2._results._found_result(result_name)

            if size == 4:
                #C:\NASA\m4\formats\git\examples\x33_blog\blog_materials\pressure_vessel_fem1_sim1-my_limit_load.bdf
                #GRID*                  1               0-1.792343211E+003.4539249252E+00+
                #*       3.6302008282E-01               0
                #GRID*                  2               0-1.792313297E+003.3970433059E+00+
                #*       7.2220293673E-01               0
                #GRID*                  3               0-1.792337984E+003.3029109367E+00+
                #*       1.0733895770E+00               0
                #GRID*                  4               0-1.792341147E+003.1725968908E+00+
                #*       1.4128346960E+00               0
                #GRID*                  5               0-1.792333634E+003.0076863630E+00+
                #*       1.7364594482E+00               0

                assert nvalues4 % nnodes == 0, nvalues4 % nnodes
                #print(size, numwide, nvalues4, len(data))

                if numwide == 4:
                    assert numwide == 4, numwide
                    ntotal = 16 # 4*4
                    assert len(data) % ntotal == 0, len(data) % ntotal
                    structi = Struct(endian + b'i3f')  # 16
                    for j in range(nnodes):
                        out = structi.unpack(data[i:i+ntotal])
                        cd, x, y, z = out
                        outs = f'nid={j+1} cd={cd} xyz=({x:g},{y:g},{z:g})'
                        #print(outs)
                        i += ntotal
                    # [cd, x, y, z]
                    cd = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 4)[:, :-3]
                    xyz = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 4)[:, -3:]

                elif numwide == 12:
                    ntotal = 48 # (9+3)*4
                    structi = Struct(endian + b'6i3d')
                    # len(data) // 4 / ngrids = 4
                    # 30768 // 4 / 1923
                    nrows = get_table_size_from_ncolumns('BGPDT', nvalues4, 12)
                    for j in range(nnodes):
                        edata = data[i:i+ntotal]
                        #op2_reader.show_data(edata, 'ifqd')
                        out = structi.unpack(edata)

                        #op2_reader.show_data(data[i:i+ntotal], types='qd')
                        #print(out)
                        cd, sil, nid, sixtyone, ps, zero_c, x, y, z = out
                        outs = f'cd={cd} sil={sil} nid={nid} sixtyone={sixtyone} ps={ps} zero_c={zero_c} xyz=({x:g},{y:g},{z:g})'
                        assert nid > 0, outs
                        #print(outs)
                        #assert zero_a in [0, 225, 362, 499], out
                        #assert zero_c == 0, outs
                        #11: C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\gluedg01f.op2
                        assert sixtyone in [11, 12, 61], f'sixtyone={sixtyone} outs={outs}'
                        i += ntotal

                    # [cd, x, y, z]
                    #op2_reader.show_data(data, types='if', endian=None, force=False)
                    #[0, num, nid, 61, 0, 0, x, _, y, _, z, _]
                    ints = np.frombuffer(data, op2.idtype).reshape(nrows, 12).copy()
                    floats = np.frombuffer(data, 'float64').reshape(nrows, 6).copy()
                    cd = ints[:, 0]
                    xyz = floats[:, 1:]
                    #print(ints[:-6])
                    #print(xyz)
                else:  # pragma: no cover
                    raise NotImplementedError((size, numwide))
                op2.op2_results.bgpdt = BGPDT(cd, xyz)
            else:
                #bad = []
                #nvalues = len(data) // 4
                #for i in [2, 3, 6]: # 2-16 checked
                    #if nvalues % i != 0:
                        #bad.append(i)
                #if bad:
                    #print(nvalues, bad)
                #op2_reader.show_data(data, types='ifqd')
                #print(nvalues4)
                # 112 / 28.
                #nrows9 = nvalues4 // 18
                #nrows4 = nvalues4 // 8
                #nvalues_9 = nvalues4 % 18
                #nvalues_4 = nvalues4 % 8
                #print('nvalues_9', nvalues_9, nvalues_9)
                #print('nvalues_4', nvalues_4, nvalues_4)
                if numwide == 18: # nvalues_9 == 0 and nvalues_4 != 0:
                    #assert nvalues4 % 18 == 0, nvalues4 % 18
                    #print(nrows)
                    #i = np.arange(0, nrows)
                    #i1 = i * 4
                    #i2 = (i + 1) * 4
                    #zero_a, sil, nid, sixtyone, ps, zero_c, x, y, z
                    ints = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 9)[:, :-3]
                    floats = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 9)[:, -3:]
                    nid = ints[:, 2]
                    cd = ints[:, [4, 5]]
                    xyz = floats

                    #GRID           1       0    27.5    20.0     0.0       0
                    #GRID           2       0 25.30329 25.30329     0.0       0
                    #GRID           3       0    20.0    27.5     0.0       0
                    #GRID           4       0 14.69671 25.30329     0.0       0

                    #GRID           1       0    27.5    20.0     0.0       0
                    #GRID           2       025.3032925.30329     0.0       0
                    #GRID           3       0    20.0    27.5     0.0       0
                    #GRID           4       014.6967125.30329     0.0       0
                    #GRID           5       0    12.5    20.0     0.0       0
                    #GRID           6       014.6967114.69671     0.0       0
                    #GRID           7       0    20.0    12.5     0.0       0
                    #GRID           8       025.3032914.69671     0.0       0
                    #GRID           9       0 22.3122 7.30685     0.0       0
                    #GRID          10       017.8372532.93714     0.0       0
                    #GRID          11       07.22962217.75039     0.0       0
                    #GRID          12       032.86626 17.7885     0.0       0
                    #GRID          13       0 30.56189.511215     0.0       0
                    #GRID          14       024.9442132.70425     0.0       0
                    #GRID          15       032.4224632.42249     0.0       0
                    #GRID          16       032.6858324.92913     0.0       0
                    #GRID          17       09.68575330.94779     0.0       0
                    #GRID          18       07.90302224.50037     0.0       0
                    #GRID          19       09.0897749.672205     0.0       0
                    #GRID          20       015.524677.918942     0.0       0
                    #GRID          21       0     0.0    32.0     0.0       0
                    #GRID          22       0     0.0    24.0     0.0       0
                    #GRID          23       0     0.0    16.0     0.0       0
                    #GRID          24       0     0.0     8.0     0.0       0
                    #GRID          25       0     8.0     0.0     0.0       0
                    #GRID          26       0    16.0     0.0     0.0       0
                    #GRID          27       0    24.0     0.0     0.0       0
                    #GRID          28       0    32.0     0.0     0.0       0
                    #GRID          29       0    32.0    40.0     0.0       0
                    #GRID          30       0    24.0    40.0     0.0       0
                    #GRID          31       0    16.0    40.0     0.0       0
                    #GRID          32       0     8.0    40.0     0.0       0
                    #GRID          33       0    40.0     8.0     0.0       0
                    #GRID          34       0    40.0    16.0     0.0       0
                    #GRID          35       0    40.0    24.0     0.0       0
                    #GRID          36       0    40.0    32.0     0.0       0
                    #GRID          37       0     0.0     0.0     0.0       0
                    #GRID          38       0     0.0    40.0     0.0       0
                    #GRID          39       0    40.0     0.0     0.0       0
                    #GRID          40       0    40.0    40.0     0.0       0
                    #GRID          41       0    12.5    20.0    20.1       0
                    #(0, 1, 1, 61, 0, 0, 27.5, 20.0, 0.0)
                    #(0, 7, 2, 61, 0, 0, 25.30329, 25.30329, 0.0)
                    #(0, 13, 3, 61, 0, 0, 20.0, 27.5, 0.0)
                    #(0, 19, 4, 61, 0, 0, 14.69671, 25.30329, 0.0)
                    #(0, 25, 5, 61, 0, 0, 12.5, 20.0, 0.0)
                    #(0, 31, 6, 61, 0, 0, 14.69671, 14.69671, 0.0)
                    #(0, 37, 7, 61, 0, 0, 20.0, 12.5, 0.0)
                    #(0, 43, 8, 61, 0, 0, 25.30329, 14.69671, 0.0)
                    #(0, 49, 9, 61, 0, 0, 22.3122, 7.30685, 0.0)
                    #(0, 55, 10, 61, 0, 0, 17.83725, 32.93714, 0.0)
                    #(0, 61, 11, 61, 0, 0, 7.229622, 17.75039, 0.0)
                    #(0, 67, 12, 61, 0, 0, 32.86626, 17.7885, 0.0)
                    #(0, 73, 13, 61, 0, 0, 30.5618, 9.511215, 0.0)
                    #(0, 79, 14, 61, 0, 0, 24.94421, 32.70425, 0.0)
                    #(0, 85, 15, 61, 0, 0, 32.42246, 32.42249, 0.0)
                    #(0, 91, 16, 61, 0, 0, 32.68583, 24.92913, 0.0)
                    #(0, 97, 17, 61, 0, 0, 9.685753, 30.94779, 0.0)
                    #(0, 103, 18, 61, 0, 0, 7.903022, 24.50037, 0.0)
                    #(0, 109, 19, 61, 0, 0, 9.089774, 9.672205, 0.0)
                    #(0, 115, 20, 61, 0, 0, 15.52467, 7.918942, 0.0)
                    #(0, 121, 21, 61, 0, 0, 0.0, 32.0, 0.0)
                    #(0, 127, 22, 61, 0, 0, 0.0, 24.0, 0.0)
                    #(0, 133, 23, 61, 0, 0, 0.0, 16.0, 0.0)
                    #(0, 139, 24, 61, 0, 0, 0.0, 8.0, 0.0)
                    #(0, 145, 25, 61, 0, 0, 8.0, 0.0, 0.0)
                    #(0, 151, 26, 61, 0, 0, 16.0, 0.0, 0.0)
                    #(0, 157, 27, 61, 0, 0, 24.0, 0.0, 0.0)
                    #(0, 163, 28, 61, 0, 0, 32.0, 0.0, 0.0)
                    #(0, 169, 29, 61, 0, 0, 32.0, 40.0, 0.0)
                    #(0, 175, 30, 61, 0, 0, 24.0, 40.0, 0.0)
                    #(0, 181, 31, 61, 0, 0, 16.0, 40.0, 0.0)
                    #(0, 187, 32, 61, 0, 0, 8.0, 40.0, 0.0)
                    #(0, 193, 33, 61, 0, 0, 40.0, 8.0, 0.0)
                    #(0, 199, 34, 61, 0, 0, 40.0, 16.0, 0.0)
                    #(0, 205, 35, 61, 0, 0, 40.0, 24.0, 0.0)
                    #(0, 211, 36, 61, 0, 0, 40.0, 32.0, 0.0)
                    #(0, 217, 37, 61, 0, 0, 0.0, 0.0, 0.0)
                    #(0, 223, 38, 61, 0, 0, 0.0, 40.0, 0.0)
                    #(0, 229, 39, 61, 0, 0, 40.0, 0.0, 0.0)
                    #(0, 235, 40, 61, 0, 0, 40.0, 40.0, 0.0)
                    #(0, 241, 41, 61, 0, 0, 12.5, 20.0, 20.1)
                    #doubles (float64) = (0, 7, 2, 61, 0.0, 0.0, 25.30329, 25.30329, 0.0)
                    #long long (int64) = (0, 7, 2, 61, 0, 0, 25.30329, 25.30329, 0)
                    #op2_reader.show_data(data[:80], types='qd')
                    ntotal = 72 # 9*8
                    structi = Struct(endian + b'6q3d')
                    for j in range(nnodes):
                        out = structi.unpack(data[i:i+ntotal])
                        #op2_reader.show_data(data[i:i+ntotal], types='qd')
                        #print(out)
                        cd, sil, nid, sixtyone, ps, zero_c, x, y, z = out
                        outs = f'cd={cd} sil={sil} nid={nid} sixtyone={sixtyone} ps={ps} zero_c={zero_c} xyz=({x:g},{y:g},{z:g})'
                        assert nid > 0, outs

                        assert zero_c == 0, outs
                        assert sixtyone == 61, outs
                        i += ntotal


                # elif nvalues_9 != 0 and nvalues_4 == 0:
                elif numwide == 8:
                    #assert numwide == 8, numwide
                    ntotal = 32 # 24+4 = 28
                    assert len(data) % ntotal == 0, len(data) % ntotal
                    structi = Struct(endian + b'q3d')  # 28
                    for j in range(nnodes):
                        out = structi.unpack(data[i:i+ntotal])
                        #op2_reader.show_data(data[i:i+ntotal], types='ifqd')
                        #print(out)
                        cd, x, y, z = out
                        i += ntotal
                    # [cd, x, y, z]
                    cd = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 4)[:, :-3]
                    xyz = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 4)[:, -3:]
                else:
                    raise RuntimeError((size, numwide))
                #cd = ints[::7]
                # [cd, x, _, y, _, z, _]

                #print(ints)
                #print(floats)
                #print(nrows*7, len(floats))
                #ints = ints[2:].reshape(nrows, 7)
                #floats = floats[2:].reshape(nrows, 7)
                #for inti, floati in zip(ints, floats):
                    #print(inti[:-3], floats[-3:])
                #print(xyz, xyz.shape)
                op2.op2_results.bgpdt = BGPDT(cd, xyz)
            #print('cd = %s' % cd.tolist())
            #print('xyz:\n%s' % xyz)

    op2_reader.read_3_markers([-4, 1, 0])
    marker = op2_reader.get_nmarkers(1, rewind=True)[0]
    if marker == 0:
        op2_reader.read_markers([0])
        return

    ## TODO: why is this needed??? (it is, but dmap is not clear)
    data = op2_reader._read_record()
    #op2_reader.show_data(data, types='i')
    isubtable = -5
    while 1:
        op2_reader.read_3_markers([isubtable, 1, 0])
        marker = op2_reader.get_nmarkers(1, rewind=True)[0]
        if marker == 0:
            break
        data = op2_reader._read_record()
        isubtable -= 1

    op2_reader.read_markers([0])

def _get_gpdt_nnodes_numwide(size: int, ndata: int,
                             header_ints,
                             log: SimpleLogger) -> tuple[int, int]:  # pragma: no cover
    """
    size=4; ndata=1120 [102  40   0   0   0   0   0] -> nnodes=40
      1120=(7*4)*4*10 -> 7 words; nnodes=40
      1120=28*40
    """
    unused_table_id, nnodes_nbytes, nwords, *other = header_ints
    #print("header_ints =", header_ints)
    nvalues = ndata // 4
    #nnodes_nbytes,
    if nwords == 0:
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\c402pen12f.op2
        # size=8 ndata=1568 nnodes=28
        #    1568 = 2* 4^2 * 7^2
        #
        # [102  18   0   0   0   0   0]
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\rbarthm1.op2
        # size=4; ndata=168
        #   168 = 4 * 2 * 3*7
        #if self.op2.table_name in [b'GPDT', b'GPDTS']:
            #is_nodes = True
            ##self.size = 4 168
            #self.log.warning(f'size={self.size} ndata={ndata}')
        #else:  # pragma: no cover
        try:
            is_nodes, numwide, nnodes = _get_gpdt_numwide_from_nodes_null_nwords(size, nnodes_nbytes, ndata)
        except Exception:
            is_nodes = True
            numwide = 0
            raise
        log.debug(f'ndata={ndata} numwide={numwide} nnodes={nnodes}; is_nodes={is_nodes}')

        if is_nodes:
            nnodes = nnodes_nbytes
            #nnodes = ndata // nwords // 4
            #assert nnodes == 3, nnodes
            numwide = nvalues // nnodes
        if numwide not in [7, 10, 14]:
            if size == 8:
                numwide = 14
            else:
                numwide = 7
            nnodes = ndata // (4 * numwide)
            assert ndata % (4 * numwide) == 0
            if numwide not in [7, 10, 14]:
                raise RuntimeError(f'numwide={numwide} must be 7, 10, or 14')
    else:
        nnodes = nnodes_nbytes
        numwide = nvalues // nnodes
        assert ndata % 4 == 0

    #print(f'nnodes={nnodes} numwide={numwide} error={nvalues % nnodes}')
    return nnodes, numwide


def _get_gpdt_numwide_from_nodes_null_nwords(size: bool, nbytes: int,  # pragma: no cover
                                             ndata: int) -> tuple[bool, int, int]:
    """
    size=4 ndata=392 nnodes_nbytes=28
      392 = 4 * 2. * 7^2
    """
    is_nodes = False
    if nbytes % 2: # or (size == 4 and nnodes_nbytes not in [28, 40]):
        is_nodes = True
        #self.log.warning('is_nodes = True')
        numwide = 0
        nnodes = 0
        return is_nodes, numwide, nnodes

    nnodes = ndata // nbytes
    #print('size =', size, ndata, nbytes)
    assert ndata % nbytes == 0

    nvalues = ndata // 4
    assert ndata % 4 == 0
    numwide = nvalues // nnodes
    assert nvalues % nnodes == 0, f'size={size} ndata={ndata} nvalues={nvalues} nnodes={nnodes}'
    #if (self.size == 4 and nnodes_nbytes in [28, 40]):
    if size == 4:
        assert ndata == numwide * nnodes * 4, f'size=4 ndata={ndata} numwide={numwide} nnodes={nnodes}'
        if numwide not in [7, 10]:
            is_nodes = True
    elif size == 8:
        assert ndata == numwide * nnodes * 4, f'size=8 ndata={ndata} numwide={numwide} nnodes={nnodes}'
        assert numwide == 14, numwide
        #if numwide not in [7, 10]:
            #is_nodes = True
    return is_nodes, numwide, nnodes

def _read_gpdt_8_14(op2: OP2, data: bytes, nnodes: int) -> tuple[np.ndarray, np.ndarray]:
    """
    (nid, cp, x, y, z, cd, 0)
    """
    endian = op2._endian

    i = 0
    ntotal = 7 * 8  #  6 ints, 3 floats
    # (nid, cp, x, y, z, cd, 0)
    structi = Struct(endian + b'qq 3d qq')

    for j in range(nnodes):
        edata = data[i:i+ntotal]
        out = structi.unpack(edata)
        #print(out)
        i += ntotal

    #   0   1  2  3  4  5   6
    # nid, cp, x, y, z, cd, ps
    iints = [0, 1, 5, 6]
    ifloats = [2, 3, 4]

    ints = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()[:, iints]
    floats = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz


def _read_gpdt_4_7(op2: OP2, data: bytes,
                   nnodes: int) -> tuple[np.ndarray, np.ndarray]:
    """
    (nid, cp, x, y, z, cd, 0)
    """
    i = 0
    ntotal = 7 * 4  #  6 ints, 3 floats
    # (nid, cp, x, y, z, cd, 0)
    endian = op2._endian
    structi = Struct(endian + b'ii 3f ii')

    #ntotal = 16
    #structi = Struct(endian + b'')
    #self.show_data(data, types='if')
    for j in range(nnodes):
        edata = data[i:i+ntotal]
        #self.show_data(edata, types='if')
        out = structi.unpack(edata)
        #print(out)
        nid, zero_a, x, y, z, cd, zero_b = out
        outs = f'nid={nid} zero_a={zero_a} xyz=({x},{y},{z}) cd={cd} zero_b={zero_b}'
        assert nid > 0, outs
        #assert zero_a == 0, (nid, zero_a, x, y, z, cd, zero_b)
        #assert zero_b == 0, (nid, zero_a, x, y, z, cd, zero_b)
        i += ntotal

    #   0   1  2  3  4  5   6
    # nid, cp, x, y, z, cd, ps
    iints = [0, 1, 5, 6]
    ifloats = [2, 3, 4]

    ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 7).copy()[:, iints]
    floats = np.frombuffer(data, op2.idtype).reshape(nnodes, 7).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz

def _read_gpdt_4_10(op2: OP2, data: bytes,
                    nnodes: int) -> tuple[np.ndarray, np.ndarray]:
    """
    (nid, 0, x, y, z, cd, 0)
    """
    i = 0
    ntotal = 40  #  6 ints, 3 floats
    endian = op2._endian
    # (nid, 0, x, y, z, cd, 0)
    structi = Struct(endian + b'2i 3d 2i')
    for j in range(nnodes):
        edata = data[i:i+ntotal]
        #self.show_data(edata, types='ifqd')
        out = structi.unpack(edata)
        nid, zero_a, x, y, z, cd, zero_b = out
        outs = f'nid={nid} zero_a={zero_a} xyz=({x},{y},{z}) cd={cd} zero_b={zero_b}'
        #assert zero_a == 0, (nid, zero_a, x, y, z, cd, zero_b)
        #assert zero_b == 0, (nid, zero_a, x, y, z, cd, zero_b)
        i += ntotal
        #print(outs)
        assert nid > 0, outs

    #   0   1    2  3  4  5  6  7  8   9
    # nid, zero, x, _, y, _, z, _, cd, zero
    iints = [0, 1, 8, 9]

    # 0  1  2  3  4
    # a, x, y, z, b
    ifloats = [1, 2, 3]

    ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 10).copy()[:, iints]
    floats = np.frombuffer(data, 'float64').reshape(nnodes, 5).copy()[:, ifloats]
    nid_cp_cd_ps = ints
    xyz = floats
    return nid_cp_cd_ps, xyz


def _get_gpdt_nnodes2(ndata: int, header_ints: list[int],
                      size: int) -> tuple[int, int]:
    unused_table_id, nnodes_nbytes, nwords, *other = header_ints
    nvalues = ndata // 4
    assert nvalues > 0
    assert ndata % 4 == 0
    try:
        # assume nodes
        nnodes = nnodes_nbytes
        numwide = nvalues // nnodes
        assert nvalues % nnodes == 0
        if size == 4:
            assert numwide in [7, 10], numwide
        else:
            assert numwide == 14, numwide
    except AssertionError:
        # calculate the bytes
        if size == 4:
            numwide = 7
        elif numwide == 8:
            numwide = 14
        nnodes = nvalues // numwide
    assert ndata == nnodes * numwide * 4
    return nnodes, numwide

def get_table_size_from_ncolumns(table_name: str, nvalues: int,
                                 ncolumns: int) -> int:
    nrows = nvalues // ncolumns
    if nvalues % ncolumns != 0:
        msg = 'table=%s: nrows=nvalues/ncolumns=%s/%s=%s; nrows=%s must be an int' % (
            table_name, nvalues, ncolumns, nrows, nvalues / ncolumns)
        raise RuntimeError(msg)
    return nrows
