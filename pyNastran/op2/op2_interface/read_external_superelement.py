"""
- read_cmodeext(op2_reader)
- read_cmodeext_helper(op2_reader)
"""
from __future__ import annotations
from struct import unpack, Struct # , error as struct_error
from typing import TYPE_CHECKING

import numpy as np

#from cpylog import SimpleLogger
#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.f06.errors import FatalError
#from pyNastran.op2.errors import FortranMarkerError, SortCodeError, EmptyRecordError

from pyNastran.op2.op2_interface.utils import (
    mapfmt, reshape_bytes_block,
)

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2


def read_cmodext(op2_reader: OP2Reader) -> None:
    r"""
    fails if a streaming block???:
     - nx_spike\mnf16_0.op2

    """
    op2: OP2 = op2_reader.op2
    size = op2_reader.size

    op2.table_name = op2_reader._read_table_name(rewind=False)
    #op2_reader.log.debug('table_name = %r' % op2.table_name)
    if op2.is_debug_file:
        op2.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)

    op2_reader.read_markers([-1])
    if op2.is_debug_file:
        op2.binary_debug.write('---markers = [-1]---\n')
    data = op2_reader._read_record()
    out = unpack(mapfmt(b'7i', size), data)
    op2_reader.log.debug(str(out))
    # (101, 500, 37232, 2, 1, 18355, 158)
    #op2_reader.show_data(data, types='sqd')

    markers = op2_reader.get_nmarkers(1, rewind=True)
    if op2.is_debug_file:
        op2.binary_debug.write('---marker0 = %s---\n' % markers)

    marker = -2
    markers = op2_reader.read_markers([marker, 1, 0])

    data = op2_reader._read_record()
    if size == 4:
        unused_table_name, oneseventy_a, oneseventy_b = unpack('8sii', data)
    else:
        table_name, oneseventy_a, oneseventy_b = unpack('16sqq', data)
        unused_table_name = reshape_bytes_block(table_name)
    assert oneseventy_a == 170, oneseventy_a
    assert oneseventy_b == 170, oneseventy_b
    #print('170*4 =', 170*4)
    #op2_reader.show_data(data)
    marker -= 1
    n = op2.n
    while marker < 0:
        #marker = op2_reader._read_cmodext_helper(marker) # -3
        marker, is_done = read_cmodext_helper(op2_reader, marker)
        #print(f'--end marker={marker}')
        #print(f'--end marker2={marker2}')
        if is_done:
            break
        marker2 = op2_reader.get_nmarkers(1, rewind=True)[0]
        n = op2.n
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    #marker = op2_reader._read_cmodext_helper(marker); print(marker)
    print(f'table end; {marker2-1}')
    #marker = op2_reader._read_cmodext_helper(marker, debug=True)
    #op2.show_ndata(200)
    #sss

def read_cmodext_helper(op2_reader: OP2Reader, marker_orig,
                        debug: bool=False) -> tuple[int, bool]:
    r"""

    64-bit:
      C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\mbdrecvr_c_0.op2
      C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cntlmtl05_0.op2
    """
    op2: OP2 = op2_reader.op2
    size = op2_reader.size
    log: SimpleLogger = op2.log

    marker = marker_orig
    #markers = op2_reader.read_nmarkers([marker, 1, 1]) # -3
    debug = False
    if debug:
        op2.show_ndata(100)
    markers = op2_reader.get_nmarkers(3, rewind=False)
    assert markers == [marker_orig, 1, 1], markers
    #print('markers =', markers)

    #marker = op2_reader.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
    val_old = 0
    if debug:
        print('-----------------------------')
    #i = 0
    #icheck = 7
    #factor = op2_reader.factor
    if size == 4:
        expected_marker = 3
        sint = 'i'
        ifs = 'ifs'
        #dq = 'if'
        structi = Struct(b'i')
        struct_qd = Struct(b'i f')
        struct_q2d = Struct(b'i 2f')
        struct_q3d = Struct(b'i 3f')
        struct_d4q = Struct(b'i 4f')
    else:
        expected_marker = 3
        sint = 'q'
        ifs = 'qds'
        #dq = 'dq'
        structi = Struct(b'q')
        struct_qd = Struct(b'q d')
        #struct_dq = Struct(b'd q')
        struct_q2d = Struct(b'q 2d')
        struct_q3d = Struct(b'q 3d')
        struct_d4q = Struct(b'd 4q')

    marker = op2_reader.get_nmarkers(1, rewind=True, macro_rewind=True)[0]
    #print('AAAAA', marker)
    if marker is None:
        asdf
    elif marker == 3:
        pass
    elif marker < 0:
        log.warning('finished CMODEXT :)!')
        done = True
        return None, done
    else:
        raise RuntimeError(marker)

    i = 0
    while 1:
        #print('------------------------------')
        #print('i = %i' % i)
        marker = op2_reader.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        if marker != expected_marker:
            log.info(f'CMODEXT i={i} marker={marker}')

        if size == 4:
            assert marker in [expected_marker], marker
        else:
            #if marker == -4:
                #op2_reader.show(200, types='isq', endian=None, force=False)
            assert marker in [1, 2, 3, 10, 13, 16], marker
        data = op2_reader._read_block()
        ndata = len(data)
        nvalues = ndata // size
        if nvalues == 2:
            #op2_reader.show_data(data, types=dq)
            #print(ndata)
            #op2_reader.show_data(data, types=ifs)
            a, b = struct_qd.unpack(data)
            print(f' 2 {i}: {a:-4d} {b:.6g}')
        elif nvalues == 3:
            #op2_reader.show_data(data, types=dq)
            #print(ndata)
            #op2_reader.show_data(data, types=ifs)
            a, b, c = struct_q2d.unpack(data)
            print(f' 3 {i}: {a:-4d} {b:.6g}')
        elif nvalues == 4:
            #op2_reader.show_data(data, types=dq)
            #print(ndata)
            #op2_reader.show_data(data, types=ifs)
            a, b, c, d = struct_q3d.unpack(data)
            #print(f'4 {i}: {a:-4d} {b:10g} {c:10g} {d:10g}')
        elif nvalues == 5:
            #op2_reader.show_data(data, types=dq)
            #print(ndata)
            out = struct_d4q.unpack(data)
            print('5 {i}:', out)
        #elif ndata == 24 * factor:
            #op2_reader.show_data(data, types=dq)
            #print(ndata)
            #op2_reader.show_data(data, types=ifs)
            #a, b, c, d = struct_q2d.unpack(data)
        elif ndata == 11:
            # ints    = (697, 0, 0, 1072693248, -427255041, -1128310260, -1065349346, 1019092170, -458662240, 1015082339, -1959380180, 1018327314, 863132413, -1129941844, 477962639, 992946870, -935437031, 992276971, -1286826876, -1158808745, 1765409376, 986144275)
            # floats  = (9.767050296343975e-43, 0.0, 0.0, 1.875, -3.2255050747205135e+23, -0.023358367383480072, -4.001845359802246, 0.023207087069749832, -2.4994545586342025e+22, 0.015738194808363914, -3.509644095891333e-32, 0.02178243175148964, 5.642776912395675e-08, -0.02031930536031723, 8.375405145620162e-22, 0.0026728338561952114, -194932.390625, 0.002516860840842128, -4.76325254794574e-08, -0.0018156570149585605, 1.4054463629294865e+25, 0.001521053141914308)
            # doubles (float64) = (3.444e-321, 1.0, -4.350930636102116e-16, 4.1789445167220847e-16, 2.936419423485894e-17, 2.559298496168014e-16, -1.5581808467238972e-16, 1.2890303437976399e-23, 8.662650619686277e-24, -7.750115449925671e-25, 1.5100882399748836e-25)
            # long long (int64) = (697, 4607182418800017408, -4846055662573544705, 4376967544989290270, 4359745452588490400, 4373682512589110060, -4853063265498801411, 4264674333793526159, 4261797142378470681, -4977045659085663100, 4235457412028039776)

            #ints    = (697, 0, 0, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z)
            #floats  = (697, 0.0,
                       #0.0, 1.875,
                       #-3.2255050747205135e+23, -0.023358367383480072, -4.001845359802246, 0.023207087069749832, -2.4994545586342025e+22, 0.015738194808363914, -3.509644095891333e-32, 0.02178243175148964, 5.642776912395675e-08, -0.02031930536031723, 8.375405145620162e-22, 0.0026728338561952114, -194932.390625, 0.002516860840842128, -4.76325254794574e-08, -0.0018156570149585605, 1.4054463629294865e+25, 0.001521053141914308)
            #doubles (float64) = (697, 1.0, y, y, y, y, y, y, y, y, y)
            #long long (int64) = (697, 1.0, x, x, x, x, x, x, x, x, x)
            op2_reader.show_data(data, types='ifsqd')
            #op2_reader.show_data(data[4:], types='qd')
            raise RuntimeError(f'marker={marker} ndata={ndata}; nvalues={nvalues}')
        elif nvalues == 17:
            inti = structi.unpack(data[:size])
            floats = np.frombuffer(data, dtype=op2.fdtype8)[1:]
            print(inti, floats)
        else:
            op2_reader.show_data(data, types=ifs)
            raise RuntimeError(f'marker={marker} ndata={ndata}; nvalues={nvalues}')
            #sdf
            #ints = [
                #1387, 0, -1574726656, -1076024976, 12958534, -1079706775, -1204798216, -1076795484, -674762419, -1125074255, -1234714250, 1025640630, 367990681, 1024314941, -1752243687, 1021642278,
                #-3, 1017741311,
                #-3, 1019878559,
                #-2, 1019293439,
                #-5, -1134034945,
                #-3, 1017733119,
                #-3, 1017241599]
            #floats = [
                #1387, 0.0, -2.2168969812013176e-18, -1.7278270721435547, 1.815877379410093e-38, -1.2889224290847778, -4.201845149509609e-05, -1.6359753608703613, -439678385782784.0, -0.029385896399617195, -3.45342914442881e-06, 0.03955908864736557, 9.656856178702571e-26, 0.034620512276887894, -9.233609993079022e-25, 0.02795703336596489,
                #nan, 0.02069091610610485,
                #nan, 0.024671850726008415,
                #nan, 0.023581979796290398,
                #nan, -0.014160155318677425,
                #nan, 0.02067565731704235,
                #nan, 0.01976012997329235]
            #op2_reader.show_data(data, types='ifsqd')
            #ints = np.frombuffer(data, dtype='int32').tolist()
            #floats = np.frombuffer(data, dtype='float32').tolist()
            #print(ints)
            #print(floats)
            #np.frombuffer(data, dtype='ifsdq')
            print('------------------------------')
            op2_reader.show_data(data, types=ifs)
        #else:
            #op2_reader.show_data(data, types=ifs)

        val = unpack(sint, data[:size])[0]
        if debug:
            print('val=%s delta=%s' % (val, val - val_old))
            #op2_reader.show_data(data, types=ifs)
            #op2_reader.show_data(data[4:], types=ifs)
        assert len(data) > size
        #print('i=%s val=%s delta=%s' % (i, val, val - val_old))
        val_old = val

        marker2 = op2_reader.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
        #print(marker2)
        if marker2 < 0:
            op2_reader.log.warning(f'breaking marker2={marker2}')
            op2_reader.show(300, types=ifs)
            break
        #if marker2 == 696:
            #op2_reader.log.warning('breaking 696')
            #break
        i += 1
    if debug:
        print('----------------------------------------')

    marker = op2_reader.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
    #if debug:
    print('****marker  = %s' % marker)
    #print('****marker2 = %s' % marker2)
    #assert marker == 696, marker
    #data = op2_reader._read_block()
    #op2_reader.show_data(data)

    #marker = op2_reader.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
    assert marker == (marker_orig - 1), f'marker={marker} marker_orig={marker_orig}'

    if debug:
        op2.show_ndata(200)
    done = False
    return marker, done

    #data = op2_reader._read_record()
    #marker -= 1
    #op2.show_ndata(100)

    ##marker -= 1
    ##marker_end = op2.get_marker1(rewind=False)

