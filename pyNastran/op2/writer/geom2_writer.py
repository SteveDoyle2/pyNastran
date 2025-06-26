from __future__ import annotations

#import struct
from collections import defaultdict
from struct import pack, Struct
from typing import BinaryIO, TYPE_CHECKING

from pyNastran.op2.errors import SixtyFourBitError
from .geom1_writer import write_geom_header, close_geom_table
integer_types = int
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom


def write_geom2(op2_file: BinaryIO, op2_ascii,
                obj, endian=b'<',
                nastran_format: str='nx'):
    if not hasattr(obj, 'elements'):
        return
    #if not hasattr(obj, 'nodes'):
        #return

    nspoints = len(obj.spoints)
    nplotels = len(obj.plotels)
    nelements = len(obj.elements)
    nmic = len(obj.micpnt)
    nbc = len(obj.bcs)
    nall = nelements + nplotels + nspoints + nmic + nbc
    if nall == 0:
        return
    cards_written = {}
    write_geom_header(b'GEOM2', op2_file, op2_ascii)
    itable = -3

    log = obj.log
    # etypes_to_skip = [
    #     'CHBDYE', 'CBEND',
    #     #'CHBDYP',
    # ]
    out = defaultdict(list)
    for eid, element in obj.elements.items():
        out[element.type].append(eid)
    for eid, element in obj.masses.items():  # CONM2
        out[element.type].append(eid)
    for eid, element in obj.micpnt.items():
        out[element.type].append(eid)
    for bc_id, bcs in obj.bcs.items():
        log.error(str(bcs))
        for bc in bcs:
            out[bc.type].append(bc)

    if nspoints:
        out['SPOINT'] = list(obj.spoints.keys())
    if nplotels:
        out['PLOTEL'] = list(obj.plotels.keys())

    # elements with fixed lengths
    geom2_key_mapper = {
        # key, spack, nfields
        'CHBDYP' : ((10908, 109, 407), b'12i 3f', 15),
        'CHBDYG' : ((10808, 108, 406), b'16i', 16),
        'PLOTEL' : ((5201, 52, 11), b'3i', 3),
        'CTUBE' : ((3701, 37, 49), b'4i', 4),
        'CSHEAR' : ((3101, 31, 61), b'6i', 6),
        'CQUAD4' : ((2958, 51, 177), b'6iffii4f', 14),
        'CTRIA3' : ((5959, 59, 282), b'5iff3i3f', 13),
        'CQUADR' : ((8009, 80, 367), b'6iffii4f', 14),  # same as CQUAD4
        'CTRIAR' : ((9200, 92, 385), b'5iff3i3f', 13),  # same as CTRIA3
        'CQUAD8' : ((4701, 47, 326), b'10i 6f i', 17),  # current; not 2001
        'CTRIA6' : ((4801, 48, 327), b'8i 5f i', 14),  # current; not 2001
        'CTRIAX' : ((10108, 101, 512), b'9i', 9),
        'CTRIAX6' : ((6108, 61, 107), b'8i f ii', 11),
        'CQUAD' : ((9108, 91, 507), b'11i', 11),
        'CQUADX' : ((9008, 90, 508), b'11i', 11),  # same as CQUAD
        'CROD' : ((3001, 30, 48), b'4i', 4),
        'CONROD' : ((1601, 16, 47), b'4i4f', 8),

        'CDAMP1' : ((201, 2, 69), b'6i', 6),
        'CDAMP2' : ((301, 3, 70), b'if4i', 6),
        'CDAMP3' : ((401, 4, 71), b'4i', 4),
        'CDAMP4' : ((501, 5, 72), b'ifii', 4),
        'CDAMP5' : ((10608, 106, 404), b'ifii', 4),

        'CELAS1' : ((601, 6, 73), b'6i', 6),
        'CELAS2' : ((701, 7, 74), b'if4iff', 8),
        'CELAS3' : ((801, 8, 75), b'4i', 4),
        'CELAS4' : ((901, 9, 76), b'ifii', 4),

        'CMASS1' : ((1001, 10, 65), b'6i', 6),
        'CMASS2' : ((1101, 11, 66), b'if4i', 6),
        'CMASS3' : ((1201, 12, 67), b'4i', 4),
        'CMASS4' : ((1301, 13, 68), b'ifii', 4),

        'CVISC' : ((3901, 39, 50), b'4i', 4),
        'CTRAX3' : ((6111, 61, 996), b'5if', 6),
        'CQUADX4' : ((6112, 61, 997), b'6if', 7),
        'CQUADX8' : ((6114, 61, 999), b'10if', 11),
        'CTRAX6' : ((6113, 61, 998), b'8if', 9),

        # masses :
        'CONM1' : ((1401, 14, 63), b'3i 21f', 24),
        'CONM2' : ((1501, 15, 64), b'3i 10f', 13),
        'CPLSTS3': ((8801, 88, 984), b'6i f 4i i3f i', 16),
        'CPLSTS4': ((8401, 84, 985), b'6i f 4i i4f', 16),
        'CPLSTS6': ((1801, 18, 986), b'2i 8i fi 4f 4i 4f', 24),
        'CPLSTS8': ((3601, 36, 987), b'2i 8i fi 4f 4i 4f', 24),

        'CPLSTN6': ((5801, 58, 982), b'2i 6i f 7i', 16),
        'CPLSTN8': ((7201, 72, 983), b'2i 8i f 5i', 16),

        'RADBC': ((12801, 128, 417), b'ifii', 4),
    }
    allowed_etypes = list(geom2_key_mapper) + list(GEOM2_MAP)
    for name, eids in sorted(out.items()):
        nelements = len(eids)
        # if name in etypes_to_skip:

        eid0 = eids[0]
        if isinstance(eid0, integer_types):
            max_eid_id = max(eids)
            if max_eid_id > 99999999:
                raise SixtyFourBitError(f'64-bit OP2 writing is not supported; {name}: max eid={max_eid_id}')
        #if max_nid > 99999999:
            #raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max SPC nid={max_nid}')

        #if nelements == 0:
            #continue
        #if name not in etypes:
            #obj.log.warning('skipping GEOM2-%s' % name)
            #continue

        cards_written[name] = nelements
        if name in GEOM2_MAP:
            func = GEOM2_MAP[name]
            itable = func(obj, name, eids, nelements, itable,
                          op2_file, op2_ascii, endian)
            continue
        #elif name == 'CFAST':
            #_write_cfast(obj, name, eids, nelements, itable, op2_file, op2_ascii, endian,
                         #nastran_format=nastran_format)

        if name in geom2_key_mapper:
            key, spacki, nfields = geom2_key_mapper[name]
            spack = Struct(endian + spacki)
            #print(name, spacki)
        elif name == 'CBUSH':
            key = (2608, 26, 60)
            spack = None
            nfields = 14
        elif name == 'CBUSH1D':
            key = (5608, 56, 218)
            spack = Struct(endian + b'8i')
            nfields = 8
        elif name == 'CGAP':
            key = (1908, 19, 104)
            spack = None
            nfields = 9
        elif name == 'SPOINT':
            key = (5551, 49, 105)
            spack = None
            nfields = 1
        # -------------------
        # masses
        #elif name == 'CONM2':
            #key = (1501, 15, 64)
            #spack = None
            #nfields = 13

        # -------------------
        else:
            log.warning(f'skipping GEOM2-{name}')
            del cards_written[name]
            continue
        #else:  # pragma: no cover
            #raise NotImplementedError(name)

        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        nbytes = _write_intermediate_block(name, key, nfields, nelements,
                                           op2_file, op2_ascii)

        try:
            write_card(name, eids, spack, obj, op2_file, op2_ascii, endian)
        except Exception:
            obj.log.error('failed GEOM2-%s' % name)
            raise
        itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2_file, op2_ascii, itable)
    obj.log.debug(str(cards_written))

    #-------------------------------------

def _write_intermediate_block(name: str, key: tuple[int, int, int],
                              nfields: int, nelements: int,
                              op2_file: BinaryIO, op2_ascii):
    """writes the start of the geometry block; goes in the middle of the writer"""
    nvalues = nfields * nelements + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2_file.write(pack('3i', *[4, nvalues, 4]))
    op2_file.write(pack('i', nbytes)) #values, nbtyes))

    op2_file.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes

def _write_end_block(nbytes, itable, op2_file, op2_ascii):
    """closes off the geometry block"""
    op2_file.write(pack('i', nbytes))

    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2_file.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
    return itable

def write_cbeam(obj, name, eids, nelements, itable, op2_file, op2_ascii, endian):
    """writes the CBEAM"""
    key = (5408, 54, 261)
    #spack = None
    nfields = 18
    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2_file, op2_ascii)

    s1 = Struct(endian + b'6i3f3i6f')
    s3 = Struct(endian + b'12i6f')
    for eid in sorted(eids):
        elem = obj.elements[eid]
        ga, gb = elem.node_ids
        pid = elem.pid

        # per DMAP: F = FE bit-wise AND with 3
        #f = fe & 3
        w1a, w2a, w3a = elem.wa
        w1b, w2b, w3b = elem.wb
        pa = elem.pa
        pb = elem.pb
        sa = elem.sa
        sb = elem.sb
        if elem.g0 is None:
            x1, x2, x3 = elem.x
            fe = 0
            #(eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
             #pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
            data = [
                eid, pid, ga, gb, sa, sb, x1, x2, x3, fe,
                pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            op2_file.write(s1.pack(*data))
        else:
            fe = 2
            g0 = elem.g0
            #(eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
            # pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
            data = [
                eid, pid, ga, gb, sa, sb, g0, 0, 0, fe,
                pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            op2_file.write(s3.pack(*data))

    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    return itable

def write_cbar(obj, name, eids, nelements, itable, op2_file, op2_ascii, endian):
    """writes the CBAR"""
    key = (2408, 24, 180)
    #spack = None
    nfields = 16
    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2_file, op2_ascii)

    s1 = Struct(endian + b'4i3f3i6f')
    s3 = Struct(endian + b'7ii2i6f')
    for eid in sorted(eids):
        elem = obj.elements[eid]
        ga, gb = elem.node_ids
        pid = elem.pid

        # per DMAP: F = FE bit-wise AND with 3
        #f = fe & 3
        w1a, w2a, w3a = elem.wa
        w1b, w2b, w3b = elem.wb
        pa = elem.pa
        pb = elem.pb
        if elem.g0 is None:
            x1, x2, x3 = elem.x
            fe = 0
            #(eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
             #w1a, w2a, w3a, w1b, w2b, w3b) = out; fe=0
             #(eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
              #w1a, w2a, w3a, w1b, w2b, w3b) = out; fe=1
            data = [
                eid, pid, ga, gb, x1, x2, x3, fe, pa, pb,
                w1a, w2a, w3a, w1b, w2b, w3b, ]
            assert None not in data, 'CBAR-1; data=%s' % (data)
            #print('CBAR data1 =', data)
            op2_file.write(s1.pack(*data))
        else:
            fe = 2
            g0 = elem.g0
            #(eid, pid, ga, gb, g0, junk, junk, _f, pa,
             #pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
            data = [
                eid, pid, ga, gb, g0, 0, 0, fe, pa, pb,
                w1a, w2a, w3a, w1b, w2b, w3b]
            assert None not in data, 'CBAR-1; data=%s' % (data)
            #print('CBAR data2 =', data)
            op2_file.write(s3.pack(*data))

        #if f == 0:
            #out = s1.unpack(edata)
            #(eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
             #w1a, w2a, w3a, w1b, w2b, w3b) = out
            #data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
        #elif f == 1:
            #out = s2.unpack(edata)
            #(eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
             #w1a, w2a, w3a, w1b, w2b, w3b) = out
            #data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
        #elif f == 2:
            #out = s3.unpack(edata)
            #(eid, pid, ga, gb, g0, junk, junk, _f, pa,
             #pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
            #data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                        #w2a, w3a, w1b, w2b, w3b], [f, g0]]
        #else:
            #raise RuntimeError('invalid f value...f=%s' % (f))
        op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, ga, gb))

    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    return itable

def write_solid(model, name: str, eids: np.ndarray, nelements: int, itable: int,
                op2_file: BinaryIO, op2_ascii, endian: bytes) -> int:
    """writes the solid elements"""
    if name == 'CTETRA':
        key = (5508, 55, 217)
        nnodes = 10
        # 12 = eid, pid, n1, n2, n3, n4, ..., n10
    elif name == 'CHEXA':
        key = (7308, 73, 253)
        nnodes = 20
    elif name == 'CPENTA':
        key = (4108, 41, 280)
        nnodes = 15
    elif name == 'CPYRAM':
        key = (17200, 172, 1000)
        # it's 13, but there's a 14th node just because...
        nnodes = 14
    else:  # pragma: no cover
        raise NotImplementedError(name)
    nfields = nnodes + 2
    spack = Struct(endian + b'%ii' % nfields)

    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2_file, op2_ascii)
    for eid in sorted(eids):
        elem = model.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        if None in nids:
            nids = [nid if nid is not None else 0 for nid in nids]
        nnids = len(nids)
        if nnids < nnodes:
            nids2 = [0] * (nnodes - nnids)
            data = [eid, pid] + nids + nids2
        else:
            data = [eid, pid] + nids
        #print(name, data)
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        op2_file.write(spack.pack(*data))

    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    return itable

def write_chbdyp(eids, spack, obj, op2_file, op2_ascii):
    surface_type_str_to_int = {
        'POINT' : 1,
        'LINE' : 2,
        'ELCYL' : 6,
        'FTUBE' : 7,
        'TUBE' : 10,
        #'AREA3' : ,
    }
    for eid in sorted(eids):
        elem = obj.elements[eid]
        pid = elem.pid
        #print(elem.get_stats())
        surface_type_int = surface_type_str_to_int[elem.surface_type]
        #(eid, pid, Type, iviewf, iviewb, g1, g2, g0, radmidf, radmidb,
         #dislin, ce, e1, e2, e3) = out
        nids = elem.node_ids

        g2 = 0 if elem.g2 is None else elem.g2
        dislin = 0 if elem.gmid is None else elem.gmid
        g0 = 0 if elem.g0 is None else elem.g0
        e1 = 0. if elem.e1 is None else elem.e1
        e2 = 0. if elem.e2 is None else elem.e2
        e3 = 0. if elem.e3 is None else elem.e3
        data = (eid, pid, surface_type_int, elem.iview_front, elem.iview_back,
                elem.g1, g2, g0, elem.rad_mid_front, elem.rad_mid_back,
                dislin, elem.ce, e1, e2, e3)
        #data = [eid, 0, surface_type_int,
                #elem.iview_front, elem.iview_back,
                #elem.rad_mid_front, elem.rad_mid_back, 0] + all_nids
        assert None not in data, data
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        op2_file.write(spack.pack(*data))

def write_chbdyg(eids, spack, obj, op2_file, op2_ascii):
    surface_type_str_to_int = {
        'REV' : 3,
        'AREA3' : 4,
        'AREA4' : 5,
        'AREA6' : 8,
        'AREA8' : 9,
    }
    for eid in sorted(eids):
        elem = obj.elements[eid]
        #print(elem.get_stats())
        nids = elem.node_ids
        #if None in nids:
            #nids = [nid if nid is not None else 0 for nid in nids]
        all_nids = [0] * 8
        nnodes = len(nids)
        all_nids[:nnodes] = nids
        assert None not in nids, nids
        surface_type_int = surface_type_str_to_int[elem.surface_type]
        #(eid, unused_blank, Type, iviewf, iviewb, radmidf, radmidb, unused_blank2,
         #g1, g2, g3, g4, g5, g6, g7, g8) = out
        data = [eid, 0, surface_type_int,
                elem.iview_front, elem.iview_back,
                elem.rad_mid_front, elem.rad_mid_back, 0] + all_nids
        assert None not in data, data
        op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
        op2_file.write(spack.pack(*data))

def write_cbush(eids, spack, obj, op2_file, op2_ascii, endian):
    spacki = Struct(endian + b'4i iii i ifi3f')
    spackf = Struct(endian + b'4i fff i ifi3f')
    for eid in sorted(eids):
        elem = obj.elements[eid]

        pid = elem.pid
        ga, gb = elem.node_ids
        s = elem.s
        s1, s2, s3 = elem.si
        cid = elem.cid
        ocid = elem.ocid
        if cid is None:
            cid = -1

        # not 100%
        s1 = 0.0 if s1 is None else s1
        s2 = 0.0 if s2 is None else s2
        s3 = 0.0 if s3 is None else s3

        if elem.x[0] is None and elem.g0 is None:
            # Use Element CID below for orientation
            f = -1
            data = [eid, pid, ga, gb, 0, 0, 0,
                    f, cid, s, ocid, s1, s2, s3]
            assert None not in data, ('CBUSH-1 %s' %
                                      (data))
            op2_file.write(spacki.pack(*data))
        elif elem.x[0] is not None:
            f = 0
            x1, x2, x3 = elem.x
            data = [eid, pid, ga, gb, x1, x2, x3,
                    f, cid, s, ocid, s1, s2, s3]
            assert None not in data, 'CBUSH-2 %s x=%s' % (data, elem.x)
            op2_file.write(spackf.pack(*data))
        elif elem.g0 is not None:
            f = 2
            g0 = elem.g0
            data = [eid, pid, ga, gb, g0, 0, 0,
                    f, cid, s, ocid, s1, s2, s3]
            assert None not in data, 'CBUSH-3 %s' % (data)
            op2_file.write(spacki.pack(*data))
        else:
            raise RuntimeError('invalid CBUSH')

def write_cbush1d(eids, spack, obj, op2_file, op2_ascii, endian):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        #(eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
        g1, g2 = elem.node_ids
        cid = elem.cid
        if cid is None:
            cid = -1
        data = [eid, elem.pid, g1, g2, cid, 0, 0, 0]
        op2_file.write(spack.pack(*data))

def write_cgap(eids, spack, obj, op2_file, op2_ascii, endian):
    structf = Struct(endian + b'4i3fii')
    structi = Struct(endian + b'4i3iii')
    for eid in sorted(eids):
        elem = obj.elements[eid]
        #(eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
        pid = elem.pid
        ga, gb = elem.node_ids
        cid = elem.cid
        #print(elem.get_stats())
        if cid is None:
            cid = -1

        if elem.x[0] is not None and elem.g0 is None:
            f = 1
            x1, x2, x3 = elem.x
            data = [eid, pid, ga, gb, x1, x2, x3, f, cid]
            op2_file.write(structf.pack(*data))
        elif elem.x[0] is None and elem.g0 is None:
            f = 1
            data = [eid, pid, ga, gb, 1., 0., 0., f, cid]
            op2_file.write(structf.pack(*data))
        elif elem.x[0] is not None:
            f = 1
            x1, x2, x3 = elem.x
            data = [eid, pid, ga, gb, x1, x2, x3, f, cid]
            #print('CGAP x; x=%s data=%s' % (elem.x, data))
            op2_file.write(structf.pack(*data))
        else:
            f = 2
            g0 = elem.g0
            data = [eid, pid, ga, gb, g0, 0, 0, f, cid]
            print('CGAP g0; x=%s gab0=%s data=%s' % (g0, [ga, gb, g0], data))
            op2_file.write(structi.pack(*data))

def write_cfast(model: OP2Geom,
                 name, eids, nelements, itable, op2_file, op2_ascii, endian,
                 nastran_format='nx') -> int:
    # MSC 2005r2 -> 2016
    # CFAST(9801,98,506) - the marker for Record 34
    # gs     : 1
    # ida    : 10
    # idb    : 12
    # pid    : 15
    # xs     : None
    # ys     : None
    # zs     : None
    formati = 0
    typei = 0
    cid = 0
    g_upper = [0] * 8
    g_lower = g_upper
    tavg = 0.0
    blank1 = 0
    blank2 = 0
    tmin = 0.0
    #                   up low thickness
    fmt = endian + b'8i 8i 8i f 2i f'
    spack = Struct(fmt)
    # 1 EID       I Element identification number
    # 2 PID       I Property identification number
    # 3 GS        I Spot weld master node identification numberGS
    # 4 FORMAT(C) I Connection format (0=gridid)
    # 5 GA        I ID of GA
    # 6 GB        I ID of GB
    # 7 TYPE      I Types of upper and lower elements for FORM="GRIDID"
    # 8 CID       I C
    # 9 GUPPER(8) I Grid IDs of the upper shell
    # FORMAT =0 GRIDID of GBI
    #   17 GLOWER(8) I
    # FORMAT =1 ALIGN (not used)
    #   17 GLOWER(8) I
    # FORMAT =2 ELEMID (not used)
    #   17 GLOWER(8) I
    # FORMAT =9 ELPAT for xyz
    #   17 XYZ(3) RS
    #   20 UNDEF(5) none Not used
    # FORMAT =10 PARTPAT for xyz
    #   17 XYZ(3) RS
    #   20 UNDEF(5) none Not used
    # End FORMAT
    # 25 TAVG RS Average shell thickness
    # 26 UNDEF(2) none Not used
    # 28 TMIN RS Minimum shell thickness

    nfields = 16
    key = ()
    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2_file, op2_ascii)
    for eid in sorted(eids):
        elem = model.elements[eid]
        print(elem.get_stats())
        data_in = [eid, elem.pid, elem.gs, formati, elem.ida, elem.idb, typei, cid] + g_upper + [tavg, blank1, blank2, tmin]
        #op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
        op2_file.write(spack.pack(*data_in))
    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    asdf
    return itable

def write_card(name, eids, spack, obj, op2_file, op2_ascii, endian):
    """writes the GEOM2 elements"""
    op2_ascii.write('GEOM2-%s\n' % name)

    eid0 = eids[0]
    if isinstance(eid0, integer_types):
        eid_max = max(eids)
        if eid_max > 99999999:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; {name} max(eid)={eid_max}')

    if name == 'CHBDYP':
        write_chbdyp(eids, spack, obj, op2_file, op2_ascii)
    elif name == 'CHBDYG':
        write_chbdyg(eids, spack, obj, op2_file, op2_ascii)

    elif name == 'PLOTEL':
        for eid in sorted(eids):
            elem = obj.plotels[eid]
            nids = elem.node_ids
            #(eid, n1, n2) = out
            data = [eid] + nids
            op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
            op2_file.write(spack.pack(*data))
    elif name == 'CBUSH':
        write_cbush(eids, spack, obj, op2_file, op2_ascii, endian)

    elif name == 'CBUSH1D':
        write_cbush1d(eids, spack, obj, op2_file, op2_ascii, endian)
    elif name == 'CGAP':
        write_cgap(eids, spack, obj, op2_file, op2_ascii, endian)

    elif name in ['CQUAD4', 'CQUADR']:
        write_cquad4_cquadr(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CQUAD8':  # current; not 2001
        write_cquad8(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CTRIA6':  # current; not 2001
        write_ctria6(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CPLSTS3':
        write_cplsts3(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CPLSTS4':
        write_cplsts4(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CPLSTS6':
        write_cplsts6(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CPLSTS8':
        write_cplsts8(eids, spack, obj, op2_file, op2_ascii, name)

    elif name == 'CPLSTN6':
        write_cplstn6(eids, spack, obj, op2_file, op2_ascii, name)
    elif name == 'CPLSTN8':
        write_cplstn8(eids, spack, obj, op2_file, op2_ascii, name)

    elif name == 'CTRIAX':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = [nid if nid is not None else 0
                    for nid in elem.node_ids]
            pid = elem.pid
            #eid, pid, n1, n2, n3, n4, n5, n6, unused_undef1 = data
            data = [eid, pid] + nids + [0]
            assert None not in data, '%s data=%s' % (name, data)
            #print('  CTRIAX eid=%s mid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2_file.write(spack.pack(*data))

    elif name == 'CTRIAX6':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = [nid if nid is not None else 0
                    for nid in elem.node_ids]
            mid = elem.mid
            #eid, mid, n1, n2, n3, n4, n5, n6, theta, unused_undef1, unused_undef2 = data
            data = [eid, mid] + nids + [elem.theta, 0, 0]
            assert None not in data, '%s data=%s' % (name, data)
            #print('  CTRIAX6 eid=%s mid=%s nids=%s data=%s\n' % (eid, mid, str(nids), data[6:]))
            op2_ascii.write('  eid=%s mid=%s nids=%s\n' % (eid, mid, str(nids)))
            op2_file.write(spack.pack(*data))

    elif name in ['CQUAD', 'CQUADX']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = [nid if nid is not None else 0
                    for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, n9) = out
            data = [eid, pid] + nids
            assert None not in data, '%s data=%s' % (name, data)
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2_file.write(spack.pack(*data))

    elif name in ['CTRIA3', 'CTRIAR']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            theta = get_theta_from_theta_mcid(elem.theta_mcid)
            t1 = elem.T1 if elem.T1 is not None else -1.
            t2 = elem.T2 if elem.T2 is not None else -1.
            t3 = elem.T3 if elem.T3 is not None else -1.
            #eid, pid, n1, n2, n3, theta_mcid, zoffs, blank1, blank2, tflag, t1, t2, t3
            data = [eid, pid] + nids + [theta, elem.zoffset, 0, 0,
                                        elem.tflag, t1, t2, t3]
            assert elem.tflag in [0, 1], elem.get_stats()
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2_file.write(spack.pack(*data))
    elif name in ['CTRAX3', 'CTRAX6', 'CQUADX4', 'CQUADX8']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            data = [eid, pid] + nids + [elem.theta]
            assert None not in data, '  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta)
            #print('  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta))
            op2_ascii.write('  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta))
            op2_file.write(spack.pack(*data))


    elif name in ['CROD', 'CTUBE', 'CVISC', 'CSHEAR']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            data = [eid, pid] + nids
            #print(data)
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2_file.write(spack.pack(*data))
    elif name == 'CONROD':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            #(eid, n1, n2, mid, a, j, c, nsm) = out
            data = [eid] + nids + [elem.mid, elem.A, elem.j, elem.c, elem.nsm]
            op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
            op2_file.write(spack.pack(*data))
    elif name in ['CELAS1', 'CDAMP1']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, g1, g2, c1, c2)
            data = [eid, pid, n1, n2, elem.c1, elem.c2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'CELAS2':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, k, g1, g2, c1, c2, ge, s) = out
            c2 = elem.c2 if elem.c2 is not None else 0
            data = [eid, elem.k, n1, n2, elem.c1, c2, elem.ge, elem.s]
            #print('CELAS2', data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name in ['CELAS3', 'CDAMP3', 'CDAMP5']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, s1, s2) = out
            data = [eid, pid, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'CELAS4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, k, s1, s2) = out
            data = [eid, elem.k, n1, n2]
            #print(data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'CDAMP2':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, bdamp, g1, g2, c1, c2) = out
            c1 = elem.c1 if elem.c1 is not None else 0
            c2 = elem.c2 if elem.c2 is not None else 0
            data = [eid, elem.b, n1, n2, c1, c2]
            #print(name, data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'CDAMP4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, b, s1, s2) = out
            data = [eid, elem.b, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2_file.write(spack.pack(*data))

    elif name == 'CMASS1':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            pid = elem.pid
            #(eid, pid, g1, g2, c1, c2)
            gc = []
            gc1, gc2 = elem.grid_component()
            data = [eid, pid, gc1[0], gc2[0], gc1[1], gc2[1]]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s gc=[%s, %s]\n' % (eid, pid, gc1, gc2))
            op2_file.write(spack.pack(*data))
    elif name == 'CMASS2':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            #(eid, mass, g1, g2, c1, c2) = out
            gc1, gc2 = elem.grid_component()
            data = [eid, elem.mass, gc1[0], gc2[0], gc1[1], gc2[1]]
            assert None not in data, data
            #print(name, data)
            op2_ascii.write(f'  eid={eid} data={data}\n')
            op2_file.write(spack.pack(*data))
    elif name == 'CMASS3':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, g1, g2)
            data = [eid, pid, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'CMASS4':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, b, s1, s2) = out
            data = [eid, elem.mass, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2_file.write(spack.pack(*data))
    elif name == 'SPOINT':
        nids = eids
        nids.sort()
        spack = Struct('%ii' % len(nids))
        op2_ascii.write('  spoints%s\n' % str(nids))
        op2_file.write(spack.pack(*nids))
    #--------------------------------------
    # masses
    elif name == 'CONM1':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            mass_matrix = elem.mass_matrix
            m1  = mass_matrix[0, 0]   # M11
            m2a = mass_matrix[1, 0]   # M21
            m2b = mass_matrix[1, 1]   # M22
            m3a = mass_matrix[2, 0]   # M31
            m3b = mass_matrix[2, 1]   # M32
            m3c = mass_matrix[2, 2]   # M33
            m4a = mass_matrix[3, 0]   # M41
            m4b = mass_matrix[3, 1]   # M42
            m4c = mass_matrix[3, 2]   # M43
            m4d = mass_matrix[3, 3]   # M44
            m5a = mass_matrix[4, 0]   # M51
            m5b = mass_matrix[4, 1]   # M52
            m5c = mass_matrix[4, 2]   # M53
            m5d = mass_matrix[4, 3]   # M54
            m5e = mass_matrix[4, 4]   # M55
            m6a = mass_matrix[5, 0]   # M61
            m6b = mass_matrix[5, 1]   # M62
            m6c = mass_matrix[5, 2]   # M63
            m6d = mass_matrix[5, 3]   # M64
            m6e = mass_matrix[5, 4]   # M65
            m6f = mass_matrix[5, 5]   # M66

            data = (eid, elem.nid, elem.cid,
                    m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
                    m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f)
            assert None not in data, '%s data=%s' % (name, data)
            op2_ascii.write('  CONM2 eid=%s data=%s\n' % (eid, str(data)))
            op2_file.write(spack.pack(*data))

    elif name == 'CONM2':
        for eid in sorted(eids):
            elem = obj.masses[eid]
            # 3i 10f
            #(eid, g, cid, m, x1, x2, x3, i1, i2a, i2b, i3a, i3b, i3c) = out
            data = [eid, elem.nid, elem.cid, elem.mass] + list(elem.X) + list(elem.I)
            assert None not in data, '%s data=%s' % (name, data)
            #assert len(data) == 13, data
            op2_ascii.write('  CONM2 eid=%s data=%s\n' % (eid, str(data)))
            op2_file.write(spack.pack(*data))
    else:  # pragma: no cover
        raise NotImplementedError(name)


def write_cquad4_cquadr(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        # (eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
        # t1, t2, t3, t4) = out
        theta = get_theta_from_theta_mcid(elem.theta_mcid)
        tflag = elem.tflag
        # if tflag is None:
        # tflag =
        t1 = elem.T1 if elem.T1 is not None else -1.
        t2 = elem.T2 if elem.T2 is not None else -1.
        t3 = elem.T3 if elem.T3 is not None else -1.
        t4 = elem.T4 if elem.T4 is not None else -1.
        # assert t4 == -1.0, elem.T4
        data = [eid, pid] + nids + [theta, elem.zoffset, 0,
                                    tflag, t1, t2, t3, t4]
        assert tflag in [0, 1], elem.get_stats()
        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        # 6i ff ii 4f
        op2_file.write(spack.pack(*data))

def write_cquad8(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = [nid if nid is not None else 0
                for nid in elem.node_ids]
        pid = elem.pid
         #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
          #t3, t4, theta, zoffs, tflag) = out # current
        #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
        #t1, t2, t3, t4, theta, zoffs) = out  # cquad8; 2001
        theta = get_theta_from_theta_mcid(elem.theta_mcid)
        tflag = elem.tflag if elem.tflag is not None else 0
        t1 = elem.T1 if elem.T1 is not None else -1.
        t2 = elem.T2 if elem.T2 is not None else -1.
        t3 = elem.T3 if elem.T3 is not None else -1.
        t4 = elem.T4 if elem.T4 is not None else -1.
        data = [eid, pid] + nids + [t1, t2, t3, t4,
                                    theta, elem.zoffset, tflag]
        assert None not in data, '%s data=%s' % (name, data)
        assert isinstance(elem.tflag, int), elem.get_stats()
        assert elem.tflag in [-1, 0, 1], elem.get_stats()
        #print('  CQUAD8 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        op2_file.write(spack.pack(*data))

def write_ctria6(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = [nid if nid is not None else 0
                for nid in elem.node_ids]
        pid = elem.pid
        #(eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3, tflag) = out
        theta = get_theta_from_theta_mcid(elem.theta_mcid)
        t1 = elem.T1 if elem.T1 is not None else -1.
        t2 = elem.T2 if elem.T2 is not None else -1.
        t3 = elem.T3 if elem.T3 is not None else -1.
        data = [eid, pid] + nids + [t1, t2, t3,
                                    theta, elem.zoffset, elem.tflag]
        assert None not in data, '%s data=%s' % (name, data)
        assert elem.tflag in [-1, 0, 1], elem.get_stats()
        #print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        op2_file.write(spack.pack(*data))

def write_cplsts3(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        tflag = elem.tflag
        # if tflag is None:
        # tflag =
        t1 = elem.T1 if elem.T1 is not None else -1.
        t2 = elem.T2 if elem.T2 is not None else -1.
        t3 = elem.T3 if elem.T3 is not None else -1.
        # t4 = elem.T4 if elem.T4 is not None else -1.
        # assert t4 == -1.0, elem.T4
        # (eid, pid, n1, n2, n3, undef6, theta,
        #  undef8, undef9, undef10, undef11,
        # tflag, t1, t2, t3, undef16) = out
        data = [eid, pid] + nids + [0, theta, 0, 0, 0, 0,
                                    tflag, t1, t2, t3, 0]
        assert tflag in [0, 1], elem.get_stats()
        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        # 6i ff ii 4f
        op2_file.write(spack.pack(*data))


def write_cplsts4(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        tflag = elem.tflag
        # if tflag is None:
        # tflag =
        t1 = elem.T1 if elem.T1 is not None else -1.
        t2 = elem.T2 if elem.T2 is not None else -1.
        t3 = elem.T3 if elem.T3 is not None else -1.
        t4 = elem.T4 if elem.T4 is not None else -1.
        # assert t4 == -1.0, elem.T4
        # (eid, pid, n1, n2, n3, n4, theta,
        # undef8, undef9, undef10, undef11,
        #  tflag, t1, t2, t3, t4) = out
        data = [eid, pid] + nids + [theta, 0, 0, 0, 0,
                                    tflag, t1, t2, t3, t4]
        assert tflag in [0, 1], elem.get_stats()
        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        # 6i ff ii 4f
        op2_file.write(spack.pack(*data))

def write_cplsts6(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        tflag = 0 #elem.tflag

        t1 = t2 = t3 = -1.
        t4 = t5 = t6 = -1.
        # t1 = elem.T1 if elem.T1 is not None else -1.
        # t2 = elem.T2 if elem.T2 is not None else -1.
        # t3 = elem.T3 if elem.T3 is not None else -1.
        # t4 = elem.T4 if elem.T4 is not None else -1.
        # (eid, pid, n1, n2, n3, n4, n5, n6, undef7, undef8, theta,
        #  tflag, t1, t2, t3, zero1,
        #  zero2, zero3, zero4, zero5, t4, t5, t6, undef8b) = out
        #nids = [n1, n2, n3, n4, n5, n6]
        #undef = (zero1, zero2, zero3, zero4, zero5, undef7, undef8, undef8b)
        #b'2i 8i fi 4f 4i 4f'
        data = [eid, pid] + nids + [
            0, 0,
            theta, tflag,
            t1, t2, t3, 0.,
            0, 0, 0, 0,
            t4, t5, t6, 0.,
        ]
        assert tflag in [0, 1], elem.get_stats()
        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        print(spack.format)
        op2_file.write(spack.pack(*data))

def write_cplsts8(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        tflag = elem.tflag
        # t1 = elem.T1 if elem.T1 is not None else -1.
        # t2 = elem.T2 if elem.T2 is not None else -1.
        # t3 = elem.T3 if elem.T3 is not None else -1.
        # t4 = elem.T4 if elem.T4 is not None else -1.
        t1 = t2 = t3 = t4 = -1.
        t5 = t6 = t7 = t8 = -1.
        #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta,
        # tflag, t1, t2, t3, t4, zero1, zero2, zero3, zero4, t5, t6, t7, t8) = out
        data = [eid, pid] + nids + [theta, tflag,
                                    t1, t2, t3, t4,
                                    0, 0, 0, 0,
                                    t5, t6, t7, t8]
        assert tflag in [0, 1], elem.get_stats()
        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        op2_file.write(spack.pack(*data))

def write_cplstn6(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        #(eid, pid, n1, n2, n3, n4, n5, n6, theta, *undef) = out
        data = [eid, pid] + nids + [theta, 0, 0, 0, 0, 0, 0, 0]

        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        # 6i ff ii 4f
        op2_file.write(spack.pack(*data))

def write_cplstn8(eids, spack, obj, op2_file, op2_ascii, name: str):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        nids = elem.node_ids
        pid = elem.pid
        theta = elem.theta
        # t1 = elem.T1 if elem.T1 is not None else -1.
        # t2 = elem.T2 if elem.T2 is not None else -1.
        # t3 = elem.T3 if elem.T3 is not None else -1.
        # t4 = elem.T4 if elem.T4 is not None else -1.
        # assert t4 == -1.0, elem.T4
        #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, theta, *undef) = out
        data = [eid, pid] + nids + [theta, 0, 0, 0, 0, 0]

        # print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
        op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
        assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
        # 6i ff ii 4f
        op2_file.write(spack.pack(*data))


def get_theta_from_theta_mcid(theta_mcid: int | float) -> float:
    """the theta/mcid field is stored in a strange way"""
    if isinstance(theta_mcid, integer_types):
        theta = 512. * (theta_mcid + 1)
    else:
        theta = theta_mcid
    return theta

def write_radbc(model, name: str, elems, nelements: int, itable: int,
                op2_file: BinaryIO, op2_ascii, endian: bytes) -> int:
    """
    RADBC(12801,128,417)

    Word Name Type Description
    1 EID      I Element identification number
    2 FAMB    RS Radiation view factor between the face and the ambient point
    3 CNTRLND  I Control point for radiation boundary condition
    4 NODAMB   I
    """
    structi = Struct(endian + b'ifii')

    nelements = 0
    for elem in elems:
        nelements += len(elem.eids)

    nfields = 4
    key = (12801, 128, 417)
    nbytes = _write_intermediate_block(name, key, nfields, nelements,
                                       op2_file, op2_ascii)
    for elem in elems:
        famb = elem.famb
        cntrlnd = elem.cntrlnd
        nodamb = elem.nodamb
        for eid in elem.eids:
            data = [eid, famb, cntrlnd, nodamb]
            # print(data) # 'ifii'

            # print('  RADBC eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
            # op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            assert None not in data, data
            op2_file.write(structi.pack(*data))
    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    return itable


def write_micpnt(model, name: str, eids: np.ndarray, nelements: int, itable: int,
                 op2_file: BinaryIO, op2_ascii, endian: bytes) -> int:
    """
    RECORD – MICPNT(2801,28,630)
    Word Name Type Description
    1 EID I Element identification number
    2 GID I Fluid grid identification number
    3 DESC(12) CHAR4 Description - 48 characters maximum
    """
    key = (2801, 28, 630)
    nfields = 14
    struc = Struct(endian + b'2i 48s')
    nbytes = _write_intermediate_block(name, key, nfields, nelements,
                                       op2_file, op2_ascii)
    for eid in sorted(eids):
        elem = model.micpnt[eid]
        # eid, node_id, name_bytes = out
        name = elem.name
        name_bytes = (f'{name:<48}').encode('latin1')
        data = [eid, elem.nid, name_bytes]
        op2_ascii.write(f'  eid={eid} nid={elem.nid} nids={name}\n')
        datai = struc.pack(*data)
        assert len(datai) == nfields*4, len(datai)
        op2_file.write(struc.pack(*data))
    itable = _write_end_block(nbytes, itable, op2_file, op2_ascii)
    return itable

    op2: OP2Geom = self.op2
    #size = self.size
    struc = Struct(op2._endian + b'2i 48s')
    ntotal = 8 + 48
    nelements = (len(data) - n) // ntotal
    assert (len(data) - n) % ntotal == 0
    assert nelements > 0
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]  # 4*4
        out = struc.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('  MICPNT=%s\n' % str(out))
        #(eid,n1,n2) = out
        eid, node_id, name_bytes = out
        name = name_bytes.decode('latin1').rstrip()

GEOM2_MAP = {
    'CBAR': write_cbar,
    'CBEAM': write_cbeam,
    'CTETRA': write_solid,
    'CHEXA': write_solid,
    'CPENTA': write_solid,
    'CPYRAM': write_solid,
    'MICPNT': write_micpnt,
    'RADBC': write_radbc,
}
