from collections import defaultdict
from struct import pack, Struct

from pyNastran.op2.errors import SixtyFourBitError
from .geom1_writer import write_geom_header, close_geom_table
integer_types = int

def write_geom2(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'elements'):
        return
    #if not hasattr(obj, 'nodes'):
        #return
    nspoints = len(obj.spoints)
    nplotels = len(obj.plotels)
    nelements = len(obj.elements)
    if nelements == 0 and nplotels == 0 and nspoints == 0:
        return
    write_geom_header(b'GEOM2', op2, op2_ascii)
    itable = -3

    #etypes = [
        #'CROD', 'CONROD',
        #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        #'CTRIA3', 'CQUAD4',
        #'CTETRA', 'CHEXA', 'CPENTA',
    #]
    etypes_to_skip = [
        'CHBDYE', 'CBEND',
        #'CHBDYP',
    ]
    out = defaultdict(list)
    for eid, element in obj.elements.items():
        out[element.type].append(eid)
    if nspoints:
        out['SPOINT'] = list(obj.spoints.keys())
    if nplotels:
        out['PLOTEL'] = list(obj.plotels.keys())

    # elements with fixed lengths
    mapper = {
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

        'CVISC' : ((3901, 39, 50), b'4i', 4),
        'CTRAX3' : ((6111, 61, 996), b'5if', 6),
        'CQUADX4' : ((6112, 61, 997), b'6if', 7),
        'CQUADX8' : ((6114, 61, 999), b'10if', 11),
        'CTRAX6' : ((6113, 61, 998), b'8if', 9),
    }
    for name, eids in sorted(out.items()):
        nelements = len(eids)
        if name in etypes_to_skip:
            obj.log.warning('skipping GEOM2-%s' % name)
            continue

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

        if name in ['CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM']:
            itable = _write_solid(obj, name, eids, nelements, itable, op2, op2_ascii, endian)
            continue

        if name in mapper:
            key, spacki, nfields = mapper[name]
            spack = Struct(endian + spacki)
            #print(name, spacki)
        elif name == 'CBAR':
            itable = _write_cbar(obj, name, eids, nelements, itable, op2, op2_ascii, endian)
            continue
        elif name == 'CBEAM':
            itable = _write_cbeam(obj, name, eids, nelements, itable, op2, op2_ascii, endian)
            continue
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
        else:
            obj.log.warning('skipping %s' % name)
            continue
        #else:  # pragma: no cover
            #raise NotImplementedError(name)

        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        nbytes = _write_intermediate_block(name, key, nfields, nelements, op2, op2_ascii)

        try:
            write_card(name, eids, spack, obj, op2, op2_ascii, endian)
        except:
            obj.log.error('failed GEOM2-%s' % name)
            raise
        itable = _write_end_block(nbytes, itable, op2, op2_ascii)

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)

    #-------------------------------------

def _write_intermediate_block(name, key, nfields, nelements, op2, op2_ascii):
    """writes the start of the geometry block; goes in the middle of the writer"""
    nvalues = nfields * nelements + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes

def _write_end_block(nbytes, itable, op2, op2_ascii):
    """closes off the geometry block"""
    op2.write(pack('i', nbytes))

    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
    return itable

def _write_cbeam(obj, name, eids, nelements, itable, op2, op2_ascii, endian):
    """writes the CBEAM"""
    key = (5408, 54, 261)
    #spack = None
    nfields = 18
    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2, op2_ascii)

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
            op2.write(s1.pack(*data))
        else:
            fe = 2
            g0 = elem.g0
            #(eid, pid, ga, gb, sa, sb, g0, xxa, xxb, fe,
            # pa, pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
            data = [
                eid, pid, ga, gb, sa, sb, g0, 0, 0, fe,
                pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
            op2.write(s3.pack(*data))

    itable = _write_end_block(nbytes, itable, op2, op2_ascii)
    return itable

def _write_cbar(obj, name, eids, nelements, itable, op2, op2_ascii, endian):
    """writes the CBAR"""
    key = (2408, 24, 180)
    #spack = None
    nfields = 16
    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2, op2_ascii)

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
            op2.write(s1.pack(*data))
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
            op2.write(s3.pack(*data))

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

    itable = _write_end_block(nbytes, itable, op2, op2_ascii)
    return itable

def _write_solid(model, name, eids, nelements, itable, op2, op2_ascii, endian):
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
    spack = Struct(endian + b'%ii' % (nfields))

    nbytes = _write_intermediate_block(name, key, nfields, nelements, op2, op2_ascii)
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
        op2.write(spack.pack(*data))

    itable = _write_end_block(nbytes, itable, op2, op2_ascii)
    return itable

def _write_chbdyp(eids, spack, obj, op2, op2_ascii):
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
        op2.write(spack.pack(*data))

def _write_chbdyg(eids, spack, obj, op2, op2_ascii):
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
        op2.write(spack.pack(*data))

def _write_cbush(eids, spack, obj, op2, op2_ascii, endian):
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
            assert None not in data, 'CBUSH-1 %s' % (data)
            op2.write(spacki.pack(*data))
        elif elem.x[0] is not None:
            f = 0
            x1, x2, x3 = elem.x
            data = [eid, pid, ga, gb, x1, x2, x3,
                    f, cid, s, ocid, s1, s2, s3]
            assert None not in data, 'CBUSH-2 %s x=%s' % (data, elem.x)
            op2.write(spackf.pack(*data))
        elif elem.g0 is not None:
            f = 2
            g0 = elem.g0
            data = [eid, pid, ga, gb, g0, 0, 0,
                    f, cid, s, ocid, s1, s2, s3]
            assert None not in data, 'CBUSH-3 %s' % (data)
            op2.write(spacki.pack(*data))
        else:
            raise RuntimeError('invalid CBUSH')

def _write_cbush1d(eids, spack, obj, op2, op2_ascii, endian):
    for eid in sorted(eids):
        elem = obj.elements[eid]
        #(eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
        g1, g2 = elem.node_ids
        cid = elem.cid
        if cid is None:
            cid = -1
        data = [eid, elem.pid, g1, g2, cid, 0, 0, 0]
        op2.write(spack.pack(*data))

def _write_cgap(eids, spack, obj, op2, op2_ascii, endian):
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
            op2.write(structf.pack(*data))
        elif elem.x[0] is None and elem.g0 is None:
            f = 1
            data = [eid, pid, ga, gb, 1., 0., 0., f, cid]
            op2.write(structf.pack(*data))
        elif elem.x[0] is not None:
            f = 1
            x1, x2, x3 = elem.x
            data = [eid, pid, ga, gb, x1, x2, x3, f, cid]
            #print('CGAP x; x=%s data=%s' % (elem.x, data))
            op2.write(structf.pack(*data))
        else:
            f = 2
            g0 = elem.g0
            data = [eid, pid, ga, gb, g0, 0, 0, f, cid]
            print('CGAP g0; x=%s gab0=%s data=%s' % (g0, [ga, gb, g0], data))
            op2.write(structi.pack(*data))

def write_card(name, eids, spack, obj, op2, op2_ascii, endian):
    """writes the GEOM2 elements"""
    op2_ascii.write('GEOM2-%s\n' % name)

    eid_max = max(eids)
    if eid_max > 99999999:
        raise SixtyFourBitError(f'64-bit OP2 writing is not supported; {name} max(eid)={eid_max}')

    if name == 'CHBDYP':
        _write_chbdyp(eids, spack, obj, op2, op2_ascii)
    elif name == 'CHBDYG':
        _write_chbdyg(eids, spack, obj, op2, op2_ascii)

    elif name == 'PLOTEL':
        for eid in sorted(eids):
            elem = obj.plotels[eid]
            nids = elem.node_ids
            #(eid, n1, n2) = out
            data = [eid] + nids
            op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CBUSH':
        _write_cbush(eids, spack, obj, op2, op2_ascii, endian)

    elif name == 'CBUSH1D':
        _write_cbush1d(eids, spack, obj, op2, op2_ascii, endian)
    elif name == 'CGAP':
        _write_cgap(eids, spack, obj, op2, op2_ascii, endian)

    elif name in ['CQUAD4', 'CQUADR']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            #(eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
             #t1, t2, t3, t4) = out
            theta = get_theta_from_theta_mcid(elem.theta_mcid)
            tflag = elem.tflag
            #if tflag is None:
                #tflag =
            t1 = elem.T1 if elem.T1 is not None else -1.
            t2 = elem.T2 if elem.T2 is not None else -1.
            t3 = elem.T3 if elem.T3 is not None else -1.
            t4 = elem.T4 if elem.T4 is not None else -1.
            data = [eid, pid] + nids + [theta, elem.zoffset, 0,
                                        tflag, t1, t2, t3, t4]
            assert tflag in [0, 1], elem.get_stats()
            #print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            assert None not in data, '  %s eid=%s pid=%s nids=%s\n%s' % (name, eid, pid, str(nids), data)
            op2.write(spack.pack(*data))
    elif name == 'CQUAD8':  # current; not 2001
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
            op2.write(spack.pack(*data))
    elif name == 'CTRIA6':  # current; not 2001
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
            op2.write(spack.pack(*data))
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
            op2.write(spack.pack(*data))

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
            op2.write(spack.pack(*data))

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
            op2.write(spack.pack(*data))

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
            op2.write(spack.pack(*data))
    elif name in ['CTRAX3', 'CTRAX6', 'CQUADX4', 'CQUADX8']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            data = [eid, pid] + nids + [elem.theta]
            assert None not in data, '  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta)
            #print('  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta))
            op2_ascii.write('  eid=%s pid=%s nids=%s theta=%r\n' % (eid, pid, str(nids), elem.theta))
            op2.write(spack.pack(*data))


    elif name in ['CROD', 'CTUBE', 'CVISC', 'CSHEAR']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            data = [eid, pid] + nids
            #print(data)
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CONROD':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            #(eid, n1, n2, mid, a, j, c, nsm) = out
            data = [eid] + nids + [elem.mid, elem.A, elem.j, elem.c, elem.nsm]
            op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
            op2.write(spack.pack(*data))
    elif name in ['CELAS1', 'CDAMP1']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, g1, g2, c1, c2)
            data = [eid, pid, n1, n2, elem.c1, elem.c2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CELAS2':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, k, g1, g2, c1, c2, ge, s) = out
            c2 = elem.c2 if elem.c2 is not None else 0
            data = [eid, elem.k, n1, n2, elem.c1, c2, elem.ge, elem.s]
            #print('CELAS2', data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name in ['CELAS3', 'CDAMP3', 'CDAMP5']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            pid = elem.pid
            #(eid, pid, s1, s2) = out
            data = [eid, pid, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CELAS4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, k, s1, s2) = out
            data = [eid, elem.k, n1, n2]
            #print(data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
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
            op2.write(spack.pack(*data))
    elif name == 'CDAMP4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = [nid if nid else 0 for nid in elem.node_ids]
            #(eid, b, s1, s2) = out
            data = [eid, elem.b, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'SPOINT':
        nids = eids
        nids.sort()
        spack = Struct('%ii' % len(nids))
        op2_ascii.write('  spoints%s\n' % str(nids))
        op2.write(spack.pack(*nids))
    else:  # pragma: no cover
        raise NotImplementedError(name)

def get_theta_from_theta_mcid(theta_mcid):
    """the theta/mcid field is stored in a strange way"""
    if isinstance(theta_mcid, integer_types):
        theta = 512. * (theta_mcid + 1)
    else:
        theta = theta_mcid
    return theta
