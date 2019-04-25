from __future__ import print_function, absolute_import
from collections import defaultdict
from struct import pack, Struct

from .geom1 import write_geom_header, close_geom_table
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
        'CHBDYE', 'CHBDYG', 'CHBDYP', 'CBEND',
    ]
    out = defaultdict(list)
    for eid, element in obj.elements.items():
        out[element.type].append(eid)
    if nspoints:
        out['SPOINT'] = list(obj.spoints.keys())
    if nplotels:
        out['PLOTEL'] = list(obj.plotels.keys())

    for name, eids in sorted(out.items()):
        nelements = len(eids)
        if name in etypes_to_skip:
            obj.log.warning('skipping GEOM2-%s' % name)
            continue

        #if nelements == 0:
            #continue
        #if name not in etypes:
            #obj.log.warning('skipping GEOM2-%s' % name)
            #continue

        if name in ['CTETRA', 'CHEXA', 'CPENTA']:
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
            else:  # pragma: no cover
                raise NotImplementedError(name)
            nfields = nnodes + 2
            spack = Struct(endian + b'%ii' % (nfields))

        elif name == 'PLOTEL':
            key = (5201, 52, 11)
            spack = Struct(endian + b'3i')
            nfields = 3
        elif name == 'CTUBE':
            key = (3701, 37, 49)
            spack = Struct(endian + b'4i')
            nfields = 4
        elif name == 'CBAR':
            key = (2408, 24, 180)
            spack = None
            nfields = 16
        elif name == 'CBEAM':
            key = (5408, 54, 261)
            spack = None
            nfields = 18
        elif name == 'CSHEAR':
            key = (3101, 31, 61)
            spack = Struct(endian + b'6i')
            nfields = 6
        elif name == 'CQUAD4':
            key = (2958, 51, 177)
            spack = Struct(endian + b'6iffii4f')
            nfields = 14
        elif name == 'CTRIA3':
            key = (5959, 59, 282)
            spack = Struct(endian + b'5iff3i3f')
            nfields = 13
        elif name == 'CQUADR':  # same as CQUAD4
            key = (8009, 80, 367)
            spack = Struct(endian + b'6iffii4f')
            nfields = 14
        elif name == 'CTRIAR':  # same as CTRIA3
            key = (9200, 92, 385)
            spack = Struct(endian + b'5iff3i3f')
            nfields = 13
        elif name == 'CQUAD8':  # current; not 2001
            key = (4701, 47, 326)
            spack = Struct(endian + b'10i 6f i')
            nfields = 17
        elif name == 'CTRIA6':  # current; not 2001
            key = (4801, 48, 327)
            spack = Struct(endian + b'8i 5f i')
            nfields = 14
        elif name == 'CTRIAX6':
            key = (6108, 61, 107)
            spack = Struct(endian + b'8i f ii')
            nfields = 11
        elif name == 'CQUAD':
            key = (9108, 91, 507)
            spack = Struct(endian + b'11i')
        elif name == 'CQUADX':  # same as CQUAD
            key = (9008, 90, 508)
            spack = Struct(endian + b'11i')
        elif name == 'CROD':
            key = (3001, 30, 48)
            spack = Struct(endian + b'4i')
            nfields = 4
        elif name == 'CONROD':
            key = (1601, 16, 47)
            spack = Struct(endian + b'4i4f')
            nfields = 8
        elif name == 'CELAS1':
            spack = Struct(endian + b'6i')
            key = (601, 6, 73)
            nfields = 6
        elif name == 'CELAS2':
            key = (701, 7, 74)
            spack = Struct(endian + b'if4iff')
            nfields = 8
        elif name == 'CELAS3':
            key = (801, 8, 75)
            spack = Struct(endian + b'4i')
            nfields = 4
        elif name == 'CELAS4':
            key = (901, 9, 76)
            spack = Struct(endian + b'ifii')
            nfields = 4
        elif name == 'CDAMP1':
            key = (201, 2, 69)
            spack = Struct(endian + b'6i')
            nfields = 6
        elif name == 'CDAMP2':
            key = (301, 3, 70)
            spack = Struct(endian + b'if4i')
            nfields = 6
        elif name == 'CDAMP3':
            key = (401, 4, 71)
            spack = Struct(endian + b'4i')
            nfields = 4
        elif name == 'CDAMP4':
            key = (501, 5, 72)
            spack = Struct(endian + b'ifii')
            nfields = 4
        elif name == 'CDAMP5':
            key = (10608, 106, 404)
            spack = Struct(endian + b'ifii')
            nfields = 4
        elif name == 'CVISC':
            key = (3901, 39, 50)
            spack = Struct(endian + b'4i')
            nfields = 4
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
        else:  # pragma: no cover
            raise NotImplementedError(name)

        #if self.is_debug_file:
            #self.binary_debug.write('ndata=%s\n' % (nelements * 44))

        nvalues = nfields * nelements + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(pack('3i', *[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

        try:
            write_card(name, eids, spack, obj, op2, op2_ascii, endian)
        except:
            obj.log.error('failed GEOM2-%s' % name)
            raise
        op2.write(pack('i', nbytes))

        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)

    #-------------------------------------

def write_card(name, eids, spack, obj, op2, op2_ascii, endian):
    op2_ascii.write('GEOM2-%s\n' % name)
    if name in ['CTETRA', 'CHEXA', 'CPENTA']:
        if name == 'CTETRA':
            nnodes = 10
        elif name == 'CHEXA':
            nnodes = 20
        elif name == 'CPENTA':
            nnodes = 15
        else:  # pragma: no cover
            raise NotImplementedError(name)

        for eid in sorted(eids):
            elem = obj.elements[eid]
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
    elif name == 'PLOTEL':
        for eid in sorted(eids):
            elem = obj.plotels[eid]
            nids = elem.node_ids
            #(eid, n1, n2) = out
            data = [eid] + nids
            op2_ascii.write('  eid=%s nids=%s\n' % (eid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CBUSH':
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
                raise RuntimeError('invalid CBBUSH')
    elif name == 'CBUSH1D':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            #(eid, pid, g1, g2, cid, unused_a, unused_b, unused_c) = out
            g1, g2 = elem.node_ids
            cid = elem.cid
            if cid is None:
                cid = -1
            data = [eid, elem.pid, g1, g2, cid, 0, 0, 0]
            op2.write(spack.pack(*data))
    elif name == 'CGAP':
        structf = Struct(endian + b'4i3fii')
        structi = Struct(endian + b'4i3iii')
        for eid in sorted(eids):
            elem = obj.elements[eid]
            #(eid, pid, ga, gb, x1, x2, x3, f, cid) = out  # f=0,1
            pid = elem.pid
            ga, gb = elem.node_ids
            cid = elem.cid
            if cid is None:
                cid = -1

            if elem.x[0] is not None:
                f = 1
                x1, x2, x3 = elem.x
                data = [eid, pid, ga, gb, x1, x2, x3, f, cid]
                #print('CGAP x; x=%s data=%s' % (elem.x, data))
                op2.write(structf.pack(*data))
            else:
                f = 2
                g0 = elem.g0
                data = [eid, pid, ga, gb, g0, 0, 0, f, cid]
                #print('CGAP g0; x=%s data=%s' % (g0, data))
                op2.write(structi.pack(*data))

    elif name == 'CBAR':
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
    elif name == 'CBEAM':
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

    elif name in ['CQUAD4', 'CQUADR']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            #print(elem.get_stats())
            #(eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
             #t1, t2, t3, t4) = out
            theta = get_theta_from_theta_mcid(elem.theta_mcid)
            tflag = elem.tflag
            #if tflag is None:
                #tflag =
            data = [eid, pid] + nids + [theta, elem.zoffset, 0,
                                        tflag, elem.T1, elem.T2, elem.T3, elem.T4]
            assert tflag in [0, 1], elem.get_stats()
            #print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CQUAD8':  # current; not 2001
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = [nid if nid is not None else 0
                    for nid in elem.node_ids]
            pid = elem.pid
            #print(elem.get_stats())
             #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8, t1, t2,
              #t3, t4, theta, zoffs, tflag) = out # current
            #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
            #t1, t2, t3, t4, theta, zoffs) = out  # cquad8; 2001
            theta = get_theta_from_theta_mcid(elem.theta_mcid)
            tflag = elem.tflag if elem.tflag is not None else 0
            #t1 = elem.T1 if elem.T1 is not None else -1.
            #t2 = elem.T2 if elem.T2 is not None else -1.
            #t3 = elem.T3 if elem.T3 is not None else -1.
            #t4 = elem.T4 if elem.T4 is not None else -1.
            data = [eid, pid] + nids + [elem.T1, elem.T2, elem.T3, elem.T4,
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
            #print(elem.get_stats())
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
    elif name == 'CTRIAX6':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = [nid if nid is not None else 0
                    for nid in elem.node_ids]
            mid = elem.mid
            #print(elem.get_stats())
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
            #print(elem.get_stats())
            theta = get_theta_from_theta_mcid(elem.theta_mcid)
            #eid, pid, n1, n2, n3, theta_mcid, zoffs, blank1, blank2, tflag, t1, t2, t3
            data = [eid, pid] + nids + [theta, elem.zoffset, 0, 0,
                                        elem.tflag, elem.T1, elem.T2, elem.T3]
            assert elem.tflag in [0, 1], elem.get_stats()
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
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
            n1, n2 = elem.node_ids
            pid = elem.pid
            #(eid, pid, g1, g2, c1, c2)
            if n2 is None:
                n2 = 0
            data = [eid, pid, n1, n2, elem.c1, elem.c2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CELAS2':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            #(eid, k, g1, g2, c1, c2, ge, s) = out
            if n1 is None:
                n1 = 0
            if n2 is None:
                n2 = 0
            c2 = elem.c2 if elem.c2 is not None else 0
            data = [eid, elem.k, n1, n2, elem.c1, c2, elem.ge, elem.s]
            #print('CELAS2', data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name in ['CELAS3', 'CDAMP3', 'CDAMP5']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            pid = elem.pid
            if n1 is None:
                n1 = 0
            if n2 is None:
                n2 = 0
            #(eid, pid, s1, s2) = out
            data = [eid, pid, n1, n2]
            #print(name, data)
            op2_ascii.write('  eid=%s pid=%s nids=[%s, %s]\n' % (eid, pid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CELAS4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            if n1 is None:
                n1 = 0
            if n2 is None:
                n2 = 0
            #(eid, k, s1, s2) = out
            data = [eid, elem.k, n1, n2]
            #print(data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CDAMP2':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            #(eid, bdamp, g1, g2, c1, c2) = out
            if n2 is None:
                n2 = 0
            c1 = elem.c1 if elem.c1 is not None else 0
            c2 = elem.c2 if elem.c2 is not None else 0
            data = [eid, elem.b, n1, n2, c1, c2]
            #print(name, data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'CDAMP4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            if n1 is None:
                n1 = 0
            if n2 is None:
                n2 = 0
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
