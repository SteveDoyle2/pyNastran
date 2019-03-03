from __future__ import print_function, absolute_import
from struct import pack, Struct

from .geom1 import write_geom_header, close_geom_table

def write_geom2(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'elements'):
        return
    #if not hasattr(obj, 'nodes'):
        #return
    nspoints = len(obj.spoints)
    nelements = len(obj.elements)
    if nelements == 0 and nspoints == 0:
        return
    write_geom_header(b'GEOM2', op2, op2_ascii)
    itable = -3

    etypes = [
        'CROD', 'CONROD',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        'CTRIA3', 'CQUAD4',
        'CTETRA', 'CHEXA', 'CPENTA',
    ]
    out = obj.get_card_ids_by_card_types(etypes)
    if nspoints:
        print(obj.spoints)
        out['SPOINT'] = obj.spoints
        aaa
    for name, eids in out.items():
        nelements = len(eids)
        if nelements == 0:
            continue

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

        elif name == 'CQUAD4':
            key = (2958, 51, 177)
            spack = Struct(endian + b'6iffii4f')
            nfields = 14
            nnodes = 4
        elif name == 'CTRIA3':
            key = (5959, 59, 282)
            spack = Struct(endian + b'5iff3i3f')
            nfields = 13
            nnodes = 3
        elif name == 'CROD':
            spack = Struct(endian + b'4i')
            key = (3001, 30, 48)
            nfields = 4
        elif name == 'CONROD':
            spack = Struct(endian + b'4i4f')
            key = (1601, 16, 47)
            nfields = 8
        elif name == 'CELAS1':
            spack = Struct(endian + b'6i')
            key = (601, 6, 73)
            nfields = 6
        elif name == 'CELAS2':
            spack = Struct(endian + b'if4iff')
            key = (701, 7, 74)
            nfields = 8
        elif name == 'CELAS3':
            spack = Struct(endian + b'4i')
            key = (801, 8, 75)
            nfields = 4
        elif name == 'CELAS4':
            spack = Struct(endian + b'ifii')
            key = (901, 9, 76)
            nfields = 4

        elif name == 'CDAMP1':
            spack = Struct(endian + b'6i')
            key = (201, 2, 69)
            nfields = 6
        elif name == 'CDAMP2':
            spack = Struct(endian + b'if4i')
            key = (301, 3, 70)
            nfields = 6
        elif name == 'CDAMP3':
            spack = Struct(endian + b'4i')
            key = (401, 4, 71)
            nfields = 4
        elif name == 'CDAMP4':
            spack = Struct(endian + b'ifii')
            key = (501, 5, 72)
            nfields = 4
        elif name == 'SPOINT':
            spack = None
            key = (5551, 49, 105)
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
            write_card(name, eids, spack, obj, op2, op2_ascii)
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
    close_geom_table(op2, op2_ascii, itable, include_last=False)

    #-------------------------------------

def write_card(name, eids, spack, obj, op2, op2_ascii):
    if name in ['CTETRA', 'CHEXA', 'CPENTA']:
        if name == 'CTETRA':
            nnodes = 10
        elif name == 'CHEXA':
            nnodes = 20
        elif name == 'CPENTA':
            nnodes = 15

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
    elif name == 'CQUAD4':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            #print(elem.get_stats())
            #(eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
             #t1, t2, t3, t4) = out
            theta = 0.0
            data = [eid, pid] + nids + [theta, elem.zoffset, 0,
                                        elem.tflag, elem.T1, elem.T2, elem.T3, elem.T4]
            assert elem.tflag in [0, 1], elem.get_stats()
            #print('  CQUAD4 eid=%s pid=%s nids=%s data=%s\n' % (eid, pid, str(nids), data[6:]))
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CTRIA3':
        for eid in sorted(eids):
            elem = obj.elements[eid]
            nids = elem.node_ids
            pid = elem.pid
            #print(elem.get_stats())
            theta = 0.0
            #eid, pid, n1, n2, n3, theta_mcid, zoffs, blank1, blank2, tflag, t1, t2, t3
            data = [eid, pid] + nids + [theta, elem.zoffset, 0, 0,
                                        elem.tflag, elem.T1, elem.T2, elem.T3]
            assert elem.tflag in [0, 1], elem.get_stats()
            op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
            op2.write(spack.pack(*data))
    elif name == 'CROD':
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
            print(data)
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
            print(data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name in ['CELAS3', 'CDAMP3']:
        for eid in sorted(eids):
            elem = obj.elements[eid]
            n1, n2 = elem.node_ids
            pid = elem.pid
            if n1 is None:
                n1 = 0
            #(eid, pid, s1, s2) = out
            data = [eid, pid, n1, n2]
            print(data)
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
            print(data)
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
            print(data)
            op2_ascii.write('  eid=%s nids=[%s, %s]\n' % (eid, n1, n2))
            op2.write(spack.pack(*data))
    elif name == 'SPOINT':
        nids = eids
        nids.sort()
        spack = Struct('%ii' % len(nids))
        spack.pack(*nids)
    else:  # pragma: no cover
        raise NotImplementedError(name)
