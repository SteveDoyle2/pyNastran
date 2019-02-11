from __future__ import absolute_import
from struct import pack, Struct

from .geom1 import close_geom_table

def write_geom2(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'elements'):
        return
    #if not hasattr(obj, 'nodes'):
        #return
    nelements = len(obj.elements)
    if nelements == 0:
        return
    data = [
        4, 2, 4,
        #4, 2,4,
        8, b'GEOM2   ', 8,
        4, -1, 4,
        #4, 1, 4,
        #4, 0, 4,
    ]
    op2.write(pack('4i 8s i 3i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        4, 7, 4,
        28, 1, 2, 3, 4, 5, 6, 7, 28,
    ]
    op2.write(pack('3i 9i', *data))
    op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    data = [
        4, -2, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        #4, 0, 4,
        4, 2, 4,
        8, 1, 2, 8,
    ]
    op2.write(pack('3i 4i', *data))
    op2_ascii.write(str(data) + '\n')
    #data = [8, 1, 2, 8]
    #op2.write(pack('4i', *data))
    #-------------------------------------

    data = [
        4, -3, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
    itable = -3

    etypes = [
        'CROD', 'CONROD',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        'CTRIA3', 'CQUAD4',
        'CTETRA', 'CHEXA', 'CPENTA',
    ]
    out = obj.get_card_ids_by_card_types(etypes)
    for name, eids in out.items():
        nelements = len(eids)
        if nelements == 0:
            continue

        if name == 'CTETRA':
            key = (5508, 55, 217)
            nnodes = 10
            # 12 = eid, pid, n1, n2, n3, n4, ..., n10
        elif name == 'CHEXA':
            key = (7308, 73, 253)
            nnodes = 20
        else:  # pragma: no cover
            raise NotImplementedError(name)

        # add eid, pid
        nfields = nnodes + 2

        spack = Struct(endian + b'%ii' % (nfields))

        nvalues = nfields * nelements + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(pack('3i', *[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

        if name in ['CTETRA', 'CHEXA', 'CPENTA']:
            for eid in sorted(eids):
                elem = obj.elements[eid]
                nids = elem.node_ids
                pid = elem.pid
                assert None not in nids, nids
                nnids = len(nids)
                if nnids < nnodes:
                    nids2 = [0] * (nnodes - nnids)
                    data = [eid, pid] + nids + nids2
                else:
                    data = [eid, pid] + nids
                #print(name, data)
                op2_ascii.write('  eid=%s pid=%s nids=%s\n' % (eid, pid, str(nids)))
                op2.write(spack.pack(*data))
        else:  # pragma: no cover
            raise NotImplementedError(name)
        op2.write(pack('i', nbytes))
        itable -= 1

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)

    #-------------------------------------

