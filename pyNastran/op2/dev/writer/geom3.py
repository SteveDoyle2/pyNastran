from __future__ import absolute_import
from struct import pack, Struct
from collections import defaultdict

from .geom1 import write_geom_header, close_geom_table

def write_geom3(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'loads') and not hasattr(obj, 'load_combinations'):
        return
    loads_by_type = defaultdict(list)
    for load_id, loads in obj.loads.items():
        for load in loads:
            loads_by_type[load.type].append(load)

    # return if no supported cards are found
    supported_cards = ['FORCE', 'GRAV']
    for load_type in sorted(loads_by_type.keys()):
        if load_type in supported_cards:
            break
        obj.log.warning('skipping GEOM3-%s' % load_type)
    else:
        return

    write_geom_header(b'GEOM3', op2, op2_ascii)
    itable = -3
    for load_type, loads in sorted(loads_by_type.items()):
        nloads = len(loads)
        if load_type == 'FORCE':
            key = (4201, 42, 18)
            nfields = 7
            spack = Struct(b'3i 4f')
            nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
            for load in loads:
                data = [load.sid, load.node_id, load.Cid(), load.mag] + list(load.xyz)
                op2_ascii.write('  FORCE data=%s\n' % str(data))
                op2.write(spack.pack(*data))
        elif load_type == 'GRAV':
            key = (4401, 44, 26)
            nfields = 7
            #(sid, cid, a, n1, n2, n3, mb) = out
            spack = Struct(b'ii4fi')
            nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
            for load in loads:
                data = [load.sid, load.Cid(), load.scale] + list(load.N) + [load.mb]
                op2_ascii.write('  GRAV data=%s\n' % str(data))
                op2.write(spack.pack(*data))
        else:
            load0 = loads[0]
            raise NotImplementedError(load0)

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


def write_header(name, nfields, nloads, key, op2, op2_ascii):
    nvalues = nfields * nloads + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes
