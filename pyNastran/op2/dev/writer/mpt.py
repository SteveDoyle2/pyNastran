from __future__ import absolute_import
from struct import pack, Struct

from .geom1 import init_table, close_geom_table

def write_mpt(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'materials'):
        return
    #if not hasattr(obj, 'nodes'):
        #return
    nmaterials = len(obj.materials)
    if nmaterials == 0:
        return
    data = init_table(b'MPT')
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

    mtypes = [
        'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9'
    ]
    out = obj.get_card_ids_by_card_types(mtypes)
    for name, mids in out.items():
        nmaterials = len(mids)
        if nmaterials == 0:
            continue

        if name == 'MAT1':
            key = (103, 1, 77)
            nfields = 12
            spack = Struct(endian + b'i10fi')
            # mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
        else:  # pragma: no cover
            raise NotImplementedError(name)

        nvalues = nfields * nmaterials + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(pack('3i', *[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

        if name == 'MAT1':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                #print(mat.get_stats())
                #mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
                data = [mid, mat.e, mat.g, mat.nu, mat.rho, mat.a, mat.tref,
                        mat.ge, mat.St, mat.Sc, mat.Ss, mat.mcsid]

                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        else:  # pragma: no cover
            raise NotImplementedError(name)
        op2.write(pack('i', nbytes))
        itable -= 1

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)
    #-------------------------------------

