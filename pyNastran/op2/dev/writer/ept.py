from __future__ import absolute_import
from struct import pack, Struct

from .geom1 import init_table, close_geom_table

def write_ept(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'properties'):
        return
    #if not hasattr(obj, 'nodes'):
        #return
    nproperties = len(obj.properties)
    if nproperties == 0:
        return
    data = init_table(b'EPT')
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

    ptypes = [
        'PSOLID', 'PSHELL', 'PCOMP', 'PROD',
    ]
    out = obj.get_card_ids_by_card_types(ptypes)
    for name, pids in out.items():
        nproperties = len(pids)
        if nproperties == 0:
            continue

        if name == 'PSOLID':
            key = (2402, 24, 281)
            nfields = 7
            spack = Struct(endian + b'6i4s')
        else:  # pragma: no cover
            raise NotImplementedError(name)

        nvalues = nfields * nproperties + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(pack('3i', *[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

        if name == 'PSOLID':
            #pid = data[0]
            #mid = data[1]
            #cordm = data[2]
            #integ = data[3]
            #stress = data[4]
            #isop = data[5]
            #fctn = data[6].decode('latin1')
            for pid in sorted(pids):
                prop = obj.properties[pid]
                mid = prop.mid
                cordm = prop.cordm
                #print(prop.get_stats())

                ## stress : int, string, or blank
                ##    blank/GRID
                ##    1-GAUSS
                if prop.stress == 'GRID':
                    stress = 0
                elif prop.stress == 'GAUSS':
                    stress = 1
                else:
                    raise NotImplementedError('prop.stress=%s and must be [0, 1]' % prop.stress)

                if prop.integ == 'BUBBLE':
                    integ = 0
                elif prop.integ == 'GAUSS':
                    integ = 1
                elif prop.integ == 'TWO':
                    integ = 2
                elif prop.integ == 'THREE':
                    integ = 3
                else:
                    raise NotImplementedError('prop.stress=%s and must be [0, 1, 2, 3]' % prop.integ)

                if prop.isop == 'REDUCED':
                    isop = 0
                elif prop.isop == 'FULL':
                    isop = 1
                else:
                    raise NotImplementedError('isop=%s and must be [0, 1]' % prop.isop)

                if prop.fctn == 'SMECH':
                    fctn = b'SMEC'
                elif prop.fctn == 'PFLUID':
                    fctn = b'PFLU'
                else:
                    raise NotImplementedError('PSOLID; fctn=%r' % prop.fctn)


                data = [pid, mid, cordm, integ, stress, isop, fctn]
                print(data)
                op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, mid, data[2:]))
                op2.write(spack.pack(*data))
        else:  # pragma: no cover
            raise NotImplementedError(name)
        op2.write(pack('i', nbytes))
        itable -= 1

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)
    #-------------------------------------

