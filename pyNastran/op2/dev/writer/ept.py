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
    struct_3i = Struct(endian + b'3i')
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
        elif name == 'PSHELL':
            key = (2302, 23, 283)
            nfields = 11
            spack = Struct(endian + b'iififi4fi')
        elif name == 'PCOMP':
            itable = write_pcomp(name, pids, itable, op2, op2_ascii, obj, endian=endian)
            continue
        elif name == 'PROD':
            key = (902, 9, 29)
            nfields = 6
            spack = Struct(endian + b'2i4f')
        else:  # pragma: no cover
            raise NotImplementedError(name)

        nvalues = nfields * nproperties + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(struct_3i.pack(*[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(struct_3i.pack(*key))
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
                op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, mid, data[2:]))
                op2.write(spack.pack(*data))
        elif name == 'PSHELL':
            for pid in sorted(pids):
                #(pid, mid1, t, mid2, bk, mid3, ts, nsm, z1, z2, mid4) = out
                prop = obj.properties[pid]
                data = [pid, prop.mid1, prop.t, prop.mid2, prop.twelveIt3, prop.mid3,
                        prop.tst, prop.nsm, prop.z1, prop.z2, prop.mid4]
                op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid1, data[2:]))
                op2.write(spack.pack(*data))
        elif name == 'PROD':
            for pid in sorted(pids):
                #(pid, mid, a, j, c, nsm) = out
                prop = obj.properties[pid]
                #print(prop.get_stats())
                data = [pid, prop.mid, prop.A, prop.j, prop.c, prop.nsm]
                op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
                op2.write(spack.pack(*data))
        else:  # pragma: no cover
            raise NotImplementedError(name)
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

def write_pcomp(name, pids, itable, op2, op2_ascii, obj, endian=b'<'):
    key = (2706, 27, 287)

    nproperties = len(pids)
    nlayers = 0
    for pid in sorted(pids):
        #(pid, mid, a, j, c, nsm) = out
        prop = obj.properties[pid]
        #(pid, nlayers, z0, nsm, sb, ft, Tref, ge) = out # 8
        nlayers += prop.nplies

    nvalues = 8 * nproperties + (4  * nlayers) + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))

    #is_symmetrical = 'NO'
    #if nlayers < 0:
        #is_symmetrical = 'SYM'
        #nlayers = abs(nlayers)
    #assert nlayers > 0, out

    s1 = Struct(endian + b'2i3fi2f')
    s2 = Struct(endian + b'i2fi')
    for pid in sorted(pids):
        prop = obj.properties[pid]
        #print(prop.get_stats())

        if prop.ft is None:
            ft = 0
        elif prop.ft == 'HILL':
            ft = 1
        elif prop.ft == 'HOFF':
            ft = 2
        elif prop.ft == 'TSAI':
            ft = 3
        elif prop.ft == 'STRN':
            ft = 4
        else:
            raise RuntimeError('unsupported ft.  pid=%s ft=%r.'
                               '\nPCOMP = %s' % (pid, prop.ft, data))

        #is_symmetric = True
        symmetric_factor = 1
        data = [pid, symmetric_factor * prop.nplies, prop.z0, prop.nsm, prop.sb, ft, prop.tref, prop.ge]
        op2.write(s1.pack(*data))
        op2_ascii.write(str(data) + '\n')

        for (mid, t, theta, sout) in zip(prop.mids, prop.thicknesses, prop.thetas, prop.souts):
            if sout == 'NO':
                sout = 0
            elif sout == 'YES':
                sout = 1
            else:
                raise RuntimeError('unsupported sout.  sout=%r and must be 0 or 1.'
                                   '\nPCOMP = %s' % (sout, data))
            data2 = [mid, t, theta, sout]
            op2.write(s2.pack(*data2))
            op2_ascii.write(str(data2) + '\n')


    #for ilayer in range(nlayers):
        #(mid, t, theta, sout) = s2.unpack(data[n:n+16])
        #mids.append(mid)
        #T.append(t)
        #thetas.append(theta)
        #souts.append(sout)
        #n += 16

    #(mid, t, theta, sout) = s2.unpack(data[n:n+16])

    #data_in = [
        #pid, z0, nsm, sb, ft, Tref, ge,
        #is_symmetrical, mids, T, thetas, souts]

    op2.write(pack('i', nbytes))
    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

    return itable
