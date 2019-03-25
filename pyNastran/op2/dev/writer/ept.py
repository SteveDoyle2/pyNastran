from __future__ import absolute_import, print_function
from collections import defaultdict
from struct import pack, Struct

from .geom1 import write_geom_header, close_geom_table

def write_ept(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'properties'):
        return

    out = defaultdict(list)
    for pid, phbdy in obj.phbdys.items():
        out[phbdy.type].append(pid)
    for pid, pelast in obj.pelast.items():
        out[pelast.type].append(pid)

    #if not hasattr(obj, 'nodes'):
        #return
    nproperties = len(obj.properties)
    if nproperties == 0:
        return
    write_geom_header(b'EPT', op2, op2_ascii, endian=endian)
    struct_3i = Struct(endian + b'3i')

    itable = -3

    #ptypes = [
        #'PSOLID', 'PSHELL', 'PCOMP', 'PROD',
    #]
    #out = obj.get_card_ids_by_card_types(ptypes)
    for pid, prop in obj.properties.items():
        out[prop.type].append(pid)

    skip_properties = [
        'PBEND', 'PBUSH', 'PBUSH1D',
        'PMASS', 'PBCOMP',

        # thermal
        'PHBDY',

    ]
    for name, pids in out.items():
        nproperties = len(pids)
        if nproperties == 0:
            continue
        elif name in skip_properties:
            obj.log.warning('skipping EPT-%s' % name)
            continue

        if name == 'PBARL':
            key = (9102, 91, 52)
            continue
        elif name == 'PCOMP':
            itable = write_pcomp(name, pids, itable, op2, op2_ascii, obj, endian=endian)
            continue
        elif name == 'PELAS':
            key = (302, 3, 46)
            nfields = 4
            spack = Struct(endian + b'i3f')
        elif name == 'PDAMP':
            key = (202, 2, 45)
            nfields = 2
            spack = Struct(endian + b'if')
        elif name == 'PVISC':
            key = (1802, 18, 31)
            nfields = 3
            spack = Struct(endian + b'i2f')
        elif name == 'PROD':
            key = (902, 9, 29)
            nfields = 6
            spack = Struct(endian + b'2i4f')
        elif name == 'PTUBE':
            key = (1602, 16, 30)
            nfields = 5
            spack = Struct(endian + b'2i3f')
        elif name == 'PGAP':
            key = (3201, 32, 55)
            nfields = 11
            spack = Struct(endian + b'i10f')
        elif name == 'PBAR':
            key = (52, 20, 181)
            nfields = 19
            spack = Struct(endian + b'2i17f')
        elif name == 'PBEAM':
            key = (5402, 54, 262)
            nfields = 197 # 5+16*12
            spack = None
        elif name == 'PSHEAR':
            key = (1002, 10, 42)
            nfields = 6
            spack = Struct(endian + b'2i4f')
        elif name == 'PSHELL':
            key = (2302, 23, 283)
            nfields = 11
            spack = Struct(endian + b'iififi4fi')
        elif name == 'PLPLANE':
            key = (4606, 46, 375)
            nfields = 11
            spack = Struct(endian + b'3i 4s f 6i')
        elif name == 'PSOLID':
            key = (2402, 24, 281)
            nfields = 7
            spack = Struct(endian + b'6i4s')
        elif name == 'PLSOLID':
            key = (4706, 47, 376)
            nfields = 7
            spack = Struct(endian + b'2i 4s 4i')
        elif name == 'PMASS':
            key = (402, 4, 44)
            nfields = 2
            spack = Struct(endian + b'if')
        elif name == 'PELAST':
            key = (1302, 13, 34)
            nfields = 4
            spack = Struct(endian + b'4i')
        else:  # pragma: no cover
            raise NotImplementedError(name)

        nvalues = nfields * nproperties + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2.write(struct_3i.pack(*[4, nvalues, 4]))
        op2.write(pack('i', nbytes)) #values, nbtyes))

        op2.write(struct_3i.pack(*key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

        try:
            write_card(op2, op2_ascii, obj, name, pids, spack, endian)
        except:
            obj.log.error('failed EPT-%s' % name)
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

def write_card(op2, op2_ascii, obj, name, pids, spack, endian):
    op2_ascii.write('EPT-%s\n' % name)
    if name == 'PELAS':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, k, ge, s) = out
            data = [pid, prop.k, prop.ge, prop.s]
            op2.write(spack.pack(*data))
    elif name == 'PDAMP':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, b) = out
            data = [pid, prop.b]
            op2.write(spack.pack(*data))
    elif name == 'PVISC':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, ce, cr) = out
            data = [pid, prop.ce, prop.cr]
            op2.write(spack.pack(*data))
    elif name == 'PTUBE':
        #.. todo:: OD2 only exists for heat transfer...
        #          how do i know if there's heat transfer?
        #          I could store all the tubes and add them later,
        #          but what about themal/non-thermal subcases?
        #.. warning:: assuming OD2 is not written (only done for thermal)
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, OD, t, nsm) = out
            data = [pid, prop.mid, prop.OD1, prop.t, prop.nsm]
            op2.write(spack.pack(*data))

    elif name == 'PGAP':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            data = [pid, prop.u0, prop.f0, prop.ka, prop.kb, prop.kt,
                    prop.mu1, prop.mu2, prop.tmax, prop.mar, prop.trmin]
            assert None not in data, data
            op2_ascii.write('  pid=%s data=%s' % (pid, data[1:]))
            op2.write(spack.pack(*data))

    elif name == 'PBAR':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, a, I1, I2, J, nsm, fe, c1, c2, d1, d2,
             #e1, e2, f1, f2, k1, k2, I12) = out
            fe = 0
            k1 = prop.k1 if prop.k1 is not None else 0
            k2 = prop.k2 if prop.k2 is not None else 0
            data = [
                pid, prop.mid, prop.A, prop.i1, prop.i2, prop.j, prop.nsm,
                fe, prop.c1, prop.c2, prop.d1, prop.d2, prop.e1, prop.e2,
                prop.f1, prop.f2, k1, k2, prop.i12]
            assert None not in data, data
            op2.write(spack.pack(*data))
    elif name == 'PBEAM':  # probably wrong stations
        struct1 = Struct(endian + b'4if')
        struct2 = Struct(endian + b'16f')
        struct3 = Struct(endian + b'16f')
        for pid in sorted(pids):
            nfieldsi = 0
            prop = obj.properties[pid]
            nsegments = len(prop.xxb)
            if nsegments == 2:
                ccf = 1 # True
            elif nsegments > 2:
                ccf = 0 # False
            else:
                raise NotImplementedError(nsegments)

            #(pid, mid, nsegments, ccf, x) = data_in
            x = 0.
            data = [pid, prop.mid, nsegments, ccf, x]
            nfieldsi += len(data)
            op2.write(struct1.pack(*data))

            nzero_segments = nsegments - 1
            j = 0
            for i in range(11):
                if i > nzero_segments:
                    data = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                            0., 0., 0., 0., 0., 0.]
                else:
                    #(soi, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                     #d1, d2, e1, e2, f1, f2) = pack
                    xxb = prop.xxb[j]
                    so_str = prop.so[j]
                    a = prop.A[j]
                    i1 = prop.i1[j]
                    i2 = prop.i2[j]
                    i12 = prop.i12[j]
                    nsm = prop.nsm[j]
                    c1 = prop.c1[j]
                    c2 = prop.c2[j]
                    d1 = prop.d1[j]
                    d2 = prop.d2[j]
                    e1 = prop.e1[j]
                    e2 = prop.e2[j]
                    f1 = prop.f1[j]
                    f2 = prop.f2[j]

                    if so_str == 'NO':
                        soi = 0.0
                    elif so_str == 'YES':
                        soi = 1.0
                    else:
                        try:
                            soi = float(so_str)
                        except ValueError:  # pragma: no cover
                            print(prop.get_stats())
                            raise NotImplementedError('SO=%s SO%s=%r' % (prop.so, j, so_str))
                    j += 1
                    data = [
                        soi, xxb, a, i1, i2, i12, j, nsm, c1, c2,
                        d1, d2, e1, e2, f1, f2]
                    #so_str = str(soi)
                assert None not in data, data
                nfieldsi += len(data)
                op2.write(struct2.pack(*data))

            #self.log.info('PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s' % tuple(data_in))

            # Constant cross-section flag: 1=yes and 0=no
            # what is 2?
            #if ccf not in [0, 1, 2]:
                #msg = ('  PBEAM pid=%s mid=%s nsegments=%s ccf=%s x=%s; '
                       #'ccf must be in [0, 1, 2]\n' % tuple(data_in))
                #raise ValueError(msg)

            #(k1, k2, s1, s2, nsia, nsib, cwa, cwb, # 8
             #m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = endpack # 8 -> 16
            data = [
                prop.k1, prop.k2, prop.s1, prop.s2,
                prop.nsia, prop.nsib,
                prop.cwa, prop.cwb,
                prop.m1a, prop.m2a, prop.m1b, prop.m2b,
                prop.n1a, prop.n2a, prop.n1b, prop.n2b,
            ]
            #k1 / k2 : float; default=1.
            #s1 / s2 : float; default=0.
            #nsia / nsia : float; default=0. / nsia
            #cwa / cwb : float; default=0. / cwa
            #m1a / m2a : float; default=0. / m1a
            #m1b / m2b : float; default=0. / m1b
            #n1a / n2a : float; default=0. / n1a
            #n1b / n2b : float; default=0. / n1b
            assert None not in data, data
            nfieldsi += len(data)
            op2.write(struct3.pack(*data))
            assert nfieldsi == 197, nfieldsi

    elif name == 'PSOLID':
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
    elif name == 'PSHEAR':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, t, nsm, f1, f2) = out
            data = [pid, prop.mid, prop.t, prop.nsm, prop.f1, prop.f2]
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2.write(spack.pack(*data))
    elif name == 'PSHELL':
        for pid in sorted(pids):
            #(pid, mid1, t, mid2, bk, mid3, ts, nsm, z1, z2, mid4) = out
            prop = obj.properties[pid]
            mid2 = prop.mid2 if prop.mid2 is not None else 0
            data = [pid, prop.mid1, prop.t, mid2, prop.twelveIt3, prop.mid3,
                    prop.tst, prop.nsm, prop.z1, prop.z2, prop.mid4]
            #print('PSHELL', data)
            #print(prop.mid1, mid2, prop.mid3, prop.mid4)

            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid1, data[2:]))
            op2.write(spack.pack(*data))
    elif name == 'PLPLANE':
        #NX 10
        #1 PID     I Property identification number
        #2 MID     I Material identification number
        #3 CID     I Coordinate system identification number
        #4 STR CHAR4 Location of stress and strain output
        #5 T      RS Default membrane thickness for Ti on the connection entry
        #6 CSOPT  I  Reserved for coordinate system definition of plane
        #7 UNDEF(5) None

        #MSC 2016
        #PID       I Property identification number
        #2 MID     I Material identification number
        #3 CID     I Coordinate system identification number
        #4 STR CHAR4 Location of stress and strain output
        #5 UNDEF(7 ) none Not used

        #.. warning:: CSOPT ad T are not supported
        for pid in sorted(pids):
            prop = obj.properties[pid]
            location = prop.stress_strain_output_location.encode('latin1')

            # MSC
            data = [pid, prop.mid, prop.cid, location,
                    0, 0, #prop.t, prop.csopt, # unsupported NX
                    0, 0, 0, 0, 0]
            #print(name, data)
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2.write(spack.pack(*data))
    elif name == 'PROD':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, a, j, c, nsm) = out
            data = [pid, prop.mid, prop.A, prop.j, prop.c, prop.nsm]
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2.write(spack.pack(*data))
    elif name == 'PLSOLID':
        #MSC 2016
        #1 PID I Property identification number
        #2 MID I Material identification number
        #3 STR CHAR4 Location of stress and strain output
        #4 UNDEF(4 ) none Not used

        #NX 10
        #1 PID I Property identification number
        #2 MID I Material identification number
        #3 STR CHAR4 Location of stress and strain output
        #4 CSOPT I Reserved for coordinate system definition of plane
        #5 UNDEF(3) None

        #.. warning:: CSOPT is not supported

        for pid in sorted(pids):
            prop = obj.properties[pid]
            location = prop.stress_strain.encode('latin1')

            # MSC
            #pid, mid, location, csopt, null_a, null_b, null_c = out
            csopt = 0
            data = [pid, prop.mid, location, csopt, 0, 0, 0]
            #print(name, data)
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2.write(spack.pack(*data))
    elif name == 'PMASS':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            data = [pid, prop.mass]
            op2_ascii.write('  pid=%s mass=%s\n' % (pid, prop.mass))
            op2.write(spack.pack(*data))
    elif name == 'PELAST':
        for pid in sorted(pids):
            prop = obj.pelast[pid]
            data = [pid, prop.tkid, prop.tgeid, prop.tknid]
            op2_ascii.write('  pid=%s tables=%s\n' % (pid, data[1:]))
            op2.write(spack.pack(*data))
    else:  # pragma: no cover
        raise NotImplementedError(name)


def write_pbarl(name, pids, itable, op2, op2_ascii, obj, endian=b'<'):
    nproperties = len(pids)
    for pid in sorted(pids):
        prop = obj.properties[pid]
    return itable

def write_pcomp(name, pids, itable, op2, op2_ascii, obj, endian=b'<'):
    """writes the PCOMP"""
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
                               '\nPCOMP = %s' % (pid, prop.ft, prop))

        #is_symmetric = True
        symmetric_factor = 1
        data = [pid, symmetric_factor * prop.nplies, prop.z0,
                prop.nsm, prop.sb, ft, prop.tref, prop.ge]
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
