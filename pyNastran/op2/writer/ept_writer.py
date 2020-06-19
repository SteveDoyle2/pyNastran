from collections import defaultdict
from struct import pack, Struct

from .geom1_writer import write_geom_header, close_geom_table

def write_ept(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'properties'):
        return

    out = defaultdict(list)
    for pid, phbdy in obj.phbdys.items():
        out[phbdy.type].append(pid)
    for pid, pelast in obj.pelast.items():
        out[pelast.type].append(pid)
    for pid, pdampt in obj.pdampt.items():
        out[pdampt.type].append(pid)
    for pid, pbusht in obj.pbusht.items():
        out[pbusht.type].append(pid)

    #if not hasattr(obj, 'nodes'):
        #return
    nproperties = len(obj.properties) + len(obj.properties_mass) + len(out)
    if nproperties == 0:
        return
    write_geom_header(b'EPT', op2, op2_ascii, endian=endian)
    struct_3i = Struct(endian + b'3i')

    itable = -3

    #ptypes = [
        #'PSOLID', 'PSHELL', 'PCOMP', 'PROD',

        # thermal
        #'PHBDY',
    #]
    #out = obj.get_card_ids_by_card_types(ptypes)
    for pid, prop in obj.properties.items():
        out[prop.type].append(pid)
    for pid, prop in obj.properties_mass.items():
        out[prop.type].append(pid)

    skip_properties = [
        'PBEND', #'PBUSH1D',
        'PBCOMP',
        #'PGAP',
        #'PCOMPG',

    ]
    nastran_format = 'nx'
    for name, pids in out.items():
        nproperties = len(pids)
        if nproperties == 0:  # pragma: no cover
            continue
        if name in skip_properties:  # pragma: no cover
            obj.log.warning('skipping EPT-%s' % name)
            continue
        #obj.log.debug('writing EPT-%s' % name)

        #print('EPT', itable, name)
        if name == 'PBARL':
            itable = write_pbarl(name, pids, itable, op2, op2_ascii, obj, endian=endian)
            continue
        elif name == 'PCOMP':
            itable = write_pcomp(name, pids, itable, op2, op2_ascii, obj, endian=endian)
            continue
        elif name == 'PCOMPG':
            itable = write_pcompg(name, pids, itable, op2, op2_ascii, obj, endian=endian)
            continue
        elif name == 'PBUSH':
            itable = write_pbush(name, pids, itable, op2, op2_ascii, obj, endian=endian,
                                 nastran_format=nastran_format)
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
            key = (2102, 21, 121)
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
        elif name == 'PHBDY':
            key = (2802, 28, 236)
            nfields = 4
            spack = Struct(endian + b'i3f')
        elif name == 'PBUSH1D':
            key = (3101, 31, 219)
            nfields = 38
            spack = Struct(endian + b'i 6f i 4f 24i 2f')
        elif name == 'PIHEX':
            obj.log.warning('skipping PIHEX')
            continue
        else:
            obj.log.warning(f'skipping {name}')
            continue
        #else:  # pragma: no cover
            #raise NotImplementedError(name)

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
    elif name == 'PBUSH1D':
        _write_pbush1d(name, obj, pids, spack, op2, op2_ascii, endian)
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
        _write_pbeam(name, obj, pids, op2, op2_ascii, endian)

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
            if prop.stress == 'GRID' or prop.stress is None:
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
            elif  prop.integ is None:  # TODO: not sure
                integ = 0
            else:
                raise NotImplementedError('prop.integ=%s and must be [0, 1, 2, 3]' % prop.integ)

            if prop.isop == 'REDUCED':
                isop = 0
            elif prop.isop == 'FULL':
                isop = 1
            elif  prop.isop is None:  # TODO: not sure
                isop = 0
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
            mid1 = 0 if prop.mid1 is None else prop.mid1
            mid2 = 0 if prop.mid2 is None else prop.mid2
            mid3 = 0 if prop.mid3 is None else prop.mid3
            mid4 = 0 if prop.mid4 is None else prop.mid4
            data = [pid, mid1, prop.t, mid2, prop.twelveIt3, mid3,
                    prop.tst, prop.nsm, prop.z1, prop.z2, mid4]
            #print('PSHELL', data)
            #print(prop.mid1, mid2, prop.mid3, prop.mid4)

            op2_ascii.write(f'  {name} pid={pid} mid1={mid1} data={data[2:]}\n')
            assert None not in data, f'  {name} pid={pid} mid1={mid1} data={data[2:]}'
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
            prop = obj.properties_mass[pid]
            data = [pid, prop.mass]
            op2_ascii.write('  pid=%s mass=%s\n' % (pid, prop.mass))
            op2.write(spack.pack(*data))
    elif name == 'PELAST':
        for pid in sorted(pids):
            prop = obj.pelast[pid]
            data = [pid, prop.tkid, prop.tgeid, prop.tknid]
            op2_ascii.write('  pid=%s tables=%s\n' % (pid, data[1:]))
            op2.write(spack.pack(*data))
    elif name == 'PHBDY':
        for pid in sorted(pids):
            prop = obj.phbdys[pid]
            af = 0.0 if prop.af is None else prop.af
            d1 = 0.0 if prop.d1 is None else prop.d1
            d2 = 0.0 if prop.d2 is None else prop.d2
            data = [pid, af, d1, d2]
            op2_ascii.write('  pid=%s [af,d1,d2]=%s\n' % (pid, data[1:]))
            #print('  pid=%s [af,d1,d2]=%s\n' % (pid, data[1:]))
            op2.write(spack.pack(*data))
            #(pid, af, d1, d2) = out
    else:  # pragma: no cover
        raise NotImplementedError(name)


def write_pbush(name, pids, itable, op2, op2_ascii, obj, endian=b'<',
                nastran_format='nx'):
    """writes the PBUSH"""

    # TODO: there are 3 different types of PBUSH cards...
    key = (1402, 14, 37) # 23 fields
    fmt = endian + b'i22f'
    struct1 = Struct(fmt)

    nproperties = len(pids)

    nvalues = 23 * nproperties + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))

    for pid in sorted(pids):
        prop = obj.properties[pid]
        (b1, b2, b3, b4, b5, b6) = (0., 0., 0., 0., 0., 0.)
        (k1, k2, k3, k4, k5, k6) = (0., 0., 0., 0., 0., 0.)
        (g1, g2, g3, g4, g5, g6) = (0., 0., 0., 0., 0., 0.)
        if prop.Bi:
            (b1, b2, b3, b4, b5, b6) = [
                bi if bi is not None else 0.0 for bi in prop.Bi]  # damping
        if prop.Ki:
            (k1, k2, k3, k4, k5, k6) = [
                ki if ki is not None else 0.0 for ki in prop.Ki] # stiffness
        if prop.GEi:
            (g1, g2, g3, g4, g5, g6) = [
                gi if gi is not None else 0.0 for gi in prop.GEi] # ???
        # TODO: not 100%
        sa = prop.sa if prop.sa is not None else 0.
        st = prop.st if prop.st is not None else 0.
        ea = prop.ea if prop.ea is not None else 0.
        et = prop.et if prop.et is not None else 0.
        data_in = [pid,
                   k1, k2, k3, k4, k5, k6,
                   b1, b2, b3, b4, b5, b6,
                   g1, g2, g3, g4, g5, g6,
                   sa, st, ea, et]

        assert len(data_in) == 23, data_in
        assert None not in data_in
        op2.write(struct1.pack(*data_in))
        op2_ascii.write(str(data_in) + '\n')

    op2.write(pack('i', nbytes))
    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
    return itable
    #"""PBUSH"""
    #ntotal = 92  # 23*4
    #

    #nentries = ndata // ntotal
    #assert nentries > 0, 'table={self.table_name} len={ndata - n}'
    #assert ndata % ntotal == 0, f'table={self.table_name} leftover=({ndata}%{ntotal}={ndata % ntotal}'

    #props = []
    #for unused_i in range(nentries):
        #edata = data[n:n+92]
        #out = struct1.unpack(edata)
        ## = out
        #prop = PBUSH.add_op2_data(out)
        #props.append(prop)
        #n += ntotal
    #return n, props

def write_pbarl(name, pids, itable, op2, op2_ascii, obj, endian=b'<'):
    """writes the PBARL"""
    key = (9102, 91, 52)
    fmt0 = endian + b'2i8s8sf'

    ndims = 0
    nproperties = len(pids)
    for pid in sorted(pids):
        prop = obj.properties[pid]
        ndim = len(prop.dim)
        ndims += ndim

    nvalues = 8 * nproperties + ndims + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))

    for pid in sorted(pids):
        prop = obj.properties[pid]

        # value is the first term in dim
        #(pid, mid, group, beam_type, value) = out
        bar_type = ('%-8s' % prop.beam_type).encode('ascii')
        group = ('%-8s' % prop.group).encode('ascii')
        data_in = [prop.pid, prop.mid, group, bar_type]

        ndim = len(prop.dim)
        fmti = b'%ifi' % ndim
        struct1 = Struct(fmt0 + fmti)
        data_in += prop.dim
        data_in.append(prop.nsm)
        data_in.append(-1)
        op2.write(struct1.pack(*data_in))
        op2_ascii.write(str(data_in) + '\n')

    op2.write(pack('i', nbytes))
    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
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

        # prop.nplies is total
        nplies = len(prop.mids)
        nlayers += nplies

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
            raise RuntimeError(f'unsupported ft.  pid={pid} ft={prop.ft!r}.'
                               f'\nPCOMP = {prop}')

        #is_symmetric = True
        # nplies is total
        nplies = len(prop.mids)
        symmetric_factor = 1
        if nplies != prop.nplies:
            assert nplies * 2 == prop.nplies
            symmetric_factor = -1
        data = [pid, symmetric_factor * nplies, prop.z0,
                prop.nsm, prop.sb, ft, prop.tref, prop.ge]
        op2.write(s1.pack(*data))
        op2_ascii.write(str(data) + '\n')

        for (mid, t, theta, sout) in zip(prop.mids, prop.thicknesses, prop.thetas, prop.souts):
            if sout == 'NO':
                sout = 0
            elif sout == 'YES':
                sout = 1
            else:
                raise RuntimeError(f'unsupported sout.  sout={sout!r} and must be 0 or 1.'
                                   f'\nPCOMP = {data}')
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

def write_pcompg(name, pids, itable, op2, op2_ascii, obj, endian=b'<'):
    """writes the PCOMPG"""
    key = (15006, 150, 604)

    nproperties = len(pids)
    nlayers = 0
    for pid in pids:
        #(pid, mid, a, j, c, nsm) = out
        prop = obj.properties[pid]
        #(pid, nlayers, z0, nsm, sb, ft, Tref, ge) = out # 8
        nlayers += prop.nplies

    # we add a layer for the (-1, -1, -1, -1, -1) at the end of each property
    nlayers += nproperties


    nvalues = 8 * nproperties + (5 * nlayers) + 3 # +3 comes from the keys
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

    s1 = Struct(endian + b'2i 3f i 2f')
    s2 = Struct(endian + b'ii 2f i')
    struct_i5 = Struct(endian + b'5i')

    lam_map = {
        None : 0,
    }
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
            raise RuntimeError(f'unsupported ft.  pid={pid} ft={prop.ft!r}.'
                               f'\nPCOMP = {prop}')

        #is_symmetric = True
        #symmetric_factor = 1
        lam = lam_map[prop.lam]
        #(pid, lam_int, z0, nsm, sb, ft_int, tref, ge) = out
        data = [pid, lam, prop.z0,
                prop.nsm, prop.sb, ft, prop.tref, prop.ge]
        op2.write(s1.pack(*data))
        op2_ascii.write(str(data) + '\n')

        for (glply, mid, t, theta, sout) in zip(prop.global_ply_ids, prop.mids, prop.thicknesses, prop.thetas, prop.souts):
            if sout == 'NO':
                sout = 0
            elif sout == 'YES':
                sout = 1
            else:
                raise RuntimeError(f'unsupported sout.  sout={sout!r} and must be 0 or 1.'
                                   f'\nPCOMPG = {data}')
            data2 = [glply, mid, t, theta, sout]
            op2.write(s2.pack(*data2))
            op2_ascii.write(str(data2) + '\n')
        data2 = [-1, -1, -1, -1, -1]
        op2.write(struct_i5.pack(*data2))
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

def _write_pbeam(name, model, pids, op2, op2_ascii, endian):
    struct1 = Struct(endian + b'4if')
    struct2 = Struct(endian + b'16f')
    struct3 = Struct(endian + b'16f')
    for pid in sorted(pids):
        nfieldsi = 0
        prop = model.properties[pid]
        nsegments = len(prop.xxb)
        if nsegments == 1:
            #  ccf = constant cross section flag
            ccf = 1 # True
            #print(prop.get_stats())
        elif nsegments == 2:
            #  ccf = constant cross section flag
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
                elif so_str == 'YESA':  # TODO: not sure...
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
    return

def _write_pbush1d(name, model, pids, spack, op2, op2_ascii, endian):
    type_map = {
        None : 0,  # NULL
        'EQUAT' : 1,
        'TABLE' : 2,
    }

    alpha = 0.
    for pid in sorted(pids):
        prop = model.properties[pid]
        #print(prop.get_stats())
        #(pid, k, c, m, unused_alpha, sa, se,
         #typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud,
         #types, idts, idcs, idtdus, idcdus,
         #typed, idtd, idcd, idtdvd, idcdvd,
         #typeg, idtg, idcg, idtdug, idcdug, idtdvg, idcdvg,
         #typef, idtf, idcf,
         #unused_ut, unused_uc) = out

        typea = cvt = cvc = expvt = expvc = idtsu = idtcu = idtsud = idcsud = 0
        if 'SHOCKA' in prop.vars:
            #optional_vars['SHOCKA'] = [typea_str, cvt, cvc, expvt, expvc,
                                       #idts, idets, idtcu, idtsud, idcsud]
            #shock_cvc : None
            #shock_cvt : 1000.0
            #shock_exp_vc : 1.0
            #shock_exp_vt : 1.0
            #shock_idecs : None
            #shock_idecsd : None
            #shock_idets : None
            #shock_idetsd : None
            #shock_idts : None
            #shock_type : 'TABLE'
            typea = type_map[prop.shock_type]
            cvt = prop.shock_cvc if prop.shock_cvc is not None else 0
            cvc = prop.shock_cvt if prop.shock_cvt is not None else 0
            expvt = prop.shock_exp_vc if prop.shock_exp_vc is not None else 0
            expvc = prop.shock_exp_vt if prop.shock_exp_vt is not None else 0
            idtsu = prop.shock_idts if prop.shock_idts is not None else 0
            idtcu = prop.shock_idecs if prop.shock_idecs is not None else 0
            idtsud = prop.shock_idetsd if prop.shock_idetsd is not None else 0
            idcsud = prop.shock_idecsd if prop.shock_idecsd is not None else 0

        types = idts = idcs = idtdus = idcdus = 0
        if 'SPRING' in prop.vars:
            #optional_vars['SPRING'] = [types_str, idts, idcs, idtdus, idcdus]
            #spring_idc : None
            #spring_idcdu : None
            #spring_idt : 205
            #spring_idtc : None
            #spring_idtdu : None
            #spring_type : 'TABLE'
            types = type_map[prop.spring_type]
            idts = prop.spring_idt if prop.spring_idt is not None else 0
            idcs = prop.spring_idc if prop.spring_idc is not None else 0
            idtdus = prop.spring_idtdu if prop.spring_idtdu is not None else 0
            idcdus = prop.spring_idcdu if prop.spring_idcdu is not None else 0

        typed = idtd = idcd = idtdvd = idcdvd = 0
        if 'DAMPER' in prop.vars:
            #optional_vars['DAMPER'] = [typed_str, idtd, idcd, idtdvd, idcdvd]
            #damper_idc : None
            #damper_idcdv : None
            #damper_idt : 206
            #damper_idtdv : None
            #damper_type : 'TABLE'
            typed = type_map[prop.damper_type]
            idtd = prop.damper_idt if prop.damper_idt is not None else 0
            idcd = prop.damper_idc if prop.damper_idc is not None else 0
            idtdvd = prop.damper_idtdv if prop.damper_idtdv is not None else 0
            idcdvd = prop.damper_idcdv if prop.damper_idcdv is not None else 0

        typeg = idtg = idcg = idtdug = idcdug = idtdvg = idcdvg = 0
        if 'GENER' in prop.vars:
            #typeg = type_map[typeg_str]
            #optional_vars['GENER'] = [idtg, idcg, idtdug, idcdug, idtdvg, idcdvg]
            gener

        typef = idtf = idcf = 0  #type_map[typef_str]  # FUSE...what is this???
        ut = uc = 0
        data = [pid, prop.k, prop.c, prop.m, alpha, prop.sa, prop.se,
                typea, cvt, cvc, expvt, expvc, idtsu, idtcu, idtsud, idcsud,
                types, idts, idcs, idtdus, idcdus,
                typed, idtd, idcd, idtdvd, idcdvd,
                typeg, idtg, idcg, idtdug, idcdug, idtdvg, idcdvg,
                typef, idtf, idcf,
                ut, uc]
        assert len(data) == 38, len(data)
        op2.write(spack.pack(*data))
    #ntotal = 152  # 38*4
    #struct1 = Struct(self._endian + b'i 6f i 4f 24i 2f')
    #nentries = (len(data) - n) // ntotal
    #for unused_i in range(nentries):
        #edata = data[n:n+152]
        #out = struct1.unpack(edata)
        ##  test_op2_other_05
        ##pbush1d, 204, 1.e+5, 1000., , , , , , +pb1
        ##+pb1, spring, table, 205, , , , , , +pb2
        ##+pb2, damper, table, 206
        ##pid=204 k=100000.0 c=1000.0 m=0.0 sa=nan se=nan


        #msg = f'PBUSH1D pid={pid} k={k} c={c} m={m} sa={sa} se={se}'
        #optional_vars = {}
        #typea_str = type_map[typea]
        #types_str = type_map[types]
        #typed_str = type_map[typed]
        #unused_typeg_str = type_map[typeg]
        #unused_typef_str = type_map[typef]

        #if min([typea, types, typed, typeg, typef]) < 0:
            #raise RuntimeError(f'typea={typea} types={types} typed={typed} typeg={typeg} typef={typef}')
    return
