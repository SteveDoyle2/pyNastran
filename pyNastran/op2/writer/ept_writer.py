from collections import defaultdict
from typing import BinaryIO
from struct import pack, Struct
import numpy as np

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import _write_spcadd
from pyNastran.op2.op2_interface.op2_reader import mapfmt


def write_ept(op2_file, op2_ascii, obj, endian=b'<',
              nastran_format: str='nx', size: int=4) -> None:
    if not hasattr(obj, 'properties'):
        return

    out = defaultdict(list)
    for pid, card in obj.convection_properties.items():
        out[card.type].append(card)
    for pid, phbdy in obj.phbdys.items():
        out[phbdy.type].append(pid)
    for pid, pelast in obj.pelast.items():
        out[pelast.type].append(pid)
    for pid, pdampt in obj.pdampt.items():
        out[pdampt.type].append(pid)
    for pid, pbusht in obj.pbusht.items():
        out[pbusht.type].append(pid)
    for nsm_id, nsmadds in obj.nsmadds.items():
        for card in nsmadds:
            out[card.type].append(card)

    #if not hasattr(obj, 'nodes'):
        #return
    nproperties = (len(obj.properties) + len(obj.properties_mass) +
                   len(obj.nsms)) + len(out)
    if nproperties == 0:
        return
    write_geom_header(b'EPT', op2_file, op2_ascii, endian=endian)
    # struct_3i = Struct(endian + b'3i')

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
    for pid, nsms in obj.nsms.items():
        for nsm in nsms:
            out[nsm.type].append(nsm)
    # for pid, prop in obj.nsmadds.items():
    #     out[prop.type].append(pid)

    skip_properties = [
        'PBEND', #'PBUSH1D',
        'PBCOMP',
        #'PGAP',
        #'PCOMPG',
    ]
    for name, pids in out.items():
        nproperties = len(pids)
        if nproperties == 0:  # pragma: no cover
            continue
        if name in skip_properties:  # pragma: no cover
            obj.log.warning(f'skipping EPT-{name}')
            continue
        #obj.log.debug('writing EPT-%s' % name)

        #print('EPT', itable, name)
        log = obj.log
        if name in EPT_MAP:
            func = EPT_MAP[name]
            # log.warning(f'reading EPT-{name}')
            itable = func(name, pids, itable, op2_file, op2_ascii, obj, endian=endian,
                          nastran_format=nastran_format)
            # itable = write_pfast(name, pids, itable, op2_file, op2_ascii, obj, endian=endian,
            #                      nastran_format=nastran_format)
            continue

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
        else:
            obj.log.warning(f'skipping {name}')
            continue
        #else:  # pragma: no cover
            #raise NotImplementedError(name)

        # doesn't include the key
        nvalues = nfields * nproperties

        nbytes = _write_table_header(
            op2_file, op2_ascii, name, key, nvalues, size)

        try:
            write_card(op2_file, op2_ascii, obj, name, pids, spack, endian)
        except Exception:
            obj.log.error('failed EPT-%s' % name)
            raise
        op2_file.write(pack('i', nbytes))
        itable -= 1

        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2_file.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2_file, op2_ascii, itable)
    #-------------------------------------

def write_card(op2_file, op2_ascii, obj, name, pids, spack, endian):
    op2_ascii.write('EPT-%s\n' % name)
    if name == 'PVISC':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, ce, cr) = out
            data = [pid, prop.ce, prop.cr]
            op2_file.write(spack.pack(*data))
    elif name == 'PBUSH1D':
        _write_pbush1d(name, obj, pids, spack, op2_file, op2_ascii, endian)
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
            op2_file.write(spack.pack(*data))

    elif name == 'PGAP':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            data = [pid, prop.u0, prop.f0, prop.ka, prop.kb, prop.kt,
                    prop.mu1, prop.mu2, prop.tmax, prop.mar, prop.trmin]
            assert None not in data, data
            op2_ascii.write('  pid=%s data=%s' % (pid, data[1:]))
            op2_file.write(spack.pack(*data))

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
            op2_file.write(spack.pack(*data))
    elif name == 'PBEAM':  # probably wrong stations
        _write_pbeam(name, obj, pids, op2_file, op2_ascii, endian)

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
            else:  # pragma: no cover
                raise NotImplementedError('prop.stress=%s and must be [0, 1]' % prop.stress)

            if prop.integ == 'BUBBLE':
                integ = 0
            elif prop.integ == 'GAUSS':
                integ = 1
            elif prop.integ == 'TWO':
                integ = 2
            elif prop.integ == 'THREE':
                integ = 3
            elif prop.integ is None:  # TODO: not sure
                integ = 0
            else:  # pragma: no cover
                raise NotImplementedError('prop.integ=%s and must be [0, 1, 2, 3]' % prop.integ)

            if prop.isop == 'REDUCED':
                isop = 0
            elif prop.isop == 'FULL':
                isop = 1
            elif prop.isop is None:  # TODO: not sure
                isop = 0
            elif prop.isop == 2: # 'TWO':
                isop = 2
            else:  # pragma: no cover
                raise NotImplementedError('isop=%r and must be [0, 1, 2]' % prop.isop)

            if prop.fctn == 'SMECH':
                fctn = b'SMEC'
            elif prop.fctn == 'PFLUID':
                fctn = b'PFLU'
            elif prop.fctn == 'FFLUID':
                fctn = b'FFLU'
            else:  # pragma: no cover
                raise NotImplementedError('PSOLID; fctn=%r' % prop.fctn)

            data = [pid, mid, cordm, integ, stress, isop, fctn]
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, mid, data[2:]))
            op2_file.write(spack.pack(*data))
    elif name == 'PSHEAR':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, t, nsm, f1, f2) = out
            data = [pid, prop.mid, prop.t, prop.nsm, prop.f1, prop.f2]
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2_file.write(spack.pack(*data))
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
            op2_file.write(spack.pack(*data))
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
            op2_file.write(spack.pack(*data))
    elif name == 'PROD':
        for pid in sorted(pids):
            prop = obj.properties[pid]
            #(pid, mid, a, j, c, nsm) = out
            data = [pid, prop.mid, prop.A, prop.j, prop.c, prop.nsm]
            op2_ascii.write('  pid=%s mid=%s data=%s\n' % (pid, prop.mid, data[2:]))
            op2_file.write(spack.pack(*data))
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
            op2_file.write(spack.pack(*data))
    elif name == 'PMASS':
        for pid in sorted(pids):
            prop = obj.properties_mass[pid]
            data = [pid, prop.mass]
            op2_ascii.write('  pid=%s mass=%s\n' % (pid, prop.mass))
            op2_file.write(spack.pack(*data))
    elif name == 'PELAST':
        for pid in sorted(pids):
            prop = obj.pelast[pid]
            data = [pid, prop.tkid, prop.tgeid, prop.tknid]
            op2_ascii.write('  pid=%s tables=%s\n' % (pid, data[1:]))
            op2_file.write(spack.pack(*data))
    elif name == 'PHBDY':
        for pid in sorted(pids):
            prop = obj.phbdys[pid]
            af = 0.0 if prop.af is None else prop.af
            d1 = 0.0 if prop.d1 is None else prop.d1
            d2 = 0.0 if prop.d2 is None else prop.d2
            data = [pid, af, d1, d2]
            op2_ascii.write('  pid=%s [af,d1,d2]=%s\n' % (pid, data[1:]))
            #print('  pid=%s [af,d1,d2]=%s\n' % (pid, data[1:]))
            op2_file.write(spack.pack(*data))
            #(pid, af, d1, d2) = out
    else:  # pragma: no cover
        raise NotImplementedError(name)


def write_pbush(name: str, pids: np.ndarray, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx') -> int:
    """writes the PBUSH"""
    size = 4
    key = (1402, 14, 37)
    is_mass = False
    is_alpha = False
    if nastran_format == 'nx':
        nfields = 23
        fmt = mapfmt(endian + b'i22f', size)
    elif nastran_format == 'msc':
        # TODO: there are 3 different types of PBUSH cards...
        # if ndata == 23:
        #     (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
        #      g1, g2, g3, g4, g5, g6, sa, st, ea, et) = data
        #     mass = 0.
        # elif ndata == 24:
        #     (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
        #      g1, g2, g3, g4, g5, g6, sa, st, ea, et,
        #      mass) = data
        # elif ndata == 27:
        #     (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
        #      g1, g2, g3, g4, g5, g6, sa, st, ea, et,
        #      mass, alpha, tref, coinl) = data
        for pid in sorted(pids):
            prop = obj.properties[pid]
            if prop.mass is not None:
                is_mass = True
            if prop.alpha is not None:
                is_mass = True
                is_alpha = True
                break
        if is_alpha:
            nfields = 27
            fmt = mapfmt(endian + b'i22f 4f', size)
        elif is_mass:
            nfields = 24
            fmt = mapfmt(endian + b'i22f f', size)
        else:
            nfields = 23
            fmt = mapfmt(endian + b'i22f', size)
    else:
        raise RuntimeError(f'EPT writer: PBUSH nastran_format={nastran_format}')


    struct1 = Struct(fmt)
    nproperties = len(pids)

    # nvalues = nfields * nproperties + 3 # +3 comes from the keys
    nvalues = nfields * nproperties

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nvalues, size)
    if 0:
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

    for pid in sorted(pids):
        prop = obj.properties[pid]
        (b1, b2, b3, b4, b5, b6) = (0., 0., 0., 0., 0., 0.)
        (k1, k2, k3, k4, k5, k6) = (0., 0., 0., 0., 0., 0.)
        (g1, g2, g3, g4, g5, g6) = (0., 0., 0., 0., 0., 0.)
        if prop.b:
            (b1, b2, b3, b4, b5, b6) = [
                bi if bi is not None else 0.0 for bi in prop.b]  # damping
        if prop.k:
            (k1, k2, k3, k4, k5, k6) = [
                ki if ki is not None else 0.0 for ki in prop.k] # stiffness
        if prop.ge:
            (g1, g2, g3, g4, g5, g6) = [
                gi if gi is not None else 0.0 for gi in prop.ge] # ???

        sa = prop.sa if prop.sa is not None else 0.
        st = prop.st if prop.st is not None else 0.
        ea = prop.ea if prop.ea is not None else 0.
        et = prop.et if prop.et is not None else 0.
        data_in = [pid,
                   k1, k2, k3, k4, k5, k6,
                   b1, b2, b3, b4, b5, b6,
                   g1, g2, g3, g4, g5, g6,
                   sa, st, ea, et]
        if is_alpha:
            data_in.extend([prop.alpha, prop.tref, prop.coincident_length])
            assert len(data_in) == 27, data_in
        elif is_mass:
            data_in.append(prop.mass)
            assert len(data_in) == 24, data_in
        else:
            pass
            #data_in = [pid,
                       #k1, k2, k3, k4, k5, k6,
                       #b1, b2, b3, b4, b5, b6,
                       #g1, g2, g3, g4, g5, g6,
                       #sa, st, ea, et]
            assert len(data_in) == 23, data_in

        assert None not in data_in
        op2_file.write(struct1.pack(*data_in))
        op2_ascii.write(str(data_in) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_pfast(name: str, pids: np.ndarray, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx') -> int:
    """writes the PFAST"""
    # TODO: there are 2 different types of PFAST cards...

    nproperties = len(pids)
    if nastran_format == 'nx':
        key = (3601, 36, 55)
        fmt = endian + b'ifii 8f'
        nvalues_per_property = 12
    elif nastran_format == 'msc':
        key = (13501, 135, 510)
        nvalues_per_property = 25
        fmt = endian + b'2if 5i 2f2i2f 3i 2i 6f'
    else:  # pragma: no cover
        raise NotImplementedError(nastran_format)
    nvalues = nvalues_per_property * nproperties + 3 # +3 comes from the keys
    struct1 = Struct(fmt)
    nbytes = nvalues * 4
    op2_file.write(pack('3i', *[4, nvalues, 4]))
    op2_file.write(pack('i', nbytes)) #values, nbtyes))

    op2_file.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    for pid in sorted(pids):
        prop = obj.properties[pid]
        #print(prop.get_stats())
        if nastran_format == 'nx':
            data_in = [pid, prop.d, prop.mcid, prop.mflag,
                       prop.kt1, prop.kt2, prop.kt3,
                       prop.kr1, prop.kr2, prop.kr3,
                       prop.mass, prop.ge]
        elif nastran_format == 'msc':
            #4 CONNBEH   I Connection behavior (0=FF/F, 1=FR, 10=RF/R, 11=RR)
            #5 CONNTYPE  I Connection type (0=clamp, 1=hinge, 2=bolt)
            #6 EXTCON    I External constraint flag (0=off, 1=on)
            #7 CONDTYPE  I Condition type (0=rigid, 1=equivalent)
            #8 WELDTYPE  I Weld type (0=spot weld, 1=but seam, 2=T-seam)
            connbeh = 0
            conntype = 0
            weldtype = 0
            extcon = 0
            condtype = 0

            #9 MINLEN   RS Minimum length of spot weld
            #10 MAXLEN  RS Maximum length of spot weld
            #11 GMCHK    I Perform geometry check
            #12 SPCGS    I SPC the master grid GS
            gmcheck = 0
            spcgs = 0
            minlen = 0.0
            maxlen = 0.0

            blank1 = 0
            blank2 = 0
            blank3 = 0
            #pid mid  D    con  con  ext  cond weld min max  chk  spc  cmass ge  und  und  und  mcid mfag kt1      kt2       kt3       kr1    kr2      kr3
            mcid0 = prop.mcid
            data_in = [
                # this 0 value is really weird...
                pid, 0, prop.d, connbeh, conntype, extcon,
                condtype, weldtype, minlen, maxlen,
                gmcheck, spcgs, prop.mass, prop.ge,
                blank1, blank2, blank3, prop.mcid, prop.mflag,
                prop.kt1, prop.kt2, prop.kt3,
                prop.kr1, prop.kr2, prop.kr3]
        else:  # pragma: no cover
            raise NotImplementedError(nastran_format)

        assert None not in data_in, data_in
        op2_file.write(struct1.pack(*data_in))
        op2_ascii.write(str(data_in) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_nsmadd(card_type: str, cards: list, itable: int,
                 op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                 nastran_format: str='nx', size: int=4) -> int:
    ncards = len(cards)
    nbytes = _write_spcadd(card_type, cards, ncards, op2_file, op2_ascii,
                           endian)
    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_nsm(name: str, nsms: list, itable: int,
              op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
              nastran_format: str='nx', size: int=4) -> int:
    """
    NX 2019.2
    RECORD – NSML(3201,32,55)

    Word Name Type Description
    1 SID         I Set identification number
    2 PROP(2) CHAR4 Set of properties or elements
    4 ORIGIN      I  Entry origin
    5 ID          I  Property or element identification number
    6 VALUE      RS Nonstructural mass value
    Words 5 through 6 repeat until End of Record

      ints    = (3, ELEMENT, 0,    200, 0.7, -1,
                 4, PSHELL,  0,   6401, 4.2, -1)
      floats  = (3, ELEMENT, 0.0,  200, 0.7, -1,
                 4, PSHELL,  0.0, 6401, 4.2, -1)

    id     : 10
    ids    : [10]
    nsm_type : 'PSHELL'
    sid    : 3000
    value  : 1.0
    """
    key = (3201, 32, 55)

    nfieldsi = 0
    for nsm in nsms:
        nfieldsi += 5 + 2 * len(nsm.ids)
        assert len(nsm.ids) == 1, nsm.get_stats()

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)

    if 0:
        nvalues = nfieldsi + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

    # nfields = 0
    for nsm in nsms:
        # print(nsm.get_stats())
        # sid    : 1000
        # nsm_type : 'ELEMENT'
        # value  : 1.0
        # ids    : [1]
        # if nsm.nsm_type == 'ELEMENT':
        #     nsm_type_bytes = b'ELEM'
        # else:
        nsm_type8 = f'{nsm.nsm_type:<8}'
        nsm_type_bytes = nsm_type8.encode('latin1')

        # ints = (3, ELEMENT, 0, 200, 0.7, -1,
        #         4, PSHELL, 0, 6401, 4.2, -1)
        fmt = b'i8si'
        data = [nsm.sid, nsm_type_bytes, 0, ]
        value = nsm.value
        for idi in nsm.ids:
            data.extend([idi, value])
            fmt += b'if'
        fmt += b'i'
        data.append(-1)

        struct1 = Struct(fmt)
        op2_file.write(struct1.pack(*data))
        op2_ascii.write(str(data) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_nsml(name: str, nsms: list, itable: int,
               op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
               nastran_format: str='nx', size: int=4) -> int:
    """
    NX 2019.2
    RECORD – NSML(3501, 35, 994)

    Defines a set of lumped nonstructural mass by ID.
    Word Name Type Description
    1 SID         I Set identification number
    2 PROP(2) CHAR4 Set of properties or elements
    4 ID          I Property of element identification number
    5 VALUE      RS Lumped nonstructural mass value
    Words 4 and 5 repeat until -1 occurs

      ints    = (3, ELEMENT, 0,    200, 0.7, -1,
                 4, PSHELL,  0,   6401, 4.2, -1)
      floats  = (3, ELEMENT, 0.0,  200, 0.7, -1,
                 4, PSHELL,  0.0, 6401, 4.2, -1)

    id     : 10
    ids    : [10]
    nsm_type : 'PSHELL'
    sid    : 3000
    value  : 1.0
    """
    key = (3501, 35, 994)

    nfieldsi = 0
    for nsm in nsms:
        nfieldsi += 5 + 2 * len(nsm.ids)
        assert len(nsm.ids) == 1, nsm.get_stats()

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)
    if 0:
        nvalues = nfieldsi + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

    # nfields = 0
    for nsm in nsms:
        # print(nsm.get_stats())
        # sid    : 1000
        # nsm_type : 'ELEMENT'
        # value  : 1.0
        # ids    : [1]
        # if nsm.nsm_type == 'ELEMENT':
        #     nsm_type_bytes = b'ELEM'
        # else:
        nsm_type8 = f'{nsm.nsm_type:<8}'
        nsm_type_bytes = nsm_type8.encode('latin1')

        # ints = (3, ELEMENT, 0, 200, 0.7, -1,
        #         4, PSHELL, 0, 6401, 4.2, -1)
        fmt = b'i8si'
        data = [nsm.sid, nsm_type_bytes, 0, ]
        value = nsm.value
        for idi in nsm.ids:
            data.extend([idi, value])
            fmt += b'if'
        fmt += b'i'
        data.append(-1)

        struct1 = Struct(fmt)
        op2_file.write(struct1.pack(*data))
        op2_ascii.write(str(data) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_nsm1(name: str, nsms: list, itable: int,
               op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
               nastran_format: str='nx', size: int=4) -> int:
    """
    Writes the NX cards:
      NSM1(3301, 33, 992)

    Defines the properties of a nonstructural mass.
     Word Name Type Description
     1 SID      I Set identification number
     2 PROP CHAR4 Set of properties
     3 TYPE CHAR4 Set of elements
     4 ORIGIN   I Entry origin
     5 VALUE   RS Nonstructural mass value
     6 SPECOPT  I Specification option
    """
    key = (3301, 33, 992)
    # fmt0 = endian + b'2i8s 8sf'

    nfieldsi = 0
    for nsm in nsms:
        if nsm.ids == ['ALL']:  # SPECOPT=2
            #  1     2      3      4         5       6      7      8   9
            #[sid, 'ELEM, 'ENT ', value, SPECOPT, 'ALL ', '    ', -1, -2]
            #fmt += b'i8sf i8s'  # 7 fields here; fmt0

            #
            #data.extend([2, b'ALL     '], -1)
            nfieldsi += 9
        else:
            nidsi = len(nsm.ids) + 4  # 1 flag; -1, -2
            nfieldsi += 3 + nidsi

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)
    if 0:
        nvalues = nfieldsi + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

    # nfields = 0
    # all_data = []
    # minus2_bytes = pack('i', -2)
    for nsm in nsms:
        # sid    : 1000
        # nsm_type : 'ELEMENT'
        # value  : 1.0
        # ids    : [1]
        # if nsm.nsm_type == 'ELEMENT':
        #     nsm_type_bytes = b'ELEM'
        # else:
        nsm_type8 = f'{nsm.nsm_type:<8}'
        nsm_type_bytes = nsm_type8.encode('latin1')

        data = [nsm.sid, nsm_type_bytes, 0, nsm.value]
        if nsm.ids == ['ALL']:  # SPECOPT=2
            fmt = b'i8sif i8s i'
            data.extend([2, b'ALL     ', -1])
        else:
            # SPECOPT=1 By IDs
            #   6 ID I Property of element identification number
            #   Word 6 repeats until -1 occurs
            nidsi = len(nsm.ids) + 2  # 3=1 for specopt, 2 for -1/-1
            fmt = f'i8sif {nidsi}i'.encode('latin1')
            assert isinstance(nsm.ids, list), nsm.ids
            data.append(1)
            data.extend(nsm.ids)
            data.extend([-1])

        # fmt += b'i'
        # data.append(-2)  # end of table
        # all_data.extend(data)
        struct1 = Struct(fmt)

        # print(fmt, nfieldsi)
        # print(data)
        op2_file.write(struct1.pack(*data))
        op2_ascii.write(str(data) + '\n')

    # print(all_data)
    # assert len(all_data) == nfieldsi, (nfieldsi, len(all_data), all_data)
    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable


def write_nsml1(name: str, nsms: list, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx', size: int=4) -> int:
    """
    Writes the NX cards:
      NSML1(3701, 37, 995)
    lumped nonstructural mass entries by VALUE, ID list.

    Word Name Type Description
    1 SID      I Set identification number
    2 PROP CHAR4 Set of properties
    3 TYPE CHAR4 Set of elements
    4 VALUE   RS Lumped nonstructural mass value
    5 SPECOPT  I Specification option
    SPECOPT=1 By IDs
      6 ID I Property of element identification number
      Word 6 repeats until -1 occurs
    SPECOPT=2 All
      6 ALL(2) CHAR4 Keyword ALL
      Words 6 and 7 repeat until -1 occurs
    SPECOPT=3 Thru range
      6 ID1         I Starting identification number
      7 THRU(2) CHAR4 Keyword THRU
      9 ID2         I Ending identification number
      Words 6 through 9 repeat until -1 occurs
    SPECOPT=4 Thru range with by
      6 ID1         I Starting identification number
      7 THRU(2) CHAR4 Keyword THRU
      9 ID2         I Ending identification number
      10 BY(2)  CHAR4 Keyword BY
      12 N I Increment
      Words 6 through 12 repeat until -1 occurs

    data = (1, ELEMENT, 466.2,
            3, 249311, THRU, 250189, -1,
            3, 250656, THRU, 251905, -1,
            3, 270705, THRU, 275998, -1,
            3, 332687, THRU, 334734, -1,
            -2,)
    [1000, b'ELEMENT ', 1.0,
        1, 42, -1,
     1001, b'ELEMENT ', 1.0, 1, 42, -1, -1]
    """
    key = (3701, 37, 995)
    # fmt0 = endian + b'2i8s 8sf'

    nfieldsi = 0
    for nsm in nsms:
        if nsm.ids == ['ALL']:  # SPECOPT=2
            #  1     2      3      4         5       6      7      8   9
            #[sid, 'ELEM, 'ENT ', value, SPECOPT, 'ALL ', '    ', -1, -2]
            #fmt += b'i8sf i8s'  # 7 fields here; fmt0

            #
            #data.extend([2, b'ALL     '], -1)
            nfieldsi += 9
        else:
            nidsi = len(nsm.ids) + 4  # 1 flag; -1, -2
            nfieldsi += 3 + nidsi

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)
    if 0:
        nvalues = nfieldsi + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
        op2_ascii.write('%s %s\n' % (name, str(key)))

    # nfields = 0
    # all_data = []
    # minus2_bytes = pack('i', -2)
    for nsm in nsms:
        # sid    : 1000
        # nsm_type : 'ELEMENT'
        # value  : 1.0
        # ids    : [1]
        # if nsm.nsm_type == 'ELEMENT':
        #     nsm_type_bytes = b'ELEM'
        # else:
        nsm_type8 = f'{nsm.nsm_type:<8}'
        nsm_type_bytes = nsm_type8.encode('latin1')

        data = [nsm.sid, nsm_type_bytes, nsm.value]
        if nsm.ids == ['ALL']:  # SPECOPT=2
            fmt = b'i8sf i8s i'
            data.extend([2, b'ALL     ', -1])
        else:
            # SPECOPT=1 By IDs
            #   6 ID I Property of element identification number
            #   Word 6 repeats until -1 occurs
            nidsi = len(nsm.ids) + 2  # 3=1 for specopt, 2 for -1/-1
            fmt = f'i8sf {nidsi}i'.encode('latin1')
            assert isinstance(nsm.ids, list), nsm.ids
            data.append(1)
            data.extend(nsm.ids)
            data.extend([-1])

        fmt += b'i'
        data.append(-2)  # end of table
        # all_data.extend(data)
        struct1 = Struct(fmt)
        op2_file.write(struct1.pack(*data))
        op2_ascii.write(str(data) + '\n')

    # print(all_data)
    # assert len(all_data) == nfieldsi, (nfieldsi, len(all_data), all_data)
    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_pbarl(name: str, pids: np.ndarray, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx', size: int=4) -> int:
    """writes the PBARL"""
    key = (9102, 91, 52)
    fmt0 = endian + b'2i8s8sf'

    ndims = 0
    nproperties = len(pids)
    for pid in sorted(pids):
        prop = obj.properties[pid]
        ndim = len(prop.dim)
        ndims += ndim

    nvalues = 8 * nproperties + ndims
    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nvalues, size)
    if 0:
        nvalues = 8 * nproperties + ndims + 3 # +3 comes from the keys
        nbytes = nvalues * 4
        op2_file.write(pack('3i', *[4, nvalues, 4]))
        op2_file.write(pack('i', nbytes)) #values, nbtyes))

        op2_file.write(pack('3i', *key))
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
        op2_file.write(struct1.pack(*data_in))
        op2_ascii.write(str(data_in) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_pcomp(name: str, pids: np.ndarray, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx') -> int:
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
    op2_file.write(pack('3i', *[4, nvalues, 4]))
    op2_file.write(pack('i', nbytes)) #values, nbtyes))

    op2_file.write(pack('3i', *key))
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
        elif prop.ft == 'HFAI':  # secret MSC
            ft = 5
        elif prop.ft == 'HTAP':  # secret MSC
            ft = 6
        elif prop.ft == 'HFAB':  # secret MSC
            ft = 7
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
        op2_file.write(s1.pack(*data))
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
            op2_file.write(s2.pack(*data2))
            op2_ascii.write(str(data2) + '\n')


    #data_in = [
        #pid, z0, nsm, sb, ft, Tref, ge,
        #is_symmetrical, mids, T, thetas, souts]

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_pcompg(name: str, pids: np.ndarray, itable: int,
                 op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                 nastran_format: str='nx') -> int:
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
    op2_file.write(pack('3i', *[4, nvalues, 4]))
    op2_file.write(pack('i', nbytes)) #values, nbtyes))

    op2_file.write(pack('3i', *key))
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
    ft_to_int_map = {
        None: 0,
        'HILL': 1,
        'HOFF': 2,
        'TSAI': 3,
        'STRN': 4,
    }

    for pid in sorted(pids):
        prop = obj.properties[pid]
        #print(prop.get_stats())

        try:
            ft = ft_to_int_map[prop.ft]
        except KeyError:
            raise KeyError(f'unsupported ft.  pid={pid} ft={prop.ft!r}.'
                           f'\nPCOMP = {prop}')

        #is_symmetric = True
        #symmetric_factor = 1
        lam = lam_map[prop.lam]
        #(pid, lam_int, z0, nsm, sb, ft_int, tref, ge) = out
        data = [pid, lam, prop.z0,
                prop.nsm, prop.sb, ft, prop.tref, prop.ge]
        op2_file.write(s1.pack(*data))
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
            op2_file.write(s2.pack(*data2))
            op2_ascii.write(str(data2) + '\n')
        data2 = [-1, -1, -1, -1, -1]
        op2_file.write(struct_i5.pack(*data2))
        op2_ascii.write(str(data2) + '\n')

    #data_in = [
        #pid, z0, nsm, sb, ft, Tref, ge,
        #is_symmetrical, mids, T, thetas, souts]

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable

def write_pconv(name: str, props: list, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx', size: int=4) -> int:
    """
    PCONV(11001,110,411)- MSC/NX

    """
    key = (11001, 110, 411)

    if nastran_format == 'nx':
        nfieldsi = len(props) * 4
    else:
        # msc
        nfieldsi = len(props) * 14

    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)

    # nfields = 0

    if nastran_format == 'nx':
        struct_3if = Struct(endian + b'3if')
        for prop in props:
            mid = 0 if prop.mid is None else prop.mid
            data = [prop.pconid, mid, prop.form, prop.expf]
            op2_file.write(struct_3if.pack(*data))
            op2_ascii.write(str(data) + '\n')
    else:
        # msc
        struct1 = Struct(endian + b'3if 4i fii 3f')
        for prop in props:
            mid = 0 if prop.mid is None else prop.mid
            tid = 0 if prop.tid is None else prop.tid
            e1 = 0.0 if prop.e[0] is None else prop.e[0]
            e2 = 0.0 if prop.e[1] is None else prop.e[1]
            e3 = 0.0 if prop.e[2] is None else prop.e[2]
            chlen = 0.0 if prop.chlen is None else prop.chlen
            gidin = 0 if prop.gidin is None else prop.gidin
            ftype = 0 if prop.ftype is None else prop.ftype
            ce = 0 if prop.ce is None else prop.ce
            data = [prop.pconid, mid, prop.form, prop.expf,
                    ftype, tid, 0, 0,
                    chlen, gidin, ce,
                    e1, e2, e3]
            assert None not in data, (data, data.index(None))
            op2_file.write(struct1.pack(*data))
            op2_ascii.write(str(data) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable


def _write_pbeam(name, model, pids, op2_file, op2_ascii, endian):
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
        op2_file.write(struct1.pack(*data))

        nzero_segments = nsegments - 1
        #print(f'nsegments={nsegments} nzero_segments={nzero_segments}')
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
            op2_file.write(struct2.pack(*data))

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
        op2_file.write(struct3.pack(*data))
        assert nfieldsi == 197, nfieldsi
    return

def _write_pbush1d(name, model, pids, spack, op2_file, op2_ascii, endian):
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
        op2_file.write(spack.pack(*data))
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

def write_pelas(name: str, props: list, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx', size: int=4) -> int:
    key = (302, 3, 46)

    nfieldsi = len(props) * 4
    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)

    spack = Struct(endian + b'i3f')
    for pid in sorted(props):
        prop = obj.properties[pid]
        # (pid, k, ge, s) = out
        data = [pid, prop.k, prop.ge, prop.s]
        # assert None not in data, (data, data.index(None))
        op2_file.write(spack.pack(*data))
        op2_ascii.write(str(data) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable


def write_pdamp(name: str, props: list, itable: int,
                op2_file: BinaryIO, op2_ascii, obj, endian: bytes=b'<',
                nastran_format: str='nx', size: int=4) -> int:
    key = (202, 2, 45)
    spack = Struct(endian + b'if')

    nfieldsi = len(props) * 2
    nbytes = _write_table_header(
        op2_file, op2_ascii, name, key, nfieldsi, size)

    for pid in sorted(props):
        prop = obj.properties[pid]
        # (pid, b) = out
        data = [pid, prop.b]
        # assert None not in data, (data, data.index(None))
        op2_file.write(spack.pack(*data))
        op2_ascii.write(str(data) + '\n')

    itable = _write_table_footer(op2_file, op2_ascii, nbytes, itable)
    return itable


def _write_table_header(op2_file: BinaryIO, op2_ascii,
                        name: str, key: tuple[int, int, int],
                        nfields: int, size: int) -> int:
    nvalues = nfields + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2_file.write(pack('3i', *[4, nvalues, 4]))
    op2_file.write(pack('i', nbytes)) #values, nbtyes))

    op2_file.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes


def _write_table_footer(op2_file: BinaryIO, op2_ascii,
                        nbytes: int, itable: int) -> int:
    op2_file.write(pack('i', nbytes))
    itable -= 1
    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2_file.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')
    return itable

EPT_MAP = {
    'NSMADD': write_nsmadd,
    'NSM': write_nsm,
    'NSM1': write_nsm1,

    'NSML': write_nsml,
    'NSML1': write_nsml1,

    'PBARL': write_pbarl,
    'PCOMP': write_pcomp,
    'PCOMPG': write_pcompg,
    'PBUSH': write_pbush,
    'PFAST': write_pfast,
    'PCONV': write_pconv,
    'PELAS': write_pelas,
    'PDAMP': write_pdamp,
}
