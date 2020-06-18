"""writes the MPT/MPTS table"""
from collections import defaultdict
from struct import pack, Struct

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header

def write_mpt(op2, op2_ascii, model, endian=b'<'):
    """writes the MPT/MPTS table"""
    if not hasattr(model, 'materials'):
        return
    nmaterials = len(model.materials) + len(model.thermal_materials)

    # the code will crash if a something is in the out dict, but
    # not handled properly
    materials_to_skip = [
        'MAT11',
        'MATS3', 'MATS8',
        'MATT3', 'MATT9',
        #  other
        'NLPARM', 'NLPCI', 'TSTEPNL', 'MAT3D',
        'RADBC',
    ]
    # these are specifically for material dependency objects (e.g., model.MATS1)
    # in general, materials (e.g., MAT1) just need to be removed from the
    # skip list to become active.  If it's a material dependency type,  it must
    # be in this list.
    mtypes = [
        'MATS1', 'MATS3', 'MATS8',
        'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
    ]
    out = defaultdict(list)
    for mtype in mtypes:
        mat_dict = getattr(model, mtype)
        if len(mat_dict) == 0:
            continue
        if mtype in materials_to_skip:
            model.log.warning(f'skipping MPT-{mtype}')
            continue
        out[mtype] = list(mat_dict.keys())
    if len(model.nlparms):
        out['NLPARM'] = model.nlparms.keys()
    if len(model.nlpcis):
        out['NLPCI'] = model.nlpcis.keys()

    if len(model.bcs):
        #'bcs' : ['CONV', 'CONVM', 'RADBC', 'RADM', 'TEMPBC'],
        for key, objs in model.bcs.items():
            if isinstance(objs, list):
                for i, obji in enumerate(objs):
                    if obji.type == 'RADBC':
                        out[obji.type].append((key, i))
            else:
                raise TypeError(objs)  # weird...was this happening or was I messing with things wrong...
                #print('obj...')
                #print(typei)
                #print(key, objs, type(objs))
            #print(out['RADBC'])
            del objs

    if len(model.tstepnls):
        # strangely TSTEPs are stored in the DYNAMIC table, but TSTEPNLs are
        # stored in the MPT table
        #'tstepnls' : ['TSTEPNL', 'TSTEP1'],
        out['TSTEPNL'] = model.tstepnls.keys()


    if nmaterials == 0 and len(out) == 0:
        return
    write_geom_header(b'MPT', op2, op2_ascii, endian=endian)
    itable = -3

    for mid, mat in model.materials.items():
        out[mat.type].append(mat.mid)
    for mid, mat in model.thermal_materials.items():
        out[mat.type].append(mat.mid)


    for name, mids in sorted(out.items()):
        #model.log.debug('MPT %s %s' % (name, mids))
        nmaterials = len(mids)
        if name in materials_to_skip:
            model.log.warning(f'skipping MPT-{name}')
            continue

        #if nmaterials == 0:
            #continue
        if name == 'MAT1':
            nbytes = _write_mat1(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT2':
            nbytes = _write_mat2(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT3':
            nbytes = _write_mat3(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT4':
            nbytes = _write_mat4(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT5':
            nbytes = _write_mat5(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT8':
            nbytes = _write_mat8(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT9':
            nbytes = _write_mat9(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MAT10':
            nbytes = _write_mat10(model, name, mids, nmaterials, op2, op2_ascii, endian)

        elif name == 'MATS1':
            nbytes = _write_mats1(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MATT1':
            nbytes = _write_matt1(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MATT2':
            nbytes = _write_matt2(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MATT4':
            nbytes = _write_matt4(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MATT5':
            nbytes = _write_matt5(model, name, mids, nmaterials, op2, op2_ascii, endian)
        elif name == 'MATT8':
            nbytes = _write_matt8(model, name, mids, nmaterials, op2, op2_ascii, endian)

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
    close_geom_table(op2, op2_ascii, itable)
    #-------------------------------------

def _write_mat1(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT1"""
    key = (103, 1, 77)
    nfields = 12
    spack = Struct(endian + b'i10fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
        data = [mid, mat.e, mat.g, mat.nu, mat.rho, mat.a, mat.tref,
                mat.ge, mat.St, mat.Sc, mat.Ss, mat.mcsid]
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat2(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT2"""
    key = (203, 2, 78)
    nfields = 17
    spack = Struct(endian + b'i15fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
         #tref, ge, St, Sc, Ss, mcsid) = out
        data = [
            mid, mat.G11, mat.G12, mat.G13, mat.G22, mat.G23, mat.G33,
            mat.rho,
            0.0 if mat.a1 is None else mat.a1,
            0.0 if mat.a2 is None else mat.a2,
            0.0 if mat.a3 is None else mat.a3,
            mat.tref, mat.ge,
            0.0 if mat.St is None else mat.St,
            0.0 if mat.Sc is None else mat.Sc,
            0.0 if mat.Ss is None else mat.Ss,
            0 if mat.mcsid is None else mat.mcsid]

        #print(data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        assert None not in data, 'MAT2 %s' % data
        #print('MAT2', data)
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat3(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT3"""
    key = (1403, 14, 122)
    nfields = 16
    spack = Struct(endian + b'i8fi5fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
         #blank, ax, ath, az, tref, ge, blank) = out

        data = [
            mid,
            mat.ex, mat.eth, mat.ez, mat.nuxth, mat.nuthz, mat.nuzx,
            mat.rho,
            0.0 if mat.gzx is None else mat.gzx,
            0, mat.ax, mat.ath, mat.az, mat.tref,
            mat.ge, 0]
        assert None not in data, 'MAT3 %s' % data
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat4(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT4"""
    key = (2103, 21, 234)
    nfields = 11
    spack = Struct(endian + b'i10f')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.thermal_materials[mid]
        #(mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
        data = [
            mid,
            mat.k, mat.cp, mat.rho,
            #  not sure
            0.0 if mat.H is None else mat.H,
            0.0 if mat.mu is None else mat.mu,
            mat.hgen,
            #  not sure
            0.0 if mat.ref_enthalpy is None else mat.ref_enthalpy,
            0.0 if mat.tch is None else mat.tch,
            0.0 if mat.tdelta is None else mat.tdelta,
            0.0 if mat.qlat is None else mat.qlat,
        ]
        #print('MAT4 -', data)
        #print(data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat5(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT5"""
    key = (2203, 22, 235)
    nfields = 10
    spack = Struct(endian + b'i9f')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.thermal_materials[mid]
        #(mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
        data = [
            mid,
            mat.kxx, mat.kxy, mat.kxz, mat.kyy, mat.kyz, mat.kzz, # 6
            mat.cp, mat.rho, mat.hgen,
        ]
        #print('MAT5 -', data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat8(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT8"""
    key = (2503, 25, 288)
    nfields = 19
    spack = Struct(endian + b'i18f')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
        # tref, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
        data = [mid, mat.e11, mat.e22, mat.nu12, mat.g12, mat.g1z, mat.g2z,
                mat.rho, mat.a1, mat.a2, mat.tref,
                mat.Xt, mat.Xc, mat.Yt, mat.Yc, mat.S,
                mat.ge, mat.F12, mat.strn]

        #print('MAT8 -', data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat9(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT9"""
    key = (2603, 26, 300)
    nfields = 35
    spack = Struct(endian + b'i 30f 4i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
         #g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
         #rho, a1, a2, a3, a4, a5, a6, tref, ge,
         #blank1, blank2, blank3, blank4) = out
        data = [
            mid,
            mat.G11, mat.G12, mat.G13, mat.G14, mat.G15, mat.G16,
            mat.G22, mat.G23, mat.G24, mat.G25, mat.G26,
            mat.G33, mat.G34, mat.G35, mat.G36,
            mat.G44, mat.G45, mat.G46,
            mat.G55, mat.G56, mat.G66,
            mat.rho] + mat.A + [mat.tref, mat.ge, 0, 0, 0, 0,]
        #assert len(mat.A) == 6, mat.A
        assert len(data) == nfields

        #print('MAT9 -', data, len(data))
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mat10(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MAT10"""
    key = (2801, 28, 365)
    nfields = 5
    spack = Struct(endian + b'i4f')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, bulk, rho, c, ge) = out
        data = [mid, mat.bulk, mat.rho, mat.c, mat.ge]
        assert len(data) == nfields

        #print('MAT10 -', data, len(data))
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_mats1(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATS1"""
    key = (503, 5, 90)
    nfields = 11
    spack = Struct(endian + b'3ifiiff3i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATS1[mid]
        #(mid, tid, Type, h, yf, hr, limit1, limit2, a, bmat, c) = out
        a = 0
        bmat = 0
        c = 0

        if mat.Type == 'NLELAST':
            Type = 1
        elif mat.Type == 'PLASTIC':
            Type = 2
        elif mat.Type == 'PLSTRN':
            Type = 3
        else:
            raise RuntimeError(f'Invalid Type:  Type={mat.Type}; must be 1=NLELAST '
                               '2=PLASTIC or 3=PLSTRN')
        data = [mid,
                # not sure
                0 if mat.tid is None else mat.tid,
                Type,
                0.0 if mat.h is None else mat.h,
                0 if mat.yf is None else mat.yf,
                0 if mat.hr is None else mat.hr,
                0.0 if mat.limit1 is None else mat.limit1,
                0.0 if mat.limit2 is None else mat.limit2,
                a, bmat, c]
        assert None not in data, f'MATS1 {data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_matt1(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATT1"""
    key = (703, 7, 91)
    nfields = 12
    spack = Struct(endian + b'12i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT1[mid]
        #(mid, E_table, G_table, nu_table, rho_table, A_table, dunno_a, ge_table,
         #st_table, sc_table, ss_table, dunno_b) = data

        e_table = mat.e_table if mat.e_table is not None else 0
        g_table = mat.g_table if mat.g_table is not None else 0
        nu_table = mat.nu_table if mat.nu_table is not None else 0
        rho_table = mat.rho_table if mat.rho_table is not None else 0
        a_table = mat.a_table if mat.a_table is not None else 0
        ge_table = mat.ge_table if mat.ge_table is not None else 0
        st_table = mat.st_table if mat.st_table is not None else 0
        sc_table = mat.sc_table if mat.sc_table is not None else 0
        ss_table = mat.ss_table if mat.ss_table is not None else 0

        e_table = -e_table + 100000000 if e_table < 0 else e_table
        g_table = -g_table + 100000000 if g_table < 0 else g_table
        nu_table = -nu_table + 100000000 if nu_table < 0 else nu_table
        rho_table = -rho_table + 100000000 if rho_table < 0 else rho_table
        a_table = -a_table + 100000000 if a_table < 0 else a_table
        ge_table = -ge_table + 100000000 if ge_table < 0 else ge_table
        st_table = -st_table + 100000000 if st_table < 0 else st_table
        sc_table = -sc_table + 100000000 if sc_table < 0 else sc_table
        ss_table = -ss_table + 100000000 if ss_table < 0 else ss_table

        data = [mid, e_table, g_table, nu_table, rho_table,
                a_table, 0, ge_table, st_table, sc_table,
                ss_table, 0]
        assert None not in data, f'MATT1 data={data}'
        assert min(data) == 0, f'MATT1 data={data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_matt2(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATT2"""
    #Record - MATT2(803,8,102)
    #Word Name Type Description
    #1 MID I Material identification number
    #2 TID(15) I TABLEMi entry identification numbers
    #17 UNDEF None
    key = (803, 8, 102)
    nfields = 17
    spack = Struct(endian + b'17i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT2[mid]
        data = [
            mat.mid,
            0 if mat.g11_table is None else mat.g11_table,
            0 if mat.g12_table is None else mat.g12_table,
            0 if mat.g13_table is None else mat.g13_table,

            0 if mat.g22_table is None else mat.g22_table,
            0 if mat.g23_table is None else mat.g23_table,
            0 if mat.g33_table is None else mat.g33_table,
            0 if mat.rho_table is None else mat.rho_table,

            0 if mat.a1_table is None else mat.a1_table,
            0 if mat.a2_table is None else mat.a2_table,
            0 if mat.a3_table is None else mat.a3_table,
            0,

            0 if mat.ge_table is None else mat.ge_table,
            0 if mat.st_table is None else mat.st_table,
            0 if mat.sc_table is None else mat.sc_table,
            0 if mat.ss_table is None else mat.ss_table,
            0,
        ]
        assert None not in data, f'MATT2 data={data}'
        assert min(data) == 0, f'MATT2 data={data}'
        #print('MATT2', data, len(data))
        assert len(data) == nfields, len(data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_matt4(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATT4"""
    key = (2303, 23, 237)
    nfields = 7
    spack = Struct(endian + b'7i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT4[mid]
        k_table = mat.k_table if mat.k_table is not None else 0
        cp_table = mat.cp_table if mat.cp_table is not None else 0
        h_table = mat.h_table if mat.h_table is not None else 0
        mu_table = mat.mu_table if mat.mu_table is not None else 0
        hgen_table = mat.hgen_table if mat.hgen_table is not None else 0

        #(mid, tk, tcp, null, th, tmu, thgen) = out
        data = [mid, k_table, cp_table, 0, h_table, mu_table, hgen_table]
        assert None not in data, f'MATT4 data={data}'
        assert min(data) >= 0, f'MATT4 data={data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_matt5(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATT5"""
    key = (2403, 24, 238)
    nfields = 10
    spack = Struct(endian + b'10i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT5[mid]
        kxx_table = mat.kxx_table if mat.kxx_table is not None else 0
        kxy_table = mat.kxy_table if mat.kxy_table is not None else 0
        kxz_table = mat.kxz_table if mat.kxz_table is not None else 0
        kyy_table = mat.kyy_table if mat.kyy_table is not None else 0
        kyz_table = mat.kyz_table if mat.kyz_table is not None else 0
        kzz_table = mat.kzz_table if mat.kzz_table is not None else 0

        cp_table = mat.cp_table if mat.cp_table is not None else 0
        #h_table = mat.h_table if mat.h_table is not None else 0
        #mu_table = mat.mu_table if mat.mu_table is not None else 0
        hgen_table = mat.hgen_table if mat.hgen_table is not None else 0


        #(mid, tk1, tk2, tk3, tk4, tk5, tk6, tcp, null, thgen) = out
        data = [mid, kxx_table, kxy_table, kxz_table, kyy_table, kyz_table, kzz_table,
                cp_table, 0, hgen_table]
        assert None not in data, f'MATT5 data={data}'
        assert min(data) == 0, f'MATT5 data={data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes

def _write_matt8(model, name, mids, nmaterials, op2, op2_ascii, endian):
    """writes the MATT8"""
    key = (903, 9, 336)
    nfields = 19
    # nx
    # Record - MATT8(903,9,336)
    # Word Name Type Description
    # 1 MID I
    # 2 TID(9) I TABLEMi entry identification numbers
    # 11 UNDEF None
    # 12 TID(7) I TABLEMi entry identification numbers
    # 19 UNDEF None
    spack = Struct(endian + b'19i')
    nbytes = write_header(name, nfields, nmaterials, key, op2, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT8[mid]
        #mat.uncross_reference()
        #print(mat.get_stats())
        data = [
            mat.mid,
            mat.e1_table, mat.e2_table, mat.nu12_table,
            mat.g12_table, mat.g1z_table, mat.g2z_table,
            mat.rho_table,
            mat.a1_table, mat.a2_table,
            0,
            mat.xt_table, mat.xc_table,
            mat.yt_table, mat.yc_table,
            mat.s_table,
            mat.ge_table,
            mat.f12_table,
            0,
        ]
        assert None not in data, f'MATT8 data={data}'
        assert min(data) == 0, f'MATT8 data={data}'
        #k_table = mat.k_table if mat.k_table is not None else 0
        #cp_table = mat.cp_table if mat.cp_table is not None else 0
        #h_table = mat.h_table if mat.h_table is not None else 0
        #mu_table = mat.mu_table if mat.mu_table is not None else 0
        #hgen_table = mat.hgen_table if mat.hgen_table is not None else 0

        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2.write(spack.pack(*data))
    return nbytes
