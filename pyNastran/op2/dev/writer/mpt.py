from collections import defaultdict
from struct import pack, Struct

from .geom1 import write_geom_header, close_geom_table

def write_mpt(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'materials'):
        return
    nmaterials = len(obj.materials) + len(obj.thermal_materials)
    mtypes = [
        'MATS1', 'MATS3', 'MATS8',
        'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
    ]
    out = defaultdict(list)
    for mtype in mtypes:
        mat_dict = getattr(obj, mtype)
        if len(mat_dict) == 0:
            continue
        out[mtype] = list(mat_dict.keys())
    if len(obj.nlparms):
        out['NLPARM'] = obj.nlparms.keys()

    if nmaterials == 0 and len(out) == 0:
        return
    write_geom_header(b'MPT', op2, op2_ascii, endian=endian)
    itable = -3

    materials_to_skip = ['NLPARM']
    #mtypes = [
        #'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9',
    #]
    #out = obj.get_card_ids_by_card_types(mtypes)
    for mid, mat in obj.materials.items():
        #if mat.type in mtypes:
        out[mat.type].append(mat.mid)
    for mid, mat in obj.thermal_materials.items():
        #if mat.type in mtypes:
        out[mat.type].append(mat.mid)


    for name, mids in sorted(out.items()):
        #print('MPT', name, mids)
        nmaterials = len(mids)
        if name in materials_to_skip:
            obj.log.warning('skipping MPT-%s' % name)
            continue

        #if nmaterials == 0:
            #continue

        if name == 'MAT1':
            key = (103, 1, 77)
            nfields = 12
            spack = Struct(endian + b'i10fi')
            # mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
        elif name == 'MAT2':
            key = (203, 2, 78)
            nfields = 17
            spack = Struct(endian + b'i15fi')
            #(mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
             #tref, ge, St, Sc, Ss, mcsid) = out
        elif name == 'MAT3':
            key = (1403, 14, 122)
            nfields = 16
            spack = Struct(endian + b'i8fi5fi')
        elif name == 'MAT4':
            key = (2103, 21, 234)
            nfields = 11
            spack = Struct(endian + b'i10f')
        elif name == 'MAT5':
            key = (2203, 22, 235)
            nfields = 10
            spack = Struct(endian + b'i9f')
        elif name == 'MAT8':
            key = (2503, 25, 288)
            nfields = 19
            spack = Struct(endian + b'i18f')
            #(mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
            # tref, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
        elif name == 'MAT9':
            key = (2603, 26, 300)
            nfields = 35
            spack = Struct(endian + b'i 30f 4i')
        elif name == 'MAT10':
            key = (2801, 28, 365)
            nfields = 5
            spack = Struct(endian + b'i4f')
        elif name == 'MATS1':
            key = (503, 5, 90)
            nfields = 11
            spack = Struct(endian + b'3ifiiff3i')
        elif name == 'MATT1':
            key = (703, 7, 91)
            nfields = 12
            spack = Struct(endian + b'12i')
        elif name == 'MATT4':
            key = (2303, 23, 237)
            nfields = 7
            spack = Struct(endian + b'7i')
        elif name == 'MATT5':
            key = (2403, 24, 238)
            nfields = 10
            spack = Struct(endian + b'10i')
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
                #mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
                data = [mid, mat.e, mat.g, mat.nu, mat.rho, mat.a, mat.tref,
                        mat.ge, mat.St, mat.Sc, mat.Ss, mat.mcsid]

                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT2':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                #print(mat.get_stats())
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
                op2.write(spack.pack(*data))
        elif name == 'MAT3':
            for mid in sorted(mids):
                mat = obj.materials[mid]
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
                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT4':
            for mid in sorted(mids):
                mat = obj.thermal_materials[mid]
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
        elif name == 'MAT5':
            for mid in sorted(mids):
                mat = obj.thermal_materials[mid]
                #(mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
                data = [
                    mid,
                    mat.kxx, mat.kxy, mat.kxz, mat.kyy, mat.kyz, mat.kzz, # 6
                    mat.cp, mat.rho, mat.hgen,
                ]
                #print('MAT5 -', data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))

        elif name == 'MAT8':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                data = [mid, mat.e11, mat.e22, mat.nu12, mat.g12, mat.g1z, mat.g2z,
                        mat.rho, mat.a1, mat.a2, mat.tref,
                        mat.Xt, mat.Xc, mat.Yt, mat.Yc, mat.S,
                        mat.ge, mat.F12, mat.strn]

                #print('MAT8 -', data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT9':
            for mid in sorted(mids):
                mat = obj.materials[mid]
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
                    mat.rho] + mat.A + [mat.tref, mat.ge,
                    0, 0, 0, 0,
                ]
                #assert len(mat.A) == 6, mat.A
                assert len(data) == nfields

                #print('MAT9 -', data, len(data))
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT10':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                #(mid, bulk, rho, c, ge) = out
                data = [mid, mat.bulk, mat.rho, mat.c, mat.ge]
                assert len(data) == nfields

                #print('MAT10 -', data, len(data))
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MATS1':
            for mid in sorted(mids):
                mat = obj.MATS1[mid]
                #(mid, tid, Type, h, yf, hr, limit1, limit2, a, bmat, c) = out
                a = 0
                bmat = 0
                c = 0

                if mat.Type == 'NLELAST':
                    Type = 1
                elif mat.Type == 'PLASTIC':
                    Type = 2
                else:
                    raise RuntimeError('Invalid Type:  Type=%s; must be 1=NLELAST '
                                       'or 2=PLASTIC' % (mat.Type))
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
                assert None not in data, 'MATS1 %s' % data
                assert len(data) == nfields
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MATT1':
            for mid in sorted(mids):
                mat = obj.MATT1[mid]
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
                assert None not in data, 'MATT1 %s' % data
                assert len(data) == nfields
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MATT4':
            for mid in sorted(mids):
                mat = obj.MATT4[mid]
                k_table = mat.k_table if mat.k_table is not None else 0
                cp_table = mat.cp_table if mat.cp_table is not None else 0
                h_table = mat.h_table if mat.h_table is not None else 0
                mu_table = mat.mu_table if mat.mu_table is not None else 0
                hgen_table = mat.hgen_table if mat.hgen_table is not None else 0

                #(mid, tk, tcp, null, th, tmu, thgen) = out
                data = [mid, k_table, cp_table, 0, h_table, mu_table, hgen_table]
                assert None not in data, 'MATT4 %s' % data
                assert len(data) == nfields
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MATT5':
            for mid in sorted(mids):
                mat = obj.MATT5[mid]
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
                data = [mid, kxx_table, kxy_table, kxz_table, kyy_table, kyy_table, kzz_table,
                        cp_table, 0, hgen_table]
                assert None not in data, 'MATT4 %s' % data
                assert len(data) == nfields
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
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
    close_geom_table(op2, op2_ascii, itable)
    #-------------------------------------
