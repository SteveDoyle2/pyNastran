from __future__ import absolute_import
from struct import pack, Struct

from .geom1 import write_geom_header, close_geom_table

def write_mpt(op2, op2_ascii, obj, endian=b'<'):
    if not hasattr(obj, 'materials'):
        return
    nmaterials = len(obj.materials)
    if nmaterials == 0:
        return
    write_geom_header(b'MPT', op2, op2_ascii, endian=endian)
    itable = -3

    mtypes = [
        'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9'
    ]
    #out = obj.get_card_ids_by_card_types(mtypes)
    from collections import defaultdict
    out = defaultdict(list)
    for mid, mat in obj.materials.items():
        if mat.type in mtypes:
            out[mat.type].append(mat.mid)

    for name, mids in out.items():
        nmaterials = len(mids)
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
                data = [mid, mat.G11, mat.G12, mat.G13, mat.G22, mat.G23, mat.G33,
                        mat.rho, mat.a1, mat.a2, mat.a3, mat.tref,
                        mat.ge, mat.St, mat.Sc, mat.Ss, mat.mcsid]

                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT3':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                #(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
                 #blank, ax, ath, az, tref, ge, blank) = out

                data = [
                    mid,
                    mat.ex, mat.eth, mat.ez, mat.nuxth, mat.nuthz, mat.nuzx,
                    mat.rho, mat.gzx, 0, mat.ax, mat.ath, mat.az, mat.tref,
                    mat.ge, 0]

                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))
        elif name == 'MAT4':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                #(mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
                print(mat)
                data = [
                    mid,
                    mat.k, mat.cp, mat.rho, mat.H, mat.mu, mat.hgen,
                    mat.ref_enthalpy, mat.tch, mat.tdelta, mat.qlat,
                ]
                #print(data)
                op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
                op2.write(spack.pack(*data))

        elif name == 'MAT8':
            for mid in sorted(mids):
                mat = obj.materials[mid]
                data = [mid, mat.e11, mat.e22, mat.nu12, mat.g12, mat.g1z, mat.g2z,
                        mat.rho, mat.a1, mat.a2, mat.tref,
                        mat.Xt, mat.Xc, mat.Yt, mat.Yc, mat.S,
                        mat.ge, mat.F12, mat.strn]

                print('MAT8', data)
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

                print('MAT9', data, len(data))
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

