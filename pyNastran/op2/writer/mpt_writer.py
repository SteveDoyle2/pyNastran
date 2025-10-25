"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import BinaryIO, TYPE_CHECKING

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header
from pyNastran.bdf.cards.dynamic import (
    TSTEPNL, NLPARM, NLPARM_CONV_MAP, NLPARM_KMETHOD_MAP, NLPARM_INT_OUT_MAP,
) # TSTEP
CONV_NLPARM_MAP = {value: key for key, value in NLPARM_CONV_MAP.items()}
KMETHOD_NLPARM_MAP = {value: key for key, value in NLPARM_KMETHOD_MAP.items()}
INT_OUT_NLPARM_MAP = {value: key for key, value in NLPARM_INT_OUT_MAP.items()}
if TYPE_CHECKING:
    from pyNastran.op2.op2_geom import BDF, OP2Geom

def write_mpt(op2_file, op2_ascii, model, endian=b'<',
              nastran_format: str='nx'):
    """writes the MPT/MPTS table"""
    if not hasattr(model, 'materials'):
        return
    nmaterials = len(model.materials) + len(model.thermal_materials)
    #print(f'nmaterials={nmaterials}')

    # the code will crash if a something is in the out dict, but
    # not handled properly
    materials_to_skip = [
        'MATS3', 'MATS8',
        'MATT3', 'MATT9',
        #  other
        'NLPCI', 'MAT3D',
        'MATPOR',
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

    card_name_attrs = [
        ('NLPARM', model.nlparms),
        ('NLPCI', model.nlpcis),

        # strangely TSTEPs are stored in the DYNAMIC table, but TSTEPNLs are
        # stored in the MPT table
        #'tstepnls' : ['TSTEPNL', 'TSTEP1'],
        ('TSTEPNL', model.tstepnls),
    ]
    for card_name, data_dict in card_name_attrs:
        keysi = list(data_dict.keys())
        out[card_name] = keysi

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

    if nmaterials == 0 and len(out) == 0:
        return
    write_geom_header(b'MPT', op2_file, op2_ascii, endian=endian)
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
        func = MATERIAL_MAP[name]
        nbytes = func(model, name, mids, nmaterials,
                      op2_file, op2_ascii, endian)

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

def write_tstepnl(model: BDF, name: str, tstepnl_ids, nmaterials,
                  op2_file, op2_ascii, endian):
    """
    TSTEPNL(3103,31,337)

    NX 2019.2
    23 KDAMP    I Flags to include differential stiffness to form structural damping
    24 KUPDATE  I Method for dynamic matrix update
    25 KUSTEP   I Number of iterations before the stiffness update
    26 TINTOPT  I Time integration method
    27 GAMMA   RS Amplitude decay factor for 2nd order transient integration

    """
    key = (3103, 31, 337)
    nfields = 27
    spack = Struct(endian + b'iif5i3f3if3i4f 4if')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for tstepnl_id in sorted(tstepnl_ids):
        print('tstepnl_id =', tstepnl_id)
        tstep = model.tstepnls[tstepnl_id]
        #(sid, ndt, dt, no, method_int, kstep, max_iter, conv_int, eps_u, eps_p, eps_w,
        #     max_div, max_qn, max_ls, fstress, max_bisect,
         #    adjust, mstep, rb, max_r, utol, rtol_b,
        #     kdamp, kupdate, kustep, tintopt, gamma)
        method = tstep.method
        if method == 'AUTO':
            method_int = 1
        elif method == 'TSTEP':
            method_int = 2
        elif method == 'ADAPT':
            method_int = 3
        else:  # pragma: no cover
            raise NotImplementedError('tstepnl=%s method=%r' % (tstep.sid, method))

        conv = tstep.conv
        if conv == 'W':
            conv_int = 1
        elif conv == 'P':  # guess based on format
            conv_int = 2
        elif conv == 'PW':
            conv_int = 3
        elif conv == 'U':
            conv_int = 4
        elif conv == 'UPW':
            conv_int = 7
        #elif conv_int == 3:
            #conv = 'ADAPT'
        else:  # pragma: no cover
            raise NotImplementedError('tstepnl=%s conv=%r' % (tstep.sid, conv))

        #print(tstep.get_stats())
        kstep = 0 if tstep.kstep is None else tstep.kstep
        #min_iter = 0 if tstep.min_iter == None else tstep.min_iter
        mstep = 0 if tstep.mstep is None else tstep.mstep

        kdamp = 0
        kupdate = 0
        kustep = 0
        tintopt = 0
        gamma = 0.0
        data = [tstep.sid, tstep.ndt, tstep.dt, tstep.no, method_int,
                kstep, tstep.max_iter, conv_int,
                tstep.eps_u, tstep.eps_p, tstep.eps_w,
                tstep.max_div, tstep.max_qn, tstep.max_ls,
                tstep.fstress, tstep.max_bisect,
                tstep.adjust, mstep, tstep.rb, tstep.max_r,
                tstep.utol, tstep.rtol_b,
                # nx 2019
                kdamp, kupdate, kustep,
                tintopt, gamma]
        assert None not in data, data
        op2_ascii.write('  tstepnl_id=%s data=%s\n' % (tstepnl_id, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_nlparm(model: OP2Geom, name: str,
                 mids: list[int], nmaterials: int,
                 op2_file: BinaryIO, op2_ascii, endian: bytes):
    r"""
    NLPARM(3003,30,286) - record 27

    NX 2019.2
    Word Name Type Description
    1 SID       I Set identification number
    2 NINC      I Number of increments
    3 DT       RS Incremental time interval for creep analysis
    4 KMETHOD   I Method for controlling stiffness updates
    5 KSTEP     I Number of iterations before the stiffness update
    6 MAXITER   I Limit on number of iterations for each load increment
    7 CONV      I Flags to select convergence criteria
    8 INTOUT    I Intermediate output flag
    9 EPSU     RS Error tolerance for displacement U criterion
    10 EPSP    RS Error tolerance for displacement P criterion
    11 EPSW    RS Error tolerance for displacement W criterion
    12 MAXDIV   I Limit on probable divergence conditions
    13 MAXQN    I Maximum number of quasi-Newton correction vectors
    14 MAXLS    I Maximum number of line searches
    15 FSTRESS RS Fraction of effective stress
    16 LSTOL   RS Line search tolerance
    17 MAXBIS   I Maximum number of bisections
    18 MAXR    RS Maximum ratio for the adjusted arc-length increment
    19 RTOLB   RS Maximum value of incremental rotation

    ndata = 80:
              sid nic dt   km ks max con int  epu   epp   epw   mx mx  mx fstr  lso  mx mx    rtolb
    ints    = (1, 10, 0,   1, 5, 25, -1, 0,   0.01, 0.01, 0.01, 3, 25, 4, 0.20, 0.5, 5, 20.0, 20.0, 0)
    floats  = (1, 10, 0.0, 1, 5, 25, -1, 0.0, 0.01, 0.01, 0.01, 3, 25, 4, 0.20, 0.5, 5, 20.0, 20.0, 0.0)

    # C:\MSC.Software\msc_nastran_runs\lcdf07a.op2
    ints    = (10000001, 1, 0,   1, 500, 25, 14, 0,   0.01, 0.01, 0.01, 5, 25, 0,   0.2, 0.5, 5, 20.0, 20.0, 0)
    floats  = (10000001, 1, 0.0, 1, 500, 25, 14, 0.0, 0.01, 0.01, 0.01, 5, 25, 0.0, 0.2, 0.5, 5, 20.0, 20.0, 0.0)

    """
    key = (3003, 30, 286)
    nfields = 19
    #s = Struct(mapfmt(op2._endian + b'iif5i3f3iffiff', self.size))
    spack = Struct(endian + b'iif5i3f3iffiff')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for nlparm_id in sorted(mids):
        nlparm = model.nlparms[nlparm_id]
        ninc = 0 if nlparm.ninc is None else nlparm.ninc

        conv_int = CONV_NLPARM_MAP[nlparm.conv]
        kmethod_int = KMETHOD_NLPARM_MAP[nlparm.kmethod]
        int_out_int = INT_OUT_NLPARM_MAP[nlparm.int_out]
        data = [nlparm_id, ninc, nlparm.dt, kmethod_int, nlparm.kstep,
                nlparm.max_iter, conv_int, int_out_int,
                nlparm.eps_u, nlparm.eps_p,
                nlparm.eps_w, nlparm.max_div, nlparm.max_qn, nlparm.max_ls,
                nlparm.fstress, nlparm.ls_tol, nlparm.max_bisect, nlparm.max_r,
                nlparm.rtol_b]
        assert None not in data, data
        op2_ascii.write(f'  NLPARM {nlparm_id}=%s data={data[1:]}\n')
        op2_file.write(spack.pack(*data))
    return nbytes


def write_mat1(model: OP2Geom, name: str,
               mids: list[int], nmaterials: int,
               op2_file: BinaryIO, op2_ascii, endian: bytes):
    """writes the MAT1"""
    key = (103, 1, 77)
    nfields = 12
    spack = Struct(endian + b'i10fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid
        data = [mid, mat.e, mat.g, mat.nu, mat.rho, mat.a, mat.tref,
                mat.ge, mat.St, mat.Sc, mat.Ss, mat.mcsid]
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat2(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT2"""
    key = (203, 2, 78)
    nfields = 17
    spack = Struct(endian + b'i15fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat3(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT3"""
    key = (1403, 14, 122)
    nfields = 16
    spack = Struct(endian + b'i8fi5fi')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat4(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT4"""
    key = (2103, 21, 234)
    nfields = 11
    spack = Struct(endian + b'i10f')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat5(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT5"""
    key = (2203, 22, 235)
    nfields = 10
    spack = Struct(endian + b'i9f')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes


def write_mat8(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT8"""
    key = (2503, 25, 288)
    nfields = 19
    spack = Struct(endian + b'i18f')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes


def write_mat9(model: BDF, name, mids, nmaterials,
               op2_file, op2_ascii, endian):
    """writes the MAT9"""
    key = (2603, 26, 300)
    nfields = 35
    spack = Struct(endian + b'i 30f 4i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat10(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MAT10"""
    key = (2801, 28, 365)
    nfields = 5
    spack = Struct(endian + b'i4f')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, bulk, rho, c, ge) = out
        data = [mid, mat.bulk, mat.rho, mat.c, mat.ge]
        assert len(data) == nfields

        #print('MAT10 -', data, len(data))
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mat11(model: BDF, name: str, mids: list[int], nmaterials: int,
                op2_file, op2_ascii, endian):
    """writes the MAT11"""
    key = (2903,29,371)
    nfields = 32
    spack = Struct(endian + b'i 15f 16i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.materials[mid]
        #(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
             #rho, a1, a2, a3, tref, ge, ???) = out
        data = [mid, mat.e1, mat.e2, mat.e3,
                mat.nu12, mat.nu13, mat.nu23,
                mat.g12, mat.g13, mat.g23,
                mat.rho, mat.a1, mat.a2, mat.a3, mat.tref, mat.ge] + [0] * 16
        assert len(data) == nfields, (len(data), nfields)

        #print('MAT10 -', data, len(data))
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_mats1(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MATS1"""

    # check if strmeas entry should be written
    # strategy: check random MATS1 card. if STRMEAS field is blank, do not write it for backwards compatibility
    test_mat = model.MATS1[mids[0]]
    write_strmeas = False if test_mat.strmeas is None else True

    strmeas_map = {
        None : 0,  # NULL
        'UNDEF' : 1,
        'ENG' : 2,
        'TRUE': 3,
        'CAUCHY': 4,
    }

    key = (503, 5, 90)
    nfields = 12 if write_strmeas else 11
    spack = Struct(endian + b'3ifiiff4i') if write_strmeas else Struct(endian + b'3ifiiff3i')

    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATS1[mid]
        #(mid, tid, Type, h, yf, hr, limit1, limit2, a, bmat, c) = out
        a = 0
        bmat = 0
        c = 0

        if mat.nl_type == 'NLELAST':
            nl_type_int = 1
        elif mat.nl_type == 'PLASTIC':
            nl_type_int = 2
        elif mat.nl_type == 'PLSTRN':
            nl_type_int = 3
        else:  # pragma: no cover
            raise RuntimeError(f'Invalid Type:  Type={mat.nl_type}; must be 1=NLELAST '
                               '2=PLASTIC or 3=PLSTRN')

        if write_strmeas:
            data = [mid,
                    # not sure
                    mat.tid,
                    nl_type_int,
                    0.0 if mat.h is None else mat.h,
                    0 if mat.yf is None else mat.yf,
                    0 if mat.hr is None else mat.hr,
                    0.0 if mat.limit1 is None else mat.limit1,
                    0.0 if mat.limit2 is None else mat.limit2,
                    strmeas_map[mat.strmeas],

                    a, bmat, c]
        else:
            data = [mid,
                    # not sure
                    mat.tid,
                    nl_type_int,
                    0.0 if mat.h is None else mat.h,
                    0 if mat.yf is None else mat.yf,
                    0 if mat.hr is None else mat.hr,
                    0.0 if mat.limit1 is None else mat.limit1,
                    0.0 if mat.limit2 is None else mat.limit2,

                    a, bmat, c]

        assert None not in data, f'MATS1 {data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_matt1(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MATT1"""
    key = (703, 7, 91)
    nfields = 12
    spack = Struct(endian + b'12i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT1[mid]
        #(mid, E_table, G_table, nu_table, rho_table, A_table, dunno_a, ge_table,
         #st_table, sc_table, ss_table, dunno_b) = data

        e_table = mat.e_table
        g_table = mat.g_table
        nu_table = mat.nu_table
        rho_table = mat.rho_table
        a_table = mat.a_table
        ge_table = mat.ge_table
        st_table = mat.st_table
        sc_table = mat.sc_table
        ss_table = mat.ss_table

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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_matt2(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MATT2"""
    #Record - MATT2(803,8,102)
    #Word Name Type Description
    #1 MID I Material identification number
    #2 TID(15) I TABLEMi entry identification numbers
    #17 UNDEF None
    key = (803, 8, 102)
    nfields = 17
    spack = Struct(endian + b'17i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT2[mid]
        data = [
            mat.mid,
            mat.g11_table, mat.g12_table, mat.g13_table,
            mat.g22_table, mat.g23_table, mat.g33_table,
            mat.rho_table,

            mat.a1_table, mat.a2_table, mat.a3_table,
            0,
            mat.ge_table, mat.st_table, mat.sc_table, mat.ss_table,
            0,
        ]
        assert None not in data, f'MATT2 data={data}'
        assert min(data) == 0, f'MATT2 data={data}'
        #print('MATT2', data, len(data))
        assert len(data) == nfields, len(data)
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_matt4(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MATT4"""
    key = (2303, 23, 237)
    nfields = 7
    spack = Struct(endian + b'7i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT4[mid]
        k_table = mat.k_table
        cp_table = mat.cp_table
        h_table = mat.h_table
        mu_table = mat.mu_table
        hgen_table = mat.hgen_table

        #(mid, tk, tcp, null, th, tmu, thgen) = out
        data = [mid, k_table, cp_table, 0, h_table, mu_table, hgen_table]
        assert None not in data, f'MATT4 data={data}'
        assert min(data) >= 0, f'MATT4 data={data}'
        assert len(data) == nfields
        op2_ascii.write('  mid=%s data=%s\n' % (mid, data[1:]))
        op2_file.write(spack.pack(*data))
    return nbytes

def write_matt5(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
    """writes the MATT5"""
    key = (2403, 24, 238)
    nfields = 10
    spack = Struct(endian + b'10i')
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT5[mid]
        kxx_table = mat.kxx_table
        kxy_table = mat.kxy_table
        kxz_table = mat.kxz_table
        kyy_table = mat.kyy_table
        kyz_table = mat.kyz_table
        kzz_table = mat.kzz_table

        cp_table = mat.cp_table
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
        op2_file.write(spack.pack(*data))
    return nbytes

def write_matt8(model: BDF, name, mids, nmaterials,
                op2_file, op2_ascii, endian):
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
    nbytes = write_header(name, nfields, nmaterials, key, op2_file, op2_ascii)
    for mid in sorted(mids):
        mat = model.MATT8[mid]
        #mat.uncross_reference()
        #print(mat.get_stats())
        data = [
            mat.mid,
            mat.e1_table, mat.e2_table, mat.nu12_table,
            mat.g12_table, mat.g1z_table, mat.g2z_table,
            mat.rho_table, mat.a1_table, mat.a2_table,
            0,
            mat.xt_table, mat.xc_table,
            mat.yt_table, mat.yc_table,
            mat.s_table, mat.ge_table, mat.f12_table,
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
        op2_file.write(spack.pack(*data))
    return nbytes


MATERIAL_MAP = {
    'MAT1': write_mat1,
    'MAT2': write_mat2,
    'MAT3': write_mat3,
    'MAT4': write_mat4,
    'MAT5': write_mat5,
    'MAT8': write_mat8,
    'MAT9': write_mat9,
    'MAT10': write_mat10,
    'MAT11': write_mat11,
    'MATS1': write_mats1,
    'MATT1': write_matt1,
    'MATT2': write_matt2,
    'MATT4': write_matt4,
    'MATT5': write_matt5,
    'MATT8': write_matt8,
    'TSTEPNL': write_tstepnl,
    'NLPARM': write_nlparm,
}
