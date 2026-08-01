from __future__ import annotations
import io
import copy
from typing import cast, TYPE_CHECKING
import numpy as np

from pyNastran.utils import PathLike
from pyNastran.utils.numpy_utils import float_types
from pyNastran.bdf.cards.aero.aero import (
    CAERO1, SPLINE1, SPLINE2, SPLINE3,)

# from .zona_cards.zaero_sets import (
#     SETADD)
from pyNastran.bdf.cards.aero.zaero_cards.atm import (
    ATMOS, FIXMATM, FIXHATM, FIXMACH, FIXMDEN)
from pyNastran.bdf.cards.aero.zaero_cards.spline import (
    SPLINE1_ZAERO, SPLINE2_ZAERO,
    #SPLINE2_ZAERO, SPLINE3_ZAERO,
)
from pyNastran.bdf.cards.aero.zaero_cards.geometry import (
    PANLST1, PANLST2, #PANLST3, SEGMESH,
    CAERO7, AESURFZ, # BODY7, PAFOIL7, PAFOIL8, AESLINK,
)
from pyNastran.bdf.cards.aero.zaero_cards.flutter import (
    # FLUTTER_ZAERO,
    MKAEROZ)
from pyNastran.bdf.cards.aero.zaero_cards.trim import (
    TRIM_ZAERO, TRIMVAR, TRIMLNK,)
from pyNastran.bdf.cards.aero.zaero_cards.manuever import (
    ACTU) # MLOADS, LOADMOD, RBRED,)
from pyNastran.bdf.cards.aero.zaero_cards.cards import (
    # MLDPRNT, MLDSTAT, MINSTAT, MLDTRIM, MLDCOMD, MLDTIME,
    AEROZ, ACOORD, #ATTACH,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, SPLINE1, AEFACT, AELIST, AEROS, FLUTTER, FLFACT


def nastran_to_zaero(bdf_filename: PathLike | BDF,
                     zaero_inp_filename: PathLike='',
                     length_unit: str='IN',
                     mass_unit: str='SLIN',) -> BDF:
    from pyNastran.bdf.bdf import BDF, read_bdf
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename)

    bdf_file = io.StringIO()
    modelz = BDF(mode='zaero')
    zaero = modelz.zaero
    zaero_add = zaero._add_methods

    # print(aeros.get_stats())
    naeroz, coords = _convert_aeros(
        model, modelz, mass_unit, length_unit)
    modelz._write_coords(bdf_file)

    ncaero7 = _convert_caeros(model, modelz)
    modelz._write_coords(bdf_file)
    ncord2r, nspline1, nspline2, nspline3 = _convert_splines(
        model, modelz)

    for aefact_id, aefact in model.aefacts.items():
        modelz.add_aefact(aefact_id, aefact.fractions*100,
                          comment=aefact.comment)

    # aesurf_names = []
    aesurf_dict, naesurfz, npanlst2, nactu = _convert_aesurf(
        model, modelz)

    trimvar_dict, nxyz_root, pqr_dot_root = _convert_aestats(model)

    ntrimvar = _make_trimvars(trimvar_dict, modelz)
    ntrim, nmkaeroz_trim = _convert_trim_mkaeroz(
        model, modelz,
        nxyz_root, pqr_dot_root,
        trimvar_dict, aesurf_dict,
    )
    modelz._write_coords(bdf_file)
    nflutter, nmkaeroz_flutter = _convert_flutter(
        model, modelz, mass_unit, length_unit)
    modelz._write_coords(bdf_file)
    card_count_dict = {
        'CORD2R': ncord2r,
        'AEROZ': naeroz,
        'FLUTTER': nflutter,
        'MKAEROZ': nmkaeroz_trim + nmkaeroz_flutter,
        'CAERO7': ncaero7,
        'TRIMVAR': ntrimvar,
        'AESURFZ': naesurfz,
        'PANLST2': npanlst2,
        'ACTU': nactu,
        'TRIM': ntrim,
        'SPLINE1': nspline1,
        'SPLINE2': nspline2,
        'SPLINE3': nspline3,
    }

    for key, ncard in card_count_dict.items():
        if ncard:
            modelz.card_count[key] = ncard

    if zaero_inp_filename != '':
        modelz.write_bdf(zaero_inp_filename)
    return modelz


def _convert_flutter(model: BDF, modelz: BDF,
                     mass_unit: str, length_unit: str,
                     ) -> tuple[int, int]:
    print_flag = 0
    flutter_id = 0

    flt_id = 0
    mkaerozs = []
    add_methods = modelz.zaero._add_methods
    if len(model.flutters) == 0:
        mkaeroz_id = 99
        for mkaero in model.mkaeros:
            # machs: array([0.8])
            # reduced_freqs: array([0.1, 0.215, 0.464, 1.])
            mach = mkaero.machs[0]
            freqs = mkaero.reduced_freqs
            filename = f'MK_{mach:g}.out'
            mkaeroz = MKAEROZ(
                mkaeroz_id, mach, flt_id, filename, print_flag,
                freqs, method=0, save='SAVE')
            add_methods.add_mkaeroz_object(mkaeroz)
            break

    for sid, flutter in model.flutters.items():
        mkaeroz_id = sid + 100
        mkaerozs.append(mkaeroz_id)
        #   density : 76
        #   epsilon : 0.001
        #   headers : ['density', 'mach', 'velocity']
        #   imethod : 'L'
        #   mach   : 77
        #   method : 'PKNL'
        #   nvalue : None
        #   omax   : None
        #   reduced_freq_velocity : 78
        mach_ref = flutter.mach_ref
        density_ref = flutter.density_ref
        is_mach_constant = np.allclose(mach_ref.max(), mach_ref.min())
        is_density_constant = np.allclose(density_ref.max(), density_ref.min())
        method = flutter.method

        machs = mach_ref.factors
        mach = machs[0]
        if model.mkaeros:
            for makero in mkaeros:
                freqs = mkaero.reduced_freqs.tolist()
                break
        else:
            freqs = [1.0]
        filename = f'MK{mkaeroz_id}_{mach:g}.out'
        mkaeroz = MKAEROZ(
            mkaeroz_id, mach, flt_id, filename, print_flag,
            freqs, method=0, save='SAVE')
        add_methods.add_mkaeroz_object(mkaeroz)

        if method in {'PKNL', 'PKNLS'}:
            velocity_ref = flutter.reduced_freq_velocity_ref
            if not is_mach_constant:
                model.log.warning(f'assuming mach is constant...')
            if is_mach_constant or 1:
                velocity = velocity_ref.factors.tolist()
                rho = density_ref.factors.tolist()
                atm = FIXMACH(sid, mkaeroz_id, mass_unit, length_unit,
                              flutter_id, print_flag, velocity, rho)
            else:  # pragma: no cover
                print(flutter.get_stats())
                raise RuntimeError((method, is_mach_constant, is_density_constant))
        elif method == 'KE':
            model.log.warning(f'skipping flutter card because method={method!r}...')
            continue
        else:  # pragma: no cover
            print(flutter.get_stats())
            raise RuntimeError((method, is_mach_constant, is_density_constant))
        add_methods.add_flutter_table_object(atm)
        # modelz.zaero.flutter_table[sid] = atm

    nflutter = len(modelz.zaero.flutter_table)
    nmkaeroz2 = len(modelz.zaero.mkaeroz)
    return nflutter, nmkaeroz2


def _convert_trim_mkaeroz(
        model: BDF,
        modelz: BDF,
        nxyz_root: list[str],
        pqr_dot_root: list[str],
        trimvar_dict: dict[str, int],
        aesurf_dict: dict[int, str],
        ) -> tuple[int, int]:
    zaero = modelz.zaero
    zaero_add = zaero._add_methods

    nmkaeroz = 0
    ntrim = 0

    mass = 1.
    mass_inertia = np.ones(6)
    if len(model.trims):
        dcg = [0.1, 0.2, 0.3]
        from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
        mass, cg, mass_inertia = mass_properties(model)

    # weight = 1000.
    # inertia = [0.4, 0.5, 0.6,
    #            0.7, 0.8, 0.9]

    true_g = 'TRUE'
    wtmass = model.wtmass
    weight = mass * wtmass
    inertia = mass_inertia * wtmass

    if 'AUNITS' in model.params:
        param = model.params['AUNITS']
        # print(param.get_stats())  # check param.value is correct
        if param.values[0] != 1.0:
            true_g = 'G'

    loadset = 0
    trimobj_id = 0
    trimcon_id = 0
    for trim_id, trim in model.trims.items():
        trim.cross_reference(model)
        mkaeroz_id = trim_id + 1
        flt_id = trim_id + 2

        nxyz = copy.deepcopy(nxyz_root)
        pqr_dot = copy.deepcopy(pqr_dot_root)
        trimvar_ids = []
        uxs = []

        print_flag = 0
        freqs = [0.]
        filename = ''
        mkaeroz = MKAEROZ(
            mkaeroz_id, trim.mach, flt_id,
             filename, print_flag, freqs,
             method=0, save=None)
        zaero.mkaeroz[mkaeroz_id] = mkaeroz
        nmkaeroz += 1

        commenti = ''
        for label, ux in zip(trim.labels, trim.uxs):
            if label in ['URDD1', 'URDD2', 'URDD3']:
                iurddt = int(label[-1]) - 1  # 1->0
                nxyz[iurddt] = ux
            elif label in ['URDD4', 'URDD5', 'URDD6']:
                iurddr = int(label[-1]) - 4  # 4->0
                pqr_dot[iurddr] = ux
            elif label in trimvar_dict:
                # no control surfaces
                aestat_id = trimvar_dict[label]
                trimvar_ids.append(aestat_id)
                # trimvar_dict[aestat_id] = label
                uxs.append(ux)
                commenti += f' {label} = {ux} (TRIMVAR={aestat_id})\n'

        lower = ''
        upper = ''
        trimlnk = ''
        dmi = None
        sym = 'SYM'
        comment = (
            f' weight={weight} (assumed)\n'
            f' dref={dcg} (assumed)\n'
            f' inertia={inertia} (assumed)\n'
            f' nxyz={nxyz} g\n'
            f' pqr_dot={pqr_dot} rad/s^2\n'
            f'{commenti}'
        )
        for label, aesurf_id in aesurf_dict.items():
            if label not in trimvar_dict:
                trimvar_ids.append(aesurf_id)
                comment += f' {label} = FREE (AESURFZ={aesurf_id})\n'
                uxs.append('FREE')
                trimvar = TRIMVAR(aesurf_id, label,
                                  lower, upper, trimlnk,
                                  dmi, sym)
                zaero_add.add_trimvar_object(trimvar)

        trimz = TRIM_ZAERO(
            trim_id, mkaeroz_id, trim.q,
            trimobj_id, trimcon_id,
            weight, dcg, inertia,
            true_g, nxyz, pqr_dot,
            loadset, trimvar_ids, uxs,
            wtmass=wtmass, comment=comment)
        modelz.trims[trim_id] = trimz
        ntrim += 1
    return ntrim, nmkaeroz


def _convert_aeros(model: BDF,
                   modelz: BDF,
                   mass_unit: str,
                   length_unit: str) -> tuple[int, list]:
    aeros: AEROS = model.aeros
    coords = {}
    if aeros is None:
        cref = 1.1
        bref = 1.2
        sref = 1.3
        acsid = 0
        rcsid = 0
    else:
        aeros.cross_reference(model)
        cref = aeros.cref
        bref = aeros.bref
        sref = aeros.sref
        acsid = aeros.acsid
        rcsid = aeros.rcsid
    aeroz = AEROZ(
        mass_unit, length_unit,
        cref, bref, sref,
        acsid=acsid, rcsid=rcsid, xyz_ref=None,
        sym_xz='NO',  # sym_xy='NO',
    )
    naeroz = 1
    modelz.aeros = aeroz
    if aeros is None:
        return naeroz, coords

    coords[acsid] = aeros.acsid_ref
    coords[rcsid] = aeros.rcsid_ref
    modelz.coords = coords
    bdf_file = io.StringIO()
    modelz._write_coords(bdf_file)
    # (aeroz.acsid, aeroz.bref, aeroz.sref, aeroz.sym_xy, aeroz.sym_xz))
    return naeroz, coords


def _make_trimvars(trimvar_dict: dict[str, int],
                   modelz: BDF) -> int:
    zaero = modelz.zaero
    zaero_add = zaero._add_methods
    lower = ''
    upper = ''
    trimlnk = ''
    dmi = None
    sym = 'SYM'

    ntrimvar = 0
    # if len(trimvar_dict):
    #     modelz.log.debug(f'trimvar_dict = {trimvar_dict}')
    for label, aestat_id in trimvar_dict.items():
        trimvar = TRIMVAR(aestat_id, label,
                          lower, upper, trimlnk,
                          dmi, sym)
        zaero_add.add_trimvar_object(trimvar)
        ntrimvar += 1
    return ntrimvar


def _convert_aestats(model: BDF) -> tuple[
                                      dict[str, int],
                                      list[str],
                                      list[str], ]:
    trimvar_dict = {}
    nxyz_root = ['NONE', 'NONE', 'NONE']
    pqr_dot_root = ['NONE', 'NONE', 'NONE']
    for aestat_id, aestat in model.aestats.items():
        label = aestat.label
        if label in ['URDD1', 'URDD2', 'URDD3']:
            iurdd = int(label[-1]) - 1  # 1->0
            nxyz_root[iurdd] = 'FREE'
        elif label in ['URDD4', 'URDD5', 'URDD6']:
            iurdd = int(label[-1]) - 4  # 4->0
            pqr_dot_root[iurdd] = 'FREE'
        else:
            # ALPHA, BETA
            trimvar_dict[label] = aestat_id

    if len(trimvar_dict):
        model.log.debug(f'trimvar_dict = {trimvar_dict}')
    return trimvar_dict, nxyz_root, pqr_dot_root


def _convert_aesurf(model: BDF,
                    modelz: BDF) -> tuple[dict, int, int, int]:
    zaero = modelz.zaero
    zaero_add = zaero._add_methods
    naesurfz = 0
    npanlst2 = 0
    nactu = 0

    # TODO: only does cid1
    surface_type = 'SYM' # assumed
    setg = 0
    aesurf_dict = {}
    for aesurf_id, aesurf in model.aesurf.items():
        aesurf.cross_reference(model)
        # aesurf_names.append('FREE')
        actu_id = aesurf_id
        panlst_id = aesurf.aelist_id1

        label = aesurf.label
        aesurf_dict[label] = aesurf_id
        cid = aesurf.cid1
        cid_ref = aesurf.cid1_ref
        cid_ref.comment = f'AESURF label={label!r}'
        modelz.coords[cid] = cid_ref
        aesurfz = AESURFZ(label, surface_type, cid, panlst_id, setg, actu_id)
        modelz.aesurf[label] = aesurfz
        naesurfz += 1

        # bag of panels
        aelist: AELIST = aesurf.aelist_id1_ref
        boxes = aelist.elements
        macro_id = panlst_id
        panlst = PANLST2(panlst_id, macro_id, boxes)
        zaero_add.add_panlst_object(panlst)
        npanlst2 += 1

        # actuator
        actu = ACTU(actu_id, None, None, None)
        zaero_add.add_actu_object(actu)
        nactu += 1
    return aesurf_dict, naesurfz, npanlst2, nactu


def _get_panlst_id_from_caero(caero: CAERO7) -> int:
    assert caero.type in {'CAERO1', 'CAERO7'}, str(caero)
    panlst_id = caero.eid
    return panlst_id


def _convert_caeros(model: BDF,
                    modelz: BDF) -> int:
    icaero = 11
    ncaero7 = 0

    # just map all the coords
    # overwrite if necessary
    #modelz.nodes = model.nodes
    for cid, coord in model.coords.items():
        modelz.coords[cid] = copy.deepcopy(coord)

    xref_obj = model.xref_obj
    xref_obj.cross_reference_nodes()
    xref_obj.cross_reference_coordinates()

    log = model.log
    if len(model.caeros):
        assert len(model.paeros), model.paeros

    for caero_id, caero in model.caeros.items():
        assert caero.type == 'CAERO1', caero
        caero.cross_reference(model)

        label = f'caero{icaero}'
        #label = ''
        comment = f'{caero.comment}\n{str(caero)}'

        # Nastran: xyz location frame
        # ZAero: spline plane
        cp_ref = caero.cp_ref
        assert cp_ref is not None, cp_ref
        #if cp_ref.cid != 0:
        #    log.warning('p1, p2 must be in cid=0')
        if cp_ref.cid > 0:
            # modelz.zaero.acoord
            cid = cp_ref.cid
            origin = cp_ref.origin
            dpitch = np.arccos(np.dot(cp_ref.i, [1., 0., 0.]))
            droll = np.arccos(np.dot(cp_ref.j, [0., 1., 0.]))
            delta_pitch = np.degrees(dpitch).round(3)
            theta_roll = np.degrees(droll).round(3)
            assert isinstance(delta_pitch, float_types), delta_pitch
            assert isinstance(theta_roll, float_types), theta_roll
            acoord = ACOORD(
                cid, origin, delta=delta_pitch, theta=theta_roll,
                comment=comment)
            #modelz.coords[caero.cp] = cp_ref
            modelz.coords[caero.cp] = acoord

        # CAERO1 is element-based, while CAERO7 is node-based
        nspan = 0 if caero.nspan == 0 else caero.nspan + 1
        nchord = 0 if caero.nchord == 0 else caero.nchord + 1
        assert nchord > 0, caero.get_stats()
        caero7 = CAERO7(caero.eid, label,
                        caero.p1, caero.x12, caero.p4, caero.x43,
                        caero.cp, nspan, nchord,
                        caero.lspan,
                        comment=comment)
        if caero.lspan:
            assert caero.lspan_ref.type == 'AEFACT', caero.lspan_ref
        if caero.lchord:
            log.warning('lchord is not suppported')
        panlst_id = _get_panlst_id_from_caero(caero)
        modelz.caeros[caero_id] = caero7
        ncaero7 += 1
        icaero += 1
    return ncaero7


def _convert_splines(model: BDF,
                     modelz: BDF) -> tuple[int, int, int, int]:
    cp = 950
    ncord2r = 0
    nspline1 = 0
    nspline2 = 0
    nspline3 = 0
    zaero_add = modelz.zaero._add_methods
    for spline_id, spline in model.splines.items():
        caero_id = spline.caero
        caero = modelz.caeros[caero_id]
        nboxes = caero.nboxes
        assert caero.type == 'CAERO7', caero.get_stats()
        # spline.cross_reference(model)
        # model_str = ''
        spline_type = spline.type
        if spline_type == 'SPLINE1':
            panlst_id = _get_panlst_id_from_caero(caero)

            comment = str(spline.comment)
            model_str = caero.label

            spline = cast(SPLINE1, spline)
            setg = int(spline.setg)
            set_card = model.sets[setg]
            set_ids = set_card.ids

            coord_id, ncord2r = _fit_plate_spline(model, modelz, set_ids, cp, ncord2r, spline_id)

            # cp2 = None
            modelz.sets[setg] = set_card
            splinez = SPLINE1_ZAERO(
                spline_id, panlst_id, setg, model=model_str, cp=coord_id,
                dz=spline.dz, eps=0.01, comment=comment)
            modelz.splines[spline_id] = splinez
            nspline1 += 1
            cp += 1
        elif spline_type == 'SPLINE2':
            caero_id = spline.caero
            caero = modelz.caeros[caero_id]
            panlst_id = _get_panlst_id_from_caero(caero)

            comment = str(spline.comment)
            # model_str = caero.label
            model_str = ''

            spline = cast(SPLINE2, spline)
            setg = int(spline.setg)
            set_card = model.sets[setg]
            set_ids = set_card.ids

            #coord_id, ncord2r = _fit_plate_spline(model, modelz, set_ids, cp, ncord2r, spline_id)
            coord_id = spline.cid
            modelz.sets[setg] = set_card
            splinez = SPLINE2_ZAERO(
                spline_id, panlst_id,
                setg, model=model_str, dz=spline.dz,
                eps=0.01, cp=coord_id, curvature=None,
                comment=comment)
            #print(splinez.get_stats())

            modelz.splines[spline_id] = splinez
            nspline2 += 1
            cp += 1
        else:  # pragma: no cover
            raise NotImplementedError(f'spline_type={spline_type!r} not supported')

        macro_id = caero.eid
        # panlst = PANLST1(panlst_id, macro_id,
        #                  box1=panlst_id, box2=panlst_id+nboxes)
        boxes = np.arange(panlst_id, panlst_id + nboxes)
        assert len(boxes) == nboxes
        nchord = len(caero.lchord_ref.elements) if caero.nchord == 0 else caero.nchord
        nspan = len(caero.lspan_ref.elements) if caero.nspan == 0 else caero.nspan
        # print(f'caero7={caero.eid} nspan={nspan} nchord={nchord}; nboxes={len(boxes)}')
        panlst = PANLST2(panlst_id, macro_id, boxes)
        zaero_add.add_panlst_object(panlst)
    return ncord2r, nspline1, nspline2, nspline3


def _fit_plate_spline(model: BDF, modelz: BDF,
                      set_ids: list[int],
                      cp: int, ncord2r: int,
                      spline_id: int) -> tuple[int, int]:
    spline = model.splines[spline_id]
    spline_type = spline.type
    if spline_id in model.coords:
        modelz.coords[spline_id] = model.coords[spline_id]
        coord_id = spline_id
    elif len(model.nodes):
        model.log.warning(f'No {spline_type}={spline_id} coord (cid={spline_id}) defined. Attempting to fit a plane...')

        xyz_list = [model.nodes[nid].get_position() for nid in set_ids]
        xyzs = np.vstack(xyz_list)
        origin, zaxis, xzplane = fit_plane_to_point_cloud(xyzs)
        origin = np.zeros(3)
        modelz.add_cord2r(cp, origin, origin + zaxis, origin + xzplane)
        ncord2r += 1
        coord_id = cp
    else:
        msg = f'No {spline_type}={spline_id} coord (cid={spline_id}) defined and no spline nodes were defined.'
        model.log.error(msg)
        raise RuntimeError(msg)
    return coord_id, ncord2r


def fit_plane_to_point_cloud(point_cloud: np.ndarray,
                             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fits a best-fit plane to a 3D point cloud using SVD.

    Parameters
    ----------
    point_cloud : (nnode, 3) float np.ndarray
        An array of shape (N, 3) representing N points (x, y, z).

    Returns
    -------
    centroid : (3, ) float np.ndarray
        centroid of the plane
    zaxis : (3, ) float np.ndarray
        zaxis of the plane
    xzplane : (3, ) float np.ndarray
        xzplane of the plane

    Note:
    normal_vector : (nnode, 3) float np.ndarray
        The normal_vector is a unit vector.
        The plane equation is normal[0]*x + normal[1]*y + normal[2]*z + d = 0,
        where d = -np.dot(normal, centroid).
    """
    # 1. Calculate the centroid (center of mass) of the points
    centroid = np.mean(point_cloud, axis=0)

    # 2. Translate points to the origin
    # This centers the data around (0, 0, 0)
    points_centered = point_cloud - centroid

    # 3. Apply Singular Value Decomposition (SVD)
    # The normal vector of the best-fit plane corresponds to the left singular vector
    # (or right singular vector, depending on implementation) associated with the smallest singular value.
    # In numpy's SVD (np.linalg.svd), the right singular vectors (vh) are returned.
    # The last row of vh is the normal vector.
    _, _, vh = np.linalg.svd(points_centered)

    # The normal vector is the last row of vh
    normal_vector = vh[2, :]
    orthogonal_v = find_orthogonal_vector_cross(normal_vector)
    zaxis = normal_vector
    xzplane = orthogonal_v
    return centroid, zaxis, xzplane


def find_orthogonal_vector_cross(v: np.ndarray) -> np.ndarray:
    # Convert list to numpy array if it isn't already
    v = np.array(v)

    # Define a "helper" vector that is not parallel to v
    # A common way to pick one is to find the smallest component of v,
    # set it to 1, and the other two to 0, or by some more robust logic.
    # For general robustness, we can check which standard axis is least aligned.

    if np.abs(v[0]) < np.abs(v[1]) and np.abs(v[0]) < np.abs(v[2]):
        helper = np.array([1, 0, 0])
    elif np.abs(v[1]) < np.abs(v[2]):
        helper = np.array([0, 1, 0])
    else:
        helper = np.array([0, 0, 1])

    # Calculate the cross product, which is orthogonal to both v and helper
    orthogonal_v = np.cross(v, helper)

    # Normalize the resulting vector (optional, but often useful)
    norm_orthogonal_v = np.linalg.norm(orthogonal_v)
    if norm_orthogonal_v == 0:
        # This case is extremely rare and only happens if the helper
        # vector was perfectly parallel. The logic above prevents this.
        pass
    return orthogonal_v
