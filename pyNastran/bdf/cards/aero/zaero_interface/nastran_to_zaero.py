from __future__ import annotations
import copy
from typing import cast, TYPE_CHECKING
import numpy as np

from pyNastran.utils import PathLike
from pyNastran.bdf.cards.aero.aero import (
    CAERO1, SPLINE1, SPLINE2, SPLINE3,)

# from .zona_cards.zaero_sets import (
#     SETADD)
from pyNastran.bdf.cards.aero.zaero_cards.atm import (
    ATMOS, FIXMATM, FIXHATM, FIXMACH, FIXMDEN)
from pyNastran.bdf.cards.aero.zaero_cards.spline import (
    SPLINE1_ZAERO,
    #SPLINE2_ZAERO, SPLINE3_ZAERO,
)
from pyNastran.bdf.cards.aero.zaero_cards.geometry import (
    PANLST2, #PANLST1, PANLST3, SEGMESH,
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
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF, SPLINE1, AEFACT, AELIST, AEROS


def nastran_to_zaero(bdf_filename: PathLike | BDF,
                     zaero_inp_filename: PathLike='',
                     length_unit: str='IN',
                     mass_unit: str='SLIN',) -> BDF:
    from pyNastran.bdf.bdf import BDF, read_bdf
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename)

    model2 = BDF(mode='zaero')
    zaero = model2.zaero
    zaero_add = zaero._add_methods

    # print(aeros.get_stats())
    naeroz, coords = _convert_aeros(
        model, model2, mass_unit, length_unit)

    ncaero7 = _convert_caeros(model, model2)
    ncord2r, nspline1, nspline2, nspline3 = _convert_splines(
        model, model2)

    for aefact_id, aefact in model.aefacts.items():
        model2.add_aefact(aefact_id, aefact.fractions*100,
                          comment=aefact.comment)

    # aesurf_names = []
    aesurf_dict, naesurfz, npanlst2, nactu = _convert_aesurf(
        model, model2)

    trimvar_dict, nxyz_root, pqr_dot_root = _convert_aestats(model)

    ntrimvar = _make_trimvars(trimvar_dict, model2)
    ntrim, nmkaeroz_trim = _convert_trim_mkaeroz(
        model, model2,
        nxyz_root, pqr_dot_root,
        trimvar_dict, aesurf_dict,
    )
    nflutter, nmkaeroz_flutter = _convert_flutter(
        model, model2, mass_unit, length_unit)
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
            model2.card_count[key] = ncard

    if zaero_inp_filename != '':
        model2.write_bdf(zaero_inp_filename)
    return model2


def _convert_flutter(model: BDF, model2: BDF,
                     mass_unit: str, length_unit: str,
                     ) -> tuple[int, int]:
    print_flag = 0
    flutter_id = 0

    flt_id = 0
    mkaerozs = []
    add_methods = model2.zaero._add_methods
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
        # model2.zaero.flutter_table[sid] = atm

    nflutter = len(model2.zaero.flutter_table)
    nmkaeroz2 = len(model2.zaero.mkaeroz)
    return nflutter, nmkaeroz2

def _convert_trim_mkaeroz(
        model: BDF,
        model2: BDF,
        nxyz_root: list[str],
        pqr_dot_root: list[str],
        trimvar_dict: dict[str, int],
        aesurf_dict: dict[int, str],
        ) -> tuple[int, int]:
    zaero = model2.zaero
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
            else:
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
        model2.trims[trim_id] = trimz
        ntrim += 1
    return ntrim, nmkaeroz

def _convert_aeros(model: BDF,
                   model2: BDF,
                   mass_unit: str,
                   length_unit: str) -> None:
    aeros: AEROS = model.aeros
    coords = {}
    if aeros is None:
        cref = 1.1
        bref = 1.2
        sref = 1.3
        acsid = 0
        rcsid = 0
    else:
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
    model2.aeros = aeroz
    if aeros is None:
        return naeroz, coords

    coords[acsid] = aeros.acsid_ref
    coords[rcsid] = aeros.rcsid_ref
    model2.coords = coords
    # (aeroz.acsid, aeroz.bref, aeroz.sref, aeroz.sym_xy, aeroz.sym_xz))
    return naeroz, coords


def _make_trimvars(trimvar_dict: dict[str, int],
                   model2: BDF) -> int:
    zaero = model2.zaero
    zaero_add = zaero._add_methods
    lower = ''
    upper = ''
    trimlnk = ''
    dmi = None
    sym = 'SYM'

    ntrimvar = 0
    if len(trimvar_dict):
        model2.log.debug(f'trimvar_dict = {trimvar_dict}')
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
    return trimvar_dict, nxyz_root, pqr_dot_root


def _convert_aesurf(model: BDF,
                    model2: BDF) -> tuple[dict, int, int, int]:
    zaero = model2.zaero
    zaero_add = zaero._add_methods
    naesurfz = 0
    npanlst2 = 0
    nactu = 0

    # TODO: only does cid1
    surface_type = 'SYM' # assumed
    setg = 0
    aesurf_dict = {}
    for aesurf_id, aesurf in model.aesurf.items():
        # aesurf_names.append('FREE')
        actu_id = aesurf_id
        panlst_id = aesurf.aelist_id1

        label = aesurf.label
        aesurf_dict[label] = aesurf_id
        cid = aesurf.cid1
        cid_ref = aesurf.cid1_ref
        cid_ref.comment = f'AESURF label={label!r}'
        model2.coords[cid] = cid_ref
        aesurfz = AESURFZ(label, surface_type, cid, panlst_id, setg, actu_id)
        model2.aesurf[label] = aesurfz
        naesurfz += 1

        # bag of panels
        aelist: AELIST = aesurf.aelist_id1_ref
        boxes = aelist.elements
        panlst = PANLST2(panlst_id, panlst_id, boxes)
        zaero_add.add_panlst_object(panlst)
        npanlst2 += 1

        # actuator
        actu = ACTU(actu_id, None, None, None)
        zaero_add.add_actu_object(actu)
        nactu += 1
    return aesurf_dict, naesurfz, npanlst2, nactu


def _get_panlst_id_from_caero(caero: CAERO1) -> int:
    panlst_id = caero.eid + 1
    return panlst_id


def _convert_caeros(model: BDF,
                    model2: BDF) -> int:
    icaero = 11
    ncaero7 = 0
    for caero_id, caero in model.caeros.items():
        label = f'caero{icaero}'
        comment = f'{caero.comment}\n{str(caero)}'
        model2.coords[caero.cp] = caero.cp_ref
        caero7 = CAERO7(caero.eid, label,
                        caero.p1, caero.x12, caero.p4, caero.x43,
                        caero.cp, caero.nspan, caero.nchord,
                        caero.lspan,
                        comment=comment)
        if caero.lspan:
            assert caero.lspan_ref.type == 'AEFACT', caero.lspan_ref
        if caero.lchord:
            model2.log.warning('lchord is not suppported')
        panlst_id = _get_panlst_id_from_caero(caero)
        model2.caeros[caero_id] = caero7
        ncaero7 += 1
        icaero += 1
    return ncaero7

def _convert_splines(model: BDF,
                     model2: BDF) -> tuple[int, int, int, int]:
    cp = 950
    ncord2r = 0
    nspline1 = 0
    nspline2 = 0
    nspline3 = 0
    for spline_id, spline in model.splines.items():
        assert spline.type == 'SPLINE1', spline
        # model_str = ''
        if spline.type == 'SPLINE1':
            caero_id = spline.caero
            caero = model2.caeros[caero_id]
            panlst_id = _get_panlst_id_from_caero(caero)

            comment = str(spline.comment)
            model_str = caero.label

            spline = cast(SPLINE1, spline)
            setg = int(spline.setg)
            set_card = model.sets[setg]
            set_ids = set_card.ids

            xyz_list = [model.nodes[nid].get_position() for nid in set_ids]
            xyzs = np.vstack(xyz_list)
            origin, zaxis, xzplane = fit_plane_to_point_cloud(xyzs)
            origin = np.zeros(3)
            model2.add_cord2r(cp, origin, origin+zaxis, origin+xzplane)
            ncord2r += 1

            cp2 = None
            model2.sets[setg] = set_card
            splinez = SPLINE1_ZAERO(
                spline_id, panlst_id, setg, model_str, cp2,
                spline.dz, eps=0.01, comment=comment)
            model2.splines[spline_id] = splinez
            nspline1 += 1
            cp += 1
    return ncord2r, nspline1, nspline2, nspline3


def fit_plane_to_point_cloud(point_cloud: np.ndarray,) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fits a best-fit plane to a 3D point cloud using SVD.

    Args:
        point_cloud (numpy.ndarray): An array of shape (N, 3) representing N points (x, y, z).

    Returns:
        tuple: (normal_vector, centroid) of the best-fit plane.
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
