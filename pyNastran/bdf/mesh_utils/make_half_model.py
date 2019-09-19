"""
This file defines:
  - model = make_half_model(
        bdf_filename, plane='xz', zero_tol=1e-12,
        log=None, debug=True)

"""
from cpylog import get_logger2
from pyNastran.utils import deprecated
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
from pyNastran.bdf.mesh_utils.mirror_mesh import _plane_to_iy


def make_symmetric_model(bdf_filename, plane: str='xz',
                         zero_tol: float=1e-12, log=None, debug: bool=True):
    log = get_logger2(log=log, debug=debug, encoding='utf-8')
    deprecated('make_symmetric_model', 'make_half_model', '1.3', levels=[0, 1, 2])
    return make_half_model(bdf_filename, plane=plane, zero_tol=zero_tol, log=log, debug=debug)

def make_half_model(bdf_filename, plane: str='xz',
                         zero_tol: float=1e-12, log=None, debug: bool=True):
    """
    Makes a half/symmetric model from a full model

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z
    zaero_tol : float; default=1e-12
        the symmetry plane tolerance

    Returns
    -------
    model : BDF()
        BDF : the BDF model object

    ## TODO: doesn't handle elements straddling the centerline

    """
    model = get_bdf_model(bdf_filename, xref=True, log=log, debug=debug)
    iy, plane = _plane_to_iy(plane)
    nids_to_remove = []
    eids_to_remove = []
    caero_ids_to_remove = []
    zero = -zero_tol
    for eid, elem in model.elements.items():
        xyz = elem.Centroid()

        if xyz[iy] < zero:
            eids_to_remove.append(eid)

    for nid, node in model.nodes.items():
        xyz = node.get_position()
        if xyz[iy] < zero:
            nids_to_remove.append(nid)

    for nid in nids_to_remove:
        del model.nodes[nid]

    for eid in eids_to_remove:
        del model.elements[eid]

    for caero_id, caero in model.caeros.items():
        if caero.type == 'CAERO1':
            p1, p2, p3, p4 = caero.get_points()
            #print(caero)
            if p1[iy] <= zero and p4[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove.append(caero_id)
            elif p1[iy] < zero:
                p1[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
            elif p4[iy] < zero:
                p4[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
        elif caero.type == 'CAERO2':
            # TODO: a CAERO2 can't be half symmetric...can it?
            # TODO: it can be skewed though...
            p1, p2 = caero.get_points()
            if p1[iy] <= zero and p2[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove.append(caero_id)
        else:  # pragma: no cover
            raise NotImplementedError(caero)

    for caero_id in caero_ids_to_remove:
        del model.caeros[caero_id]

    #print('nids_to_remove =', nids_to_remove)
    for unused_spline_id, spline in model.splines.items():
        caero = spline.caero
        #setg = spline.setg
        #print('caero = ', caero)
        nids = spline.setg_ref.ids  # list
        #spline.uncross_reference()

        #i = 0
        nids = list(set(nids) - set(nids_to_remove))
        nids.sort()
        spline.setg_ref.ids_ref = None
        spline.setg_ref.ids = nids

    plane_to_labels_keep_map = {
        'yz' : ['URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW'], # yz
        'xz' : ['URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA'], # xz plane
        'xy' : ['URDD1', 'URDD2', 'URDD6', 'ROLL'], # xy plane
    }

    all_labels = {
        'URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW',
        'URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA',
        'URDD1', 'URDD2', 'URDD6', 'ROLL',
    }
    labels_to_keep = plane_to_labels_keep_map[plane]
    labels_to_remove = [label for label in all_labels if label not in labels_to_keep]

    #print('labels_to_remove =', labels_to_remove)
    for aestat_id in list(model.aestats.keys()):
        aestat = model.aestats[aestat_id]
        if aestat.label in labels_to_remove:
            del model.aestats[aestat_id]

    for unused_trim_id, trim in model.trims.items():
        labels = trim.labels
        ilabels_to_remove = [labels.index(label) for label in labels_to_remove
                             if label in labels]
        #print("ilabels_to_remove =", ilabels_to_remove)
        trim.uxz = [trim.uxs[ilabel] for ilabel in ilabels_to_remove]
        trim.labels = [trim.labels[ilabel] for ilabel in ilabels_to_remove]
    return model
