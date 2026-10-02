import os
from typing import cast, Any, Optional
import numpy as np

from cpylog import SimpleLogger
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
from pyNastran.converters.fluent.fluent import read_fluent, Fluent, filter_by_region
from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot
from .utils import log_range

def get_aero_model(aero_filename: PathLike, aero_format: str,
                   aero_xyz_scale: float=1.0,
                   xyz_units_in: str='???',
                   xyz_units_out: str='???',
                   stop_on_failure: bool=True,
                   qinf: float=1.0,
                   reference_point: Optional[np.ndarray]=None,
                   pressure_units: str='',
                   regions_to_include=None,
                   log: Optional[SimpleLogger]=None) -> tuple[Any, list[str]]:
    """
    Loads (or accepts) an aero surface model and scales its xyz to
    structure units.

    Parameters
    ----------
    aero_filename : PathLike | Cart3D | Fluent
        a file path, or an already-loaded model.  Tecplot must be a path.
    aero_format : str
        'cart3d', 'fluent', 'tecplot'
    aero_xyz_scale : float; default=1.0
        xyz_structure = aero_xyz_scale * xyz_aero, applied IN PLACE.
          39.3701    aero m  -> structure in
          0.0393701  aero mm -> structure in
          1000.      aero m  -> structure mm
        Passing a loaded model scales that same object again, so call
        this once per model.  If you then hand the model to
        ``pressure_map``, pass ``aero_xyz_scale=1.0`` there.
    xyz_units_in / xyz_units_out : str; default='???'
        label for the scaled (structure) length unit; log-only, does not
        scale anything
    stop_on_failure : bool; default=True
        raise on an unknown aero_format (else return None, [])
    qinf : float; default=1.0
        dynamic pressure for the logged force/moment sum (structure
        pressure units).  1.0 logs F/q and M/q.
    reference_point : (3,) float ndarray; default=None -> [0, 0, 0]
        moment reference point in the SCALED (structure) units, basic
    pressure_units : str; default=''
        label for qinf; log only

    Before returning, logs the aero force/moment sum over every aero
    element (F = sum(qinf*Cp*A*n), M = sum((c - reference_point) x F))
    using the scaled xyz.  It's skipped with a warning if the model has
    no Cp ('Cp' for Cart3D, 'Pressure Coefficient' for Fluent) or is
    Tecplot.  Log only; the model isn't changed.

    Returns
    -------
    model : Cart3D | Fluent | Tecplot
        the (scaled) aero model; the same object if one was passed in
    variables : list[str]
        result names (e.g. ['Pressure Coefficient'] for Fluent)

    Cp, not pressure, is read from the model, so no pressure conversion
    happens here; dimensional pressure is qinf*Cp in ``pressure_map``.
    """
    if log is None:
        log = SimpleLogger(level='debug')
    # log.debug(f'aero_filename={aero_filename}')
    # log.info(f'get_aero_model')
    log.info(f'get_aero_model (aero_xyz_scale={aero_xyz_scale}):')
    if isinstance(aero_filename, PathLike):
        assert os.path.exists(aero_filename), print_bad_path(aero_filename)

    # if regions_to_include is None:
    #     regions_to_include = []
    # if regions_to_remove is None:
    #     regions_to_remove = []

    if aero_xyz_scale == 1.0 and not isinstance(aero_filename, PathLike):
        if isinstance(aero_filename, Fluent):
            model = aero_filename
            variables = model.titles[1:]
            # xyz = model.xyz
            log_aero_force_moment(
                model, 'fluent', qinf=qinf, reference_point=reference_point,
                xyz_units=xyz_units_out, pressure_units=pressure_units)
            return model, variables

    aero_format = aero_format.lower()
    if aero_format == 'cart3d':
        if isinstance(aero_filename, Cart3D):
            model = aero_filename
        else:
            model: Cart3D = read_cart3d(aero_filename, log=log)
        log = model.log
        variables = list(model.loads)
        log.info(f'aero xyz range:')
        log_range(log.debug, model.points, xyz_units_in)  # Cart3D has no .xyz
        model.points *= aero_xyz_scale
        xyz = model.points
    # elif aero_format == 'Fund3D':
    #     return None, []
    # elif aero_format == 'Fluent Press':
    #     pass
    elif aero_format == 'tecplot':
        model: Tecplot = read_tecplot(aero_filename, log=log)
        log = model.log
        #print(model.object_stats())
        variables = model.result_variables
        model.xyz *= aero_xyz_scale
        log.info(f'aero xyz range:')
        log_range(log.debug, model.xyz, xyz_units_in)
        xyz = model.xyz
    elif aero_format == 'fluent':
        if isinstance(aero_filename, Fluent):
            model = aero_filename
        else:
            model: Fluent = read_fluent(aero_filename, debug='debug', log=log)
        log = model.log
        # print(model.object_stats())
        variables = model.titles[1:]
        log.info(f'aero xyz range:')
        log_range(log.debug, model.xyz, xyz_units_in)
        model.xyz *= aero_xyz_scale
        xyz = model.xyz
    else:  # pragma: no cover
        raise NotImplementedError(aero_format)

    log.info(f'aero xyz range (aero_xyz_scale={aero_xyz_scale}):')
    log_range(log.debug, xyz, xyz_units_out)
    log_aero_force_moment(
        model, aero_format, qinf=qinf, reference_point=reference_point,
        xyz_units=xyz_units_out, pressure_units=pressure_units)
    return model, variables


def log_aero_force_moment(model: Cart3D | Tecplot | Fluent,
                          aero_format: str,
                          qinf: float=1.0,
                          reference_point: Optional[np.ndarray]=None,
                          xyz_units: str='',
                          pressure_units: str='') -> Optional[tuple[np.ndarray, np.ndarray]]:
    """
    Logs the aero force/moment over every aero element.

    F = sum(qinf*Cp_i*A_i*n_i);  M = sum((c_i - reference_point) x F_i)

    Uses the model's current xyz (so the scaled units) and the element
    normals as the mesh defines them (no flipping).

    Returns
    -------
    (force, moment) : (3,) float ndarrays, or None if skipped
    """
    log = model.log
    aero_format = aero_format.lower()
    if reference_point is None:
        reference_point = np.zeros(3)
    reference_point = np.asarray(reference_point, dtype='float64')

    if aero_format == 'cart3d':
        cp_name = 'Cp'
        has_cp = cp_name in model.loads
    elif aero_format == 'fluent':
        cp_name = 'Pressure Coefficient'
        has_cp = cp_name in list(model.titles[1:])
    else:
        log.warning(f'aero force/moment sum: not supported for aero_format={aero_format!r}; skipping')
        return None
    if not has_cp:
        log.warning(f'aero force/moment sum: no {cp_name!r} result; skipping')
        return None

    aero_dict = get_aero_pressure_centroid(model, aero_format, 'pressure')
    force_i = (qinf * aero_dict['Cp_centroid'] * aero_dict['area'])[:, np.newaxis] * aero_dict['normal']
    moment_i = np.cross(aero_dict['centroid'] - reference_point, force_i)
    force = force_i.sum(axis=0)
    moment = moment_i.sum(axis=0)
    area = aero_dict['area'].sum()

    log.info(f'aero force/moment sum (all {len(force_i)} aero elements):')
    log.info(f'  qinf={qinf} {pressure_units}; reference_point={reference_point.tolist()} {xyz_units}; '
             f'area={area:g} {xyz_units}^2')
    log.info(f'  F = [{force[0]:g}, {force[1]:g}, {force[2]:g}]')
    log.info(f'  M = [{moment[0]:g}, {moment[1]:g}, {moment[2]:g}]')
    return force, moment


def get_aero_pressure_centroid(aero_model: Cart3D | Tecplot | Fluent,
                               aero_format: str,
                               map_type: str,
                               variable: str='Cp',
                               regions_to_include=None,
                               regions_to_remove=None,
                               idtype: str='int32',
                               fdtype: str='float64') -> dict[str, np.ndarray]:
    """
    Gets per-element Cp, area, centroid, and unit normal.

    All lengths come from the model's current xyz, so they are in
    whatever units ``get_aero_model`` scaled to (area in length**2).
    Cp is returned as-is (no qinf applied).

    variable: str
        cart3d: 'Cp' only; variable is ignored
        tecplot: use variable
        fun3d:        ???
        fluent press: ???
        fluent vrt:   ???
    """
    log = aero_model.log
    if regions_to_include is None:
        regions_to_include = []
    if regions_to_remove is None:
        regions_to_remove = []

    log = aero_model.log
    aero_format = aero_format.lower()
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format

    quad_nodes = np.zeros((0, 4), dtype=idtype)
    quad_area = np.array([], dtype=fdtype)
    quad_centroid = np.zeros((0, 3), dtype=fdtype)
    quad_normal = np.zeros((0, 3), dtype=fdtype)
    quad_Cp_centroid = np.array([], dtype=fdtype)
    if aero_format == 'cart3d':
        aero_model = cast(Cart3D, aero_model)
        tri_nodes = aero_model.elements
        xyz_nodal = aero_model.nodes
        node_id = np.arange(len(aero_model.nodes), dtype=idtype)
        ntri = len(tri_nodes)
        xyz1 = xyz_nodal[tri_nodes[:, 0], :]
        xyz2 = xyz_nodal[tri_nodes[:, 1], :]
        xyz3 = xyz_nodal[tri_nodes[:, 2], :]
        tri_centroid = (xyz1 + xyz2 + xyz3) / 3.
        tri_normal = np.cross(xyz2-xyz1, xyz3-xyz1, axis=1)
        normi =np.linalg.norm(tri_normal, axis=1)
        assert tri_normal.shape == xyz1.shape
        assert len(normi) == ntri, (len(normi), ntri)
        tri_normal /= normi[:, np.newaxis]
        tri_area = 0.5 * normi
        tri_Cp_centroid = aero_model.loads['Cp']
        assert len(tri_area) == ntri
        assert len(tri_Cp_centroid) == ntri
        centroid = tri_centroid
        area = tri_area
        normal = tri_normal
        Cp_centroid = tri_Cp_centroid
    elif aero_format == 'tecplot':
        aero_model = cast(Tecplot, aero_model)
        #raise RuntimeError(aero_model)
    elif aero_format == 'fluent':
        aero_model = cast(Fluent, aero_model)
        xyz_nodal = aero_model.xyz

        # if len(regions_to_include) > 0 or len(regions_to_remove) > 0:
        log.debug(f'regions_to_remove={regions_to_remove} regions_to_include={regions_to_include}')
        element_id, tris, quads, quad_results, tri_results = filter_by_region(
            aero_model, regions_to_remove, regions_to_include)
        tri_eids = tris[:, 0]
        quad_eids = quads[:, 0]

        tri_nodes = tris[:, 2:]
        quad_nodes = quads[:, 2:]

        assert len(tri_results) == len(tri_nodes), (len(tri_results), len(tri_nodes))
        (
            tri_area, quad_area,
            tri_centroid, quad_centroid,
            tri_normal, quad_normal,
        ) = aero_model.get_area_centroid_normal_from_nodes(tri_nodes, quad_nodes)

        node_id = aero_model.node_id
        assert len(tri_results) == len(tri_area), (len(tri_results), len(tri_area))

        normal = np.vstack([tri_normal, quad_normal])
        centroid = np.vstack([tri_centroid, quad_centroid])

        titles = list(aero_model.titles[1:])
        name = 'Pressure Coefficient'
        try:
            iresult = titles.index(name)
        except ValueError:
            log.error(f'cant find {name!r} in titles={list(aero_model.titles)}')
            raise
        # all combined arrays are stacked [tri, quad] so row i of
        # area/Cp_centroid/centroid/normal is the same element
        area = np.hstack([tri_area, quad_area])
        quad_Cp_centroid = quad_results[:, iresult]
        tri_Cp_centroid = tri_results[:, iresult]
        Cp_centroid = np.hstack([
            tri_Cp_centroid,
            quad_Cp_centroid,
        ])
    else:  # pragma: no cover
        raise RuntimeError(aero_format)

    assert xyz_nodal.shape[1] == 3, xyz_nodal
    if len(tri_nodes) + len(quad_nodes) == 0:
        raise RuntimeError(
            f'no aero elements (ntris=0, nquads=0); '
            f'regions_to_include={regions_to_include} '
            f'regions_to_remove={regions_to_remove}')
    log.info(f'aero: nnodes={len(node_id)} nxyz={len(xyz_nodal)} '
             f'ntris={len(tri_nodes)} nquads={len(quad_nodes)} ncentroid={len(centroid)} narea={len(area)}')

    if len(Cp_centroid):
        msg = f'aero Cp: min={Cp_centroid.min():g} max={Cp_centroid.max():g}'
        if np.abs(Cp_centroid).max() > 2:
            log.warning(msg)
        else:
            log.info(msg)

    ntri = len(tri_nodes)
    nquad = len(quad_nodes)
    nelem = ntri + nquad
    assert len(area) == nelem, (len(area), nelem)
    assert len(Cp_centroid) == nelem, (len(Cp_centroid), nelem)
    assert len(centroid) == nelem, (len(centroid), nelem)
    assert len(normal) == nelem, (len(normal), nelem)

    assert len(tri_area) == ntri, (len(tri_area), ntri)
    assert len(tri_centroid) == ntri, (len(tri_centroid), ntri)
    assert len(tri_normal) == ntri, (len(tri_normal), ntri)
    assert len(tri_Cp_centroid) == ntri, (len(tri_Cp_centroid), ntri)

    if ntri:
        assert np.abs(tri_normal).max() <= 1.001, (tri_normal.min(), tri_normal.max())
    if nquad:
        assert np.abs(quad_normal).max() <= 1.001, (quad_normal.min(), quad_normal.max())

    aero_dict = {
        'node_id': node_id,
        'xyz_nodal': xyz_nodal,
        'Cp_centroid': Cp_centroid,
        'centroid': centroid,
        'area': area,
        'normal': normal,

        'tri_nodes': tri_nodes,
        'tri_centroid': tri_centroid,
        'tri_area': tri_area,
        'tri_Cp_centroid': tri_Cp_centroid,
        'tri_normal': tri_normal,

        'quad_nodes': quad_nodes,
        'quad_centroid': quad_centroid,
        'quad_area': quad_area,
        'quad_Cp_centroid': quad_Cp_centroid,
        'quad_normal': quad_normal,
    }
    # aero_Cp_centroidal = out_dict['aero_Cp_centroidal']
    # aero_area = out_dict['aero_area']
    return aero_dict
