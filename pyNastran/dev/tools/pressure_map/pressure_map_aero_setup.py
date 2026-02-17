import os
from typing import cast, Any
import numpy as np

from pyNastran.utils import PathLike, print_bad_path

from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
from pyNastran.converters.fluent.fluent import read_fluent, Fluent, filter_by_region
from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot


def get_aero_model(aero_filename: PathLike, aero_format: str,
                   aero_xyz_scale: float=1.0,
                   xyz_units: str='???',
                   stop_on_failure: bool=True) -> tuple[Any, list[str]]:
    if isinstance(aero_filename, PathLike):
        assert os.path.exists(aero_filename), print_bad_path(aero_filename)

    # if regions_to_include is None:
    #     regions_to_include = []
    # if regions_to_remove is None:
    #     regions_to_remove = []

    aero_format = aero_format.lower()
    if aero_format == 'cart3d':
        if isinstance(aero_filename, Cart3D):
            model = aero_filename
        else:
            model: Cart3D = read_cart3d(aero_filename)
        log = model.log
        variables = list(model.loads)
        model.points *= aero_xyz_scale
        xyz = model.points
    # elif aero_format == 'Fund3D':
    #     return None, []
    # elif aero_format == 'Fluent Press':
    #     pass
    elif aero_format == 'tecplot':
        model: Tecplot = read_tecplot(aero_filename)
        log = model.log
        #print(model.object_stats())
        variables = model.result_variables
        model.xyz *= aero_xyz_scale
        xyz = model.xyz
    elif aero_format == 'fluent':
        if isinstance(aero_filename, Fluent):
            model = aero_filename
        else:
            model: Fluent = read_fluent(aero_filename, debug='debug')
        log = model.log
        # print(model.object_stats())
        variables = model.titles[1:]
        model.xyz *= aero_xyz_scale
        xyz = model.xyz
    else:  # pragma: no cover
        if stop_on_failure:
            raise NotImplementedError(aero_format)
        else:
            return None, []

    xyz_min = xyz.min(axis=0)
    xyz_max = xyz.max(axis=0)
    dxyz = xyz_max - xyz_min
    log.info(f'aero xyz range (aero_xyz_scale={aero_xyz_scale}):')
    log.info(f'    xyz_min  ({xyz_units}) = {xyz_min}')
    log.info(f'    xyz_max  ({xyz_units}) = {xyz_max}')
    log.info(f'    dxyz     ({xyz_units}) = {dxyz}')

    return model, variables


def get_aero_pressure_centroid(aero_model: Cart3D | Tecplot | Fluent,
                               aero_format: str,
                               map_type: str,
                               variable: str='Cp',
                               regions_to_include=None,
                               regions_to_remove=None,
                               idtype: str='int32',
                               fdtype: str='float64') -> dict[str, np.ndarray]:
    """
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
        area = np.hstack([quad_area, tri_area])
        quad_Cp_centroid = quad_results[:, iresult]
        tri_Cp_centroid = tri_results[:, iresult]
        Cp_centroid = np.hstack([
            quad_Cp_centroid,
            tri_Cp_centroid,
        ])
    else:  # pragma: no cover
        raise RuntimeError(aero_format)

    assert xyz_nodal.shape[1] == 3, xyz_nodal
    log.info(f'aero: nnodes={len(node_id)} nxyz={len(xyz_nodal)} '
             f'ntris={len(tri_nodes)} nquads={len(quad_nodes)} ncentroid={len(centroid)} narea={len(area)}')

    ntri = len(tri_nodes)
    nquad = len(quad_nodes)

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
