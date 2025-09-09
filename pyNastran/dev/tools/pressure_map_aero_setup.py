import os
from typing import cast, Any
import numpy as np

from pyNastran.utils import PathLike, print_bad_path

from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
from pyNastran.converters.fluent.fluent import read_fluent, Fluent, filter_by_region
from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot


def get_aero_model(aero_filename: PathLike, aero_format: str,
                   aero_xyz_scale: float=1.0,
                   xyz_units: float='in',
                   stop_on_failure: bool=True) -> tuple[Any, list[str]]:
    assert os.path.exists(aero_filename), print_bad_path(aero_filename)

    aero_format = aero_format.lower()
    if aero_format == 'cart3d':
        model: Cart3D = read_cart3d(aero_filename)
        log = model.log
        variables = list(model.loads)
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
    elif aero_format == 'fluent':
        model: Fluent = read_fluent(aero_filename, debug='debug')
        log = model.log
        # print(model.object_stats())
        variables = model.titles[1:]
        model.xyz *= aero_xyz_scale
    else:  # pragma: no cover
        if stop_on_failure:
            raise NotImplementedError(aero_format)
        else:
            return None, []

    xyz_min = model.xyz.min(axis=0)
    xyz_max = model.xyz.max(axis=0)
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
                               regions_to_remove=None) -> dict[str, np.ndarray]:
    """
    variable: str
        cart3d: 'Cp' only; variable is ignored
        tecplot: use variable
        fun3d:        ???
        fluent press: ???
        fluent vrt:   ???
    """
    if regions_to_include is None:
        regions_to_include = []
    if regions_to_remove is None:
        regions_to_remove = []

    log = aero_model.log
    aero_format = aero_format.lower()
    if aero_format == 'cart3d':
        aero_model = cast(Cart3D, aero_model)
        assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
        aero_elements = aero_model.elements
        aero_xyz_nodal = aero_model.nodes
        xyz1 = aero_xyz_nodal[aero_elements[:, 0], :]
        xyz2 = aero_xyz_nodal[aero_elements[:, 1], :]
        xyz3 = aero_xyz_nodal[aero_elements[:, 2], :]
        #aero_xyz_centroid = (xyz1 + xyz2 + xyz3) / 3
        aero_centroid = (xyz1 + xyz2 + xyz3) / 3.
        aero_normal = np.cross(xyz2-xyz1, xyz3-xyz1)
        aero_area = 0.5 * np.linalg.norm(aero_normal)
        aero_Cp_centroidal = aero_model.loads['Cp']
        assert len(aero_Cp_centroidal) == len(aero_elements)
        assert len(aero_Cp_centroidal) == len(aero_xyz)
    elif aero_format == 'tecplot':
        aero_model = cast(Tecplot, aero_model)
        #raise RuntimeError(aero_model)
    elif aero_format == 'fluent':
        aero_model = cast(Fluent, aero_model)
        (
            tri_area, quad_area,
            tri_centroid, quad_centroid,
            tri_normal, quad_normal,
        ) = aero_model.get_area_centroid_normal(aero_model.tris, aero_model.quads)

        # tri_nodes = tris[:, 2:]
        # quad_nodes = quads[:, 2:]
        aero_xyz_nodal = aero_model.xyz
        if len(regions_to_include) > 0 or len(regions_to_remove) > 0:
            element_id, tris, quads, quad_results, tri_results = filter_by_region(
                aero_model, regions_to_remove, regions_to_include)
        else:
            log.info(f'regions_to_include={regions_to_include}')
            log.info(f'regions_to_remove={regions_to_remove}')
            tri_eids = aero_model.tris[:, 0]
            tris = aero_model.tris[:, 2:]
            quads = aero_model.quads[:, 2:]

            quad_eids = aero_model.tris[:, 0]
            element_id = np.unique(np.hstack([quad_eids, tri_eids]))
            tri_results = aero_model.results[:len(tris)]
            quad_results = aero_model.results[:len(quads)]
            assert len(aero_model.result_element_id) > 0, aero_model.result_element_id

        node_id = aero_model.node_id
        itri = np.searchsorted(node_id, tris)
        iquad = np.searchsorted(node_id, quads)

        xyz1 = aero_xyz_nodal[itri[:, 0], :]
        xyz2 = aero_xyz_nodal[itri[:, 1], :]
        xyz3 = aero_xyz_nodal[itri[:, 2], :]
        tri_centroid = (xyz1 + xyz2 + xyz3) / 3.
        tri_normal = np.cross(xyz2 - xyz1, xyz3 - xyz1)

        xyz1 = aero_xyz_nodal[iquad[:, 0], :]
        xyz2 = aero_xyz_nodal[iquad[:, 1], :]
        xyz3 = aero_xyz_nodal[iquad[:, 2], :]
        xyz4 = aero_xyz_nodal[iquad[:, 3], :]
        quad_centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
        quad_normal = np.cross(xyz3 - xyz1, xyz4 - xyz2)

        aero_normal = np.vstack([tri_normal, quad_normal])
        aero_centroid = np.vstack([tri_centroid, quad_centroid])

        titles = list(aero_model.titles[1:])
        iresult = titles.index('Pressure Coefficient')
        aero_area = np.hstack([quad_area, tri_area])
        aero_Cp_centroidal = np.hstack([
            quad_results[:, iresult],
            tri_results[:, iresult],
        ])
        #aero_model = read_fluent(aero_filename)
    else:  # pragma: no cover
        raise RuntimeError(aero_format)

    assert aero_xyz_nodal.shape[1] == 3, aero_xyz_nodal
    aero_dict = {
        'xyz_nodal': aero_xyz_nodal,
        'Cp_centroidal': aero_Cp_centroidal,
        'centroid': aero_centroid,
        'area': aero_area,
        'normal': aero_normal,
    }
    # aero_Cp_centroidal = out_dict['aero_Cp_centroidal']
    # aero_area = out_dict['aero_area']
    return aero_dict
