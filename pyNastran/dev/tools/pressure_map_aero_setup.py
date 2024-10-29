import os
from typing import cast, Any
import numpy as np

from pyNastran.utils import PathLike, print_bad_path

from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
from pyNastran.converters.fluent.fluent import read_fluent, Fluent
from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot


def get_aero_model(aero_filename: str, aero_format: str,
                   stop_on_failure=True) -> tuple[Any, list[str]]:
    assert os.path.exists(aero_filename), print_bad_path(aero_filename)
    if aero_format == 'Cart3D':
        model = read_cart3d(aero_filename)
        variables = list(model.loads)
    # elif aero_format == 'Fund3D':
    #     return None, []
    # elif aero_format == 'Fluent Press':
    #     pass
    elif aero_format == 'Tecplot':
        model: Tecplot = read_tecplot(aero_filename)
        #print(model.object_stats())
        variables = model.result_variables
    elif aero_format == 'Fluent':
        model: Fluent = read_fluent(aero_filename)
        # print(model.object_stats())
        variables = model.titles[1:]
    else:  # pragma: no cover
        if stop_on_failure:
            raise NotImplementedError(aero_format)
        else:
            return None, []
    return model, variables

def get_aero_pressure_centroid(aero_model: Cart3D | Tecplot | Fluent,
                               aero_format: str,
                               map_type: str,
                               variable: str='Cp') -> dict[str, np.ndarray]:
    """
    variable: str
        cart3d: 'Cp' only; variable is ignored
        tecplot: use variable
        fun3d:        ???
        fluent press: ???
        fluent vrt:   ???
    """
    if aero_format == 'cart3d':
        aero_model = cast(Cart3D, aero_model)
        assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
        aero_elements = aero_model.elements
        aero_xyz_nodal = aero_model.nodes
        xyz1 = aero_xyz_nodal[aero_elements[:, 0], :]
        xyz2 = aero_xyz_nodal[aero_elements[:, 1], :]
        xyz3 = aero_xyz_nodal[aero_elements[:, 2], :]
        #aero_xyz_centroid = (xyz1 + xyz2 + xyz3) / 3
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
        ) = aero_model.get_area_centroid_normal()
        aero_xyz_nodal = aero_model.nodes
        regions_to_remove = []
        regions_to_include = []
        element_id, tris, quads, quad_results, tri_results = aero_model.filter_by_region(
            regions_to_remove, regions_to_include)
        iresult = aero_model.titles[1:].index('Pressure Coefficient')
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
        'xyz_nodal' : aero_xyz_nodal,
        'Cp_centroidal' : aero_Cp_centroidal,
        'area' : aero_area,
    }
    # aero_Cp_centroidal = out_dict['aero_Cp_centroidal']
    # aero_area = out_dict['aero_area']
    return aero_dict
