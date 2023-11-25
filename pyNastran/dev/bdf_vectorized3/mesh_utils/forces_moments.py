import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF

def get_temperatures_array(model: BDF, loadcase_id: int,
                           fdtype: str='float32') -> np.ndarray:
    """gets a static temperature array; supports GRID/SPOINT"""
    ngrid = len(model.grid)
    nspoint = len(model.spoint)
    ndof = ngrid + nspoint

    default_temp = 0.
    if loadcase_id in model.tempd.load_id:
        tempd = model.tempd.slice_card_by_id(loadcase_id)
        default_temp = tempd.temperature[0]
    temperature = np.full(ndof, default_temp, dtype=fdtype)

    if len(model.temp) and loadcase_id in model.temp:
        temp = model.temp.slice_card_by_id(loadcase_id)

        if ngrid:
            grid_id = model.grid.node_id
            grid_temps = temperature[:ngrid]

            # A+B
            ugrid_id = np.intersect1d(grid_id, temp.node_id)
            # now that we know the common set, apply grid loads to only temp grids
            if len(ugrid_id):
                igrid = np.searchsorted(grid_id, ugrid_id)
                itemp = np.searchsorted(temp.node_id, ugrid_id)
                grid_temps[igrid] = temp.temperature[itemp]
        if nspoint:
            spoint_id = model.spoint.ids
            uspoint_id = np.intersect1d(spoint_id, temp.node_id)
            spoint_temps = temperature[ngrid:]

            # now that we know the common set, apply grid loads to only temp grids
            if len(uspoint_id):
                ispoint = np.searchsorted(spoint_id, uspoint_id)
                itemp = np.searchsorted(temp.node_id, uspoint_id)
                spoint_temps[ispoint] = temp.temperature[itemp]
    return temperature
