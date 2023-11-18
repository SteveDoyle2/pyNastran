import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF

def get_temperatures_array(model: BDF, loadcase_id: int, fdtype: str='float32') -> np.ndarray:
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
        grid_id = model.grid.node_id
        grid_temps = temperature[:ngrid]
        igrid = np.searchsorted(grid_id, temp.node_id)
        grid_temps[igrid] = temp.temperature
    return temperature
