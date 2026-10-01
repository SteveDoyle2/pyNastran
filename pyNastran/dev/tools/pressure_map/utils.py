import numpy as np


def log_range(log,
              xyz: np.ndarray,
              xyz_units: str) -> None:
    xyz_min = xyz.min(axis=0)
    xyz_max = xyz.max(axis=0)
    dxyz = xyz_max - xyz_min
    log(f'    xyz_min  ({xyz_units}) = {xyz_min}')
    log(f'    xyz_max  ({xyz_units}) = {xyz_max}')
    log(f'    dxyz     ({xyz_units}) = {dxyz}')
