from typing import Union
from .static_loads import (
    FORCE, FORCE1, FORCE2,
    MOMENT, MOMENT1, MOMENT2,
    GRAV, ACCEL, ACCEL1,
    SLOAD, TEMP, TEMPD,
    RFORCE, RFORCE1)
from .static_pressure_loads import PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1

Loads = Union[FORCE, FORCE1, FORCE2,
              MOMENT, MOMENT1, MOMENT2,
              PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
              GRAV, ACCEL, ACCEL1,
              SLOAD, TEMP, TEMPD,
              RFORCE, RFORCE1]
