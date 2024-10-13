from typing import Callable, Union
from vtkmodules.vtkRenderingCore import vtkActor
from vtkmodules.vtkRenderingLOD import vtkLODActor

# fmt_name, major_name, geom_wildcard, geom_func, res_wildcard, _resfunc
Format = tuple[str, str, str, Callable, str, Callable]

ColorFloat = tuple[float, float, float]
ColorInt = tuple[int, int, int]

Actor = vtkActor | vtkLODActor
