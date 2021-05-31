from __future__ import annotations
from typing import Optional, List, TYPE_CHECKING
if TYPE_CHECKING:
    import numpy as np
    from pyNastran.gui.dev.gui2.gui2 import MainWindow2

from pyNastran.gui.qt_files.colors import BLACK_FLOAT


class ScalarBar:
    def __init__(self):
        pass
    def Modified(self):
        pass
    def VisibilityOn(self):
        pass
    def VisibilityOff(self):
        pass


class VtkInterface:
    def __init__(self, gui: MainWindow2):
        self.gui = gui
        self.scalar_bar_actor = ScalarBar()

    @property
    def log(self):
        return self.gui.log

    def set_quad_grid(self, box_name: str,
                      nodes: np.ndarray, elements: np.ndarray,
                      color: Optional[List[float]]=None,
                      line_width: float=1, opacity: float=1.) -> None:
        if color is None:
            color = BLACK_FLOAT
        self.log.warning('set_quad_grid')

    def create_global_axes(self, dim_max: float) -> None:
        self.log.warning('create_global_axes')
