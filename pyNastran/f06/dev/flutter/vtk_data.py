from typing import Any

DT_MS_DEFAULT = 100
NPHASE_DEFAULT = 10
DT_MS_MIN = 100
DT_MS_MAX = 5000


class VtkData:
    def __init__(self):
        self.nphase = NPHASE_DEFAULT
        self.icase = 0
        self.ncase = 1
        self.animate = True
        self.dt_ms = DT_MS_DEFAULT
        self.global_scale_factor = 1.0
        self.point_size = 2

    def apply_settings(self, data: dict[str, Any]) -> None:
        apply_vtk_settings(self, data)

    def to_json(self) -> dict[str, Any]:
        return {
            'nphase': self.nphase,
            'icase': self.icase,
            'ncase': self.ncase,
            'animate': self.animate,
            'dt_ms': self.dt_ms,
        }

def apply_vtk_settings(self, data: dict[str, Any]):
    vtk_data = data.get('vtk', {})
    self.dt_ms = DT_MS_DEFAULT if 'dt_ms' not in vtk_data else int(vtk_data['dt_ms'])
    self.nphase = NPHASE_DEFAULT if 'nphase' not in vtk_data else int(vtk_data['nphase'])
    self.icase = 0 if 'icase' not in vtk_data else int(vtk_data['icase'])
    self.animate = True
