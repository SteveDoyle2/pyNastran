from collections import defaultdict
import numpy as np
TrimVariable = tuple[int, str, str, float, str]
ControllerState = dict[str, float]
TrimVariables = dict[str, TrimVariable]


class MonitorLoads:
    def __init__(self, name: np.ndarray,
                 comp: np.ndarray,
                 classi: np.ndarray,
                 label: np.ndarray,
                 cp: np.ndarray,
                 xyz: np.ndarray,
                 coefficient: np.ndarray,
                 cd: np.ndarray):
        self.name = name
        self.comp = comp
        self.group = classi
        self.label = label
        self.cp = cp
        self.xyz = xyz
        self.coefficient = coefficient
        self.cd = cd

    def __repr__(self) -> str:
        msg = f'MonitorLoads; n={len(self.name)}'
        return msg


class AeroPressure:
    def __init__(self, subcase: int, title: str, subtitle: str,
                 mach: float, q: float,
                 cref: float, bref: float, sref: float,
                 nodes: np.ndarray,
                 cp: np.ndarray, pressure: np.ndarray):
        self.subcase = subcase
        self.title = title
        self.subtitle = subtitle
        self.mach = mach
        self.q = q
        self.cref = cref
        self.bref = bref
        self.sref = sref

        self.nodes = nodes
        self.pressure = pressure
        self.cp = cp

    @classmethod
    def from_f06(self, subcase: int,
                 nodes: np.ndarray,
                 cp_pressure: np.ndarray):
        title = ''
        subtitle = ''
        mach = np.nan
        q = np.nan
        cref = np.nan
        bref = np.nan
        sref = np.nan
        apress = AeroPressure(
            subcase, title, subtitle,
            mach, q,
            cref, bref, sref,
            nodes,
            cp_pressure[:, 0], cp_pressure[:, 1])
        return apress

    def get_element_pressure(self, nid_to_eid_map: dict[int, list[int]]):
        """
                           AERODYNAMIC PRES.       AERODYNAMIC
        GRID   LABEL          COEFFICIENTS           PRESSURES             ELEMENT
        156     LS            6.035214E-01         9.052821E-01            900014
        """
        element_pressures = defaultdict(list)
        nids = self.nodes
        nnids = len(nids)
        pressure = self.pressure
        assert pressure.shape == (nnids,)
        for nid, cp in zip(nids, self.cp):
            # print(nid, cp, q)
            eids = nid_to_eid_map[nid]
            for eid in eids:
                element_pressures[eid].append(cp)
        return dict(element_pressures)

    def __repr__(self) -> str:
        nnids = len(self.nodes)
        msg = (
            'AeroPressure:\n'
            f' - nodes ({nnids})\n'
            f' - cp ({nnids})\n'
            f' - pressure ({nnids})'
        )
        return msg

class AeroForce:
    def __init__(self, subcase: int, title: str, subtitle: str,
                 mach: float, q: float,
                 cref: float, bref: float, sref: float,
                 nodes: np.ndarray,
                 force: np.ndarray):
        self.subcase = subcase
        self.title = title
        self.subtitle = subtitle
        self.mach = mach
        self.q = q
        self.cref = cref
        self.bref = bref
        self.sref = sref

        self.nodes = nodes
        self.force = force

    @classmethod
    def from_f06(self, subcase: int,
                 nodes: np.ndarray,
                 force: np.ndarray):
        title = ''
        subtitle = ''
        mach = np.nan
        q = np.nan
        cref = np.nan
        bref = np.nan
        sref = np.nan
        aforce = AeroForce(
            subcase, title, subtitle,
            mach, q,
            cref, bref, sref,
            nodes, force)
        return aforce

    def __repr__(self) -> str:
        nnids = len(self.nodes)
        msg = (
            'AeroForce:\n'
            f' - nodes ({nnids})\n'
            f' - force ({nnids}, 6)'
        )
        return msg


class TrimResults:
    def __init__(self):
        self.metadata = {}

        # aero_pressure[subcase] = AeroPressure(...)
        self.aero_pressure: dict[int, AeroPressure] = {}

        # aero_force[subcase] = AeroForce(...)
        self.aero_force: dict[int, AeroForce] = {}
        self.structural_monitor_loads: dict[int, MonitorLoads] = {}

        # controller_state = {'ALPHA': 2.0}
        self.controller_state: dict[int, ControllerState] = {}

        # trim_variables[name] = [idi, Type, trim_status, ux, ux_unit]
        self.trim_variables: dict[int, TrimVariables] = {}  # TODO: not supported

    def __repr__(self) -> str:
        msg = (
            'TrimResults:'
        )
        if len(self.aero_force):
            keys = [str(case) for case in self.aero_force]
            msg += '\n  aero_force keys:\n - ' + '\n   - '.join(keys)
        else:
            msg += '\n  len(aero_force) = 0'

        if len(self.aero_pressure):
            keys = [str(case) for case in self.aero_pressure]
            msg += '\n  aero_pressure keys:\n - ' + '\n   - '.join(keys)
        else:
            msg += '\n  len(aero_pressure) = 0'
        return msg
