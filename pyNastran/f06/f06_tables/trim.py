from collections import defaultdict
import numpy as np
from pyNastran.utils.numpy_utils import integer_types
#TrimVariable = tuple[int, str, str, float, str]
ControllerState = dict[str, float]
#TrimVariables = dict[str, TrimVariable]


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
                 force: np.ndarray,
                 label: np.ndarray):
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
        self.label = label

    @classmethod
    def from_f06(self, subcase: int,
                 nodes: np.ndarray,
                 force: np.ndarray,
                 label: np.ndarray):
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
            nodes, force, label)
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


# trim_variables = TrimVariables(
#     mach, q, cref, bref, sref,
#     name_type_status, trim_values_array, derivatives_array,
#     subcase=subcase_id, title=title,
#     subtitle=subtitle, label=label)
class TrimVariables:
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 name_type_status: np.ndarray, data: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        self.nonlinear_factor = None
        self.is_complex = False
        self.is_real = False

        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

        self.title = title
        self.subtitle = subtitle
        self.label = label
        assert isinstance(title, str), title
        assert isinstance(subtitle, str), subtitle
        assert isinstance(label, str), label

        self.name_type_status = name_type_status
        self.data = data

    def __eq__(self, other) -> bool:
        return True

    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  variables[{self.subcase}]:\n'
        # self.name_type_status = name_type_status
        # self.data = data

        if short:
            nnames = len(self.name_type_status)
            msg += f'    name_type_status: {self.name_type_status.tolist()}; n={nnames}\n'
            msg += f'    data.shape = {str(self.data.shape)}\n'
        else:
            for (name, typei, status), data in zip(self.name_type_status, self.data):
                msg += f'   {name:<8} {typei:<8} {status:<8}: {data}\n'
        return msg

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        msg = (
            '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
            '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
            '\n')

        for (name, typei, status), data in zip(self.name_type_status, self.data):
            msg += f'   {name:<8} {typei:<8} {status:<8}: {data}\n'
        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1


class TrimDerivatives:
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 names: np.ndarray, derivatives: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        self.nonlinear_factor = None
        self.is_complex = False
        self.is_real = False

        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

        self.title = title
        self.subtitle = subtitle
        self.label = label
        assert isinstance(title, str), title
        assert isinstance(subtitle, str), subtitle
        assert isinstance(label, str), label

        self.names = names
        self.data = derivatives

    def __eq__(self, other) -> bool:
        return True
    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  derivatives[{self.subcase}]:\n'
        coeffs = ['Cx', 'Cy', 'Cz', 'Cmx', 'Cmy', 'Cmz']
        headers = ['rigid_unsplined', 'rigid_splined', 'elastic_unsplined', 'elastic_splined',
                   'inertial_restrained', 'inertial_unrestrained']

        if short:
            nnames = len(self.names)
            msg += f'    names   = {self.names.tolist()}; n={nnames}\n'
            msg += f'    headers = {headers}\n'
            msg += f'    coeffs  = {coeffs}\n'
            msg += f'    data.shape = ({nnames}, 6, 6)\n'
        else:
            for name, derivs in zip(self.names, self.data):
                msg += f'    {name}:\n'
                for coeff, line in zip(coeffs, derivs):
                    msg += f'     {coeff}: {line}\n'
        return msg

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        msg = (
            '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
            '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
            '\n')

        coeffs = ['Cx', 'Cy', 'Cz', 'Cmx', 'Cmy', 'Cmz']
        for name, derivs in zip(self.names, self.data):
            # msg += f'    {name}:\n'
            name_str = name
            for coeff, line in zip(coeffs, derivs):
                # msg += f'     {coeff}: {line}\n'
                rigid_unsplined, rigid_splined, elastic_unsplined, elastic_splined, inertial_restrained, inertial_unrestrained = line
                try:
                    e2r_unsplined = elastic_unsplined / rigid_unsplined
                except FloatingPointError:
                    e2r_unsplined = np.nan
                try:
                    e2r_splined = elastic_splined / rigid_splined
                except FloatingPointError:
                    e2r_splined = np.nan
                msg += (
                    f'    {name_str:<14}   {coeff:9s}   '
                    f'{rigid_unsplined:>13.6E}  {rigid_splined:>13.6E}  '
                    f'{elastic_unsplined:>13.6E}  {elastic_splined:>13.6E}  '
                    f'{inertial_restrained:>13.6E}  {inertial_unrestrained:>13.6e}  '
                    f'{e2r_unsplined:>13.6E}  {e2r_splined:>13.6E}  \n')
                name_str = ''
            msg += '\n'
        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1

    def __repr__(self):
        return f'TrimDerivatives(subcase={self.subcase}, title, subtitle)'


class ControlSurfacePostiionHingeMoment:
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 names: np.ndarray,
                 trim_values: np.ndarray,
                 hinge_moments: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        self.nonlinear_factor = None
        self.is_complex = False
        self.is_real = False

        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

        self.title = title
        self.subtitle = subtitle
        self.label = label
        assert isinstance(title, str), title
        assert isinstance(subtitle, str), subtitle
        assert isinstance(label, str), label

        self.names = names
        self.trim_values = trim_values
        self.data = hinge_moments

    def __eq__(self, other) -> bool:
        return True
    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  control_surface_position_hinge_moment[{self.subcase:d}]:\n'
        coeffs = ['Cx', 'Cy', 'Cz', 'Cmx', 'Cmy', 'Cmz']
        nnames = len(self.names)

        if short or 1:
            nnames = len(self.names)
            msg += f'    names       = {self.names.tolist()}; n={nnames:d}\n'
            msg += f'    trim_values = {self.trim_values.tolist()}\n'
        else:
            for name, derivs in zip(self.names, self.data):
                msg += f'    {name}:\n'
                for coeff, line in zip(coeffs, derivs):
                    msg += f'     {coeff}: {line}\n'
        return msg

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        # msg = (
        #     '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
        #     '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
        #     '\n')
        msg = 'name, trim_value, position_lower, position_upper, hinge_moment, hinge_moment_lower, hinge_moment_upper\n'
        for name, trim_value, values in zip(self.names, self.trim_values, self.data):
            # values = [-1.571e+00  1.571e+00  1.673e+06 - 1.000e+10  1.000e+10]
            (position_lower, position_upper, hinge_moment,
             hinge_moment_lower, hinge_moment_upper) = values
            msg += (
                f'    {name:<14}   {trim_value:>13.6E}   '
                f'{position_lower:>13.6E}  {position_upper:>13.6E}  {hinge_moment:>13.6E}  '
                f'{hinge_moment_lower:>13.6E}  {hinge_moment_upper:>13.6E}\n')
            msg += '\n'
        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1

    def __repr__(self) -> str:
        out = f'ControlSurfacePostiionHingeMoment(subcase={self.subcase!r}, '\
              f'title={self.title!r}, subtitle={self.subtitle!r})'
        return out
