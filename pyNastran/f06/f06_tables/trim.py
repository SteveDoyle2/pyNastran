from collections import defaultdict
import io
import numpy as np
from pyNastran.utils.numpy_utils import integer_types
#TrimVariable = tuple[int, str, str, float, str]
ControllerState = dict[str, float]
#TrimVariables = dict[str, TrimVariable]

class Statics:
    def __init__(self, title: str, subtitle: str, label: str):
        self.nonlinear_factor = None
        self.is_complex = False
        self.is_real = False

        assert isinstance(title, str), title
        assert isinstance(subtitle, str), subtitle
        assert isinstance(label, str), label
        self.title = title
        self.subtitle = subtitle
        self.label = label


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


class AeroPressure(Statics):
    def __init__(self, subcase: int, title: str, subtitle: str,
                 mach: float, q: float,
                 cref: float, bref: float, sref: float,
                 nodes: np.ndarray,
                 cp: np.ndarray, pressure: np.ndarray):
        label = ''
        super().__init__(title, subtitle, label)
        self.subcase = subcase
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

    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  aero_pressure[{self.subcase}]:\n'

        # if short:
        msg += f'    nodes:     n={len(self.nodes)}\n'
        msg += f'    pressure:  {str(self.pressure.shape)}\n'
        msg += f'    cp:        {str(self.cp.shape)}\n'
        # else:
        #     for (name, typei, status, units), data in zip(self.name_type_status_units, self.data):
        #         msg += f'   {name:<8} {typei:<8} {status:<8}: {data} {units}\n'
        return msg

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

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        # msg = (
        #     '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
        #     '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
        #     '\n')
        msg = '     Node    PRESSURE     CP\n'
        for (nid, press, cp) in zip(self.nodes, self.pressure, self.cp):
            msg += f' {nid:>10d}  {press:13.6e}  {cp:13.6e}\n'
        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1

    def __repr__(self) -> str:
        nnids = len(self.nodes)
        msg = (
            'AeroPressure:\n'
            f' - nodes ({nnids})\n'
            f' - cp ({nnids})\n'
            f' - pressure ({nnids})'
        )
        return msg

class AeroForce(Statics):
    def __init__(self, subcase: int, title: str, subtitle: str,
                 mach: float, q: float,
                 cref: float, bref: float, sref: float,
                 nodes: np.ndarray,
                 force: np.ndarray,
                 force_label: np.ndarray):
        label = ''
        super().__init__(title, subtitle, label)
        self.subcase = subcase
        self.mach = mach
        self.q = q
        self.cref = cref
        self.bref = bref
        self.sref = sref

        self.nodes = nodes
        self.force = force
        self.force_label = force_label

    @classmethod
    def from_f06(self, subcase: int,
                 nodes: np.ndarray,
                 force: np.ndarray,
                 force_label: np.ndarray):
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
            nodes, force, force_label)
        return aforce

    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  aero_force[{self.subcase}]:\n'

        # if short:
        msg += f'    nodes:        n={len(self.nodes)}\n'
        msg += f'    force:        {str(self.force.shape)}\n'
        msg += f'    force_label:  {str(self.force_label.shape)}\n'
        # else:
        #     for (name, typei, status, units), data in zip(self.name_type_status_units, self.data):
        #         msg += f'   {name:<8} {typei:<8} {status:<8}: {data} {units}\n'
        return msg

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        # msg = (
        #     '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
        #     '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
        #     '\n')
        msg = '     Node          FX          FY          FZ          MX          MY          MZ     CP\n'
        for (nid, force, label) in zip(self.nodes, self.force, self.force_label):
            # print(nid, force, label)
            fx, fy, fz, mx, my, mz = force
            msg += f' {nid:>10d}  {fx:13.6e}  {fy:13.6e}  {fz:13.6e}  '\
                   f'{mx:13.6e}  {my:13.6e}  {mz:13.6e}  {label:10s}\n'
        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1

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
class TrimVariables(Statics):
    """

                                   A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S
                             CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC
                                               MACH = 7.890000E-01                Q = 1.500000E+00
                             CHORD = 1.3110E+02           SPAN = 2.5564E+03            AREA = 7.3400E+05

              TRIM ALGORITHM USED: LINEAR TRIM SOLUTION WITHOUT REDUNDANT CONTROL SURFACES.


                                                          AEROELASTIC TRIM VARIABLES

                                      ID     LABEL                 TYPE        TRIM STATUS      VALUE OF UX

                                             INTERCEPT          RIGID BODY           FIXED      1.000000E+00
                                     500     ANGLEA             RIGID BODY            FREE      1.047265E-01  RADIANS
                                     501     PITCH              RIGID BODY           FIXED      0.000000E+00  NONDIMEN. RATE
                                     502     URDD3              RIGID BODY           FIXED      2.500000E+00  LOAD FACTOR
                                     503     URDD5              RIGID BODY           FIXED      0.000000E+00  LOAD FACTOR
                                  110000     TFLAP         CONTROL SURFACE            FREE     -4.541439E-01  RADIANS


                                               CONTROL SURFACE POSITION AND HINGE MOMENT RESULTS

                                ACTIVE LIMITS ARE FLAGGED WITH AN (A),  VIOLATED LIMITS ARE FLAGGED WITH A (V).

                                                    POSITION                                         HINGE MOMENT
              CONTROL SURFACE      LOWER LIMIT       VALUE         UPPER LIMIT        LOWER LIMIT       VALUE         UPPER LIMIT
                  TFLAP          -1.570796E+00   -4.541439E-01    1.570796E+00            N/A        1.672770E+06         N/A
    1    N2A STATIC AEROELASTIC AND FLUTTER MODEL                                 MARCH  10, 2025  SIMCENTER NASTRAN  8/25/23   PAGE    24

    """
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 name_type_status_units: np.ndarray, data: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        super().__init__(title, subtitle, label)
        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

        self.name_type_status_units = name_type_status_units
        self.data = data

    def __eq__(self, other) -> bool:
        return True

    def get_stats(self, short: bool=False) -> str:
        msg = ''
        msg += f'  variables[{self.subcase}]:\n'
        # self.name_type_status = name_type_status
        # self.data = data

        if short:
            nnames = len(self.name_type_status_units)
            msg += f'    name_type_status_units: {self.name_type_status_units.tolist()}; n={nnames}\n'
            msg += f'    data.shape = {str(self.data.shape)}\n'
        else:
            for (name, typei, status, units), data in zip(self.name_type_status_units, self.data):
                msg += f'   {name:<8} {typei:<8} {status:<8}: {data} {units}\n'
        return msg

    def print_f06(self) -> str:
        f06_file = io.StringIO()
        header = []
        page_stamp = '%d'
        out = self.write_f06(f06_file, header, page_stamp)
        return f06_file.getvalue()[:-1]

    def write_f06(self, f06_file, header: list[str], page_stamp: str, page_num: int=1,
                  is_mag_phase: bool=False, is_sort1: bool=True):
        f06_file.write(''.join(header))
        msg = (
            '    TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL                  ELASTIC/RIGID\n'
            '                                   UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED    UNSPLINED  SPLINED\n'
            '\n')

        msg = (
            ''
            '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S\n'
            '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC\n'
            f'                                           MACH = {self.mach:13.6E}                Q = {self.q:13.6E}\n'
            f'                         CHORD = {self.chord:7.4E}           SPAN = {self.span:7.4E}            AREA = {self.sref:7.4E}\n'
            '\n'
            '          TRIM ALGORITHM USED: LINEAR TRIM SOLUTION WITHOUT REDUNDANT CONTROL SURFACES.                                           \n'
            '\n'
            '\n'
            '                                                      AEROELASTIC TRIM VARIABLES'
            '\n'
            '                                  ID     LABEL                 TYPE        TRIM STATUS      VALUE OF UX\n'
            '\n'
            # '                                         INTERCEPT          RIGID BODY           FIXED      1.000000E+00\n'
            # '                                 500     ANGLEA             RIGID BODY            FREE      1.047265E-01  RADIANS         \n'
            # '                                 501     PITCH              RIGID BODY           FIXED      0.000000E+00  NONDIMEN. RATE  \n'
            # '                                 502     URDD3              RIGID BODY           FIXED      2.500000E+00  LOAD FACTOR     \n'
            # '                                 503     URDD5              RIGID BODY           FIXED      0.000000E+00  LOAD FACTOR     \n'
            # '                              110000     TFLAP         CONTROL SURFACE            FREE     -4.541439E-01  RADIANS         \n'
            # '\n'
            # '\n'
            # '                                           CONTROL SURFACE POSITION AND HINGE MOMENT RESULTS\n'
            # '\n'
            # '                            ACTIVE LIMITS ARE FLAGGED WITH AN (A),  VIOLATED LIMITS ARE FLAGGED WITH A (V).\n'
            # '\n'
            # '                                                POSITION                                         HINGE MOMENT\n'
            # '          CONTROL SURFACE      LOWER LIMIT       VALUE         UPPER LIMIT        LOWER LIMIT       VALUE         UPPER LIMIT\n'
            # '              TFLAP          -1.570796E+00   -4.541439E-01    1.570796E+00            N/A        1.672770E+06         N/A     \n'
        )
        # for (name, typei, status) in self.name_type_status:
        for (name, typei, status, units), data in zip(self.name_type_status_units.tolist(), self.data):
            # assert name == 'INTERCEPT', f'name={name!r}'
            # assert typei.upper() == 'RIGID BODY', f'typei={typei!r}'
            # assert status.upper() == 'FIXED', f'status={status!r}'
            # print(f'name={name}; data={data}')
            # assert len(data) == 1, data
            value = data
            extra = ''
            if units in 'RADIANS':
                extra = f'  {np.degrees(value):13.6E}  DEGREES'
            elif units == 'NONDIMEN. RATE':
                extra = f'  {np.degrees(value):13.6E}  DEGREES/SEC'

            msg += f'                                        {name.upper():<10}    {typei.upper():>16}      {status.upper():>10}     {value:13.6E}  {units.upper():<14s}{extra}\n'

        f06_file.write(msg)
        f06_file.write(page_stamp % page_num)
        return page_num + 1


class TrimDerivatives(Statics):
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 names: np.ndarray, derivatives: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        super().__init__(title, subtitle, label)

        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

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


class ControlSurfacePostiionHingeMoment(Statics):
    def __init__(self, mach: float, q: float,
                 chord: float, span: float, sref: float,
                 names: np.ndarray,
                 trim_values: np.ndarray,
                 hinge_moments: np.ndarray,
                 subcase: int=1, title: str='', subtitle: str='', label: str=''):
        super().__init__(title, subtitle, label)

        self.mach = mach
        self.q = q
        self.chord = chord
        self.span = span
        self.sref = sref
        self.subcase = subcase
        assert isinstance(self.subcase, integer_types), self.subcase

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
