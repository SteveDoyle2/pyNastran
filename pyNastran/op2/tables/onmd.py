from __future__ import annotations
from collections import defaultdict
import numpy as np
from typing import TYPE_CHECKING
#from pyNastran.op2.op2_interface.op2_reader import mapfmt
#from pyNastran.op2.op2_interface.op2_common import OP2Common
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class NormalizedMassDensity:
    def __init__(self, isubcase: int,
                 approach_code: int, analysis_code: int, device_code: int, num_wide: int,
                 dcycle: int, robj: float, rcon: float,
                 title: str, subtitle: str, label: str):
        self.isubcase = isubcase

        self.approach_code = approach_code
        self.analysis_code = analysis_code
        self.device_code = device_code
        self.dcycle = dcycle
        self.robj = robj
        self.rcon = rcon

        self.title = title
        self.subtitle = subtitle
        self.label = label

        ## number of words per entry in record
        self.num_wide = num_wide

        self.eids = np.array([])
        self.data = np.array([])

    def get_stats(self, key='', short=True):
        #key2 = f'[{key!r}]'
        if short:
            msg = f'NormalizedMassDensity(isubcase={self.isubcase}, dcycle={self.dcycle:g}, robj={self.robj:g}, rcon={self.rcon:g}); neids={len(self.eids)}\n'
            #msg = (f'NormalizedMassDensity{key2}: ref_point=%s mass=%g; '
                   #'[reference_point, M0, S, mass, cg, IS, IQ, Q]\n' % (
                       #self.reference_point, self.mass.max()))
        else:
            msg = f'NormalizedMassDensity(isubcase={self.isubcase}, dcycle={self.dcycle:g}, robj={self.robj:g}, rcon={self.rcon:g}); neids={len(self.eids)}\n'
            #msg = (
                #f'GridPointWeight{key2}:'
                #'  reference_point=%s\n'
                #'  mass=[%10g %10g %10g]\n'
                #'  cg  =[%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n\n'

                #'  IS  =[%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n\n'

                #'  IQ  =[%10g %10s %10s]\n'
                #'       [%10s %10g %10s]\n'
                #'       [%10s %10s %10g]\n\n'

                #'  Q  = [%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n'
                #'       [%10g %10g %10g]\n' % (
                    #self.reference_point, self.mass[0], self.mass[1], self.mass[2],
                    #self.cg[0, 0], self.cg[0, 1], self.cg[0, 2],
                    #self.cg[1, 0], self.cg[1, 1], self.cg[1, 2],
                    #self.cg[2, 0], self.cg[2, 1], self.cg[2, 2],

                    #self.IS[0, 0], self.IS[0, 1], self.IS[0, 2],
                    #self.IS[1, 0], self.IS[1, 1], self.IS[1, 2],
                    #self.IS[2, 0], self.IS[2, 1], self.IS[2, 2],

                    #self.IQ[0], '', '',
                    #'', self.IQ[1], '',
                    #'', '', self.IQ[2],

                    #self.Q[0, 0], self.Q[0, 1], self.Q[0, 2],
                    #self.Q[1, 0], self.Q[1, 1], self.Q[1, 2],
                    #self.Q[2, 0], self.Q[2, 1], self.Q[2, 2],
                    #)
            #)
        return msg

    #def __repr__(self) -> str:
        #return ''


class ONMD:
    def __init__(self, op2: OP2):
        self.op2 = op2
        self.obj = None

    def _read_onmd_3(self, data: bytes, ndata: int):
        """
        reads table 3 (the header table)

        Word Name Type Description
        1 ACODE(C)     I     Device code + 10*Approach code
        2 TCODE(C)     I     Table code 92
        3 UNDEF(3)     None
        6 DCYCLE       I     Design cycle number
        7 ROBJ         RS    Objective value
        8 RCON         RS    Critical constraint value
        9 UNDEF        None
        10 NUMWDE      I     Number of words per entry in DATA record (always 2)
        11 UNDEF(40)   None
        51 TITLE(32)   CHAR4 Title character string (TITLE)
        83 SUBTITL(32) CHAR4 Subtitle character string (SUBTITLE)
        115 LABEL(32)  CHAR4 LABEL character string (LABEL)
        """
        op2 = self.op2
        op2.to_nx('; found ONMD (normalized mass density) table')
        #self.log.info('OUG table 3')
        #self.show_data(data, types='ifs')
        #self.log.info('----------------------------')
        #self.show_ndata(400, types='if')
        #self._set_times_dtype()
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = True
        op2.is_table_2 = False
        assert data is not None
        op2.parse_approach_code(data)  # field 3
        op2.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', 'dcycle', 'robj',
            'rcon', 'num_wide', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']
        #op2.log.warning(f"approach_code={op2.approach_code}; analysis_code={op2.analysis_code}")

        #6 DCYCLE       I     Design cycle number
        #7 ROBJ         RS    Objective value
        #8 RCON         RS    Critical constraint value

        op2.dcycle = op2.add_data_parameter(data, 'dcycle', b'i', 7-1, False)
        op2.robj = op2.add_data_parameter(data, 'robj', b'f', 8-1, False)
        op2.rcon = op2.add_data_parameter(data, 'rcon', b'f', 9-1, False)

        ## number of words per entry in record
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)
        assert op2.num_wide == 2, op2.num_wide

        analysis_code = op2.analysis_code
        #if analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            #op2.setNullNonlinearFactor()
        #elif analysis_code == 2:  # real eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            #op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            ## float - C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mftank.op2
            ##self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'i', 7, False)  # nope...
            #op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False) # radians
            #op2.reader_oug.update_mode_cycle('mode_cycle')
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        ##elif analysis_code == 3: # differential stiffness
            ##self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            ##self.data_code['lsdvmn'] = self.lsdvmn
        ##elif analysis_code == 4: # differential stiffness
            ##self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        #elif analysis_code == 5:   # frequency
            ## frequency
            #op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        #elif analysis_code == 6:  # transient
            ## time step
            #op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['dt'])
        #elif analysis_code == 7:  # pre-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif analysis_code == 8:  # post-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        #elif analysis_code == 9:  # complex eigenvalues
            ## mode number
            #op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            #op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        #elif analysis_code == 10:  # nonlinear statics
            ## load step
            #op2.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lftsfq'])
        #elif analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
        #elif analysis_code == 13:
        # 11 Geometric nonlinear statics
        # 12 Response Simulation
        # 13 Response Simulation Random
        # 14 Response Simulation Shock Response Spectrum
        # 15 Response Simulation Quasi-Static Event
        if analysis_code == 16:  # topology optimization; not defined in DMAP
            pass
        else:  # pragma: no cover
            msg = f'invalid analysis_code...analysis_code={analysis_code}\ndata={op2.data_code}'
            raise RuntimeError(msg)

        #print self.code_information()
        #op2._fix_oug_format_code()
        #op2._parse_thermal_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()
        self.obj = NormalizedMassDensity(
            op2.isubcase,
            op2.approach_code, analysis_code, op2.device_code, op2.num_wide,
            op2.dcycle, op2.robj, op2.rcon,
            op2.title, op2.subtitle, op2.label)

        #op2._correct_eigenvalue()
        #print(op2.code_information())

    def _read_onmd_4(self, data: bytes, ndata: int) -> int:
        """
        Word Name Type Description
        1 EKEY  I Device code + 10*Element ID
        2 VALUE RS Scalar value for element
        """
        if not data:
            return ndata
        op2 = self.op2

        fdata = np.frombuffer(data, dtype=op2.fdtype8)
        idata = np.frombuffer(data, dtype=op2.idtype8)
        ndata = len(idata)
        #op2.log.warning(f'ndata={ndata}')
        eids = idata[::2] // 10
        data = fdata[1::2]
        #print(f'eids = {eids}')
        #print(f'data = {data}')
        self.obj.eids = eids
        self.obj.data = data

        key = op2.isubcase
        responses = op2.op2_results.responses
        if responses.normalized_mass_density is None:
            responses.normalized_mass_density = defaultdict(list)
        responses.normalized_mass_density[key].append(self.obj)
        return ndata
