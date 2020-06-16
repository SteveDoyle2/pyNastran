import warnings
import numpy as np

class Responses:
    """Defines SOL 200 responses"""
    def __init__(self):
        self.convergence_data = None
        self.desvars = None
        self.dscmcol = None
        self.weight_response = None
        self.displacement_response = None
        self.stress_response = None
        self.strain_response = None
        self.force_response = None
        self.composite_stress_response = None
        self.composite_strain_response = None
        self.flutter_response = None
        self.fractional_mass_response = None

    def get_stats(self, short=False) -> str:
        objects = [
            self.convergence_data,
            self.desvars,
            self.dscmcol,
            self.displacement_response,
            self.weight_response,
            self.stress_response,
            self.strain_response,
            self.force_response,
            self.composite_stress_response,
            self.composite_strain_response,
            self.flutter_response,
            self.fractional_mass_response,
        ]
        msg = []
        for obj in objects:
            if obj is not None:
                msg += obj.get_stats(short=short) + '\n'
        return msg

    def get_table_types(self):
        tables = [
            'convergence_data', 'desvars', 'dscmcol',
            'displacement_response', 'weight_response',
            'stress_response', 'strain_response', 'force_response',
            'composite_stress_response', 'composite_strain_response', 'flutter_response',
            'fractional_mass_response',
        ]
        return ['responses.' + table for table in tables
                #if getattr(self, table) is not None
                ]

class WeightResponse:
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []

    def add_from_op2(self, out, log):
        """
        Weight Response

        1 IRID         I Internal response identification number
        2 RID          I External response identification number
        3 TYPE(C)      I Response type
        4 LABEL(2) CHAR4 Label
        6 REGION       I Region identifier
        7 SCID         I Subcase identification number
        8 UNDEF(2)     I Not used
        10 SEID        I Superelement identification number or ALL
        11 UNDEF(2)    I Not used

        13 UNDEF       I Not used
        14 TYFLG       I Flag to indicate how response is referenced
        15 SEID        I Superelement identificaiton number

        ---> 3i 8s 7i 3i

        #                             -----  WEIGHT RESPONSE  -----
        # ---------------------------------------------------------------------------------
        #  INTERNAL  DRESP1  RESPONSE  ROW  COLUMN  LOWER     INPUT      OUTPUT     UPPER
        #     ID       ID     LABEL     ID    ID    BOUND     VALUE       VALUE     BOUND
        # ---------------------------------------------------------------------------------
        #       1       1    WEIGHT     3     3       N/A   2.9861E+05  2.9852E+05   N/A

        # ?  ?     ?  LABEL?      ?  ?     ROW_ID? COL_ID? ?  ?  ?          ?          ?  ?
        #(1, 1,    1, 'WEIGHT  ', 0, 1011, 3,      3,      0, 0, 0,         0,         0, 0)
        #(1, 1000, 1, 'W       ', 0, 1,    3,      3,      0, 0, 0,         0,         0, 0)
        #
        # per dev forum; 538976288 is probably just '    '
        #(1, 1,    1, 'WEIGHT  ', 0, 1,    3,      3,      0, 0, 538976288, 538976288, 0, 0)
        """
        # F:\work\pyNastran\examples\Dropbox\move_tpl\mbcgen.op2
        # (1,  15, 1, 'W       ', -1,    1, 3, 3, 0, 0, 0, 0, 0, 10)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\i2002.op2
        # (20, 15, 1, 'W       ', -1,    1, 3, 3, 0, 0, 0, 0, 0, 20)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\edr2n.op2
        # (1, 201, 1, 'WEIGHT  ', -1,    0, 3, 3, 0, 0, 0, 0, 0, 100)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\ss200m2.op2
        # (1,   1, 1, 'W       ', 0, 30001, 3, 3, 0, 0, 0, 0, 0, 300)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\mcso43.op2
        # (1,  10, 1, 'W       ', 0,     1, 3, 3, 0, 0, 0, 0, 1, 0)
        #--------------------------
        # common
        #
        # per the R1TAB DMAP page:
        #   all indicies are downshift by 1
        #   indices above out[3] are off by +2 because of the 2 field response_label
        internal_id = out[0]
        dresp_id = out[1]
        response_type = out[2]
        response_label = out[3].strip()
        # -1 for 2 field wide response_label
        region = out[4]
        subcase = out[5]
        type_flag = out[12]  # no meaning per MSC DMAP 2005
        seid = out[13]

        #--------------------------------------------------
        #row_id = out[4]

        # these should be blank?
        row_id = out[6]
        column_id = out[7]
        seid_weight = out[8]

        if np.abs(out[8:-1]).sum() != 0.0:
            msg = 'WEIGHT response sum error; out=%s 8=%s' % (out, out[8:-1])
            log.warning(msg)
        assert seid == out[-1]
        #assert out[-1] in [0, 1, 2, 3, 4, 5, 10, 20, 100, 200, 300], out
        #dunno_8 = out[8]
        #dunno_9 = out[9]
        #dunno_10 = out[10]
        #dunno_11 = out[11]
        #dunno_12 = out[12]
        #msg = ('WEIGHT - response_type=%r response_label=%r row_id=%r column_id=%r '
            #'6=%r 7=%r 8=%r 9=%r 10=%r 11=%r 12=%r seid=%r' % (
            #response_type, response_label, row_id, column_id,
            #dunno_6, dunno_7, dunno_8, dunno_9, dunno_10, dunno_11, dunno_12, seid))
        #out = unpack(self._endian + b'iii 8s iiff f fffff', data)
        #msg = 'WEIGHT - label=%r region=%s subcase=%s row_id=%r column_id=%r' % (
            #response_label, region, subcase, row_id, column_id)
        self.append(internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    row_id, column_id)
        #print(msg)
        #self.log.debug(msg)

    def append(self, internal_id, dresp_id, response_label, region,
               subcase, type_flag, seid,
               row_id, column_id):
        self.internal_id.append(internal_id)
        #self.response_label.append(response_label)
        #self.subcase.append(subcase)
        #self.mode.append(mode)
        #self.mach.append(mach)
        #self.velocity.append(velocity)
        #self.density.append(density)
        #self.flutter_id.append(flutter_id)
        #self.subcase.append(subcase)
        self._n += 1

    def __repr__(self):
        msg = 'responses.WeightResponse()\n'
        msg += '  n=%s\n' % self.n
        #msg += '  velocity=%s\n' % (velocity)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.weight_response (%s)\n' % (self.n)
        return self.__repr__() + '\n'


class PropertyResponse:
    """common class for StressResponse, StrainResponse, and ForceResponse"""
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.response_label = []
        self.subcase = []
        self.item_code = []
        self.pid = []

    def append(self, internal_id, dresp_id, response_label, region, subcase, type_flag, seid,
               item_code, pid):
        assert isinstance(response_label, str), response_label
        self.internal_id.append(internal_id)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.item_code.append(item_code)
        self.pid.append(pid)
        self._n += 1

    def __repr__(self):
        name = self.__class__.__name__
        msg = 'responses.%s()\n' % name
        msg += '  internal_id = %s\n' % np.array(self.internal_id)
        msg += '  response_label = %s\n' % np.array(self.response_label)
        msg += '  subcase = %s\n' % np.array(self.subcase)
        msg += '  item_code = %s\n' % np.array(self.item_code)
        msg += '  pid = %s\n' % np.array(self.pid)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)\n' % (self.name, self.n)
        return self.__repr__() + '\n'

class FractionalMassResponse:
    name = 'fractional_mass'
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.dresp_id = []
        self.region = []
        self.response_label = []
        self.subcase = []
        self.type_flag = []
        self.seid = []

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)\n' % (self.name, self.n)
        return self.__repr__() + '\n'

    def append(self, internal_id, dresp_id, response_label, region,
               subcase, type_flag, seid):
        self.internal_id.append(internal_id)
        self.dresp_id.append(dresp_id)
        self.region.append(region)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.type_flag.append(type_flag)
        self.seid.append(seid)
        self._n += 1

    def __repr__(self):
        msg = 'FractionalMassResponse()\n'
        msg += '  response_label=%s\n' % np.array(self.response_label)
        msg += '  internal_id=%s\n' % np.array(self.internal_id)
        msg += '  region=%s\n' % np.array(self.region)
        msg += '  subcase=%s\n' % np.array(self.subcase)
        msg += '  type_flag=%s\n' % np.array(self.type_flag)
        msg += '  seid=%s\n' % np.array(self.seid)
        return msg


class DisplacementResponse:
    name = 'displacement'
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.dresp_id = []
        self.region = []
        self.response_label = []
        self.subcase = []
        self.nid = []
        self.component = []
        self.type_flag = []
        self.seid = []

    def append(self, internal_id, dresp_id, response_label, region,
               subcase, type_flag, seid,
               nid, component):
        self.internal_id.append(internal_id)
        self.dresp_id.append(dresp_id)
        self.region.append(region)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.nid.append(nid)
        self.component.append(component)
        self.type_flag.append(type_flag)
        self.seid.append(seid)
        self._n += 1

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)\n' % (self.name, self.n)
        return self.__repr__() + '\n'

    def __repr__(self):
        msg = 'DisplacementResponse()\n'
        msg += '  response_label=%s\n' % np.array(self.response_label)
        msg += '  nid=%s\n' % np.array(self.nid)
        msg += '  component=%s\n' % np.array(self.component)
        msg += '  internal_id=%s\n' % np.array(self.internal_id)
        msg += '  subcase=%s\n' % np.array(self.subcase)
        msg += '  type_flag=%s\n' % np.array(self.type_flag)
        msg += '  seid=%s\n' % np.array(self.seid)
        return msg

class ForceResponse(PropertyResponse):
    name = 'force'
class StressResponse(PropertyResponse):
    name = 'stress'
class StrainResponse(PropertyResponse):
    name = 'strain'

class FlutterResponse:
    name = 'flutter'
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.response_label = []
        self.subcase = []
        self.mode = []
        self.mach = []
        self.velocity = []
        self.density = []
        self.flutter_id = []
        self.subcase = []

    def append(self, internal_id, dresp_id, response_label, region, subcase, type_flag, seid,
               mode, mach, velocity, density, flutter_id):
        self.internal_id.append(internal_id)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.mode.append(mode)
        self.mach.append(mach)
        self.velocity.append(velocity)
        self.density.append(density)
        self.flutter_id.append(flutter_id)
        self._n += 1

    def __repr__(self):
        msg = 'FlutterResponse()\n'
        msg += '  velocity=%s\n' % np.array(self.velocity)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)\n' % (self.name, self.n)
        return self.__repr__() + '\n'


class DSCMCOL:
    """
    '                                   -----   IDENTIFICATION OF COLUMNS IN THE DESIGN SENSITIVITY  -----'
    '                                   -----     MATRIX THAT ARE ASSOCIATED WITH DRESP1  ENTRIES    -----'
    ''
    ''
    '             -----  WEIGHT/VOLUME RESPONSES  -----'
    '          ------------------------------------------'
    '            COLUMN         DRESP1         RESPONSE'
    '              NO.         ENTRY ID          TYPE  '
    '          ------------------------------------------'
    '                1               1        WEIGHT  '
    ''
    ''
    '             -----  STATICS RESPONSES  -----'
    '          ------------------------------------------------------------------------------------------------------------------------'
    '            COLUMN         DRESP1         RESPONSE        GRID/ELM          VIEW         COMPONENT           SUB             PLY  '
    '              NO.         ENTRY ID          TYPE             ID            ELM ID            NO.             CASE             NO. '
    #'            COLUMN         DRESP1         RESPONSE        GRID/ELM          VIEW         COMPONENT           SUB  '
    #'              NO.         ENTRY ID          TYPE             ID            ELM ID            NO.             CASE '
    '          ------------------------------------------------------------------------------------------------------------------------'
    '                2               2        STRESS                 1                               5               1'
    '                3               2        STRESS                 3                               5               1'
    '                4           10501        DISP                 100                               1               1'
    ''
    ''
    #'                                                 ---- RETAINED DRESP2 RESPONSES ----'
    #''
    #'     ----------------------------------------------------------------------------------------------------------'
    #'         INTERNAL      DRESP2      RESPONSE     EQUATION         LOWER                             UPPER   '
    #'            ID           ID         LABEL          ID            BOUND            VALUE            BOUND   '
    #'     ----------------------------------------------------------------------------------------------------------'
    #'                3          105     DISPMAG             3      -1.0000E+20       9.4677E-05       2.5000E-04'
    '                                   -----   IDENTIFICATION OF COLUMNS IN THE DESIGN SENSITIVITY  -----'
    '                                   -----     MATRIX THAT ARE ASSOCIATED WITH DRESP2  ENTRIES    -----'
    ''
    ''
    '          ----------------------------------------------------------'
    '            COLUMN         DRESP2           SUB             FREQ/ '
    '              NO.         ENTRY ID          CASE            TIME  '
    '          ----------------------------------------------------------'
    '               96             103               0         0.00000'
    '               97             104               0         0.00000'
    '              100             105               1'


    """
    def __init__(self, responses):
        """internal_response_id = iresponse + 1 = column in DSCM2"""
        self.responses = responses

    #def write_f06(self):
        #msg = [
            #'                                   -----   IDENTIFICATION OF COLUMNS IN THE DESIGN SENSITIVITY  -----'
            #'                                   -----     MATRIX THAT ARE ASSOCIATED WITH DRESP1  ENTRIES    -----'
            #''
            #''
            #'             -----  WEIGHT/VOLUME RESPONSES  -----'
            #'          ------------------------------------------'
            #'            COLUMN         DRESP1         RESPONSE'
            #'              NO.         ENTRY ID          TYPE  '
            #'          ------------------------------------------'
            #'                1               1        WEIGHT  '
            #''
            #''
            #'             -----  STATICS RESPONSES  -----'
            #'          ------------------------------------------------------------------------------------------------------------------------'
            #'            COLUMN         DRESP1         RESPONSE        GRID/ELM          VIEW         COMPONENT           SUB  '
            #'              NO.         ENTRY ID          TYPE             ID            ELM ID            NO.             CASE '
            #'          ------------------------------------------------------------------------------------------------------------------------'
            #'                2               2        STRESS                 1                               5               1'
            #'                3               2        STRESS                 3                               5               1'
            #''
            #''
            #'             -----  EIGENVALUE RESPONSES  -----'
            #'          --------------------------------------------------------------------------'
            #'            COLUMN         DRESP1         RESPONSE          MODE            SUB  '
            #'              NO.         ENTRY ID          TYPE             NO.            CASE '
            #'          --------------------------------------------------------------------------'
            #'                2               4        FREQ                   1               1'
        #]
        #return msg

    def get_responses_by_group(self):
        response_groups_order = [
            'weight_volume', 'static', 'eigenvalue',
            '???', 'psd',
            '2',
        ]
        responses_groups = {key: [] for key in response_groups_order}


        response_name_to_group = {
            'weight': 'weight_volume',
            'volume': 'weight_volume',

            'static stress': 'static',
            'static displacement': 'static',
            'static strain': 'static',
            'composite failure': 'static',
            'composite strain': 'static',

            'psd displacement': 'psd',
            'psd acceleration': 'psd',

            'normal modes': 'eigenvalue',

            'equivalent radiated power': '???',
            'frequency response displacement': '???',
            'frequency response stress?': '???',
            'ceig': '???',
        }
        response_name_to_f06_response_type = {
            # weight/volume
            'weight' : 'WEIGHT',
            'volume': 'VOLUME',
            'normal modes': 'FREQ',

            # static
            'static stress': 'STRESS',
            'static strain': 'STRAIN',
            'static displacement': 'DISP',
            'composite failure': 'CFAILURE',
            'composite strain': 'CSTRAIN',

            # psd
            'psd displacement': 'DISP',
            'psd acceleration': 'ACCE',

            # ???
            'frequency response displacement': '???',
            'frequency response stress?': '???',
            'ceig': '???',
        }

        for i, respi in self.responses.items():
            response_number = respi['response_number']
            if 'name' in respi:
                name = respi['name']
                group = response_name_to_group[name]
                responses_groups[group].append(i)
                assert response_number == 1, respi
            else:
                assert response_number == 2, respi
                responses_groups['2'].append(i)
        #print(f'responses_groups = {responses_groups}')

        msg = ''
        #'                                   -----   IDENTIFICATION OF COLUMNS IN THE DESIGN SENSITIVITY  -----'
        #'                                   -----     MATRIX THAT ARE ASSOCIATED WITH DRESP1  ENTRIES    -----'
        #''
        #''
        missing_keys = set(list(responses_groups.keys())) - set(response_groups_order)
        assert len(missing_keys) == 0, missing_keys
        for group_name in response_groups_order:
            ids = responses_groups[group_name]
            if len(ids) == 0:
                continue

            if group_name == 'weight_volume':
                msg += (
                    '             -----  WEIGHT/VOLUME RESPONSES  -----\n'
                    '          ------------------------------------------\n'
                    '            COLUMN         DRESP1         RESPONSE\n'
                    '              NO.         ENTRY ID          TYPE  \n'
                    '          ------------------------------------------\n')
                for i in ids:
                    respi = self.responses[i]
                    external_id = respi['external_response_id']
                    name = respi['name']
                    response_type = response_name_to_f06_response_type[name]
                    msg += f'         {i+1:8d}        {external_id:8d}        {response_type:8s}\n'
            elif group_name == 'static':
                msg += self._write_static(ids, response_name_to_f06_response_type)
            elif group_name == 'eigenvalue':
                msg += (
                    '             -----  EIGENVALUE RESPONSES  -----\n'
                    '          --------------------------------------------------------------------------\n'
                    '            COLUMN         DRESP1         RESPONSE          MODE            SUB  \n'
                    '              NO.         ENTRY ID          TYPE             NO.            CASE \n'
                    '          --------------------------------------------------------------------------\n'
                )
                for i in ids:
                    respi = self.responses[i]
                    external_id = respi['external_response_id']
                    name = respi['name']
                    mode_num = respi['mode_num']
                    subcase = respi['subcase']
                    response_type = response_name_to_f06_response_type[name]
                    msg += f'         {i+1:8d}        {external_id:8d}    {response_type:8s}            {mode_num:8d}        {subcase:8d}\n'
            elif group_name == 'psd':
                msg += self._write_psd(ids, response_name_to_f06_response_type)
            elif group_name == '2':
                msg += self._write_dresp2(ids)
            else:
                warnings.warn(f'skipping DSCMCOL group_name={group_name}')
                for i in ids:
                    respi = self.responses[i]
                    warnings.warn(str(respi))
                continue
            msg += '\n\n'
        str(msg)

    def _write_static(self, ids, response_name_to_f06_response_type):
        msg = ''
        is_composite = False
        for i in ids:
            respi = self.responses[i]
            name = respi['name']
            if name in ['composite failure', 'composite strain']:
                is_composite = True
                break

        msg += (
            '             -----  STATICS RESPONSES  -----\n'
            '          ------------------------------------------------------------------------------------------------------------------------\n'
        )
        if is_composite:
            msg += (
                '            COLUMN         DRESP1         RESPONSE        GRID/ELM          VIEW         COMPONENT           SUB             PLY  \n'
                '              NO.         ENTRY ID          TYPE             ID            ELM ID            NO.             CASE             NO. \n'
            )
        else:
            msg += (
                '            COLUMN         DRESP1         RESPONSE        GRID/ELM          VIEW         COMPONENT           SUB  \n'
                '              NO.         ENTRY ID          TYPE             ID            ELM ID            NO.             CASE \n'
                #'          ------------------------------------------------------------------------------------------------------------------------\n'
                #'                2               2        STRESS                 1                               5               1'
            )
        msg += '          ------------------------------------------------------------------------------------------------------------------------\n'

        for i in ids:
            respi = self.responses[i]
            external_id = respi['external_response_id']
            name = respi['name']
            subcase = respi['subcase']
            comp = respi['component']
            if name in ['static stress', 'static strain']:
                eid = respi['eid']
                response_type = response_name_to_f06_response_type[name]
                msg += f'         {i+1:8d}        {external_id:8d}        {response_type:8s}        {eid:8d}                        {comp:8d}        {subcase:8d}\n'
            elif name in ['static displacement']:
                grid = respi['grid']
                response_type = response_name_to_f06_response_type[name]
                msg += f'         {i+1:8d}        {external_id:8d}        {response_type:8s}        {grid:8d}                        {comp:8d}        {subcase:8d}\n'
            elif name in ['composite failure', 'composite strain']:
                eid = respi['eid']
                ply = respi['ply']
                response_type = response_name_to_f06_response_type[name]
                msg += f'         {i+1:8d}        {external_id:8d}        {response_type:8s}        {eid:8d}                0       {comp:8d}        {subcase:8d}        {ply:8d}\n'
                #msg += f'         {i+1:8d}        {external_id:8d}      {response_type:8s}          {eid:8d}                        {comp:8d}        {subcase:8d}\n'
            else:
                raise NotImplementedError(respi)
            return msg

    def _write_psd(self, ids, response_name_to_f06_response_type):
        msg = ''
        #is_composite = False
        #for i in ids:
            #respi = self.responses[i]
            #name = respi['name']
            #if name in ['composite failure', 'composite strain']:
                #is_composite = True
                #break

        msg += (
            '             -----  PSD RESPONSES  -----\n'
        )
        msg += (
            '          ------------------------------------------------------------------------------------------------------------------------\n'
            '            COLUMN         DRESP1         RESPONSE        GRID/ELM        RANDPS         COMPONENT           SUB             PLY  \n'
            '              NO.         ENTRY ID          TYPE             ID            ID                NO.             CASE             NO. \n'
            '          ------------------------------------------------------------------------------------------------------------------------\n'
        )

        for i in ids:
            respi = self.responses[i]
            external_id = respi['external_response_id']
            name = respi['name']
            subcase = respi['subcase']
            comp = respi['component']
            if name in ['psd displacement', 'psd acceleration']:
                grid = respi['grid']
                randps = respi['randps']
                response_type = response_name_to_f06_response_type[name]
                msg += f'         {i+1:8d}        {external_id:8d}        {response_type:8s}        {grid:8d}        {randps:8d}        {comp:8d}        {subcase:8d}\n'
            else:
                raise RuntimeError(respi)
        return msg

    def _write_dresp2(self, ids):
        msg = (
            '                                   -----   IDENTIFICATION OF COLUMNS IN THE DESIGN SENSITIVITY  -----\n'
            '                                   -----     MATRIX THAT ARE ASSOCIATED WITH DRESP2  ENTRIES    -----\n'
            '\n'
            '\n'
            '          ----------------------------------------------------------\n'
            '            COLUMN         DRESP2           SUB             FREQ/ \n'
            '              NO.         ENTRY ID          CASE            TIME  \n'
            '          ----------------------------------------------------------\n'
        )
        #'               96             103               0         0.00000'
        #'               97             104               0         0.00000'
        #'              100             105               1'
        for i in ids:
            #{'iresponse': 95, 'response_number': 2, 'internal_response_id': 1, 'external_response_id': 103,
             #'subcase': 0, 'dflag': 0, 'freq': 0.0, 'seid': 0}
            respi = self.responses[i]
            external_id = respi['external_response_id']
            freq = respi['freq']
            subcase = respi['subcase']
            dflag = respi['dflag']
            if dflag:
                msg += f'         {i+1:8d}        {external_id:8d}        {subcase:8d}\n'
            else:
                msg += f'         {i+1:8d}        {external_id:8d}        {subcase:8d}        {freq:8f}\n'
        return msg

    def __repr__(self):
        self.get_responses_by_group()
        nresponses = len(self.responses)
        internal_ids = np.zeros(nresponses, dtype='int32')
        external_ids = np.zeros(nresponses, dtype='int32')
        response_types = -np.full(nresponses, -1, dtype='int32')
        names_per_response = {}
        #internal_ids = [resp['internal_id'] for resp in self.responses.values()]
        #external_ids = [resp['external_id'] for resp in self.responses.values()]
        #response_types = [resp['response_type'] for resp in self.responses.values()]
        is_dresp2 = False
        for i, respi in self.responses.items():
            internal_ids[i] = respi['internal_response_id']
            external_ids[i] = respi['external_response_id']
            try:
                response_type = respi['response_type']
                response_types[i] = response_type
                name = respi['name']
                if (response_type, name) not in names_per_response:
                    names_per_response[(response_type, name)] = list(respi.keys())
            except KeyError:
                is_dresp2 = True
        msg = 'DSCMCOL()\n'
        msg += f'  nresponses: {nresponses}\n'
        msg += f'  internal_ids: {internal_ids}\n'
        msg += f'  external_ids: {external_ids}\n'
        msg += f'  response_types: {response_types}\n'
        if names_per_response:
            msg += '  response_type names:\n'
            for (response_type, name), names in names_per_response.items():
                names.remove('internal_response_id')
                names.remove('external_response_id')
                names2 = str(names).replace("'", '')
                msg += f'    {response_type:2d}: name={name!r}; names={names2}\n'
        return msg

    def get_stats(self, short=False):
        if short:
            return f'responses.dscmcol ({len(self.responses)})\n'
        return self.__repr__() + '\n'


class Desvars:
    def __init__(self, desvars):
        #internal_id, desvar_id, label, lower, upper, delxv, dunno = desvar
        ndesvars = len(desvars)
        self.internal_id = np.zeros(ndesvars, dtype='int32')
        self.desvar_id = np.zeros(ndesvars, dtype='int32')
        self.label = np.zeros(ndesvars, dtype='|U8')
        self.lower = np.zeros(ndesvars, dtype='float32')
        self.upper = np.zeros(ndesvars, dtype='float32')
        self.delxv = np.zeros(ndesvars, dtype='float32')
        self.dunno = np.zeros(ndesvars, dtype='float32')

        encoding = 'ascii'
        for i, (internal_id, desvar_id, label, lower, upper, delxv, dunno) in enumerate(desvars):
            #print((internal_id, desvar_id, label, lower, upper, delxv, dunno))
            self.internal_id[i] = internal_id
            self.desvar_id[i] = desvar_id
            self.label[i] = label.decode(encoding).strip()
            self.lower[i] = lower
            self.upper[i] = upper
            self.delxv[i] = delxv
            self.dunno[i] = dunno

    def __repr__(self):
        msg = 'Desvars()\n'
        #msg += '  shape=(%s, %s)\n' % (self.n, self.ndesign_variables)
        msg += '  internal_id = %s\n' % self.internal_id
        msg += '  desvar_id = %s\n' % self.desvar_id
        msg += '  label = %s\n' % self.label
        msg += '  lower = %s\n' % self.lower
        msg += '  upper = %s\n' % self.upper
        msg += '  delxv = %s\n' % self.delxv
        msg += '  dunno = %s\n' % self.dunno
        return msg

    def get_stats(self, short=False):
        if short:
            return f'responses.desvars ({len(self.internal_id)})\n'
        return self.__repr__() + '\n'


class Convergence:
    def __init__(self, ndesign_variables):
        self.n = 1
        self._n = 0
        self.is_built = False
        self.ndesign_variables = ndesign_variables
        self.design_iter = []
        self.iconvergence = []  #       1-soft, 2-hard
        self.conv_result = []   # 0-no, 1-soft, 2-hard
        self.obj_initial = []
        self.obj_final = []
        self.constraint_max = []
        self.row_constraint_max = []
        self.desvar_values = []

    def get_convergence(self):
        iconvergence = self.iconvergence[-1] # hard, soft (2)
        conv_result = self.conv_result[-1] # no, hard, soft (3)
        # 2*3 = 6
        if (conv_result, iconvergence) == ('hard', 'hard'):  # sure about this
            convergence = 'HARD'
        elif (conv_result, iconvergence) == ('no', 'hard'):  # pretty sure about this...
            convergence = 'MAX DESIGN CYCLES'
        elif (conv_result, iconvergence) == ('no', 'soft'):  # is this possible???
            convergence = 'MAX DESIGN CYCLES'

        # are these possible???
        # we'll just assume SOFT convergence
        if (conv_result, iconvergence) == ('hard', 'soft'):  # pretty sure about this...
            convergence = 'SOFT'
        elif (conv_result, iconvergence) == ('soft', 'hard'): # is this possible???
            convergence = 'SOFT'
        elif (conv_result, iconvergence) == ('soft', 'soft'): # is this possible???
            convergence = 'SOFT'
        else:  # pragma: no cover
            msg = 'conv_result=%r iconvergence=%r' % (conv_result, iconvergence)
            raise NotImplementedError(msg)
        return convergence, (conv_result, iconvergence)

    def append(self, design_iter, iconvergence, conv_result, obj_initial, obj_final,
               constraint_max, row_constraint_max, desvar_values):
        if not self.is_built:
            self.design_iter = np.zeros(self.n, dtype='int32')
            self.iconvergence = np.zeros(self.n, dtype=object)
            self.conv_result = np.zeros(self.n, dtype=object)
            self.obj_initial = np.zeros(self.n, dtype='float32')
            self.obj_final = np.zeros(self.n, dtype='float32')
            self.constraint_max = np.zeros(self.n, dtype='float32')
            self.row_constraint_max = np.zeros(self.n, dtype='int32')
            self.desvar_values = np.zeros((self.n, self.ndesign_variables), dtype='float32')
            self.is_built = True

        n = self._n
        self.design_iter[n] = design_iter
        self.iconvergence[n] = iconvergence
        self.conv_result[n] = conv_result
        self.obj_initial[n] = obj_initial
        self.obj_final[n] = obj_final
        self.constraint_max[n] = constraint_max
        self.row_constraint_max[n] = row_constraint_max
        self.desvar_values[n, :] = desvar_values
        self._n += 1

    def __repr__(self):
        msg = 'Convergence()\n'
        msg += '  design_iter = %s\n' % self.design_iter
        msg += '  icovergence = %s\n' % self.iconvergence
        msg += '  conv_result = %s\n' % self.conv_result
        msg += '  obj_initial = %s\n' % self.obj_initial
        msg += '  obj_final = %s\n' % self.obj_final
        msg += '  constraint_max = %s\n' % self.constraint_max
        msg += '  row_constraint_max = %s\n' % self.row_constraint_max
        msg += f'  desvar_values; shape=({self.n}, {self.ndesign_variables}):\n{self.desvar_values}'
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.convergence_data (%s, %s)\n' % (self.n, self.ndesign_variables)
        return self.__repr__() + '\n'
