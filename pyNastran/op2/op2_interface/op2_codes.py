from typing import List, Tuple, Union, Optional
from pyNastran.op2.op2_interface.nx_tables import NX_ELEMENTS, NX_TABLE_CONTENT
from pyNastran.op2.op2_interface.msc_tables import MSC_ELEMENTS, MSC_TABLE_CONTENT

# strings
SORT1_TABLES = [b'OSTRMS1C', b'OSTNO1C', b'OES1X', b'OSTR1X',
                b'OESRMS2', b'OESNO2', b'OESXRMS1',
                b'OES1C', b'OSTR1C',
                'OES1C', 'OSTR1C', 'OESNLXR']
SORT2_TABLES = [b'OUGPSD2', b'OUGATO2', b'OESCP',
                b'OES2C', b'OSTR2C',
                b'OFMPF2M', b'OLMPF2M', b'OPMPF2M', b'OSMPF2M', b'OGPMPF2M',
                'OFMPF2M', 'OLMPF2M', 'OPMPF2M', 'OSMPF2M', 'OGPMPF2M',
                'OES2C', 'OSTR2C']
NO_SORT_METHOD = [b'QHHA']

NASA95_ELEMENTS = {
    #             OES       OEF
    1 : 'CROD  ',  # done   done
    2 : 'C.....',
    3 : 'CTUBE ',  # done   done
    4 : 'CSHEAR',  # done   done
    5 : 'CTWIST',
    6 : 'CTRIA1',
    7 : 'CTRBSC',
    8 : 'CTRPLT',
    9 : 'CTRMEM',
    10 : 'CONROD', # done   done
    11 : 'ELAS1',  # done   done
    12 : 'ELAS2',  # done   done
    13 : 'ELAS3',  # done   done
    14 : 'ELAS4',  # done   done
    15 : 'CQDPLT',
    16 : 'CQDMEM',
    17 : 'CTRIA2',
    18 : 'CQUAD2',
    19 : 'CQUAD1',
    20 : 'CDAMP1', #        done
    21 : 'CDAMP2', #        done
    22 : 'CDAMP3', #        done
    23 : 'CDAMP4', #        done
    24 : 'CVISC',  #        done
    25 : 'CMASS1',
    26 : 'CMASS2',
    27 : 'CMASS3',
    28 : 'CMASS4',
    29 : 'CONM1',
    30 : 'CONM2',
    31 : 'PLOTEL',
    32 : 'C.....',
    33 : 'C.....',
    34 : 'CBAR',    # done   done
    35 : 'CCONE',   #        done <--- CCONEAX
    36 : 'CTRIARG',
    37 : 'CTRAPRG',
    38 : 'CTORDRG',
    39 : 'CTETRA',
    40 : 'CWEDGE',
    41 : 'CHEXA1',
    42 : 'CHEXA2',
    43 : 'CFLUID2',
    44 : 'CFLUID3',
    45 : 'CFLUID4',
    46 : 'CFLMASS',
    47 : 'CAXIF2',
    48 : 'CAXIF3',
    49 : 'CAXIF4',
    50 : 'CSLOT3',
    51 : 'CSLOT4',
    52 : 'CHBDY',
    53 : 'CDUM1',
    54 : 'CDUM2',
    55 : 'CDUM3',
    56 : 'CDUM4',
    57 : 'CDUM5',
    58 : 'CDUM6',
    59 : 'CDUM7',
    60 : 'CDUM8',
    61 : 'CDUM9',
    62 : 'CQDMEM1',
    63 : 'CQDMEM2',
    64 : 'CQUAD4',
    65 : 'CIHEX1',
    66 : 'CIHEX2',
    67 : 'CIHEX3',
    68 : 'CQUADTS',
    69 : 'CTRIATS',
    70 : 'CTRIAAX',
    71 : 'CTRAPAX',
    72 : 'CAERO1',
    73 : 'CTRIM6',
    74 : 'CTRPLT1',
    75 : 'CTRSHL',
    76 : 'CFHEX1',
    77 : 'CFHEX2',
    78 : 'CFTETRA',
    79 : 'CFWEDGE',
    80 : 'CIS2D8',
    81 : 'CELBOW',
    82 : 'CFTUBE',
    83 : 'CTRIA3',  # done   done
}
ANALYSIS_CODE_MAP = {
    1 : "Statics",
    2 : "Normal modes or buckling (real eigenvalues)",
    3 : "Differential Stiffness 0 - obsolete",
    4 : "Differential Stiffness 1 - obsolete",
    5 : "Frequency",
    6 : "Transient",
    7 : "Pre-buckling",
    8 : "Post-buckling",
    9 : "Complex eigenvalues",
    10 : "Nonlinear statics",
    11 : "Geometric nonlinear statics",
}

DEVICE_CODE_MAP = {
    1 : "Print",
    2 : "Plot",
    3 : "Print and Plot",
    4 : "Punch",
    5 : "Print and Punch",
    6 : "Plot and Punch",
    7 : "Print, Plot, and Punch",
}

THERMAL_MAP = {
    0 : 'isHeatTransfer = False',
    1 : 'isHeatTransfer = True',
    2 : 'Scaled response spectra ABS',
    #3 : 'Scaled response spectra SRSS',
    4 : 'Scaled response spectra SRSS', # NRL???
    5 : 'Scaled response spectra NRLO',
    #6 :
    #7 :
    8 : 'Scaled response spectra NRL',
}

TABLE_CODE_MAP = {
    2 : "OPG - Load vector",
    3 : "OQG - SPC Force vector",
    #4: 'OEF - Element force/flux',
    #5: 'OES - Element Stress/Strain',
    6 : "LAMA - Eigenvalue summary",
    7 : "OUG - Eigenvector",
    8 : "none - Grid point singularity table (obsolete)",
    9 : 'OEIGS - Eigenvalue analysis summary',
    10 : "OUG - Velocity vector",
    11 : "OUG - Acceleration vector",
    12 : "OPG - Nonlinear force vector",
    13 : "OGPWG - Grid point weight generator",
    14 : "OUG - Eigenvector (solution set)",
    15 : "OUG - Displacement vector (solution set)",
    16 : "OUG - Velocity vector (solution set)",
    17 : "OUG - Acceleration vector (solution set)",
    18 : "OEE - Element strain energy",
    19 : "OGF - Grid point force balance",
    20 : "OES - Stresses at grid points (from the CURV module)",
    21 : "OES - Strain/curvature at grid points",
    22 : "OELOF1 - Element internal forces and moments",
    23 : "OELOP1 - Summation of element oriented forces on adjacent elements",
    24 : "OEP - Element pressures",
    25 : "OEF - Composite failure indicies",
    26 : "OGS - Grid point stresses (surface)",
    27 : "OGS - Grid point stresses (volume -- direct)",
    28 : "OGS - Grid point stresses (volume -- principal)",
    29 : "OGS - Element stress discontinuities (surface)",
    30 : "OGS - Element stress discontinuities (volume -- direct)",
    31 : "OGS - Element stress discontinuities (volume -- principal)",
    32 : "OGS - Grid point stress discontinuities (surface)",
    33 : "OGS - Grid point stress discontinuities (volume -- direct)",
    34 : "OGS - Grid point stress discontinuities (volume -- principal)",
    35 : "OGS - Grid point stress discontinuities (plane strain)",
    36 : "OEE - Element kinetic energy",
    37 : "OEE - Element energy loss",
    38 : "OMM - Max/Min summary",
    39 : "OQG - MPC Forces",
    40 : "OGPKE - Grip point kinetic energy",
}

GEOM_TABLES = { # no analysis code
    'CASECC', 'EDOM', 'VIEWTB', 'AXIC',
    'GEOM1', 'GEOM2', 'GEOM3', 'GEOM4', 'DYNAMIC', 'CONTACT',
    'GEOM1S', 'GEOM2S', 'GEOM3S', 'GEOM4S', 'DYNAMICS', 'CONTACTS',
    'GEOM1N',
    'PVT', 'PVT0',
    'EPT', 'MPT', 'DIT', 'EDT',
    'EPTS', 'MPTS', 'DITS', 'EDTS',
    'DBCOPT', 'DSCMCOL', 'DESCYC', 'R1TABRG',
}

def get_sort_method_from_table_name(table_name: bytes) -> int:
    """helper method"""
    if table_name in SORT1_TABLES:
        sort_method = 1
    elif table_name in SORT2_TABLES:
        sort_method = 2
    elif table_name in NO_SORT_METHOD:
        table_name_str = table_name.decode('utf8')
        table_num = table_name_str[-1]
        sort_method = -1
        #raise ValueError('%r is not a table' % table_name_str)
    else:
        table_name_str = table_name.decode('utf8')
        table_num = table_name_str[-1]
        try:
            sort_method = int(table_num)
        except ValueError:
            print(f'error determining sort_method: table_name={table_name_str} type={type(table_name)}')
            raise
    return sort_method


class Op2Codes:
    def __init__(self):
        pass

    def set_table_type(self) -> None:
        if self.is_msc:
            self.element_mapper = MSC_ELEMENTS
        elif self.is_nasa95:
            self.element_mapper = NASA95_ELEMENTS
        else:  # default
            self.element_mapper = NX_ELEMENTS

    def get_element_type(self, elem_code: int) -> str:
        self.set_table_type()
        try:
            etype = self.element_mapper[elem_code]
        except TypeError:
            print('elem_code=%r' % elem_code)
            raise
        return etype

    def print_table_code(self, table_code: int) -> str:
        #table_code_content = table_code % 1000
        #data_format = table_code / 1000
        msg = ''
        #msg += 'table_code_content=%s data_format=%s\n' %(table_code_content, data_format)

        table = get_table_from_table_code(table_code, self.table_name_str, is_msc=self.is_msc)
        if self.is_msc:
            msg += 'n=%s msc table=%s-%s' % (self.n, self.table_name, table)
        else:
            msg += 'n=%s nx table=%s-%s' % (self.n, self.table_name, table)
        return msg

    def approach_code_str(self, approach_code: int) -> str:
        """TODO: not done

        The approach code is the type of solution:
          1 : statics
          2 : modes
        ...

        This is not the same thing as SOL 101.  SOL 101 (linear statics)
        and SOL 144 (aero-statics) are both statics.

        """
        return ''

    def code_information(self, include_time: bool=True) -> str:
        """
        prints the general table information
        DMAP - page 60-63

        """
        device_code = self.device_code
        #analysis_code = self.analysis_code
        #table_code = self.table_code
        sort_code = self.sort_code

        format_code = None
        if hasattr(self, 'format_code'):
            format_code = self.format_code

        s_code = None
        if hasattr(self, 's_code'):
            s_code = self.s_code

        thermal = None
        if hasattr(self, 'thermal'):
            thermal = self.thermal

        s_word = ''
        stress_word = ''
        if hasattr(self, 'stress_bits'):
            if self.is_stress:
                stress_word = 'Stress'
            else:
                stress_word = 'Strain'
            s_word = get_scode_word(s_code, self.stress_bits)

        element_type = None
        if hasattr(self, 'element_type'):
            element_type = self.element_type

        format_word = get_format_word(format_code)

        if self.sort_bits[0] == 0:
            sort_word1 = 'Real'
        else:
            sort_word1 = 'Real/Imaginary'
        if self.sort_bits[1] == 0:
            sort_word2 = 'Sort1'
        else:
            sort_word2 = 'Sort2'
        if self.sort_bits[2] == 0:
            sort_word3 = 'Sorted Responses'
        else:
            sort_word3 = 'Random Responses'

        #if   self.sort_code==0: sortWord = 'Real'
        #elif self.sort_code==1: sortWord = 'Real/Imaginary'
        #elif self.sort_code==2: sortWord = 'Random Responses'
        #else:
            #sortWord = '???'
            #msg = 'unsupported sort_code:  sort_code=%s\n' %(sort_code)
            #print msg
            #raise RuntimeError(msg)

        try:
            thermal_word = THERMAL_MAP[thermal]
        except KeyError:
            thermal_word = '???'
            #msg = 'unsupported thermal:  thermal=%s\n' %(thermal)
            #raise ValueError(msg)

        analysis = '???'
        analysis_code = None
        is_geom_table = self.table_name_str in GEOM_TABLES
        if not is_geom_table:
            if hasattr(self, 'analysis_code'):
                analysis_code = self.analysis_code
                try:
                    analysis = ANALYSIS_CODE_MAP[analysis_code]
                except KeyError:
                    pass
            else:
                self.log.warning('%s has no analysis code' % self.table_name_str)
                raise  RuntimeError('%s has no analysis code' % self.table_name_str)

        try:
            device = DEVICE_CODE_MAP[self.device_code]
        except KeyError:
            device = '???'

        force_flux = self.get_force_flux(thermal)
        disp_temp = self.get_disp_temp(thermal)

        self_table_code = None
        table_code = None
        table = None
        if not is_geom_table:
            self_table_code = self.table_code
            table_code, table = self.get_table_code_name(disp_temp, force_flux, stress_word)

        msg = '--Table3Data--\n\n'
        msg += "  device_code   = %-3s %s\n" % (self.device_code, device)
        msg += "  analysis_code = %-3s %s\n" % (analysis_code, analysis)
        msg += "  table_code    = %-3s %s-%s\n" % (self_table_code, self.table_name_str, table)
        msg += "  format_code   = %-3s %s\n" % (format_code, format_word)

        msg += (
            f'  sort_method   = {self.sort_method}\n'
            f'  sort_code     = {self.sort_code}\n'
        )
        msg += "    sort_bits   = (%s, %s, %s)\n" % tuple(self.sort_bits)
        msg += "    data_format = %-3s %s\n" % (self.sort_bits[0], sort_word1)
        msg += "    sort_type   = %-3s %s\n" % (self.sort_bits[1], sort_word2)
        msg += "    is_random   = %-3s %s\n" % (self.sort_bits[2], sort_word3)

        random_code = self.random_code if hasattr(self, 'random_code') else 0
        msg += "  random_code   = %-3s\n" % (random_code)

        if element_type is not None:
            if isinstance(element_type, str):
                etype = element_type
            else:
                etype = self.get_element_type(element_type)
            msg += "  element_type  = %-3s %s\n" % (element_type, etype)

        if s_word:  # stress code
            msg += "  s_code        = %-3s %s\n" % (s_code, s_word)
        if thermal is not None:
            msg += "  thermal       = %-3s %s\n" % (thermal, thermal_word)
            if hasattr(self, 'thermal_bits'):
                msg += "  thermal_bits  = %s\n" % str(self.thermal_bits)

        if hasattr(self, 'num_wide'):
            msg += "  num_wide      = %-3s\n" % self.num_wide
        if hasattr(self, 'isubcase'):
            msg += "  isubcase      = %-3s\n" % self.isubcase
        else:
            msg += "  ID            = %-3s\n" % self.ID

        dt_names = [
            'dt', 'time', 'mode', 'eign', 'cycle', 'mode2',
            'freq', 'lsdvmn', 'eigr', 'eigi', 'lftsfq']
        for name in dt_names:
            if hasattr(self, name):
                dvalue = getattr(self, name)
                msg += "  %-6s        = %s\n" % (name, dvalue)


        if self.is_msc:
            msg += '  MSC Nastran\n'
        elif self.is_nasa95:
            msg += '  NASA 95 Nastran\n'
        else:
            msg += '  NX Nastran\n'
        #print msg
        if hasattr(self, 'format_code'):
            assert isinstance(self.format_code, int), type(self.format_code)
        return msg

    def get_force_flux(self, thermal: Optional[int]=None) -> str:
        """is this a force or flux table"""
        if thermal == 0:
            force_flux = 'Force'
        elif thermal == 1:
            force_flux = 'Flux'
        else:
            force_flux = f'Force (or Flux); thermal={thermal!r}'
        return force_flux

    def get_disp_temp(self, thermal: int) -> str:
        """is displacement (static) or thermal analysis being done"""
        if thermal == 0:
            disp_temp = 'Displacement'
        elif thermal == 1:
            disp_temp = 'Temperature'
        #elif thermal is None:
            #raise RuntimeError('thermal_code is not specified; thermal_code=None')
        else:
            disp_temp = f'Displacement/Temperature; thermal={thermal!r}'
        return disp_temp

    def get_table_code_name(self, disp_temp: str='', force_flux: str='',
                            stress_word: str='') -> Tuple[int, str]:
        """gets the name of the active table"""
        table = '???'
        table_code = self.table_code
        table_code = _adjust_table_code(self.table_code)

        if table_code == 1:
            table = f'OUG - {disp_temp} vector/scalar'
        elif table_code == 4:
            table = f'OEF - Element {force_flux}'
        elif table_code == 5:
            table = f'OES - Element {stress_word}'
        else:
            #try:
            table = get_table_from_table_code(table_code, self.table_name_str, is_msc=self.is_msc)
            #except KeyError:
                #table = '%s - Unknown' % self.table_name

        return table_code, table

    @property
    def table_name_str(self) -> str:
        """
        Converts the table_name from bytes/str to a str

        Returns
        -------
        table_name_str : str
            the table name as a string

        ..note :: Refers to bytes/str in the Python 3 sense.
        """
        table_name = self.table_name
        if isinstance(table_name, bytes):
            table_name = self.table_name.decode(self._encoding)
        return table_name

    #----
    def is_thermal(self) -> Union[bool, str]:
        """is this result thermal solution?"""
        if self.thermal == 0:
            return False
        elif self.thermal == 1:
            return True
        return '???'

    #----
    # format_code 3
    def is_magnitude_phase(self) -> bool:
        if self.format_code == 3:
            return True
        return False

    def is_sort1_new(self) -> bool: # pragma: no cover
        #is_sort1_table = self.is_sort1
        table_name = self.table_name_str
        if table_name in SORT1_TABLES:
            is_sort1_table = True
        elif table_name in SORT2_TABLES:
            is_sort1_table = False
        else:
            try:
                sort_method, is_real, is_random = self._table_specs()
                is_sort1_table = (sort_method == 1)
            except AssertionError:
                try:
                    is_sort1_table = int(table_name[-1]) == 1
                except ValueError:
                    raise ValueError(f'is this SORT1/2?  table_name={table_name!r}')
        return is_sort1_table

    @property
    def is_sort1(self) -> bool:
        #is_sort1_table = self.is_sort1
        try:
            sort_method, is_real, is_random = self._table_specs()
            is_sort1_table = (sort_method == 1)
        except AssertionError:
            table_name = self.table_name_str
            if table_name in SORT1_TABLES:
                is_sort1_table = True
            elif table_name in SORT2_TABLES:
                is_sort1_table = False
            else:
                try:
                    is_sort1_table = int(table_name[-1]) == 1
                except ValueError:
                    raise ValueError(f'is this SORT1/2?  table_name={table_name!r}')
        return is_sort1_table

    @property
    def is_sort2(self) -> bool:
        #return not self.is_sort1
        try:
            sort_method, is_real, is_random = self._table_specs()
            is_sort2_table = (sort_method == 2)
        except AssertionError:
            table_name = self.table_name_str
            if table_name in SORT2_TABLES:
                is_sort2_table = True
            elif table_name in SORT1_TABLES:
                is_sort2_table = False
            else:
                try:
                    is_sort2_table = int(table_name[-1]) == 2
                except ValueError:
                    raise ValueError(f'is this SORT1/2?  table_name={table_name!r}')
        return is_sort2_table

    def update_t_code(self) -> None:
        """
        Value Sort type Data format Random
        ===== ========= =========== ======
        0     SORT1     Real        No
        1     SORT1     Complex     No
        2     SORT2     Real        No
        3     SORT2     Complex     No
        4     SORT1     Real        Yes
        5     SORT2     Complex?    Yes
        6     SORT2     Real        Yes

        table_code%1000 = function3()

        SPCForce = table_code % 1000 (function 3)

        """
        is_complex, is_sort2, is_random = self.sort_bits
        map_sort_bits = {
            # is_complex, is_sort2, is_random
            (0, 0, 0) : 0,
            (1, 0, 0) : 1,

            (0, 1, 0) : 2,
            (1, 1, 0) : 3,

            # random
            (0, 0, 1) : 4,
            (1, 1, 1) : 5, # not 100%
            (0, 1, 1) : 6,
        }
        unused_t_code = map_sort_bits[(is_complex, is_sort2, is_random)]

    def _table_specs(self) -> Tuple[int, bool, bool]:
        """
        +-------+-----------+-------------+----------+
        | Value | Sort Type | Data Format | Random ? |
        +-------+-----------+-------------+----------+
        |   0   |   SORT1   |    Real     |   No     |
        +-------+-----------+-------------+----------+
        |   1   |   SORT1   |    Complex  |   No     |
        +-------+-----------+-------------+----------+
        |   2   |   SORT2   |    Real     |   No     |
        +-------+-----------+-------------+----------+
        |   3   |   SORT2   |    Complex  |   No     |
        +-------+-----------+-------------+----------+
        |   4   |   SORT1   |    Real     |   Yes    |
        +-------+-----------+-------------+----------+
        |   5   |   SORT1   |    Real     |   ???    |
        +-------+-----------+-------------+----------+
        |   6   |   SORT2   |    Real     |   Yes    |
        +-------+-----------+-------------+----------+
        |   7   |    ???    |    ???      |   ???    |
        +-------+-----------+-------------+----------+

        +-----+-------------+---------+
        | Bit |     0       |    1    |
        +-----+-------------+---------+
        |  0  | Not Random  | Random  |
        |  1  | SORT1       | SORT2   |
        |  2  | Real        | Complex |
        +-----+-------------+---------+
        """
        #tcode = self.table_code // 1000
        table_code = self.tCode
        tcode = self.sort_code
        sort_code = tcode
        #if self.table_name_str == 'OQGPSD2':
            #print(self.code_information())
            #print('table_name=%s tCode=%s sort_code=%s self.sort_bits=%s' % (self.table_name_str, self.tCode, sort_code, self.sort_bits))
        assert sort_code in [0, 1, 2, 3, 4, 5, 6], 'sort_code=%s\n%s' % (sort_code, self.code_information())
        try:
            sort_method, is_real, is_random = determine_sort_bits_meaning(table_code, sort_code, self.sort_bits)
        except AssertionError:
            #print(self.code_information())
            raise
        return sort_method, is_real, is_random

    #----
    # sort_code
    # disabled 11/2015
    #def isSortedResponse(self):
        #if self.sort_bits[0] == 0:
            #return True
        #return False

    # disabled 11/2015
    #def isRandomResponse(self):
        #return not self.isSortedResponse()

    #----
    # combos
    #def isRealOrRandom(self):  # been broken for a long time
        #return self.isReal() or self.isRandom()

    #def isRealImaginaryOrMagnitudePhase(self):  # been broken for a long time
        #return self.is_real_imaginary or self.MagnitudePhase()

    #----

SCODE_MAP = {
    # word, stress_bits_expected
    0 : ('Coordinate Element - Stress Max Shear (Octahedral)',       [0, 0, 0, 0, 0]),
    14: ('Coordinate Element - Strain Fiber Max Shear (Octahedral)', [0, 1, 1, 1, 0]),

    1: ('Coordinate Element - Stress von Mises',                         [0, 0, 0, 0, 1]),
    10: ('Coordinate Element - Strain Curvature Max Shear (Octahedral)', [0, 1, 0, 1, 0]),

    11: ('Coordinate Element - Strain Curvature von Mises', [0, 1, 0, 1, 1]),
    15: ('Coordinate Element - Strain Fiber von Mises',     [0, 1, 1, 1, 1]),

    16: ('Coordinate Material - Stress Max Shear (Octahedral)', [1, 0, 0, 0, 0]),
    17: ('Coordinate Material - Stress von Mises',              [1, 0, 0, 0, 1]),

    26: ('Coordinate Material - Strain Curvature Max Shear', [1, 1, 0, 1, 0]),
    30: ('Coordinate Material - Strain Fiber Max Shear (Octahedral)', [1, 1, 1, 1, 0]),

    27: ('Coordinate Material - Strain Curvature von Mises', (1, 1, 0, 1, 1)),
    31: ('Coordinate Material - Strain Fiber von Mises', (1, 1, 1, 1, 1)),
}

def get_scode_word(s_code: int, stress_bits: List[int]) -> str:
    try:
        s_word, stress_bits_expected = SCODE_MAP[s_code]
    except KeyError:
        #s_word = 'Stress or Strain - UNDEFINED'
        s_word = '???'
    return s_word

def determine_sort_bits_meaning(table_code: int, sort_code: int,
                                sort_bits: Tuple[int, int, int]) -> Tuple[int, bool, bool]:
    """
    Value Sort type Data format Random
    ===== ========= =========== ======
    0     SORT1     Real        No
    1     SORT1     Complex     No
    2     SORT2     Real        No
    3     SORT2     Complex     No
    4     SORT1     Real        Yes
    5     SORT2     ???         Yes
    6     SORT2     Real        Yes

    table_code%1000 = function3()

    SPCForce = table_code % 1000 (function 3)

    """
    sort_method = 1
    is_real = True
    is_random = False
    # old
    #if sort_code in [2, 3, 5, 6]:
        #sort_method = 2
    #if sort_code in [1, 3]:
        #is_real = False
    #if sort_code in [4, 5, 6]:
        #is_random = True

    # new
    if sort_code in [2, 3, 5, 6]:
        sort_method = 2
    if sort_code in [1, 3]:
        is_real = False
    if sort_code in [4, 5, 6]:
        is_random = True

    try:
        if is_random:
            assert sort_bits[0] == 1, 'should be RANDOM; sort_bits=%s; sort_code=%s' % (sort_bits, sort_code)
        else:
            assert sort_bits[0] == 0, 'should be NOT RANDOM; sort_bits=%s; sort_code=%s' % (sort_bits, sort_code)

        if sort_method == 1:
            assert sort_bits[1] == 0, 'should be SORT1; sort_bits=%s; sort_code=%s' % (sort_bits, sort_code)
        else:
            assert sort_bits[1] == 1, 'should be SORT2; sort_bits=%s; sort_code=%s' % (sort_bits, sort_code)

        if is_real:
            assert sort_bits[2] == 0, 'should be REAL; sort_bits=%s; sort_code=%s; table_code=%s table_code%%1000=%s' % (sort_bits, sort_code, table_code, table_code % 1000)
        else:
            assert sort_bits[2] == 1, 'should be IMAG; sort_bits=%s; sort_code=%s; table_code=%s table_code%%1000=%s' % (sort_bits, sort_code, table_code, table_code % 1000)
    except AssertionError:
        #print('sort_method=%r; is_real=%r is_random=%r' % (sort_method, is_real, is_random))
        raise
    return sort_method, is_real, is_random

def get_table_from_table_code(table_code: int, table_name: str, is_msc: bool=True) -> str:
    """translates that a key of say 1 is the 'OUG - Displacement vector' table"""
    try:
        if is_msc:
            table = MSC_TABLE_CONTENT[table_code]
        else:
            table = NX_TABLE_CONTENT[table_code]
    except:
        print(f'count not determine the table description for {table_name}')
        raise

    #table = TABLE_CODE_MAP[table_code]
    return table

def _adjust_table_code(table_code: int) -> int:
    """
    table code handler for alternate table code form:
       501 -> 1
       510 -> 10
       511 -> 11
       601 -> 1
       etc.
    """
    if table_code in [501, 510, 511]:
        table_code -= 500
    elif table_code in [601, 610, 611]:
        table_code -= 600
    elif table_code in [701, 710, 711]:
        table_code -= 700
    elif table_code in [801, 810, 811]:
        table_code -= 800
    elif table_code in [901, 910, 911]:
        table_code -= 900
    return table_code

def get_format_word(format_code: int) -> str:
    format_word = '???'
    if format_code == 1:
        format_word = 'Real'
    elif format_code == 2:
        format_word = 'Real/Imaginary'
    elif format_code == 3:
        format_word = 'Magnitude/Phase'
    else:
        format_word = '\n%18s1 - Real\n' % ''
        format_word += '%18s2 - Real/Imaginary\n' % ''
        format_word += '%18s3 - Magnitude/Phase\n' % ''
        #msg = 'unsupported format_code:  format_code=%s\n' % format_code
        #raise InvalidFormatCodeError(msg)
    return format_word
