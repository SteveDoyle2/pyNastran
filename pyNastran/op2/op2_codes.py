from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types

class Op2Codes(object):
    def __init__(self):
        pass

    def set_table_type(self):
        msc_elements, nx_elements = self._get_element_mappers()
        if self.is_msc:
            self.element_mapper = msc_elements
        else:
            self.element_mapper = nx_elements

    def get_element_type(self, elem_code):
        msc_elements, nx_elements = self._get_element_mappers()
        if self.is_msc:
            self.element_mapper = msc_elements
        else:
            self.element_mapper = nx_elements
        try:
            etype = self.element_mapper[elem_code]
        except TypeError:
            print('elem_code=%r' % elem_code)
            raise
        return etype

    def _get_element_mappers(self):
        msc_elements = {
            None: '',
            0: 'GRID',
            1: 'CROD',
            2: 'CBEAM',
            3: 'CTUBE',
            4: 'CSHEAR',
            5: 'FORMON12-FORCEi/MOMENTi follower',
            6: 'FORCE',
            7: 'PLOAD4',
            8: 'PLOADX1',
            9: 'PLOAD/PLOAD2',
            10: 'CONROD',
            11: 'CELAS1',
            12: 'CELAS2',
            13: 'CELAS3',
            14: 'CELAS4',
            15: 'AEROT3',
            16: 'AEROBEAM',
            17: 'CTRIAX',
            18: 'CQUADX',
            19: 'CQUAD1-old',
            20: 'CDAMP1',
            21: 'CDAMP2',
            22: 'CDAMP3',
            23: 'CDAMP4',
            24: 'CVISC',
            25: 'CMASS1',
            26: 'CMASS2',
            27: 'CMASS3',
            28: 'MASS4',
            29: 'CONM1',
            30: 'CONM2',
            31: 'PLOTEL',
            32: None,
            33: 'CQUAD4',
            34: 'CBAR-34',
            35: 'CCONEAX',
            36: 'CTRIARG-old',
            37: 'CTRAPRG-old',
            38: 'CGAP',
            39: 'CTETRA',
            40: 'CBUS1D',
            41: 'CHEXA1-old',
            42: 'CHEXA2-old',
            43: 'CFLUID2',
            44: 'CFLUID3',
            45: 'CFLUID4',
            46: 'FLMASS',
            47: 'AXIF2',
            48: 'AXIF3',
            49: 'AXIF4',
            50: 'SLOT3',
            51: 'SLOT4',
            52: 'CHBDYG/CHBDYP',
            53: 'CTRIAX6',
            54: 'TRIM6-old',
            55: 'CDUM3',
            56: 'CDUM4',
            57: 'CDUM5',
            58: 'CDUM6',
            59: 'CDUM7',
            60: 'CDUM8/CCRAC2D',
            61: 'CDUM9/CCRAC3D',
            62: 'CQDMEM1',
            63: 'IFQUAD',
            64: 'CQUAD8',
            65: 'IFHEX',
            66: 'IFPENT',
            67: 'CHEXA',
            68: 'CPENTA',
            69: 'CBEND',
            70: 'CTRIAR',
            71: None,
            72: 'AEROQ4',
            73: 'IFQDX',
            74: 'CTRIA3',
            75: 'CTRIA6',
            76: 'HEXPR',
            77: 'PENPR',
            78: 'TETPR',
            79: None,
            80: None,
            81: None,
            82: 'QUADR',
            83: 'HACAB',
            84: 'HACBR',
            85: 'TETRANL-nonlinear',
            86: 'GAPNL-nonlinear',
            87: 'TUBENL-nonlinear',
            88: 'TRIA3NL-nonlinear',
            89: 'RODNL-nonlinear',
            90: 'QUAD4NL-nonlinear',
            91: 'PENTANL-nonlinear',
            92: 'CONRODNL-nonlinear',
            93: 'HEXANL-nonlinear',
            94: 'BEAMNL-nonlinear',
            95: 'QUAD4LC-composite',
            96: 'QUAD8LC-composite',
            97: 'TRIA3LC-composite',
            98: 'TRIA6LC-composite',
            99: None,
            100: 'BAR-100',
            101: 'AABSF',
            102: 'CBUSH',
            103: 'QUADP',
            104: 'TRIAP',
            105: 'BEAMP',
            106: 'DAMP5',
            107: 'CHBDYE',
            108: 'CHBDYG',
            109: 'CHBDYP',
            110: 'CONV',
            111: 'CONVM',
            112: 'QBDY3',
            113: 'QVECT',
            114: 'QVOL',
            115: 'RADBC',
            116: 'SLIF1D',
            117: 'WELDC',
            118: 'WELDP',
            119: 'SEAM',
            120: 'GENEL',
            121: 'DMIG',
            122: 'DIEL-old',
            123: 'HEXAE-old',
            124: 'IND-old',
            125: 'LINE-old',
            126: 'FASTP',
            127: 'CQUAD-old',
            128: 'CQUADX-old',
            129: 'RELUC-old',
            130: 'RES-old',
            131: 'CTETRAE-old',
            132: 'CTRIA-old',
            133: 'CTRIAX-old',
            134: 'LINEOB-old',
            135: 'LINXOB-old',
            136: 'QUADOB-old',
            137: 'TRIAOB-old',
            138: 'LINEX-old',
            139: 'QUAD4FD',
            140: 'HEXA8FD',
            141: 'HEXAP',
            142: 'PENTAP',
            143: 'TETRAP',
            144: 'QUAD144',
            145: 'VUHEXA',
            146: 'VUPENTA',
            147: 'VUTETRA',
            148: 'HEXAM-old',
            149: 'PENTAM-old',
            150: 'TETRAM-old',
            151: 'QUAD-old',
            152: 'TRIAM-old',
            153: 'QUADXM-old',
            154: 'TRIAXM-old',
            155: 'RADINT',
            156: 'BUSH2D',
            157: 'LINEPW-old',
            158: 'QUADOBM-old',
            159: 'SEAMP',
            160: 'PENTA6FD',
            161: 'TETRA4FD',
            162: 'TRIA3FD',
            163: 'HEXAFD',
            164: 'QUADFD',
            165: 'PENTAFD',
            166: 'TETRAFD',
            167: 'TRIAFD',
            168: 'TRIAX3FD',
            169: 'TRIAXFD',
            170: 'QUADX4FD',
            171: 'QUADXFD',
            172: 'QUADRNL',
            173: 'TRIARNL',
            174: 'LINEOBM-old',
            175: 'LINXOBM-old',
            176: 'QUADWGM-old',
            177: 'TRIAWGM-old',
            178: 'QUADIB-old',
            179: 'TRIAIB-old',
            180: 'LINEIB-old',
            181: 'LINEXIB-old',
            182: 'QUADIBM-old',
            183: 'TRIAIBM-old',
            184: 'CBEAM3',
            185: 'LINXIBM-old',
            186: 'QUADPWM-old',
            187: 'TRIAPWM-old',
            188: 'LINEPWM-old',
            189: 'VUQUAD',
            190: 'VUTRIA',
            191: 'VUBEAM',
            192: 'CVINT',
            193: 'QUADFR-old',
            194: 'TRIAFR-old',
            195: 'LINEFR-old',
            196: 'LINXFR-old',
            197: 'SFINT',
            198: 'CNVPEL',
            199: 'VUHBDY',
            200: 'WELD',
            201: 'QUAD4FD',
            202: 'HEXA8FD',
            203: 'SLIF1D',
            204: 'PENTA6FD',
            205: 'TETRA4FD',
            206: 'TRIA3FD',
            207: 'HEXAFD',
            208: 'QUADFD',
            209: 'PENTAFD',
            210: 'TETRAFD',
            211: 'TRIAFD',
            212: 'TRIAX3FD',
            213: 'TRIAXFD',
            214: 'QUADX4FD',
            215: 'QUADXFD',
            216: 'TETRA4FD',
            217: 'TRIA3FD',
            218: 'HEXAFD',
            219: 'QUADFD',
            220: 'PENTAFD',
            221: 'TETRAFD',
            222: 'TRIAX3FD',
            223: 'QUADXFD',
            224: 'ELAS1',
            225: 'ELAS3',
            226: 'BUSH',
            227: 'RBAR',
            228: 'RBE1',
            229: 'RBE3',
            230: 'RJOINT',
            231: 'RROD',
            232: 'QUADRLC',
            233: 'TRIARLC',
            234: None,
            235: 'CQUADR',
            236: 'CTRIAR',
            237 : 'CTRIAR',
            238 : 'CBAR',
            239 : 'CEAM',
            240 : 'CBAR',
            241 : 'AXISYM',
        }
        nx_elements = {
            0 : 'GRID',
            1 : 'CROD',
            2 : 'CBEAM',
            3 : 'CTUBE',
            4 : 'CSHEAR',
            5 : 'FORCEi/MOMENTi follower',
            6 : 'CTRIA1-old',
            7 : 'PLOAD4 follower',
            8 : 'FLOAD/PLOAD2 follower',
            9 : 'PLOAD/PLOAD2',
            10 : 'CONROD',
            11 : 'CELAS1',
            12 : 'CELAS2',
            13 : 'CELAS3',
            14 : 'CELAS4',
            15 : 'AEROT3',
            16 : 'AEROBEAM',
            17 : 'CTRIA2-old',
            18 : 'CQUAD2-old',
            19 : 'CQUAD1-old',
            20 : 'CDAMP1',
            21 : 'CDAMP2',
            22 : 'CDAMP3',
            23 : 'CDAMP4',
            24 : 'CVISC',
            25 : 'CMASS1',
            26 : 'CMASS2',
            27 : 'CMASS3',
            28 : 'CMASS4',
            29 : 'CONM1',

            30 : 'CONM2',
            31 : 'PLOTEL',
            32 : None,
            33 : 'CQUAD4',
            34 : 'CBAR-34',
            35 : 'CCONE',
            36 : 'CTRIARG-old',
            37 : 'CTRAPRG-old',
            38 : 'CGAP',
            39 : 'CTETRA',

            40 : 'CBUSH1D',
            41 : 'CHEXA1-old',
            42 : 'CHEXA2-old',
            43 : 'CFLUID2',
            44 : 'CFLUID3',
            45 : 'CFLUID4',
            46 : 'FLMASS',
            47 : 'AXIF2',
            48 : 'AXIF3',
            49 : 'AXIF4',

            50 : 'SLOT3',
            51 : 'SLOT4',
            52 : 'CHBDYG/CHBDYP',
            53 : 'CTRIAX6',
            54 : 'TRIM6-old',
            55 : 'CDUM3',
            56 : 'CDUM4',
            57 : 'CDUM5',
            58 : 'CDUM6',
            59 : 'CDUM7',

            60 : 'CDUM8',
            61 : 'CDUM9',
            62 : 'CQDMEM1-old',
            63 : 'CQDMEM2-old',
            64 : 'CQUAD8',
            65 : 'CHEX8-old',
            66 : 'CHEX20-old',
            67 : 'CHEXA',
            68 : 'CPENTA',
            69 : 'CBEND',

            70 : 'CTRIAR',
            71 : '',
            72 : 'AEROQ4',
            73 : 'CFTUBE-old',
            74 : 'CTRIA3',
            75 : 'CTRIA6',
            76 : 'CHEXPR',
            77 : 'CPENPR',
            78 : 'CTETPR',
            79 : 'CPYRAM',
            80 : '',
            81 : '',
            82 : 'CQUADR',
            83 : 'HACAB',
            84 : 'HACBR',
            85 : 'TETRA-nonlinear',
            86 : 'GAP-nonlinear',
            87 : 'TUBE-nonlinear',
            88 : 'TRIA3-nonlinear',
            89 : 'ROD-nonlinear',
            90 : 'QUAD4-nonlinear',
            91 : 'PENTA-nonlinear',
            92 : 'CONROD-nonlinear',
            93 : 'HEXA-nonlinear',
            94 : 'BEAM-nonlinear',
            95 : 'QUAD4-nonlinear',
            96 : 'QUAD8-nonlinear',
            97 : 'TRIA3-nonlinear',
            98 : 'TRIA6-nonlinear',
            99 : '',
            100 : 'CBAR-100',
            101 : 'AABSF',
            102 : 'CBUSH',
            103 : 'CQUADP',
            104 : 'CTRIAP',
            105 : 'CBEAMP',
            106 : 'CDAMP5',
            107 : 'CHBDYE',
            108 : 'CHBDYG',
            109 : 'CHBDYP',

            110 : 'CONV',
            111 : 'CONVM',
            112 : 'QBDY3',
            113 : 'QVECT',
            114 : 'QVOL',
            115 : 'RADBC',
            116 : 'SLIF1D',
            117 : 'CWELDC',  # unlisted in the main table, used in OEF table
            118 : 'CWELDP',  # unlisted in the main table, used in OEF table
            119 : 'CFAST',  # unlisted in the main table, used in OEF table
            120 : '',

            121 : '',
            122 : '',
            123 : '',
            124 : '',
            125 : '',
            126 : '',
            127 : 'CQUAD',  # unlisted in the main table, used in OEF table
            128 : 'CQUADX',  # unlisted in the main table, used in OEF table
            129 : 'RELUC',  # unlisted in the main table, used in OEF table
            130 : 'RES',  # unlisted in the main table, used in OEF table

            131 : 'TETRAE',  # unlisted in the main table, used in OEF table
            132 : 'CTRIA',  # unlisted in the main table, used in OEF table
            133 : 'CTRIAX',  # unlisted in the main table, used in OEF table
            134 : 'LINEOB',  # unlisted in the main table, used in OEF table
            135 : 'LINXOB',  # unlisted in the main table, used in OEF table
            136 : 'QUADOB',  # unlisted in the main table, used in OEF table
            137 : 'TRIAOB',  # unlisted in the main table, used in OEF table
            138 : 'LINX',  # unlisted in the main table, used in OEF table
            139 : 'CQUAD4FD',
            140 : 'CHEXA8FD',

            141 : 'CHEXAP',
            142 : 'CPENTAP',
            143 : 'CTETRAP',
            144 : 'CQUAD144',
            145 : 'VUHEXA',
            146 : 'VUPENTA',
            147 : 'VUTETRA',
            148 : 'HEXAM',  # unlisted in the main table, used in OEF table
            149 : 'PENTAM',  # unlisted in the main table, used in OEF table
            150 : 'TETRAM',  # unlisted in the main table, used in OEF table

            151 : 'QUADM',  # unlisted in the main table, used in OEF table
            152 : 'TRIAM',  # unlisted in the main table, used in OEF table
            153 : 'QUADXM',  # unlisted in the main table, used in OEF table
            154 : 'TRIAXM',  # unlisted in the main table, used in OEF table
            155 : 'QUADPW',  # unlisted in the main table, used in OEF table
            156 : 'TRIAPW',  # unlisted in the main table, used in OEF table
            157 : 'LINEPW',  # unlisted in the main table, used in OEF table
            158 : 'QUADOBM',  # unlisted in the main table, used in OEF table
            159 : 'TRIAOBM',  # unlisted in the main table, used in OEF table
            160 : 'CPENTA6FD',

            161 : 'CTETRA4FD',
            162 : 'CTRIA3FD',
            163 : 'CHEXAFD',
            164 : 'CQUADFD',
            165 : 'CPENTAFD',
            166 : 'CTETRAFD',
            167 : 'CTRIAFD',
            168 : 'CTRIAX3FD',
            169 : 'CTRIAXFD',
            170 : 'CQUADX4FD',

            171 : 'CQUADXFD',
            172 : '',
            173 : '',
            174 : 'LINEOBM',  # unlisted in the main table, used in OEF table
            175 : 'LINXOBM',  # unlisted in the main table, used in OEF table
            176 : 'QUADWGM',  # unlisted in the main table, used in OEF table
            177 : 'TRIAWGM',  # unlisted in the main table, used in OEF table
            178 : 'QUADIB',  # unlisted in the main table, used in OEF table
            179 : 'TRIAIB',  # unlisted in the main table, used in OEF table
            180 : 'LINEIB',  # unlisted in the main table, used in OEF table

            181 : 'LINXIB',  # unlisted in the main table, used in OEF table
            182 : 'QUADIBM',  # unlisted in the main table, used in OEF table
            183 : 'TRIAIBM',  # unlisted in the main table, used in OEF table
            184 : 'LINEIBM',  # unlisted in the main table, used in OEF table
            185 : 'LINXIBM',  # unlisted in the main table, used in OEF table
            186 : 'QUADPWM',  # unlisted in the main table, used in OEF table
            187 : 'TRIAPWM',  # unlisted in the main table, used in OEF table
            188 : 'LINEPWM',  # unlisted in the main table, used in OEF table
            189 : 'VUQUAD',
            190 : 'VUTRIA',

            191 : 'VUBEAM',
            192 : 'CVINT',
            193 : 'QUADFR',  # unlisted in the main table, used in OEF table
            194 : 'TRIAFR',  # unlisted in the main table, used in OEF table
            195 : 'LINEFR',  # unlisted in the main table, used in OEF table
            196 : 'LINXFR',  # unlisted in the main table, used in OEF table
            197 : 'SFINT',  # TODO: GMINTS-OEF??
            198 : 'CNVPEL',
            199 : 'VUHBDY',
            200 : 'CWELD',

            201 : 'CQUAD4FD',
            202 : 'CHEXA8FD',
            203 : 'SLIF1D',
            204 : 'CPENTA6FD',
            205 : 'CTETRA4FD',
            206 : 'CTRIA3FD',
            207 : 'CHEXAFD',
            208 : 'CQUADFD',
            209 : 'CPENTAFD',
            210 : 'CTETRAFD',

            211 : 'CTRIAFD',
            212 : 'CTRIAX3FD',
            213 : 'CTRIAXFD',
            214 : 'CQUADX4FD',
            215 : 'CQUADXFD',
            216 : 'CTETRA4FD',
            217 : 'CTRIA3FD',
            218 : 'CHEXAFD',
            219 : 'CQUADFD',
            220 : 'CPENTAFD',

            221 : 'CTETRAFD',
            222 : 'CTRIAX3FD',
            223 : 'CQUADXFD',
            224 : 'CELAS1',
            225 : 'CELAS3',
            226 : 'CBUSH',
            227 : 'CTRIAR',
            228 : 'CQUADR',
            229 : '',
            230 : '',

            231 : '',
            232 : 'CQUADR-composite',
            233 : 'CTRIAR-composite',
            234 : '',
            235 : '',
            236 : '',
            237 : '',
            238 : '',
            239 : '',

            # per nx 8.5 QRG
            240 : 'CTRIA6',
            241 : 'CQUAD8',
            242 : 'CTRAX3',
            243 : 'CQUADX4',
            244 : 'CTRAX6',
            245 : 'CQUADX8',
            246 : 'CTRAX3',
            247 : 'CQUADX4',
            248 : 'CTRAX6',
            249 : 'CQUADX8',
            255 : 'CPYRAM',
            256 : 'CPYRAM',
            257 : 'CPYRAMFD',
            258 : 'CPYRAMFD',
            259 : 'CTRAX3FD',
            261 : 'CTRAX3FD',
            262 : 'CQUADX4FD',
            263 : 'CTRAX6FD',
            264 : 'QUADX3FD',
            266 : 'QUADX8FD',
            271 : 'CPLSTN3',
            272 : 'CPLSTN4',
            273 : 'CPLSTN6',
            274 : 'CPLSTN3',
            275 : 'CPLSTS3',
            276 : 'CPLSTS4',
            277 : 'CPLSTS6',
            278 : 'CPLSTS8',
            281 : 'CPLSTN3',
            282 : 'CPLSTN4',
            283 : 'CPLSTN6',
            284 : 'CPLSTN8',
            285 : 'CPLSTS3',
            287 : 'CPLSTS6',
            288 : 'CPLSTS8',
            289 : 'CPLSTN4',
            290 : 'CPLSTS4',
            291 : 'CPLSTN3',
            292 : 'CPLSTN4',
            293 : 'CPLSTN6',
            294 : 'CPLSTS8',
            295 : 'CPLSTS3',
            296 : 'CPLSTS4',
            297 : 'CPLSTS6',
            298 : 'CPLSTS8',
        }
        return msc_elements, nx_elements

    def print_table_code(self, table_code):
        table_code_content = table_code % 1000
        #dataFormat = table_code/1000
        msg = ''
        #msg += 'tableCodeContent=%s dataFormat=%s\n' %(tableCodeContent,dataFormat)

        msc_table_content = {
            # dmap 2014
            0: '',
            1: 'OUG - Displacement vector',
            2: 'OPG - Load vector',
            3: 'OQG - SPC Force vector',
            4: 'OEF - Element force/flux',
            5: 'OES - Element stress/strain',
            6: 'LAMA - Eigenvalue summary',
            7: 'OUG - Eigenvector',
            8: 'Grid Point Singularity Table (obsolete)',
            9: 'OEIGS - Eigenvalue analysis summary',
            10: 'OUG - Velocity vector',
            11: 'OUG - Acceleration vector',
            12: 'OPG - Nonlinear force vector',
            13: 'OGPWG - Grid point weight generator',
            14: 'OUG - Eigenvector (solution set)',
            15: 'OUG - Displacement vector (solution set)',
            16: 'OUG - Velocity vector (solution set)',
            17: 'OUG - Acceleration vector (solutin set)',
            18: 'OEE - Element strain energy',
            19: 'OGF - Grid point force balance',
            20: 'OES - Stresses at grid points',
            21: 'OES - Strain/curvature at grid points',
            22: 'OELOF1 - Element internal forces/moments',
            23: 'OELOP1 - Summation of element oriented forces on adjacent elements',
            24: 'OEP - Element pressures',
            25: 'OEF - Composite failure indices',
            26: 'OGS - Grid point stresses (surface)',
            27: 'OGS - Grid point stresses (volume - direct)',
            28: 'OGS - Grid point stresses (volume - princial)',
            29: 'OGS - Element stress discontinuities (surface)',
            30: 'OGS - Element stress discontinuities (volume - direct)',
            31: 'OGS - Element stress discontinuities (volume - princial)',
            32: 'OGS - Grid point stress discontinuities (surface)',
            33: 'OGS - Grid point stress discontinuities (volume - direct)',
            34: 'OGS - Grid point stress discontinuities (volume - princial)',
            35: 'OGS - Grid point stresses (plane strain)',
            36: 'OEE - Element kinetic energy',
            37: 'OEE - Element energy loss',
            38: 'OMM - MaxMin summary',
            39: 'OQG - MPC forces',
            40: 'OGPKE - Grip point kinetic energy',
            #51: 'OFMPF2M - ???',
            #52: 'OSMPF2M - ???',
            #53: 'OPMPF2M - ???',
            #54: 'OLMPF2M - ???',
            #55: 'OGMPF2M - ???',
        }
        nx_table_content = {
            # nx 8.5
            0: '',
            1: 'OUG - Displacement vector',
            2: 'OPG - Load vector',
            3: 'OQG - SPC/MPC Force vector',
            4: 'OEF - Element force/flux',
            5: 'OES - Element stress/strain',
            6: 'LAMA - Eigenvalue summary',
            7: 'OUG - Eigenvector',
            8: 'Grid Point Singularity Table (obsolete)',
            9: 'OEIGS - Eigenvalue analysis summary',
            10: 'OUG - Velocity vector',
            11: 'OUG - Acceleration vector',
            12: 'OPG - Nonlinear force vector',
            13: 'OGPWG - Grid point weight generator',
            14: 'OUG - Eigenvector (solution set)',
            15: 'OUG - Displacement vector (solution set)',
            16: 'OUG - Velocity vector (solution set)',
            17: 'OUG - Acceleration vector (solutin set)',
            18: 'OEE - Element strain energy',
            19: 'OGF - Grid point force balance',
            20: 'OES - Stresses at grid points',
            21: 'OES - Strain/curvature at grid points',
            22: 'OELOF1 - Element internal forces/moments',
            23: 'OELOP1 - Summation of element oriented forces on adjacent elements',
            24: 'OEP - Element pressures',
            25: 'OEF - Composite failure indices',
            26: 'OGS - Grid point stresses (surface)',
            27: 'OGS - Grid point stresses (volume - direct)',
            28: 'OGS - Grid point stresses (volume - princial)',
            29: 'OGS - Element stress discontinuities (surface)',
            30: 'OGS - Element stress discontinuities (volume - direct)',
            31: 'OGS - Element stress discontinuities (volume - princial)',
            32: 'OGS - Grid point stress discontinuities (surface)',
            33: 'OGS - Grid point stress discontinuities (volume - direct)',
            34: 'OGS - Grid point stress discontinuities (volume - princial)',
            35: 'OGS - Grid point stresses (plane strain)',
            36: 'OEE - Element kinetic energy',
            37: 'OEE - Element energy loss',

            38 : 'OMSEC - Constant modal strain energy',
            39 : 'OMSED - Oscillating modal strain energy',
            40 : 'OMKEC - Constant modal kinetic energy',
            41 : 'OMKED - Oscillating modal kinetic energy',
            42 : 'OMECON - Constant total modal energy',
            43 : 'OMEOSC - Oscillating total modal energy',
            44 : 'OUGMC - Displacement/velocity/acceleration modal contributions',
            45 : 'OEFMC - Element force modal contributions',
            46 : 'OESMC - Element stress modal contributions',
            47 : 'OSTRMC - Element strain modal contributions',
            48 : 'OQGMC - SPC force modal contributions',
            49 : 'OUGPC - Panel contributions',
            50 : 'OUGGC - Grid contributions',
            51 : 'OUGRC - Reciprocal panel contributions',
            61 : 'OGK - Gasket Element Results',
            62 : 'OBC - Contact Pressure and Traction',
            63 : 'OQG - Contact Force Results',
            64 : 'OSPDSI - Contact Separation Distance - Initial',
            65 : 'OSPDS - Contact Separation Distance',
            66 : 'OBG - Glue force results (normal and in-plane tractions)',
            67 : 'OQG - Glue force resutls',
            68 : '??? - Tosca normalized material properties',
        }
        if self.is_msc:
            msg += 'n=%s table=%s-%s' % (self.n, self.table_name,
                                         msc_table_content[table_code_content])
        else:
            msg += 'n=%s table=%s-%s' % (self.n, self.table_name,
                                         nx_table_content[table_code_content])

        return msg

    def approach_code_str(self, approach_code):
        return ''

    def code_information(self, include_time=True):
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
            if self.is_stress():
                stress_word = 'Stress'
            else:
                stress_word = 'Strain'
            s_word = get_scode_word(s_code, self.stress_bits)

        element_type = None
        if hasattr(self, 'element_type'):
            element_type = self.element_type

        format_word = '???'
        if format_code == 1:
            format_word = "Real"
        elif format_code == 2:
            format_word = "Real/Imaginary"
        elif format_code == 3:
            format_word = "Magnitude/Phase"
        else:
            format_word = '\n%18s1 - Real\n' % ''
            format_word += '%18s2 - Real/Imaginary\n' % ''
            format_word += '%18s3 - Magnitude/Phase\n' % ''
            #msg = 'unsupported format_code:  format_code=%s\n' % format_code
            #raise InvalidFormatCodeError(msg)

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

        if thermal == 0:
            thermal_word = 'isHeatTransfer = False'
        elif thermal == 1:
            thermal_word = 'isHeatTransfer = True'
        elif thermal == 2:
            thermal_word = 'Scaled response spectra ABS'
        #elif thermal == 3:
            #thermal_word = 'Scaled response spectra SRSS'
        elif thermal == 4:
            thermal_word = 'Scaled response spectra SRSS' # NRL???
        elif thermal == 5:
            thermal_word = 'Scaled response spectra NRLO'
        elif thermal == 8:
            thermal_word = 'Scaled response spectra NRL'
        else:
            thermal_word = '???'
            #msg = 'unsupported thermal:  thermal=%s\n' %(thermal)
            #raise ValueError(msg)

        analysis = '???'
        if self.analysis_code == 1:
            analysis = "Statics"
        elif self.analysis_code == 2:
            analysis = "Normal modes or buckling (real eigenvalues)"
        elif self.analysis_code == 3:
            analysis = "Differential Stiffness 0 - obsolete"
        elif self.analysis_code == 4:
            analysis = "Differential Stiffness 1 - obsolete"
        elif self.analysis_code == 5:
            analysis = "Frequency"
        elif self.analysis_code == 6:
            analysis = "Transient"
        elif self.analysis_code == 7:
            analysis = "Pre-buckling"
        elif self.analysis_code == 8:
            analysis = "Post-buckling"
        elif self.analysis_code == 9:
            analysis = "Complex eigenvalues"
        elif self.analysis_code == 10:
            analysis = "Nonlinear statics"
        elif self.analysis_code == 11:
            analysis = "Geometric nonlinear statics"

        device = '???'
        if self.device_code == 1:
            device = "Print"
        elif self.device_code == 2:
            device = "Plot"
        elif self.device_code == 3:
            device = "Print and Plot"
        elif self.device_code == 4:
            device = "Punch"
        elif self.device_code == 5:
            device = "Print and Punch"
        elif self.device_code == 6:
            device = "Plot and Punch"
        elif self.device_code == 7:
            device = "Print, Plot, and Punch"

        if thermal == 0:
            force_flux = 'Force'
        elif thermal == 1:
            force_flux = 'Flux'
        else:
            force_flux = 'Force (or Flux); thermal=%r' % thermal

        if thermal == 0:
            disp_temp = 'Displacement'
        elif thermal == 1:
            disp_temp = 'Temperature'
        #elif thermal is None:
            #raise RuntimeError('thermal_code is not specified; thermal_code=None')
        else:
            disp_temp = 'Displacement/Temperature; thermal=%r' % thermal

        table = '???'
        table_code = self.table_code
        if table_code in [601, 610, 611]:
            table_code -= 600
        if table_code == 1:
            table = "OUG - %s vector/scalar" % disp_temp
        elif table_code == 2:
            table = "OPG - Load vector"
        elif table_code == 3:
            table = "OQG - SPC Force vector"
        elif table_code == 4:
            table = "OEF - Element %s" % force_flux
        elif table_code == 5:
            table = "OES - Element %s" % stress_word
        elif table_code == 6:
            table = "LAMA - Eigenvalue summary"
        elif table_code == 7:
            table = "OUG - Eigenvector"
        elif table_code == 8:
            table = "none - Grid point singularity table (obsolete)"
        elif table_code == 9:
            table = "OEIGS - Eigenvalue analysis summary"
        elif table_code == 10:
            table = "OUG - Velocity vector"
        elif table_code == 11:
            table = "OUG - Acceleration vector"
        elif table_code == 12:
            table = "OPG - Nonlinear force vector"
        elif table_code == 13:
            table = "OGPWG - Grid point weight generator"
        elif table_code == 14:
            table = "OUG - Eigenvector (solution set)"
        elif table_code == 15:
            table = "OUG - Displacement vector (solution set)"
        elif table_code == 16:
            table = "OUG - Velocity vector (solution set)"
        elif table_code == 17:
            table = "OUG - Acceleration vector (solution set)"
        elif table_code == 18:
            table = "OEE - Element strain energy"
        elif table_code == 19:
            table = "OGF - Grid point force balance"
        elif table_code == 20:
            table = "OES - Stresses at grid points (from the CURV module)"
        elif table_code == 21:
            table = "OES - Strain/curvature at grid points"
        elif table_code == 22:
            table = "OELOF1 - Element internal forces and moments"
        elif table_code == 23:
            table = "OELOP1 - Summation of element oriented forces on adjacent elements"
        elif table_code == 24:
            table = "OEP - Element pressures"
        elif table_code == 25:
            table = "OEF - Composite failure indicies"
        elif table_code == 26:
            table = "OGS - Grid point stresses (surface)"
        elif table_code == 27:
            table = "OGS - Grid point stresses (volume -- direct)"
        elif table_code == 28:
            table = "OGS - Grid point stresses (volume -- principal)"
        elif table_code == 29:
            table = "OGS - Element stress discontinuities (surface)"
        elif table_code == 30:
            table = "OGS - Element stress discontinuities (volume -- direct)"
        elif table_code == 31:
            table = "OGS - Element stress discontinuities (volume -- principal)"
        elif table_code == 32:
            table = "OGS - Grid point stress discontinuities (surface)"
        elif table_code == 33:
            table = "OGS - Grid point stress discontinuities (volume -- direct)"
        elif table_code == 34:
            table = "OGS - Grid point stress discontinuities (volume -- principal)"
        elif table_code == 35:
            table = "OGS - Grid point stress discontinuities (plane strain)"
        elif table_code == 36:
            table = "OEE - Element kinetic energy"
        elif table_code == 37:
            table = "OEE - Element energy loss"
        elif table_code == 38:
            table = "OMM - Max/Min summary"
        elif table_code == 39:
            table = "OQG - MPC Forces"
        elif table_code == 40:
            table = "OGPKE - Grip point kinetic energy"
        else:
            table = '%s - Unknown' % self.table_name
        msg = '--Table3Data--\n\n'
        msg += "  device_code   = %-3s %s\n" % (self.device_code, device)
        msg += "  analysis_code = %-3s %s\n" % (self.analysis_code, analysis)
        msg += "  table_code    = %-3s %s-%s\n" % (self.table_code, self.table_name, table)
        msg += "  format_code   = %-3s %s\n" % (format_code, format_word)

        msg += "  dataFormat    = %-3s %s\n" % (self.sort_bits[0], sort_word1)
        msg += "  sortType      = %-3s %s\n" % (self.sort_bits[1], sort_word2)
        msg += "  isRandom      = %-3s %s\n" % (self.sort_bits[2], sort_word3)
        random_code = self.random_code if hasattr(self, 'random_code') else 0
        msg += "  random_code   = %-3s\n" % (random_code)

        if element_type is not None:
            if isinstance(element_type, string_types):
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

        if self.is_msc:
            msg += '  MSC Nastran\n'
        else:
            msg += '  NX Nastran\n'
        #print msg
        assert isinstance(self.format_code, int), type(self.format_code)
        return msg

    #----
    def is_thermal(self):
        if self.thermal == 0:
            return False
        elif self.thermal == 1:
            return True
        return '???'

    #----
    # format_code 3
    def is_magnitude_phase(self):
        if self.format_code == 3:
            return True
        return False


    def is_sort1(self):
        sort_method, is_real, is_random = self._table_specs()
        return True if sort_method == 1 else False

    def is_sort2(self):
        sort_method, is_real, is_random = self._table_specs()
        return True if sort_method == 2 else False

    def _table_specs(self):
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
        tcode = self.sort_code
        sort_method = 1
        is_real = True
        is_random = False
        assert tcode in [0, 1, 2, 3, 4, 5, 6], 'tcode=%s\n%s' % (tcode, self.code_information())
        if tcode in [2, 3, 5, 6]:
            sort_method = 2
        if tcode in [1, 3]:
            is_real = False
        if tcode in [4, 5, 6]:
            is_random = True

        if is_random:
            assert self.sort_bits[0] == 1, 'should be RANDOM; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)
        else:
            assert self.sort_bits[0] == 0, 'should be NOT RANDOM; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)

        if sort_method == 1:
            assert self.sort_bits[1] == 0, 'should be SORT1; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)
        else:
            assert self.sort_bits[1] == 1, 'should be SORT2; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)

        if is_real:
            assert self.sort_bits[2] == 0, 'should be REAL; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)
        else:
            assert self.sort_bits[2] == 1, 'should be IMAG; sort_bits=%s; tcode=%s' % (self.sort_bits, tcode)
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

    def _set_op2_date(self, month, day, year):
        self.date = (month, day, 2000 + year)
        return self.date

def get_scode_word(s_code, stress_bits):
    if s_code == 0:
        s_word = 'Coordinate Element - Stress Max Shear (Octahedral)'
        assert stress_bits == [0, 0, 0, 0, 0], stress_bits
    elif s_code == 14:
        s_word = 'Coordinate Element - Strain Fiber Max Shear (Octahedral)'
        assert stress_bits == [0, 1, 1, 1, 0], stress_bits

    elif s_code == 1:
        s_word = 'Coordinate Element - Stress von Mises'
        assert stress_bits == [0, 0, 0, 0, 1], stress_bits
    elif s_code == 10:
        s_word = 'Coordinate Element - Strain Curvature Max Shear (Octahedral)'
        assert stress_bits == [0, 1, 0, 1, 0], stress_bits

    elif s_code == 11:
        s_word = 'Coordinate Element - Strain Curvature von Mises'
        assert stress_bits == [0, 1, 0, 1, 1], stress_bits
    elif s_code == 15:
        s_word = 'Coordinate Element - Strain Fiber von Mises'
        assert stress_bits == [0, 1, 1, 1, 1], stress_bits

    elif s_code == 16:
        s_word = 'Coordinate Material - Stress Max Shear (Octahedral)'
        assert stress_bits == [1, 0, 0, 0, 0], stress_bits
    elif s_code == 17:
        s_word = 'Coordinate Material - Stress von Mises'
        assert stress_bits == [1, 0, 0, 0, 1], stress_bits

    elif s_code == 26:
        s_word = 'Coordinate Material - Strain Curvature Max Shear'
        assert stress_bits == [1, 1, 0, 1, 0], stress_bits
    elif s_code == 30:
        s_word = 'Coordinate Material - Strain Fiber Max Shear (Octahedral)'
        assert stress_bits == [1, 1, 1, 1, 0], stress_bits

    elif s_code == 27:
        s_word = 'Coordinate Material - Strain Curvature von Mises'
        assert stress_bits == (1, 1, 0, 1, 1), stress_bits
    elif s_code == 31:
        s_word = 'Coordinate Material - Strain Fiber von Mises'
        assert stress_bits == (1, 1, 1, 1, 1), stress_bits
    else:
        #s_word = 'Stress or Strain - UNDEFINED'
        s_word = '???'
    return s_word