from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)


class Op2Codes(object):
    def __init__(self):
        pass

    def get_element_type(self, eCode):
        elements = {
            None: '',
            0: 'GRID',
            1: 'ROD',
            2: 'BEAM',
            3: 'TUBE',
            4: 'SHEAR',
            5: 'FORMON12',
            6: 'FORCE',
            7: 'PLOAD4',
            8: 'PLOADX1',
            9: 'PLOAD/PLOAD2',
            10: 'CONROD',
            11: 'ELAS1',
            12: 'ELAS2',
            13: 'ELAS3',
            14: 'ELAS4',
            15: 'AEROT3',
            16: 'AEROBEAM',
            17: None,
            18: None,
            19: None,
            20: 'DAMP1',
            21: 'DAMP2',
            22: 'DAMP3',
            23: 'DAMP4',
            24: 'VISC',
            25: 'MASS1',
            26: 'MASS2',
            27: 'MASS3',
            28: 'MASS4',
            29: 'CONM1',
            30: 'CONM2',
            31: 'PLOTEL',
            32: None,
            33: 'QUAD4',
            34: 'BAR',
            35: 'CON',
            36: None,
            37: None,
            38: 'GAP',
            39: 'TETRA',
            40: 'BUS1D',
            41: None,
            42: None,
            43: 'FLUID2',
            44: 'FLUID3',
            45: 'FLUID4',
            46: 'FLMASS',
            47: 'AXIF2',
            48: 'AXIF3',
            49: 'AXIF4',
            50: 'SLOT3',
            51: 'SLOT4',
            52: 'HBDY',
            53: 'TRIAX6',
            54: None,
            55: 'DUM3',
            56: 'DUM4',
            57: 'DUM5',
            58: 'DUM6',
            59: 'DUM7',
            60: 'DUM8',
            61: 'DUM9',
            62: None,
            63: None,
            64: 'QUAD8',
            65: None,
            66: None,
            67: 'HEXA',
            68: 'PENTA',
            69: 'BEND',
            70: 'TRIAR',
            71: None,
            72: 'AEROQ4',
            73: None,
            74: 'TRIA3',
            75: 'TRIA6',
            76: 'HEXPR',
            77: 'PENPR',
            78: 'TETPR',
            79: None,
            80: None,
            81: None,
            82: 'QUADR',
            83: 'HACAB',
            84: 'HACBR',
            85: 'TETRANL',
            86: 'GAPNL',
            87: 'TUBENL',
            88: 'TRIA3NL',
            89: 'RODNL',
            90: 'QUAD4NL',
            91: 'PENTANL',
            92: 'CONRODNL',
            93: 'HEXANL',
            94: 'BEAMNL',
            95: 'QUAD4LC',
            96: 'QUAD8LC',
            97: 'TRIA3LC',
            98: 'TRIA6LC',
            99: None,
            100: 'BARS',
            101: 'AABSF',
            102: 'BUSH',
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
            119: 'GENEL',
            120: 'DMIG',
            121: None,
            122: None,
            123: None,
            124: None,
            125: None,
            126: None,
            127: None,
            128: None,
            129: None,
            130: None,
            131: None,
            132: None,
            133: None,
            134: None,
            135: None,
            136: None,
            137: None,
            138: None,
            139: 'QUAD4FD',
            140: 'HEXA8FD',
            141: 'HEXAP',
            142: 'PENTAP',
            143: 'TETRAP',
            144: 'QUAD144',
            145: 'VUHEXA',
            146: 'VUPENTA',
            147: 'VUTETRA',
            148: None,
            149: None,
            150: None,
            151: None,
            152: None,
            153: None,
            154: None,
            155: None,
            156: None,
            157: None,
            158: None,
            159: None,
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
            174: None,
            175: None,
            176: None,
            177: None,
            178: None,
            179: None,
            180: None,
            181: None,
            182: None,
            183: None,
            184: None,
            185: None,
            186: None,
            187: None,
            188: None,
            189: 'VUQUAD',
            190: 'VUTRIA',
            191: 'VUBEAM',
            192: 'CVINT',
            193: None,
            194: None,
            195: None,
            196: None,
            197: 'SFINT',
            198: 'CNVPEL',
            199: 'VUHBDY',
            200: 'WELD',
            201: 'QUAD4FD',
            202: 'HEXA8FD',
            203: 'SLIF1D?',
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
            234: '???',
            235: 'CQUADR',  # was blank in DMAP, found reference in OEF table
            236: 'CTRIAR',  # was blank in DMAP, found reference in OEF table
        }
        return elements[eCode]  # +'_'+str(eCode)

    def print_table_code(self, table_code):
        tableCodeContent = table_code % 1000
        #dataFormat = table_code/1000
        msg = ''
        #msg += 'tableCodeContent=%s dataFormat=%s\n' %(tableCodeContent,dataFormat)

        tableContent = {
            0: '',
            1: 'OUG - Displacement vector',
            2: 'OPG - Load vector',
            3: 'OQG - SPC Force vector',
            4: 'OEF - Element force/flux',
            5: 'OES - Element stress/strain',
            6: 'LAMA - Eigenvalue summary',
            7: 'OUG - Eigenvector',
            8: None,
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
            35: 'OGS - Grid point stresses (plane stress)',
            36: 'OEE - Element kinetic energy',
            37: 'OEE - Element energy loss',
            38: 'OMM - MaxMin summary',
            39: 'OQG - MPC forces',
            40: 'OGPKE - Grip point kinetic energy',
            51: 'OFMPF2M - ???',
            52: 'OSMPF2M - ???',
            53: 'OPMPF2M - ???',
            54: 'OLMPF2M - ???',
            55: 'OGMPF2M - ???',
        }
        msg += 'n=%s table=%s-%s' % (self.n, self.table_name,
                                     tableContent[tableCodeContent])
        return msg

    def code_information(self):
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

        stressWord = ''
        if hasattr(self, 'stress_bits'):
            if self.isStress():
                stressWord = 'Stress'
            else:
                stressWord = 'Strain'

        element_type = None
        if hasattr(self, 'element_type'):
            element_type = self.element_type

        sWord = ''
        if s_code == 0:
            sWord += 'Coordinate Element - Stress Max Shear (Octahedral)'
        elif s_code == 14:
            sWord += 'Coordinate Element - Strain Fiber Max Shear (Octahedral)'

        elif s_code == 1:
            sWord += 'Coordinate Element - Stress von Mises'
        elif s_code == 10:
            sWord += 'Coordinate Element - Strain Curvature Max Shear (Octahedral)'

        elif s_code == 11:
            sWord += 'Coordinate Element - Strain Curvature von Mises'
        elif s_code == 15:
            sWord += 'Coordinate Element - Strain Fiber von Mises'

        elif s_code == 16:
            sWord += 'Coordinate Material - Stress Max Shear (Octahedral)'
        elif s_code == 17:
            sWord += 'Coordinate Material - Stress von Mises'

        elif s_code == 26:
            sWord += 'Coordinate Material - Strain Curvature Max Shear'
        elif s_code == 30:
            sWord += 'Coordinate Material - Strain Fiber Max Shear (Octahedral)'

        elif s_code == 27:
            sWord += 'Coordinate Material - Strain Curvature von Mises'
        elif s_code == 31:
            sWord += 'Coordinate Material - Strain Fiber von Mises'
        else:
            #sWord = 'Stress or Strain - UNDEFINED'
            sWord = ''

        formatWord = '???'
        if format_code == 1:
            formatWord = "Real"
        elif format_code == 2:
            formatWord = "Real/Imaginary"
        elif format_code == 3:
            formatWord = "Magnitude/Phase"
        else:
            formatWord = '\n%18s1 - Real\n%18s2-Real/Imaginary\n%18s3-Magnitude/Phase\n' % ('', '', '')
            #msg = 'unsupported format_code:  format_code=%s\n' % format_code
            #raise InvalidFormatCodeError(msg)

        if   self.sort_bits[0] == 0:
            sortWord1 = 'Sort1'
        else:
            sortWord1 = 'Sort2'
        if   self.sort_bits[1] == 0:
            sortWord2 = 'Real'
        else:
            sortWord2 = 'Real/Imaginary'
        if   self.sort_bits[2] == 0:
            sortWord3 = 'Sorted Responses'
        else:
            sortWord3 = 'Random Responses'

        #if   self.sort_code==0: sortWord = 'Real'
        #elif self.sort_code==1: sortWord = 'Real/Imaginary'
        #elif self.sort_code==2: sortWord = 'Random Responses'
        #else:
            #sortWord = '???'
            #msg = 'unsupported sort_code:  sort_code=%s\n' %(sort_code)
            #print msg
            #raise RuntimeError(msg)

        if thermal == 0:
            thermalWord = 'isHeatTransfer = False'
        elif thermal == 1:
            thermalWord = 'isHeatTransfer = True'
        elif thermal == 2:
            thermalWord = 'Scaled response spectra ABS'
        elif thermal == 3:
            thermalWord = 'Scaled response spectra SRSS'
        elif thermal == 4:
            thermalWord = 'Scaled response spectra NRL'
        elif thermal == 5:
            thermalWord = 'Scaled response spectra NRLO'
        else:
            thermalWord = '???'
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
            ForceFlux = 'Force'
        elif thermal == 1:
            ForceFlux = 'Flux'
        else:
            ForceFlux = 'Force (or Flux); thermal=%r' % thermal

        if thermal == 0:
            DispTemp = 'Displacement'
        elif thermal == 1:
            DispTemp = 'Temperature'
        #elif thermal is None:
            #raise RuntimeError('thermal_code is not specified; thermal_code=None')
        else:
            DispTemp = 'Displacement/Temperature; thermal=%r' % thermal

        table = '???'
        if self.table_code == 1:
            table = "OUG - %s vector/scalar" % DispTemp
        elif self.table_code == 2:
            table = "OPG - Load vector"
        elif self.table_code == 3:
            table = "OQG - SPC Force vector"
        elif self.table_code == 4:
            table = "OEF - Element %s" % ForceFlux
        elif self.table_code == 5:
            table = "OES - Element %s" % stressWord
        elif self.table_code == 6:
            table = "LAMA - Eigenvalue summary"
        elif self.table_code == 7:
            table = "OUG - Eigenvector"
        elif self.table_code == 8:
            table = "none - Grid point singularity table (obsolete)"
        elif self.table_code == 9:
            table = "OEIGS - Eigenvalue analysis summary"
        elif self.table_code == 10:
            table = "OUG - Velocity vector"
        elif self.table_code == 11:
            table = "OUG - Acceleration vector"
        elif self.table_code == 12:
            table = "OPG - Nonlinear force vector"
        elif self.table_code == 13:
            table = "OGPWG - Grid point weight generator"
        elif self.table_code == 14:
            table = "OUG - Eigenvector (solution set)"
        elif self.table_code == 15:
            table = "OUG - Displacement vector (solution set)"
        elif self.table_code == 16:
            table = "OUG - Velocity vector (solution set)"
        elif self.table_code == 17:
            table = "OUG - Acceleration vector (solution set)"
        elif self.table_code == 18:
            table = "OEE - Element strain energy"
        elif self.table_code == 19:
            table = "OGF - Grid point force balance"
        elif self.table_code == 20:
            table = "OES - Stresses at grid points (from the CURV module)"
        elif self.table_code == 21:
            table = "OES - Strain/curvature at grid points"
        elif self.table_code == 22:
            table = "OELOF1 - Element internal forces and moments"
        elif self.table_code == 23:
            table = "OELOP1 - Summation of element oriented forces on adjacent elements"
        elif self.table_code == 24:
            table = "OEP - Element pressures"
        elif self.table_code == 25:
            table = "OEF - Composite failure indicies"
        elif self.table_code == 26:
            table = "OGS - Grid point stresses (surface)"
        elif self.table_code == 27:
            table = "OGS - Grid point stresses (volume -- direct)"
        elif self.table_code == 28:
            table = "OGS - Grid point stresses (volume -- principal)"
        elif self.table_code == 29:
            table = "OGS - Element stress discontinuities (surface)"
        elif self.table_code == 30:
            table = "OGS - Element stress discontinuities (volume -- direct)"
        elif self.table_code == 31:
            table = "OGS - Element stress discontinuities (volume -- principal)"
        elif self.table_code == 32:
            table = "OGS - Grid point stress discontinuities (surface)"
        elif self.table_code == 33:
            table = "OGS - Grid point stress discontinuities (volume -- direct)"
        elif self.table_code == 34:
            table = "OGS - Grid point stress discontinuities (volume -- principal)"
        elif self.table_code == 35:
            table = "OGS - Grid point stress discontinuities (plane strain)"
        elif self.table_code == 36:
            table = "OEE - Element kinetic energy"
        elif self.table_code == 37:
            table = "OEE - Element energy loss"
        elif self.table_code == 38:
            table = "OMM - Max/Min summary"
        elif self.table_code == 39:
            table = "OQG - MPC Forces"
        elif self.table_code == 40:
            table = "OGPKE - Grip point kinetic energy"

        msg = '--Table3Data--\n\n'
        msg += "  device_code   = %-3s %s\n" % (self.device_code, device)
        msg += "  analysis_code = %-3s %s\n" % (self.analysis_code, analysis)
        msg += "  table_code    = %-3s %s-%s\n" % (self.table_code, self.table_name, table)
        msg += "  format_code   = %-3s %s\n" % (format_code, formatWord)

        msg += "  sortType      = %-3s %s\n" % (self.sort_bits[0], sortWord1)
        msg += "  dataFormat    = %-3s %s\n" % (self.sort_bits[1], sortWord2)
        msg += "  isRandom      = %-3s %s\n" % (self.sort_bits[2], sortWord3)

        if element_type is not None:
            msg += "  element_type  = %-3s %s\n" % (element_type, self.get_element_type(element_type))
        if sWord:  # stress code
            msg += "  s_code        = %-3s %s\n" % (s_code, sWord)
        if thermal is not None:
            msg += "  thermal       = %-3s %s\n" % (thermal, thermalWord)

        if hasattr(self, 'num_wide'):
            msg += "  num_wide      = %-3s\n" % self.num_wide
        if hasattr(self, 'isubcase'):
            msg += "  isubcase      = %-3s\n" % self.isubcase
        else:
            msg += "  ID            = %-3s\n" % self.ID
        #print msg
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

    #----
    # sort_code 0
    def is_sort1(self):
        if self.sort_bits[0] == 0:
            return True
        return False

    def is_sort2(self):
        return not self.is_sort1()

    #----
    # sort_code 1
    def isReal(self):  # format_code=1, this one is tricky b/c you can overwrite the Real code
        #if self.format_code==1:
        #    return True
        if self.sort_bits[1] == 0:
            return True
        return False

    def is_real_imaginary(self):  # format_code=2...does that dominate?
        return not self.isReal()

    #----
    # sort_code 2
    def isSortedResponse(self):
        if self.sort_bits[2] == 0:
            return True
        return False

    def isRandomResponse(self):
        return not self.isSortedResponse()

    #----
    # combos
    #def isRealOrRandom(self):  # been broken for a long time
        #asfd
        #return self.isReal() or self.isRandom()

    def isRealImaginaryOrMagnitudePhase(self):
        return self.is_real_imaginary or self.MagnitudePhase()

    #----
    def isStress(self):
        if self.stress_bits[1] == 0:
            return True
        return False

    def _set_op2_date(self, month, day, year):
        self.date = (month, day, 2000 + year)
        return self.date
