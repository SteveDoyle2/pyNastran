from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)


class Op2Codes(object):
    def ElementType(self, eCode):
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

    def printTableCode(self, tableCode):
        tableCodeContent = tableCode % 1000
        #dataFormat = tableCode/1000
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
        msg += 'n=%s table=%s-%s' % (self.n, self.tableName,
                                     tableContent[tableCodeContent])
        return msg

    def codeInformation(self):
        """
        prints the general table information
        DMAP - page 60-63
        """
        deviceCode = self.deviceCode
        #analysisCode = self.analysisCode
        #tableCode    = self.tableCode
        sortCode = self.sortCode

        formatCode = None
        if hasattr(self, 'formatCode'):
            formatCode = self.formatCode

        sCode = None
        if hasattr(self, 'sCode'):
            sCode = self.sCode

        thermal = None
        if hasattr(self, 'thermal'):
            thermal = self.thermal

        stressWord = ''
        if hasattr(self, 'stressBits'):
            if self.isStress():
                stressWord = 'Stress'
            else:
                stressWord = 'Strain'
            ###
        ###

        elementType = None
        if hasattr(self, 'elementType'):
            elementType = self.elementType

        sWord = ''
        if(sCode == 0):
            sWord += 'Coordinate Element - Stress Max Shear (Octahedral)'
        elif(sCode == 14):
            sWord += 'Coordinate Element - Strain Fiber Max Shear (Octahedral)'

        elif(sCode == 1):
            sWord += 'Coordinate Element - Stress von Mises'
        elif(sCode == 10):
            sWord += 'Coordinate Element - Strain Curvature Max Shear (Octahedral)'

        elif(sCode == 11):
            sWord += 'Coordinate Element - Strain Curvature von Mises'
        elif(sCode == 15):
            sWord += 'Coordinate Element - Strain Fiber von Mises'

        elif(sCode == 16):
            sWord += 'Coordinate Material - Stress Max Shear (Octahedral)'
        elif(sCode == 17):
            sWord += 'Coordinate Material - Stress von Mises'

        elif(sCode == 26):
            sWord += 'Coordinate Material - Strain Curvature Max Shear'
        elif(sCode == 30):
            sWord += 'Coordinate Material - Strain Fiber Max Shear (Octahedral)'

        elif(sCode == 27):
            sWord += 'Coordinate Material - Strain Curvature von Mises'
        elif(sCode == 31):
            sWord += 'Coordinate Material - Strain Fiber von Mises'
        else:
            #sWord = 'Stress or Strain - UNDEFINED'
            sWord = ''

        formatWord = '???'
        if(formatCode == 1):
            formatWord = "Real"
        elif(formatCode == 2):
            formatWord = "Real/Imaginary"
        elif(formatCode == 3):
            formatWord = "Magnitude/Phase"
        else:
            formatWord = '???'
            #msg = 'unsupported formatCode:  formatCode=%s\n' %(formatCode)
            #raise InvalidFormatCodeError(msg)

        if   self.sortBits[0] == 0:
            sortWord1 = 'Sort1'
        else:
            sortWord1 = 'Sort2'
        if   self.sortBits[1] == 0:
            sortWord2 = 'Real'
        else:
            sortWord2 = 'Real/Imaginary'
        if   self.sortBits[2] == 0:
            sortWord3 = 'Sorted Responses'
        else:
            sortWord3 = 'Random Responses'

        #if(  self.sortCode==0): sortWord = 'Real'
        #elif(self.sortCode==1): sortWord = 'Real/Imaginary'
        #elif(self.sortCode==2): sortWord = 'Random Responses'
        #else:
            #sortWord = '???'
            #msg = 'unsupported sortCode:  sortCode=%s\n' %(sortCode)
            #print msg
            #raise RuntimeError(msg)

        if(thermal == 0):
            thermalWord = 'isHeatTransfer = False'
        elif(thermal == 1):
            thermalWord = 'isHeatTransfer = True'
        elif(thermal == 2):
            thermalWord = 'Scaled response spectra ABS'
        elif(thermal == 3):
            thermalWord = 'Scaled response spectra SRSS'
        elif(thermal == 4):
            thermalWord = 'Scaled response spectra NRL'
        elif(thermal == 5):
            thermalWord = 'Scaled response spectra NRLO'
        else:
            thermalWord = '???'
            #msg = 'unsupported thermal:  thermal=%s\n' %(thermal)
            #raise ValueError(msg)

        analysis = '???'
        if(self.analysisCode == 1):
            analysis = "Statics"
        elif(self.analysisCode == 2):
            analysis = "Normal modes or buckling (real eigenvalues)"
        elif(self.analysisCode == 3):
            analysis = "Differential Stiffness 0 - obsolete"
        elif(self.analysisCode == 4):
            analysis = "Differential Stiffness 1 - obsolete"
        elif(self.analysisCode == 5):
            analysis = "Frequency"
        elif(self.analysisCode == 6):
            analysis = "Transient"
        elif(self.analysisCode == 7):
            analysis = "Pre-buckling"
        elif(self.analysisCode == 8):
            analysis = "Post-buckling"
        elif(self.analysisCode == 9):
            analysis = "Complex eigenvalues"
        elif(self.analysisCode == 10):
            analysis = "Nonlinear statics"
        elif(self.analysisCode == 11):
            analysis = "Geometric nonlinear statics"

        device = '???'
        if(self.deviceCode == 1):
            device = "Print"
        elif(self.deviceCode == 2):
            device = "Plot"
        elif(self.deviceCode == 3):
            device = "Print and Plot"
        elif(self.deviceCode == 4):
            device = "Punch"
        elif(self.deviceCode == 5):
            device = "Print and Punch"
        elif(self.deviceCode == 6):
            device = "Plot and Punch"
        elif(self.deviceCode == 7):
            device = "Print, Plot, and Punch"

        if thermal == 0:
            ForceFlux = 'Force'
        elif thermal == 1:
            ForceFlux = 'Flux'
        else:
            ForceFlux = 'Force (or Flux)'

        if thermal == 0:
            DispTemp = 'Displacement'
        elif thermal == 1:
            DispTemp = 'Temperature'
        else:
            DispTemp = 'Displacement/Temperature'

        table = '???'
        if(self.tableCode == 1):
            table = "OUG - %s vector/scalar" % (DispTemp)
        elif(self.tableCode == 2):
            table = "OPG - Load vector"
        elif(self.tableCode == 3):
            table = "OQG - SPC Force vector"
        elif(self.tableCode == 4):
            table = "OEF - Element %s" % (ForceFlux)
        elif(self.tableCode == 5):
            table = "OES - Element %s" % (stressWord)
        elif(self.tableCode == 6):
            table = "LAMA - Eigenvalue summary"
        elif(self.tableCode == 7):
            table = "OUG - Eigenvector"
        elif(self.tableCode == 8):
            table = "none - Grid point singularity table (obsolete)"
        elif(self.tableCode == 9):
            table = "OEIGS - Eigenvalue analysis summary"
        elif(self.tableCode == 10):
            table = "OUG - Velocity vector"
        elif(self.tableCode == 11):
            table = "OUG - Acceleration vector"
        elif(self.tableCode == 12):
            table = "OPG - Nonlinear force vector"
        elif(self.tableCode == 13):
            table = "OGPWG - Grid point weight generator"
        elif(self.tableCode == 14):
            table = "OUG - Eigenvector (solution set)"
        elif(self.tableCode == 15):
            table = "OUG - Displacement vector (solution set)"
        elif(self.tableCode == 16):
            table = "OUG - Velocity vector (solution set)"
        elif(self.tableCode == 17):
            table = "OUG - Acceleration vector (solution set)"
        elif(self.tableCode == 18):
            table = "OEE - Element strain energy"
        elif(self.tableCode == 19):
            table = "OGF - Grid point force balance"
        elif(self.tableCode == 20):
            table = "OES - Stresses at grid points (from the CURV module)"
        elif(self.tableCode == 21):
            table = "OES - Strain/curvature at grid points"
        elif(self.tableCode == 22):
            table = "OELOF1 - Element internal forces and moments"
        elif(self.tableCode == 23):
            table = "OELOP1 - Summation of element oriented forces on adjacent elements"
        elif(self.tableCode == 24):
            table = "OEP - Element pressures"
        elif(self.tableCode == 25):
            table = "OEF - Composite failure indicies"
        elif(self.tableCode == 26):
            table = "OGS - Grid point stresses (surface)"
        elif(self.tableCode == 27):
            table = "OGS - Grid point stresses (volume -- direct)"
        elif(self.tableCode == 28):
            table = "OGS - Grid point stresses (volume -- principal)"
        elif(self.tableCode == 29):
            table = "OGS - Element stress discontinuities (surface)"
        elif(self.tableCode == 30):
            table = "OGS - Element stress discontinuities (volume -- direct)"
        elif(self.tableCode == 31):
            table = "OGS - Element stress discontinuities (volume -- principal)"
        elif(self.tableCode == 32):
            table = "OGS - Grid point stress discontinuities (surface)"
        elif(self.tableCode == 33):
            table = "OGS - Grid point stress discontinuities (volume -- direct)"
        elif(self.tableCode == 34):
            table = "OGS - Grid point stress discontinuities (volume -- principal)"
        elif(self.tableCode == 35):
            table = "OGS - Grid point stress discontinuities (plane strain)"
        elif(self.tableCode == 36):
            table = "OEE - Element kinetic energy"
        elif(self.tableCode == 37):
            table = "OEE - Element energy loss"
        elif(self.tableCode == 38):
            table = "OMM - Max/Min summary"
        elif(self.tableCode == 39):
            table = "OQG - MPC Forces"
        elif(self.tableCode == 40):
            table = "OGPKE - Grip point kinetic energy"

        msg = '--Table3Data--\n\n'
        msg += "  deviceCode   = %-3s %s\n" % (self.deviceCode, device)
        msg += "  analysisCode = %-3s %s\n" % (self.analysisCode, analysis)
        msg += "  tableCode    = %-3s %s-%s\n" % (
            self.tableCode, self.tableName, table)
        msg += "  formatCode   = %-3s %s\n" % (formatCode, formatWord)

        msg += "  sortType     = %-3s %s\n" % (self.sortBits[0], sortWord1)
        msg += "  dataFormat   = %-3s %s\n" % (self.sortBits[1], sortWord2)
        msg += "  isRandom     = %-3s %s\n" % (self.sortBits[2], sortWord3)

        if elementType is not None:
            msg += "  elementType  = %-3s %s\n" % (
                elementType, self.ElementType(elementType))
        if sWord:  # stress code
            msg += "  sCode        = %-3s %s\n" % (sCode, sWord)
        if thermal is not None:
            msg += "  thermal      = %-3s %s\n" % (thermal, thermalWord)
        msg += "  numWide      = %-3s\n" % (self.numWide)
        msg += "  iSubcase     = %-3s\n" % (self.iSubcase)
        #print msg
        return msg

    #----
    def isThermal(self):
        if self.thermal == 0:
            return False
        elif self.thermal == 1:
            return True
        return '???'

    #----
    # formatCode 3
    def isMagnitudePhase(self):
        if self.formatCode == 3:
            return True
        return False

    #----
    # sortCode 0
    def isSort1(self):
        if self.sortBits[0] == 0:
            return True
        return False

    def isSort2(self):
        return not(self.isSort1())

    #----
    # sortCode 1
    def isReal(self):  # formatCode=1, this one is tricky b/c you can overwrite the Real code
        #if self.formatCode==1:
        #    return True
        if self.sortBits[1] == 0:
            return True
        return False

    def isRealImaginary(self):  # formatCode=2...does that dominate?
        return not(self.isReal())

    #----
    # sortCode 2
    def isSortedResponse(self):
        if self.sortBits[2] == 0:
            return True
        return False

    def isRandomResponse(self):
        return not(self.isSortedResponse())

    #----
    # combos
    def isRealOrRandom(self):
        return self.isReal() or self.isRandom()

    def isRealImaginaryOrMagnitudePhase(self):
        return self.isRealImaginary or self.MagnitudePhase()

    #----
    def isStress(self):
        if self.stressBits[1] == 0:
            return True
        return False
