from struct import Struct, unpack

from pyNastran import isRelease
from pyNastran.op2.op2_helper import polar_to_real_imag


class OP2Common(object):
    def __init__(self):
        self.show_table3_map = [
            #'OUGV1',
            #'OEF1X',
            #'OES1X1',
        ]
        self.show_table4_map = [
            #'OUGV1',
            #'OEF1X',
            #'OES1X1',
        ]

        self.element_mapper = {
            # rods
            1 : 'CROD',
            3 : 'CTUBE',
            10: 'CONROD',

            # beam
            2 : 'CBEAM',

            # springs
            11: 'CELAS1',
            12: 'CELAS2',
            13: 'CELAS3',
            14: 'CELAS4',

            # bar
            34 : 'CBAR',

            # solids
            39: 'CTETRA',
            67: 'CHEXA',
            68: 'CPENTA',

            # triax
            53 : 'CTRIAX6',

            # centroidal plate
            33: 'CQUAD4-centroidal', # 1 node
            74: 'CTRIA3-centroidal', # 1 node

            # bilinear plate
            144: 'CQUAD4-bilinear', # 5 nodes

            # composite plates
            95 : 'CQUAD4-composite',
            96 : 'CQUAD8-composite',
            97 : 'CTRIA3-composite',
            98 : 'CTRIA6-composite',

            # unorganized
            None: '',
            0: 'GRID',
            #1: 'ROD',
            #2: 'BEAM',
            #3: 'TUBE',
            4: 'SHEAR',
            5: 'FORMON12',
            6: 'FORCE',
            7: 'PLOAD4',
            8: 'PLOADX1',
            9: 'PLOAD/PLOAD2',
            #10: 'CONROD',
            #11: 'ELAS1',
            #12: 'ELAS2',
            #13: 'ELAS3',
            #14: 'ELAS4',
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
            #33: 'QUAD4',
            #34: 'BAR',
            35: 'CON',
            36: None,
            37: None,
            38: 'GAP',
            #39: 'TETRA',
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
            #53: 'TRIAX6',
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
            #67: 'HEXA',
            #68: 'PENTA',
            69: 'BEND',
            70: 'TRIAR',
            71: None,
            72: 'AEROQ4',
            73: None,
            #74: 'TRIA3',
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
            #95: 'QUAD4LC',
            #96: 'QUAD8LC',
            #97: 'TRIA3LC',
            #98: 'TRIA6LC',
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
            #144: 'QUAD144',
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
        pass

    def add_data_parameter(self, data, var_name, Type, field_num,
                           applyNonlinearFactor=True, fixDeviceCode=False, add_to_dict=True):

        datai = data[4*(field_num-1) : 4*(field_num)]
        assert len(datai) == 4, len(datai)
        value, = unpack(Type, datai)
        #print "%-12s = %r" % (var_name, value)
        if self.debug:
            self.binary_debug.write('  %-12s = %r\n' % (var_name, value))
        #setattr(self, var_name, value)  # set the parameter to the local namespace
        self.words[field_num-1] = var_name
        return value

    def apply_data_code_value(self, name, value):
        pass
    def setNullNonlinearFactor(self):
        pass

    def read_title(self, data):
        assert len(data) == 584, len(data)
        Title, subtitle, label = unpack(b'128s128s128s', data[200:])  # titleSubtitleLabel

        self.Title = Title.strip()
        #: the subtitle of the subcase
        self.subtitle = subtitle.strip()
        #: the label of the subcase
        self.label = label.strip()
        if self.debug:
            self.binary_debug.write('  Title    = %r\n' % self.Title)
            self.binary_debug.write('  subtitle = %r\n' % self.subtitle)
            self.binary_debug.write('  label    = %r\n' % self.label)

    def write_debug_bits(self):
        if self.debug:
            msg = ''
            for i, param in enumerate(self.words):
                if param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '
            if hasattr(self, 'format_code'):
                self.binary_debug.write('  sort_bits[0] = %i -> is_sort1 =%s\n' % (self.sort_bits[0], self.is_sort1() ))
                self.binary_debug.write('  sort_bits[1] = %i -> is_real  =%s vs real/imag\n' % (self.sort_bits[1], self.is_real()   ))
                self.binary_debug.write('  sort_bits[2] = %i -> is_random=%s vs mag/phase\n' % (self.sort_bits[2], self.is_random() ))
                if self.is_complex():
                    self.binary_debug.write('  format_code  = %i -> is_mag_phase=%s vs is_real_imag\n' % (self.format_code, self.is_mag_phase() ))
                else:
                    self.binary_debug.write('  format_code  = %i\n' % self.format_code)
            self.binary_debug.write('  recordi = [%s]\n\n' % msg)

    def read_oug_table(self, data, result_name, real_obj, complex_obj,
                       thermal_real_obj, node_elem):
        assert real_obj is None
        assert complex_obj is None

        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                # real_obj
                self.read_real_table(data, result_name, node_elem)
            else:
                # thermal_real_obj
                self.read_real_table(data, result_name, node_elem)
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                # complex_obj
                self.read_complex_table(data, result_name, node_elem)
            else:
                self.not_implemented_or_skip('thermal=%r' % self.thermal)
        else:
            msg = 'only num_wide=8 or 14 is allowed  num_wide=%s' % self.num_wide
            self.not_implemented_or_skip(msg)

    def read_real_table(self, data, result_name, flag):
        #return
        if self.debug4():
            self.binary_debug.write('  read_real_table\n')
        assert flag in ['node', 'elem'], flag
        format1 = '2i6f' # 8

        ntotal = 32 # 8 * 4
        nnodes = len(data) // ntotal > 0
        assert len(data) % ntotal == 0

        n = 0
        s = Struct(format1)
        for inode in xrange(nnodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid_device, gridType, tx, ty, tz, rx, ry, rz) = out

            eid = (eid_device - self.device_code) // 10
            if self.debug4():
                self.binary_debug.write('  %s=%i; %s\n' % (flag, eid_device, str(out)))
            #print "eType=%s" %(eType)

            dataIn = [eid, gridType, tx, ty, tz, rx, ry, rz]
            #self.obj.add(dt, dataIn)
            n += ntotal

    def read_complex_table(self, data, result_name, flag):
        #return
        if self.debug4():
            self.binary_debug.write('  read_complex_table\n')
        assert flag in ['node', 'elem'], flag

        format1 = '2i12f'
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 56  # 14 * 4
        nnodes = len(data) // ntotal
        s = Struct(format1)
        for inode in xrange(nnodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)

            (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
             txi, tyi, tzi, rxi, ryi, rzi) = out
            eid = (eid_device - self.device_code) // 10

            if self.debug4():
                self.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
            if is_magnitude_phase:
                tx = polar_to_real_imag(txr, txi)
                rx = polar_to_real_imag(rxr, rxi)
                ty = polar_to_real_imag(tyr, tyi)
                ry = polar_to_real_imag(ryr, ryi)
                tz = polar_to_real_imag(tzr, tzi)
                rz = polar_to_real_imag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                rx = complex(rxr, rxi)
                ty = complex(tyr, tyi)
                ry = complex(ryr, ryi)
                tz = complex(tzr, tzi)
                rz = complex(rzr, rzi)

            dataIn = [eid, grid_type, tx, ty, tz, rx, ry, rz]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            #self.obj.add(dt, dataIn)
            n += ntotal

    def not_implemented_or_skip(self, msg=''):
        if isRelease:
            pass
        else:
            raise NotImplementedError('table_name=%s table_code=%s %s' % (self.table_name, self.table_code, msg))

    def parse_approach_code(self, data):
        (aCode, tCode, int3, isubcase) = unpack(b'4i', data[:16])
        self.aCode = aCode
        self.tCode = tCode
        self.int3 = int3

        #: the local subcase ID
        self.isubcase = isubcase
        #print("isubcase = %s" %(isubcase))
        #self.subcases.add(self.isubcase)  # set notation

        #: the type of result being processed
        self.table_code = tCode % 1000

        #: used to create sort_bits
        self.sort_code = tCode // 1000
        self.sort_code2 = ((tCode // 1000) + 2) // 2

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        self.device_code = aCode % 10

        #: what solution was run (e.g. Static/Transient/Modal)
        self.analysis_code = (aCode - self.device_code) // 10

        #print('parse_approach_code - aCode=%s tCode=%s int3=%s isubcase=%s' % (aCode, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (self.analysis_code, self.device_code, self.table_code, self.sort_code))
        if self.debug3():
            self.binary_debug.write('  table_code    = %r\n' % self.table_code)
            self.binary_debug.write('  sort_code     = %r\n' % self.sort_code)
            self.binary_debug.write('  sort_code2    = %r\n' % self.sort_code2)
            self.binary_debug.write('  device_code   = %r\n' % self.device_code)
            self.binary_debug.write('  analysis_code = %r\n' % self.analysis_code)

        if self.device_code == 3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results..."
            #                 "if there's a crash, try plot only\n")
            self.device_code = 1

            #self.log.info('The op2 may be inconsistent...')
            #self.log.info('  print and plot can cause bad results...'
            #              'if there's a crash, try plot only')

        self._parse_sort_code()

    #===================================
    def _parse_sort_code(self):
        """
        +-------------------------+
        | sort_code  | sort_bits  |
        +=========================+
        | 0          | [0, 0, 0]  |
        +------------+------------+
        | 1          | [0, 0, 1]  |
        +------------+------------+
        | 2          | [0, 1, 0]  |
        +------------+------------+
        | 3          | [0, 1, 1]  |
        +------------+------------+
        | ...        | ...        |
        +------------+------------+
        | 7          | [1, 1, 1]  |
        +------------+------------+

        sort_code = 0 -> sort_bits = [0,0,0]
        sort_code = 1 -> sort_bits = [0,0,1]
        sort_code = 2 -> sort_bits = [0,1,0]
        sort_code = 3 -> sort_bits = [0,1,1]
        etc.
        sort_code = 7 -> sort_bits = [1,1,1]

        sort_bits[0] = 0 -> is_sort1=True  isSort2=False
        sort_bits[1] = 0 -> isReal=True   isReal/Imaginary=False
        sort_bits[2] = 0 -> isSorted=True isRandom=False
        """
        bits = [0, 0, 0]
        sort_code = self.sort_code

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        if self.sort_code not in [0, 1]:
            msg = 'Invalid sort_code=%s sort_code2=%s' %(self.sort_code, self.sort_code2)
            raise RuntimeError(msg)

        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1
        #: the bytes describe the SORT information
        self.sort_bits = bits
        #print "sort_bits =", bits
        #self.data_code['sort_bits'] = self.sort_bits

    def is_sort1(self):
        if self.sort_bits[0] == 0:
            return True
        return False

    def is_sort2(self):
        return not(self.is_sort1())

    def is_real(self):
        return not(self.is_complex())

    def is_complex(self):
        if self.sort_bits[1] == 1:
            return True
        return False

    def is_random(self):
        if self.sort_bits[1] == 1:
            return True
        return False

    def is_mag_phase(self):
        assert self.format_code in [0, 1], self.format_code
        return bool(self.format_code)

    def is_magnitude_phase(self):
        if self.format_code == 3:
            return True
        return False

    def debug3(self):
        return True
        if self.debug and self.table_name in self.show_table3_map:
            return True
        return False

    def debug4(self):
        return True
        if self.debug and self.table_name in self.show_table4_map:
            return True
        return False