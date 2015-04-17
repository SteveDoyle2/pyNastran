from __future__ import print_function
from six import string_types
from six.moves import range
import copy
from struct import Struct, unpack

from pyNastran import is_release
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.utils import object_attributes

from pyNastran.f06.f06Writer import F06Writer
from pyNastran.op2.op2Codes import Op2Codes


class SortCodeError(RuntimeError):
    pass

class OP2Common(Op2Codes, F06Writer):
    def __init__(self):
        Op2Codes.__init__(self)
        F06Writer.__init__(self)

        #: flag for vectorization
        #: 0 - no vectorization
        #: 1 -   first pass
        #: 2 -   second pass
        self.read_mode = None

        #: the results
        self.result_names = set([])
        #: bool
        self.is_vectorized = None

        #: the storage dictionary that is passed to OP2 objects (e.g. DisplacementObject)
        #: the key-value pairs are extracted and used to generate dynamic self
        #: variables for the OP2 objects
        self.data_code = {}

        #: current subcase ID
        #: non-transient (SOL101) cases have isubcase set to None
        #: transient (or frequency/modal) cases have isubcase set to a int/float value
        self.isubcase = None

        #: the corresponding piece to isubcase
        #: used only for SORT2 (not supported)
        self.ID = None

        #: should the op2 debugging file be written
        self.debug = False

        #: op2 debug file or None (for self.debug=False)
        self.binary_debug = None

        #: the list of "words" on a subtable 3
        self.words = []

        #: The current table_name (e.g. OES1)
        #: None indicates no table_name has been read
        self.table_name = None

        # the date stamp used in the F06
        self.date = (1, 1, 2000)

        #: set of all the subcases that have been found
        self.subcases = set()

        #: the list/set/tuple of times/modes/frequencies that should be read
        #: currently unused
        self.expected_times = None

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
            33: 'CQUAD4', # 1 node
            74: 'CTRIA3', # 1 node

            # bilinear plate
            144: 'CQUAD4', # 5 nodes

            # composite plates
            95 : 'CQUAD4',
            96 : 'CQUAD8',
            97 : 'CTRIA3',
            98 : 'CTRIA6',

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
            64: 'CQUAD8',
            65: None,
            66: None,
            #67: 'HEXA',
            #68: 'PENTA',
            69: 'CBEND',
            70: 'CTRIAR',
            71: None,
            72: 'AEROQ4',
            73: None,
            #74: 'TRIA3',
            75: 'CTRIA6',
            76: 'HEXPR',
            77: 'PENPR',
            78: 'TETPR',
            79: None,
            80: None,
            81: None,
            82: 'CQUADR',
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
        if self.debug:
            self.binary_debug.write('  %-12s = %r\n' % (var_name, value))
        #setattr(self, var_name, value)  # set the parameter to the local namespace

        if applyNonlinearFactor:
            self.nonlinear_factor = value
            self.data_code['nonlinear_factor'] = value
            self.data_code['name'] = var_name

        if add_to_dict:
            self.data_code[var_name] = value

        try:
            self.words[field_num - 1] = var_name
        except IndexError:
            msg = 'Trying to set word, but...len(words)=%s ifield=%s' % (len(self.words), field_num - 1)
            raise IndexError(msg)
        return value

    def apply_data_code_value(self, name, value):
        self.data_code[name] = value

    def setNullNonlinearFactor(self):
        self.nonlinear_factor = None
        self.data_code['nonlinear_factor'] = None

    def _read_title_helper(self, data):
        assert len(data) == 584, len(data)
        # titleSubtitleLabel
        Title, subtitle, label = unpack(b'128s128s128s', data[200:])

        self.Title = Title.strip()

        #: the subtitle of the subcase
        self.subtitle = subtitle.strip()
        self.data_code['subtitle'] = self.subtitle

        #: the label of the subcase
        self.label = label.strip()
        self.data_code['label'] = self.label
        self.data_code['Title'] = self.Title

        if self.debug:
            self.binary_debug.write('  Title    = %r\n' % self.Title)
            self.binary_debug.write('  subtitle = %r\n' % self.subtitle)
            self.binary_debug.write('  label    = %r\n' % self.label)

    def _read_title(self, data):
        self._read_title_helper(data)

        if hasattr(self, 'isubcase'):
            if self.isubcase not in self.iSubcaseNameMap:
                self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.label]
        else:
            raise  RuntimeError('isubcase is not defined')

        if hasattr(self, 'subtitle') and hasattr(self, 'label'):
            if (self.isubcase, self.subtitle) not in self.labels:
                self.subtitles[self.isubcase].append(self.subtitle)
                self.labels[(self.isubcase, self.subtitle)] = self.label

    def _write_debug_bits(self):
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

    def _read_geom_4(self, mapper, data):
        if not self.make_geom:
            return len(data)
        n = 0
        keys = unpack('3i', data[n:n+12])
        n += 12
        if len(data) == 12:
            #print('*self.istream = %s' % self.istream)
            #print('self.isubtable = %s' % self.isubtable)
            self.istream -= 1
            self.isubtable_old = self.isubtable
            return n

        #print('is_start_of_subtable=%s' % self.is_start_of_subtable)
        #print('self.istream = %s' % self.istream)
        if not hasattr(self, 'isubtable_old'):
            self.isubtable_old = None
        elif self.isubtable_old > self.isubtable:
            self.isubtable_old = None

        #self.binary_debug.write('isubtable=%s isubtable_old=%s\n' % (self.isubtable, self.isubtable_old))
        #ni = self.f.tell() - len(data) + 12
        #self.binary_debug.write('**:  f.tell()=%s; n=%s:%s\n\n' % (self.f.tell(), ni, self.n))

        # we're only going to use the keys if istream=0 (so the beginning of the record)
        if self.istream == 0 and keys in mapper:
            pass
        elif self.isubtable_old == self.isubtable:
            # we didn't increment the record, so we fix the n+=12 statement we called before
            # then we toss the keys and use the old geom_keys
            n = 0
            keys = self.geom_keys
        else:
            raise NotImplementedError('keys=%s not found - %s; istream=%s; isubtable=%s isubtable_old=%s' % (str(keys), self.table_name, self.istream, self.isubtable, self.isubtable_old))

        name, func = mapper[keys]
        self.binary_debug.write('  found keys=%s -> name=%-6s - %s\n' % (str(keys), name, self.table_name))
        #print("  found keys=(%5s,%4s,%4s) name=%-6s - %s" % (keys[0], keys[1], keys[2], name, self.table_name))

        n = func(data, n)  # gets all the grid/mat cards
        assert n != None, name

        self.geom_keys = keys
        self.is_start_of_subtable = False
        self.isubtable_old = self.isubtable

        #assert n == len(data), 'n=%s len(data)=%s' % (n, len(data))
        return n

    def _read_table(self, data, result_name, storage_obj,
                    real_obj, complex_obj,
                    real_vector, complex_vector,
                    node_elem, random_code=None, is_cid=False):

        assert isinstance(result_name, string_types), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #assert real_obj is None
        #assert complex_obj is None
        #assert thermal_real_obj is None

        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real_obj
            assert real_obj is not None
            nnodes = len(data) // 32  # 8*4
            auto_return = self._create_table_object(result_name, nnodes, storage_obj, real_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return len(data)
            n = self._read_real_table(data, result_name, node_elem, is_cid=is_cid)
        elif self.format_code in [1, 2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            # complex_obj
            assert complex_obj is not None
            nnodes = len(data) // 56  # 14*4
            self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_object(result_name, nnodes, storage_obj, complex_obj, complex_vector)
            if auto_return:
                return len(data)
            n = self._read_complex_table(data, result_name, node_elem)
        elif self.format_code in [2, 3] and self.num_wide == 8:  # real mag/phase - this is confusing...I think there is no phase?
            # real_obj
            assert real_obj is not None
            nnodes = len(data) // 32  # 8*4
            auto_return = self._create_table_object(result_name, nnodes, storage_obj, real_obj, real_vector)
            if auto_return:
                return len(data)
            n = self._read_real_table(data, result_name, node_elem)
        else:
            msg = 'only num_wide=8 or 14 is allowed;  num_wide=%s' % self.num_wide
            n = self._not_implemented_or_skip(data, msg)
        #else:
        #msg = 'invalid random_code=%s num_wide=%s' % (random_code, self.num_wide)
        #n = self._not_implemented_or_skip(data, msg)
        return n

    def _read_real_table(self, data, result_name, flag, is_cid=False):
        if self.debug4():
            self.binary_debug.write('  _read_real_table\n')
        assert flag in ['node', 'elem'], flag
        n = 0
        ntotal = 32 # 8 * 4
        dt = self.nonlinear_factor
        assert self.obj is not None

        obj = self.obj
        format1 = '2i6f' # 8

        nnodes = len(data) // ntotal

        assert nnodes > 0
        #assert len(data) % ntotal == 0
        s = Struct(format1)
        for inode in range(nnodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out

            eid = (eid_device - self.device_code) // 10
            if self.debug4():
                self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
            obj.add(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
            n += ntotal
        return n

    def _read_complex_table(self, data, result_name, flag):
        if self.debug4():
            self.binary_debug.write('  _read_complex_table\n')
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor

        format1 = '2i12f'
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 56  # 14 * 4

        obj = self.obj
        nnodes = len(data) // ntotal
        s = Struct(format1)

        assert self.obj is not None
        assert nnodes > 0
        #assert len(data) % ntotal == 0

        if self.debug4():
            self.binary_debug.write('  nnodes=%i\n' % (nnodes))
        for inode in range(nnodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)

            (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
             txi, tyi, tzi, rxi, ryi, rzi) = out
            eid = (eid_device - self.device_code) // 10

            if self.debug4():
                self.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
            if is_magnitude_phase:
                tx = polar_to_real_imag(txr, txi)
                ty = polar_to_real_imag(tyr, tyi)
                tz = polar_to_real_imag(tzr, tzi)

                rx = polar_to_real_imag(rxr, rxi)
                ry = polar_to_real_imag(ryr, ryi)
                rz = polar_to_real_imag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
                tz = complex(tzr, tzi)

                rx = complex(rxr, rxi)
                ry = complex(ryr, ryi)
                rz = complex(rzr, rzi)

            obj.add(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
            n += ntotal
        return n

    def _read_complex_table2(self, data, result_name, flag):
        #return
        if self.debug4():
            self.binary_debug.write('  _read_complex_table\n')
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor

        format1 = 'fi12f'
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 56  # 14 * 4
        nnodes = len(data) // ntotal
        s = Struct(format1)

        assert self.obj is not None
        assert nnodes > 0
        #assert len(data) % ntotal == 0

        for inode in range(nnodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)

            (freq, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
             txi, tyi, tzi, rxi, ryi, rzi) = out

            if self.debug4():
                self.binary_debug.write('  %s=%i %s\n' % ('freq', freq, str(out)))
            if is_magnitude_phase:
                tx = polar_to_real_imag(txr, txi)
                ty = polar_to_real_imag(tyr, tyi)
                tz = polar_to_real_imag(tzr, tzi)

                rx = polar_to_real_imag(rxr, rxi)
                ry = polar_to_real_imag(ryr, ryi)
                rz = polar_to_real_imag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
                tz = complex(tzr, tzi)

                rx = complex(rxr, rxi)
                ry = complex(ryr, ryi)
                rz = complex(rzr, rzi)

            data_in = [freq, grid_type, tx, ty, tz, rx, ry, rz]
            self.obj.add(dt, data_in)
            n += ntotal
        return n

    def create_transient_object(self, storageObj, classObj, is_cid=False, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        :param storageName:  the name of the dictionary to store the object in (e.g. 'displacements')
        :param classObj:    the class object to instantiate
        :param debug:       developer debug

        .. note:: dt can also be load_step depending on the class
        """
        assert not isinstance(classObj, string_types), 'classObj=%r' % classObj
        if debug:
            print("create Transient Object")
            print("***NF = %s" % self.nonlinear_factor)
        #if not hasattr(self, storageName):
            #attrs =  object_attributes(obj, mode="public")
            #msg = 'storage_obj=%r does not exist.\n' % storageObj
            #msg += 'Attributes = [%s]' , ', %s'.join(attrs)
            #raise RuntimeError(msg)
        #storageObj = getattr(self, storageName)
        #assert classObj is not None, 'name=%r has no associated classObject' % storageName
        self.data_code['table_name'] = self.table_name
        assert self.log is not None

        code = self._get_code()
        if hasattr(self, 'isubcase'):
            if self.code in storageObj:
                self.obj = storageObj[code]
                self.obj.update_data_code(copy.deepcopy(self.data_code))
            else:
                classObj.is_cid = is_cid
                self.obj = classObj(self.data_code, self.is_sort1(), self.isubcase, self.nonlinear_factor)
            storageObj[code] = self.obj
        else:
            if code in storageObj:
                self.obj = storageObj[code]
            else:
                storageObj[code] = self.obj

    def _get_code(self):
        code = self.isubcase
        code = (self.isubcase, self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code

    def _not_implemented_or_skip(self, data, msg=''):
        if is_release:
            return len(data)
        else:
            raise NotImplementedError('table_name=%s table_code=%s %s\n%s' % (self.table_name, self.table_code, msg, self.code_information()))

    def parse_approach_code(self, data):
        (approach_code, tCode, int3, isubcase) = unpack(b'4i', data[:16])
        self.approach_code = approach_code
        self.tCode = tCode
        self.int3 = int3

        #: the local subcase ID
        self.isubcase = isubcase
        self.data_code['isubcase'] = self.isubcase
        #self.subcases.add(self.isubcase)  # set notation

        #: the type of result being processed
        self.table_code = tCode % 1000
        self.data_code['table_code'] = self.table_code

        #: used to create sort_bits
        self.sort_code = tCode // 1000

        #: .. todo::  sort_code2 seems to be unused...
        self.sort_code2 = ((tCode // 1000) + 2) // 2
        self.data_code['sort_code'] = self.sort_code

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        self.device_code = approach_code % 10
        self.data_code['device_code'] = self.device_code

        #: what solution was run (e.g. Static/Transient/Modal)
        self.analysis_code = (approach_code - self.device_code) // 10
        self.data_code['analysis_code'] = self.analysis_code

        #print('parse_approach_code - approach_code=%s tCode=%s int3=%s isubcase=%s' % (approach_code, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (self.analysis_code, self.device_code, self.table_code, self.sort_code))
        if self.device_code == 3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results..."
            #                 "if there's a crash, try plot only\n")
            self.device_code = 1

            #self.log.info('The op2 may be inconsistent...')
            #self.log.info('  print and plot can cause bad results...'
            #              'if there's a crash, try plot only')
            self.data_code['device_code'] = self.device_code

        if self.debug3():
            self.binary_debug.write('  table_name    = %r\n' % self.table_name)
            self.binary_debug.write('  approach_code = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode         = %r\n' % self.tCode)
            self.binary_debug.write('  table_code    = %r\n' % self.table_code)
            self.binary_debug.write('  sort_code     = %r\n' % self.sort_code)
            self.binary_debug.write('  sort_code2    = %r\n' % self.sort_code2)
            self.binary_debug.write('  device_code   = %r\n' % self.device_code)
            self.binary_debug.write('  analysis_code = %r\n' % self.analysis_code)

        self._parse_sort_code()

    def _parse_sort_code(self):
        """
        +------------+------------+
        | sort_code  | sort_bits  |
        +============+============+
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

        ::
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
            msg = 'Invalid sort_code=%s sort_code2=%s' % (self.sort_code, self.sort_code2)
            raise SortCodeError(msg)
        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1
        #: the bytes describe the SORT information
        self.sort_bits = bits
        #self.data_code['sort_bits'] = self.sort_bits

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
        |   5   |   SORT2   |    Real     |   Yes    |
        +-------+-----------+-------------+----------+
        """
        tcode = self.table_code // 1000
        sort_method = 1
        is_real = True
        is_random = False
        assert tcode in [0, 1, 2, 3, 4, 5], tcode
        if tcode in [2, 3, 5]:
            sort_method = 2
        if tcode in [1, 3]:
            is_real = False
        if tcode in [4, 5]:
            is_random = True
        return sort_method, is_real, is_random

    def is_sort1(self):
        sort_method, is_real, is_random = self._table_specs()
        return True if sort_method == 1 else False
        #if self.sort_bits[0] == 0:
        #    return True
        #return False

    def is_sort2(self):
        sort_method, is_real, is_random = self._table_specs()
        return True if sort_method == 2 else False
        #return not(self.is_sort1())

    def is_real(self):
        sort_method, is_real, is_random = self._table_specs()
        return is_real
        #return not(self.is_complex())

    def is_complex(self):
        sort_method, is_real, is_random = self._table_specs()
        return not(is_real)
        #if self.sort_bits[1] == 1:
        #    return True
        #return False

    def is_random(self):
        sort_method, is_real, is_random = self._table_specs()
        return is_random
        #if self.sort_bits[1] == 1:
        #    return True
        #return False

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

    def isStress(self):
        if self.stress_bits[1] == 0:
            return True
        return False

    def _create_table_object(self, result_name,  nnodes,
                             slot, slot_object, slot_vector, is_cid=False):
        assert isinstance(result_name, string_types), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        is_vectorized = self.is_vectorized
        if is_vectorized and slot_vector is None:
            is_vectorized = False

        if is_vectorized:
            if self.read_mode == 1:
                self.create_transient_object(slot, slot_vector, is_cid=is_cid)
                self.result_names.add(result_name)
                self.obj._nnodes += nnodes
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                self.obj = slot[self.code]
                #self.obj.update_data_code(self.data_code)
                self.obj.build()
        else:  # not vectorized
            self.result_names.add(result_name)
            if self.read_mode == 1:
                auto_return = True
                return auto_return
            # pass = 0/2
            self.create_transient_object(slot, slot_object)
        return auto_return
