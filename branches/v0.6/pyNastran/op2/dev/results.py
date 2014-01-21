from struct import unpack


class Results(object):
    def __init__(self):
        self.element_mapper = {
            # rods
            1: 'ROD',
            3 : 'CTUBE',
            10: 'CONROD',

            # beam
            2 : 'CBEAM',

            # springs
            11: 'CELAS1',
            12: 'CELAS2',
            13: 'CELAS3',
            
            # bar
            34 : 'CBAR',
            
            # solids
            39: 'TETRA',
            67: 'HEXA',
            68: 'PENTA',
            
            # triax
            53 : 'CTRIAX6',

            # centroidal plate
            33: 'QUAD4', # 1 node
            74: 'TRIA3', # 1 node

            # bilinear plate
            144: 'QUAD144', # 5 nodes

            # composite plates
            95 : 'CQUAD4',
            96 : 'CQUAD8',
            97 : 'CTRIA3',
            98 : 'CTRIA6',

            102: 'CBUSH',
        }
        pass

    def apply_data_code_value(self, name, value):
        pass
    def setNullNonlinearFactor(self):
        pass

    def read_title(self, data):
        assert len(data) == 584, len(data)
        Title, subtitle, label = unpack('128s128s128s', data[200:])  # titleSubtitleLabel

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

    def read_table(self, data, result_name, flag):
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                #obj = self.create_transient_object(self.spcForces, SPCForcesObject)
                self.read_real_table(data, result_name, flag)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()

    def read_oug_table(self, data, result_name, real_obj, complex_obj, node_elem):
        assert real_obj is None
        assert complex_obj is None

        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                self.read_real_table(data, result_name, node_elem)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                self.read_complex_table(data, result_name, node_elem)
            else:
                self.not_implemented_or_skip()
        else:
            msg = 'only num_wide=8 or 14 is allowed  num_wide=%s' % self.num_wide
            self.not_implemented_or_skip(msg)

    def read_real_table(self, data, result_name, flag):
        if self.debug:
            self.binary_debug.write('  read_real_table\n')
        assert flag in ['node', 'elem'], flag
        format1 = '2i6f' # 8

        ntotal = 32 # 8 * 4
        nnodes = len(data) // ntotal > 0
        assert len(data) % ntotal == 0

        n = 0
        for inode in xrange(nnodes):
            eData = data[n:n+32]
            out = unpack(format1, eData)
            (eid_device, gridType, tx, ty, tz, rx, ry, rz) = out
            eid = (eid_device - self.device_code) // 10
            #print "eType=%s" %(eType)

            dataIn = [eid, gridType, tx, ty, tz, rx, ry, rz]
            #self.obj.add(dt, dataIn)
            n += ntotal

    def read_complex_table(self, data, result_name, flag):
        if self.debug:
            self.binary_debug.write('  read_real_table\n')
        assert flag in ['node', 'elem'], flag

        format1 = '2i12f'
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 56  # 14 * 4
        nnodes = len(data) // ntotal
        for inode in xrange(nnodes):
            edata = data[n:n+ntotal]

            out = unpack(format1, edata)
            if self.debug:
                self.binary_debug.write('read_complex_table - %s\n' % str(out))
            (eid_device, gridType, txr, tyr, tzr, rxr, ryr, rzr,
                                   txi, tyi, tzi, rxi, ryi, rzi) = out

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

            eid = (eid_device - self.device_code) // 10

            dataIn = [eid, gridType, tx, ty, tz, rx, ry, rz]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            #self.obj.add(dt, dataIn)
            n += ntotal

    def not_implemented_or_skip(self, msg):
        if 1:
            raise NotImplementedError('table_name=%s table_code=%s' % (self.table_name, self.table_code))
        else:
            pass

