from struct import unpack


from numpy import array

from pyNastran.op2.dev.op2_common import OP2Common


class OGPWG(object):
    def __init__(self):
        #OP2Common.__init__(self)
        pass

    def readTable_OGPWG(self, table_name):
        self.readRecordTable_OGPWG(table_name)

    def readRecordTable_OGPWG(self, expectedTableName):
        """
        .. note:: assumes self.iTableMap has already been set
        """
        table_name = self.read_table_name(rewind=False)  # GEOM1
        self._table_init(table_name)

        #print "*table_name = |%r|" %(table_name)

        self.read_markers([-1, 7])
        fields = self.read_int_block()
        #print "fields = ",fields

        self.read_markers([-2, 1, 0])  # 2
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words,buffer_words*4
        word = self.read_string_block()

        iTable = -3
        while 1:  # TODO could this cause an infinite loop...i dont this so...
            (table_name, isNextTable, isNextSubTable,
                isFileDone) = self.readGeomSubTable_OGPWG(iTable)

            if self.checkForNextTable() or isFileDone:
                #sys.exit('end of geom1')
                return
            iTable -= 1
        sys.exit('end of %s-this should never happen...' % expectedTableName)

    def readGeomSubTable_OGPWG(self, iTable):
        i = 0
        isNextTable = False
        isNextSubTable = False
        self.read_markers([iTable, 1, 0])

        table_name = self.read_table_name(rewind=True, stopOnFailure=False)
        if table_name:
            return table_name, isNextTable, isNextSubTable, False

        data = b''
        isTableActive = False
        while isNextSubTable == False and isNextTable == False:
            marker = self.get_marker()
            if marker < 0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data += self.read_block()

            #print "iTable=%s lenGeomData=%s" %(iTable,len(data))
            if iTable == -3:
                self._read_ogpwg_3(data)
            elif iTable == -4:
                self._read_ogpwg_4(data)
            else:
                raise RuntimeError('bad iTable')

            isNextTable = self.checkForNextTable()
            isNextSubTable = self.checkForNextSubTable(iTable - 1)
            isFileDone = self.checkFileDone(iTable - 1)
            if isFileDone:
                isNextTable = True

            i += 1
            isTableActive = True

        return (table_name, isNextTable, isNextSubTable, isFileDone)

    def _read_ogpwg_3(self, data):
        """
        Grid Point Weight Generator
        ..todo:: find the reference_point...
        """
        #self.show_data(data)
        self.words = [
                 'aCode',       'tCode',    '???',     'isubcase',
                 '???',         '???',      '???',          '???',
                 '???',         'num_wide', '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', 'Title', 'subtitle', 'label']
        #self.print_block(data)

        self.parse_approach_code(data)
        self.add_data_parameter(data, 'reference_point', 'i', 3)

        #if self.debug3():
            #self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            #self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            #self.binary_debug.write('  isubcase = %r\n' % self.isubcase)

        #self.read_title(data)
        #self.write_debug_bits()
        self.data = b''

    def _read_ogpwg_4(self, data):
        """
        Grid Point Weight Generator
        """
        MO = array(unpack('36f', data[:4*36]))
        MO = MO.reshape(6,6)

        S = array(unpack('9f', data[4*36:4*(36+9)]))
        S = S.reshape(3,3)

        mxyz = array(unpack('12f', data[4*(36+9):4*(36+9+12)]))
        mxyz = mxyz.reshape(3,4)
        mass = mxyz[:, 0]
        cg = mxyz[:, 1:]

        IS = array(unpack('9f', data[4*(36+9+12):4*(36+9+12+9)]))
        IS = IS.reshape(3,3)

        IQ = array(unpack('3f', data[4*(36+9+12+9):4*(36+9+12+9+3)]))

        Q = array(unpack('9f', data[4*(36+9+12+9+3):4*(36+9+12+9+3+9)]))
        Q = Q.reshape(3,3)

        self.grid_point_weight.set_grid_point_weight(self.reference_point,
            MO, S, mass, cg, IS, IQ, Q)
        del self.reference_point