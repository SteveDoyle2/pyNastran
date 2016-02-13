from struct import unpack


from numpy import array

from pyNastran.op2.op2_common import OP2Common


class OGPWG(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogpwg_3(self, data, ndata):
        """
        Grid Point Weight Generator
        .. todo:: find the reference_point...
        """
        #self.show_data(data)
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', 'num_wide', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)
        self.reference_point = self.add_data_parameter(data, 'reference_point', 'i', 3, add_to_dict=False)

        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)

        self._read_title(data)
        self._write_debug_bits()

    def _read_ogpwg_4(self, data, ndata):
        """
        Grid Point Weight Generator
        """
        if self.read_mode == 1:
            return ndata
        MO = array(unpack('36f', data[:4*36]))
        MO = MO.reshape(6, 6)

        S = array(unpack('9f', data[4*36:4*(36+9)]))
        S = S.reshape(3, 3)

        mxyz = array(unpack('12f', data[4*(36+9):4*(36+9+12)]))
        mxyz = mxyz.reshape(3, 4)
        mass = mxyz[:, 0]
        cg = mxyz[:, 1:]

        IS = array(unpack('9f', data[4*(36+9+12):4*(36+9+12+9)]))
        IS = IS.reshape(3, 3)

        IQ = array(unpack('3f', data[4*(36+9+12+9):4*(36+9+12+9+3)]))

        Q = array(unpack('9f', data[4*(36+9+12+9+3):4*(36+9+12+9+3+9)]))
        Q = Q.reshape(3, 3)

        self.grid_point_weight.set_grid_point_weight(
            self.reference_point,
            MO, S, mass, cg, IS, IQ, Q)
        #del self.reference_point
        return ndata
