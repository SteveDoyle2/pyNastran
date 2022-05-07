from struct import unpack
from numpy import array
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.result_objects.grid_point_weight import GridPointWeight


class OGPWG(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogpwg_3(self, data: bytes, ndata: int):
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
        self.reference_point = self.add_data_parameter(data, 'reference_point', b'i', 3, add_to_dict=False)

        #13 OGPWG Grid Point Weight Generator
        #approach_code   = 1
        #tCode           = 13
        #isubcase        = 0
        #reference_point = 0
        #print('  approach_code   = %r' % self.approach_code)
        #print('  tCode           = %r' % self.tCode)
        #print('  isubcase        = %r' % self.isubcase)
        #print('  reference_point = %r' % self.reference_point)
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)

        self._read_title(data)
        #print('title = %r' % self.title)
        #print('subtitle = %r' % self.subtitle)
        #print('label = %r' % self.label)
        self._write_debug_bits()

    def _read_ogpwg_4(self, data: bytes, ndata: int):
        """
        Grid Point Weight Generator
        """
        if self.read_mode == 1:
            return ndata
        #print('  num_wide = %r' % self.num_wide)
        size = self.size
        fmt_mo = mapfmt(b'36f', size)
        fmt_s = mapfmt(b'9f', size)
        fmt_mxyz = mapfmt(b'12f', size)
        fmt_iq = mapfmt(b'3f', size)
        assert ndata == 78 * size, ndata
        MO = array(unpack(fmt_mo, data[:36*size]))
        MO = MO.reshape(6, 6)

        S = array(unpack(fmt_s, data[36*size:(36+9)*size]))
        S = S.reshape(3, 3)

        mxyz = array(unpack(fmt_mxyz, data[(36+9)*size:(36+9+12)*size]))
        mxyz = mxyz.reshape(3, 4)
        mass = mxyz[:, 0]
        cg = mxyz[:, 1:]

        IS = array(unpack(fmt_s, data[(36+9+12)*size:(36+9+12+9)*size]))
        IS = IS.reshape(3, 3)

        IQ = array(unpack(fmt_iq, data[(36+9+12+9)*size:(36+9+12+9+3)*size]))

        Q = array(unpack(fmt_s, data[(36+9+12+9+3)*size:(36+9+12+9+3+9)*size]))
        Q = Q.reshape(3, 3)

        #print(self.object_attributes())
        #print(self._count)
        #print(self.title)
        #print(self.subtitle)
        #print(self.label)
        #print(self.pval_step)
        #print(self.superelement_adaptivity_index)
        weight = GridPointWeight(
            self.reference_point,
            MO, S, mass, cg, IS, IQ, Q,
            approach_code=self.approach_code, table_code=self.table_code,
            title=self.title, subtitle=self.subtitle, label=self.label,
            superelement_adaptivity_index=self.superelement_adaptivity_index,
        )
        str(weight)
        self.grid_point_weight[self.superelement_adaptivity_index] = weight
        #del self.reference_point
        return ndata
