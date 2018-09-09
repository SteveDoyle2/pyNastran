from numpy import array

from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray
from pyNastran.f06.f06_formatting import write_floats_13e


class ComplexEigenvectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        words = []
        #freq = self.eigrs[i]
        #freq = 0.0
        #words.append('%16s = %12E\n' % ('EIGENVALUE', freq))
        #words.append('%%16s = %%12E\n')
        has_cycle = False
        if hasattr(self, 'cycles'):
            has_cycle = True
        #if has_cycle:
        #words.append('%s = %s          C O M P L E X   E I G E N V E C T O R   N O. %s\n' % ('%16s', '%12E', '%10i'))
            #msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
        #else:
        words.append('      COMPLEX EIGENVALUE = %s, %s\n' % ('%12E', '%12E'))
        words.append('                                       C O M P L E X   E I G E N V E C T O R   NO. %s\n' % '%10i')
        #msg.append('                                         C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % (imode))

        #msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)


class RealEigenvectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt, f06_flag=False):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

        if f06_flag:
            self.itime = -1
            self._times = []
            self.data = []
            self.node_gridtype = []

            self.data_code = data_code
            self.apply_data_code()
            self.set_data_members()

    def get_phi(self):
        """
        gets the eigenvector matrix

        Returns
        -------
        phi : (ndof, nmodes)
            the eigenvector matrix

        TODO: doesn't consider SPOINTs/EPOINTs
        """
        nmodes, nnodes = self.data.shape[:2]
        ndof = nnodes * 6
        phi_transpose = self.data.reshape(nmodes, ndof)
        return phi_transpose.T

    @classmethod
    def phi_to_data(self, phi):
        """(ndof, nmodes) -> (nmodes, nnodes, 6)"""
        ndof, nmodes = phi.shape
        nnodes = ndof // 6
        assert ndof % 6 == 0
        phi2 = phi.T  # nmodes, ndof
        data = phi2.reshape(nmodes, nnodes, 6)
        return data

    def set_phi(self, phi):
        """(ndof, nmodes) -> (nmodes, nnodes, 6)"""
        self.data = self.phi_to_data(phi)

    def build_f06_vectorization(self):
        self.data = array(self.data, dtype='float32')
        self.node_gridtype = array(self.node_gridtype)
        self._times = array(self._times)

    def add_new_f06_mode(self, imode):
        self._times.append(imode)
        self.itime += 1
        self.ntimes += 1
        self.data.append([])
        self.node_gridtype = []

    def update_f06_mode(self, data_code, imode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode >= 0.
        self.data_code = data_code
        self.apply_data_code()
        #self.caseVal = imode
        self.add_new_f06_mode(imode)
        self.set_data_members()

    def read_f06_data(self, data_code, lines):
        """
        it is assumed that all data coming in is correct, so...

           [node_id, grid_type, t1, t2, t3, r1, r2, r3]
           [100, 'G', 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] #  valid

           [101, 'S', 1.0, 2.0] #  invalid
           [101, 'S', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
           [102, 'S', 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] #  valid
        """
        imode = data_code['mode']
        if imode not in self._times:
            self.update_f06_mode(data_code, imode)

        itime = self.itime
        datai = self.data[itime]

        for line in lines:
            if len(line) != 8:
                msg = 'invalid length; even spoints must be in \n'
                msg += '[nid, type, t1, t2, t3, r1, r2, t3] format\nline=%s\n' % line
                msg += 'expected length=8; length=%s' % len(line)
                raise RuntimeError(msg)
            (node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.node_gridtype.append([node_id, self.cast_grid_type(grid_type)])
            datai.append([t1, t2, t3, r1, r2, r3])
        assert self.eigrs[-1] == data_code['eigr'], 'eigrs=%s\ndata_code[eigrs]=%s' %(self.eigrs, data_code['eigr'])

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #if self.nonlinear_factor not in (None, np.nan):
            #return self._write_f06_transient(header, page_stamp, page_num, f06_file,
                                             #is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        # modes get added
        words = '                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' \
                '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n'

        #if not len(header) >= 3:
            #header.append('')
        for itime in range(self.ntimes):
            node = self.node_gridtype[:, 0]
            gridtype = self.node_gridtype[:, 1]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            dt = self._times[itime]
            #if isinstance(dt, float):
                #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #else:
                #header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + [words % dt]))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write(
                    '%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        node_id, sgridtype, dx, dy, dz, rx, ry, rz))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
