from __future__ import print_function
from six import b
from six.moves import range
from struct import Struct

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)


class LAMA(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_complex_eigenvalue_3(self, data, ndata):
        """parses the Complex Eigenvalues Table 3 Data"""
        #raise NotImplementedError(self.table_name)
        self.words = [
            'aCode',       'tCode',    '???', 'isubcase',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        three = self.parse_approach_code(data)
        self.six = self.add_data_parameter(data, 'six', 'i', 10, False)  # seven
        self._read_title(data)

    def _read_buckling_eigenvalue_3(self, data, ndata):
        """parses the Buckling Eigenvalues Table 3 Data"""
        #print(self.show_data(data))
        #self._read_title_helper(data)

        self.words = [
            'aCode',       'tCode',    '???', 'isubcase',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        three = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.resFlag = self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.fldFlag = self.add_data_parameter(data, 'fldFlag', 'i', 12, False)

        self._read_title(data)

    def _read_complex_eigenvalue_4(self, data, ndata):
        """parses the Complex Eigenvalues Table 4 Data"""
        if self.read_mode == 1:
            return ndata

        ntotal = 4 * 6
        nmodes = ndata // ntotal
        n = 0
        #assert self.isubcase != 0, self.isubcase
        clama = ComplexEigenvalues(11)
        self.eigenvalues[self.title] = clama
        #self.eigenvalues[self.isubcase] = lama
        s = Struct(b(self._endian + 'ii4f'))
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            (imode, order, eigr, eigc, freq, damping) = out # CLAMA
            #print('imode=%s order=%s eigr=%s eigc=%s freq=%s damping=%s' %
                  #(imode, order, eigr, eigc, freq, damping))
            clama.add_f06_line(out)
            n += ntotal
        assert n == ndata, 'clama length error'
        return n

    def _read_buckling_eigenvalue_4(self, data, ndata):
        """parses the Buckling Eigenvalues Table 4 Data"""
        # BLAMA - Buckling eigenvalue summary table
        # CLAMA - Complex eigenvalue summary table
        # LAMA - Normal modes eigenvalue summary table
        if self.read_mode == 1:
            return ndata

        ntotal = 4 * 7
        nmodes = ndata // ntotal
        n = 0
        #assert self.isubcase != 0, self.isubcase
        blama = BucklingEigenvalues(11)
        self.eigenvalues[self.title] = blama
        #self.eigenvalues[self.isubcase] = lama
        s = Struct(b(self._endian + 'ii5f'))
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            (imode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(rootNum, extractOrder, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.add_f06_line(out)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data, ndata):
        """parses the Real Eigenvalues Table 3 Data"""
        self.words = [
            'aCode',       'tCode',    '???', 'isubcase',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???',         '???',      '???',          '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        three = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.resFlag = self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.fldFlag = self.add_data_parameter(data, 'fldFlag', 'i', 12, False)
        self.title = None

        #print(self.data_code)
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code
        #self.add_data_parameter(data,'num_wide',     'i',10,False)  ## number of words per entry in record; .. note:: is this needed for this table ???

        #if self.analysis_code == 2: # sort2
            #self.lsdvmn = self.get_values(data,'i',5)

        #print("*isubcase=%s" % self.isubcase)
        #print("analysis_code=%s table_code=%s thermal=%s" % (self.analysis_code, self.table_code, self.thermal))

        #self.print_block(data)
        self._read_title(data)

    def _read_real_eigenvalue_4(self, data, ndata):
        """parses the Real Eigenvalues Table 4 Data"""
        if self.read_mode == 1:
            return ndata
        #self.show_data(data)
        nModes = ndata // 28
        n = 0
        ntotal = 28
        #assert self.isubcase != 0, self.isubcase
        lama = RealEigenvalues(self.title)
        self.eigenvalues[self.title] = lama
        s = Struct('ii5f')
        for i in range(nModes):
            edata = data[n:n+28]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigen, omega, freq, mass, stiff) = out
            (imode, extract_order, eigenvalue, radian, cycle, genM, genK) = out
            lama.add_f06_line(out)
            n += ntotal
        return n
