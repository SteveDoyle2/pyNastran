from six.moves import range
from struct import Struct

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)


class LAMA(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_complex_eigenvalue_3(self, data):
        raise NotImplementedError(self.table_name)

    def _read_buckling_eigenvalue_3(self, data):
        #print self.show_data(data)
        #self._read_title_helper(data)

        self.words = ['aCode',       'tCode',    '???', 'isubcase',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', '???', '???', '???']
        #self.show_data(data)

        (three) = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.resFlag = self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.fldFlag = self.add_data_parameter(data, 'fldFlag', 'i', 12, False)

        #print(self.data_code)
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code
        #self.add_data_parameter(data,'num_wide',     'i',10,False)  ## number of words per entry in record; .. note:: is this needed for this table ???

        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print("*isubcase=%s"%(self.isubcase)()
        #print("analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code self.table_code, self.thermal))

        #self.print_block(data)
        self._read_title(data)

        #assert self.isubcase != 0, self.isubcase
        #print(self.Title)
        #raise NotImplementedError(self.table_name)

    def _read_complex_eigenvalue_4(self, data):
        if self.read_mode == 1:
            return len(data)
        raise NotImplementedError(self.table_name)

    def _read_buckling_eigenvalue_4(self, data):
        # BLAMA - Buckling eigenvalue summary table
        # CLAMA - Complex eigenvalue summary table
        # LAMA - Normal modes eigenvalue summary table
        if self.read_mode == 1:
            return len(data)

        msg = '_read_buckling_eigenvalue_4'
        #return self._not_implemented_or_skip(data, msg)  # TODO: implement buckling eigenvalues

        ntotal = 4 * 7
        nModes = len(data) // ntotal
        n = 0
        #assert self.isubcase != 0, self.isubcase
        blama = BucklingEigenvalues(11)
        self.eigenvalues[self.Title] = blama
        #self.eigenvalues[self.isubcase] = lama
        s = Struct(b'ii5f')
        for i in range(nModes):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if self.debug4():
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            (iMode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(rootNum, extractOrder, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.addF06Line(out)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data):
        self.words = ['aCode',       'tCode',    '???', 'isubcase',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', '???', '???', '???']
        #self.show_data(data)

        (three) = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.resFlag = self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.fldFlag = self.add_data_parameter(data, 'fldFlag', 'i', 12, False)
        self.Title = None

        #print(self.data_code)
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code
        #self.add_data_parameter(data,'num_wide',     'i',10,False)  ## number of words per entry in record; .. note:: is this needed for this table ???

        #if self.analysis_code == 2: # sort2
            #self.lsdvmn = self.get_values(data,'i',5)

        #print("*isubcase=%s" % self.isubcase)
        #print("analysis_code=%s table_code=%s thermal=%s" % (self.analysis_code, self.table_code, self.thermal))

        #self.print_block(data)
        self._read_title(data)

    def _read_real_eigenvalue_4(self, data):
        if self.read_mode == 1:
            return len(data)
        #self.show_data(data)
        nModes = len(data) // 28
        n = 0
        ntotal = 28
        #assert self.isubcase != 0, self.isubcase
        lama = RealEigenvalues(self.Title)
        self.eigenvalues[self.Title] = lama
        s = Struct('ii5f')
        for i in range(nModes):
            edata = data[n:n+28]
            out = s.unpack(edata)
            if self.debug4():
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(iMode, order, eigen, omega, freq, mass, stiff) = out
            (modeNum, extractOrder, eigenvalue, radian, cycle, genM, genK) = out
            lama.addF06Line(out)
            n += ntotal
        return n
