"""
defines the LAMA class to read:
 - RealEigenvalues
 - ComplexEigenvalues
 - BucklingEigenvalues

from the OP2
"""
from __future__ import print_function
from struct import Struct

#from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)

from pyNastran.op2.tables.table_reader import GenericTableReader


class LAMA(GenericTableReader):
    def __init__(self, op2_reader, op2):
        GenericTableReader.__init__(self, op2_reader, op2)

    #---------------------------------------------------------------------------
    def _read_complex_eigenvalue_3(self, data, ndata):
        """parses the Complex Eigenvalues Table 3 Data"""
        #raise NotImplementedError(self.table_name)
        op2_reader = self.op2_reader
        op2_reader.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)
        op2_reader.six = self.add_data_parameter(data, 'six', b'i', 10, False)  # seven
        self._read_title(data)

    def _read_buckling_eigenvalue_3(self, data, ndata):
        """parses the Buckling Eigenvalues Table 3 Data"""
        #print(self.show_data(data))
        #self._read_title_helper(data)
        op2_reader = self.op2_reader
        op2_reader.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #op2_reader.show_data(data)

        unused_three = op2_reader.parse_approach_code(data)

        op2_reader.seven = op2_reader.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #: residual vector augmentation flag
        op2_reader.residual_flag = op2_reader.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #: fluid modes flag
        op2_reader.fluid_flag = op2_reader.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        self._read_title(data)

    def _read_complex_eigenvalue_4(self, data, ndata):
        """parses the Complex Eigenvalues Table 4 Data"""
        if self.read_mode == 1:
            return ndata
        op2 = self.op2
        op2_reader = self.op2_reader

        ntotal = 4 * 6
        nmodes = ndata // ntotal
        n = 0
        #assert op2.isubcase != 0, op2.isubcase
        clama = ComplexEigenvalues(op2_reader.title, nmodes)
        op2.eigenvalues[op2_reader.title] = clama
        #op2.eigenvalues[op2.isubcase] = lama
        structi = Struct(self._endian + b'ii4f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigr, eigc, freq, damping) = out # CLAMA
            #print('imode=%s order=%s eigr=%s eigc=%s freq=%s damping=%s' %
                  #(imode, order, eigr, eigc, freq, damping))
            clama.add_f06_line(out, i)
            n += ntotal
        assert n == ndata, 'clama length error'
        return n

    def _read_buckling_eigenvalue_4(self, data, ndata):
        """parses the Buckling Eigenvalues Table 4 Data"""
        # BLAMA - Buckling eigenvalue summary table
        # CLAMA - Complex eigenvalue summary table
        # LAMA - Normal modes eigenvalue summary table
        op2_reader = self.op2_reader
        op2 = self.op2
        if op2_reader.read_mode == 1:
            return ndata

        ntotal = 4 * 7
        nmodes = ndata // ntotal
        n = 0
        #assert op2.isubcase != 0, op2.isubcase
        blama = BucklingEigenvalues(op2_reader.title, nmodes)
        op2.eigenvalues[op2_reader.title] = blama
        #op2.eigenvalues[self.isubcase] = lama
        structi = Struct(self._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(mode_num, extract_order, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(root_num, extract_order, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.add_f06_line(out, i)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data, ndata):
        """parses the Real Eigenvalues Table 3 Data"""
        op2_reader = self.op2_reader
        op2_reader.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)

        op2_reader.seven = self.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        ## residual vector augmentation flag
        op2_reader.residual_flag = self.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        ## fluid modes flag
        op2_reader.fluid_flag = self.add_data_parameter(data, 'fluid_flag', b'i', 12, False)
        op2_reader.title = None

        #print(op2.data_code)
        #op2.add_data_parameter(data,'format_code',  'i',9,False)   ## format code

        #: number of words per entry in record;
        #.. todo:: is this needed for this table ???
        #op2.add_data_parameter(data,'num_wide',     'i',10,False)

        #if self.analysis_code == 2: # sort2
            #op2.lsdvmn = op2.get_values(data,'i',5)

        #print("*isubcase=%s" % op2.isubcase)
        #print("analysis_code=%s table_code=%s thermal=%s" % (
            #op2.analysis_code, op2.table_code, op2.thermal))

        #op2.print_block(data)
        self._read_title(data)

    def _read_real_eigenvalue_4(self, data, ndata):
        """parses the Real Eigenvalues Table 4 Data"""
        op2 = self.op2
        if self.read_mode == 1:
            return ndata

        op2_reader = self.op2_reader
        #self.show_data(data)
        nmodes = ndata // 28
        n = 0
        ntotal = 28
        #assert op2.isubcase != 0, op2.isubcase
        lama = RealEigenvalues(op2.title, nmodes=nmodes)
        op2.eigenvalues[op2_reader.title] = lama
        structi = Struct(self._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+28]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = out
            lama.add_f06_line(out, i)
            n += ntotal
        return n
