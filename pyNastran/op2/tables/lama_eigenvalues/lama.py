"""
defines the LAMA class to read:
 - RealEigenvalues
 - ComplexEigenvalues
 - BucklingEigenvalues

from the OP2

"""
from struct import Struct

from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)


class LAMA(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_complex_eigenvalue_3(self, data, ndata):
        """parses the Complex Eigenvalues Table 3 Data"""
        #raise NotImplementedError(self.table_name)
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)
        self.six = self.add_data_parameter(data, 'six', b'i', 10, False)  # seven
        self._read_title(data)

    def _read_buckling_eigenvalue_3(self, data, ndata):
        """parses the Buckling Eigenvalues Table 3 Data"""
        #print(self.show_data(data))
        #self._read_title_helper(data)

        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #: residual vector augmentation flag
        self.residual_flag = self.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #: fluid modes flag
        self.fluid_flag = self.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        self._read_title(data)

    def _read_complex_eigenvalue_4(self, data, ndata):
        """parses the Complex Eigenvalues Table 4 Data"""
        if self.read_mode == 1:
            return ndata

        ntotal = 24 # 4 * 6
        nmodes = ndata // ntotal
        n = 0
        #assert self.isubcase != 0, self.isubcase
        clama = ComplexEigenvalues(self.title, self.table_name, nmodes)
        #assert self.title not in self.eigenvalues, f'table={self.table_name_str} title={self.title} optimization_count={self._count}'
        self.eigenvalues[self.title] = clama
        #self.eigenvalues[self.isubcase] = clama
        structi = Struct(self._endian + b'ii4f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigr, eigc, freq, damping) = out # CLAMA
            #print('imode=%s order=%s eigr=%s eigc=%s freq=%s damping=%s' %
                  #(imode, order, eigr, eigc, freq, damping))
            clama.add_op2_line(out, i)
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

        ntotal = 28 # 4 * 7
        nmodes = ndata // ntotal
        n = 0
        #assert self.isubcase != 0, self.isubcase
        blama = BucklingEigenvalues(self.title, self.table_name, nmodes)
        #assert self.title not in self.eigenvalues, f'table={self.table_name_str} title={self.title} optimization_count={self._count}'
        self.eigenvalues[self.title] = blama
        #self.eigenvalues[self.isubcase] = lama
        structi = Struct(self._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(mode_num, extract_order, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(root_num, extract_order, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.add_op2_line(out, i)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data, ndata):
        """parses the Real Eigenvalues Table 3 Data"""
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)

        self.seven = self.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.residual_flag = self.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        ## fluid modes flag
        self.fluid_flag = self.add_data_parameter(data, 'fluid_flag', b'i', 12, False)
        self.title = None

        #print(self.data_code)
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code

        #: number of words per entry in record;
        #.. todo:: is this needed for this table ???
        #self.add_data_parameter(data,'num_wide',     'i',10,False)

        #if self.analysis_code == 2: # sort2
            #self.lsdvmn = self.get_values(data,'i',5)

        #print("*isubcase=%s" % self.isubcase)
        #print("analysis_code=%s table_code=%s thermal=%s" % (
            #self.analysis_code, self.table_code, self.thermal))

        #self.print_block(data)
        self._read_title(data)

    def _read_real_eigenvalue_4(self, data, ndata):
        """parses the Real Eigenvalues Table 4 Data"""
        if self.read_mode == 1:
            return ndata
        nmodes = ndata // 28
        n = 0
        ntotal = 28
        #assert self.isubcase != 0, self.isubcase
        lama = RealEigenvalues(self.title, self.table_name, nmodes=nmodes)

        if self.table_name in [b'LAMA', b'LAMAS']:
            result_name = 'eigenvalues'
        elif self.table_name == b'LAMAF':
            result_name = 'eigenvalues_fluid'
        else:  # pragma: no cover
            raise NotImplementedError(self.table_name)
        slot = getattr(self, result_name)
        #assert self.title not in slot, f'{result_name}: table={self.table_name_str} title={self.title!r} optimization_count={self._count}'
        slot[self.title] = lama

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
