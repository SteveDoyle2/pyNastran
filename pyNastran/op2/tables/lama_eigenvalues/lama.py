"""
defines the LAMA class to read:
 - RealEigenvalues
 - ComplexEigenvalues
 - BucklingEigenvalues

from the OP2

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from struct import Struct

from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


class LAMA:
    def __init__(self, op2: OP2):
        self.op2 = op2

    def _read_complex_eigenvalue_3(self, data: bytes, ndata: int):
        """parses the Complex Eigenvalues Table 3 Data"""
        #raise NotImplementedError(self.table_name)
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = op2.parse_approach_code(data)
        op2.six = op2.add_data_parameter(data, 'six', b'i', 10, False)  # seven
        op2._read_title(data)

    def _read_buckling_eigenvalue_3(self, data: bytes, ndata: int):
        """parses the Buckling Eigenvalues Table 3 Data"""
        op2 = self.op2

        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']

        unused_three = op2.parse_approach_code(data)

        op2.seven = op2.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #: residual vector augmentation flag
        op2.residual_flag = op2.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #: fluid modes flag
        op2.fluid_flag = op2.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        op2._read_title(data)

    def _read_complex_eigenvalue_4(self, data: bytes, ndata: int):
        """parses the Complex Eigenvalues Table 4 Data"""
        op2: OP2 = self.op2
        if op2.read_mode == 1:
            return ndata

        n = 0
        ntotal = 24 * op2.factor # 4 * 6
        nmodes = ndata // ntotal
        if op2.table_name in [b'CLAMA']:
            result_name = 'eigenvalues'
        #elif op2.table_name == b'LAMAF':
            #result_name = 'eigenvalues_fluid'
        else:  # pragma: no cover
            raise NotImplementedError(op2.table_name)

        #assert self.isubcase != 0, self.isubcase
        clama = ComplexEigenvalues(op2.title, op2.table_name, nmodes)
        #assert self.title not in self.eigenvalues, f'table={self.table_name_str} title={self.title} optimization_count={self._count}'
        op2.eigenvalues[op2.title] = clama
        #self.eigenvalues[self.isubcase] = clama
        fmt = mapfmt(op2._endian + b'ii4f', op2.size)
        structi = Struct(fmt)
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigr, eigc, freq, damping) = out # CLAMA
            #print('imode=%s order=%s eigr=%s eigc=%s freq=%s damping=%s' %
                  #(imode, order, eigr, eigc, freq, damping))
            clama.add_op2_line(out, i)
            n += ntotal
        assert n == ndata, 'clama length error'
        return n

    def _read_buckling_eigenvalue_4(self, data: bytes, ndata: int):
        """parses the Buckling Eigenvalues Table 4 Data"""
        op2: OP2 = self.op2
        # BLAMA - Buckling eigenvalue summary table
        # CLAMA - Complex eigenvalue summary table
        # LAMA - Normal modes eigenvalue summary table
        if op2.read_mode == 1:
            return ndata

        n = 0
        ntotal = 28 * op2.factor # 4 * 7
        nmodes = ndata // ntotal
        #assert self.isubcase != 0, self.isubcase
        blama = BucklingEigenvalues(op2.title, op2.table_name, nmodes)

        if op2.table_name in [b'BLAMA']:
            result_name = 'eigenvalues'
        #elif op2.table_name == b'LAMAF':
            #result_name = 'eigenvalues_fluid'
        else:  # pragma: no cover
            raise NotImplementedError(op2.table_name)

        #assert self.title not in self.eigenvalues, f'table={self.table_name_str} title={self.title} optimization_count={self._count}'
        op2.eigenvalues[op2.title] = blama
        #self.eigenvalues[self.isubcase] = lama
        structi = Struct(op2._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(mode_num, extract_order, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(root_num, extract_order, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.add_op2_line(out, i)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data: bytes, ndata: int):
        """parses the Real Eigenvalues Table 3 Data"""
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = op2.parse_approach_code(data)

        op2.seven = op2.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        ## residual vector augmentation flag
        op2.residual_flag = op2.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        ## fluid modes flag
        op2.fluid_flag = op2.add_data_parameter(data, 'fluid_flag', b'i', 12, False)
        op2.title = None

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
        op2._read_title(data)

    def _read_real_eigenvalue_4(self, data: bytes, ndata: int):
        """parses the Real Eigenvalues Table 4 Data"""
        op2: OP2 = self.op2
        if op2.read_mode == 1:
            return ndata

        n = 0
        ntotal = 28 * op2.factor
        nmodes = ndata // ntotal
        #assert self.isubcase != 0, self.isubcase
        lama = RealEigenvalues(op2.title, op2.table_name, nmodes=nmodes)

        if op2.table_name in [b'LAMA', b'LAMAS']:
            result_name = 'eigenvalues'
        elif op2.table_name == b'LAMAF':
            result_name = 'eigenvalues_fluid'
        else:  # pragma: no cover
            raise NotImplementedError(op2.table_name)
        slot = getattr(op2, result_name)
        #assert self.title not in slot, f'{result_name}: table={self.table_name_str} title={self.title!r} optimization_count={self._count}'
        slot[op2.title] = lama

        structi = Struct(op2._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = out
            lama.add_f06_line(out, i)
            n += ntotal
        return n
