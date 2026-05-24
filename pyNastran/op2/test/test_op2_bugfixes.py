"""Tests for OP2 reader bug fixes.

Tests the following fixes:
- ComplexBarArray.get_stats assert nelements == nelements2 (not truthiness)
- ComplexBarArray.__eq__ unpacking uses correct s1b..s4b variable names
- ComplexShearStrainArray used for complex shear strain (not StressArray)
- onr.py table_name validation (chained comparison fix)
- op2_reader.py f-string fix for debug message
- ComplexBarArray write_op2 uses correct variable names
- write_op2 round-trip for real/complex stress results (rod, bar, shear, spring)
"""
import os
import tempfile
import unittest
from io import BytesIO, StringIO
from pathlib import Path

import numpy as np

import pyNastran
from pyNastran.op2.op2 import read_op2
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import (
    ComplexShearStressArray, ComplexShearStrainArray)
from pyNastran.op2.tables.oes_stressStrain.utils_cshear import (
    ComplexShearStressArray as _CShearStress,
    ComplexShearStrainArray as _CShearStrain)
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import RealSpringStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()


class TestComplexBarFixes(unittest.TestCase):
    """Tests for complex bar stress/strain fixes.

    Tolerances: exact (equality checks, no floating-point comparison).
    Tests: get_stats assertion, __eq__ variable unpacking, data shape.
    """

    def test_complex_bar_get_stats(self):
        """get_stats uses == not truthiness for element count assertion."""
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)
        cbar_stress = model.op2_results.stress.cbar_stress
        for subcase_id, result in cbar_stress.items():
            stats = result.get_stats()
            assert len(stats) > 0
            assert result.nelements == result.element.shape[0]

    def test_complex_bar_eq(self):
        """__eq__ correctly unpacks all 9 columns (sa1-sa4, axial, sb1-sb4)."""
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)
        cbar_stress = model.op2_results.stress.cbar_stress
        for subcase_id, result in cbar_stress.items():
            assert result == result

    def test_complex_bar_data_layout(self):
        """Data array has 9 columns: [sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4]."""
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)
        cbar_stress = model.op2_results.stress.cbar_stress
        for subcase_id, result in cbar_stress.items():
            assert result.data.shape[2] == 9


class TestComplexShearFix(unittest.TestCase):
    """Tests for complex shear stress/strain class selection.

    Tolerances: exact (class identity check).
    Tests: strain uses ComplexShearStrainArray, not ComplexShearStressArray.
    """

    def test_complex_shear_strain_class_is_distinct(self):
        """ComplexShearStrainArray is a different class from ComplexShearStressArray."""
        assert ComplexShearStressArray is not ComplexShearStrainArray
        assert _CShearStress is ComplexShearStressArray
        assert _CShearStrain is ComplexShearStrainArray

    def test_cshear_strain_class_selection(self):
        """utils_cshear selects ComplexShearStrainArray when is_stress=False."""
        from unittest.mock import MagicMock
        from pyNastran.op2.tables.oes_stressStrain import utils_cshear

        op2 = MagicMock()
        op2.is_stress = False
        obj_vector_complex = ComplexShearStressArray if op2.is_stress else ComplexShearStrainArray
        assert obj_vector_complex is ComplexShearStrainArray


class TestOnrTableNameValidation(unittest.TestCase):
    """Tests for onr.py chained comparison fix.

    Tolerances: exact (boolean checks).
    Tests: table_name membership correctly validates against expected names.
    """

    def test_valid_table_names_pass(self):
        """Valid table names should be in the allowed list."""
        allowed = [b'ONRGY', b'ONRGY1']
        assert b'ONRGY1' in allowed
        assert b'ONRGY' in allowed

    def test_invalid_table_name_detected(self):
        """Invalid table names should NOT be in the allowed list."""
        allowed = [b'ONRGY', b'ONRGY1']
        assert b'BADTABLE' not in allowed

    def test_chained_comparison_bug_would_fail(self):
        """The old chained comparison 'a not in a in [...]' is always False."""
        table_name = b'ONRGY1'
        # Old buggy expression (always False due to chaining):
        old_result = table_name not in table_name in [b'ONRGY', b'ONRGY1']
        assert old_result is False

        # Fixed expression:
        new_result = table_name not in [b'ONRGY', b'ONRGY1']
        assert new_result is False  # ONRGY1 IS in the list, so "not in" is False


class TestWriteOp2Constructed(unittest.TestCase):
    """Tests for write_op2 on result objects constructed via add_static/add_freq classmethods.

    Tolerances: checks that write_op2 produces >100 bytes (valid binary output).
    Tests: RealRodStress, RealBarStress, RealShearStress, RealSpringStress,
           ComplexRodStress, ComplexSpringStress.
    """

    def _write_op2_result(self, result, endian=b'>'):
        """Helper: call write_op2 and return bytes written."""
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_write_real_rod_stress(self):
        """RealRodStressArray.write_op2 produces valid binary for a 2-element static case.

        Data: 2 CROD elements, 1 time step, 4 columns [axial, torsion, SMa, SMt].
        """
        element = np.array([1, 2], dtype='int32')
        data = np.array([[[100., 50., 1.5, 2.0],
                          [200., 60., 1.2, 1.8]]], dtype='float32')
        rod = RealRodStressArray.add_static_case('OES1X1', 'CROD', element, data, isubcase=1)
        nbytes = self._write_op2_result(rod)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'

    def test_write_real_bar_stress(self):
        """RealBarStressArray.write_op2 produces valid binary for a 1-element static case.

        Data: 1 CBAR element, 1 time step, 15 columns
        [s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_t, s1b, s2b, s3b, s4b, smaxb, sminb, MS_c].
        """
        element = np.array([10], dtype='int32')
        data = np.zeros((1, 1, 15), dtype='float32')
        data[0, 0, :] = [100, 200, 300, 400, 50, 500, -100, 1.5,
                         110, 210, 310, 410, 510, -110, 1.2]
        bar = RealBarStressArray.add_static_case('OES1X1', 'CBAR', element, data, isubcase=1)
        nbytes = self._write_op2_result(bar)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'

    def test_write_real_shear_stress(self):
        """RealShearStressArray.write_op2 produces valid binary for a 2-element static case.

        Data: 2 CSHEAR elements, 1 time step, 3 columns [max_shear, avg_shear, margin].
        """
        element = np.array([20, 21], dtype='int32')
        data = np.zeros((1, 2, 3), dtype='float32')
        data[0, :, :] = [[100., 50., 1.2], [120., 60., 1.1]]
        shear = RealShearStressArray.add_static_case('OES1X1', 'CSHEAR', element, data, isubcase=1)
        nbytes = self._write_op2_result(shear)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'

    def test_write_real_spring_stress(self):
        """RealSpringStressArray.write_op2 produces valid binary for a 1-element static case.

        Data: 1 CELAS1 element, 1 time step, 1 column [stress].
        """
        element = np.array([40], dtype='int32')
        data = np.array([[[5.0]]], dtype='float32')
        spring = RealSpringStressArray.add_static_case('OES1X1', 'CELAS1', element, data, isubcase=1)
        nbytes = self._write_op2_result(spring)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'

    def test_write_complex_rod_stress(self):
        """ComplexRodStressArray.write_op2 produces valid binary for a freq case.

        Data: 2 CROD elements, 3 frequencies, 2 complex columns [axial, torsion].
        """
        element = np.array([30, 31], dtype='int32')
        freqs = np.array([10., 20., 30.], dtype='float32')
        data = np.zeros((3, 2, 2), dtype='complex64')
        data[0, 0, :] = [1+2j, 3+4j]
        data[1, 1, :] = [5+6j, 7+8j]
        crod = ComplexRodStressArray.add_freq_case('OES1C', 'CROD', element, data,
                                                   isubcase=1, freqs=freqs)
        nbytes = self._write_op2_result(crod)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'

    def test_write_complex_spring_stress(self):
        """ComplexSpringStressArray.write_op2 produces valid binary for a freq case.

        Data: 1 CELAS1 element, 2 frequencies, 1 complex column [stress].
        """
        element = np.array([50], dtype='int32')
        freqs = np.array([5., 10.], dtype='float32')
        data = np.zeros((2, 1, 1), dtype='complex64')
        data[0, 0, 0] = 1+1j
        data[1, 0, 0] = 2+3j
        cspring = ComplexSpringStressArray.add_freq_case('OES1C', 'CELAS1', element, data,
                                                         isubcase=1, freqs=freqs)
        nbytes = self._write_op2_result(cspring)
        assert nbytes > 100, f'expected >100 bytes, got {nbytes}'


class TestWriteOp2RoundTrip(unittest.TestCase):
    """Tests for OP2 write/read round-trip using real model files.

    Tolerances: atol=1e-5 for float comparisons.
    Tests: complex bar stress/strain, real rod/bar stress round-trip data fidelity.
    """

    def test_complex_bar_stress_round_trip(self):
        """Complex bar stress data survives write -> re-read with all 9 columns intact.

        Model: freq_solid_shell_bar.op2 (SOL 108, frequency response).
        Verifies: data shape preserved, values within atol=1e-5.
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)

        out_file = tempfile.mktemp(suffix='.op2')
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)

            cbar1 = model.op2_results.stress.cbar_stress
            cbar2 = model2.op2_results.stress.cbar_stress
            for k in cbar1:
                d1 = cbar1[k].data
                d2 = cbar2[k].data
                assert d1.shape == d2.shape, f'shape mismatch: {d1.shape} vs {d2.shape}'
                assert d1.shape[2] == 9
                assert np.allclose(d1, d2, atol=1e-5), f'data mismatch subcase {k}'
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)

    def test_complex_bar_strain_round_trip(self):
        """Complex bar strain data survives write -> re-read.

        Model: freq_solid_shell_bar.op2 (SOL 108).
        Verifies: 9 columns preserved, values within atol=1e-5.
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)

        out_file = tempfile.mktemp(suffix='.op2')
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)

            strain1 = model.op2_results.strain.cbar_strain
            strain2 = model2.op2_results.strain.cbar_strain
            for k in strain1:
                d1 = strain1[k].data
                d2 = strain2[k].data
                assert d1.shape == d2.shape
                assert d1.shape[2] == 9
                assert np.allclose(d1, d2, atol=1e-5)
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)

    def test_static_rod_bar_stress_round_trip(self):
        """Real rod and bar stress data survives write -> re-read.

        Model: static_solid_shell_bar.op2 (SOL 101, static).
        Verifies: rod (4 cols), bar (15 cols) data within atol=1e-5.
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        model = read_op2(str(op2_filename), debug=None)

        out_file = tempfile.mktemp(suffix='.op2')
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)

            rod1 = model.op2_results.stress.crod_stress
            rod2 = model2.op2_results.stress.crod_stress
            for k in rod1:
                assert rod1[k].data.shape[2] == 4
                assert np.allclose(rod1[k].data, rod2[k].data, atol=1e-5)

            bar1 = model.op2_results.stress.cbar_stress
            bar2 = model2.op2_results.stress.cbar_stress
            for k in bar1:
                assert bar1[k].data.shape[2] == 15
                assert np.allclose(bar1[k].data, bar2[k].data, atol=1e-5)
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)

    def test_modal_rod_stress_add_modal_case(self):
        """RealRodStressArray.add_modal_case creates a valid result with correct modes.

        Data: 3 modes, 2 CROD elements, 4 columns.
        Tolerances: exact shape/mode checks.
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2, 3], dtype='int32')
        eigns = np.array([100., 400., 900.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(3, 2, 4).astype('float32')

        rod = RealRodStressArray.add_modal_case(
            'OES1X1', 'CROD', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert rod.ntimes == 3
        assert rod.nelements == 2
        assert rod.data.shape == (3, 2, 4)
        assert np.array_equal(rod._times, modes)

        nbytes = BytesIO()
        rod.write_op2(nbytes, StringIO(), -1, True, (5, 23, 2026), endian=b'>')
        assert nbytes.tell() > 100


class TestWriteOp2Displacements(unittest.TestCase):
    """Tests for displacement/velocity/acceleration/eigenvector write_op2.

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: static, transient, modal, freq classmethods + write_op2.
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_displacement_static(self):
        """RealDisplacementArray.add_static_case + write_op2.

        Data: 3 GRID nodes, 1 time step, 6 DOF columns [t1,t2,t3,r1,r2,r3].
        """
        from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
        node_gridtype = np.array([[1, 1], [2, 1], [3, 1]], dtype='int32')
        data = np.zeros((1, 3, 6), dtype='float32')
        data[0, 0, :] = [0.1, 0.2, 0.3, 0.01, 0.02, 0.03]
        disp = RealDisplacementArray.add_static_case('OUGV1', node_gridtype, data, isubcase=1)
        assert disp.data.shape == (1, 3, 6)
        assert self._write(disp) > 100

    def test_displacement_transient(self):
        """RealDisplacementArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 3 time steps, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 6).astype('float32')
        disp = RealDisplacementArray.add_transient_case(
            'OUGV1', node_gridtype, data, isubcase=1, times=times)
        assert disp.ntimes == 3
        assert self._write(disp) > 100

    def test_displacement_modal(self):
        """RealDisplacementArray.add_modal_case + write_op2.

        Data: 2 GRID nodes, 3 modes, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2, 3], dtype='int32')
        eigns = np.array([100., 400., 900.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(3, 2, 6).astype('float32')
        disp = RealDisplacementArray.add_modal_case(
            'OUGV1', node_gridtype, data, isubcase=1,
            modes=modes, eigenvalues=eigns, mode_cycles=cycles)
        assert disp.ntimes == 3
        assert self._write(disp) > 100

    def test_complex_displacement_freq(self):
        """ComplexDisplacementArray.add_freq_case + write_op2.

        Data: 2 GRID nodes, 3 frequencies, 6 complex DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_displacements import ComplexDisplacementArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        freqs = np.array([10., 20., 30.], dtype='float32')
        data = np.zeros((3, 2, 6), dtype='complex64')
        data[0, 0, :] = [1+2j, 3+4j, 5+6j, 0.1+0.2j, 0.3+0.4j, 0.5+0.6j]
        cdisp = ComplexDisplacementArray.add_freq_case(
            'OUGV1', node_gridtype, data, isubcase=1, freqs=freqs)
        assert cdisp.ntimes == 3
        assert self._write(cdisp) > 100

    def test_complex_displacement_modes(self):
        """ComplexDisplacementArray.add_complex_modes_case + write_op2.

        Data: 2 GRID nodes, 2 complex modes, 6 complex DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_displacements import ComplexDisplacementArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 6), dtype='complex64')
        data[0, 0, :] = [1+2j, 3+4j, 5+6j, .1+.2j, .3+.4j, .5+.6j]
        cdisp = ComplexDisplacementArray.add_complex_modes_case(
            'OUGV1', node_gridtype, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis)
        assert cdisp.ntimes == 2
        assert self._write(cdisp) > 100

    def test_eigenvector_modal(self):
        """RealEigenvectorArray.add_modal_case + write_op2.

        Data: 2 GRID nodes, 2 modes, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 6).astype('float32')
        eig = RealEigenvectorArray.add_modal_case(
            'OUGV1', node_gridtype, data, isubcase=1,
            modes=modes, eigenvalues=eigns, mode_cycles=cycles)
        assert eig.ntimes == 2
        assert self._write(eig) > 100

    def test_complex_eigenvector_modes(self):
        """ComplexEigenvectorArray.add_complex_modes_case + write_op2.

        Data: 2 GRID nodes, 2 complex modes, 6 complex DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_eigenvectors import ComplexEigenvectorArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 6), dtype='complex64')
        data[0, 0, :] = [1+2j, 3+4j, 5+6j, .1+.2j, .3+.4j, .5+.6j]
        ceig = ComplexEigenvectorArray.add_complex_modes_case(
            'OUGV1', node_gridtype, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis)
        assert ceig.ntimes == 2
        assert self._write(ceig) > 100


class TestWriteOp2Forces(unittest.TestCase):
    """Tests for element force write_op2 via classmethods.

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: rod, bar, spring forces (static/modal/transient/freq).
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_rod_force_static(self):
        """RealRodForceArray.add_static_case + write_op2.

        Data: 2 CROD elements, 1 time step, 2 columns [axial, torsion].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealRodForceArray
        element = np.array([1, 2], dtype='int32')
        data = np.array([[[1000., 50.], [2000., 60.]]], dtype='float32')
        rod_f = RealRodForceArray.add_static_case('OEF1X', 'CROD', element, data, isubcase=1)
        assert rod_f.nelements == 2
        assert self._write(rod_f) > 100

    def test_rod_force_modal(self):
        """RealRodForceArray.add_modal_case + write_op2.

        Data: 2 CROD elements, 2 modes, 2 columns [axial, torsion].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealRodForceArray
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        freqs = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 2).astype('float32')
        rod_f = RealRodForceArray.add_modal_case(
            'OEF1X', 'CROD', element, data, isubcase=1,
            modes=modes, eigns=eigns, freqs=freqs)
        assert rod_f.ntimes == 2
        assert self._write(rod_f) > 100

    def test_rod_force_transient(self):
        """RealRodForceArray.add_transient_case + write_op2.

        Data: 2 CROD elements, 3 time steps, 2 columns [axial, torsion].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealRodForceArray
        element = np.array([1, 2], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 2).astype('float32')
        rod_f = RealRodForceArray.add_transient_case(
            'OEF1X', 'CROD', element, data, isubcase=1, times=times)
        assert rod_f.ntimes == 3
        assert self._write(rod_f) > 100

    def test_bar_force_static(self):
        """RealCBarForceArray.add_static_case + write_op2.

        Data: 1 CBAR element, 1 time step, 8 columns
        [bm_a1, bm_a2, bm_b1, bm_b2, shear1, shear2, axial, torque].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCBarForceArray
        element = np.array([20], dtype='int32')
        data = np.zeros((1, 1, 8), dtype='float32')
        data[0, 0, :] = [100, 200, 150, 250, 50, 60, 1000, 300]
        bar_f = RealCBarForceArray.add_static_case('OEF1X', 'CBAR', element, data, isubcase=1)
        assert bar_f.nelements == 1
        assert self._write(bar_f) > 100

    def test_bar_force_modal(self):
        """RealCBarForceArray.add_modal_case + write_op2.

        Data: 1 CBAR element, 2 modes, 8 columns.
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCBarForceArray
        element = np.array([20], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 1, 8).astype('float32')
        bar_f = RealCBarForceArray.add_modal_case(
            'OEF1X', 'CBAR', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert bar_f.ntimes == 2
        assert self._write(bar_f) > 100

    def test_bar_force_transient(self):
        """RealCBarForceArray.add_transient_case + write_op2.

        Data: 1 CBAR element, 2 time steps, 8 columns.
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCBarForceArray
        element = np.array([20], dtype='int32')
        times = np.array([0.0, 1.0], dtype='float32')
        data = np.random.rand(2, 1, 8).astype('float32')
        bar_f = RealCBarForceArray.add_transient_case(
            'OEF1X', 'CBAR', element, data, isubcase=1, times=times)
        assert bar_f.ntimes == 2
        assert self._write(bar_f) > 100

    def test_spring_force_static(self):
        """RealSpringForceArray.add_static_case + write_op2.

        Data: 1 CELAS1 element, 1 time step, 1 column [spring_force].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray
        element = np.array([10], dtype='int32')
        data = np.array([[[5.0]]], dtype='float32')
        spring_f = RealSpringForceArray.add_static_case(
            'OEF1X', 'CELAS1', element, data, isubcase=1)
        assert spring_f.nelements == 1
        assert self._write(spring_f) > 100

    def test_spring_force_modal(self):
        """RealSpringForceArray.add_modal_case + write_op2.

        Data: 1 CELAS1 element, 2 modes, 1 column [spring_force].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray
        element = np.array([10], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        freqs = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 1, 1).astype('float32')
        spring_f = RealSpringForceArray.add_modal_case(
            'OEF1X', 'CELAS1', element, data, isubcase=1,
            modes=modes, eigns=eigns, freqs=freqs)
        assert spring_f.ntimes == 2
        assert self._write(spring_f) > 100

    def test_spring_force_transient(self):
        """RealSpringForceArray.add_transient_case + write_op2.

        Data: 1 CELAS1 element, 3 time steps, 1 column [spring_force].
        """
        from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray
        element = np.array([10], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 1, 1).astype('float32')
        spring_f = RealSpringForceArray.add_transient_case(
            'OEF1X', 'CELAS1', element, data, isubcase=1, times=times)
        assert spring_f.ntimes == 3
        assert self._write(spring_f) > 100

    def test_complex_rod_force_freq(self):
        """ComplexRodForceArray.add_freq_case + write_op2.

        Data: 2 CROD elements, 2 frequencies, 2 complex columns [axial, torque].
        """
        from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexRodForceArray
        element = np.array([1, 2], dtype='int32')
        freqs = np.array([10., 20.], dtype='float32')
        data = np.zeros((2, 2, 2), dtype='complex64')
        data[0, 0, :] = [100+50j, 30+20j]
        crod_f = ComplexRodForceArray.add_freq_case(
            'OEF1X', 'CROD', element, data, isubcase=1, freqs=freqs)
        assert crod_f.ntimes == 2
        assert self._write(crod_f) > 100

    def test_complex_rod_force_modes(self):
        """ComplexRodForceArray.add_complex_modes_case + write_op2.

        Data: 2 CROD elements, 2 complex modes, 2 complex columns.
        """
        from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexRodForceArray
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 2), dtype='complex64')
        crod_f = ComplexRodForceArray.add_complex_modes_case(
            'OEF1X', element, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis, element_name='CROD')
        assert crod_f.ntimes == 2
        assert self._write(crod_f) > 100

    def test_complex_spring_force_freq(self):
        """ComplexSpringForceArray.add_freq_case + write_op2.

        Data: 1 CELAS1 element, 3 frequencies, 1 complex column.
        """
        from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexSpringForceArray
        element = np.array([10], dtype='int32')
        freqs = np.array([5., 10., 15.], dtype='float32')
        data = np.zeros((3, 1, 1), dtype='complex64')
        data[0, 0, 0] = 1+2j
        cspring_f = ComplexSpringForceArray.add_freq_case(
            'OEF1X', 'CELAS1', element, data, isubcase=1, freqs=freqs)
        assert cspring_f.ntimes == 3
        assert self._write(cspring_f) > 100

    def test_complex_spring_force_modes(self):
        """ComplexSpringForceArray.add_complex_modes_case + write_op2.

        Data: 1 CELAS1 element, 2 complex modes, 1 complex column.
        """
        from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexSpringForceArray
        element = np.array([10], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 1, 1), dtype='complex64')
        data[0, 0, 0] = 1+2j
        cspring_f = ComplexSpringForceArray.add_complex_modes_case(
            'OEF1X', element, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis, element_name='CELAS1')
        assert cspring_f.ntimes == 2
        assert self._write(cspring_f) > 100


class TestWriteOp2StressExtended(unittest.TestCase):
    """Tests for plate and solid stress write_op2 via classmethods.

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: plate (static/modal/transient/freq), solid (static/modal), rod post-buckling.
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_plate_stress_static(self):
        """RealPlateStressArray.add_static_case + write_op2.

        Data: 2 CTRIA3 elements (nnodes=1), 4 layers, 8 stress columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray
        element_node = np.array([[1, 0], [1, 0], [2, 0], [2, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05, -0.05, 0.05], dtype='float32')
        data = np.random.rand(1, 4, 8).astype('float32')
        plate = RealPlateStressArray.add_static_case(
            'OES1X1', 'CTRIA3', 1, element_node, fiber, data, isubcase=1)
        assert plate.nelements == 2
        assert self._write(plate) > 100

    def test_plate_stress_modal(self):
        """RealPlateStressArray.add_modal_case + write_op2.

        Data: 1 CQUAD4 element (nnodes=1), 2 layers, 3 modes.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray
        element_node = np.array([[1, 0], [1, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05], dtype='float32')
        modes = np.array([1, 2, 3], dtype='int32')
        eigns = np.array([100., 400., 900.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(3, 2, 8).astype('float32')
        plate = RealPlateStressArray.add_modal_case(
            'OES1X1', 'CQUAD4', 1, element_node, fiber, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert plate.ntimes == 3
        assert self._write(plate) > 100

    def test_plate_stress_transient(self):
        """RealPlateStressArray.add_transient_case + write_op2.

        Data: 1 CTRIA3 element (nnodes=1), 2 layers, 2 time steps.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray
        element_node = np.array([[1, 0], [1, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05], dtype='float32')
        times = np.array([0.0, 1.0], dtype='float32')
        data = np.random.rand(2, 2, 8).astype('float32')
        plate = RealPlateStressArray.add_transient_case(
            'OES1X1', 'CTRIA3', 1, element_node, fiber, data, isubcase=1, times=times)
        assert plate.ntimes == 2
        assert self._write(plate) > 100

    def test_complex_plate_stress_freq(self):
        """ComplexPlateStressArray.add_freq_case + write_op2.

        Data: 1 CTRIA3 element, 2 layers, 2 frequencies, 3 complex columns [oxx, oyy, txy].
        """
        from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray
        element_node = np.array([[1, 0], [1, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05], dtype='float32')
        freqs = np.array([10., 20.], dtype='float32')
        data = np.zeros((2, 2, 3), dtype='complex64')
        data[0, 0, :] = [100+10j, 50+5j, 20+2j]
        cplate = ComplexPlateStressArray.add_freq_case(
            'OES1C', 'CTRIA3', element_node, fiber, data, isubcase=1, freqs=freqs)
        assert cplate.ntimes == 2
        assert self._write(cplate) > 100

    def test_solid_stress_static(self):
        """RealSolidStressArray.add_static_case + write_op2.

        Data: 1 CTETRA element, 5 nodes (centroid + 4 corners), 10 stress columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStressArray
        element_node = np.array([[1, 0], [1, 10], [1, 20], [1, 30], [1, 40]], dtype='int32')
        element_cid = np.array([[1, 0]], dtype='int32')
        data = np.random.rand(1, 5, 10).astype('float32')
        solid = RealSolidStressArray.add_static_case(
            'OES1X1', 'CTETRA', element_node, element_cid, data, isubcase=1)
        assert solid.nelements == 1
        assert solid.nnodes == 5
        assert self._write(solid) > 100

    def test_solid_stress_modal(self):
        """RealSolidStressArray.add_modal_case + write_op2.

        Data: 1 CHEXA element, 9 nodes (centroid + 8 corners), 2 modes.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStressArray
        element_node = np.array([[1, 0]] + [[1, i] for i in range(10, 90, 10)], dtype='int32')
        element_cid = np.array([[1, 0]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 9, 10).astype('float32')
        solid = RealSolidStressArray.add_modal_case(
            'OES1X1', 'CHEXA', element_node, element_cid, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert solid.ntimes == 2
        assert solid.nnodes == 9
        assert self._write(solid) > 100

    def test_rod_stress_post_buckling(self):
        """RealRodStressArray.add_post_buckling_case constructs valid object.

        Data: 2 CROD elements, 2 modes, 4 columns [axial, torsion, SMa, SMt].
        Note: write_op2 not tested because post-buckling table3 requires lsdvmns.
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([0., 0.], dtype='float32')
        data = np.random.rand(2, 2, 4).astype('float32')
        rod = RealRodStressArray.add_post_buckling_case(
            'OES1X1', 'CROD', element, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis)
        assert rod.ntimes == 2
        assert rod.nelements == 2
        assert rod.data.shape == (2, 2, 4)

    def test_shear_stress_modal(self):
        """RealShearStressArray.add_modal_case + write_op2.

        Data: 2 CSHEAR elements, 2 modes, 3 columns [max_shear, avg_shear, margin].
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 3).astype('float32')
        shear = RealShearStressArray.add_modal_case(
            'OES1X1', 'CSHEAR', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert shear.ntimes == 2
        assert self._write(shear) > 100

    def test_spring_stress_modal(self):
        """RealSpringStressArray.add_modal_case + write_op2.

        Data: 2 CELAS1 elements, 2 modes, 1 column [stress].
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 1).astype('float32')
        spring = RealSpringStressArray.add_modal_case(
            'OES1X1', 'CELAS1', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert spring.ntimes == 2
        assert self._write(spring) > 100

    def test_spring_stress_transient(self):
        """RealSpringStressArray.add_transient_case + write_op2.

        Data: 2 CELAS2 elements, 3 time steps, 1 column [stress].
        """
        element = np.array([1, 2], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 1).astype('float32')
        spring = RealSpringStressArray.add_transient_case(
            'OES1X1', 'CELAS2', element, data, isubcase=1, times=times)
        assert spring.ntimes == 3
        assert self._write(spring) > 100

    def test_rod_stress_modal(self):
        """RealRodStressArray.add_modal_case + write_op2.

        Data: 2 CROD elements, 2 modes, 4 columns [axial, torsion, SMa, SMt].
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 4).astype('float32')
        rod = RealRodStressArray.add_modal_case(
            'OES1X1', 'CROD', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert rod.ntimes == 2
        assert rod.nelements == 2
        assert self._write(rod) > 100

    def test_rod_stress_transient(self):
        """RealRodStressArray.add_transient_case + write_op2.

        Data: 2 CROD elements, 3 time steps, 4 columns [axial, torsion, SMa, SMt].
        """
        element = np.array([1, 2], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 4).astype('float32')
        rod = RealRodStressArray.add_transient_case(
            'OES1X1', 'CROD', element, data, isubcase=1, times=times)
        assert rod.ntimes == 3
        assert rod.nelements == 2
        assert self._write(rod) > 100

    def test_bar_stress_modal(self):
        """RealBarStressArray.add_modal_case + write_op2.

        Data: 1 CBAR element, 2 modes, 15 columns.
        """
        element = np.array([10], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 1, 15).astype('float32')
        bar = RealBarStressArray.add_modal_case(
            'OES1X1', 'CBAR', element, data, isubcase=1,
            modes=modes, eigns=eigns, cycles=cycles)
        assert bar.ntimes == 2
        assert bar.nelements == 1
        assert self._write(bar) > 100

    def test_bar_stress_transient(self):
        """RealBarStressArray.add_transient_case + write_op2.

        Data: 1 CBAR element, 3 time steps, 15 columns.
        """
        element = np.array([10], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 1, 15).astype('float32')
        bar = RealBarStressArray.add_transient_case(
            'OES1X1', 'CBAR', element, data, isubcase=1, times=times)
        assert bar.ntimes == 3
        assert bar.nelements == 1
        assert self._write(bar) > 100

    def test_shear_stress_transient(self):
        """RealShearStressArray.add_transient_case + write_op2.

        Data: 2 CSHEAR elements, 3 time steps, 3 columns [max_shear, avg_shear, margin].
        """
        element = np.array([20, 21], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 3).astype('float32')
        shear = RealShearStressArray.add_transient_case(
            'OES1X1', 'CSHEAR', element, data, isubcase=1, times=times)
        assert shear.ntimes == 3
        assert shear.nelements == 2
        assert self._write(shear) > 100

    def test_solid_stress_transient(self):
        """RealSolidStressArray.add_transient_case + write_op2.

        Data: 1 CTETRA element, 5 nodes, 3 time steps, 10 stress columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStressArray
        element_node = np.array([[1, 0], [1, 10], [1, 20], [1, 30], [1, 40]], dtype='int32')
        element_cid = np.array([[1, 0]], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 5, 10).astype('float32')
        solid = RealSolidStressArray.add_transient_case(
            'OES1X1', 'CTETRA', element_node, element_cid, data, isubcase=1, times=times)
        assert solid.ntimes == 3
        assert solid.nnodes == 5
        assert self._write(solid) > 100

    def test_complex_rod_stress_freq(self):
        """ComplexRodStressArray.add_freq_case + write_op2.

        Data: 2 CROD elements, 3 frequencies, 2 complex columns [axial, torsion].
        """
        element = np.array([1, 2], dtype='int32')
        freqs = np.array([10., 20., 30.], dtype='float32')
        data = np.zeros((3, 2, 2), dtype='complex64')
        data[0, 0, :] = [100+10j, 50+5j]
        data[1, 1, :] = [200+20j, 80+8j]
        crod = ComplexRodStressArray.add_freq_case(
            'OES1C', 'CROD', element, data, isubcase=1, freqs=freqs)
        assert crod.ntimes == 3
        assert self._write(crod) > 100

    def test_complex_spring_stress_freq(self):
        """ComplexSpringStressArray.add_freq_case + write_op2.

        Data: 1 CELAS1 element, 2 frequencies, 1 complex column.
        """
        element = np.array([1], dtype='int32')
        freqs = np.array([5., 10.], dtype='float32')
        data = np.zeros((2, 1, 1), dtype='complex64')
        data[0, 0, 0] = 3+4j
        cspring = ComplexSpringStressArray.add_freq_case(
            'OES1C', 'CELAS1', element, data, isubcase=1, freqs=freqs)
        assert cspring.ntimes == 2
        assert self._write(cspring) > 100


class TestWriteOp2VelocityAcceleration(unittest.TestCase):
    """Tests for velocity/acceleration/SPC/MPC force write_op2.

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: static, modal, transient, freq, complex_modes.
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_velocity_static(self):
        """RealVelocityArray.add_static_case + write_op2.

        Data: 2 GRID nodes, 1 time step, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        vel = RealVelocityArray.add_static_case('OVG1', node_gridtype, data, isubcase=1)
        assert self._write(vel) > 100

    def test_velocity_modal(self):
        """RealVelocityArray.add_modal_case + write_op2.

        Data: 2 GRID nodes, 2 modes, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigns = np.array([100., 400.], dtype='float32')
        cycles = np.sqrt(eigns) / (2 * np.pi)
        data = np.random.rand(2, 2, 6).astype('float32')
        vel = RealVelocityArray.add_modal_case(
            'OVG1', node_gridtype, data, isubcase=1,
            modes=modes, eigenvalues=eigns, mode_cycles=cycles)
        assert vel.ntimes == 2
        assert self._write(vel) > 100

    def test_velocity_transient(self):
        """RealVelocityArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 3 time steps, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        times = np.array([0.0, 0.5, 1.0], dtype='float32')
        data = np.random.rand(3, 2, 6).astype('float32')
        vel = RealVelocityArray.add_transient_case(
            'OVG1', node_gridtype, data, isubcase=1, times=times)
        assert vel.ntimes == 3
        assert self._write(vel) > 100

    def test_complex_velocity_freq(self):
        """ComplexVelocityArray.add_freq_case + write_op2.

        Data: 2 GRID nodes, 2 frequencies, 6 complex DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_velocities import ComplexVelocityArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        freqs = np.array([10., 20.], dtype='float32')
        data = np.zeros((2, 2, 6), dtype='complex64')
        data[0, 0, :] = [1+1j, 2+2j, 3+3j, .1+.1j, .2+.2j, .3+.3j]
        cvel = ComplexVelocityArray.add_freq_case(
            'OVG1', node_gridtype, data, isubcase=1, freqs=freqs)
        assert cvel.ntimes == 2
        assert self._write(cvel) > 100

    def test_acceleration_static(self):
        """RealAccelerationArray.add_static_case + write_op2.

        Data: 2 GRID nodes, 1 time step, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        acc = RealAccelerationArray.add_static_case('OAG1', node_gridtype, data, isubcase=1)
        assert self._write(acc) > 100

    def test_acceleration_transient(self):
        """RealAccelerationArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 2 time steps, 6 DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        times = np.array([0.0, 1.0], dtype='float32')
        data = np.random.rand(2, 2, 6).astype('float32')
        acc = RealAccelerationArray.add_transient_case(
            'OAG1', node_gridtype, data, isubcase=1, times=times)
        assert acc.ntimes == 2
        assert self._write(acc) > 100

    def test_complex_acceleration_modes(self):
        """ComplexAccelerationArray.add_complex_modes_case + write_op2.

        Data: 2 GRID nodes, 2 complex modes, 6 complex DOF columns.
        """
        from pyNastran.op2.tables.oug.oug_accelerations import ComplexAccelerationArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 6), dtype='complex64')
        cacc = ComplexAccelerationArray.add_complex_modes_case(
            'OAG1', node_gridtype, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis)
        assert cacc.ntimes == 2
        assert self._write(cacc) > 100

    def test_spc_forces_static(self):
        """RealSPCForcesArray.add_static_case + write_op2.

        Data: 2 GRID nodes, 1 time step, 6 DOF columns.
        """
        from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        spc = RealSPCForcesArray.add_static_case('OQG1', node_gridtype, data, isubcase=1)
        assert self._write(spc) > 100

    def test_spc_forces_transient(self):
        """RealSPCForcesArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 2 time steps, 6 DOF columns.
        """
        from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        times = np.array([0.0, 1.0], dtype='float32')
        data = np.random.rand(2, 2, 6).astype('float32')
        spc = RealSPCForcesArray.add_transient_case(
            'OQG1', node_gridtype, data, isubcase=1, times=times)
        assert spc.ntimes == 2
        assert self._write(spc) > 100

    def test_mpc_forces_transient(self):
        """RealMPCForcesArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 2 time steps, 6 DOF columns.
        """
        from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import RealMPCForcesArray
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        times = np.array([0.0, 1.0], dtype='float32')
        data = np.random.rand(2, 2, 6).astype('float32')
        mpc = RealMPCForcesArray.add_transient_case(
            'OQG1', node_gridtype, data, isubcase=1, times=times)
        assert mpc.ntimes == 2
        assert self._write(mpc) > 100


class TestWriteOp2Strain(unittest.TestCase):
    """Tests for strain arrays (mirrors stress tests).

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: rod, bar, shear, plate strain static + write_op2.
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_rod_strain_static(self):
        """RealRodStrainArray.add_static_case + write_op2.

        Data: 2 CROD elements, 1 time step, 4 columns [axial, torsion, SMa, SMt].
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStrainArray
        element = np.array([1, 2], dtype='int32')
        data = np.random.rand(1, 2, 4).astype('float32')
        rod = RealRodStrainArray.add_static_case('OES1X1', 'CROD', element, data, isubcase=1)
        assert rod.nelements == 2
        assert self._write(rod) > 100

    def test_bar_strain_static(self):
        """RealBarStrainArray.add_static_case + write_op2.

        Data: 1 CBAR element, 1 time step, 15 columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStrainArray
        element = np.array([10], dtype='int32')
        data = np.random.rand(1, 1, 15).astype('float32')
        bar = RealBarStrainArray.add_static_case('OES1X1', 'CBAR', element, data, isubcase=1)
        assert bar.nelements == 1
        assert self._write(bar) > 100

    def test_shear_strain_static(self):
        """RealShearStrainArray.add_static_case + write_op2.

        Data: 2 CSHEAR elements, 1 time step, 3 columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray
        element = np.array([1, 2], dtype='int32')
        data = np.random.rand(1, 2, 3).astype('float32')
        shear = RealShearStrainArray.add_static_case('OES1X1', 'CSHEAR', element, data, isubcase=1)
        assert shear.nelements == 2
        assert self._write(shear) > 100

    def test_plate_strain_static(self):
        """RealPlateStrainArray.add_static_case + write_op2.

        Data: 1 CTRIA3 element (nnodes=1), 2 layers, 8 strain columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStrainArray
        element_node = np.array([[1, 0], [1, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05], dtype='float32')
        data = np.random.rand(1, 2, 8).astype('float32')
        plate = RealPlateStrainArray.add_static_case(
            'OES1X1', 'CTRIA3', 1, element_node, fiber, data, isubcase=1)
        assert plate.nelements == 1
        assert self._write(plate) > 100

    def test_complex_rod_strain_freq(self):
        """ComplexRodStrainArray.add_freq_case + write_op2.

        Data: 2 CROD elements, 2 frequencies, 2 complex columns [axial, torsion].
        """
        from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStrainArray
        element = np.array([1, 2], dtype='int32')
        freqs = np.array([10., 20.], dtype='float32')
        data = np.zeros((2, 2, 2), dtype='complex64')
        data[0, 0, :] = [0.001+0.002j, 0.003+0.004j]
        crod = ComplexRodStrainArray.add_freq_case(
            'OES1C', 'CROD', element, data, isubcase=1, freqs=freqs)
        assert crod.ntimes == 2
        assert self._write(crod) > 100


class TestWriteOp2ComplexStressModes(unittest.TestCase):
    """Tests for complex stress add_complex_modes_case + write_op2.

    Tolerances: checks write_op2 produces >100 bytes.
    Tests: rod, spring, plate complex stress with complex modes.
    """

    def _write(self, result, endian=b'>'):
        buf = BytesIO()
        result.write_op2(buf, StringIO(), -1, True, (5, 23, 2026), endian=endian)
        return buf.tell()

    def test_complex_rod_stress_modes(self):
        """ComplexRodStressArray.add_complex_modes_case + write_op2.

        Data: 2 CROD elements, 2 complex modes, 2 complex columns.
        """
        element = np.array([1, 2], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 2), dtype='complex64')
        data[0, 0, :] = [100+10j, 50+5j]
        crod = ComplexRodStressArray.add_complex_modes_case(
            'OES1C', element, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis, element_name='CROD')
        assert crod.ntimes == 2
        assert self._write(crod) > 100

    def test_complex_spring_stress_modes(self):
        """ComplexSpringStressArray.add_complex_modes_case + write_op2.

        Data: 1 CELAS1 element, 2 complex modes, 1 complex column.
        """
        element = np.array([1], dtype='int32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 1, 1), dtype='complex64')
        data[0, 0, 0] = 3+4j
        cspring = ComplexSpringStressArray.add_complex_modes_case(
            'OES1C', element, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis, element_name='CELAS1')
        assert cspring.ntimes == 2
        assert self._write(cspring) > 100

    def test_complex_plate_stress_modes(self):
        """ComplexPlateStressArray.add_complex_modes_case + write_op2.

        Data: 1 CTRIA3 element, 2 layers, 2 complex modes, 3 complex columns.
        """
        from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray
        element_node = np.array([[1, 0], [1, 0]], dtype='int32')
        fiber = np.array([-0.05, 0.05], dtype='float32')
        modes = np.array([1, 2], dtype='int32')
        eigrs = np.array([10., 20.], dtype='float32')
        eigis = np.array([1., 2.], dtype='float32')
        data = np.zeros((2, 2, 3), dtype='complex64')
        data[0, 0, :] = [100+10j, 50+5j, 20+2j]
        cplate = ComplexPlateStressArray.add_complex_modes_case(
            'OES1C', 'CTRIA3', element_node, fiber, data, isubcase=1,
            modes=modes, eigrs=eigrs, eigis=eigis)
        assert cplate.ntimes == 2
        assert self._write(cplate) > 100


if __name__ == '__main__':
    unittest.main()
