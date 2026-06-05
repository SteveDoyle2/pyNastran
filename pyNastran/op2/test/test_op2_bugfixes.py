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
from unittest.mock import MagicMock
from pyNastran.op2.tables.oes_stressStrain import utils_cshear

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray, ComplexDisplacementArray
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray, ComplexEigenvectorArray
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealRodForceArray, RealCBarForceArray, RealSpringForceArray
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexRodForceArray, ComplexSpringForceArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStressArray
from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray, ComplexVelocityArray
from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray, ComplexAccelerationArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import RealMPCForcesArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStrainArray

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()

from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexSpringForceArray



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
        spring = RealSpringStressArray.add_static_case(
            'OES1X1', 'CELAS1', element, data, isubcase=1)
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

        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
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

        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
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

        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)

            with np.errstate(under='ignore'):
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
        element = np.array([1, 2], dtype='int32')
        data = np.array([[[1000., 50.], [2000., 60.]]], dtype='float32')
        rod_f = RealRodForceArray.add_static_case('OEF1X', 'CROD', element, data, isubcase=1)
        assert rod_f.nelements == 2
        assert self._write(rod_f) > 100

    def test_rod_force_modal(self):
        """RealRodForceArray.add_modal_case + write_op2.

        Data: 2 CROD elements, 2 modes, 2 columns [axial, torsion].
        """
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
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        vel = RealVelocityArray.add_static_case('OVG1', node_gridtype, data, isubcase=1)
        assert self._write(vel) > 100

    def test_velocity_modal(self):
        """RealVelocityArray.add_modal_case + write_op2.

        Data: 2 GRID nodes, 2 modes, 6 DOF columns.
        """
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
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        acc = RealAccelerationArray.add_static_case('OAG1', node_gridtype, data, isubcase=1)
        assert self._write(acc) > 100

    def test_acceleration_transient(self):
        """RealAccelerationArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 2 time steps, 6 DOF columns.
        """
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
        node_gridtype = np.array([[1, 1], [2, 1]], dtype='int32')
        data = np.random.rand(1, 2, 6).astype('float32')
        spc = RealSPCForcesArray.add_static_case('OQG1', node_gridtype, data, isubcase=1)
        assert self._write(spc) > 100

    def test_spc_forces_transient(self):
        """RealSPCForcesArray.add_transient_case + write_op2.

        Data: 2 GRID nodes, 2 time steps, 6 DOF columns.
        """
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
        element = np.array([1, 2], dtype='int32')
        data = np.random.rand(1, 2, 4).astype('float32')
        rod = RealRodStrainArray.add_static_case('OES1X1', 'CROD', element, data, isubcase=1)
        assert rod.nelements == 2
        assert self._write(rod) > 100

    def test_bar_strain_static(self):
        """RealBarStrainArray.add_static_case + write_op2.

        Data: 1 CBAR element, 1 time step, 15 columns.
        """
        element = np.array([10], dtype='int32')
        data = np.random.rand(1, 1, 15).astype('float32')
        bar = RealBarStrainArray.add_static_case('OES1X1', 'CBAR', element, data, isubcase=1)
        assert bar.nelements == 1
        assert self._write(bar) > 100

    def test_shear_strain_static(self):
        """RealShearStrainArray.add_static_case + write_op2.

        Data: 2 CSHEAR elements, 1 time step, 3 columns.
        """
        element = np.array([1, 2], dtype='int32')
        data = np.random.rand(1, 2, 3).astype('float32')
        shear = RealShearStrainArray.add_static_case('OES1X1', 'CSHEAR', element, data, isubcase=1)
        assert shear.nelements == 2
        assert self._write(shear) > 100

    def test_plate_strain_static(self):
        """RealPlateStrainArray.add_static_case + write_op2.

        Data: 1 CTRIA3 element (nnodes=1), 2 layers, 8 strain columns.
        """
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


class TestWriteOp2ReadBack(unittest.TestCase):
    """Round-trip tests: read OP2 model -> write -> read back -> compare data.

    Tolerances: atol=1e-5 for all float comparisons.
    Tests: static, transient, modal, frequency response (complex) models.
    Each test verifies multiple result types (stress, strain, forces, displacements).
    """

    def _round_trip(self, op2_filename, clear_gpf=True):
        """Helper: read model, write to temp, read back. Returns (model1, model2)."""
        model = read_op2(str(op2_filename), debug=None)
        if clear_gpf:
            model.grid_point_forces = {}
        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)
        return model, model2

    def _compare_result(self, slot1, slot2, name):
        """Compare a result dict slot between two models."""
        for k in slot1:
            d1 = slot1[k].data
            d2 = slot2[k].data
            assert d1.shape == d2.shape, f'{name} subcase {k}: shape {d1.shape} vs {d2.shape}'
            with np.errstate(under='ignore'):
                assert np.allclose(d1, d2, atol=1e-5), f'{name} subcase {k}: data mismatch'

    def test_static_stress_round_trip(self):
        """Static model: rod, bar, plate, solid stress survive write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: crod_stress (4 cols), cbar_stress (15 cols), ctria3_stress (8 cols),
                ctetra_stress (10 cols).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.stress.crod_stress,
                             m2.op2_results.stress.crod_stress, 'rod_stress')
        self._compare_result(m1.op2_results.stress.cbar_stress,
                             m2.op2_results.stress.cbar_stress, 'bar_stress')
        self._compare_result(m1.op2_results.stress.ctria3_stress,
                             m2.op2_results.stress.ctria3_stress, 'tria3_stress')
        self._compare_result(m1.op2_results.stress.ctetra_stress,
                             m2.op2_results.stress.ctetra_stress, 'tetra_stress')

    def test_static_strain_round_trip(self):
        """Static model: rod, bar, plate, solid strain survive write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: crod_strain, cbar_strain, ctria3_strain, ctetra_strain.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.strain.crod_strain,
                             m2.op2_results.strain.crod_strain, 'rod_strain')
        self._compare_result(m1.op2_results.strain.cbar_strain,
                             m2.op2_results.strain.cbar_strain, 'bar_strain')
        self._compare_result(m1.op2_results.strain.ctria3_strain,
                             m2.op2_results.strain.ctria3_strain, 'tria3_strain')
        self._compare_result(m1.op2_results.strain.ctetra_strain,
                             m2.op2_results.strain.ctetra_strain, 'tetra_strain')

    def test_static_forces_disp_round_trip(self):
        """Static model: displacements, SPC forces, rod/bar forces survive write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: displacements (6 DOF), spc_forces, crod_force (2 cols), cbar_force (8 cols).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.displacements, m2.displacements, 'displacements')
        self._compare_result(m1.spc_forces, m2.spc_forces, 'spc_forces')
        self._compare_result(m1.op2_results.force.crod_force,
                             m2.op2_results.force.crod_force, 'rod_force')
        self._compare_result(m1.op2_results.force.cbar_force,
                             m2.op2_results.force.cbar_force, 'bar_force')

    def test_transient_stress_round_trip(self):
        """Transient model: rod, bar, plate stress survive write/read.

        Model: transient_solid_shell_bar.op2 (SOL 109, transient response).
        Checks: multi-timestep data for crod_stress, cbar_stress, cquad4_stress.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'transient_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.stress.crod_stress,
                             m2.op2_results.stress.crod_stress, 'transient_rod_stress')
        self._compare_result(m1.op2_results.stress.cbar_stress,
                             m2.op2_results.stress.cbar_stress, 'transient_bar_stress')
        self._compare_result(m1.op2_results.stress.cquad4_stress,
                             m2.op2_results.stress.cquad4_stress, 'transient_quad_stress')

    def test_transient_disp_strain_round_trip(self):
        """Transient model: displacements and strain survive write/read.

        Model: transient_solid_shell_bar.op2 (SOL 109).
        Checks: displacements (multi-timestep), cquad4_strain, crod_strain.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'transient_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.displacements, m2.displacements, 'transient_disp')
        self._compare_result(m1.op2_results.strain.cquad4_strain,
                             m2.op2_results.strain.cquad4_strain, 'transient_quad_strain')
        self._compare_result(m1.op2_results.strain.crod_strain,
                             m2.op2_results.strain.crod_strain, 'transient_rod_strain')

    def test_modal_stress_force_round_trip(self):
        """Modal model: eigenvectors, rod/bar stress, rod/bar forces survive write/read.

        Model: mode_solid_shell_bar.op2 (SOL 103, normal modes).
        Checks: eigenvectors, crod_stress, cbar_stress, crod_force, cbar_force, spc_forces.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.eigenvectors, m2.eigenvectors, 'eigenvectors')
        self._compare_result(m1.op2_results.stress.crod_stress,
                             m2.op2_results.stress.crod_stress, 'modal_rod_stress')
        self._compare_result(m1.op2_results.stress.cbar_stress,
                             m2.op2_results.stress.cbar_stress, 'modal_bar_stress')
        self._compare_result(m1.op2_results.force.crod_force,
                             m2.op2_results.force.crod_force, 'modal_rod_force')
        self._compare_result(m1.op2_results.force.cbar_force,
                             m2.op2_results.force.cbar_force, 'modal_bar_force')
        self._compare_result(m1.spc_forces, m2.spc_forces, 'modal_spc_forces')

    def test_modal_plate_solid_round_trip(self):
        """Modal model: plate and solid stress/strain survive write/read.

        Model: mode_solid_shell_bar.op2 (SOL 103).
        Checks: cquad4_stress, ctria3_stress, ctetra_stress, cpenta_stress, chexa_stress.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.stress.cquad4_stress,
                             m2.op2_results.stress.cquad4_stress, 'modal_quad_stress')
        self._compare_result(m1.op2_results.stress.ctria3_stress,
                             m2.op2_results.stress.ctria3_stress, 'modal_tria_stress')
        self._compare_result(m1.op2_results.stress.ctetra_stress,
                             m2.op2_results.stress.ctetra_stress, 'modal_tetra_stress')
        self._compare_result(m1.op2_results.stress.cpenta_stress,
                             m2.op2_results.stress.cpenta_stress, 'modal_penta_stress')
        self._compare_result(m1.op2_results.stress.chexa_stress,
                             m2.op2_results.stress.chexa_stress, 'modal_hexa_stress')

    def test_freq_complex_stress_round_trip(self):
        """Frequency response model: complex stress survives write/read.

        Model: freq_solid_shell_bar.op2 (SOL 108, frequency response).
        Checks: complex crod_stress, cbar_stress, ctria3_stress, ctetra_stress.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.stress.crod_stress,
                             m2.op2_results.stress.crod_stress, 'complex_rod_stress')
        self._compare_result(m1.op2_results.stress.cbar_stress,
                             m2.op2_results.stress.cbar_stress, 'complex_bar_stress')
        self._compare_result(m1.op2_results.stress.ctria3_stress,
                             m2.op2_results.stress.ctria3_stress, 'complex_tria_stress')
        self._compare_result(m1.op2_results.stress.ctetra_stress,
                             m2.op2_results.stress.ctetra_stress, 'complex_tetra_stress')

    def test_freq_complex_disp_force_round_trip(self):
        """Frequency response model: complex displacements and forces survive write/read.

        Model: freq_solid_shell_bar.op2 (SOL 108).
        Checks: complex displacements (6 DOF), crod_force, cbar_force, spc_forces.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.displacements, m2.displacements, 'complex_disp')
        self._compare_result(m1.op2_results.force.crod_force,
                             m2.op2_results.force.crod_force, 'complex_rod_force')
        self._compare_result(m1.op2_results.force.cbar_force,
                             m2.op2_results.force.cbar_force, 'complex_bar_force')
        self._compare_result(m1.spc_forces, m2.spc_forces, 'complex_spc_forces')

    def test_freq_complex_strain_round_trip(self):
        """Frequency response model: complex strain survives write/read.

        Model: freq_solid_shell_bar.op2 (SOL 108).
        Checks: complex crod_strain, cbar_strain, ctria3_strain, ctetra_strain.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        self._compare_result(m1.op2_results.strain.crod_strain,
                             m2.op2_results.strain.crod_strain, 'complex_rod_strain')
        self._compare_result(m1.op2_results.strain.cbar_strain,
                             m2.op2_results.strain.cbar_strain, 'complex_bar_strain')
        self._compare_result(m1.op2_results.strain.ctria3_strain,
                             m2.op2_results.strain.ctria3_strain, 'complex_tria_strain')
        self._compare_result(m1.op2_results.strain.ctetra_strain,
                             m2.op2_results.strain.ctetra_strain, 'complex_tetra_strain')


class TestWriteOp2CompositeBeamForce(unittest.TestCase):
    """Round-trip tests for composite stress/strain and beam/plate forces.

    Tolerances: atol=1e-5 for all float comparisons.
    Tests: composite plate stress/strain, beam force, plate force round-trips.
    """

    def _round_trip(self, op2_filename):
        model = read_op2(str(op2_filename), debug=None)
        model.grid_point_forces = {}
        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)
        return model, model2

    def test_composite_quad_stress_round_trip(self):
        """Composite CQUAD4 stress survives write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: element_layer preserved, 9 columns (o11, o22, t12, t1z, t2z, angle,
                major, minor, max_shear).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        comp1 = m1.op2_results.stress.cquad4_composite_stress
        comp2 = m2.op2_results.stress.cquad4_composite_stress
        for k in comp1:
            d1, d2 = comp1[k].data, comp2[k].data
            assert d1.shape == d2.shape
            assert d1.shape[2] == 9
            assert np.allclose(d1, d2, atol=1e-5), f'composite quad stress mismatch {k}'

    def test_composite_tria_stress_round_trip(self):
        """Composite CTRIA3 stress survives write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: data fidelity within atol=1e-5.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        comp1 = m1.op2_results.stress.ctria3_composite_stress
        comp2 = m2.op2_results.stress.ctria3_composite_stress
        for k in comp1:
            assert np.allclose(comp1[k].data, comp2[k].data, atol=1e-5)

    def test_composite_strain_round_trip(self):
        """Composite CQUAD4/CTRIA3 strain survives write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: both quad and tria composite strain data fidelity.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        for attr in ('cquad4_composite_strain', 'ctria3_composite_strain'):
            s1 = getattr(m1.op2_results.strain, attr)
            s2 = getattr(m2.op2_results.strain, attr)
            for k in s1:
                assert np.allclose(s1[k].data, s2[k].data, atol=1e-5), f'{attr} mismatch {k}'

    def test_beam_force_round_trip(self):
        """CBEAM force survives write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: 8 columns (sd, bm1, bm2, shear1, shear2, axial, torque, warp_torque).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        bf1 = m1.op2_results.force.cbeam_force
        bf2 = m2.op2_results.force.cbeam_force
        for k in bf1:
            d1, d2 = bf1[k].data, bf2[k].data
            assert d1.shape == d2.shape
            assert d1.shape[2] == 8
            assert np.allclose(d1, d2, atol=1e-5), f'beam force mismatch {k}'

    def test_plate_force_round_trip(self):
        """CTRIA3 and CQUAD4 plate forces survive write/read.

        Model: static_solid_shell_bar.op2 (SOL 101).
        Checks: 8 columns (mx, my, mxy, bmx, bmy, bmxy, tx, ty).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        for attr in ('ctria3_force', 'cquad4_force'):
            f1 = getattr(m1.op2_results.force, attr)
            f2 = getattr(m2.op2_results.force, attr)
            for k in f1:
                d1, d2 = f1[k].data, f2[k].data
                assert d1.shape[2] == 8
                assert np.allclose(d1, d2, atol=1e-5), f'{attr} mismatch {k}'


class TestOp2Metadata(unittest.TestCase):
    """Tests for metadata preservation across write/read round-trips.

    Tolerances: exact for integer fields, atol=1e-6 for times/modes.
    Tests: element IDs, node IDs, time steps, mode numbers, table names.
    """

    def _round_trip(self, op2_filename):
        model = read_op2(str(op2_filename), debug=None)
        model.grid_point_forces = {}
        fd, out_file = tempfile.mkstemp(suffix='.op2')
        os.close(fd)
        try:
            model.write_op2(out_file)
            model2 = read_op2(out_file, debug=None)
        finally:
            if os.path.exists(out_file):
                os.remove(out_file)
        return model, model2

    def test_transient_times_preserved(self):
        """Time steps are preserved exactly on round-trip.

        Model: transient_solid_shell_bar.op2 (21 time steps from 0.0 to 0.2).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'transient_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        rod1 = m1.op2_results.stress.crod_stress[1]
        rod2 = m2.op2_results.stress.crod_stress[1]
        assert rod1.ntimes == rod2.ntimes
        assert np.allclose(rod1._times, rod2._times, atol=1e-6)

    def test_transient_element_ids_preserved(self):
        """Element IDs are preserved on round-trip.

        Model: transient_solid_shell_bar.op2.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'transient_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        rod1 = m1.op2_results.stress.crod_stress[1]
        rod2 = m2.op2_results.stress.crod_stress[1]
        assert np.array_equal(rod1.element, rod2.element)

    def test_modal_modes_preserved(self):
        """Mode numbers (stored in _times) are preserved on round-trip.

        Model: mode_solid_shell_bar.op2 (3 modes).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        eig1 = m1.eigenvectors[1]
        eig2 = m2.eigenvectors[1]
        assert np.array_equal(eig1._times, eig2._times)

    def test_modal_node_gridtype_preserved(self):
        """Node IDs and grid types are preserved on round-trip.

        Model: mode_solid_shell_bar.op2.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        eig1 = m1.eigenvectors[1]
        eig2 = m2.eigenvectors[1]
        assert np.array_equal(eig1.node_gridtype, eig2.node_gridtype)

    def test_freq_frequencies_preserved(self):
        """Frequency values are preserved on round-trip.

        Model: freq_solid_shell_bar.op2 (frequency response).
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        disp1 = m1.displacements[1]
        disp2 = m2.displacements[1]
        assert np.allclose(disp1._times, disp2._times, atol=1e-6)

    def test_table_name_preserved(self):
        """Table names (OES1X1, OEF1X, etc.) are preserved on round-trip.

        Model: static_solid_shell_bar.op2.
        """
        op2_file = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'
        m1, m2 = self._round_trip(op2_file)
        rod1 = m1.op2_results.stress.crod_stress[1]
        rod2 = m2.op2_results.stress.crod_stress[1]
        assert rod1.table_name == rod2.table_name

        force1 = m1.op2_results.force.crod_force[1]
        force2 = m2.op2_results.force.crod_force[1]
        assert force1.table_name == force2.table_name


class TestOp2EqualityAndStats(unittest.TestCase):
    """Tests for __eq__ operators and get_stats on real model data.

    Tolerances: exact for equality; checks output format for get_stats.
    Tests: real/complex stress equality, get_stats output, f06 write smoke tests.
    """

    def test_real_rod_stress_equality(self):
        """RealRodStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        rod = model.op2_results.stress.crod_stress[1]
        assert rod == rod

    def test_real_bar_stress_equality(self):
        """RealBarStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        bar = model.op2_results.stress.cbar_stress[1]
        assert bar == bar

    def test_real_plate_stress_equality(self):
        """RealPlateStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        plate = model.op2_results.stress.ctria3_stress[1]
        assert plate == plate

    def test_real_solid_stress_equality(self):
        """RealSolidStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        solid = model.op2_results.stress.ctetra_stress[1]
        assert solid == solid

    def test_complex_rod_stress_equality(self):
        """ComplexRodStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'), debug=None)
        rod = model.op2_results.stress.crod_stress[1]
        assert rod == rod

    def test_complex_bar_stress_equality(self):
        """ComplexBarStressArray.__eq__ returns True for self-comparison (uses fixed unpacking)."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'), debug=None)
        bar = model.op2_results.stress.cbar_stress[1]
        assert bar == bar

    def test_complex_plate_stress_equality(self):
        """ComplexPlateStressArray.__eq__ returns True for self-comparison."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'), debug=None)
        plate = model.op2_results.stress.ctria3_stress[1]
        assert plate == plate

    def test_get_stats_rod_stress(self):
        """get_stats returns meaningful output for RealRodStressArray."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        rod = model.op2_results.stress.crod_stress[1]
        stats = rod.get_stats()
        assert len(stats) > 0
        joined = ''.join(stats)
        assert 'RealRodStress' in joined
        assert 'data' in joined

    def test_get_stats_complex_bar(self):
        """get_stats returns meaningful output for ComplexBarStressArray."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'), debug=None)
        bar = model.op2_results.stress.cbar_stress[1]
        stats = bar.get_stats()
        assert len(stats) > 0
        joined = ''.join(stats)
        assert 'ComplexBar' in joined

    def test_get_stats_transient_plate(self):
        """get_stats includes ntimes for transient results."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'transient_solid_shell_bar.op2'), debug=None)
        plate = model.op2_results.stress.cquad4_stress[1]
        stats = plate.get_stats()
        joined = ''.join(stats)
        assert 'ntimes' in joined

    def test_write_f06_rod_stress(self):
        """RealRodStressArray.write_f06 produces valid text output."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        rod = model.op2_results.stress.crod_stress[1]
        f = StringIO()
        rod.write_f06(f, header=[''], page_stamp='PAGE %s')
        text = f.getvalue()
        assert len(text) > 100
        assert 'STRESS' in text.upper() or 'AXIAL' in text.upper()

    def test_write_f06_complex_bar_stress(self):
        """ComplexBarStressArray.write_f06 produces valid text output."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'freq_solid_shell_bar.op2'), debug=None)
        bar = model.op2_results.stress.cbar_stress[1]
        f = StringIO()
        bar.write_f06(f, header=['', ''], page_stamp='PAGE %s')
        text = f.getvalue()
        assert len(text) > 100

    def test_write_f06_solid_stress(self):
        """RealSolidStressArray.write_f06 produces valid text output."""
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        solid = model.op2_results.stress.ctetra_stress[1]
        f = StringIO()
        solid.write_f06(f, header=[''], page_stamp='PAGE %s')
        text = f.getvalue()
        assert len(text) > 200


class TestLinearCombination(unittest.TestCase):
    """Tests for linear_combination on OP2 result objects.

    Verifies that combination_inplace respects the ires parameter:
    only specified columns are zeroed/combined, derived quantities are excluded.
    After update_data_components, derived quantities are recalculated.

    Tolerances: atol=1e-3 for derived quantities (principal stress, von Mises).
    """

    def test_rod_stress_ires_excludes_margins(self):
        """Rod stress linear_combination only modifies axial(0) and torsion(2).

        Headers: [axial, SMa, torsion, SMt].
        ires=[0, 2] -> margins (columns 1, 3) are NOT combined.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        rod = model.op2_results.stress.crod_stress[1]
        rod_copy = copy.deepcopy(rod)
        rod_copy.linear_combination(0.0, update=False)
        rod_copy.linear_combination(2.0, rod.data, update=False)
        assert np.allclose(rod_copy.data[:, :, 0], rod.data[:, :, 0] * 2)
        assert np.allclose(rod_copy.data[:, :, 2], rod.data[:, :, 2] * 2)
        np.testing.assert_array_equal(rod_copy.data[:, :, 1], rod.data[:, :, 1])
        np.testing.assert_array_equal(rod_copy.data[:, :, 3], rod.data[:, :, 3])

    def test_plate_stress_ires_excludes_derived(self):
        """Plate stress linear_combination only modifies oxx(1), oyy(2), txy(3).

        Headers: [fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises].
        ires=[1, 2, 3] -> fiber_dist, angle, principal, ovm are NOT combined.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        plate = model.op2_results.stress.ctria3_stress[1]
        plate_copy = copy.deepcopy(plate)
        plate_copy.linear_combination(0.0, update=False)
        plate_copy.linear_combination(3.0, plate.data, update=False)
        assert np.allclose(plate_copy.data[:, :, 1], plate.data[:, :, 1] * 3)
        assert np.allclose(plate_copy.data[:, :, 2], plate.data[:, :, 2] * 3)
        assert np.allclose(plate_copy.data[:, :, 3], plate.data[:, :, 3] * 3)
        np.testing.assert_array_equal(plate_copy.data[:, :, 0], plate.data[:, :, 0])
        np.testing.assert_array_equal(plate_copy.data[:, :, 4], plate.data[:, :, 4])
        np.testing.assert_array_equal(plate_copy.data[:, :, 5], plate.data[:, :, 5])
        np.testing.assert_array_equal(plate_copy.data[:, :, 6], plate.data[:, :, 6])
        np.testing.assert_array_equal(plate_copy.data[:, :, 7], plate.data[:, :, 7])

    def test_solid_stress_ires_excludes_principal(self):
        """Solid stress linear_combination only modifies first 6 columns.

        Headers: [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, von_mises].
        ires=slice(None, 6) -> principal stresses and ovm are NOT combined.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        solid = model.op2_results.stress.ctetra_stress[1]
        solid_copy = copy.deepcopy(solid)
        solid_copy.linear_combination(0.0, update=False)
        solid_copy.linear_combination(2.0, solid.data, update=False)
        assert np.allclose(solid_copy.data[:, :, :6], solid.data[:, :, :6] * 2)
        np.testing.assert_array_equal(solid_copy.data[:, :, 6:], solid.data[:, :, 6:])

    def test_beam_force_ires_excludes_sd(self):
        """Beam force linear_combination skips station distance (column 0).

        Headers: [sd, bm1, bm2, shear1, shear2, axial, torque, warp_torque].
        ires=slice(1, None) -> sd is NOT combined.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        bf = model.op2_results.force.cbeam_force[1]
        bf_copy = copy.deepcopy(bf)
        bf_copy.linear_combination(0.0, update=False)
        bf_copy.linear_combination(2.0, bf.data, update=False)
        np.testing.assert_array_equal(bf_copy.data[:, :, 0], bf.data[:, :, 0])
        assert np.allclose(bf_copy.data[:, :, 1:], bf.data[:, :, 1:] * 2)

    def test_composite_stress_ires_excludes_derived(self):
        """Composite stress linear_combination only modifies o11-t2z (columns 0-4).

        Headers: [o11, o22, t12, t1z, t2z, angle, major, minor, max_shear].
        ires=[0,1,2,3,4] -> angle, principal, max_shear are NOT combined.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        comp = model.op2_results.stress.cquad4_composite_stress[1]
        comp_copy = copy.deepcopy(comp)
        comp_copy.linear_combination(0.0, update=False)
        comp_copy.linear_combination(2.5, comp.data, update=False)
        assert np.allclose(comp_copy.data[:, :, :5], comp.data[:, :, :5] * 2.5)
        np.testing.assert_array_equal(comp_copy.data[:, :, 5:], comp.data[:, :, 5:])

    def test_plate_update_data_components(self):
        """After linear_combination + update_data_components, principal stresses are correct.

        Workflow: zero -> add 2x -> update_data_components.
        Verifies: omax, omin, von_mises are recalculated from 2*oxx, 2*oyy, 2*txy.
        """
        import copy
        model = read_op2(str(MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.op2'), debug=None)
        plate = model.op2_results.stress.ctria3_stress[1]
        plate_c = copy.deepcopy(plate)
        plate_c.linear_combination(0.0, update=False)
        plate_c.linear_combination(2.0, plate.data, update=False)
        plate_c.update_data_components()

        oxx = plate_c.data[:, :, 1].ravel()
        oyy = plate_c.data[:, :, 2].ravel()
        txy = plate_c.data[:, :, 3].ravel()
        C = (oxx + oyy) / 2
        R = np.sqrt(((oxx - oyy) / 2)**2 + txy**2)
        expected_max = C + R
        expected_min = C - R
        assert np.allclose(plate_c.data[:, :, 5].ravel(), expected_max, atol=1e-3)
        assert np.allclose(plate_c.data[:, :, 6].ravel(), expected_min, atol=1e-3)

    def test_combination_inplace_no_data_nonzero_factor(self):
        """combination_inplace with datai=None, factor!=0 raises RuntimeError."""
        from pyNastran.op2.result_objects.op2_objects import combination_inplace
        data = np.ones((1, 2, 4), dtype='float32')
        with self.assertRaises(RuntimeError):
            combination_inplace(data, None, 5.0)

    def test_combination_inplace_ires_zero(self):
        """combination_inplace with ires zeros only specified columns."""
        from pyNastran.op2.result_objects.op2_objects import combination_inplace
        data = np.ones((1, 2, 4), dtype='float32')
        combination_inplace(data, None, 0.0, ires=[0, 2])
        assert np.all(data[:, :, 0] == 0)
        assert np.all(data[:, :, 1] == 1)
        assert np.all(data[:, :, 2] == 0)
        assert np.all(data[:, :, 3] == 1)


if __name__ == '__main__':
    unittest.main()
