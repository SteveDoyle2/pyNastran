from __future__ import print_function
from six import  iteritems
from six.moves import zip
import os
import unittest
from math import sqrt

from numpy import array, array_equiv, array_equal, allclose
from itertools import count


class DummyWriter(object):
    def __init__(self, *args):
        pass
    def open(self, *args, **kwargs):
        pass
    def close(self):
        pass
    def write(self, msg):
        raise RuntimeError('remove this')

#import sys
#sys.stdout = DummyWriter()
import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.errors import FatalError
from pyNastran.op2.op2 import OP2
from pyNastran.op4.op4 import OP4

test_path = os.path.join(pyNastran.__path__[0], 'f06', 'test', 'tests')
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')

def run_model(bdf_name=None, op2_name=None, f06_name=None,
              op4_name=None, dynamic_vars=None, f06_has_weight=True,
              vectorized=False, encoding=None):
    outputs = []
    if bdf_name:
        bdf = BDF(debug=False, log=None)
        if dynamic_vars is not None:
            bdf.set_dynamic_syntax(dynamic_vars)
        bdf.read_bdf(bdf_name, encoding=encoding)
        bdf.write_bdf(bdf_name+'.out', interspersed=False)
        bdf.write_bdf(bdf_name+'.out', interspersed=True)
        outputs.append(bdf)

    if op2_name:
        op2 = OP2(debug=False)
        op2.read_op2(op2_name, vectorized=vectorized)
        op2.write_f06(op2_name[:-4] + '.test_op2.out', quiet=True)
        outputs.append(op2)

    if f06_name:
        f06 = F06(debug=False, log=None)
        f06.read_f06(f06_name)
        f06.write_f06(f06_name[:-4] + '.test_f06.out', quiet=True)
        outputs.append(f06)
        if f06_has_weight:
            assert f06.grid_point_weight.reference_point is not None
        else:
            assert f06.grid_point_weight.reference_point is None

    if op4_name:
        op4 = OP4()
        op4.read_op4(op4_name, matrixNames=None, precision='default')
        outputs.append(op4)

    assert len(outputs) > 0
    if len(outputs) == 1: return outputs[0]
    return outputs

class TestF06(unittest.TestCase):

    def test_blade2dv_fatal_1(self):
        f06_filename = os.path.join(model_path, 'blade_2dv', 'blade_2dv.f06_fatal')
        f06 = F06(debug=False, log=None)
        with self.assertRaises(AttributeError):
            f06.readF06(f06_filename)
        with self.assertRaises(FatalError):
            f06.read_f06(f06_filename)

    def test_blade2dv_fatal_2(self):
        f06_filename = os.path.join(model_path, 'blade_2dv', 'blade_2dv.f06_fatal')
        bdf_filename = os.path.join(model_path, 'blade_2dv', 'blade_2dv.bdf')
        #bdf2 = run_model(bdfname2, dynamic_vars=dynamic_vars)
        #self.assertEqual(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)

        f06 = F06(debug=False, log=None)

        try:  # this is supposed to raise a FatalError
            f06.read_f06(f06_filename)
            assert False, 'a FATAL should have been raised'
        except FatalError:
            pass

        f06.write_f06(f06_filename + '.out', quiet=True)

        ref_point = f06.grid_point_weight.reference_point
        #print("ref_point = %s" % ref_point)
        ref_point_exact = 0
        self.assertEqual(ref_point, ref_point_exact)

        MO = f06.grid_point_weight.MO
        #print("MO %s" % MO)
        MO_exact = array(
            [[  1.22085800e-01,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  5.33146300e-01,  -1.22767700e-05],
             [  0.00000000e+00,   1.22085800e-01,   0.00000000e+00,  -5.33146300e-01,  0.00000000e+00,   1.57186600e-01],
             [  0.00000000e+00,   0.00000000e+00,   1.22085800e-01,   1.22767700e-05, -1.57186600e-01,   0.00000000e+00],
             [  0.00000000e+00,  -5.33146300e-01,   1.22767700e-05,   3.20227600e+00, -7.13340800e-06,  -6.83890800e-01],
             [  5.33146300e-01,   0.00000000e+00,  -1.57186600e-01,  -7.13340800e-06,  3.45033400e+00,  -7.35886500e-05],
             [ -1.22767700e-05,   1.57186600e-01,   0.00000000e+00,  -6.83890800e-01, -7.35886500e-05,   2.50287600e-01]])


        S = f06.grid_point_weight.S
        S_exact = array([[ 1.,  0.,  0.],
                         [ 0.,  1.,  0.],
                         [ 0.,  0.,  1.]])
        #print("S %s" % S)

        mass = f06.grid_point_weight.mass
        mass_exact = array([ 0.1220858,  0.1220858,  0.1220858])
        self.assertTrue(array_equiv(mass, mass_exact))
        #print("mass = %s" % mass)

        cg = f06.grid_point_weight.cg
        cg_exact = array(
            [[  0.00000000e+00,   1.00558600e-04,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   0.00000000e+00,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   1.00558600e-04,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00]])
        #print("cg = %s", cg)
        self.assertTrue(array_equiv(cg, cg_exact))

        IS = f06.grid_point_weight.IS
        #print("IS = %s" % IS)
        IS_exact = array([[  8.74036600e-01,  -8.67305300e-06,  -2.54028500e-03],
                          [ -8.67305300e-06,   9.19714300e-01,   1.99762300e-05],
                          [ -2.54028500e-03,   1.99762300e-05,   4.79082500e-02]])

        IQ = f06.grid_point_weight.IQ
        #print("IQ %s" % IQ)
        IQ_exact = array([[ 0.04790044, 0.9197143, 0.8740444 ]])
        msg = 'IQ=%s\nexact=%s' % (str(IQ), str(IQ_exact))
        self.assertTrue(array_equiv(IQ, IQ_exact), msg=msg)

    def test_blade2dv_fatal_3(self):
        f06_filename = os.path.join(model_path, 'blade_2dv', 'blade_2dv.f06_fatal')
        bdf_filename = os.path.join(model_path, 'blade_2dv', 'blade_2dv.bdf')
        #bdf2 = run_model(bdfname2, dynamic_vars=dynamic_vars)
        #self.assertEqual(bdf2.properties[1].t, 42., 't=%s' % bdf2.properties[1].t)
        f06 = F06(debug=False, log=None)

        # we skip the fatal by stopping after reading the matrices
        f06.stop_after_reading_grid_point_weight()
        f06.read_f06(f06_filename)

        f06.write_f06(f06_filename + '.out', quiet=True)

        ref_point = f06.grid_point_weight.reference_point
        #print("ref_point = %s" % ref_point)
        ref_point_exact = 0
        self.assertEqual(ref_point, ref_point_exact)

        MO = f06.grid_point_weight.MO
        #print("MO = %s" % MO)
        MO_exact = array(
            [[  1.22085800e-01,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  5.33146300e-01,  -1.22767700e-05],
             [  0.00000000e+00,   1.22085800e-01,   0.00000000e+00,  -5.33146300e-01,  0.00000000e+00,   1.57186600e-01],
             [  0.00000000e+00,   0.00000000e+00,   1.22085800e-01,   1.22767700e-05, -1.57186600e-01,   0.00000000e+00],
             [  0.00000000e+00,  -5.33146300e-01,   1.22767700e-05,   3.20227600e+00, -7.13340800e-06,  -6.83890800e-01],
             [  5.33146300e-01,   0.00000000e+00,  -1.57186600e-01,  -7.13340800e-06,  3.45033400e+00,  -7.35886500e-05],
             [ -1.22767700e-05,   1.57186600e-01,   0.00000000e+00,  -6.83890800e-01, -7.35886500e-05,   2.50287600e-01]])


        S = f06.grid_point_weight.S
        S_exact = array([[ 1.,  0.,  0.],
                         [ 0.,  1.,  0.],
                         [ 0.,  0.,  1.]])
        #print("S %s" % S)

        mass = f06.grid_point_weight.mass
        mass_exact = array([ 0.1220858,  0.1220858,  0.1220858])
        self.assertTrue(array_equiv(mass, mass_exact))
        #print("mass = %s" % mass)

        cg = f06.grid_point_weight.cg
        cg_exact = array(
            [[  0.00000000e+00,   1.00558600e-04,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   0.00000000e+00,   4.36698100e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  1.28751000e+00,   1.00558600e-04,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00]])
        #print("cg = %s" % cg)
        self.assertTrue(array_equiv(cg, cg_exact))

        IS = f06.grid_point_weight.IS
        #print("IS  %s" % IS)
        IS_exact = array([[  8.74036600e-01,  -8.67305300e-06,  -2.54028500e-03],
                          [ -8.67305300e-06,   9.19714300e-01,   1.99762300e-05],
                          [ -2.54028500e-03,   1.99762300e-05,   4.79082500e-02]])

        IQ = f06.grid_point_weight.IQ
        #print("IQ %s" % IQ)
        IQ_exact = array([[ 0.04790044, 0.9197143, 0.8740444 ]])
        self.assertTrue(array_equiv(IQ, IQ_exact))

    def test_complex_tets_1(self):
        bdfname = None
        bdfname = os.path.join(model_path, 'complex', 'tet10', 'Simple_Example.bdf')
        f06name = os.path.join(model_path, 'complex', 'tet10', 'simple_example.f06')
        op2name = os.path.join(model_path, 'complex', 'tet10', 'simple_example.op2')

        f06name2 = os.path.join(model_path, 'complex', 'tet10', 'simple_example.test_f06.f06')
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)

        assert len(f06.displacements) == 1, len(f06.displacements)
        assert len(f06.spc_forces) == 0, len(f06.spc_forces) # 0 is correct

        assert len(f06.ctetra_strain) == 0, len(f06.ctetra_strain)  # 1 is correct
        assert len(f06.ctetra_stress) == 0, len(f06.ctetra_stress)  # 1 is correct
        assert len(f06.cpenta_strain) == 0, len(f06.cpenta_strain)
        assert len(f06.cpenta_stress) == 0, len(f06.cpenta_stress)
        assert len(f06.chexa_strain) == 0, len(f06.chexa_strain)
        assert len(f06.chexa_stress) == 0, len(f06.chexa_stress)

    def test_beam_modes_1(self):
        bdfname = os.path.join(model_path, 'beam_modes', 'beam_modes.dat')
        op2name = os.path.join(model_path, 'beam_modes', 'beam_modes_m1.op2')
        f06name = os.path.join(model_path, 'beam_modes', 'beam_modes.f06')
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=True)

        assert op2.Title == 'SIMPLE BEAM EXAMPLE', '%r' % op2.Title
        assert f06.Title == 'SIMPLE BEAM EXAMPLE', '%r' % f06.Title

        subtitle_label= f06.iSubcaseNameMap[1]
        assert subtitle_label[0] == 'MODES', subtitle_label
        assert subtitle_label[1] == 2, subtitle_label  # 2=modal
        assert subtitle_label[2] == '', subtitle_label

        subtitle_label = op2.iSubcaseNameMap[1]
        assert subtitle_label[0] == 'MODES', subtitle_label
        assert subtitle_label[1] == 2, subtitle_label  # 2=modal
        assert subtitle_label[2] == '', subtitle_label

        assert len(f06.displacements) == 0, len(f06.displacements)
        assert len(f06.eigenvectors) == 1, len(f06.eigenvectors)

        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_beam_modes_2(self):
        bdfname = None
        op2name = os.path.join(model_path, 'beam_modes', 'beam_modes_m2.op2')
        f06name = None
        op2 = run_model(bdfname, op2name, f06name, f06_has_weight=True)

        assert op2.Title == 'SIMPLE BEAM EXAMPLE', '%r' % op2.Title
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    #@unittest.expectedFailure - fails if subtitle check is different
    def test_bar3truss_1(self):
        bdfname = None
        op2name = None
        titles = ['', '', '', 'THIS_IS_A_BAD_TITLE', '']
        subtitles = ['UNTITLED.SC4', '', 'UNTITLED.SC4', 'UNTITLED.SC4', 'UNTITLED.SC4']
        labels = ['', '', 'MYLABEL', 'MYLABEL', '']
        f06_filenames = ['no_subcase.f06', 'subcase_no_subtitle.f06',
                         'subcase_subtitle_label.f06', 'subcase_subtitle_label_title.f06',
                         'with_subcase.f06']
        for i, title, subtitle, label, f06_filename in zip(count(), titles, subtitles, labels, f06_filenames):
            f06name = os.path.join(model_path, 'bar3truss', f06_filename)

            #f06name2 = os.path.join(model_path, 'bar3truss', 'no_subcase.test_f06.f06')
            f06 = run_model(bdfname, op2name, f06name, f06_has_weight=True)

            assert f06.Title == title, 'i=%i title=%r expected=%r' % (i, f06.Title, title)
            subtitle_label = f06.iSubcaseNameMap[1]
            analysis_code = 1 # static
            assert subtitle_label[0] == subtitle, 'f06_filename=%s i=%i subtitle=%r expected=%r' % (f06_filename, i, subtitle_label[0], subtitle)
            assert subtitle_label[1] == 1, 'f06_filename=%s i=%i label=%r expected=%r' % (f06_filename, i, subtitle_label[1], analysis_code)
            assert subtitle_label[2] == label, 'f06_filename=%s i=%i label=%r expected=%r' % (f06_filename, i, subtitle_label[2], label)

            assert len(f06.displacements) == 1, len(f06.displacements)
            assert len(f06.spc_forces) == 1, len(f06.spc_forces)
            assert len(f06.applied_loads) == 0, len(f06.applied_loads)  # 0 is correct
            assert len(f06.load_vectors) == 1, len(f06.load_vectors)
            assert len(f06.spc_forces) == 1, len(f06.spc_forces)

            assert len(f06.crod_force) == 1, len(f06.crod_force)
            assert len(f06.crod_strain) == 0, len(f06.crod_strain)  # 0 is correct
            assert len(f06.crod_stress) == 1, len(f06.crod_stress)

            # 0 is correct
            assert len(f06.ctetra_strain) == 0, len(f06.ctetra_strain)
            assert len(f06.ctetra_stress) == 0, len(f06.ctetra_stress)
            assert len(f06.cpenta_strain) == 0, len(f06.cpenta_strain)
            assert len(f06.cpenta_stress) == 0, len(f06.cpenta_stress)
            assert len(f06.chexa_strain) == 0, len(f06.chexa_strain)
            assert len(f06.chexa_stress) == 0, len(f06.chexa_stress)


    def test_fsi_1(self):
        bdfname = os.path.join(model_path, 'fsi', 'fsi.bdf')
        f06name = os.path.join(model_path, 'fsi', 'fsi.f06')
        op2name = os.path.join(model_path, 'fsi', 'fsi.op2')

        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)
        assert len(f06.eigenvectors) == 1, len(f06.eigenvectors)  # 1 is correct
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)  # 1 is correct

    def test_cbush_1(self):
        bdfname = os.path.join(model_path, 'cbush', 'cbush.dat')
        f06name = os.path.join(model_path, 'cbush', 'cbush.f06')
        op2name = os.path.join(model_path, 'cbush', 'cbush.op2')

        f06name2 = os.path.join(model_path, 'cbush', 'cbush.test_f06.f06')
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)

        assert len(f06.displacements) == 1, len(f06.displacements)  # 1 is correct
        assert len(f06.spc_forces) == 1, len(f06.spc_forces)        # 1 is correct

        assert len(f06.cbush_strain) == 0, len(f06.cbush_strain)  # 1 is correct
        assert len(f06.cbush_stress) == 0, len(f06.cbush_stress)  # 1 is correct

    def test_solid_shell_bar_1(self):
        bdfname = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        f06name = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.f06')
        op2name = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.op2')

        f06name2 = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.test_f06.f06')
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)

        assert len(f06.displacements) == 1, len(f06.displacements)
        assert len(f06.spc_forces) == 1, len(f06.spc_forces)
        assert len(f06.load_vectors) == 1, len(f06.load_vectors)

        assert len(f06.grid_point_forces) == 1, len(f06.grid_point_forces)

        assert len(f06.crod_force) == 1, len(f06.crod_force)
        assert len(f06.crod_strain) == 1, len(f06.crod_strain)
        assert len(f06.crod_stress) == 1, len(f06.crod_stress)

        assert len(f06.cquad4_force) == 1, len(f06.cquad4_force)  # 1 is correct
        assert len(f06.ctria3_force) == 0, len(f06.ctria3_force)  # 1 is correct

        assert len(f06.ctria3_stress) == 1, len(f06.ctria3_stress)
        assert len(f06.ctria3_strain) == 1, len(f06.ctria3_strain)
        assert len(f06.ctria6_stress) == 0, len(f06.ctria6_stress)
        assert len(f06.ctria6_strain) == 0, len(f06.ctria6_strain)

        assert len(f06.cquad4_stress) == 1, len(f06.cquad4_stress)
        assert len(f06.cquad4_strain) == 1, len(f06.cquad4_strain)
        assert len(f06.cquad8_stress) == 0, len(f06.cquad8_stress)
        assert len(f06.cquad8_strain) == 0, len(f06.cquad8_strain)

        assert len(f06.cbeam_force) == 0, len(f06.cbeam_force)    # 1 is correct
        assert len(f06.cbeam_strain) == 0, len(f06.cbeam_strain)  # 1 is correct

        assert len(f06.cbar_strain) == 1, len(f06.cbar_strain)
        assert len(f06.cbar_stress) == 1, len(f06.cbar_stress)

        assert len(f06.ctetra_strain) == 1, len(f06.ctetra_strain)
        assert len(f06.ctetra_stress) == 1, len(f06.ctetra_stress)
        assert len(f06.cpenta_strain) == 1, len(f06.cpenta_strain)
        assert len(f06.cpenta_stress) == 1, len(f06.cpenta_stress)
        assert len(f06.chexa_strain) == 1, len(f06.chexa_strain)
        assert len(f06.chexa_stress) == 1, len(f06.chexa_stress)

    #@unittest.expectedFailure  # fails
    def test_solid_shell_bar_2(self):
        bdfname = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.bdf')
        f06name = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.f06')
        op2name = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.op2')

        f06name2 = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.test_f06.f06')
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)

        assert op2.Title == 'MSC.NASTRAN JOB', '%r' % op2.Title
        assert f06.Title == 'MSC.NASTRAN JOB', '%r' % f06.Title

        subtitle_label = f06.iSubcaseNameMap[1]
        #assert subtitle_label[0] == 'DEFAULT', subtitle_label
        #assert subtitle_label[1] == 'SUBCASE 1', subtitle_label

        subtitle_label = op2.iSubcaseNameMap[1]
        assert subtitle_label[0] == 'DEFAULT', subtitle_label
        assert subtitle_label[1] == 2, subtitle_label  # 2=modal
        assert subtitle_label[2] == 'SUBCASE 1', subtitle_label

        assert len(f06.displacements) == 0, len(f06.displacements)  # 0 is correct
        assert len(f06.spc_forces) == 1, len(f06.spc_forces)        # 1 is correct
        assert len(f06.load_vectors) == 0, len(f06.load_vectors)    # 0 is correct
        assert len(f06.grid_point_forces) == 0, len(f06.grid_point_forces)  # 1 is correct

        assert len(f06.crod_force) == 1, len(f06.crod_force)    # 1 is correct
        assert len(f06.crod_strain) == 1, len(f06.crod_strain)  # 1 is correct
        assert len(f06.crod_stress) == 1, len(f06.crod_stress)  # 1 is correct

        assert len(f06.cquad4_force) == 1, len(f06.cquad4_force)  # 1 is correct
        assert len(f06.ctria3_force) == 0, len(f06.ctria3_force)  # 1 is correct


        assert len(f06.ctria3_stress) == 1, len(f06.ctria3_stress)  # 1 is correct
        assert len(f06.ctria3_strain) == 1, len(f06.ctria3_strain)  # 1 is correct
        assert len(f06.ctria6_stress) == 0, len(f06.ctria6_stress)  # 0 is correct
        assert len(f06.ctria6_strain) == 0, len(f06.ctria6_strain)  # 0 is correct

        assert len(f06.cquad4_stress) == 1, len(f06.cquad4_stress)  # 1 is correct
        assert len(f06.cquad4_strain) == 1, len(f06.cquad4_strain)  # 1 is correct
        assert len(f06.cquad8_stress) == 0, len(f06.cquad8_stress)  # 0 is correct
        assert len(f06.cquad8_strain) == 0, len(f06.cquad8_strain)  # 0 is correct

        assert len(f06.cbeam_force) == 0, len(f06.cbeam_force)    # 1 is correct
        assert len(f06.cbeam_strain) == 0, len(f06.cbeam_strain)  # 1 is correct

        assert len(f06.cbar_strain) == 1, len(f06.cbar_strain)  # 1 is correct
        assert len(f06.cbar_stress) == 1, len(f06.cbar_stress)  # 1 is correct

        assert len(f06.ctetra_strain) == 1, len(f06.ctetra_strain)  # 1 is correct
        assert len(f06.ctetra_stress) == 1, len(f06.ctetra_stress)  # 1 is correct
        assert len(f06.cpenta_strain) == 1, len(f06.cpenta_strain)  # 1 is correct
        assert len(f06.cpenta_stress) == 1, len(f06.cpenta_stress)  # 1 is correct
        assert len(f06.chexa_strain) == 1, len(f06.chexa_strain)    # 1 is correct
        assert len(f06.chexa_stress) == 1, len(f06.chexa_stress)    # 1 is correct

    def test_failure_index(self):
        bdfname = None
        f06name1 = os.path.join(test_path, 'failure_index_test.f06')
        f06name2 = os.path.join(test_path, 'failure_index_test.test_f06.f06')
        op2name = None
        f06 = run_model(bdfname, op2name, f06name1, f06_has_weight=False)

        #      FREQUENCY =  2.100000E-01
        #                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R
        #                                                          (REAL/IMAGINARY)
        #      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
        #0           21      G      9.121415E-18   5.869901E-15  -1.456074E+02   0.0            2.171751E+02   2.631261E-13
        #                           4.977449E-19   3.364763E-15  -6.035482E+00   0.0            2.480264E+01   2.345500E-14
        isubcase = 1
        assert len(f06.cquad4_force) == 2, len(f06.cquad4_force)  # 2 is correct
        assert len(f06.ctria3_force) == 0, len(f06.ctria3_force)  # 0 is correct

        assert len(f06.displacements) == 2, len(f06.displacements)
        assert len(f06.spc_forces) == 2, len(f06.spc_forces)
        assert len(f06.grid_point_forces) == 2, len(f06.grid_point_forces)
        assert len(f06.ctria3_composite_stress) == 0, len(f06.ctria3_composite_stress)
        assert len(f06.cquad4_composite_stress) == 2, len(f06.cquad4_composite_stress)

        assert len(f06.cquad4_force) == 2, len(f06.cquad4_force)  # 2 is correct
        assert len(f06.ctria3_force) == 0, len(f06.ctria3_force)  # 0 is correct

        #disp = f06.displacements[isubcase]
        #frequency = .21
        #T3 = disp.translations[frequency][21][2]
        #self.assertEqual(T3, -1.456074E+02 + -6.035482E+00j)  # T3

    def test_plate_openmdao(self):
        bdfname = os.path.join(model_path, 'plate', 'plate_openmdao.bdf')
        f06name = os.path.join(model_path, 'plate', 'plate.f06')
        op2name = os.path.join(model_path, 'plate', 'plate.op2')
        dynamic_vars = {'t' : 42.}
        bdf, op2, f06 = run_model(bdfname, op2name, f06name, dynamic_vars=dynamic_vars, f06_has_weight=False)
        self.assertEqual(bdf.properties[1].t, 42., 't=%s' % bdf.properties[1].t)

        assert len(f06.displacements) == 1, len(f06.displacements)
        assert len(f06.spc_forces) == 1, len(f06.spc_forces)

        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.spc_forces) == 1, len(op2.spc_forces)

        assert len(op2.ctria3_stress) == 0, len(op2.ctria3_stress)
        assert len(op2.ctria3_strain) == 0, len(op2.ctria3_strain)
        assert len(op2.ctria6_stress) == 0, len(op2.ctria6_stress)
        assert len(op2.ctria6_strain) == 0, len(op2.ctria6_strain)

        assert len(op2.cquad4_stress) == 1, len(op2.cquad4_stress)
        assert len(op2.cquad4_strain) == 0, len(op2.cquad4_strain)
        assert len(op2.cquad8_stress) == 0, len(op2.cquad8_stress)
        assert len(op2.cquad8_strain) == 0, len(op2.cquad8_strain)

        assert len(f06.ctria3_stress) == 0, len(f06.ctria3_stress)
        assert len(f06.ctria3_strain) == 0, len(f06.ctria3_strain)
        assert len(f06.ctria6_stress) == 0, len(f06.ctria6_stress)
        assert len(f06.ctria6_strain) == 0, len(f06.ctria6_strain)

        assert len(f06.cquad4_stress) == 1, len(f06.cquad4_stress)
        assert len(f06.cquad4_strain) == 0, len(f06.cquad4_strain)
        assert len(f06.cquad8_stress) == 0, len(f06.cquad8_stress)
        assert len(f06.cquad8_strain) == 0, len(f06.cquad8_strain)

        dynamic_vars = {'t' : 42}
        with self.assertRaises(SyntaxError):
            bdf2 = run_model(bdfname, dynamic_vars=dynamic_vars)

        dynamic_vars = {'t' : 'asdddddddf'}
        with self.assertRaises(SyntaxError):
            bdf3 = run_model(bdfname, dynamic_vars=dynamic_vars)

    def test_complex_displacement(self):
        bdfname = None
        f06name1 = os.path.join(test_path, 'complex_displacement.f06')
        f06name2 = os.path.join(test_path, 'complex_displacement.test_f06.f06')
        op2name = None
        f06 = run_model(bdfname, op2name, f06name1, f06_has_weight=False)

        #      FREQUENCY =  2.100000E-01
        #                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R
        #                                                          (REAL/IMAGINARY)
        #      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
        #0           21      G      9.121415E-18   5.869901E-15  -1.456074E+02   0.0            2.171751E+02   2.631261E-13
        #                           4.977449E-19   3.364763E-15  -6.035482E+00   0.0            2.480264E+01   2.345500E-14
        isubcase = 1
        disp = f06.displacements[isubcase]
        frequency = .21
        T3 = disp.translations[frequency][21][2]
        self.assertTrue(allclose(T3.real, -1.456074E+02))
        self.assertTrue(allclose(T3.imag, -6.035482E+00))

        #f06.write_f06(f06name2, quiet=True)
        #os.remove(f06name2)

    def test_eigenvectors1(self):
        bdfname = None
        f06name = os.path.join(test_path, 'test_no_rotations.f06')
        op2name = None
        f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)
        subcase = 1
        eig = 1
        nid = 2
        self.assertTrue(array_equal(f06.eigenvectors[subcase].translations[eig][2], array([4., 5., 6.])),
                            msg=str(f06.eigenvectors[subcase].translations[eig][2]))
        self.assertTrue(array_equal(f06.eigenvectors[subcase].translations[eig][4], array([10., 0., 0.])),
                            msg=str(f06.eigenvectors[subcase].translations[eig][4]))


        self.assertTrue(array_equal(f06.eigenvectors[subcase].rotations[eig][2], array([0., 0., 0.])))
        self.assertTrue(array_equal(f06.eigenvectors[subcase].rotations[eig][4], array([0., 0., 0.])))

    def test_eigenvectors2(self):
        bdfname = None
        f06name = os.path.join(test_path, 'test_with_rotations.f06')
        op2name = None

        subcase = 1
        eig = 1
        f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)
        self.assertTrue(array_equal(f06.eigenvectors[subcase].translations[eig][2], array([4., 5., 6.])),
                            msg=str(f06.eigenvectors[subcase].translations[eig][2]))

        self.assertTrue(array_equal(f06.eigenvectors[subcase].translations[eig][4], array([10., 11., 12.])),
                            msg=str(f06.eigenvectors[subcase].translations[eig][4]))

        self.assertTrue(array_equal(f06.eigenvectors[subcase].rotations[eig][2], array([1., 1., 1.])))
        self.assertTrue(array_equal(f06.eigenvectors[subcase].rotations[eig][4], array([1., 1., 1.])))

    def test_plate_vonmises(self):
        bdfname = os.path.join(model_path, 'plate', 'plate.bdf')
        f06name = os.path.join(model_path, 'plate', 'plate.f06')
        op2name = os.path.join(model_path, 'plate', 'plate.op2')

        bdf, op2, f06 = run_model(bdfname, op2name, f06name, f06_has_weight=False)
        self.assertEqual(bdf.properties[1].t,  0.3, 't=%s' % bdf.properties[1].t)

        self.assertEqual(len(bdf.nodes), 36, bdf.nodes)
        self.assertEqual(len(bdf.elements), 25, bdf.elements)
        self.assertEqual(len(bdf.properties), 1, bdf.properties)
        self.assertEqual(len(bdf.materials), 1, bdf.materials)
        self.assertEqual(len(bdf.loads), 9, 'nloads=%s\n%s' % (len(bdf.loads), bdf.loads))  # FORCE, LOAD
        self.assertEqual(len(bdf.params), 2, bdf.params)
        self.assertEqual(bdf.sol, 101, bdf.sol)

        cen = 0 # 'CEN/4'
        for (loadcase, stress) in iteritems(f06.cquad4_stress):
            # stress is a PlateStressObject
            if stress.is_von_mises():
                #print("%3s %3s %6s %8s" % ('eID', 'NID', 'iLayer', 'VM_Stress'))
                #vonMises = 'VON MISES'
                for eid,ovm in sorted(iteritems(stress.ovmShear)):
                    for nid, ovmi in sorted(iteritems(ovm)):
                        for ilayer, ovmii in enumerate(ovmi):
                            print("%8s %8s %6s %8s" % (eid, nid, ilayer, ovmii))
            else:
                #print("%3s %3s %6s %8s" % ('eID', 'NID', 'iLayer', 'MaxShear'))
                #vonMises = 'MAX SHEAR'
                for eid,ovm in sorted(iteritems(stress.ovmShear)):
                    for nid, ovmi in sorted(iteritems(ovm)):
                        ovmi = ovm[nid]
                        for ilayer, ovmii in enumerate(ovmi):
                            o1 = stress.oxx[eid][nid][ilayer]
                            o2 = stress.oyy[eid][nid][ilayer]
                            t12 = stress.txy[eid][nid][ilayer]
                            ovmii2 = sqrt(o1**2 - o1*o2 + o2**2 + 3*t12**2)
                            self.assertAlmostEqual(ovmii, ovmii2, places=3)
                            #print("%3s %3s %6s %8s" % (eid, nid, ilayer, ovmii2))

        assert cen == 0, cen
        assert len(op2.cquad4_stress.keys()) == 1, op2.cquad4_stress.keys()
        assert len(f06.cquad4_stress.keys()) == 1, f06.cquad4_stress.keys()
        assert list(op2.cquad4_stress.keys()) == [1], list(op2.cquad4_stress.keys())
        assert list(f06.cquad4_stress.keys()) == [1], list(f06.cquad4_stress.keys())
        self.assertEqual(op2.cquad4_stress[1].ovmShear[25][cen][0], 276.8023376464844)
        self.assertEqual(f06.cquad4_stress[1].ovmShear[25][cen][0], 276.8023)
        self.assertTrue(op2.cquad4_stress[1].is_von_mises())
        #self.assertTrue(f06.cquad4_stress[1].is_von_mises())

        self.assertTrue(op2.cquad4_stress[1].is_fiber_distance())
        #self.assertTrue(f06.cquad4_stress[1].is_fiber_distance())
        #f06.print_stats()
        #op2.print_stats()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()