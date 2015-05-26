from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
from numpy import allclose, array
from numpy.linalg import norm

import pyNastran
from pyNastran.bdf.cards.baseCard import collapse_thru_by
from pyNastran.bdf.bdf import BDF

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'bdf', 'test')
#print("testPath = %s" % testPath)
from pyNastran.bdf.test.test_bdf import run_bdf, run_all_files_in_folder


class Tester(unittest.TestCase):

    def run_bdf(self, folder, bdf_filename, xref=False, cid=None,
                meshForm='combined', debug=False, dynamic_vars={}):
        cid = 0
        #xref = False
        return run_bdf(folder, bdf_filename, xref=xref, cid=cid, isFolder=True,
                       meshForm=meshForm, dynamic_vars=dynamic_vars, debug=debug)

    def run_all_files_in_folder(self, folder, xref=False, cid=None, debug=False):
        run_all_files_in_folder(folder, xref=xref, cid=cid, debug=debug)


class TestBDF(Tester):
    def test_bdf_01(self):
        bdf_filename = os.path.join('solid_bending', 'solid_bending.bdf')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        self.run_bdf(folder, bdf_filename)
        fem1, fem2, diffCards = self.run_bdf(folder, bdf_filename, xref=True)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 2, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 72, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 186, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        mass, cg, I = fem1.mass_properties()

        assert allclose(mass, 6.0), 'mass = %s' % mass
        cg_exact = array([0.5, 1., 1.5])
        for i, (cgi, cgie) in enumerate(zip(cg, cg_exact)):
            assert allclose(cgi, cgie), 'i=%s cg=%s' % (i, str(cg))

        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')
        mass, cg, I = fem1.mass_properties(reference_point='cg')

    def _compare_mass_cg_I(self, fem1, reference_point=None, sym_axis=None):
        num_cpus = 4
        mass1, cg1, I1 = fem1.mass_properties(reference_point=reference_point, sym_axis=sym_axis)
        mass2, cg2, I2 = fem1.mass_properties(reference_point=reference_point, sym_axis=sym_axis, num_cpus=num_cpus)

        assert allclose(mass1, mass2), 'mass1_sp=%s mass2_mp=%s' % (mass1, mass2)
        assert allclose(norm((cg1 - cg2)**2), 0.0), 'cg1-cg2=%s' % (cg1 - cg2)
        assert allclose(norm((I1  -  I2)**2), 0.0), 'I1-I2=%s' % (I1 - I2)

    def test_bdf_02(self):
        bdf_filename = os.path.join('plate_py', 'plate_py.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        self.run_bdf(folder, bdf_filename)
        fem1, fem2, diffCards = self.run_bdf(folder, bdf_filename, xref=True)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2

        for fem in [fem1, fem2]:
            assert len(fem.coords) == 3, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.nodes) == 231, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 200, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')


    def test_bdf_03(self):
        bdf_filename = os.path.join('cbush', 'cbush.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        fem1, fem2, diffCards = self.run_bdf(folder, bdf_filename)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 2, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 0, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 1, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue

        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')

        self.run_bdf(folder, bdf_filename, xref=True)

    def test_bdf_04(self):
        bdf_filename = os.path.join('beam_modes', 'beam_modes.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        fem1, fem2, diffCards = self.run_bdf(folder, bdf_filename)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 12, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 10, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.masses) == 1, 'len(masses) = %i' % len(fem.elements)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 3, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue
            assert len(fem.properties_mass) == 0, 'len(properties_mass) = %i' % len(fem.properties_mass)
        self._compare_mass_cg_I(fem1)
        #self._compare_mass_cg_I(fem1, reference_point=u'cg')

        #self.run_bdf(folder, bdf_filename, xref=True) # PBEAML is not supported

    def test_bdf_05(self):
        bdf_filename = 'testA.bdf'
        folder = os.path.abspath(os.path.join(pkg_path, 'bdf', 'test', 'unit'))
        (fem1, fem2, diffCards) = self.run_bdf(folder, bdf_filename)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2
        #self.run_bdf(folder, bdf_filename, xref=True) # PBEAML is not supported

    def test_bdf_06(self):
        bdf_filename = os.path.join('bar3truss', 'vared_bar3.bdf')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))

        dynamic_vars = {
            'bar1_a': 1.0,
            'bar2_a': 1.0,
            'bar3_a': 1.0,
            'loadx': 1.0,
            'loady': 1.0,
            'loadmag': 10000.,
            'youngs' : 1e7,
            'rho': 0.01,
        }
        fem1, fem2, diffCards = self.run_bdf(folder, bdf_filename, dynamic_vars=dynamic_vars)
        diffCards2 = list(set(diffCards))
        diffCards2.sort()
        assert len(diffCards2) == 0, diffCards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 4, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 4, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 3, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 3, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue

        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')
        self._compare_mass_cg_I(fem1, reference_point='cg')


class BaseCard_Test(Tester):
    def test_base_card_01_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        msg = 'expected=%s actual=%s' % (expected, collapse_thru_by(data))
        self.assertEqual(collapse_thru_by(data), expected, msg)

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 7, 9, 11, 12, 14, 16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2]
        expected = [1, 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 5, 7, 9, 11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4]
        expected = [1, 'THRU', 4]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3]
        expected = [1, 2, 3]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
