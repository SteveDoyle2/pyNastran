import os
import unittest

import numpy as np
from numpy import allclose
from cpylog import SimpleLogger

#import pyNastran
#from pyNastran.bdf.bdf import BDF

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')
from pyNastran.bdf.cards.elements.mass import CONM2

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, CaseControlDeck, PARAM
from pyNastran.bdf.mesh_utils.convert import convert, get_scale_factors, scale_by_terms
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh

pkg_path = pyNastran.__path__[0]

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)

class TestConvert(unittest.TestCase):
    """various BDF conversion tests"""

    def test_convert_bar(self):
        """converts a bar model"""
        log = SimpleLogger(level='warning')
        model_path = os.path.join(pkg_path, '..', 'models', 'beam_modes')
        bdf_filename = os.path.join(model_path, 'beam_modes.dat')
        bdf_filename_out = os.path.join(model_path, 'beam_modes_temp.bdf')
        bdf_filename_out2 = os.path.join(model_path, 'beam_modes_converted.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)
        #card_lines = ['EIGRL', 42, None, None, 20]
        #model.add_card(card_lines, 'EIGRL')
        #model.case_control_deck = CaseControlDeck(lines)
        model.write_bdf(bdf_filename_out)
        #units_from = ['in', 'lbm', 's']
        units_from = ['mm', 'Mg', 's']
        units_to = ['m', 'kg', 's']

        convert(model, units_to, units=units_from)
        model.write_bdf(bdf_filename_out2)
        os.remove(bdf_filename_out)
        os.remove(bdf_filename_out2)

        terms = ['F', 'P', 'V']
        scales = [1.1, 0.9, 0.8]
        scale_by_terms(bdf_filename, terms, scales, bdf_filename_out=bdf_filename_out, log=log)
        os.remove(bdf_filename_out)

    def test_convert_isat(self):
        """converts a isat model"""
        log = SimpleLogger(level='error')
        model_path = os.path.join(pkg_path, '..', 'models', 'iSat')
        bdf_filename = os.path.join(model_path, 'ISat_Dploy_Sm.dat')
        bdf_filename_out = os.path.join(model_path, 'isat.bdf')
        bdf_filename_out2 = os.path.join(model_path, 'isat_converted.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)
        #card_lines = ['EIGRL', 42, None, None, 20]
        #model.add_card(card_lines, 'EIGRL')
        #model.case_control_deck = CaseControlDeck(lines)
        model.write_bdf(bdf_filename_out)
        #units_from = ['in', 'lbm', 's']
        units_from = ['mm', 'Mg', 's']
        units_to = ['m', 'kg', 's']

        convert(model, units_to, units=units_from)
        model.write_bdf(bdf_filename_out2)
        os.remove(bdf_filename_out)
        os.remove(bdf_filename_out2)

    def test_convert_bwb(self):
        """converts a bwb model"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero.bdf')
        bdf_filename_out = os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_modes.bdf')
        bdf_filename_out2 = os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_modes_converted.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)
        model.sol = 103

        lines = [
            'ECHO = NONE',
            'SUBCASE 1',
            '    DISPLACEMENT(PLOT) = ALL',
            '    MPC = 1',
            '    SPC = 100',
            '    SUPORT1 = 1',
            '    METHOD = 42',
        ]
        card_lines = ['EIGRL', 42, None, None, 20]
        model.add_card(card_lines, 'EIGRL')
        model.case_control_deck = CaseControlDeck(lines, log=log)
        model.write_bdf(bdf_filename_out)
        units_from = ['in', 'lbm', 's']
        #units_from = ['mm', 'Mg', 's']
        units_to = ['m', 'kg', 's']

        convert(model, units_to, units=units_from)
        model.write_bdf(bdf_filename_out2)
        os.remove(bdf_filename_out)
        os.remove(bdf_filename_out2)

    def test_convert_sine(self):
        """converts a sine model"""
        log = SimpleLogger(level='error')
        model_path = os.path.join(pkg_path, '..', 'models', 'freq_sine')
        bdf_filename = os.path.join(model_path, 'good_sine.dat')
        bdf_filename_out = os.path.join(model_path, 'sine_modes.bdf')
        bdf_filename_out2 = os.path.join(model_path, 'sine_converted.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)
        model.sol = 103

        lines = [
            'ECHO = NONE',
            'SUBCASE 1',
            '    DISPLACEMENT(PLOT) = ALL',
            #'$    SPC = 100',
            '    METHOD = 42',
        ]
        card_lines = ['EIGRL', 42, None, None, 20]
        model.add_card(card_lines, 'EIGRL')
        model.case_control_deck = CaseControlDeck(lines, log=log)
        model.params['GRDPNT'] = PARAM('GRDPNT', 0)
        #del model.params['WTMASS']
        model.write_bdf(bdf_filename_out)
        #units_from = ['in', 'lbm', 's']

        units_from = ['mm', 'Mg', 's']
        units_to = ['m', 'kg', 's']

        convert(model, units_to, units=units_from)
        model.write_bdf(bdf_filename_out2)
        os.remove(bdf_filename_out)
        os.remove(bdf_filename_out2)

    def test_convert_null(self):
        """null conversions to verify any conversion is undone consistently"""
        log = SimpleLogger(level='error')
        pairs = [
            ['m', 'kg', 's'],
            ['mm', 'kg', 's'],
            ['mm', 'g', 's'],
            ['mm', 'Mg', 's'],
            ['ft', 'lbm', 's'],
            ['in', 'lbm', 's'],
            ['cm', 'kg', 's'],
        ]
        for pair in pairs:
            xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
                pair, pair, log)
            assert xyz_scale == 1., xyz_scale
            assert mass_scale == 1., mass_scale
            assert time_scale == 1., time_scale
            assert weight_scale == 1., weight_scale
            assert gravity_scale == 1., gravity_scale

        invalid_pairs = [
            ['cm', 'ton', 's'],
        ]
        pair0 = ['m', 'kg', 's']
        for pair in invalid_pairs:
            with self.assertRaises(NotImplementedError):
                xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
                    pair, pair0, log)
            with self.assertRaises(NotImplementedError):
                xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
                    pair0, pair, log)

    def test_convert_units(self):
        """tests various conversions"""
        log = SimpleLogger(level='error')
        # from -> to
        xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
            ['in', 'lbm', 's'], ['ft', 'lbm', 's'], log)
        assert xyz_scale == 1./12.
        assert mass_scale == 1.
        assert time_scale == 1.
        assert weight_scale == 1., weight_scale
        assert gravity_scale == 1./12., gravity_scale
        wtmass = 1. / (32.174 * 12.)
        wtmass_expected = 1. / (32.174)
        assert allclose(wtmass/gravity_scale, wtmass_expected), 'wtmass=%s wtmass_expected=%s' % (wtmass, wtmass_expected)

        xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
            ['mm', 'Mg', 's'], ['m', 'kg', 's'], log)
        assert xyz_scale == 1./1000.
        assert mass_scale == 1000.
        assert time_scale == 1.
        assert weight_scale == 1., weight_scale
        assert gravity_scale == 1.

        xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
            ['ft', 'lbm', 's'], ['m', 'kg', 's'], log)
        assert xyz_scale == 0.3048
        assert mass_scale == 0.45359237, mass_scale
        assert time_scale == 1.
        assert allclose(weight_scale, 4.4482216526), weight_scale
        assert allclose(gravity_scale, 1/32.174), 'gravity_scale=%s 1/expected=%s' % (gravity_scale, 1/(32.2))
        wtmass = 1. / (32.174)
        wtmass_expected = 1.
        assert allclose(wtmass/gravity_scale, wtmass_expected), 'wtmass=%s wtmass_expected=%s' % (wtmass/gravity_scale, wtmass_expected)

        # both are consistent systems, so wtmass=1.0
        xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
            ['in', 'slinch', 's'], ['in', 'slug', 's'], log)
        assert xyz_scale == 1.
        assert mass_scale == 12., mass_scale
        assert time_scale == 1.
        assert np.allclose(weight_scale, 1.), weight_scale
        assert gravity_scale == 1., gravity_scale
        wtmass = 1.
        wtmass_expected = 1.
        assert allclose(wtmass/gravity_scale, wtmass_expected), 'wtmass=%s wtmass_expected=%s' % (wtmass, wtmass_expected)

        # both are consistent systems, so wtmass=1.0
        xyz_scale, mass_scale, time_scale, weight_scale, gravity_scale = get_scale_factors(
            ['ft', 'lbm', 's'], ['ft', 'slug', 's'], log)
        assert xyz_scale == 1.
        assert np.allclose(mass_scale, 1. / 32.174), 'mass_scale=%s expected=%s' % (mass_scale, 1/32.174)
        assert time_scale == 1.
        assert np.allclose(weight_scale, 1.), 'weight_scale=%s expected=%s' % (weight_scale, 1.)
        assert gravity_scale == 1/32.174, 'gravity_scale=%s expected=%s' % (gravity_scale, 1./32.174)
        wtmass = 1. / 32.174
        wtmass_expected = 1.
        assert allclose(wtmass/gravity_scale, wtmass_expected), 'wtmass=%s wtmass_expected=%s' % (wtmass, wtmass_expected)

    def test_convert_01(self):
        """converts the CONM2s units"""
        log = SimpleLogger(level='error')
        model = BDF(log=log)
        eid = 1000
        nid = 100
        cid = 0
        mass = 247200. # kg
        X = [30.16, 0., 3.55] # m
        I11 = 1.39e7 # kg-m^2
        I22 = 3.66e7
        I33 = 4.99e7
        I13 = I12 = I23 = 0.
        I = I11, I12, I22, I13, I23, I33
        elem = CONM2(eid, nid, mass, cid=cid, X=X, I=I, comment='')
        model.masses[eid] = elem

        units_to = ['in', 'lbm', 's']
        units_from = ['m', 'kg', 's']
        convert(model, units_to, units=units_from)
        #print(model.masses[eid].write_card_16())

    def test_convert_02(self):
        """converts a full model units"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero.bdf'))
        bdf_filename_out = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero.out'))

        model = read_bdf(bdf_filename, log=log)
        units_to = ['m', 'kg', 's']
        units_from = ['in', 'lbm', 's']
        #units_to = units_from
        convert(model, units_to, units_from)
        model.write_bdf(bdf_filename_out)

        caero_bdf_filename = 'caero.bdf'
        export_caero_mesh(model, caero_bdf_filename=caero_bdf_filename)
        os.remove(bdf_filename_out)
        os.remove(caero_bdf_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
