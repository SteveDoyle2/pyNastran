"""various atmosphere tests"""
import os
import unittest
import numpy as np

from pyNastran.utils.atmosphere import (
    atm_density, atm_dynamic_pressure, atm_temperature,
    atm_pressure, atm_velocity, atm_mach, atm_equivalent_airspeed,
    atm_dynamic_viscosity_mu, atm_kinematic_viscosity_nu,
    atm_calibrated_airspeed, cas_to_mach,
    get_alt_for_density, get_alt_for_pressure,
    get_alt_for_mach_eas,
    get_alt_for_q_with_constant_mach,
    get_alt_for_eas_with_constant_mach,
    atm_unit_reynolds_number,
    make_flfacts_alt_sweep_constant_mach,
    make_flfacts_mach_sweep_constant_alt,
    make_flfacts_eas_sweep_constant_alt,
    make_flfacts_eas_sweep_constant_mach,
    sutherland_viscoscity,
    _reynolds_factor,
    create_atmosphere_table,
)

from pyNastran.utils.convert import (
    convert_length, convert_area,
    convert_density, convert_mass,
    convert_velocity, convert_force, convert_pressure,
    convert_temperature,
    _length_factor, _area_factor,
    _density_factor, _velocity_factor, _force_factor,
    _temperature_factor,
)
class TestAtmConvert(unittest.TestCase):
    """various unit conversion tests"""

    def test_length(self):
        """length checks"""
        assert np.allclose(convert_length(1., 'ft', 'ft'), 1.)
        with self.assertRaises(RuntimeError):
            convert_length(1., 'ft', 'bad')
        with self.assertRaises(RuntimeError):
            convert_length(1., 'bad', 'ft')

        assert np.allclose(_length_factor('ft', 'ft'), 1.)
        assert np.allclose(_length_factor('in', 'in'), 1.)
        assert np.allclose(_length_factor('m', 'm'), 1.)
        assert np.allclose(_length_factor('cm', 'cm'), 1.)
        assert np.allclose(_length_factor('mm', 'mm'), 1.)

        assert np.allclose(_length_factor('ft', 'in'), 12.)
        assert np.allclose(_length_factor('in', 'ft'), 1 / 12.)

        assert np.allclose(_length_factor('ft', 'ft'), 1,)
        assert np.allclose(_length_factor('ft', 'm'), 1 / 3.28084)

        assert np.allclose(_length_factor('m', 'ft'), 3.28084)
        assert np.allclose(_length_factor('m', 'cm'), 100.)
        assert np.allclose(_length_factor('m', 'mm'), 1000.)

        assert np.allclose(_length_factor('cm', 'mm'), 10.)
        assert np.allclose(_length_factor('cm', 'm'), 1 / 100.)

        assert np.allclose(_length_factor('cm', 'in'), 1/2.54), _length_factor('cm', 'in')
        assert np.allclose(_length_factor('in', 'cm'), 2.54), _length_factor('in', 'cm')

        assert np.allclose(_length_factor('mm', 'm'), 1 / 1000.)
        assert np.allclose(_length_factor('mm', 'cm'), 1 / 10.)

    def test_area(self):
        """area checks"""
        assert np.allclose(convert_area(1., 'ft^2', 'ft^2'), 1.)
        with self.assertRaises(RuntimeError):
            convert_area(1., 'ft^2', 'bad^2')
        with self.assertRaises(RuntimeError):
            convert_area(1., 'bad^2', 'ft^2')

        assert np.allclose(_area_factor('ft^2', 'ft^2'), 1.)
        assert np.allclose(_area_factor('in^2', 'in^2'), 1.)
        assert np.allclose(_area_factor('m^2', 'm^2'), 1.)
        assert np.allclose(_area_factor('cm^2', 'cm^2'), 1.)
        assert np.allclose(_area_factor('mm^2', 'mm^2'), 1.)

        assert np.allclose(_area_factor('ft^2', 'in^2'), 12.**2)
        assert np.allclose(_area_factor('in^2', 'ft^2'), 1 / 12.**2)

        assert np.allclose(_area_factor('ft^2', 'ft^2'), 1.)
        assert np.allclose(_area_factor('ft^2', 'm^2'), 1 / 3.28084**2)

        assert np.allclose(_area_factor('m^2', 'ft^2'), 3.28084**2)
        assert np.allclose(_area_factor('m^2', 'cm^2'), 100.**2)
        assert np.allclose(_area_factor('m^2', 'mm^2'), 1000.**2)

        assert np.allclose(_area_factor('cm^2', 'mm^2'), 10.**2)
        assert np.allclose(_area_factor('cm^2', 'm^2'), 1 / 100.**2)

        assert np.allclose(_area_factor('cm^2', 'in^2'), 1/2.54**2), _area_factor('cm^2', 'in^2')
        assert np.allclose(_area_factor('in^2', 'cm^2'), 2.54**2), _area_factor('in^2', 'cm^2')

        assert np.allclose(_area_factor('mm^2', 'm^2'), 1 / 1000.**2)
        assert np.allclose(_area_factor('mm^2', 'cm^2'), 1 / 10.**2)

    def test_mass(self):
        """mass checks"""
        assert np.allclose(convert_mass(1., 'slinch', 'slinch'), 1.)
        with self.assertRaises(RuntimeError):
            convert_mass(1., 'slinch', 'bad')
        with self.assertRaises(RuntimeError):
            convert_mass(1., 'bad', 'slinch')

        assert np.allclose(convert_mass(1., 'slinch', 'slug'), 12.)
        assert np.allclose(convert_mass(1., 'slug', 'slinch'), 1 / 12.)

        assert np.allclose(convert_mass(1., 'slug', 'kg'), 14.5939)
        assert np.allclose(convert_mass(1., 'kg', 'slug'), 1 / 14.5939)

        #assert np.allclose(convert_mass(1., 'Mg', 'kg'), 1000)
        #assert np.allclose(convert_mass(1., 'kg', 'Mg'), 1 / 1000)

    def test_force(self):
        """force checks"""
        assert np.allclose(convert_force(1., 'lbf', 'lbf'), 1.)
        with self.assertRaises(RuntimeError):
            convert_force(1., 'lbf', 'bad')
        with self.assertRaises(RuntimeError):
            convert_force(1., 'bad', 'lbf')

        assert np.allclose(_force_factor('lbf', 'N'), 4.44822)
        assert np.allclose(_force_factor('N', 'lbf'), 1 / 4.44822)

        assert np.allclose(_force_factor('MN', 'N'), 1000000.)
        assert np.allclose(_force_factor('N', 'MN'), 1 / 1000000.)

        assert np.allclose(_force_factor('mN', 'N'), 0.001)
        assert np.allclose(_force_factor('N', 'mN'), 1 / 0.001)

        assert np.allclose(_force_factor('cN', 'N'), 0.01)
        assert np.allclose(_force_factor('N', 'cN'), 1 / 0.01)

    def test_pressure(self):
        """pressure checks"""
        assert np.allclose(convert_pressure(1., 'psf', 'psf'), 1.)
        with self.assertRaises(RuntimeError):
            convert_pressure(1., 'psf', 'bad')
        with self.assertRaises(RuntimeError):
            convert_pressure(1., 'bad', 'psf')

        assert np.allclose(convert_pressure(1., 'bar', 'Pa'), 100_000.)
        assert np.allclose(convert_pressure(1., 'Pa', 'bar'), 1 / 100_000.)

        assert np.allclose(convert_pressure(1., 'bar', 'kPa'), 100.)
        assert np.allclose(convert_pressure(1., 'kPa', 'bar'), 1 / 100.)

        #------------
        assert np.allclose(convert_pressure(1., 'psf', 'Pa'), 47.880208)
        assert np.allclose(convert_pressure(1., 'Pa', 'psf'), 1 / 47.880208)

        assert np.allclose(convert_pressure(1., 'kPa', 'Pa'), 1000.)
        assert np.allclose(convert_pressure(1., 'Pa', 'kPa'), 1 / 1000.)

        assert np.allclose(convert_pressure(1., 'MPa', 'Pa'), 1000000.)
        assert np.allclose(convert_pressure(1., 'Pa', 'MPa'), 1 / 1000000.)

        assert np.allclose(convert_pressure(1., 'psi', 'psf'), 144.)
        assert np.allclose(convert_pressure(1., 'psf', 'psi'), 1 / 144.)

    def test_velocity(self):
        """velocity checks"""
        assert np.allclose(convert_velocity(1., 'ft/s', 'ft/s'), 1.)
        with self.assertRaises(RuntimeError):
            convert_velocity(1., 'ft/s', 'bad')
        with self.assertRaises(RuntimeError):
            convert_velocity(1., 'bad', 'ft/s')

        assert np.allclose(_velocity_factor('ft/s', 'ft/s'), 1.)
        assert np.allclose(_velocity_factor('in/s', 'in/s'), 1.)
        assert np.allclose(_velocity_factor('m/s', 'm/s'), 1.)
        assert np.allclose(_velocity_factor('cm/s', 'cm/s'), 1.)
        assert np.allclose(_velocity_factor('mm/s', 'mm/s'), 1.)

        assert np.allclose(_velocity_factor('ft/s', 'in/s'), 12.)
        assert np.allclose(_velocity_factor('in/s', 'ft/s'), 1 / 12.)

        assert np.allclose(_velocity_factor('m/s', 'ft/s'), 3.28084)
        assert np.allclose(_velocity_factor('ft/s', 'm/s'), 1 / 3.28084)

        assert np.allclose(_velocity_factor('ft/s', 'm/s'), 1 / 3.28084)
        assert np.allclose(_velocity_factor('m/s', 'ft/s'), 3.28084)

        assert np.allclose(_velocity_factor('ft/s', 'cm/s'), 1 / 3.28084 * 100)
        assert np.allclose(_velocity_factor('ft/s', 'mm/s'), 1 / 3.28084 * 1000)

        assert np.allclose(_velocity_factor('cm/s', 'm/s'), 1 / 100.)
        assert np.allclose(_velocity_factor('m/s', 'cm/s'), 100.)

        assert np.allclose(_velocity_factor('m/s', 'mm/s'), 1000.)
        assert np.allclose(_velocity_factor('mm/s', 'm/s'), 1 / 1000.)

        assert np.allclose(_velocity_factor('mm/s', 'cm/s'), 1 / 10.)
        assert np.allclose(_velocity_factor('cm/s', 'mm/s'), 10.)

        assert np.allclose(_velocity_factor('knots', 'ft/s'), 1.68781)
        assert np.allclose(_velocity_factor('ft/s', 'knots'), 1 / 1.68781)

    def test_density(self):
        """density checks"""
        assert np.allclose(convert_density(1., 'slug/ft^3', 'slug/ft^3'), 1.)
        with self.assertRaises(RuntimeError):
            convert_density(1., 'slug/ft^3', 'bad')
        with self.assertRaises(RuntimeError):
            convert_density(1., 'bad', 'slug/ft^3')

        assert np.allclose(_density_factor('slinch/in^3', 'slinch/in^3'), 1)
        assert np.allclose(_density_factor('slug/ft^3', 'slug/ft^3'), 1)
        assert np.allclose(_density_factor('kg/m^3', 'kg/m^3'), 1)
        assert np.allclose(_density_factor('g/cm^3', 'g/cm^3'), 1)
        assert np.allclose(_density_factor('Mg/mm^3', 'Mg/mm^3'), 1)

        assert np.allclose(_density_factor('slinch/in^3', 'slug/ft^3'), 12**4)
        assert np.allclose(_density_factor('slug/ft^3', 'slinch/in^3'), 1 / 12**4)

        assert np.allclose(_density_factor('slug/ft^3', 'kg/m^3'), 515.379)
        assert np.allclose(_density_factor('kg/m^3', 'slug/ft^3'), 1 / 515.379)

        assert np.allclose(_density_factor('g/cm^3', 'kg/m^3'), 1000.), 'actual=%g expected=%g' % (convert_density(1., 'g/cm^3', 'kg/m^3'), 1000.)
        assert np.allclose(_density_factor('kg/m^3', 'g/cm^3'), 1 / 1000.), 'actual=%g expected=%g' % (convert_density(1., 'kg/m^3', 'g/cm^3'), 1 / 1000.)

        assert np.allclose(_density_factor('Mg/mm^3', 'kg/m^3'), 1000.**4), 'actual=%g expected=%g' % (convert_density(1., 'Mg/mm^3', 'kg/m^3'), 1000.**4)
        assert np.allclose(_density_factor('kg/m^3', 'Mg/mm^3'), 1/1000.**4), 'actual=%g expected=%g' % (convert_density(1., 'kg/m^3', 'Mg/mm^3'), 1/1000.**4)

        assert np.allclose(_density_factor('g/cm^3', 'slug/ft^3'), 1.94032)
        assert np.allclose(_density_factor('slug/ft^3', 'g/cm^3'), 1 / 1.94032), 'actual=%g expected=%g' % (convert_density(1., 'slug/ft^3', 'g/cm^3'), 1 / 1.94032)

    def test_temperature(self):
        """density checks"""
        assert np.allclose(convert_temperature(1., 'R', 'R'), 1.)
        with self.assertRaises(RuntimeError):
            convert_temperature(1., 'F', 'bad')
        with self.assertRaises(RuntimeError):
            convert_temperature(1., 'bad', 'C')

        assert np.allclose(_temperature_factor(1., 'R', 'R'), 1)
        assert np.allclose(_temperature_factor(1., 'F', 'F'), 1)
        assert np.allclose(_temperature_factor(1., 'C', 'C'), 1)
        assert np.allclose(_temperature_factor(1., 'K', 'K'), 1)

        assert np.allclose(_temperature_factor(32., 'F', 'C'), 0.)
        assert np.allclose(_temperature_factor(0., 'C', 'F'), 32.)

        assert np.allclose(_temperature_factor(104., 'F', 'C'), 40.)
        assert np.allclose(_temperature_factor(40., 'C', 'F'), 104.)

        assert np.allclose(_temperature_factor(40., 'F', 'R'), 499.67)
        assert np.allclose(_temperature_factor(499.67, 'R', 'F'), 40.)

        assert np.allclose(_temperature_factor(20., 'C', 'K'), 293.15)
        assert np.allclose(_temperature_factor(293.15, 'K', 'C'), 20.)

class TestAtm(unittest.TestCase):
    """various atmosphere tests"""
    def test_calibrated_airpseed(self):
        alt = 30000
        mach = 0.8
        cas = atm_calibrated_airspeed(alt, mach, alt_units='ft', cas_units='knots')
        self.assertAlmostEqual(cas, 303.990024, delta=0.001), cas

        mach2 = cas_to_mach(alt, cas, alt_units='ft', cas_units='knots')
        self.assertAlmostEqual(mach2, mach, delta=0.01)

    def test_temperature(self):
        """tests temperature at various altitudes"""
        with self.assertRaises(RuntimeError):
            atm_temperature(alt=0., alt_units='km', temperature_units='R')
        with self.assertRaises(RuntimeError):
            atm_temperature(alt=0., alt_units='ft', temperature_units='F')

        self.assertEqual(atm_temperature(alt=10 *1000.), 482.40003999999999)
        self.assertEqual(atm_temperature(alt=60 *1000.), 389.988)
        self.assertEqual(atm_temperature(alt=120*1000.), 451.2655824328092)
        self.assertEqual(atm_temperature(alt=165*1000.), 508.78800000000001)
        self.assertEqual(atm_temperature(alt=230*1000.), 394.18835930326833)
        self.assertEqual(atm_temperature(alt=270*1000.), 354.34800000000001)
        self.assertEqual(atm_temperature(alt=350*1000.), 354.34800000000001)

        self.assertEqual(atm_temperature(alt=0, alt_units='ft', temperature_units='K'),
                         287.77777777777777)

    def test_pressure(self):
        """tests pressure at various altitudes"""
        self.assertEqual(atm_pressure(alt=0., alt_units='ft', pressure_units='psf'),
                         2116.2247459927403)
        self.assertEqual(atm_pressure(alt=0., alt_units='m', pressure_units='psf'),
                         2116.2247459927403)

        units = {
            'alt_units': 'ft',
            'pressure_units': 'psf',
        }
        self.assertEqual(atm_pressure(alt=10*1000., **units), 1456.3074319943232)
        self.assertEqual(atm_pressure(alt=60*1000., **units), 151.20878913237249)
        self.assertEqual(atm_pressure(alt=120*1000., **units), 9.8627955961437763)
        self.assertEqual(atm_pressure(alt=165*1000., **units), 1.7725806687593277)
        self.assertEqual(atm_pressure(alt=230*1000., **units), 0.13023784776280109)
        self.assertEqual(atm_pressure(alt=270*1000., **units), 0.017353278750799964)
        self.assertEqual(atm_pressure(alt=350*1000., **units), 0.00028114006933161638)

        units = {
            'alt_units': 'kft',
            'pressure_units': 'psf',
        }
        self.assertEqual(atm_pressure(alt=10, **units), 1456.3074319943232)
        self.assertEqual(atm_pressure(alt=60, **units), 151.20878913237249)
        self.assertEqual(atm_pressure(alt=120, **units), 9.8627955961437763)
        self.assertEqual(atm_pressure(alt=165, **units), 1.7725806687593277)
        self.assertEqual(atm_pressure(alt=230, **units), 0.13023784776280109)
        self.assertEqual(atm_pressure(alt=270, **units), 0.017353278750799964)
        self.assertEqual(atm_pressure(alt=350, **units), 0.00028114006933161638)

        self.assertEqual(
            atm_pressure(alt=0., alt_units='m', pressure_units='Pa'),
            101325.20482878872)
        self.assertAlmostEqual(
            atm_pressure(alt=0., alt_units='m', pressure_units='kPa'),
            101.32520482878872, delta=0.001)
        self.assertAlmostEqual(
            atm_pressure(alt=0., alt_units='m', pressure_units='MPa'),
            0.10132520482878872, delta=0.000001)

        with self.assertRaises(RuntimeError):
            atm_pressure(alt=0., alt_units='m', pressure_units='kbar')

    def test_viscosity(self):
        """tests dynamic viscosity at various altitudes"""
        units = {
            'alt_units': 'ft',
            'visc_units': '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu(alt=0., **units), 3.7345965612371534e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=10 *1000., **units), 3.5317481186391660e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=60 *1000., **units), 2.9702384755729678e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=120*1000., **units), 3.3485025784164385e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=165*1000., **units), 3.6827603483595828e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=230*1000., **units), 2.9969664628998927e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=270*1000., **units), 2.7383347922674784e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=350*1000., **units), 2.7383347922674784e-07)

        units = {
            'alt_units': 'kft',
            'visc_units': '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu(alt=0., **units), 3.7345965612371534e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=10., **units), 3.5317481186391660e-07)

        units = {
            'alt_units': 'm',
            'visc_units': '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu(alt=0., **units), 3.7345965612371534e-07)

    def test_density(self):
        """tests density at various altitudes"""
        self.assertEqual(atm_density(299., alt_units='ft'), 0.0023600367621627013)
        self.assertEqual(atm_density(299./1000., alt_units='kft'), 0.0023600367621627013)

        rho2 = atm_density(0., alt_units='m', density_units='kg/m^3')
        assert np.allclose(rho2, 1.22699081123), rho2

        rho1 = atm_density(0.)
        assert np.allclose(rho1, 0.00238075522), rho1
        rho2 = atm_density(0., alt_units='ft', density_units='slinch/in^3')
        assert np.allclose(rho2, 0.00238075522/12.**4), rho2

    def test_velocity(self):
        """tests velocity at various altitudes"""
        vel_fts_55_2p4 = atm_velocity(55000., 2.4)
        #vel_fts_55_2 = atm_velocity(55000., 2.0)
        vel_fts_55_2b = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='ft/s')
        vel_fts_55_2c = atm_velocity(55., 2.0, alt_units='kft', velocity_units='ft/s')
        vel_ins_55_2c = atm_velocity(55., 2.0, alt_units='kft', velocity_units='in/s')
        assert np.allclose(vel_fts_55_2p4, 2323.05516639), vel_fts_55_2p4
        assert np.allclose(vel_fts_55_2b, 1935.87930533), vel_fts_55_2b
        assert np.allclose(vel_fts_55_2c, 1935.87930533), vel_fts_55_2c
        assert np.allclose(vel_fts_55_2b, vel_ins_55_2c/12.)

        #vel_fts_55_2b   = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='ft/s')
        vel_knots_55_2b = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='knots')
        assert np.allclose(vel_knots_55_2b, 1146.977130662127), vel_knots_55_2b

        vel_ms_55_2b = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='m/s')
        assert np.allclose(vel_ms_55_2b, 590.0560122641932), vel_ms_55_2b

        vel_mms_55_2b = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='mm/s')
        assert np.allclose(vel_mms_55_2b, 590.0560122641932 * 1000), vel_mms_55_2b

        with self.assertRaises(RuntimeError):
            atm_velocity(55000., 2.0, alt_units='ft', velocity_units='bad')

        self.assertAlmostEqual(vel_fts_55_2p4, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_fts_55_2p4, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_fts_55_2p4, vel_fts_55_2p4)
        self.assertAlmostEqual(vel_fts_55_2b, 1935.8, delta=0.1)

        with self.assertRaises(TypeError):
            atm_velocity(55000, [2.4])
        #self.assertRaises(AtmosphereError, atm_velocity, 55000, None)
        #-----
        vel_a = atm_velocity(55000., 2.4, alt_units='ft', velocity_units='ft/s')
        vel_b = atm_velocity(55000., 2.4, alt_units='ft', velocity_units='ft/s')
        vel_c = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='ft/s')

        self.assertAlmostEqual(vel_a, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_b, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_a, vel_b)
        self.assertAlmostEqual(vel_c, 1935.8, delta=0.1)

        vel_a = atm_velocity(55., 2.4, alt_units='kft', velocity_units='ft/s')
        vel_b = atm_velocity(55., 2.4, alt_units='kft', velocity_units='ft/s')
        vel_c = atm_velocity(55., 2.0, alt_units='kft', velocity_units='ft/s')

        self.assertAlmostEqual(vel_a, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_b, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_a, vel_b)
        self.assertAlmostEqual(vel_c, 1935.8, delta=0.1)

    def test_mach(self):
        """
        Tests mach at various speeds.
        """
        mach_a = atm_mach(14000., 743.011549709834) #, debug=False)
        mach_b = atm_mach(14000., 2122.8901420281)
        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)

        mach_a = atm_mach(14000., 743.011549709834, alt_units='ft', velocity_units='ft/s')
        mach_b = atm_mach(14000., 2122.8901420281, alt_units='ft', velocity_units='ft/s')
        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)

        mach_a = atm_mach(14., 743.011549709834, alt_units='kft', velocity_units='ft/s')
        mach_b = atm_mach(14., 2122.8901420281, alt_units='kft', velocity_units='ft/s')
        mach_c = atm_mach(14., 628.3419, alt_units='kft', velocity_units='knots')

        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)
        self.assertAlmostEqual(mach_c, 1.0, delta=0.002)

    def test_reynolds(self):
        """tests reynolds number"""
        self.assertEqual(
            atm_unit_reynolds_number(55000., 2.4, alt_units='ft', reynolds_units='1/ft'),
            2244166.3810534105)
        self.assertEqual(
            atm_unit_reynolds_number(55000., 2.4, alt_units='ft', reynolds_units='1/in'),
            187013.8650877842)
        self.assertEqual(
            atm_unit_reynolds_number(55000., 2.4, alt_units='ft', reynolds_units='1/m'),
            7362750.594007252)

        with self.assertRaises(RuntimeError):
            atm_unit_reynolds_number(55000., 2.4, alt_units='ft', reynolds_units='1/mm')

    def test_equiv_airspeed(self):
        """tests ``atm_equivalent_airspeed``"""
        alt = 0.
        mach = 1.
        veq1 = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='ft/s')
        assert np.allclose(veq1, 1115.5461442719434)

        alt = 10000.
        veq2 = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='ft/s')
        assert np.allclose(veq2, 925.407847794747)

        veq3 = atm_equivalent_airspeed(alt/1000., mach, alt_units='kft', eas_units='ft/s')
        assert np.allclose(veq2, veq3)

        alt = 0.
        veq4 = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='m/s')
        assert np.allclose(veq4, 340.0184647740884)

        alt = 0.
        veq4 = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='mm/s')
        assert np.allclose(veq4, 340.0184647740884 * 1000)

    def test_reynolds_factor(self):
        """tests ``test_reynolds_factor``"""
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/m', reynolds_units_out='1/m'), 1)
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/ft', reynolds_units_out='1/ft'), 1)
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/in', reynolds_units_out='1/in'), 1)

        assert np.allclose(_reynolds_factor(reynolds_units_in='1/m', reynolds_units_out='1/ft'), 0.3048)
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/m', reynolds_units_out='1/in'), 0.3048/12)

        assert np.allclose(_reynolds_factor(reynolds_units_in='1/ft', reynolds_units_out='1/m'), 1/0.3048)
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/ft', reynolds_units_out='1/in'), 1/12)

        assert np.allclose(_reynolds_factor(reynolds_units_in='1/in', reynolds_units_out='1/m'), 12/0.3048)
        assert np.allclose(_reynolds_factor(reynolds_units_in='1/in', reynolds_units_out='1/ft'), 12)

    def test_atm_kinematic_viscosity_nu(self):
        """tests ``atm_kinematic_viscosity_mu``"""
        mu = atm_kinematic_viscosity_nu(10., alt_units='kft', visc_units='ft^2/s')
        self.assertEqual(mu, 0.00020075264466844382)

        mu = atm_kinematic_viscosity_nu(5000., alt_units='m', visc_units='m^2/s')
        self.assertEqual(mu, 2.204293839480214e-05)

    def test_atm_dynamic_viscosity_mu(self):
        """tests ``test_atm_dynamic_viscosity_mu``"""
        mu1 = atm_dynamic_viscosity_mu(0., alt_units='ft', visc_units='(N*s)/m^2')
        mu2 = atm_dynamic_viscosity_mu(300., alt_units='kft', visc_units='(N*s)/m^2')
        self.assertEqual(mu1, 1.7881345434714084e-05), mu1
        self.assertEqual(mu2, 1.3111218182081286e-05), mu2

    def test_sutherland_viscoscity(self):
        """temperature (input) in rankine"""
        mu1 = sutherland_viscoscity(0)
        mu2 = sutherland_viscoscity(100)
        mu3 = sutherland_viscoscity(6000)  # rankine
        self.assertEqual(mu1, 0.0), mu1
        self.assertEqual(mu2, 8.0382436e-08), mu2
        self.assertEqual(mu3, 1.7019982955940058e-06), mu3

    def test_get_alt_for_density(self):
        """tests ``get_alt_for_density``"""
        alt_targets = [
            0., 10., 20., 30., 40., 50.,
            300., 350., 400.,
        ]
        for alt_target in alt_targets:
            rho1 = atm_density(alt_target * 1000.)
            alt1 = get_alt_for_density(rho1)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            rho2 = atm_density(alt_target, alt_units='kft', density_units='kg/m^3')
            tol = 0.005  # 5 feet
            alt2 = get_alt_for_density(rho2, density_units='kg/m^3', alt_units='kft', tol=tol, nmax=50)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_get_alt_for_mach_eas(self):
        """tests ``get_alt_for_mach_eas``"""
        mach = 0.8
        eas_target = 500.  # knots
        alt_expected = 3067.215571329515  # ft
        alt = get_alt_for_mach_eas(mach, eas_target, alt_units='ft', eas_units='knots', tol=1e-12)
        assert np.allclose(alt, alt_expected)
        #print(f'alt={alt}')
        eas = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='knots')
        #print(f'eas={eas}')
        assert np.allclose(eas, eas_target)

    def test_get_alt_for_pressure(self):
        """tests ``get_alt_for_pressure``"""
        alt_targets = [
            0., 10., 20., 30., 40., 50.,
            300., 350., #400.,
        ]
        for alt_target in alt_targets:
            pressure1 = atm_pressure(alt_target*1000.)
            #alt1 = get_alt_for_pressure(pressure1, tol=5., SI=False, nmax=20)
            alt1 = get_alt_for_pressure(pressure1, nmax=20)

            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            pressure2 = atm_pressure(alt_target, alt_units='kft', pressure_units='Pa')

            tol = 0.005 # 5 feet
            alt2 = get_alt_for_pressure(pressure2, pressure_units='Pa', alt_units='kft', tol=tol)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_get_alt_for_q_with_constant_mach(self):
        """tests ``get_alt_for_q_with_constant_mach`` for various altitudes"""
        alt = get_alt_for_q_with_constant_mach(
            534.2, 0.8, pressure_units='psf', alt_units='ft', nmax=20, tol=5.)
        assert np.allclose(alt, 15064.6441438)

        tol = 0.001 # in feet
        alt_a = get_alt_for_q_with_constant_mach(
            600., 0.8, pressure_units='psf', alt_units='ft', nmax=20, tol=tol)
        alt_b = get_alt_for_q_with_constant_mach(
            400., 0.8, pressure_units='psf', alt_units='ft', nmax=20, tol=tol)

        self.assertAlmostEqual(alt_a, 12144.30, delta=1.25)  # TODO: should be 0.01
        self.assertAlmostEqual(alt_b, 22058.47, delta=2.25)  # TODO: should be 0.01

    def test_get_alt_for_q_with_constant_mach2(self):
        """tests ``get_alt_for_q_with_constant_mach``"""
        mach = 0.8
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            pressure1 = atm_dynamic_pressure(alt_target*1000., mach)
            alt1 = get_alt_for_q_with_constant_mach(pressure1, mach)
            #self.assertAlmostEqual(alt1, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            pressure2 = atm_dynamic_pressure(alt_target, mach,
                                             alt_units='kft', pressure_units='psi')

            tol = 0.005 # 5 feet
            alt2 = get_alt_for_q_with_constant_mach(
                pressure2, mach, pressure_units='psi', alt_units='kft', tol=tol)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

        #get_alt_for_q_mach(q, mach, pressure_units='psf', alt_units='ft')

    def test_get_alt_for_eas_with_constant_mach(self):
        """tests get_alt_for_q_with_constant_mach"""
        mach = 0.8
        alt_targets = [
            0., 10., 20., 30., 40., 50.,
            300., 350., 400.,
        ]
        for alt_target in alt_targets:
            veq1 = atm_equivalent_airspeed(
                alt_target*1000., mach, alt_units='ft', eas_units='ft/s')
            veq2 = atm_equivalent_airspeed(
                alt_target, mach, alt_units='kft', eas_units='knots')
            tol = 5. # 5 feet
            alt1 = get_alt_for_eas_with_constant_mach(
                veq1, mach, velocity_units='ft/s', alt_units='ft', nmax=20, tol=tol)
            tol = 5/1000. # 5 feet
            alt2 = get_alt_for_eas_with_constant_mach(
                veq2, mach, velocity_units='knots', alt_units='kft', nmax=20, tol=tol)

            assert np.allclose(alt1/1000., alt_target, atol=1e-3), 'alt1=%s alt_target=%s' % (alt1, alt_target)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_get_alt_for_eas_with_constant_mach2(self):
        """tests get_alt_for_eas_with_constant_mach"""
        eas_expected = 500.
        mach = 0.8
        alt = get_alt_for_eas_with_constant_mach(
            eas_expected, mach,
            velocity_units='ft/s', alt_units='ft', nmax=20, tol=5.)
        eas = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='ft/s')
        assert np.allclose(eas, eas_expected), 'eas=%s eas_expected=%s' % (eas, eas_expected)

        eas = atm_equivalent_airspeed(alt, mach, alt_units='ft', eas_units='knots')
        assert np.allclose(eas, 296.2418795250434), eas

    def test_atm_unit_reynolds_number(self):
        """tests ``atm_unit_reynolds_number and atm_unit_reynolds_number``"""
        mach = 0.8
        targets = [
            # (alt_kft, Re_1/ft)
            (0., 5689165.64362),
            (10., 4289977.7416),
            (20., 3169466.53065),
            (30., 2286823.00426),
            (40., 1532035.46128),
            (50., 949974.915093)]
        for alt, re_expected in targets:
            rel_a = atm_unit_reynolds_number(alt*1000., mach)
            assert np.allclose(rel_a, re_expected, atol=1e-3), 'rel_a=%s re_expected=%s' % (rel_a, re_expected)

            rel_b = atm_unit_reynolds_number(alt, mach, alt_units='kft', reynolds_units='1/ft')
            assert np.allclose(rel_b, re_expected, atol=1e-3), 'rel_b=%s re_expected=%s' % (rel_b, re_expected)

    def test_sweep_eas_mach(self):
        eass = np.linspace(0., 1000., num=101, dtype='float64')
        eass[0] = 1e-5
        #neas = len(eass)

        machi = 0.5
        #machs = np.ones(neas, dtype='float64') * mach
        rhos1, machs1, velocity1, alts = make_flfacts_eas_sweep_constant_mach(
            machi, eass, gamma=1.4,
            velocity_units='ft/s', density_units='slug/ft^3',
            alt_units='ft',
            #pressure_units='psf',
            eas_units='knots')

        rhos2, machs2, velocity2 = make_flfacts_alt_sweep_constant_mach(
            machi, alts,
            eas_limit=1000000,
            velocity_units='ft/s', density_units='slug/ft^3',
            alt_units='ft',
            #pressure_units='psf',
            eas_units='knots')
        assert np.allclose(rhos1, rhos2)
        assert np.allclose(machs1, machs2)
        #assert np.allclose(velocity1, velocity2)
        #-----------------------------------------------
        rhos1, machs1, velocity1, alts = make_flfacts_eas_sweep_constant_mach(
            machi, eass, gamma=1.4,
            velocity_units='m/s', density_units='kg/m^3',
            alt_units='m',
            #pressure_units='psf',
            eas_units='m/s')

        rhos2, machs2, velocity2 = make_flfacts_alt_sweep_constant_mach(
            machi, alts,
            eas_limit=1000000,
            velocity_units='m/s', density_units='kg/m^3',
            alt_units='ft',
            #pressure_units='psf',
            eas_units='m/s')
        #print(machs1-machs2)
        #print(velocity1-velocity2)
        #assert np.allclose(rhos1, rhos2)
        assert np.allclose(machs1, machs2)
        #assert np.allclose(velocity1, velocity2)

    def test_sweep(self):
        """tests FLFACT sweeps"""
        mach = 0.8
        alts = np.linspace(-10000, 80000.)[::-1]
        #make_flfacts_alt_sweep_constant_mach
        rho, mach, vel = make_flfacts_alt_sweep_constant_mach(
            mach, alts, eas_limit=300., alt_units='m',
            velocity_units='m/s',
            density_units='kg/m^3',
            eas_units='m/s')
        vel_expected = [224.97932124482696, 224.97932124482696, 224.97932124482696, 228.02811169475368, 231.9519431728144, 235.81049204862362, 239.6069122099433, 243.34411149101786, 247.02477773681403, 250.65140142136016, 254.22629536116307, 257.75161196708507, 261.2293584002385, 264.6614099349766, 268.04952178159016, 269.58519565506737, 269.58519565506737, 269.58519565506737, 269.0029125642663, 266.3865851329001, 263.74430516975593, 261.0752846980864, 258.37869504062513, 255.65366381416462, 252.8992716326313, 250.1145484833398, 247.2984697359336, 244.4499517374548, 241.5678469398473, 238.65093849776292, 236.0224049056773, 236.0224049056773, 236.0224049056773, 236.0224049056773, 236.0224049056773, 236.0224049056773, 236.0224049056773, 236.0224049056773, 238.67960652624149, 245.01487555752158, 251.19041366543848, 257.2177260422781, 263.1069996471146, 268.8673054431369, 274.5067624106907]
        mach_expected = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
        assert np.allclose(vel, vel_expected), vel.tolist()
        assert np.allclose(mach, mach_expected), mach.tolist()
        del rho, mach, vel, alts

        alt = 10000.
        machs = np.linspace(0., 0.8)
        #make_flfacts_mach_sweep_constant_alt
        rho, mach, vel = make_flfacts_mach_sweep_constant_alt(
            alt, machs, eas_limit=300.,
            alt_units='m',
            velocity_units='m/s',
            density_units='kg/m^3',
            eas_units='m/s')
        del rho, mach, vel, alt

        alt = 0.
        machs = np.linspace(0.6, 0.8)
        with self.assertRaises(RuntimeError):
            rho, mach, vel = make_flfacts_mach_sweep_constant_alt(
                alt, machs, eas_limit=100.,
                alt_units='m',
                velocity_units='m/s',
                density_units='kg/m^3',
                eas_units='m/s')
        del alt

        alt = 10000.
        eass = np.linspace(0., 300.)
        rho, mach, vel = make_flfacts_eas_sweep_constant_alt(
            alt, eass, alt_units='m',
            velocity_units='m/s', density_units='kg/m^3',
            eas_units='m/s')
        del rho, mach, vel, alt

    def test_table(self):
        dirname = os.path.dirname(__file__)
        quantities = ['density', 'pressure', 'speed_of_sound', 'temperature',
                      'velocity', 'dynamic_viscosity',
                      'calibrated_airspeed', 'equivalent_airspeed']
        x = create_atmosphere_table(
            quantities, 0.8,
            0., 50000.,
            dalt=1000.)
        csv_filename = os.path.join(dirname, 'atmosphere.csv')
        quantity_list = ['alt'] + quantities
        header = ','.join(quantity_list)
        np.savetxt(csv_filename, x, delimiter=',', header=header)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
