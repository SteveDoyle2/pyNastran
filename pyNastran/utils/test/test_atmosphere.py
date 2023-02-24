"""various atmosphere tests"""
import unittest
import numpy as np

from pyNastran.utils.atmosphere import (
    atm_density, atm_dynamic_pressure, atm_temperature,
    atm_pressure, atm_velocity, atm_mach, atm_equivalent_airspeed,
    atm_dynamic_viscosity_mu, atm_kinematic_viscosity_nu,
    get_alt_for_density, get_alt_for_pressure,
    get_alt_for_q_with_constant_mach,
    get_alt_for_eas_with_constant_mach,
    atm_unit_reynolds_number, atm_unit_reynolds_number2,
    make_flfacts_alt_sweep_constant_mach,
    make_flfacts_mach_sweep_constant_alt,
    make_flfacts_eas_sweep_constant_alt,
)

from pyNastran.utils.convert import (
    convert_length, convert_density, convert_mass,
    convert_velocity, convert_force, convert_pressure,
)
class TestConvert(unittest.TestCase):
    """various unit conversion tests"""

    def test_length(self):
        """length checks"""
        assert np.allclose(convert_length(1., 'ft', 'ft'), 1.)
        with self.assertRaises(RuntimeError):
            convert_length(1., 'ft', 'bad')
        with self.assertRaises(RuntimeError):
            convert_length(1., 'bad', 'ft')

        assert np.allclose(convert_length(1., 'ft', 'in'), 12.)
        assert np.allclose(convert_length(1., 'in', 'ft'), 1 / 12.)

        assert np.allclose(convert_length(1., 'm', 'ft'), 3.28084)
        assert np.allclose(convert_length(1., 'ft', 'm'), 1 / 3.28084)

        assert np.allclose(convert_length(1., 'm', 'ft'), 3.28084)
        assert np.allclose(convert_length(1., 'ft', 'm'), 1 / 3.28084)

        assert np.allclose(convert_length(1., 'm', 'cm'), 100.)
        assert np.allclose(convert_length(1., 'cm', 'm'), 1 / 100.)

        #assert np.allclose(convert_length(1., 'm', 'mm'), 1000.)
        #assert np.allclose(convert_length(1., 'mm', 'm'), 1 / 1000.)

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

    def test_force(self):
        """force checks"""
        assert np.allclose(convert_force(1., 'lbf', 'lbf'), 1.)
        with self.assertRaises(RuntimeError):
            convert_force(1., 'lbf', 'bad')
        with self.assertRaises(RuntimeError):
            convert_force(1., 'bad', 'lbf')

        assert np.allclose(convert_force(1., 'lbf', 'N'), 4.44822)
        assert np.allclose(convert_force(1., 'N', 'lbf'), 1 / 4.44822)

        assert np.allclose(convert_force(1., 'MN', 'N'), 1000000.)
        assert np.allclose(convert_force(1., 'N', 'MN'), 1 / 1000000.)

        assert np.allclose(convert_force(1., 'mN', 'N'), 0.001)
        assert np.allclose(convert_force(1., 'N', 'mN'), 1 / 0.001)

        assert np.allclose(convert_force(1., 'cN', 'N'), 0.01)
        assert np.allclose(convert_force(1., 'N', 'cN'), 1 / 0.01)

    def test_pressure(self):
        """pressure checks"""
        assert np.allclose(convert_pressure(1., 'psf', 'psf'), 1.)
        with self.assertRaises(RuntimeError):
            convert_pressure(1., 'psf', 'bad')
        with self.assertRaises(RuntimeError):
            convert_pressure(1., 'bad', 'psf')

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

        assert np.allclose(convert_velocity(1., 'ft/s', 'in/s'), 12.)
        assert np.allclose(convert_velocity(1., 'in/s', 'ft/s'), 1 / 12.)

        assert np.allclose(convert_velocity(1., 'm/s', 'ft/s'), 3.28084)
        assert np.allclose(convert_velocity(1., 'ft/s', 'm/s'), 1 / 3.28084)

        assert np.allclose(convert_velocity(1., 'm/s', 'ft/s'), 3.28084)
        assert np.allclose(convert_velocity(1., 'ft/s', 'm/s'), 1 / 3.28084)

        assert np.allclose(convert_velocity(1., 'm/s', 'cm/s'), 100.)
        assert np.allclose(convert_velocity(1., 'cm/s', 'm/s'), 1 / 100.)

        #assert np.allclose(convert_velocity(1., 'm/s', 'mm/s'), 1000.)
        #assert np.allclose(convert_velocity(1., 'mm/s', 'm/s'), 1 / 1000.)

        assert np.allclose(convert_velocity(1., 'knots', 'ft/s'), 1.68781)
        assert np.allclose(convert_velocity(1., 'ft/s', 'knots'), 1 / 1.68781)

    def test_density(self):
        """density checks"""
        assert np.allclose(convert_density(1., 'slug/ft^3', 'slug/ft^3'), 1.)
        with self.assertRaises(RuntimeError):
            convert_density(1., 'slug/ft^3', 'bad')
        with self.assertRaises(RuntimeError):
            convert_density(1., 'bad', 'slug/ft^3')

        assert np.allclose(convert_density(1., 'slinch/in^3', 'slug/ft^3'), 12**4)
        assert np.allclose(convert_density(1., 'slug/ft^3', 'slinch/in^3'), 1 / 12**4)

        assert np.allclose(convert_density(1., 'slug/ft^3', 'kg/m^3'), 515.379)
        assert np.allclose(convert_density(1., 'kg/m^3', 'slug/ft^3'), 1 / 515.379)

        assert np.allclose(convert_density(1., 'g/cm^3', 'kg/m^3'), 1000.), 'actual=%g expected=%g' % (convert_density(1., 'g/cm^3', 'kg/m^3'), 1000.)
        assert np.allclose(convert_density(1., 'kg/m^3', 'g/cm^3'), 1 / 1000.), 'actual=%g expected=%g' % (convert_density(1., 'kg/m^3', 'g/cm^3'), 1 / 1000.)

        assert np.allclose(convert_density(1., 'g/cm^3', 'slug/ft^3'), 1.94032)
        assert np.allclose(convert_density(1., 'slug/ft^3', 'g/cm^3'), 1 / 1.94032), 'actual=%g expected=%g' % (convert_density(1., 'slug/ft^3', 'g/cm^3'), 1 / 1.94032)

class TestAtm(unittest.TestCase):
    """various atmosphere tests"""
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
            'alt_units' : 'ft',
            'pressure_units' : 'psf',
        }
        self.assertEqual(atm_pressure(alt=10*1000., **units), 1456.3074319943232)
        self.assertEqual(atm_pressure(alt=60*1000., **units), 151.20878913237249)
        self.assertEqual(atm_pressure(alt=120*1000., **units), 9.8627955961437763)
        self.assertEqual(atm_pressure(alt=165*1000., **units), 1.7725806687593277)
        self.assertEqual(atm_pressure(alt=230*1000., **units), 0.13023784776280109)
        self.assertEqual(atm_pressure(alt=270*1000., **units), 0.017353278750799964)
        self.assertEqual(atm_pressure(alt=350*1000., **units), 0.00028114006933161638)

        units = {
            'alt_units' : 'kft',
            'pressure_units' : 'psf',
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
            atm_pressure(alt=0., alt_units='m', pressure_units='bar')

    def test_viscosity(self):
        """tests dynamic viscosity at various altitudes"""
        units = {
            'alt_units' : 'ft',
            'visc_units' : '(lbf*s)/ft^2',
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
            'alt_units' : 'kft',
            'visc_units' : '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu(alt=0., **units), 3.7345965612371534e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=10., **units), 3.5317481186391660e-07)

        units = {
            'alt_units' : 'm',
            'visc_units' : '(lbf*s)/ft^2',
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
        vel_ms_55_2b = atm_velocity(55000., 2.0, alt_units='ft', velocity_units='m/s')
        assert np.allclose(vel_knots_55_2b, 1146.977130662127), vel_knots_55_2b
        assert np.allclose(vel_ms_55_2b, 590.0560122641932), vel_ms_55_2b

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
        """tests atm_equivalent_airspeed"""
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

    def test_atm_kinematic_viscosity_nu(self):
        """tests ``atm_kinematic_viscosity_mu``"""
        mu = atm_kinematic_viscosity_nu(10., alt_units='kft', visc_units='ft^2/s')
        self.assertEqual(mu, 0.00020075264466844382)

        mu = atm_kinematic_viscosity_nu(5000., alt_units='m', visc_units='m^2/s')
        self.assertEqual(mu, 2.204293839480214e-05)

    def test_get_alt_for_density(self):
        """tests ``get_alt_for_density``"""
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            rho1 = atm_density(alt_target * 1000.)
            alt1 = get_alt_for_density(rho1)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            rho2 = atm_density(alt_target, alt_units='kft', density_units='kg/m^3')
            tol = 0.005 # 5 feet
            alt2 = get_alt_for_density(rho2, density_units='kg/m^3', alt_units='kft', tol=tol)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_get_alt_for_pressure(self):
        """tests ``get_alt_for_pressure``"""
        alt_targets = [0., 10., 20., 30., 40., 50.]
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
        alt_targets = [0., 10., 20., 30., 40., 50.]
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
            rel_a = atm_unit_reynolds_number2(alt*1000., mach)
            assert np.allclose(rel_a, re_expected, atol=1e-3), 'rel_a=%s re_expected=%s' % (rel_a, re_expected)

            rel_b = atm_unit_reynolds_number2(alt, mach, alt_units='kft', reynolds_units='1/ft')
            assert np.allclose(rel_b, re_expected, atol=1e-3), 'rel_b=%s re_expected=%s' % (rel_b, re_expected)

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

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
