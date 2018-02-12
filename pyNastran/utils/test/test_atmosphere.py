"""various atmosphere tests"""
import unittest
import numpy as np

from pyNastran.utils.atmosphere import (
    atm_density, atm_temperature, atm_pressure,
    atm_dynamic_viscosity_mu, atm_dynamic_pressure,
    atm_velocity, atm_mach,
    atm_unit_reynolds_number,
    atm_equivalent_airspeed,
    atm_kinematic_viscosity_nu,
    get_alt_for_density,
    get_alt_for_pressure,
    get_alt_for_q_with_constant_mach,
    get_alt_for_eas_with_constant_mach,
    atm_unit_reynolds_number as atm_unit_reynolds_number1a,
    atm_unit_reynolds_number2 as atm_unit_reynolds_number2a,
)
from pyNastran.utils.atmosphere2 import (
    atm_density as atm_density2,
    atm_dynamic_pressure as atm_dynamic_pressure2,
    atm_temperature as atm_temperature2,
    atm_pressure as atm_pressure2,
    atm_velocity as atm_velocity2,
    atm_mach as atm_mach2,
    atm_dynamic_viscosity_mu as atm_dynamic_viscosity_mu2,
    atm_kinematic_viscosity_nu as atm_kinematic_viscosity_nu2,
    atm_equivalent_airspeed as atm_equivalent_airspeed2,
    get_alt_for_density as get_alt_for_density2,
    get_alt_for_pressure as get_alt_for_pressure2,
    get_alt_for_q_with_constant_mach as get_alt_for_q_with_constant_mach2,
    get_alt_for_eas_with_constant_mach as get_alt_for_eas_with_constant_mach2,
    atm_unit_reynolds_number as atm_unit_reynolds_number1b,
    atm_unit_reynolds_number2 as atm_unit_reynolds_number2b,
    make_flfacts_alt_sweep, make_flfacts_mach_sweep,
)


class TestAtm(unittest.TestCase):
    """various atmosphere tests"""
    def test_temperature(self):
        """tests temperature at various altitudes"""
        self.assertEqual(atm_temperature(alt=10 *1000.), 482.40003999999999)
        self.assertEqual(atm_temperature(alt=60 *1000.), 389.988)
        self.assertEqual(atm_temperature(alt=120*1000.), 451.2655824328092)
        self.assertEqual(atm_temperature(alt=165*1000.), 508.78800000000001)
        self.assertEqual(atm_temperature(alt=230*1000.), 394.18835930326833)
        self.assertEqual(atm_temperature(alt=270*1000.), 354.34800000000001)
        self.assertEqual(atm_temperature(alt=350*1000.), 354.34800000000001)

        self.assertEqual(atm_temperature2(alt=10 *1000.), 482.40003999999999)
        self.assertEqual(atm_temperature2(alt=60 *1000.), 389.988)
        self.assertEqual(atm_temperature2(alt=120*1000.), 451.2655824328092)
        self.assertEqual(atm_temperature2(alt=165*1000.), 508.78800000000001)
        self.assertEqual(atm_temperature2(alt=230*1000.), 394.18835930326833)
        self.assertEqual(atm_temperature2(alt=270*1000.), 354.34800000000001)
        self.assertEqual(atm_temperature2(alt=350*1000.), 354.34800000000001)

    def test_pressure(self):
        """tests pressure at various altitudes"""
        self.assertEqual(atm_pressure(alt=10 *1000.), 1456.3074319943232)
        self.assertEqual(atm_pressure(alt=60 *1000.), 151.20878913237249)
        self.assertEqual(atm_pressure(alt=120*1000.), 9.8627955961437763)
        self.assertEqual(atm_pressure(alt=165*1000.), 1.7725806687593277)
        self.assertEqual(atm_pressure(alt=230*1000.), 0.13023784776280109)
        self.assertEqual(atm_pressure(alt=270*1000.), 0.017353278750799964)
        self.assertEqual(atm_pressure(alt=350*1000.), 0.00028114006933161638)

        self.assertEqual(atm_pressure2(alt=0., alt_units='ft', pressure_units='psf'),
                         2116.2247459927403)
        self.assertEqual(atm_pressure2(alt=0., alt_units='m', pressure_units='psf'),
                         2116.2247459927403)

        units = {
            'alt_units' : 'ft',
            'pressure_units' : 'psf',
        }
        self.assertEqual(atm_pressure2(alt=10*1000., **units), 1456.3074319943232)
        self.assertEqual(atm_pressure2(alt=60*1000., **units), 151.20878913237249)
        self.assertEqual(atm_pressure2(alt=120*1000., **units), 9.8627955961437763)
        self.assertEqual(atm_pressure2(alt=165*1000., **units), 1.7725806687593277)
        self.assertEqual(atm_pressure2(alt=230*1000., **units), 0.13023784776280109)
        self.assertEqual(atm_pressure2(alt=270*1000., **units), 0.017353278750799964)
        self.assertEqual(atm_pressure2(alt=350*1000., **units), 0.00028114006933161638)

        units = {
            'alt_units' : 'kft',
            'pressure_units' : 'psf',
        }
        self.assertEqual(atm_pressure2(alt=10, **units), 1456.3074319943232)
        self.assertEqual(atm_pressure2(alt=60, **units), 151.20878913237249)
        self.assertEqual(atm_pressure2(alt=120, **units), 9.8627955961437763)
        self.assertEqual(atm_pressure2(alt=165, **units), 1.7725806687593277)
        self.assertEqual(atm_pressure2(alt=230, **units), 0.13023784776280109)
        self.assertEqual(atm_pressure2(alt=270, **units), 0.017353278750799964)
        self.assertEqual(atm_pressure2(alt=350, **units), 0.00028114006933161638)

    def test_viscosity(self):
        """tests dynamic viscosity at various altitudes"""
        self.assertEqual(atm_dynamic_viscosity_mu(alt=10 *1000.), 3.5317481186391660e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=60 *1000.), 2.9702384755729678e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=120*1000.), 3.3485025784164385e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=165*1000.), 3.6827603483595828e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=230*1000.), 2.9969664628998927e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=270*1000.), 2.7383347922674784e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=350*1000.), 2.7383347922674784e-07)

        units = {
            'alt_units' : 'ft',
            'visc_units' : '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=0., **units), 3.7345965612371534e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=10 *1000., **units), 3.5317481186391660e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=60 *1000., **units), 2.9702384755729678e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=120*1000., **units), 3.3485025784164385e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=165*1000., **units), 3.6827603483595828e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=230*1000., **units), 2.9969664628998927e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=270*1000., **units), 2.7383347922674784e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=350*1000., **units), 2.7383347922674784e-07)

        units = {
            'alt_units' : 'kft',
            'visc_units' : '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=0., **units), 3.7345965612371534e-07)
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=10., **units), 3.5317481186391660e-07)

        units = {
            'alt_units' : 'm',
            'visc_units' : '(lbf*s)/ft^2',
        }
        self.assertEqual(atm_dynamic_viscosity_mu2(alt=0., **units), 3.7345965612371534e-07)

    def test_density(self):
        """tests density at various altitudes"""
        self.assertEqual(atm_density(299.), 0.0023600367621627013)
        self.assertEqual(atm_density2(299., alt_units='ft'), 0.0023600367621627013)
        self.assertEqual(atm_density2(299./1000., alt_units='kft'), 0.0023600367621627013)

        rho1 = atm_density(0., SI=True)
        rho2 = atm_density2(0., alt_units='m', density_units='kg/m^3')
        assert np.allclose(rho1, 1.22699081123), rho1
        assert np.allclose(rho2, 1.22699081123), rho2

        rho1 = atm_density(0.)
        assert np.allclose(rho1, 0.00238075522), rho1
        rho2 = atm_density2(0., alt_units='ft', density_units='slinch/in^3')
        assert np.allclose(rho2, 0.00238075522/12.**4), rho2

    def test_velocity(self):
        """tests velocity at various altitudes"""
        vel_fts_55_2p4 = atm_velocity(55000., 2.4)
        vel_fts_55_2p4 = atm_velocity(55000., 2.4)
        vel_fts_55_2 = atm_velocity(55000., 2.0)
        vel_fts_55_2b = atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='ft/s')
        vel_fts_55_2c = atm_velocity2(55., 2.0, alt_units='kft', velocity_units='ft/s')
        vel_ins_55_2c = atm_velocity2(55., 2.0, alt_units='kft', velocity_units='in/s')
        assert np.allclose(vel_fts_55_2, vel_fts_55_2b)
        assert np.allclose(vel_fts_55_2, vel_fts_55_2c)
        assert np.allclose(vel_fts_55_2, vel_ins_55_2c/12.)

        #vel_fts_55_2b   = atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='ft/s')
        vel_knots_55_2b = atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='knots')
        vel_ms_55_2b = atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='m/s')
        assert np.allclose(vel_knots_55_2b, 1146.977130662127), vel_knots_55_2b
        assert np.allclose(vel_ms_55_2b, 590.0560122641932), vel_ms_55_2b

        with self.assertRaises(RuntimeError):
            atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='bad')


        self.assertAlmostEqual(vel_fts_55_2p4, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_fts_55_2p4, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_fts_55_2p4, vel_fts_55_2p4)
        self.assertAlmostEqual(vel_fts_55_2, 1935.8, delta=0.1)

        with self.assertRaises(TypeError):
            atm_velocity(55000, [2.4])
        #self.assertRaises(AtmosphereError, atm_velocity, 55000, None)
        #-----
        vel_a = atm_velocity2(55000., 2.4, alt_units='ft', velocity_units='ft/s')
        vel_b = atm_velocity2(55000., 2.4, alt_units='ft', velocity_units='ft/s')
        vel_c = atm_velocity2(55000., 2.0, alt_units='ft', velocity_units='ft/s')

        self.assertAlmostEqual(vel_a, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_b, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_a, vel_b)
        self.assertAlmostEqual(vel_c, 1935.8, delta=0.1)

        vel_a = atm_velocity2(55., 2.4, alt_units='kft', velocity_units='ft/s')
        vel_b = atm_velocity2(55., 2.4, alt_units='kft', velocity_units='ft/s')
        vel_c = atm_velocity2(55., 2.0, alt_units='kft', velocity_units='ft/s')

        self.assertAlmostEqual(vel_a, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_b, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_a, vel_b)
        self.assertAlmostEqual(vel_c, 1935.8, delta=0.1)

    def test_mach(self):
        """
        @todo fix error in function related to assertAlmostEqual
        """
        mach_a = atm_mach(14000., 743.011549709834, debug=False)
        mach_b = atm_mach(14000., 2122.8901420281)
        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)

        mach_a = atm_mach2(14000., 743.011549709834, alt_units='ft', velocity_units='ft/s')
        mach_b = atm_mach2(14000., 2122.8901420281, alt_units='ft', velocity_units='ft/s')
        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)

        mach_a = atm_mach2(14., 743.011549709834, alt_units='kft', velocity_units='ft/s')
        mach_b = atm_mach2(14., 2122.8901420281, alt_units='kft', velocity_units='ft/s')
        mach_c = atm_mach2(14., 628.3419, alt_units='kft', velocity_units='knots')

        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)
        self.assertAlmostEqual(mach_c, 1.0, delta=0.002)

    def test_reynolds(self):
        """tests reynolds number"""
        reynolds = atm_unit_reynolds_number(55000., 2.4)
        self.assertEqual(reynolds, 2244166.3810534105)

    def test_equiv_airspeed(self):
        """tests atm_equivalent_airspeed"""
        alt = 0.
        mach = 1.
        veq1 = atm_equivalent_airspeed(alt, mach)
        veq2 = atm_equivalent_airspeed2(alt, mach, alt_units='ft', eas_units='ft/s')
        assert np.allclose(veq1, veq2)

        alt = 10000.
        veq1 = atm_equivalent_airspeed(alt, mach)
        veq2 = atm_equivalent_airspeed2(alt, mach, alt_units='ft', eas_units='ft/s')
        assert np.allclose(veq1, veq2)

        veq2 = atm_equivalent_airspeed2(alt/1000., mach, alt_units='kft', eas_units='ft/s')
        assert np.allclose(veq1, veq2)

        alt = 0.
        veq1 = atm_equivalent_airspeed(alt, mach, SI=True)
        veq2 = atm_equivalent_airspeed2(alt, mach, alt_units='ft', eas_units='m/s')
        assert np.allclose(veq1, veq2)

    def test_atm_kinematic_viscosity_nu(self):
        """tests atm_kinematic_viscosity_mu"""
        mu1 = atm_kinematic_viscosity_nu(10000.)
        mu2 = atm_kinematic_viscosity_nu2(10., alt_units='kft', visc_units='ft^2/s')
        self.assertEqual(mu1, mu2)

        mu1 = atm_kinematic_viscosity_nu(5000., SI=True)
        mu2 = atm_kinematic_viscosity_nu2(5000., alt_units='m', visc_units='m^2/s')
        self.assertEqual(mu1, mu2)

    def test_get_alt_for_density(self):
        """tests get_alt_for_density"""
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            rho1 = atm_density(alt_target * 1000.)
            alt1 = get_alt_for_density(rho1)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            rho2 = atm_density2(alt_target, alt_units='kft', density_units='kg/m^3')
            alt2 = get_alt_for_density2(rho2, density_units='kg/m^3', alt_units='kft')
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_get_alt_for_pressure(self):
        """tests get_alt_for_pressure"""
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            pressure1 = atm_pressure(alt_target*1000.)
            alt1 = get_alt_for_pressure(pressure1, tol=5., SI=False, nmax=20)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            pressure2 = atm_pressure2(alt_target, alt_units='kft', pressure_units='Pa')
            alt2 = get_alt_for_pressure2(pressure2, pressure_units='Pa', alt_units='kft')
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)


    def test_get_alt_for_q_with_constant_mach(self):
        """tests get_alt_for_q_with_constant_mach for various altitudes"""
        get_alt_for_q_with_constant_mach(534.2, 0.8)
        alt_a = get_alt_for_q_with_constant_mach(600., 0.8, tol=0.01)
        alt_b = get_alt_for_q_with_constant_mach(400., 0.8, tol=0.01)

        self.assertAlmostEqual(alt_a, 12144.30, delta=1.5)  # TODO: should be 0.01
        self.assertAlmostEqual(alt_b, 22058.47, delta=2.5)  # TODO: should be 0.01

    def test_get_alt_for_q_with_constant_mach2(self):
        """tests get_alt_for_q_with_constant_mach"""
        mach = 0.8
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            pressure1 = atm_dynamic_pressure(alt_target*1000., mach)
            alt1 = get_alt_for_q_with_constant_mach(pressure1, mach)
            #self.assertAlmostEqual(alt, alt_target)
            assert np.allclose(alt1, alt_target*1000, atol=1.), 'alt1=%s alt_target=%s' % (alt1, alt_target*1000)

            pressure2 = atm_dynamic_pressure2(alt_target, mach, alt_units='kft', pressure_units='psi')
            alt2 = get_alt_for_q_with_constant_mach2(pressure2, mach, pressure_units='psi', alt_units='kft')
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

        #get_alt_for_q_mach(q, mach, pressure_units='psf', alt_units='ft')

    def get_alt_for_eas_with_constant_mach(self):
        """tests get_alt_for_q_with_constant_mach"""
        mach = 0.8
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt_target in alt_targets:
            veq1 = atm_equivalent_airspeed(alt_target*1000., mach)
            veq2 = atm_equivalent_airspeed2(alt_target, mach, alt_units='kft', eas_units='knots')
            alt1 = get_alt_for_eas_with_constant_mach(veq1, mach)
            alt2 = get_alt_for_eas_with_constant_mach2(veq2, mach, velocity_units='knots', alt_units='kft')

            assert np.allclose(alt1/1000., alt_target, atol=1e-3), 'alt1=%s alt_target=%s' % (alt1, alt_target)
            assert np.allclose(alt2, alt_target, atol=1e-3), 'alt2=%s alt_target=%s' % (alt2, alt_target)

    def test_atm_unit_reynolds_number(self):
        """tests atm_unit_reynolds_number and atm_unit_reynolds_number"""
        mach = 0.8
        alt_targets = [0., 10., 20., 30., 40., 50.]
        for alt in alt_targets:
            rel_a = atm_unit_reynolds_number1a(alt*1000., mach)
            rel_b = atm_unit_reynolds_number2a(alt*1000., mach)
            assert np.allclose(rel_a, rel_b, atol=1e-3), 'rel_a=%s rel_b=%s' % (rel_a, rel_b)

            rel_c = atm_unit_reynolds_number1b(alt, mach, alt_units='kft', reynolds_units='1/ft')
            rel_d = atm_unit_reynolds_number2b(alt, mach, alt_units='kft', reynolds_units='1/ft')
            assert np.allclose(rel_c, rel_d, atol=1e-3), 'rel_c=%s rel_d=%s' % (rel_c, rel_d)

            assert np.allclose(rel_a, rel_d, atol=1e-3), 'rel_a=%s rel_d=%s' % (rel_a, rel_d)

    def test_sweep(self):
        """tests FLFACT sweeps"""
        mach = 0.8
        alts = np.linspace(-10000, 80000.)
        rho, mach, vel = make_flfacts_alt_sweep(
            mach, alts, eas_limit=300., alt_units='m',
            velocity_units='m/s',
            density_units='kg/m^3',
            eas_units='m/s')
        del rho, mach, vel, alts

        alt = 10000.
        machs = np.linspace(0., 0.8)
        rho, mach, vel = make_flfacts_mach_sweep(
            alt, machs, eas_limit=300.,
            alt_units='m',
            velocity_units='m/s',
            density_units='kg/m^3',
            eas_units='m/s')
        del rho, mach, vel, alt

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
