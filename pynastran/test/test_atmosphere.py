import unittest

from pyNastran.utils.atmosphere import (
    atm_density, atm_dynamic_pressure, atm_temperature, atm_pressure,
    atm_velocity, atm_mach, atm_dynamic_viscosity_mu, get_alt_for_mach_q,
    atm_unit_reynolds_number)
from pyNastran.utils.atmosphere2 import (
    atm_density as atm_density2,
    atm_dynamic_pressure as atm_dynamic_pressure2,
    atm_temperature as atm_temperature2,
    atm_pressure as atm_pressure2,
    atm_velocity as atm_velocity2,
    atm_mach as atm_mach2,
    #atm_dynamic_viscosity_mu, get_alt_for_mach_q,
    #atm_unit_reynolds_number
)


class Atm(unittest.TestCase):
    def test_temperature(self):
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
        self.assertEqual(atm_pressure(alt=10 *1000.), 1456.3074319943232)
        self.assertEqual(atm_pressure(alt=60 *1000.), 151.20878913237249)
        self.assertEqual(atm_pressure(alt=120*1000.), 9.8627955961437763)
        self.assertEqual(atm_pressure(alt=165*1000.), 1.7725806687593277)
        self.assertEqual(atm_pressure(alt=230*1000.), 0.13023784776280109)
        self.assertEqual(atm_pressure(alt=270*1000.), 0.017353278750799964)
        self.assertEqual(atm_pressure(alt=350*1000.), 0.00028114006933161638)

        self.assertEqual(atm_pressure2(alt=10*1000., alt_units='ft'), 1456.3074319943232)
        self.assertEqual(atm_pressure2(alt=60*1000., alt_units='ft'), 151.20878913237249)
        self.assertEqual(atm_pressure2(alt=120*1000., alt_units='ft'), 9.8627955961437763)
        self.assertEqual(atm_pressure2(alt=165*1000., alt_units='ft'), 1.7725806687593277)
        self.assertEqual(atm_pressure2(alt=230*1000., alt_units='ft'), 0.13023784776280109)
        self.assertEqual(atm_pressure2(alt=270*1000., alt_units='ft'), 0.017353278750799964)
        self.assertEqual(atm_pressure2(alt=350*1000., alt_units='ft'), 0.00028114006933161638)

        self.assertEqual(atm_pressure2(alt=10, alt_units='kft'), 1456.3074319943232)
        self.assertEqual(atm_pressure2(alt=60, alt_units='kft'), 151.20878913237249)
        self.assertEqual(atm_pressure2(alt=120, alt_units='kft'), 9.8627955961437763)
        self.assertEqual(atm_pressure2(alt=165, alt_units='kft'), 1.7725806687593277)
        self.assertEqual(atm_pressure2(alt=230, alt_units='kft'), 0.13023784776280109)
        self.assertEqual(atm_pressure2(alt=270, alt_units='kft'), 0.017353278750799964)
        self.assertEqual(atm_pressure2(alt=350, alt_units='kft'), 0.00028114006933161638)

    def test_viscosity(self):
        self.assertEqual(atm_dynamic_viscosity_mu(alt=10 *1000.), 3.5317481186391660e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=60 *1000.), 2.9702384755729678e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=120*1000.), 3.3485025784164385e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=165*1000.), 3.6827603483595828e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=230*1000.), 2.9969664628998927e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=270*1000.), 2.7383347922674784e-07)
        self.assertEqual(atm_dynamic_viscosity_mu(alt=350*1000.), 2.7383347922674784e-07)

    def test_density(self):
        self.assertEqual(atm_density(299.), 0.0023600367621627013)
        self.assertEqual(atm_density2(299., alt_units='ft'), 0.0023600367621627013)
        self.assertEqual(atm_density2(299./1000., alt_units='kft'), 0.0023600367621627013)

    def test_velocity(self):
        vel_a = atm_velocity(55000., 2.4)
        vel_b = atm_velocity(55000., 2.4)
        vel_c = atm_velocity(55000., 2.0)

        self.assertAlmostEqual(vel_a, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_b, 2323.0, delta=0.1)
        self.assertAlmostEqual(vel_a, vel_b)
        self.assertAlmostEqual(vel_c, 1935.8, delta=0.1)

        self.assertRaises(TypeError, atm_velocity, 55000, [2.4])
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

        mach_a = atm_mach2(14., 743.011549709834, alt_units='kft', velocity_units='ft/s', debug=False)
        mach_b = atm_mach2(14., 2122.8901420281, alt_units='kft', velocity_units='ft/s')
        mach_c = atm_mach2(14., 628.3419, alt_units='kft', velocity_units='knots', debug=True)

        self.assertAlmostEqual(mach_a, 0.7, delta=0.01)
        self.assertAlmostEqual(mach_b, 2.0, delta=0.002)
        self.assertAlmostEqual(mach_c, 1.0, delta=0.002)

    def test_get_q(self):
        get_alt_for_mach_q(0.8, 534.2)
        alt_a = get_alt_for_mach_q(0.8, 600., tol=0.01)
        alt_b = get_alt_for_mach_q(0.8, 400., tol=0.01)

        self.assertAlmostEqual(alt_a, 12144.30, delta=1.5)  # TODO: should be 0.01
        self.assertAlmostEqual(alt_b, 22058.47, delta=2.5)  # TODO: should be 0.01

    def test_re(self):
        re = atm_unit_reynolds_number(55000., 2.4)

        self.assertEqual(atm_unit_reynolds_number(55000., 2.4), 2244166.3810534105)

if __name__ == "__main__":  # pragma: no cover
    unittest.main()
