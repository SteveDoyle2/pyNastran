
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float

from openmdao.lib.components.nastran.nastran import NastranComponent

class Bar3Dynamic(NastranComponent):
    """ Dynamic model of three bar truss """

    lumpedmass = Float(1., nastran_card="CONM2",
                 nastran_id="21", nastran_fieldnum=4,
                 iotype='in', desc="lumped mass")

    Youngs_Modulus = Float(30000000.0,
                           nastran_card="MAT1",
                           nastran_id="5",
                           nastran_fieldnum=2,
                           iotype='in',
                           units='lb/(inch*inch)',
                           desc='Youngs Modulus')

    weight_density = Float(0.284,
                           nastran_card="MAT1",
                           nastran_id="5",
                           nastran_fieldnum=5,
                           iotype='in', units='lb/(inch**3)',
                           desc='weight density of all bars')


    bar1_area  = Float(1., nastran_card="PROD",
                       nastran_id="11", nastran_fieldnum=3,
                       low=0.0009, high=10000.,
                       iotype='in', units='inch*inch',
                       desc='Cross-sectional area for bar 1')

    bar2_area  = Float(1., nastran_card="PROD",
                       nastran_id="12", nastran_fieldnum=3,
                       low=0.0009, high=10000.,
                       iotype='in', units='inch*inch',
                       desc='Cross-sectional area for bar 2')

    bar3_area  = Float(1., nastran_card='PROD',
                       nastran_id="13", nastran_fieldnum=3,
                       low=0.0009, high=10000.,
                       iotype='in', units='inch*inch',
                       desc='Cross-sectional area for bar 3')


    frequency = Float(0., nastran_header="real eigenvalues",
                      nastran_subcase=1,
                      nastran_constraints={"MODE NO." : "1"},
                      nastran_columns=["CYCLES"],
                      iotype='out', desc='frequency')
