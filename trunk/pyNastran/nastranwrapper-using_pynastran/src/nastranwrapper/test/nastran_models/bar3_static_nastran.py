from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent
from nastranwrapper.test.nastranwrapper_test_utils import calculate_stress

class Bar3Static(NastranComponent):
    """ Model of a three bar truss - Fortran Implementation."""

    bar1_area  = Float(1., nastran_card="PROD",
                       nastran_id="11",
                       nastran_field='A',
                       low=0.0009, high=10000.,
                       iotype='in', units='inch*inch',
                       desc='Cross-sectional area for bar 1')

    bar2_area  = Float(1., nastran_card="PROD",
                       nastran_id="12",
                       nastran_field='A',
                       low=0.0009, high=10000.,
                       iotype='in', units='inch*inch',
                       desc='Cross-sectional area for bar 2')

    bar3_area  = Float(1., nastran_card='PROD',
                       nastran_id="13",
                       nastran_field='A',
                       low=0.0009, high=10000.,
                        iotype='in', units='inch*inch',
                        desc='Cross-sectional area for bar 3')

    bar1_stress = Float(0., 
                        iotype='out',
                        units='lb/(inch*inch)',
                        desc='Stress in bar 1')
    bar2_stress = Float(0., 
                        iotype='out',
                        units='lb/(inch*inch)',
                        desc='Stress in bar 2')
    bar3_stress = Float(0., 
                        iotype='out',
                        units='lb/(inch*inch)',
                        desc='Stress in bar 3')

    # def disp(op2,**keywords):
    #     d = op2.displacements[keywords['isubcase']].translations[keywords['id']][keywords['xyz']]
    #     return d

    displacement_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 1},
                               nastran_dof=0)  # 0-based
                                                                                                                            
    displacement_y_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 1},
                               nastran_dof=1)  # 0-based
                                                                                                                            
    def mass(op2):
        return op2.grid_point_weight.mass[0]


    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')

    def execute(self):
        """ Simulates the analysis of a three bar truss structure.
            Force, Stress, Displacement,Frequency and Weight are returned at
            the Bar3Truss output.
        """

        super(Bar3Static, self).execute()

        stresses = []
        for i in range(1,4):
            isubcase = 1
            axial = self.op2.rodStress[isubcase].axial[ i ]
            torsion = self.op2.rodStress[isubcase].torsion[ i  ]
            stresses.append((axial, torsion))

        [self.bar1_stress, self.bar2_stress, self.bar3_stress] = \
                          map(calculate_stress, stresses)


