"""
    ring_25dv_static_nastran.py - Ring implementation for the a sixty bar truss
    example structures problem. This openMDAO component contains a sixty bar
    truss example referenced in CometBoards
"""

from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent
from nastranwrapper.test.nastranwrapper_test_utils import calculate_stress

class RingTruss(NastranComponent):

    for i in range(1,26):
        cmd = "bar%d_init_area = 1.0" %i
        exec (cmd)
        cmd = 'bar%d_area  = Float(bar%d_init_area, nastran_card="PROD",\
                       nastran_id="%d", \
                       nastran_field="A",\
                       iotype="in", units="inch*inch",\
                       desc="Cross-sectional area for bar %d")' %(i,i,i,i)
        exec(cmd)

    # these are stresses that will be  constrained
    for i in range(1,61):
        cmd = "bar%d_stress = Float(0., iotype='out', units='lb/(inch*inch)', desc='Axial stress in element %d')" %(i,i)
        exec(cmd)

    # these are displacements that will be  constrained
    displacement1_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 10},
                               nastran_dof=0)  # 0-based

    displacement1_y_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 10},
                               nastran_dof=1)  # 0-based

    displacement2_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :4},
                               nastran_dof=0)  # 0-based

    displacement2_y_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :4},
                               nastran_dof=1)  # 0-based


    displacement3_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :19},
                               nastran_dof=0)  # 0-based

    displacement3_y_dir = Float(2.75, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :19},
                               nastran_dof=1)  # 0-based

    displacement4_x_dir = Float(2.25, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :13},
                               nastran_dof=0)  # 0-based

    displacement4_y_dir = Float(2.25, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" :13},
                               nastran_dof=1)  # 0-based

    def mass(op2):
        return op2.grid_point_weight.mass[0]

    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')


    def execute(self):
        """ Simulates the analysis of a sixty bar truss structure.
            Force, Stress, Displacement,Frequency and Weight are returned at
            the Ring output.
        """
        super(RingTruss, self).execute()

        # get values from the table with this header
        #    S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        isubcase = 1
        for i in range( len( self.op2.rodStress[isubcase].axial ) ) :
            stress = calculate_stress( ( self.op2.rodStress[isubcase].axial[ i + 1 ],
                                       self.op2.rodStress[isubcase].torsion[ i + 1 ] ) )
            cmd = "self.bar%d_stress = stress" % (i+1)
            exec(cmd)

