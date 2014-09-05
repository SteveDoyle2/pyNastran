"""
    bar25_static_nastran.py - Antanna Tower truss example structures problem. 
    This openMDAO component contains a twenty bar truss example referenced in CometBoards
"""
from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent
from nastranwrapper.test.nastranwrapper_test_utils import calculate_stress

class Bar25Static(NastranComponent):
    """ Model of a twenty bar antenna tower truss """

    for i in range(1,26):
        cmd = "bar%d_init_area = 10.0" %i
        exec (cmd)

    for i in range(1,26):
        cmd = 'bar%d_area  = Float(bar%d_init_area, nastran_card="PROD",\
                       nastran_id=%d, \
                       nastran_field="A",\
                       iotype="in", units="inch*inch",\
                       desc="Cross-sectional area for bar %d")' %(i,i,i,i)
        exec(cmd)

    # these are stresses that will be  constrained
    for i in range(1,26):
        cmd = "bar%d_stress = Float(0., iotype='out', units='lb/(inch*inch)', desc='Axial stress in element %d')" %(i,i)
        exec(cmd)

    # these are displacements that will be  constrained
    displacement1_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 1},
                               nastran_dof=0)  # 0-based

    displacement1_y_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 1},
                               nastran_dof=1)  # 0-based

    displacement2_x_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in x-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 2},
                               nastran_dof=0)  # 0-based

    displacement2_y_dir = Float(0.20, iotype='out',
                               units='inch',
                               desc='Displacement in y-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 2},
                               nastran_dof=1)  # 0-based


    def mass(op2):
        return op2.grid_point_weight.mass[0]

    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')


    def execute(self):
        """ Simulates the analysis of a ten bar truss structure.
            Force, Stress, Displacement,Frequency and Weight are returned at
            the Ring output.
        """

        super(Bar25Static, self).execute()

        # Get stresses from table with header
        #   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        isubcase = 1
        for i in range( len( self.op2.rodStress[isubcase].axial ) ) :
            stress = calculate_stress( ( self.op2.rodStress[isubcase].axial[ i + 1 ],
                                       self.op2.rodStress[isubcase].torsion[ i + 1 ] ) )
            cmd = "self.bar%d_stress = stress" % (i+1)
            exec(cmd)
 
        # groups = [range(1,2), range(2,6), range(6,10), range(10,14), \
        #           range(14,15), range(20,21), range(22,23), \
        #           range(15,16), range(18,19), range(23,24), \
        #           range(16,17), range(21,22), range(25,26), \
        #           range(17,18), range(19,20), range(24,25)] 

# def group_stresses(groups, stresses):
#     final = []
#     for group in groups:
#         final.append([])
#         for element in group:
#             final[-1].append(abs(stresses[element-1])) # stresses is zero indexed
#         # we actually just wanted the maximum
#         final[-1] = max(final[-1])
#     return final
