"""
    dome_static_nastran.py - implementation for the geodesic dome
    example structures problem referenced in CometBoards
"""
from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent
from nastranwrapper.test.nastranwrapper_test_utils import calculate_stress

class DomeStatic(NastranComponent):

    for i in range(1,157):
        cmd = "bar%d_init_area = 1.00" %i
        exec (cmd)

    for i in range(157,253):
        cmd = "tria%d_init_thickness = 1.00" %i
        exec (cmd)

    for i in range(1,157):
        cmd = 'bar%d_area  = Float(bar%d_init_area, nastran_card="PROD",\
                       nastran_id=%d, \
                       nastran_field="A",\
                       iotype="in", units="inch*inch",\
                       desc="Cross-sectional area for bar %d")' %(i,i,i,i)
        exec(cmd)

    for i in range(157,253):
        cmd = 'tria%d_thickness  = Float(tria%d_init_thickness, nastran_card="PSHELL",\
                       nastran_id=%d, \
                       nastran_field="T",\
                       iotype="in", units="inch*inch",\
                       desc="Membrane thickness for tria %d")' %(i,i,i,i)
        exec(cmd)

    # these are stresses that will be  constrained
    for i in range(1,157):
        cmd = "bar%d_stress = Float(0., iotype='out', units='lb/(inch*inch)', desc='Axial stress in element %d')" %(i,i)
        exec(cmd)
    for i in range(157,253):
        cmd = "tria%d_stress = Float(0., iotype='out', units='lb/(inch*inch)', desc='Von Mises stress in element %d')" %(i,i)
        exec(cmd)

    def mass(op2):
        return op2.grid_point_weight.mass[0]

    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')


    def execute(self):
        """ Simulates the analysis of a ten bar truss structure.
            Force, Stress, Displacement,Frequency and Weight are returned at
            the Ring output.
        """

        super(DomeStatic, self).execute()

        # get stresses from table with header
        #   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        isubcase = 1
        for i in range( len( self.op2.rodStress[isubcase].axial ) ) :
            stress = calculate_stress( ( self.op2.rodStress[isubcase].axial[ i + 1 ],
                                       self.op2.rodStress[isubcase].torsion[ i + 1 ] ) )
            cmd = "self.bar%d_stress = stress" % (i+1)
            exec(cmd)

        # get stresses from table with header
        #   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )"
        for i in range(157,253):
            biggest = self.op2.plateStress[isubcase].ovmShear[i]['CEN/3'][0]
            cmd = "self.tria%d_stress = biggest" %i
            exec(cmd)

# def group_von_mises(groups, von_mises):
#     final = []
#     for group in groups:
#         final.append([])
#         for element in group:
#             final[-1].append(abs(von_mises[element-1])) # stresses is zero indexed
#         # we actually just wanted the maximum
#         final[-1] = max(final[-1])
#     return final
# 
