"""
    comp_plate_static_nastran.py - Composite plate implementation example.
"""

from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent

class Comp_Plate(NastranComponent):
    """ Model of a composite plate """

    #initial thickness of each ply(sum of all the layers)

    thick1_init = 2.6562
    thick2_init = 2.6479
    thick3_init = 1.3309

    # design variables
    thick1 = Float(thick1_init, iotype="in", units="inch", desc="Thickness of pcomp 801")
    thick2 = Float(thick2_init, iotype="in", units="inch", desc="Thickness of pcomp 802")
    thick3 = Float(thick3_init, iotype="in", units="inch", desc="Thickness of pcomp 803")

    # outputs
    property1_max_major_strain = Float(0.0, iotype="out", desc="max major strain for pcomp 801")
    property2_max_major_strain = Float(0.0, iotype="out", desc="max major strain for pcomp 802")
    property3_max_major_strain = Float(0.0, iotype="out", desc="max major strain for pcomp 803")
    property1_max_minor_strain = Float(0.0, iotype="out", desc="max minor strain for pcomp 801")
    property2_max_minor_strain = Float(0.0, iotype="out", desc="max minor strain for pcomp 802")
    property3_max_minor_strain = Float(0.0, iotype="out", desc="max minor strain for pcomp 803")

    property1_max_major_minor_strain = Float(0.0, iotype="out", desc="max major minor strain for pcomp 801")
    property2_max_major_minor_strain = Float(0.0, iotype="out", desc="max major minor strain for pcomp 802")
    property3_max_major_minor_strain = Float(0.0, iotype="out", desc="max major minor strain for pcomp 803")

    displacement_18_z_dir = Float(1.25, iotype='out',
                               units='inch',
                               desc='Displacement in z-direction',
                               nastran_header="displacements",
                               nastran_subcase=1,
                               nastran_time_step_freq_mode=None,
                               nastran_constraints={"translations" : 18},
                               nastran_dof=2)  # 0-based

    def mass(op2):
        return op2.grid_point_weight.mass[0]

    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')

    def execute(self):

        super(Comp_Plate, self).execute()

        self.comp_elm_dict = {}
        for cquad4 in range(1,26):
            elmid = self.bdf.elements[cquad4].eid
            pid = self.bdf.elements[cquad4].pid.pid
            if pid not in self.comp_elm_dict:
                self.comp_elm_dict[ pid ] = []

            self.comp_elm_dict[pid].append( elmid )

        max_minor_strain_by_pid, max_major_strain_by_pid = self.calculate_max_strains()

        self.property1_max_major_strain = max_major_strain_by_pid[ 801 ]
        self.property2_max_major_strain = max_major_strain_by_pid[ 802 ]
        self.property3_max_major_strain = max_major_strain_by_pid[ 803 ]

        self.property1_max_minor_strain = max_minor_strain_by_pid[ 801 ]
        self.property2_max_minor_strain = max_minor_strain_by_pid[ 802 ]
        self.property3_max_minor_strain = max_minor_strain_by_pid[ 803 ]

        # Calculate the maximum strain (max(major,minor)) for each property 
        self.property1_max_major_minor_strain = max( self.property1_max_major_strain,
                                                     self.property1_max_minor_strain )
        self.property2_max_major_minor_strain = max( self.property2_max_major_strain,
                                                     self.property2_max_minor_strain )
        self.property3_max_major_minor_strain = max( self.property3_max_major_strain,
                                                     self.property3_max_minor_strain )

    def calculate_max_strains( self ):
        '''Using the data from the input.out file,
            calculate the max major and minor strains
            for each of the three properties'''

        # find the max major and minor strains for each element ID
        max_major_strain_by_elmid = {}
        max_minor_strain_by_elmid = {}

        isubcase = 1
        for elmid in self.op2.compositePlateStrain[isubcase].majorP :

            # self.op2.compositePlateStrain[1].majorP[1] returns a list like this
            #  [7.347079372266307e-05, 2.4492490410921164e-05, 0.002068278146907687, 0.006204872392117977]
            max_major_strain_by_elmid[ elmid ] = max( self.op2.compositePlateStrain[isubcase].majorP[elmid] )
            max_minor_strain_by_elmid[ elmid ] = max( self.op2.compositePlateStrain[isubcase].minorP[elmid] )

        # Find the max minor and major strains for each property
        max_major_strain_by_pid = {}
        max_minor_strain_by_pid = {}
        # Go through the dictionary of self.comp_elm_dict
        for pid in self.comp_elm_dict:
            element_ids_in_pid = self.comp_elm_dict[ pid ]

            max_major_strain_by_pid[ pid ] = max( [ max_major_strain_by_elmid[ elmid ] for elmid in element_ids_in_pid ] ) 
            max_minor_strain_by_pid[ pid ] = max( [ max_minor_strain_by_elmid[ elmid ] for elmid in element_ids_in_pid ] ) 

        return max_minor_strain_by_pid, max_major_strain_by_pid
        

    def update_hook(self):

        # We want to keep the ratios of ply thickness
        # of each ply, but we want to change them in relation
        # to the overall thickness (t1, t2, and t3)

        # We'll use NastranMaker to set the individual ply's
        # thicknesses.

        value_distribution = {"801" : [.25,.25,.25,.25],
                              "802" : [.25,.25,.25,.25],
                              "803" : [.25,.25,.25,.25]}

        
        # for each pcomp, we have to set all four plys
        for pcomp in range(1,4): # there are three pcomps
            pid = str(800 + pcomp)
            values = []
            for x in value_distribution[pid]:
                values.append(x * self.__getattribute__("thick%d" % pcomp) )

            for ply in range(4): # there are four plys
                pcomp = self.bdf.Property( int(pid) , "dummy message") 
                pcomp.plies[ply][1] = values[ply]
