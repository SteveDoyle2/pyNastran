"""
    bar3_wrap_f.py - Bar3 (Fortran implementation) for the a three bar truss
    example structures problem. This openMDAO component contains a three bar
    truss example referenced in CometBoards
"""
import os

# pylint: disable-msg=E0611,F0401
from openmdao.main.api import Assembly
from openmdao.lib.datatypes.api import Float

from bar3_static_nastran import Bar3Static
from bar3_dynamic_nastran import Bar3Dynamic

class Bar3Truss(Assembly):
    """ Model of a three bar truss - Fortran Implementation."""

    # set up interface to the framework
    # pylint: disable-msg=E1101

    def __init__(self):

        super(Bar3Truss, self).__init__()

        self.add('bar3_static', Bar3Static())
        self.add('bar3_dynamic', Bar3Dynamic())

        nastran_command = "/msc/nastran/bin/nastran"

        self.bar3_static.stdout = os.devnull
        self.bar3_static.stderr = os.devnull
        self.bar3_static.nastran_filename = "vared_bar3.bdf"
        self.bar3_static.nastran_command = nastran_command

        self.bar3_dynamic.stdout = os.devnull
        self.bar3_dynamic.stderr = os.devnull
        self.bar3_dynamic.nastran_filename = "bar3_dyn.bdf"
        self.bar3_dynamic.nastran_command = nastran_command

        self.create_passthrough('bar3_static.bar1_area')
        self.create_passthrough('bar3_static.bar2_area')
        self.create_passthrough('bar3_static.bar3_area')

        self.create_passthrough('bar3_static.load_x_dir')
        self.create_passthrough('bar3_static.load_y_dir')
        self.create_passthrough('bar3_static.loadmag')

        self.create_passthrough('bar3_static.Youngs_Modulus')
        self.create_passthrough('bar3_static.weight_density')

        self.create_passthrough('bar3_static.bar1_stress')
        self.create_passthrough('bar3_static.bar2_stress')
        self.create_passthrough('bar3_static.bar3_stress')

        self.create_passthrough('bar3_static.displacement_x_dir')
        self.create_passthrough('bar3_static.displacement_y_dir')

        self.create_passthrough('bar3_static.weight')

        self.create_passthrough('bar3_dynamic.lumpedmass')
        self.create_passthrough('bar3_dynamic.frequency')

        self.connect('bar1_area', 'bar3_dynamic.bar1_area')
        self.connect('bar2_area', 'bar3_dynamic.bar2_area')
        self.connect('bar3_area', 'bar3_dynamic.bar3_area')

        self.connect('Youngs_Modulus', 'bar3_dynamic.Youngs_Modulus')
        self.connect('weight_density', 'bar3_dynamic.weight_density')


#end Bar3Truss

if __name__ == "__main__": # pragma: no cover

    truss = Bar3Truss()

    import time
    time1 = time.time()

    truss.run()

    print " "
    print "Weight = %8.4f" % (truss.weight)
    print " "
    print "Bar Forces = %8.4f" % (truss.bar1_force), \
                       "%8.4f" % (truss.bar2_force), \
                       "%8.4f" % (truss.bar3_force)
    print " "
    print "Bar Stresses = %8.4f" % (truss.bar1_stress), \
                         "%8.4f" % (truss.bar2_stress), \
                         "%8.4f" % (truss.bar3_stress)
    print " "
    print "Displacement in x & y directions = %8.4f" \
                   % (truss.displacement_x_dir),"%8.4f" \
                   % (truss.displacement_y_dir)
    print " "
    print "Frequency = %8.4f" % (truss.frequency)
    print " "
    print "Elapsed time: ", time.time()-time1
