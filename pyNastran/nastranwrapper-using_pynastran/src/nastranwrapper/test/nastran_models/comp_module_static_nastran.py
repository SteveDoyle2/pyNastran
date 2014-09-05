"""
    comp_module_static_nastran.py
"""
from openmdao.lib.datatypes.api import Float

from nastranwrapper.nastran import NastranComponent

class Comp_Module(NastranComponent):
    """ Model of a composite model """

    def mass(op2):
        return op2.grid_point_weight.mass[0]

    weight = Float(0., nastran_func=mass, iotype='out', units='lb',
                        desc='Weight of the structure')

    def execute(self):

        super(Comp_Module, self).execute()
