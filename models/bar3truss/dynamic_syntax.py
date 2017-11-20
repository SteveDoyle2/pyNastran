"""
Demonstrates reading the OpenMDAO syntax
"""
from __future__ import print_function
import os
from numpy import sqrt, searchsorted

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

try:
    from pyNastran.f06.f06 import F06
    is_f06 = True
except ImportError:
    is_f06 = False
is_f06 = False

def calculate_stress(axial, torsion):
    """calculate VM stress"""
    sigma = 2 * axial * axial
    tau = 3 * torsion * torsion
    val = sqrt(.5 * (sigma + tau))
    return val

def main():
    """runs the dynamic syntax problem"""
    is_op2 = True
    variables = {
        'bar1_a' : 1.0,
        'bar2_a' : 2.0,
        'bar3_a' : 3.0,
        'loadx' : 50000.0,
        'loady' : 100000.0,
        'loadmag': 1.,
        'rho'   : 0.284,
        'youngs': 30000000.0,
    }
    model = BDF(debug=False)
    model.set_dynamic_syntax(variables)
    model.read_bdf('vared_bar3.bdf')
    out_bdf = 'out.bdf'
    out_f06 = 'out.f06'
    out_op2 = 'out.op2'

    if 'POST' in model.params:
        # change "PARAM, POST, 0"  to "PARAM, POST, -1"

        # option 1
        #model.params['POST'].update_field(2, -1)

        # option 2
        model.params['POST'].update_values(value1=-1, value2=None)

    model.write_bdf(out_bdf, size=16, is_double=False)
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf)

    if is_f06:
        model2 = F06()
        model2.read_f06(out_f06)
        print("")
        print("mass = %s" % model2.grid_point_weight.mass)
    else:
        model2 = None

    #========================================
    model4 = OP2()
    model4.read_op2(out_op2)

    eids = [2, 3]

    #========================================
    subcase1 = model4.crod_stress[1]
    combined_stress = calculate_stress(subcase1.data[0, :, 0], subcase1.data[0, :, 1])

    search = searchsorted(subcase1.element, eids)

    print("---op2 vectorized---")
    #[axial, torsion, SMa, SMt]
    for i, j in enumerate(search):
        eid = eids[i]
        axial = subcase1.data[0, j, 0]
        torsion = subcase1.data[0, j, 1]
        combined = combined_stress[i]
        print('axial   stress[%s] = %s' % (eid, axial))
        print('torsion stress[%s] = %s' % (eid, torsion))
        print('        stress[%s] = %s\n' % (eid, combined))



if __name__ == '__main__':  # pragma: no cover
    main()
