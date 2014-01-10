import os
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2
from math import sqrt

def calculate_stress(ax, tors):
    sigma = 2 * ax * ax
    tau = 3 * tors * tors
    val = sqrt(.5 * (sigma + tau))
    return val

def main():
    vars = {
        'bar1_a' : 1.0,
        'bar2_a' : 2.0,
        'bar3_a' : 3.0,
        'loadx' : 50000.0,
        'loady' : 100000.0,
        'loadmag': 1.,
        'rho'   : 0.284,
        'youngs': 30000000.0,

    }
    model = BDF()
    model.set_dynamic_syntax(vars)
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

    model.write_bdf(out_bdf, size=16, precision='double')
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf)

    model2 = F06(out_f06)
    if 0:
        model2.markerMap = {
            'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R': model2._grid_point_weight_generator,
        }
        model2.markers = model2.markerMap.keys()

    if 0:
        model2.stop_after_reading_grid_point_weight(stop=True)
    model2.read_f06()

    #print '\n'.join(dir(subcase1))
    print ""
    print "mass = %s" % model2.grid_point_weight.mass

    #========================================
    model3 = OP2(out_op2)
    model3.read_op2()
    #========================================
    for form, modeli in [('f06', model2), ('op2', model3)]:
        print "---%s---" % form
        subcase1 = modeli.rodStress[1]

        eid = 2
        print 'axial   stress[%s] = %s' % (eid, subcase1.axial[eid])
        print 'torsion stress[%s] = %s' % (eid, subcase1.torsion[eid])
        print '        stress[%s] = %s\n' % (eid, calculate_stress(subcase1.axial[eid], subcase1.torsion[eid]))

        eid = 3
        print 'axial   stress[%s] = %s' % (eid, subcase1.axial[eid])
        print 'torsion stress[%s] = %s' % (eid, subcase1.torsion[eid])
        print '        stress[%s] = %s\n' % (eid, calculate_stress(subcase1.axial[eid], subcase1.torsion[eid]))

if __name__ == '__main__':
    main()