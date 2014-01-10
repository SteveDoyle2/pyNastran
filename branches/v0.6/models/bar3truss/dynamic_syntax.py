import os
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06
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

    model.write_bdf(out_bdf, size=16, precision='double')
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf)

    model2 = F06(out_f06)
    model2.read_f06()

    subcase1 = model2.rodStress[1]

    #print '\n'.join(dir(subcase1))
    print ""

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