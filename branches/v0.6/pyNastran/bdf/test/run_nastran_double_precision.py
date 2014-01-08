import os
import sys

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06


if __name__ == '__main__':
    bdf_name = sys.argv[1]
    base, ext = os.path.splitext(bdf_name)

    #===========================
    # run nastran and verify the starting model is correct
    os.system('nastran scr=yes bat=no news=no %s' % bdf_name)

    f06_name = base + '.f06'
    model2 = F06(f06_name)
    model2.read_f06()

    #===========================
    # read/write the model in double precision
    out_bdf = 'out.bdf'
    model3 = BDF()
    model3.read_bdf(bdf_name)
    model3.write_bdf(out_bdf, size=16, precision='double')

    #===========================
    # run nastran again
    os.system('nastran scr=yes bat=no news=no %s' % out_bdf)
    #===========================
    # verify it's correct
    out_f06 = 'out.f06'
    model4 = F06(out_f06)
    model4.read_f06()

    print('\n\npassed!!')
