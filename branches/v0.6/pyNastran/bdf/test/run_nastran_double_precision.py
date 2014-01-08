import os
import sys

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06


if __name__ == '__main__':
    model = BDF()
    model.read_bdf(sys.argv[1])

    out_bdf = 'out.bdf'
    model.write_bdf(out_bdf, size=16, precision='double')

    os.system('nastran scr=yes bat=no news=no %s' % out_bdf)
    out_f06 = 'out.f06'
    model2 = F06(out_f06)
    model2.read_f06()
