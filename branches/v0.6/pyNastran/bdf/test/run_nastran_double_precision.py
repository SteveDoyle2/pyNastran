import os
import sys

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06, FatalError


def main(bdf_name, run_first_nastran=True):
    base, ext = os.path.splitext(bdf_name)

    print "len(sys.argv) =", len(sys.argv)
    #===========================
    if run_first_nastran:
        # run nastran and verify the starting model is correct
        os.system('nastran scr=yes bat=no news=no old=no %s' % bdf_name)

        f06_name = base + '.f06'
        try:
            model2 = F06(f06_name)
            model2.read_f06()
        except FatalError as e:
            print(e)
            #return
    else:
        pass
    #===========================
    # read/write the model in double precision
    out_bdf = 'out.bdf'
    model3 = BDF()
    model3.read_bdf(bdf_name)
    model3.write_bdf(out_bdf, size=16, precision='double')
    print "---wrote the bdf---"
    #===========================
    # run nastran again
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf)
    #===========================
    # verify it's correct
    out_f06 = 'out.f06'
    model4 = F06(out_f06)
    model4.read_f06()
    print('\n\npassed!!')


if __name__ == '__main__':
    bdf_name = sys.argv[1]

    run_first_nastran = False
    if len(sys.argv) == 2:
        run_first_nastran = True

    main(bdf_name, run_first_nastran=run_first_nastran)