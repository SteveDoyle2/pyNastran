import os
import sys

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06 import F06, FatalError
from pyNastran.op2.op2 import OP2

def update_bdf(model):
    cc = model.caseControlDeck
    #for isubcase, subcase in sorted(cc.subcases.iteritems()):
        #for param, values in subcase.params.iteritems():
            #if param in ['SPCFORCES', 'STRESS', 'DISPLACEMENT', 'STRAIN', 'MPCFORCES', 'GPFORCE', 'GPSTRESS', 'VELOCITY', 'ACCELERATION']:
                #print('values =', values)
    if 'POST' in sorted(model3.params):
        model.params['POST'].update_values(value1=-1)
    else:
        model.rejects.append(['PARAM,POST,-1'])

def main(bdf_name, run_first_nastran=True, debug=True):
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
    out_bdf_8 = base + '_8.bdf'
    out_bdf_16 = base + '_16.bdf'
    model3 = BDF()
    model3.read_bdf(bdf_name)
    update_bdf(model3)

    model3.write_bdf(out_bdf_8, size=8, precision='single')
    model3.write_bdf(out_bdf_16, size=16, precision='double')
    if debug:
        print "---wrote the bdf---"
    #===========================
    # run nastran again
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf_8)
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf_16)
    #===========================
    # verify it's correct
    out_f06_8 = base + '_8.f06'
    out_f06_16 = base + '_16.f06'
    out_op2_16 = base + '_16.op2'

    model4 = F06(out_f06_8, debug=False)
    model4.read_f06()

    model5 = F06(out_f06_16, debug=False)
    model5.read_f06()

    model6 = OP2(out_op2_16, debug=False)
    model6.read_op2()
    print('\n\npassed!!')


def cmd_line():
    bdf_name = sys.argv[1]

    run_first_nastran = False
    if len(sys.argv) == 2:
        run_first_nastran = True
    main(bdf_name, run_first_nastran=run_first_nastran)

if __name__ == '__main__':
    cmd_line()