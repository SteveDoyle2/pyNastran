from __future__ import print_function
import os
import sys

from pyNastran.f06.errors import FatalError
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import read_op2

#try:
    #from pyNastran.f06.f06 import F06
    #is_f06 = True
#except ImportError:
    #is_f06 = False

#class remove_prints(object):
    #def __init__(self):
        #pass
    #def write(self, msg):
        #if 'ntotal' in msg:
            #raise RuntimeError()
    #def flush(self):
        #pass
#sys.stdout = remove_prints()

#is_f06 = False

def update_bdf(model, post=-1):
    #cc = model.case_control_deck
    #for isubcase, subcase in sorted(cc.subcases.items()):
        #for param, values in subcase..items(params):
            #if param in ['SPCFORCES', 'STRESS', 'DISPLACEMENT', 'STRAIN', 'MPCFORCES',
                         #'GPFORCE', 'GPSTRESS', 'VELOCITY', 'ACCELERATION']:
                #print('values =', values)
    if 'POST' in sorted(model.params):
        model.params['POST'].update_values(value1=post)
    else:
        model.rejects.append(['PARAM,POST,%i' % post])

def main(bdf_filename, run_first_nastran=True, debug=True):
    base = os.path.splitext(bdf_filename)[0]

    print("len(sys.argv) = %s" % len(sys.argv))
    #===========================
    if run_first_nastran:
        # run nastran and verify the starting model is correct
        os.system('nastran scr=yes bat=no news=no old=no %s' % bdf_filename)

        #f06_filename = base + '.f06'
        op2_filename = base + '.op2'
        try:
            model2 = read_op2(op2_filename)
        except FatalError as error:
            print(error)

        #try:
            #model2 = F06()
            #model2.read_f06(f06_filename)
        #except FatalError as error:
            #print(error)
            #return
    else:
        pass
    #===========================
    # read/write the model in double precision
    out_bdf_8 = base + '_8.bdf'
    out_bdf_16 = base + '_16.bdf'
    out_bdf_16s = base + '_16s.bdf'
    model3 = read_bdf(bdf_filename, debug=debug)
    update_bdf(model3, post=-1)

    model3.write_bdf(out_bdf_8, size=8, is_double=False)
    model3.write_bdf(out_bdf_16s, size=16, is_double=False)
    model3.write_bdf(out_bdf_16, size=16, is_double=True)
    if debug:
        print("---wrote the bdf---")
    #===========================
    # run nastran again
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf_8)
    os.system('nastran scr=yes bat=no news=no old=no %s' % out_bdf_16)
    #===========================
    # verify it's correct
    #if is_f06:
        #out_f06_8 = base + '_8.f06'
        #out_f06_16 = base + '_16.f06'
        #model4 = F06(debug=False)
        #model4.read_f06(out_f06_8)

        #model5 = F06(debug=False)
        #model5.read_f06(out_f06_16)

    out_op2_8 = base + '_8.op2'
    out_op2_16 = base + '_16.op2'

    model6 = read_op2(out_op2_16)
    print('\n\npassed!!')


def cmd_line():
    bdf_filename = sys.argv[1]

    run_first_nastran = False
    if len(sys.argv) == 2:
        run_first_nastran = True
    main(bdf_filename, run_first_nastran=run_first_nastran)

if __name__ == '__main__':  # pragma: no cover
    cmd_line()
