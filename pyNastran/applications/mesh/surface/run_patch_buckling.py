import os
from copy import deepcopy
from six import iteritems, string_types
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import time

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
import subprocess

def run_bdfs():
    bdf_filenames = glob('patch_*.bdf')
    print 'Start running patch jobs.'
    for bdf_filename in bdf_filenames:
        patch_id_str = bdf_filename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        op2_filename = 'patch_%s.op2' % patch_id
        if not os.path.exists(op2_filename):
            print('working on %s' % bdf_filename)
            #cmd = 'nastran {} scr=yes bat=no mem=100MB old=no'.format(patch_file)
            subprocess.call(['nastran', bdf_filename, 'scr=yes', 'bat=no', 'old=no'])
            time.sleep(10)
            #print('cmd = %r' % cmd)
            #os.system(cmd)
    print 'Done running patch jobs.'


def get_eigs():
    patch_files = glob('patch_*.bdf')
    patch_numbers = []
    evals = []
    for pf in patch_files:
        pn = pf.split('_')[1].split('.')[0]
        patch_numbers.append(int(pn))
        model = BDF()
        model.read_bdf(pf)
        # eids = model.elements.keys()
        model2 = OP2()
        op2_path = ''.join([pf[:-4], '.op2'])
        model2.read_op2(op2_path)
        cases = model2.eigenvectors.keys()
        isubcase = cases[0]
        d = model2.eigenvectors[isubcase]
        eigrs = d._eigrs
        min_eval = eigrs.min()
        evals.append(min_eval)
        print 'Patch:', pn, 'Min eigenvalue: ', min_eval
    fig = plt.figure(figsize=(12, 9))
    plt.plot(patch_numbers, evals, label='Eigenvalues')
    plt.xticks(fontsize=12, y=-0.01)
    plt.yticks(fontsize=12, x=-0.01)
    plt.title('', fontsize=16, y=1.02)
    plt.xlabel('Patch Number', fontsize=12, y=-0.01)
    plt.ylabel('Minimum Eigenvalue', fontsize=12, x=-0.01)
    # plt.legend()
    plt.show()
    plt.close()

def main():
    #if not os.path.exists('patch_0.op2'):
    run_bdfs()
    #get_eigs()


if __name__ == '__main__':
    main()
