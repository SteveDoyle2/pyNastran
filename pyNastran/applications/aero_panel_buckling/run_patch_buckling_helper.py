import os
from copy import deepcopy
from six import iteritems, string_types
import glob
import subprocess
import time

import matplotlib.pyplot as plt
import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

def run_nastran(fname, keywords=None):
    """
    Call a nastran subprocess with the given filename

    Parameters
    -----------
    fname : string
        Filename of the Nastran .bdf file
    keywords : list of strings, optional
        Default keywords are `'mem=1024mb'`, `'old=no'`, and `'news=no'`
    """
    if keywords is None:
        keywords = ['old=no','news=no'] # ['mem=1024mb','old=no','news=no']
    return subprocess.call(['nastran', fname] + keywords)

def run_bdfs_batch(bdf_filenames, workpath='results', mem='100mb', auth=None, overwrite_op2_if_exists=False):
    #print(bdf_filenames)
    print('Start running patch jobs.')

    run_filename = os.path.join(workpath, 'run_jobs.sh')
    run_file = open(run_filename, 'wb')
    count = 1
    if not os.path.exists('linux'):
        os.makedirs('linux')
    fnames = os.listdir('linux')
    for fname in fnames:
        os.remove(os.path.join('linux', fname))

    import shutil
    for bdf_filename in bdf_filenames:
        basename = os.path.basename(bdf_filename)
        patch_id_str = basename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        op2_filename = 'patch_%s.op2' % patch_id
        if not os.path.exists(op2_filename) or overwrite_op2_if_exists:
            #print(bdf_filename)
            #shutil.copyfile(bdf_filename, os.path.join('linux', basename))
            #print('working on %s' % bdf_filename)
            #cmd = 'nastran {} scr=yes bat=no mem=100MB old=no'.format(bdf_filename)  # subprocess/os.system version
            #cmd = 'nastran results/%s scr=yes' % (basename) # shell version
            cmd = 'nastran %s scr=yes bat=no old=no' % (basename) # shell version
            for key, value in nastran_keywords:
                cmd += ' %s=%s' % (key, value)

            #subprocess.call(['nastran', bdf_filename, 'scr=yes', 'bat=no', 'old=no'])
            run_file.write('echo "%s"\n' % count)
            run_file.write('%s\n' % cmd)
            #run_file.write('sleep 5\n')
            #print('cmd = %r' % cmd)
            #run_nastran(bdf_filename)
            #os.system(cmd)
            #time.sleep(15)
            count += 1
    run_file.close()
    #shutil.copyfile('run_jobs.sh', os.path.join(work))
    #os.system('bash run_jobs.sh')
    print('Done running patch jobs.')


def run_bdfs(bdf_filenames, workpath='results', nastran_keywords=None, overwrite_op2_if_exists=False):
    #print(bdf_filenames)
    print('Start running patch jobs.')
    for bdf_filename in bdf_filenames:
        basename = os.path.basename(bdf_filename)
        patch_id_str = basename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        op2_filename = 'patch_%s.op2' % patch_id
        if not os.path.exists(op2_filename) or overwrite_op2_if_exists:
            #cmd = 'nastran %s scr=yes bat=no old=no' % (basename) # shell version

            call_args = ['nastran', bdf_filename, 'scr=yes', 'bat=no', 'old=no']
            for key, value in nastran_keywords:
                call_args.append('%s=%s' % (key, value))

            subprocess.call(call_args)
    print('Done running patch jobs.')


def load_sym_regions_map(sym_regions_filename):
    print(os.getcwd())
    with open(sym_regions_filename, 'r') as sym_regions_file:
        lines = sym_regions_file.readlines()

    region_to_symregion_map = {}
    symregion_to_region_map = {}
    for line in lines[1:]:
        sline = line.strip().split(',')
        values = [int(val) for val in sline]
        region_id = values[0]
        sym_region_id = values[1]
        region_to_symregion_map[region_id] = sym_region_id
        symregion_to_region_map[sym_region_id] = region_id
    return region_to_symregion_map, symregion_to_region_map


def get_eigs():
    patch_files = glob.glob('patch_*.bdf')
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
        print('Patch:', pn, 'Min eigenvalue: ', min_eval)
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
    bdf_filenames = glob.glob('%s/patches/patch_*.bdf' % workpath)
    run_bdfs(bdf_filenames)
    #get_eigs()


if __name__ == '__main__':
    main()
