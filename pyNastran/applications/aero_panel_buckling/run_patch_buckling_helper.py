"""
Defines:
 - run_bdfs_batch
 - run_bdfs
 - load_sym_regions_map
 - get_eigenvalues
 - get_eigs
"""
from __future__ import print_function
import os
import sys
import glob
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import read_op2
#from pyNastran.utils.nastran_utils import run_nastran

def run_bdfs_batch(bdf_filenames, workpath='results', mem='100mb', auth=None,
                   overwrite_op2_if_exists=False):
    """
    Nastran can be very temperamental.

    Sometimes, it decides to just not run a job with os.system(...)
    or subprocess.  I suspect it's an issue of throwing too many
    models at the server.  Instead, we can create a bash script,
    run that on Linux (even if we have to switch operating systems),
    manually run the the jobs, bring them back, and continue the
    analysis.
    """
    #print(bdf_filenames)
    print('Start running patch jobs.')
    nastran_keywords = {
        'mem' : mem,
        'auth' : auth,
    }

    run_filename = os.path.join(workpath, 'run_jobs.sh')
    run_file = open(run_filename, 'wb')
    count = 1
    if not os.path.exists('linux'):
        os.makedirs('linux')
    fnames = os.listdir('linux')
    for fname in fnames:
        os.remove(os.path.join('linux', fname))

    for bdf_filename in bdf_filenames:
        basename = os.path.basename(bdf_filename)
        patch_id_str = basename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        op2_filename = 'patch_%s.op2' % patch_id
        if not os.path.exists(op2_filename) or overwrite_op2_if_exists:
            #print(bdf_filename)
            #shutil.copyfile(bdf_filename, os.path.join('linux', basename))
            #print('working on %s' % bdf_filename)

            # subprocess/os.system version
            #cmd = 'nastran {} scr=yes bat=no mem=100MB old=no'.format(bdf_filename)
            #cmd = 'nastran results/%s scr=yes' % (basename) # shell version
            cmd = 'nastran %s scr=yes bat=no old=no' % (basename) # shell version
            for key, value in nastran_keywords:
                if value is None:
                    continue
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


def run_bdfs(bdf_filenames, workpath='results', nastran_keywords=None,
             overwrite_op2_if_exists=False):
    """runs a series of BDF files using subprocess"""
    assert os.path.exists(workpath), print_bad_path(workpath)

    curdir = os.getcwd()
    print('switching from %r to %r' % (curdir, workpath))
    os.chdir(workpath)

    #print(bdf_filenames)
    print('Start running patch jobs.')
    sys.stdout.flush()
    op2_filenames = []
    assert len(bdf_filenames) > 0, 'bdf_filenames=%s' % bdf_filenames
    for bdf_filename in bdf_filenames:
        basename = os.path.basename(bdf_filename)
        assert os.path.exists(bdf_filename)
        patch_id_str = basename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        op2_filename = 'patch_%s.op2' % patch_id
        if not os.path.exists(op2_filename) or overwrite_op2_if_exists:
            #cmd = 'nastran %s scr=yes bat=no old=no' % (basename) # shell version

            call_args = ['nastran', bdf_filename, 'scr=yes', 'bat=no', 'old=no']
            for key, value in nastran_keywords.items():
                if key not in ['scr', 'bat', 'old']:
                    call_args.append('%s=%s' % (key, value))
            print(call_args)
            sys.stdout.flush()
            op2_filenames.append(os.path.join(workpath, op2_filename))

            subprocess.call(call_args)
    print('Done running patch jobs.')
    print('switching from %r to %r' % (workpath, curdir))
    os.chdir(curdir)
    return op2_filenames


def load_sym_regions_map(sym_regions_filename):
    """loads the file that defines which patches are symmetric"""
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

def get_eigenvalues(op2_filename, debug=False):
    """get the buckling eigenvalues for each panel"""
    model2 = read_op2(op2_filename, debug=debug)
    cases = model2.eigenvectors.keys()
    isubcase = cases[0]
    eigenvector = model2.eigenvectors[isubcase]
    try:
        eigrs = eigenvector.eigrs
        #eigrs = d._eigrs
    except AttributeError:
        msg = '%s.object_attributes() = %s' % (
            eigenvector.class_name, str(eigenvector.object_attributes(keys_to_skip='class_name')))
        raise RuntimeError(msg)
    return eigrs

def get_eigs(debug=False):
    """get the buckling eigenvalues for each panel"""
    patch_files = glob.glob('patch_*.bdf')
    patch_numbers = []
    evals = []
    for patch_filename in patch_files:
        patch_number = patch_filename.split('_')[1].split('.')[0]
        patch_numbers.append(int(patch_number))
        #model = read_bdf(patch_filename, debug=debug)

        #eids = model.elements.keys()
        op2_filename = ''.join([patch_filename[:-4], '.op2'])
        eigenvalues = np.array(get_eigenvalues(op2_filename, debug=debug))
        i = np.where(eigenvalues > 0.0)
        min_eigenvalue = eigenvalues[i].min()
        evals.append(min_eigenvalue)
        print('Patch:%s  Min eigenvalue:%s' % (patch_number, min_eigenvalue))

    plt.figure(figsize=(12, 9))
    plt.plot(patch_numbers, evals, label='Eigenvalues')
    plt.xticks(fontsize=12, y=-0.01)
    plt.yticks(fontsize=12, x=-0.01)
    plt.title('', fontsize=16, y=1.02)
    plt.xlabel('Patch Number', fontsize=12, y=-0.01)
    plt.ylabel('Minimum Eigenvalue', fontsize=12, x=-0.01)
    # plt.legend()
    plt.show()
    plt.close()

def main():  # pragma: no cover
    """test code"""
    #if not os.path.exists('patch_0.op2'):
    workpath = ''
    bdf_filenames = glob.glob('%s/patches/patch_*.bdf' % workpath)
    run_bdfs(bdf_filenames)
    #get_eigs()


if __name__ == '__main__':  # pragma: no cover
    main()
