from __future__ import print_function
import glob
import os
from six import iteritems
from collections import defaultdict
from numpy import where, unique, array, zeros, searchsorted, log10, array_equal

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2, FatalError
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import load_sym_regions_map


def load_regions(regions_filename):
    with open(regions_filename, 'r') as regions_file:
        lines = regions_file.readlines()

    regions_to_pid_map = {}
    regions_to_eids_map = {}
    for line in lines[1:]:
        sline = line.strip().split(',')
        values = [int(val) for val in sline]
        pid = values[0]
        regions_patch_id = values[1]
        eids = values[2:]
        regions_to_eids_map[regions_patch_id] = eids
        regions_to_pid_map[regions_patch_id] = pid
    return regions_to_pid_map, regions_to_eids_map


def load_regions_and_create_eigenvalue_csv(bdf_filename, op2_filenames, regions_filename, sym_regions_filename=None):
    #patch_numbers = []
    #evals = []

    min_eigenvalue_by_patch_id = {}
    is_sym_regions = False
    if sym_regions_filename is not None:
        is_sym_regions = True
        region_to_symregion_map, symregion_to_region_map = load_sym_regions_map(sym_regions_filename)
    msg = ''
    for op2_filename in op2_filenames:
        if not os.path.exists(op2_filename):
            print(op2_filename)
            continue
        patch_id_str = op2_filename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)

        # = pf.split('_')[1].split('.')[0]
        #patch_numbers.append(patch_id)
        #model = BDF()
        #model.read_bdf(pf)
        # eids = model.elements.keys()
        model2 = OP2(debug=False)
        #op2_path = '%s_.op2' % patch_id)
        try:
            model2.read_op2(op2_filename)
        except FatalError:
            #os.remove(op2_filename)
            msg += '%s\n' % op2_filename
            continue
        cases = model2.eigenvectors.keys()
        isubcase = cases[0]
        d = model2.eigenvectors[isubcase]
        eigrs = array(d.eigrs)
        #print('eigrs =', eigrs, type(eigrs))

        i = where(eigrs > 0.0)[0]
        if len(i) == 0:
            min_eigenvalue = 0.  # TODO: no buckling eigenvalue...wat?
        else:
            min_eigenvalue = log10(eigrs[i].min())

        #evals.append(min_eval)
        print('Patch:', patch_id, 'Min eigenvalue: ', min_eigenvalue)
        min_eigenvalue_by_patch_id[patch_id] = min_eigenvalue
        if is_sym_regions:
            if patch_id in symregion_to_region_map:
                sym_patch_id = symregion_to_region_map[patch_id]
            elif patch_id in region_to_symregion_map:
                sym_patch_id = region_to_symregion_map[patch_id]
            else:
                asdf
            min_eigenvalue_by_patch_id[sym_patch_id] = min_eigenvalue
    print(msg)
    model = BDF()
    model.read_bdf(bdf_filename)
    all_eids = unique(model.elements.keys())
    neids = len(all_eids)

    eigenvalues = zeros(neids, dtype='float32')
    with open(regions_filename, 'r') as regions_file:
        lines = regions_file.readlines()

        for iline, line in enumerate(lines[1:]):
            sline = line.strip().split(',')
            #print(sline)
            values = [int(val) for val in sline]
            pid = values[0]
            regions_patch_id = values[1]
            eids = values[2:]
            i  = searchsorted(all_eids, eids) # [0] ???
            assert array_equal(all_eids[i], eids), 'iline=%s pid=%s patch_id=%s' % (iline, pid, regions_patch_id)
            if regions_patch_id not in min_eigenvalue_by_patch_id:
                print('missing pid=%s' % pid)
                continue
            #patch_id[i] = panel_id
            eigenvalues[i] = min_eigenvalue_by_patch_id[regions_patch_id]

    eigenvalue_filename = 'eigenvalues_output.csv'
    with open(eigenvalue_filename, 'w') as eigenvalue_file:
        eigenvalue_file.write('# log(Eigenvalue), eigenvalue, is_buckled\n')
        for log10_eig in eigenvalues:
            if log10_eig < 0:
                is_buckled = 1.
            else:
                is_buckled = 0
            eig = 10 ** log10_eig
            eig = min([eig, 1.1])
            eigenvalue_file.write('%f, %f, %i\n' % (log10_eig, eig, is_buckled))


def split_model_by_pid_panel(workpath='results'):
    patch_filenames = glob.glob('%s/patches/patch_*.bdf' % workpath)

    pid_panel = defaultdict(list)
    for patch_filename in patch_filenames:
        #print(patch_filename)
        basename = os.path.basename(patch_filename)
        sline = basename.split('.')[0].split('_')[1]
        #print('sline', sline)
        ipanel = int(sline)
        #print('ipanel = %s' % ipanel)

        model = BDF()
        model.read_bdf(patch_filename, xref=False)

        eids = defaultdict(list)
        for eid, elem in iteritems(model.elements):
            pid = elem.pid
            key = (pid, ipanel)
            pid_panel[key].append(eid)
            eids[pid].append(eid)

        if 0:
            for pid, eidsi in sorted(iteritems(eids)):
                pid_filename = os.path.join(workpath, 'pid_%i_ipanel_%i.txt' % (pid, ipanel))
                out = str(eidsi)
                with open(pid_filename, 'w') as pid_file:
                    pid_file.write('# PSHELL pid\n')
                    pid_file.write('# eids\n')
                    pid_file.write('%s\n' % pid)
                    pid_file.write('%s\n' % out[1:-1])

    regions_filename = 'regions.txt'
    nregions = 0
    with open(regions_filename, 'w') as regions_file:
        regions_file.write('# pid, ipanel, eids\n')
        for key, eidsi in iteritems(pid_panel):
            pid, ipanel = key
            out = str(eidsi)

            regions_file.write('%s, %s, ' % (pid, ipanel))
            regions_file.write('%s\n' % out[1:-1])
            nregions += 1
    assert nregions > 0, nregions
    # done



def main():
    """
    key=AELINK   value=8
    key=AELIST   value=4
    key=AERO     value=1
    key=AEROS    value=1
    key=AESTAT   value=11
    key=AESURF   value=4
    key=CAERO1   value=6
    key=CBAR     value=8
    key=CONM2    value=1
    key=CORD2R   value=4
    key=CTRIA3   value=120338
    key=DMI      value=2
    key=EIGRL    value=1
    key=ENDDATA  value=1
    key=GRID     value=55541
    key=MAT1     value=2
    key=PAERO1   value=1
    key=PARAM    value=2
    key=PBAR     value=1
    key=PSHELL   value=117
    key=RBE2     value=4
    key=SET1     value=6
    key=SPLINE1  value=6
    key=SUPORT   value=1
    key=TRIM     value=4
    """
    workpath = 'results'
    #split_model_by_pid_panel(workpath)

    sym_regions_filename = 'sym_regions_map.csv'
    load_regions_and_create_eigenvalue_csv(bdf_filename, op2_filenames, 'regions.txt', sym_regions_filename=sym_regions_filename)

if __name__ == '__main__':
    main()
