import glob
import os
from six import iteritems
from collections import defaultdict
from numpy import where, unique, array, zeros, searchsorted, log10, array_equal

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

def load_regions_and_create_eigenvalue_csv(regions_filename):
    op2_filenames = glob.glob('patch_*.op2')
    #patch_numbers = []
    #evals = []

    min_eigenvalue_by_patch_id = {}
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
        model2.read_op2(op2_filename)
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
        print 'Patch:', patch_id, 'Min eigenvalue: ', min_eigenvalue
        min_eigenvalue_by_patch_id[patch_id] = min_eigenvalue

    model = BDF()
    model.read_bdf('model_144.bdf')
    all_eids = unique(model.elements.keys())
    neids = len(all_eids)

    eigenvalues = zeros(neids, dtype='float32')
    with open(regions_filename, 'r') as regions_file:
        lines = regions_file.readlines()

        for line in lines[1:]:
            sline = line.strip().split(',')
            #print(sline)
            values = [int(val) for val in sline]
            pid = values[0]
            regions_patch_id = values[1]
            eids = values[2:]
            i  = searchsorted(all_eids, eids) # [0] ???
            assert array_equal(all_eids[i], eids)
            if regions_patch_id not in min_eigenvalue_by_patch_id:
                print('missing pid=%s' % pid)
                continue
            #patch_id[i] = panel_id
            eigenvalues[i] = min_eigenvalue_by_patch_id[regions_patch_id]

    eigenvalue_filename = 'eigenvalues_output.csv'
    with open(eigenvalue_filename, 'w') as eigenvalue_file:
        eigenvalue_file.write('# log(Eigenvalue)\n')
        for eig in eigenvalues:
            eigenvalue_file.write('%f\n' % eig)

def split_model_by_pid_panel():
    patch_filenames = glob.glob('patch_*.bdf')

    pid_panel = defaultdict(list)
    for patch_filename in patch_filenames:
        #print(patch_filename)
        sline = patch_filename.split('.')[0].split('_')[1]
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
        for pid, eidsi in sorted(iteritems(eids)):
            pid_filename = 'pid_%i_ipanel_%i.txt' % (pid, ipanel)
            out = str(eidsi)
            with open(pid_filename, 'w') as pid_file:
                pid_file.write('# PSHELL pid\n')
                pid_file.write('# eids\n')
                pid_file.write('%s\n' % pid)
                pid_file.write('%s\n' % out[1:-1])

    regions_filename = 'regions.txt'
    with open(regions_filename, 'w') as regions_file:
        regions_file.write('# pid, ipanel, eids\n')
        for key, eidsi in iteritems(pid_panel):
            pid, ipanel = key
            out = str(eidsi)

            regions_file.write('%s, %s, ' % (pid, ipanel))
            regions_file.write('%s\n' % out[1:-1])
    # done



if __name__ == '__main__':
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
    #split_model_by_pid_panel()
    load_regions_and_create_eigenvalue_csv('regions.txt')
