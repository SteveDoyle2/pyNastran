"""
defines various functions for splitting a model
 - load_regions(regions_filename)
 - load_regions_and_create_eigenvalue_csv(bdf_model, op2_filenames,
                                           regions_filename, sym_regions_filename=None,
                                           eig_min=-1.0, eig_max=1.0, eig_default=3.0)
 - split_model_by_pid_panel(patch_filenames, workpath='results')

"""
from __future__ import print_function
import os
from collections import defaultdict

from six import iteritems
import numpy as np
#from numpy import where, unique, array, zeros, searchsorted, log10, array_equal

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import read_op2, FatalError
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import (
    load_sym_regions_map)


def load_regions(regions_filename):
    """
    Loads a regions.csv file
    """
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


def load_regions_and_create_eigenvalue_csv(bdf_model, op2_filenames,
                                           regions_filename, sym_regions_filename=None,
                                           eig_min=-1.0, eig_max=1.0, eig_default=3.0):
    """
    loads a BDF and a series of OP2 filenames and creates an eigenvalue buckling plot

    Parameters
    ----------
    eig_min : float
        the required minimum eigenvalue
    eig_max : float
        the required maximum eigenvalue
    eig_default : float
        the default eigenvalue for cases that do not calculate eigenvectors
        because there were no eigenvalues in the range
    regions_filename : str
        path to regions.txt file
    sym_regions_filename : str; default=None -> No symmetry
        path to sym_regions.txt file

    Returns
    -------
    min_eigenvalue_by_patch_id : dict
        key : patch_id : int
            the integer patch id
        value : eigenvalue or reserve factor
           the reserve factor eigenvalue for buckling
    eigenvalues : (n, ) float ndarray
        the minimum eigenvalues

    Creates
    -------
    eigenvalues_output.csv : file
        csv of log10(eigenvalue), eigenvalue, is_buckled
    """
    bdf_model.log.info('load_regions_and_create_eigenvalue_csv')
    #patch_numbers = []
    #evals = []
    assert isinstance(bdf_model, BDF), type(bdf_model)

    min_eigenvalue_by_patch_id = {}
    is_sym_regions = False
    if sym_regions_filename is not None:
        is_sym_regions = True
        region_to_symregion_map, symregion_to_region_map = load_sym_regions_map(
            sym_regions_filename)
    msg = ''

    assert len(op2_filenames) > 0, 'op2_filenames=%s' % op2_filenames
    bdf_model.log.info('eig_min=%s eig_max=%s' % (eig_min, eig_max))
    for op2_filename in op2_filenames:
        bdf_model.log.info('op2_filename = %r' % op2_filename)
        if not os.path.exists(op2_filename):
            bdf_model.log.warning(op2_filename)
            continue
        op2_filename_base = os.path.basename(op2_filename)
        patch_id_str = op2_filename_base.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)

        sym_patch_id = None
        if is_sym_regions:
            if patch_id in symregion_to_region_map:
                sym_patch_id = symregion_to_region_map[patch_id]
            elif patch_id in region_to_symregion_map:
                sym_patch_id = region_to_symregion_map[patch_id]
            else:
                raise RuntimeError("can this happen???")

        # = pf.split('_')[1].split('.')[0]
        #patch_numbers.append(patch_id)
        #model = BDF(debug=False)
        #model.read_bdf(pf)
        # eids = model.elements.keys()
        #op2_path = '%s_.op2' % patch_id)
        try:
            model2 = read_op2(op2_filename, combine=True, log=None,
                              debug=False, mode='msc')
        except FatalError:
            #os.remove(op2_filename)
            bdf_model.log.error('fatal on %r' % op2_filename)
            msg += '%s\n' % op2_filename
            continue
        cases = model2.eigenvectors.keys()
        if len(cases) == 0:
            #assert is_sym_regions == False, is_sym_regions
            min_eigenvalue_by_patch_id[patch_id] = eig_default
            min_eigenvalue_by_patch_id[sym_patch_id] = eig_default
            continue

        isubcase = cases[0]
        eigenvector = model2.eigenvectors[isubcase]
        eigrs = np.array(eigenvector.eigrs)
        #cycles = (eigrs * 2 * np.pi) ** 2.
        #print('eigrs =', eigrs)

        #----------------------------------
        # calculate what's basically a reserve factor (RF); margin = reserve_factor - 1
        # take the minimum of the "tension"/"compression" RFs, which are
        # compared to different allowables

        # lambda > 0
        i = np.where(eigrs >= 0.0)[0]
        if len(i) == 0:
            pos_eigenvalue = eig_default  # TODO: no buckling eigenvalue...what?
            pos_reserve_factor = eig_default
        else:
            pos_eigenvalue = eigrs[i].min()
            pos_reserve_factor = pos_eigenvalue / eig_max
            #max_eigenvalue = np.log10(eigi)

        # lambda < 0
        if 0:
            j = np.where(eigrs < 0.0)[0]
            if len(j) == 0:
                neg_eigenvalue = eig_default  # TODO: no buckling eigenvalue...what?
                neg_reserve_factor = eig_default
            else:
                neg_eigenvalue = np.abs(eigrs[j]).min()
                neg_reserve_factor = neg_eigenvalue / abs(eig_min)
                #min_eigenvalue = np.log10(eigi)
        else:
            neg_reserve_factor = 10.
            neg_eigenvalue = 10.

        #evals.append(min_eval)
        bdf_model.log.info('Patch=%s  compression (lambda > 0); lambda=%.3f RF=%.3f' % (
            patch_id, pos_eigenvalue, pos_reserve_factor))
        bdf_model.log.info('Patch=%s  tension    (lambda < 0); lambda=%.3f RF=%.3f' % (
            patch_id, neg_eigenvalue, neg_reserve_factor))
        reserve_factor = min(neg_reserve_factor, pos_reserve_factor, eig_default)
        assert reserve_factor > 0.
        min_eigenvalue_by_patch_id[patch_id] = reserve_factor
        if is_sym_regions:
            min_eigenvalue_by_patch_id[sym_patch_id] = reserve_factor
    bdf_model.log.info(msg)

    bdf_model.log.info('finished parsing eigenvalues...')
    #model = BDF(debug=False)
    #model.read_bdf(bdf_filename)
    all_eids = np.unique(bdf_model.elements.keys())
    neids = len(all_eids)

    eigenvalues = np.zeros(neids, dtype='float32')
    with open(regions_filename, 'r') as regions_file:
        lines = regions_file.readlines()

        for iline, line in enumerate(lines[1:]):
            sline = line.strip().split(',')
            #print(sline)
            values = [int(val) for val in sline]
            pid = values[0]
            regions_patch_id = values[1]
            eids = values[2:]
            i = np.searchsorted(all_eids, eids) # [0] ???
            assert np.array_equal(all_eids[i], eids), 'iline=%s pid=%s patch_id=%s' % (
                iline, pid, regions_patch_id)
            if regions_patch_id not in min_eigenvalue_by_patch_id:
                bdf_model.log.info('missing pid=%s' % pid)
                continue
            #patch_id[i] = panel_id
            eigenvalues[i] = min_eigenvalue_by_patch_id[regions_patch_id]

    eigenvalue_filename = 'eigenvalues_output.csv'
    with open(eigenvalue_filename, 'w') as eigenvalue_file:
        eigenvalue_file.write('# log(Eigenvalue), eigenvalue, is_buckled\n')
        for eig in eigenvalues:
            eig = max(eig, 0.000001)
            if eig < 1.0:
                is_buckled = 1.0
            else:
                is_buckled = 0.0
            log10_eig = np.log10(eig)
            eigenvalue_file.write('%f, %f, %i\n' % (log10_eig, eig, is_buckled))
    return min_eigenvalue_by_patch_id, eigenvalues


def split_model_by_pid_panel(patch_filenames, workpath='results'):
    """
    creates a list of element ids for each patch
    """
    pid_panel = defaultdict(list)
    for patch_filename in patch_filenames:
        basename = os.path.basename(patch_filename)
        try:
            sline = basename.split('.')[0].split('_')[1]
        except:
            print('patch_filename=%r' % patch_filename)
            print('basename=%r' % basename)
            raise
        #print('sline', sline)
        ipanel = int(sline)
        #print('ipanel = %s' % ipanel)

        bdf_model = BDF(debug=False)
        bdf_model.read_bdf(patch_filename, xref=False)

        eids = defaultdict(list)
        for eid, elem in iteritems(bdf_model.elements):
            pid = elem.pid
            key = (pid, ipanel)
            pid_panel[key].append(eid)
            eids[pid].append(eid)

        #if 0:
            #for pid, eidsi in sorted(iteritems(eids)):
                #pid_filename = os.path.join(workpath, 'pid_%i_ipanel_%i.txt' % (pid, ipanel))
                #out = str(eidsi)
                #with open(pid_filename, 'w') as pid_file:
                    #pid_file.write('# PSHELL pid\n')
                    #pid_file.write('# eids\n')
                    #pid_file.write('%s\n' % pid)
                    #pid_file.write('%s\n' % out[1:-1])

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
    return regions_filename


def main():  # pragma: no cover
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

    bdf_filename = 'model_144.bdf'
    op2_filenames = [os.path.join(workpath, op2_filename)
                     for op2_filename in ['patch_1.op2', 'patch_2.op2']]
    sym_regions_filename = 'sym_regions_map.csv'
    load_regions_and_create_eigenvalue_csv(bdf_filename, op2_filenames,
                                           'regions.txt', sym_regions_filename=sym_regions_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
