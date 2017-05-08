"""
This is the main panel buckling file.  Defines:
 - run_panel_buckling
 - get_patch_filenames
 -
"""
from __future__ import print_function
import os
import glob
from six import iteritems

from pyNastran.applications.aero_panel_buckling.find_surface_panels import (
    create_plate_buckling_models, get_bdf_object, get_op2_object)
from pyNastran.applications.aero_panel_buckling.split_model import (
    split_model_by_pid_panel, load_regions, load_regions_and_create_eigenvalue_csv)
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import run_bdfs
#, load_sym_regions_map


def run_panel_buckling(bdf_filename='model_144.bdf', op2_filename='model_144.op2',
                       isubcase=1, workpath=None,
                       build_model=True, rebuild_patches=True, rebulid_buckling_bdfs=True,
                       create_regions_filename=True,
                       run_nastran=True, nastran_keywords=None, overwrite_op2_if_exists=True,
                       parse_eigenvalues=True,
                       op2_filenames=None,
                       mode='displacement',
                       eig_min=-1.0, eig_max=1.0, eig_default=3.0):
    """
    Step 1 - Setup
    ==============
    bdf_filename/bdf_model : str/BDF(); default=model_144.bdf
        str : the path to the SOL 144 model
        BDF : a BDF model
    op2_filename/op2_model : str/OP2(); default=model_144.op2
        the path to the SOL 144 solution
        OP2 : an OP2 model
    isubcase : int; default=1
        the case to analyze
    workpath : str; default=cwd/results
        the location the models will be created
    build_model : bool; default=True
        builds the edge_*.csv, patch_*.csv, patch_*.bdf files
        True : sets rebuild_patches=True and rebulid_buckling_bdfs=True
        False : set rebuild_patches and rebulid_buckling_bdfs by hand
    rebuild_patches : bool; default=True
        True : rebuilds the following files:
                - element_edges.csv
                - element_patches.csv
                - free_edge_nodes.csv
                - free_edge_nodes_xyz.csv
                - nodal_edges.csv
                - patch_edges_array.csv
                - patch_edges_xyz.csv
                - patches.csv
        False : does not create these files
    rebulid_buckling_bdfs : bool; default=True
        True : creates the edge_*.csv, patch_*.csv, patch_*.bdf files
        False : does not create these files
    create_regions_filename : bool; default=True
        True : creates regions.txt from patch_filenames
        False : assumes the file already exists
        # pid, ipanel, eids
        2, 18, 22016, 22017, 22018, 22019...
        20301, 15, 3584, 3585, 3586, ...
    mode : str; default='displacement'
        'displacement' : extract displacements from the op2 and apply
            interface displacements for each patch
        'load' : extract loads from the op2 and apply interface loads to
            each patch and constrain a boundary node
    parse_eigenvalues : bool; default=True
        calculate the outputs
        requires patch_*.op2 files exist

    Step 2 - Run Jobs
    =================
    run_nastran : bool; default=True
        runs Nastran on the patch_*.bdf files
    nastran_keywords : dict; default=None
        key : nastran key
        value : value for key

        default:
          nastran_keywords = {
              'mem' : mem='100mb',
              'scr' : 'yes',
              'bat' : 'no',
              'old' : 'no',
              'news' : 'no',
          }
        the Nastran arguments
    overwrite_op2_if_exists : default=True
        regenerates the OP2 if it already exists

    Step 3 - Post Process
    =====================
    op2_filenames : List[str]; default=None
        the list of op2 filenames generated from the previous step
    eig_min : float
        the required minimum eigenvalue
    eig_max : float
        the required maximum eigenvalue
    eig_default : float
        the default eigenvalue for cases that do not calculate eigenvectors
        because there were no eigenvalues in the range

    Returns
    -------
    min_eigenvalue_by_patch_id : dict
        key : patch_id : int
            the integer patch id
        value : eigenvalue or reserve factor
           the reserve factor eigenvalue for buckling
        None if parse_eigenvalues=False
    eigenvalues : (n, ) float ndarray
        the minimum eigenvalues
        None if parse_eigenvalues=False
    """
    if workpath is None:
        workpath = os.path.join(os.getcwd(), 'results')

    bdf_model = get_bdf_object(bdf_filename)
    if op2_filename is not None:
        op2_model = get_op2_object(op2_filename)
    else:
        op2_model = None

    bdf_model.log.info('build_model=%s rebuild_patches=%s rebulid_buckling_bdfs=%s' % (
        build_model, rebuild_patches, rebulid_buckling_bdfs))

    # isubcase=2 -> 2.5g pullup at Mach=0.85 at 11 km
    patch_dir = os.path.join(workpath, 'patches')
    if any([build_model, rebuild_patches, rebulid_buckling_bdfs]):
        if build_model:
            rebuild_patches = True
            rebulid_buckling_bdfs = True
        patch_filenames, edge_filenames = create_plate_buckling_models(
            bdf_model, op2_model, mode=mode, isubcase=isubcase, consider_pids=True,
            rebuild_patches=rebuild_patches, rebulid_buckling_bdfs=rebulid_buckling_bdfs,
            workpath=workpath)
        create_regions_filename = True

        if not rebulid_buckling_bdfs:
            return patch_filenames, edge_filenames
    else:
        patch_filenames = get_patch_filenames(patch_dir, extension='.bdf')

    if create_regions_filename:
        regions_filename = split_model_by_pid_panel(patch_filenames, workpath)

    regions_filename = os.path.join(workpath, 'regions.txt')
    #regions_filename = 'regions.txt'
    regions_to_pid_map, regions_to_eids_map = load_regions(regions_filename)

    # hardcoded...
    if 0:
        offset = 74305
        min_region_eid = {}
        sym_regions_map = {}
        min_region_eid_sym = {}
        for region_id, eids in iteritems(regions_to_eids_map):
            eid_min = min(eids)
            if eid_min < offset:
                min_region_eid[region_id] = eid_min
            else:
                min_region_eid_sym[eid_min] = region_id

        for region_id, eid_min in iteritems(min_region_eid):
            sym_region_id = min_region_eid_sym[eid_min + offset]
            sym_regions_map[region_id] = sym_region_id

    write_maps = False
    if write_maps:
        sym_regions_filename = 'sym_regions_map.csv'
        with open(sym_regions_filename, 'w') as sym_regions_file:
            sym_regions_file.write('RegionID, SymmetricalRegionID\n')
            for region_id, sym_region_id in sorted(iteritems(sym_regions_map)):
                sym_regions_file.write('%s, %s\n' % (region_id, sym_region_id))

        condensed_regions_filename = 'condensed_regions.txt'
        with open(condensed_regions_filename, 'w') as condensed_regions_file:
            condensed_regions_file.write('# pid, ipanel, eids\n')
            ipanel = 0
            for region_id, sym_region_id in sorted(iteritems(sym_regions_map)):
                pid = regions_to_pid_map[region_id]
                pid_sym = regions_to_pid_map[sym_region_id]
                assert pid == pid_sym, 'pid=%s pid_sym=%s' % (pid, pid_sym)
                eids = regions_to_eids_map[pid]
                eids_sym = regions_to_eids_map[pid_sym]
                out = str(eids + eids_sym)

                condensed_regions_file.write('%s, %s, ' % (pid, ipanel))
                condensed_regions_file.write('%s\n' % out[1:-1])
                ipanel += 1

    if run_nastran:
        bdf_model.log.info('run_nastran=True')
        op2_filenames = run_bdfs(
            patch_filenames, nastran_keywords=nastran_keywords,
            workpath=patch_dir, overwrite_op2_if_exists=overwrite_op2_if_exists)
    elif op2_filenames is None:
        bdf_model.log.info('run_nastran=False; automatically determining op2s...')
        op2_filenames = get_patch_filenames(patch_dir, extension='.op2')
        #print(op2_filenames)
        #op2_filenames = glob.glob(r'Y:\work\panel_buckling\patches\patch_*.op2')

    #sym_regions_filename = 'sym_regions_map.csv'
    sym_regions_filename = None
    split_model_by_pid_panel(patch_filenames, workpath=workpath)  # outputs regions.txt

    if parse_eigenvalues:
        bdf_model.log.info('parse_eigenvalues=True...')
        min_eigenvalue_by_patch_id, eigenvalues = load_regions_and_create_eigenvalue_csv(
            bdf_model, op2_filenames, 'regions.txt',
            eig_min=eig_min, eig_max=eig_max, eig_default=eig_default,
            sym_regions_filename=sym_regions_filename)
        return min_eigenvalue_by_patch_id, eigenvalues
    else:
        bdf_model.log.info('skipping parse_eigenvalues...')

    return None, None

def get_patch_filenames(patch_dir, extension='.bdf'):
    """
    Gets all the patch_*.bdf files

    Parameters
    ----------
    patch_dir : str
        the path to the patch_*.bdf files
    extension : str; default='.bdf'
        the file extension

    Returns
    -------
    patch_filenames : List[str]
        the BDFs to run
        [patch_0.csv, patch_1.csv, ...]
    """
    patch_filenames = []
    bdf_filenames = [fname for fname in os.listdir(patch_dir)]
    for bdf_filename2 in bdf_filenames:
        basename = os.path.basename(bdf_filename2)
        #print('basename = %r' % basename)
        #print('basename[:6] = %r' % basename[:6])
        #print('basename[-5:] = %r' % basename[-4:])
        if basename[:6] == 'patch_' and basename[-4:] == extension:
            patch_filenames.append(os.path.join(patch_dir, basename))
        #else:
            #print('***skipping %r' % basename)

    assert len(bdf_filenames) > 0, bdf_filenames

    #for bdf_filename2 in patch_filenames:
        #print(bdf_filename2)
        #basename = os.path.basename(bdf_filename2)
        #patch_id_str = basename.split('_')[1].split('.')[0]
        #patch_id = int(patch_id_str)
        #op2_filename = 'patch_%s.op2' % patch_id
        #if patch_id in sym_regions_map:
        #bdf_filenames2.append(bdf_filename2)
    return patch_filenames

if __name__ == '__main__':  # pragma: no cover
    run_panel_buckling()

