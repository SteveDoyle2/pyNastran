import os
import glob
from six import iteritems

from pyNastran.applications.aero_panel_buckling.find_surface_panels import (
    find_surface_panels, get_bdf_object, get_op2_object)
from pyNastran.applications.aero_panel_buckling.split_model import (
    split_model_by_pid_panel, load_regions, load_regions_and_create_eigenvalue_csv)
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import run_bdfs#, load_sym_regions_map


def run_panel_buckling(bdf_filename='model_144.bdf', op2_filename='model_144.op2',
                       isubcase=1, workpath=None,
                       build_model=True, rebuild_patches=True, create_regions_filename=True,
                       run_nastran=True, nastran_keywords=None, overwrite_op2_if_exists=True,
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
    build_model : bool; default=True
        builds the patch model
    isubcase : int; default=1
        the case to analyze
    workpath : str; default=cwd/results
        the location the models will be created

    Step 2 - Run Jobs
    =================
    run_nastran : bool; default=True
        runs Nastran
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
    """
    if workpath is None:
        workpath = os.path.join(os.getcwd(), 'results')

    bdf_model = get_bdf_object(bdf_filename)
    if op2_filename is not None:
        op2_model = get_op2_object(op2_filename)
    else:
        op2_model = None

    bdf_model.log.info('build_model=%s rebuild_patches=%s' % (build_model, rebuild_patches))

    # isubcase=2 -> 2.5g pullup at Mach=0.85 at 11 km
    patch_dir = os.path.join(workpath, 'patches')
    if build_model:
        patch_filenames, edge_filenames = find_surface_panels(
            bdf_model, op2_model, mode=mode, isubcase=isubcase, consider_pids=True,
            rebuild_patches=rebuild_patches, workpath=workpath)
        create_regions_filename = True
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
        print('calling run_bdfs...')
        op2_filenames = run_bdfs(
            patch_filenames, nastran_keywords=nastran_keywords,
            workpath=patch_dir, overwrite_op2_if_exists=overwrite_op2_if_exists)
    else:
        if op2_filenames is None:
            op2_filenames = get_patch_filenames(patch_dir, extension='.op2')
        #op2_filenames = glob.glob(r'Y:\work\panel_buckling\patches\patch_*.op2')

    #sym_regions_filename = 'sym_regions_map.csv'
    sym_regions_filename = None
    split_model_by_pid_panel(patch_filenames, workpath=workpath)  # outputs regions.txt

    min_eigenvalue_by_patch_id, eigenvalues = load_regions_and_create_eigenvalue_csv(
        bdf_model, op2_filenames, 'regions.txt',
        eig_min=eig_min, eig_max=eig_max, eig_default=eig_default,
        sym_regions_filename=sym_regions_filename)
    return min_eigenvalue_by_patch_id, eigenvalues

def get_patch_filenames(patch_dir, extension='.bdf'):
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

if __name__ == '__main__':
    run_panel_buckling()

