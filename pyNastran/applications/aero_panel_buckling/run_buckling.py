import os
import glob
from six import iteritems

from pyNastran.applications.aero_panel_buckling.find_surface_panels import find_surface_panels
from pyNastran.applications.aero_panel_buckling.split_model import split_model_by_pid_panel, load_regions, load_regions_and_create_eigenvalue_csv
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import run_bdfs, load_sym_regions_map


def run_panel_buckling(bdf_filename='model_144.bdf', op2_filename='model_144.op2',
                       isubcase=1, workpath=None,
                       build_model=False,
                       run_nastran=False, nastran_keywords=None, overwrite_op2_if_exists=True,
                       op2_filenames=None):
    """
    Step 1 - Setup
    ==============
    bdf_filename : str
        the path to the SOL 144 model
    op2_filename : str
        the path to the SOL 144 solution
    build_model : bool; default=False
        builds the patch model
    isubcase : int; default=1
        the case to analyze
    workpath : str; default=cwd/results
        the location the models will be created

    Step 2 - Run Jobs
    =================
    run_nastran : bool; default=False
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
    """
    if workpath is None:
        workpath = os.path.join(os.getcwd(), 'results')

    # isubcase=2 -> 2.5g pullup at Mach=0.85 at 11 km
    if build_model:
        find_surface_panels(bdf_filename, op2_filename, isubcase=isubcase, consider_pids=True,
                            rebuild_patches=False, workpath=workpath)
        split_model_by_pid_panel(workpath)

    regions_filename = os.path.join(workpath, 'regions.txt')
    #regions_filename = 'regions.txt'
    regions_to_pid_map, regions_to_eids_map = load_regions(regions_filename)

    # hardcoded...
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

    bdf_filenames2 = []
    bdf_filenames = glob.glob('%s/patches/patch_*.bdf' % workpath)
    assert len(bdf_filenames) > 0, bdf_filenames
    #print(bdf_filenames)
    for bdf_filename2 in bdf_filenames:
        basename = os.path.basename(bdf_filename2)
        patch_id_str = basename.split('_')[1].split('.')[0]
        patch_id = int(patch_id_str)
        #op2_filename = 'patch_%s.op2' % patch_id
        #if patch_id in sym_regions_map:
        bdf_filenames2.append(bdf_filename2)

    if run_nastran:
        run_bdfs(bdf_filenames2, nastran_keywords=nastran_keywords, overwrite_op2_if_exists=overwrite_op2_if_exists)

    #sym_regions_filename = 'sym_regions_map.csv'
    sym_regions_filename = None
    if op2_filenames is None:
        op2_filenames = glob.glob(r'Y:\work\panel_buckling\patches\patch_*.op2')

    split_model_by_pid_panel(workpath)
    load_regions_and_create_eigenvalue_csv(bdf_filename, op2_filenames, 'regions.txt', sym_regions_filename=sym_regions_filename)

if __name__ == '__main__':
    run_panel_buckling()

