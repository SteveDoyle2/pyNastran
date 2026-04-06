import numpy as np
from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2


def downselect(
        bdf_filename: BDF,
        op2_filename: OP2,
        eids=None, nids=None,
        subcases=None,
        percent_eids_target: float=0.90):
    """
    Analysis is done in the material coordinate frame for composites
    and in the element frame for shells/solids/bars

    def theta_func_0_90(thetas: list[float]):
        thetas_out = []
        for theta in thetas:
            if np.isclose(abs(theta), 0.) or np.isclose(abs(theta), 90.):
                thetas_out.append(theta)
        return thetas_out

    # all criterion must be satisfied
    criteria = {
        'PSHELL': {
            1: {'eid': eids, 'pid': pids,   'result': {'vm': 20000,}},
            2: {             'mid1': mid1s, 'result': {'max': 35000., 'vm', 21000}},
        },
        'PCOMP': {
            1: {'eid': eids, 'ply': [1, 2, 3], 'result': {'e1': 2500}},
            2: {'eid': eids, 'theta': [0., 90.], 'result': {'e1': 2500}},
            3: {'eid': eids, 'theta': theta_func_0_90, 'result': {'e1': 2500}},
        },
        'PCOMPG': {
            2: {'glply: [4, 5, 6], 'result': {'e1': 2500}},
        },
    }
    PSHELL:
       - stress: quantity/allowable (vm, max, min)
      #- strain: quantity/allowable
    PCOMP:
      - strain quantity/allowable (ply, e1, e2, e12)
      - strain quantity/allowable (core)
    SOLID:
      - stress quantity/allowable (max/vm)
    """
    log = SimpleLogger(level='debug')
    model = read_bdf(bdf_filename, cards_to_skip=cards_to_skip, log=log)

    element_id = np.array(list(model.elements), dtype='int32')

    model_results = read_op2(op2_filename, subcases=subcases, debug=None)
    subcase_id = 1

    eid_to_nid_map = {}
    plate_shell_eids = []
    plate_comp_eids = []
    solid_eids = []
    for eid, elem in model.elements.items():
        eid_to_nid_map[eid] = elem.nodes
        prop = elem.pid_ref
        if prop.type in {'PSHELL'}:
            plate_shell_eids.append(eid)
        elif prop.type in {'PCOMP', 'PCOMPG'}:
            plate_comp_eids.append(eid)
        elif prop.type in {'PSOLID'}:
            solid_eids.append(eid)
        elif prop.type in {'PBUSH', 'PBUSH1D', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBEND'}:
            log.warning(f'skipping {prop.type}')
        else:
            log.error(f'skipping {prop.type}')


    out = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int64')
    nid_cp_cd, unused_xyz_cid, unused_xyz_cp, unused_icd_transform, unused_icp_transform = out
    node_id = nid_cp_cd[:, 0]

    is_stress = False
    obj2 = PlateStrainStressResults2.load_from_code(
        subcase_id, model, model_results, element_id,
        is_stress, eid_to_nid_map,
        node_id=node_id,
        # is_variable_data_format=False,
        require_results=False,
    )
    # plates = PlateStrainStressResults2(
    #     subcase_id, model,
    #     node_idy,
    #     element_id,
    #     cases, #: list[RealPlateArray],
    #     result, # str
    #     title='plate',
    #     is_fiber_distance,
    #     eid_to_nid_map,
    #     # dim_max: float=1.0,
    #     data_format='%g',
    #     is_variable_data_format=False,
    #     nlabels=None, labelsize=None, ncolors=None,
    #     colormap='', set_max_min=False,
    #     uname='PlateStressStrainResults2')

    pass

def main():
    import pyNastran
    from pathlib import Path
    PKG_PATH = Path(pyNastran.__file__[0]).parent
    model_path = PKG_PATH / 'models'
    bdf_filename = model_path / 'model.bdf'
    op2_filename = model_path / 'model.op2'
    downselect(bdf_filename, op2_filename)

if __name__ == '__main__':
    main()
