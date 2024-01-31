import os
from copy import deepcopy
import numpy as np
from pyNastran.bdf.bdf import BDF, CaseControlDeck

def setup_statics(lines=None) -> tuple[BDF, CaseControlDeck]:
    if lines is None:
        lines = []
    model = BDF()
    model.sol = 101
    #model.set_as_mystran()
    lines += [
        'STRESS(PLOT,PRINT) = ALL',
        'STRAIN(PLOT,PRINT) = ALL',
        'FORCE(PLOT,PRINT) = ALL',
        'SPCFORCE(PLOT,PRINT) = ALL',
        'DISP(PLOT,PRINT) = ALL',
        'MPCFORCE(PLOT,PRINT) = ALL',
        'GPFORCE(PLOT) = ALL',
    ]
    cc = CaseControlDeck(lines)
    model.case_control_deck = cc
    return model, cc

def add_forces_by_6_cases(i, model, sid, nid, F: float):
    vector = np.zeros(3)
    if i in {0, 1, 2}:
        vector[i] = F
        model.add_force(sid, nid, 1.0, vector)
    else:
        vector[i-3] = F
        model.add_moment(sid, nid, 1.0, vector)

def build_celas1(dirname: str=''):
    model, cc = setup_statics()

    eid = 1
    pid = 2
    k = 1000.
    F = 20.
    dx = F / k
    model.add_pelas(pid, k)

    all_nids2 = []
    for i in range(6):
        i0 = (i + 1) * 10
        nid2 = i0 + 2
        all_nids2.append(nid2)

    for i in range(6):
        component = i + 1
        sid = component
        subcase = cc.create_new_subcase(sid)
        subcase.add_integer_type('LOAD', sid)
        subcase.add_integer_type('SPC', sid)

        i0 = (i + 1) * 10
        eid = component
        nid1 = i0 + 1
        nid2 = i0 + 2
        model.add_grid(nid1, [0., 0., i], comment=f'eid={eid} load_id/spc_id={sid}')
        model.add_grid(nid2, [1., 0., i])
        nids = [nid1, nid2]
        model.add_celas1(eid, pid, nids, c1=component, c2=component, comment='')
        model.add_spc1(sid, str(component), nid1)

        add_forces_by_6_cases(i, model, sid, nid2, F)
        local_nids = deepcopy(all_nids2)
        local_nids.remove(nid2)
        model.add_spc1(sid, '123456', nid2)
    bdf_filename = os.path.join(dirname, f'celas1_pelas_dir.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)

def build_celas1_ground1(dirname: str=''):
    model, cc = setup_statics()
    eid = 1
    pid = 2
    k = 1000.
    F = 20.
    dx = F / k
    model.add_pelas(pid, k)

    all_nids2 = []
    for i in range(6):
        i0 = (i + 1) * 10
        nid2 = i0 + 2
        all_nids2.append(nid2)

    for i in range(6):
        component = i + 1
        sid = component
        subcase = cc.create_new_subcase(sid)
        subcase.add_integer_type('LOAD', sid)
        #subcase.add_integer_type('SPC', sid)


        i0 = (i + 1) * 10
        eid = component
        nid1 = i0 + 1
        nid2 = i0 + 2
        #nid3 = i0 + 3
        #nid4 = i0 + 4
        model.add_grid(nid1, [0., 0., i], comment=f'eid={eid} load_id={sid}')
        model.add_grid(nid2, [1., 0., i])
        nids = [nid1, nid2]
        model.add_celas1(eid, pid, nids, c1=None, c2=component, comment='')

        add_forces_by_6_cases(i, model, sid, nid2, F)
        local_nids = deepcopy(all_nids2)
        local_nids.remove(nid2)
        model.add_spc1(sid, '123456', nid2)

    bdf_filename = os.path.join(dirname, f'celas1_ground1.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)


def build_celas1_ground2(dirname: str=''):
    model, cc = setup_statics()
    eid = 1
    pid = 2
    k = 1000.
    F = 20.
    dx = F / k
    mag = 1.0
    model.add_pelas(pid, k)

    all_nids2 = []
    for i in range(6):
        i0 = (i + 1) * 10
        nid2 = i0 + 2
        all_nids2.append(nid2)

    for i in range(6):
        component = i + 1
        sid = component
        subcase = cc.create_new_subcase(sid)
        subcase.add_integer_type('LOAD', sid)
        #subcase.add_integer_type('SPC', sid)


        i0 = (i + 1) * 10
        eid = component
        nid1 = i0 + 1
        nid2 = i0 + 2
        #nid3 = i0 + 3
        #nid4 = i0 + 4
        model.add_grid(nid1, [0., 0., i], comment=f'eid={eid} load_id={sid}')
        model.add_grid(nid2, [1., 0., i])
        nids = [nid1, nid2]
        model.add_celas1(eid, pid, nids, c1=component, c2=None, comment='')

        add_forces_by_6_cases(i, model, sid, nid1, F)
        #local_nids = deepcopy(all_nids2)
        #local_nids.remove(nid2)
        #model.add_spc1(sid, '123456', nid2)

    bdf_filename = os.path.join(dirname, f'celas1_ground2.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)

def build_cbar_pbar_ydir(dirname: str):
    model = _build_cbar_pbar_ydir()
    bdf_filename = os.path.join(dirname, f'cbar_pbar_ydir.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)

def build_cbar_pbar_zdir(dirname: str):
    model = _build_cbar_pbar_ydir()

    for elem in model.elements.values():
        elem.x = [0., 0., 1.]
    bdf_filename = os.path.join(dirname, f'cbar_pbar_zdir.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)

def build_cbar_pbar_pin_ydir(dirname: str):
    model = _build_cbar_pbar_pin_ydir()
    bdf_filename = os.path.join(dirname, f'cbar_pbar_pin_ydir.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)


def build_cbar_pbar_pin_zdir(dirname: str):
    model = _build_cbar_pbar_pin_ydir()

    for elem in model.elements.values():
        elem.x = [0., 0., 1.]
    bdf_filename = os.path.join(dirname, f'cbar_pbar_pin_zdir.bdf')
    model.write_bdf(bdf_filename, enddata=True, write_header=True)

def add_steel(model, mid):
    E = 3.0e7
    G = None
    nu = 0.3
    rho = 0.1
    model.add_mat1(mid, E, G, nu, rho=rho)

def _build_cbar_pbar_ydir():
    spc_id = 1
    lines = [f'SPC = {spc_id}']
    model, cc = setup_statics(lines)
    eid = 1
    pid = 2

    k = 1000.
    F = 20.
    dx = F / k
    mid = 10
    add_steel(model, mid)
    model.add_pbar(pid, mid,
                   A=0.3, i1=0.5, i2=0.7, i12=0.11, j=0.233, nsm=0.,
                   c1=1., c2=1.,
                   d1=1., d2=-1.,
                   e1=-1., e2=1.,
                   f1=-1., f2=-1.,
                   k1=1.e8, k2=1.e8)
    for i in range(6):
        component = i + 1
        sid = component
        subcase = cc.create_new_subcase(sid)
        subcase.add_integer_type('LOAD', sid)

        i0 = (i + 1) * 10
        eid = component
        nid1 = i0 + 1
        nid2 = i0 + 2
        model.add_grid(nid1, [0., 0., i], comment=f'eid={eid} load_id={sid}')
        model.add_grid(nid2, [1., 0., i])
        nids = [nid1, nid2]
        x = [0., 1., 0.]
        model.add_cbar(eid, pid, nids,
                 x, g0=None,
                 offt='GGG',
                 pa=0, pb=0,
                 wa=None, wb=None,
                 validate=True)
        model.add_spc1(spc_id, '123456', nid1)
        add_forces_by_6_cases(i, model, sid, nid2, F)
    return model

def _build_cbar_pbar_pin_ydir():
    spc_id = 1
    lines = [f'SPC = {spc_id}']
    model, cc = setup_statics(lines)
    eid = 1
    pid = 2

    k = 1000.
    F = 20.
    dx = F / k

    mid = 10
    add_steel(model, mid)
    model.add_pbar(pid, mid,
                   A=0.3, i1=0.5, i2=0.7, i12=0.11, j=0.233, nsm=0.,
                   c1=1., c2=1.,
                   d1=1., d2=-1.,
                   e1=-1., e2=1.,
                   f1=-1., f2=-1.,
                   k1=1.e8, k2=1.e8)

    for i in range(6):
        component = i + 1
        sid = component
        subcase = cc.create_new_subcase(sid)
        subcase.add_integer_type('LOAD', sid)
        #subcase.add_integer_type('SPC', sid)


        i0 = (i + 1) * 10
        eid1 = i0 + 1
        eid2 = i0 + 2
        eid3 = i0 + 3
        eid4 = i0 + 4
        nid1 = i0 + 1
        nid2 = i0 + 2
        nid3 = i0 + 3
        nid4 = i0 + 4
        model.add_grid(nid1, [0., 0., i], comment=f'eid={eid} load_id={sid}')
        model.add_grid(nid2, [1., 0., i])
        model.add_grid(nid3, [2., 0., i])
        model.add_grid(nid4, [3., 0., i])
        x = [0., 1., 0.]

        #  1           2          3            4
        # SPC----------o----------o----------SPC
        model.add_cbar(eid1, pid, [nid1, nid2],
                       x, g0=None,
                       offt='GGG',
                       pa=0, pb=0,
                       wa=None, wb=None,
                       validate=True)
        pa_list = ['1', '2', '3', '4', '5', '6']
        pa_list.remove(str(component))
        pa = ''.join(pa_list)

        model.add_cbar(eid2, pid, [nid2, nid3],
                       x, g0=None,
                       offt='GGG',
                       pa=pa, pb=0,
                       wa=None, wb=None,
                       validate=True)
        model.add_cbar(eid3, pid, [nid3, nid4],
                       x, g0=None,
                       offt='GGG',
                       pa=0, pb=0,
                       wa=None, wb=None,
                       validate=True)
        model.add_spc1(spc_id, '123456', [nid1, nid4])
        add_forces_by_6_cases(i, model, sid, nid3, F)
    return model

def main():
    dirname = os.path.dirname(__file__)
    build_cbar_pbar_pin_ydir(dirname)
    build_cbar_pbar_pin_zdir(dirname)

    build_celas1(dirname)
    build_celas1_ground1(dirname)
    build_celas1_ground2(dirname)
    build_cbar_pbar_ydir(dirname)
    build_cbar_pbar_zdir(dirname)
    x = 1


if __name__ == '__main__':
    main()
