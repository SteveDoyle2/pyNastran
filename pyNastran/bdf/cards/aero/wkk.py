import argparse
import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.utils import PathLike


def cmd_line_wkk(argv=None):
    """
    fname: series of:
     - GRID
     - CQUAD4
     - PLOAD2 (sid 1: force)
     - PLOAD2 (sid 2: moment)

    Parameters
    ----------
    dirname

    Returns
    -------

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("wkk", help='activates wkk')
    parser.add_argument("dirname", help='path to panel dirname')
    args = parser.parse_args()
    print(f'args = {args}')
    dirname = args.dirname
    solve_wkk_from_dirname(dirname=dirname)


def solve_wkk_from_dirname(dirname: PathLike):
    bdf_filenames = dirname.glob('*.bdf')
    solve_wkk_from_bdf_filenames(bdf_filenames)


def solve_wkk_from_bdf_filenames(bdf_filenames: list[PathLike]) -> np.ndarray:
    nfiles = len(bdf_filenames)
    forces, pressures = setup_pressures_from_bdf_filenames(bdf_filenames)
    if nfiles == 1:
        # solve the diagonal
        pass
    else:
        # A x = b
        pass
    Wkk = np.zeros((5, 5), dtype='float64')
    return Wkk

def setup_pressures_from_bdf_filenames(bdf_filenames: list[PathLike]) -> np.ndarray:
    if not isinstance(bdf_filenames, list):
        bdf_filenames = [bdf_filenames]

    force_list = []
    pressure_list = []
    for bdf_filename in bdf_filenames:
        if isinstance(bdf_filename, BDF):
            model = bdf_filename
        else:
            model = read_bdf(bdf_filename)

        # loadset_force = {}
        # loadset_pressure = {}
        nelement = len(model.elements)
        elements = np.array(list(model.elements), dtype='int32')
        elements.sort()
        force = np.zeros(nelement, dtype='float64')
        pressure = np.zeros(nelement, dtype='float64')
        # for eid, elem in sorted(model.elements.items()):
        #     area, centroid, normal = elem.AreaCentroidNormal()

        force_list.append(force)
        pressure_list.append(pressure)
        for sid, loads in sorted(model.loads.items()):
            if sid == 1:
                data = force
            else:
                assert sid == 2, sid
                data = pressure

            for load in loads:
                assert load.type == 'PLOAD2', load.get_stats()

                ieids = np.searchsorted(elements, load.eids)
                data[ieids] = load.pressure
                # for ieid, eid in zip(ieids, load.eids):
                #     print(load.get_stats())
                #     data[ieid] += load.pressure
    forces = np.column_stack(force_list, dtype='float64')
    pressures = np.column_stack(pressure_list, dtype='float64')
    return forces, pressures


def main():
    # cmd_line_wkk()
    model = BDF()
    eid = 1
    pid = 1
    mid1 = 1
    model.add_mat1(mid1, 3.0e7, None, 0.3)
    model.add_pshell(pid, mid1, t=0.1)
    model.add_cquad4(eid, pid, [1, 2, 3, 5])
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(5, [10., 1., 0.])
    model.add_pload2(1, 2.0, [1])
    model.add_pload2(2, 5.0, [1])
    model.cross_reference()
    solve_wkk_from_bdf_filenames([model])


if __name__ == '__main__':
    main()
