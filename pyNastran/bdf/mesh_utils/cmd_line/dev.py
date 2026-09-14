from __future__ import annotations
import sys
import argparse
from typing import Optional, TYPE_CHECKING
from cpylog import SimpleLogger
import numpy as np

import pyNastran
from pyNastran.bdf.mesh_utils.shift import update_nodes
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

from .utils import filter_no_args
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_bin(argv=None, quiet: bool=False) -> None:  # pragma: no cover
    """bins the model into nbins"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    if len(argv) == 1:
        sys.exit("bdf bin: use 'bdf bin -h' for help")

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf bin',
        description='Bin a BDF model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Plot z (2) as a function of y (1) in y-stepsizes of 0.1:\n'
            '  bdf bin fem.bdf 1 2 --cid 0 --step 0.1\n\n'
            'Plot z (2) as a function of y (1) with 50 bins:\n'
            '  bdf bin fem.bdf 1 2 --cid 0 --nbins 50\n'
        ),
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('AXIS1', type=int,
                        help='axis to loop over')
    parser.add_argument('AXIS2', type=int,
                        help='axis to bin')
    parser.add_argument('--cid', type=int, default=0,
                        help='the coordinate system to bin (default=0)')
    parser.add_argument('--nbins', type=int, default=10,
                        help='the number of bins (default=10)')

    args = parser.parse_args(argv[2:])
    bdf_filename = args.IN_BDF_FILENAME
    axis1 = args.AXIS1
    axis2 = args.AXIS2
    cid = args.cid
    nbins = args.nbins
    assert nbins >= 2, nbins
    if not quiet:  # pragma: no cover
        print(vars(args))

    from pyNastran.bdf.bdf import read_bdf
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    model = read_bdf(bdf_filename, log=log)
    bin_model(model, axis1, axis2, cid=cid, nbins=nbins, debug=quiet)


def bin_model(model: BDF, axis1: int, axis2: int,
              cid: int=0, nbins: int=10, debug: bool=False) -> None:
    import matplotlib.pyplot as plt
    xyz_cid = model.get_xyz_in_coord(cid=cid, fdtype='float64')
    y = xyz_cid[:, axis1]
    z = xyz_cid[:, axis2]

    plt.figure(1)
    #n, bins, patches = plt.hist( [x0,x1,x2], 10, weights=[w0, w1, w2], histtype='bar')
    ys = []
    #zs = []
    zs_min = []
    zs_max = []
    y0 = y.min()
    y1 = y.max()
    dy = (y1 - y0) / nbins
    y0i = y0
    y1i = y0 + dy
    for unused_i in range(nbins):
        j = np.where((y0i <= y) & (y <= y1i))[0]
        if not len(j):
            continue
        ys.append(y[j].mean())
        zs_min.append(z[j].min())
        zs_max.append(z[j].max())
        y0i += dy
        y1i += dy
    zs_max = np.array(zs_max)
    zs_min = np.array(zs_min)
    if not debug:  # pragma: no cover
        print('ys = %s' % ys)
        print('zs_max = %s' % zs_max)
        print('zs_min = %s' % zs_min)
    plt.plot(ys, zs_max, 'r-o', label='max')
    plt.plot(ys, zs_min, 'b-o', label='min')
    plt.plot(ys, zs_max - zs_min, 'g-o', label='delta')
    #plt.xlim([y0, y1])
    plt.xlabel('Axis %s' % axis1)
    plt.ylabel('Axis %s' % axis2)
    plt.grid(True)
    plt.legend()
    plt.show()


def cmd_line_transform(argv=None, quiet: bool=False) -> None:
    """command line interface to transform"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf transform: use 'bdf transform -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf transform',
        description='Transform a BDF model',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='transform.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF file (default=transform.bdf)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--shift', default=None, metavar='XYZ',
                        help='shift the model by XYZ (comma separated, e.g. 1.0,2.0,3.0)')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    bdf_filename_out = args.output

    dxyz = None
    if args.shift is not None:
        xyz = args.shift.split(',')
        dxyz = np.array(xyz, dtype='float64')
        assert len(dxyz) == 3, dxyz

    from pyNastran.bdf.bdf import read_bdf
    model = read_bdf(bdf_filename, punch=punch, log=log)

    nid_cp_cd, xyz_cid0, unused_xyz_cp, unused_icd_transform, unused_icp_transform = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int32')

    update_nodes_flag = False
    # we pretend to change the SPOINT location
    if dxyz is not None:
        xyz_cid0 += dxyz
        update_nodes_flag = True

    if update_nodes_flag:
        update_nodes(model, nid_cp_cd, xyz_cid0)
        model.write_bdf(bdf_filename_out)


def cmd_line_filter(argv=None, quiet: bool=False) -> None:  # pragma: no cover
    """command line interface to bdf filter"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf filter: use 'bdf filter -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf filter',
        description='Filter a BDF model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  1. remove unused cards:\n'
            '     bdf filter fem.bdf\n'
            '  2. remove GRID points and associated cards with y value < 0:\n'
            "     bdf filter fem.bdf --y '< 0.'\n"
        ),
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='filter.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF file (default=filter.bdf)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--x', default=None, metavar='YSIGN_X',
                        help="a string (e.g., '< 0.')")
    parser.add_argument('--y', default=None, metavar='YSIGN_Y',
                        help="a string (e.g., '< 0.')")
    parser.add_argument('--z', default=None, metavar='YSIGN_Z',
                        help="a string (e.g., '< 0.')")

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    bdf_filename_out = args.output

    func_map = {
        '<': np.less,
        '>': np.greater,
        '<=': np.less_equal,
        '>=': np.greater_equal,
    }
    xsign = None
    ysign = None
    zsign = None
    if args.x:
        xsign, xval_str = args.x.split(' ')
        xval = float(xval_str)
        assert xsign in ['<', '>', '<=', '>='], xsign
    if args.y:  # --y < 0
        ysign, yval_str = args.y.split(' ')
        yval = float(yval_str)
        assert ysign in ['<', '>', '<=', '>='], ysign
    if args.z:
        zsign, zval_str = args.z.split(' ')
        zval = float(zval_str)
        assert zsign in ['<', '>', '<=', '>='], zsign

    from pyNastran.bdf.bdf import read_bdf
    model = read_bdf(bdf_filename, log=log, punch=punch)

    #nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array(
        #cid=0, fdtype='float64', idtype='int32')

    eids = []
    xyz_cid0 = []
    for eid, elem in sorted(model.elements.items()):
        xyz = elem.Centroid()
        xyz_cid0.append(xyz)
        eids.append(eid)
    xyz_cid0 = np.array(xyz_cid0)
    eids = np.array(eids)

    # we pretend to change the SPOINT location
    update_nodesi = False
    # we pretend to change the SPOINT location
    iunion = None
    if xsign:
        xvals = xyz_cid0[:, 0]
        xfunc = func_map[xsign]
        ix = xfunc(xvals, xval)
        iunion = _union(xval, ix, iunion)
        update_nodesi = True
    if ysign:
        yvals = xyz_cid0[:, 1]
        yfunc = func_map[ysign]
        iy = yfunc(yvals, yval)
        iunion = _union(yval, iy, iunion)
        update_nodesi = True
    if zsign:
        zvals = xyz_cid0[:, 2]
        zfunc = func_map[zsign]
        iz = zfunc(zvals, zval)
        iunion = _union(zval, iz, iunion)
        update_nodesi = True

    if update_nodesi:
        eids_to_remove = eids[iunion]
        for eid in eids_to_remove:
            etype = model.elements[eid].type
            model._type_to_id_map[etype].remove(eid)
            del model.elements[eid]

    #update_nodes(model, nid_cp_cd, xyz_cid0)
    # unxref'd model
    remove_unused(model, remove_nids=True, remove_cids=True,
                  remove_pids=True, remove_mids=True)
    model.write_bdf(bdf_filename_out)


def _union(xval: float,
           iunion: np.ndarray,
           ix: Optional[np.ndarray]) -> np.ndarray:
    """helper method for ``filter``"""
    if xval:
        if iunion:
            ## TODO: ix can be None?
            iunion = np.union1d(iunion, ix)
        else:
            pass
    return iunion
