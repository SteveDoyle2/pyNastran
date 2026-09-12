from __future__ import annotations
import os
import sys
import argparse
import numpy as np
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args


def cmd_line_remove_unused(argv=None, quiet: bool=False) -> None:
    """command line interface to remove_unused"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf remove_unused: use 'bdf remove_unused -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf remove_unused',
        description='Remove unused cards from a BDF model',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF file')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--lax', action='store_true',
                        help='lax card parser')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    out_bdf_filename = args.output
    if out_bdf_filename is None:
        abs_name = os.path.abspath(bdf_filename)
        dirname = os.path.dirname(abs_name)
        basename = os.path.basename(abs_name)
        basename_noext = os.path.splitext(basename)[0]
        out_bdf_filename = os.path.join(dirname, f'clean_{basename}')
        dict_filename = os.path.join(dirname, f'clean_summary_{basename_noext}.out')
    else:
        dirname = os.path.dirname(out_bdf_filename)
        dict_filename = os.path.join(dirname, f'clean_summary.out')

    is_strict_card_parser = not args.lax

    from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
    from .utils_bdf import read_lax_bdf
    model = read_lax_bdf(
        bdf_filename, punch=punch, xref=False,
        is_strict_card_parser=is_strict_card_parser,
        log=log)
    #model.cross_reference()
    model, out_dict = remove_unused(
        model,
        remove_nids=True, remove_cids=True,
        remove_pids=True, remove_mids=True,
        remove_spcs=True, remove_mpcs=True,
        remove_optimization=True,
        reset_type_to_id_map=False)

    if os.path.exists(dict_filename):
        os.remove(dict_filename)

    if out_dict:
        with open(dict_filename, 'w') as dict_file:
            dict_file.write('removed:\n')
            for key, myarray in out_dict.items():
                assert isinstance(key, str), key
                assert isinstance(myarray, np.ndarray), (key, myarray)
                ids = myarray.tolist()
                dict_file.write(f'  {key} = {ids}\n')

    model.write_bdf(out_bdf_filename,
                    nodes_size=None,
                    is_double=False, interspersed=False)
    #for iply in iplies:
        #csv_filename = csv_filename_base + '_ply=%i.csv' % iply
        # export_mcids(
        #     model, csv_filename,
        #     export_xaxis=export_xaxis, export_yaxis=export_yaxis, iply=iply)
        #model.log.info('wrote %s' % csv_filename)
