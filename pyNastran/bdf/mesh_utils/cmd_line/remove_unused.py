from __future__ import annotations
import os
import sys
from docopt import docopt
import numpy as np

import pyNastran
from .utils import filter_no_args, get_bdf_filename_punch_log


def cmd_line_remove_unused(argv=None, quiet: bool=False) -> None:
    """command line interface to remove_unused"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf remove_unused IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--lax]\n'
        '  bdf remove_unused -h | --help\n'
        '  bdf remove_unused -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF file\n'
        '  --punch                             flag to identify a *.pch/*.inc file\n'
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    out_bdf_filename = data['--output']
    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)

    if out_bdf_filename is None:
        abs_name = os.path.abspath(bdf_filename)
        dirname = os.path.dirname(abs_name)
        basename = os.path.basename(abs_name)
        out_bdf_filename = os.path.join(dirname, f'clean_{basename}')
        dict_filename = os.path.join(dirname, f'clean_summary_{basename}')
    else:
        dirname = os.path.dirname(out_bdf_filename)
        dict_filename = os.path.join(dirname, f'clean_summary.out')

    is_strict_card_parser = not data['--lax']

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
