from __future__ import annotations
import os
import sys

from pygments.lexer import default

from cpylog import SimpleLogger
from .utils import add_argparse_arguments
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties_nsm
# import pyNastran


def cmd_line_mass(argv=None, quiet: bool=False) -> None:
    """bdf mass model.bdf --nsm 7100000 --obj --element_ids "7100001,7101743,1;7100001"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # positional arguments
    parent_parser.add_argument('mass', type=str)
    parent_parser.add_argument('INPUT', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('--nsm', help='NSM id', type=int, default=0)
    parent_parser.add_argument('--element_ids', help='element_ids', type=str, default='')
    #parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)
    parent_parser.add_argument('--no_prop_mass', action='store_true', help='remove property/material mass')
    parent_parser.add_argument('--obj', action='store_true', help='save the obj for reload')
    add_argparse_arguments(parent_parser, ['--punch', '--lax', '--allow_dup'])
    args = parent_parser.parse_args(args=argv[1:])
    # if not quiet:  # pragma: no cover
    #     print(args)

    bdf_filename = args.INPUT
    # bdf_filename_out = args.OUTPUT
    # is_strict_card_parser = not args.lax

    # bdf_filename_out = ''
    # if bdf_filename_out is None:
    #     bdf_filename_base, ext = os.path.splitext(bdf_filename)
    #     bdf_filename_out = f'{bdf_filename_base}.zip{ext}'
    duplicate_cards = args.allow_dup.split(',') if args.allow_dup else []
    nsm_id = args.nsm if args.nsm else None

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    from .utils_bdf import read_lax_obj, get_ids
    ('7100001,7101743,1;'
     '')
    obj_filename = os.path.splitext(bdf_filename)[0] + '.obj'
    element_ids = get_ids(args.element_ids)
    model = read_lax_obj(
        bdf_filename, obj_filename, args.obj,
        xref=False,
        is_strict_card_parser=not args.lax,
        punch=args.punch,
        duplicate_cards=duplicate_cards)

    model.rigid_elements = {}
    if len(element_ids):
        print(element_ids)
        set_eids = set(element_ids)
        eids_to_delete = set(list(model.elements)) - set_eids
        for eid in eids_to_delete:
            del model.elements[eid]
        eids_to_delete = set(list(model.masses)) - set_eids
        for eid in eids_to_delete:
            del model.masses[eid]

        # TODO: need to consider coords
        # nids_used = set()
        # for eid, elem in model.masses.items():
        #     nids_used.update(elem.node_ids)
        # for eid, elem in model.elements.items():
        #     nids_used.update(elem.nodes)
        # all_nids = set(list(model.nodes))
        # nids_to_remove = all_nids - nids_used
        # for nid in nids_to_remove:
        #     del model.nodes[nid]
        element_ids = list(model.elements)
        mass_ids = list(model.masses)
        model._reset_type_to_id_map()

    if args.no_prop_mass:
        skip_properties = {'PSOLID', 'PBUSH', 'PBUSH1D'}
        for pid, prop in model.properties.items():
            if prop.type in {'PSHELL', 'PCOMP', 'PCOMPG'}:
                assert hasattr(prop, 'nsm')
            elif prop.type in {'PBAR', 'PBARL', 'PBEND', 'PBEAM', 'PBEAML'}:
                assert hasattr(prop, 'nsm')
            elif prop.type in skip_properties:
                continue
            else:
                log.warning(f'skipping prop\n{prop.get_stats()}')
    model.cross_reference()

    # nsm_id = 7100000 # ailerons
    mass, cg, inertia = mass_properties_nsm(
        model, element_ids=element_ids, mass_ids=mass_ids, nsm_id=nsm_id, debug=True)
    print(f'mass = {mass}')
    print(f'cg = {cg}')
    print(f'inertia = {inertia}')
    # if bdf_filename_out:
    #     model.write_bdf(bdf_filename_out)
