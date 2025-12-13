from __future__ import annotations
import sys
from typing import TYPE_CHECKING
from cpylog import SimpleLogger
import numpy as np

from .utils import add_argparse_arguments
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_super(argv=None, quiet: bool=False) -> None:
    """command line interface to ``bdf super``"""
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
    parent_parser.add_argument('super', type=str)
    parent_parser.add_argument('INPUT_MAIN', help='path to main BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('INPUT_PLOTEL', help='path to plotel BDF/DAT/NAS file (as CQUAD4, CTRIA3, CBAR, etc.)', type=str)
    parent_parser.add_argument('INPUT_RETAIN', help='path to retained BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('nnode', help='starting node id', type=int)
    #parent_parser.add_argument('--output', help='path to output BDF/DAT/NAS file; default=super.bdf', type=str, default='super.bdf')

    #parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)
    parent_parser.add_argument('--plotel', action='store_true', help='make a plotel model')
    parent_parser.add_argument('--nmode', help='set the number of modes', type=int, default=0)
    add_argparse_arguments(parent_parser, ['--punch', '--lax', '--allow_dup'])
    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    log = SimpleLogger(level='debug')

    main_bdf_filename = args.INPUT_MAIN
    plotel_bdf_filename = args.INPUT_PLOTEL
    retained_bdf_filename = args.INPUT_RETAIN
    out_bdf_filename = args.output

    assert plotel_bdf_filename != main_bdf_filename
    assert plotel_bdf_filename != retained_bdf_filename
    assert plotel_bdf_filename != out_bdf_filename

    is_strict_card_parser = not args.lax
    nnode = args.nnode
    nmode = args.nmode
    punch = args.punch
    assert nmode > 0, nmode

    is_strict_card_parser = not args.lax
    from .utils_bdf import read_lax_bdf
    main_model = read_lax_bdf(
        main_bdf_filename, punch=punch, xref=False,
        is_strict_card_parser=is_strict_card_parser,
        log=log)
    plotel_model = read_lax_bdf(
        plotel_bdf_filename, punch=punch, xref=False,
        is_strict_card_parser=is_strict_card_parser,
        log=log)
    retained_model = read_lax_bdf(
        retained_bdf_filename, punch=punch, xref=False,
        is_strict_card_parser=is_strict_card_parser,
        log=log)

    # matrices_to_add = {
    #     'KAAX',
    #     'MAAX',
    #     'BAAX',
    #     'K4AAX',
    # }
    #-------------------------------------------------
    subcases = main_model.subcases
    subcase0 = subcases[0]
    subcase0.add_integer_type('SUPER', 103)
    subcase0.add_set_from_values(1001, ['KAAX'])
    subcase0.add_set_from_values(1002, ['MAAX'])
    subcase0.add_set_from_values(1003, ['BAAX'])
    subcase0.add_set_from_values(1004, ['K4AAX'])

    subcase0.add_integer_type('K2GG', 1001)
    subcase0.add_integer_type('M2GG', 1002)
    subcase0.add_integer_type('B2GG', 1003)
    subcase0.add_integer_type('K42GG', 1004)

    #  remove elements that aren't retained
    eids_to_remove = []
    rigid_eids_to_remove = []
    mass_eids_to_remove = []
    for eid, elem in main_model.elements.items():
        if eid not in retained_model.elements:
            eids_to_remove.append(eid)
    for eid, elem in main_model.rigid_elements.items():
        if eid not in retained_model.rigid_elements:
            rigid_eids_to_remove.append(eid)
    for eid, elem in main_model.masses.items():
        if eid not in retained_model.masses:
            mass_eids_to_remove.append(eid)

    for eid in eids_to_remove:
        del main_model.elements[eid]
    for eid in rigid_eids_to_remove:
        del main_model.rigid_elements[eid]
    for eid in mass_eids_to_remove:
        del main_model.masses[eid]

    #  remove nodes that aren't retained
    retained_elements_to_skip = set()
    retained_nodes = set()
    for eid, elem in retained_model.elements.items():
        if elem.type in retained_elements_to_skip:
            continue
        nodesi = elem.nodes
        retained_nodes.update(nodesi)
    for eid, elem in retained_model.masses.items():
        if elem.type in retained_elements_to_skip:
            continue
        nodesi = elem.nodes
        retained_nodes.update(nodesi)
    for eid, elem in retained_model.rigid_elements.items():
        if elem.type in retained_elements_to_skip:
            continue
        nodesi = elem.nodes
        retained_nodes.update(nodesi)

    nids_to_remove = set(list(main_model.nodes)) - retained_nodes
    for nid in nids_to_remove:
        del main_model.nodes[nid]
    #-------------------------------------------------

    max_nid_main = max(main_model.nodes)
    max_nid_retained = max(retained_model.nodes) if len(retained_model.nodes) else 0
    max_nid = max(max_nid_main, max(plotel_model.nodes), max_nid_retained)
    assert nnode > max_nid, f'nnode={nnode:d}, max_nid={max_nid} is too large'

    #model.cross_reference()
    print(plotel_model.get_bdf_stats())
    #print(retained_model.get_bdf_stats())

    #-------------------------------------------------
    from pyNastran.bdf.bdf import BDF
    super_model = BDF(log=log)

    skipped_etypes = {'CBUSH'}
    comment = f'$ plot elements per {plotel_bdf_filename}'
    for eid, elem in plotel_model.elements.items():
        etype = elem.type
        nodes = elem.nodes

        if etype in skipped_etypes:
            continue

        if etype == 'CQUAD4':
            super_model.add_plotel4(eid, nodes, comment=comment)
        elif etype == 'CTRIA3':
            super_model.add_plotel3(eid, nodes, comment=comment)
        elif etype in {'CROD', 'CBAR', 'CBEAM', 'CBEND'}:
            assert len(nodes) == 2, elem
            super_model.add_plotel(eid, nodes, comment=comment)
        elif elem.type == 'CTETRA':
            super_model.add_plottet(eid, nodes, comment=comment)
        elif elem.type == 'CPENTA':
            super_model.add_plotpen(eid, nodes, comment=comment)
        elif elem.type == 'CPYRAM':
            super_model.add_plotpyr(eid, nodes, comment=comment)
        elif elem.type == 'CHEXA':
            super_model.add_plothex(eid, nodes, comment=comment)
        else:
            log.warning(elem.rstrip())
            continue
        comment = ''
    super_model.write_bdf('super_plotel.bdf')
    #-------------------------------------------------
    retain_setup_model = BDF(log=log)

    spoint_ids = np.arange(1, nmode+1) + nnode
    assert len(spoint_ids) == nmode and nmode > 0, f'spoint_ids={spoint_ids}; len={len(spoint_ids)}; nmode={nmode}'
    retain_setup_model.add_spoint(spoint_ids, f' per nmodes={nmode}')

    qset = retain_setup_model.add_qset1(spoint_ids, '0', comment=' per spoint modes')

    fixed_nids = list(set(list(retain_setup_model.nodes)) | set(list(plotel_model.nodes)))
    fixed_nids.sort()
    fixed_components = '123456'
    assert len(fixed_nids) > 0, fixed_nids

    nfree_nodes = len(fixed_nids)
    ndof = 6 * nfree_nodes
    cset = retain_setup_model.add_cset1(
        fixed_nids, fixed_components,
        comment=f' nfree_nodes={nfree_nodes}; ndof={ndof}; free nodes per {retained_bdf_filename}')  # fixed
    #print(cset)
    #print(qset)
    #print(spoint_ids)
    #super_model.write_bdf(out_bdf_filename)
    retain_setup_model.write_bdf('super_setup.bdf')
    main_model.write_bdf('super_run.bdf')
