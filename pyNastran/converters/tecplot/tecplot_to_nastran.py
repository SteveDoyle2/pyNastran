"""
Defines:
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)

"""
from __future__ import annotations
from typing import Union, List, Optional, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone, read_tecplot
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger


def tecplot_to_nastran_filename(tecplot_filename: Union[str, Tecplot], bdf_filename: str,
                                log: Optional[SimpleLogger]=None, debug: bool=True) -> BDF:
    """Converts a Tecplot file to Nastran."""
    return tecplot_to_nastran(tecplot_filename, bdf_filename, log=log, debug=debug)


def tecplot_to_nastran(tecplot_filename: Union[str, Tecplot], bdf_filename: str,
                       log: Optional[SimpleLogger]=None, debug: bool=True) -> None:
    """Converts a Tecplot file to Nastran."""
    if isinstance(tecplot_filename, str):
        model = read_tecplot(tecplot_filename, log=log, debug=debug)
    else:
        model = tecplot_filename

    zone = model.zones[0]
    removed_nodes = False
    shell_pid = 1
    solid_pid = 2
    mid = 1
    istart = 1
    with open(bdf_filename, 'w') as bdf_file:
        bdf_file.write('$pyNastran : punch=True\n')
        for inode, node in enumerate(zone.xyz):
            card = ['GRID', inode + 1, None,] + list(node)
            bdf_file.write(print_card_8(card))

        itri = 0
        if len(zone.tri_elements):
            # tris only
            for itri, tri in enumerate(zone.tri_elements):
                card = ['CTRIA3', itri + 1, shell_pid] + list(tri)
                bdf_file.write(print_card_8(card))
            #istart += bdf_model

        if len(zone.quad_elements):
            if len(zone.tri_elements) != 0:
                # if there are tris, then we assume the quads are good
                for iquad, quad in enumerate(zone.quad_elements):
                    card = ['CQUAD4', iquad + 1, shell_pid] + list(quad)
                    bdf_file.write(print_card_8(card))
            else:
                # need to split out the CQUAD4 elements
                istart = itri + 1
                for iquad, quad in enumerate(zone.quad_elements):
                    if quad[2] == quad[3]:
                        # if it's a tri
                        card = ['CTRIA3', istart + iquad, shell_pid] + list(quad[:3])
                    else:
                        card = ['CQUAD4', istart + iquad, shell_pid] + list(quad)
                    bdf_file.write(print_card_8(card))
            istart += iquad

        if len(zone.tri_elements) + len(zone.quad_elements):
            card = ['PSHELL', shell_pid, mid, 0.1]
            bdf_file.write(print_card_8(card))

        if len(zone.tet_elements) + len(zone.hexa_elements):
            card = ['PSOLID', solid_pid, mid]
            bdf_file.write(print_card_8(card))

        if len(zone.tet_elements):
            for itet, tet in enumerate(zone.tet_elements):
                card = ['CTETRA', istart + itet, solid_pid] + list(tet)
                bdf_file.write(print_card_8(card))

        if len(zone.hexa_elements):
            # need to split out the CTETRA and CPENTA elements
            for ihex, hexa in enumerate(zone.hexa_elements):
                uhexa = np.unique(hexa)
                nnodes_unique = len(uhexa)
                nids = hexa[:nnodes_unique]
                centroid_y = zone.xyz[nids, 1].max()
                if centroid_y < 0:
                    removed_nodes = True
                    continue
                if nnodes_unique == 4:
                    card = ['CTETRA', istart + ihex, solid_pid] + list(nids)
                    assert len(card) == 7, len(card)
                elif nnodes_unique == 5:
                    card = ['CPYRAM', istart + ihex, solid_pid] + list(nids)
                    assert len(card) == 8, len(card)
                elif nnodes_unique == 6:
                    card = ['CPENTA', istart + ihex, solid_pid] + list(nids)
                    assert len(card) == 9, len(card)
                elif nnodes_unique == 8:
                    card = ['CHEXA', istart + ihex, solid_pid] + list(hexa)
                bdf_file.write(print_card_8(card))

        E = 3.0e7
        G = None
        nu = 0.3
        card = ['MAT1', mid, E, G, nu]
        bdf_file.write(print_card_8(card))

    if removed_nodes:
        bdf_model = BDF(debug=debug)
        bdf_model.read_bdf(bdf_filename)
        remove_unused(bdf_model)


def nastran_table_to_tecplot(bdf_model: BDF, case, variables: List[str]) -> Tecplot:
    """assumes only triangles"""
    xyz = []
    tris = []
    nid_map = {}
    for inid, (nid, node) in enumerate(sorted(bdf_model.nodes.items())):
        xyz.append(node.get_position())
        nid_map[nid] = inid
    for eid, elem in sorted(bdf_model.elements.items()):
        tris.append([nid_map[nid] for nid in elem.node_ids])

    tecplot_model = Tecplot(log=bdf_model.log, debug=bdf_model.debug)
    zone = Zone(bdf_model.log)
    zone.xyz = np.array(xyz, dtype='float64')
    zone.tri_elements = tris = np.array(tris, dtype='int32') + 1

    tecplot_model.title = ('%s; %s' % (case.title, case.subtitle)).strip(' ;')
    zone.headers_dict['VARIABLES'] = variables
    zone.variables = variables
    tecplot_model.zones = [zone]
    return tecplot_model

def nastran_tables_to_tecplot_filenames(tecplot_filename_base: str,
                                        bdf_model: BDF, case,
                                        variables: Optional[List[str]]=None,
                                        ivars: Optional[List[int]]=None) -> None:
    if variables is None:
        variables = case.headers
    if ivars is None:
        ivars = np.arange(0, len(variables))

    tecplot_model = nastran_table_to_tecplot(bdf_model, case, variables)
    zone = tecplot_model.zones[0]
    for itime, time in enumerate(case._times):
        if '%' in tecplot_filename_base:
            tecplot_filename = tecplot_filename_base % time
        else:
            tecplot_filename = tecplot_filename_base

        # you can't combine the two lines or it transposes it...
        nodal_results = case.data[itime, :, :]
        zone.nodal_results = nodal_results[:, ivars]
        tecplot_model.write_tecplot(
            tecplot_filename, res_types=None, adjust_nids=False)
