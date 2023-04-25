"""
Defines:
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)

"""
from __future__ import annotations
from typing import Union, Optional, TextIO, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone, read_tecplot
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.result_objects.table_object import TableArray


def tecplot_to_nastran_filename(tecplot_filename: Union[str, Tecplot], bdf_filename: str,
                                log: Optional[SimpleLogger]=None, debug: bool=True) -> BDF:
    """Converts a Tecplot file to Nastran."""
    return tecplot_to_nastran(tecplot_filename, bdf_filename, log=log, debug=debug)


def tecplot_to_nastran(tecplot_filename: Union[str, Tecplot],
                       bdf_filename: str,
                       log: Optional[SimpleLogger]=None,
                       debug: bool=True) -> Optional[BDF]:
    """Converts a Tecplot file to Nastran."""
    if isinstance(tecplot_filename, str):
        model = read_tecplot(tecplot_filename, log=log, debug=debug)
    else:
        model = tecplot_filename
    log = model.log

    removed_nodes = False
    E = 3.0e7
    G = None
    nu = 0.3
    t = 0.1

    inode = 1
    ielem = 1
    nzones = len(model.zones)
    assert nzones > 0, f'nzones={nzones}; tecplot_filename={tecplot_filename!r}'
    with open(bdf_filename, 'w') as bdf_file:
        bdf_file.write('$pyNastran : punch=True\n')
        for izone, zone in enumerate(model.zones):
            mid = izone + 1
            shell_pid = izone + 1
            solid_pid = shell_pid + nzones

            if len(zone.tri_elements) + len(zone.quad_elements):
                card = ['PSHELL', shell_pid, mid, t]
                bdf_file.write(print_card_8(card))

            if len(zone.tet_elements) + len(zone.hexa_elements):
                card = ['PSOLID', solid_pid, mid]
                bdf_file.write(print_card_8(card))

            card = ['MAT1', mid, E, G, nu]
            bdf_file.write(print_card_8(card))

        for izone, zone in enumerate(model.zones):
            inode0 = inode
            shell_pid = izone + 1
            solid_pid = shell_pid + nzones
            _write_loads(bdf_file, zone, inode=inode)

            inode = _write_grids(bdf_file, zone, inode=inode)
            ielem = _write_shells(
                bdf_file, zone, shell_pid, log,
                inode=inode0, ielem=ielem)
            ielem, removed_nodes = _write_solids(
                bdf_file, zone, solid_pid, log,
                inode=inode0, ielem=ielem, removed_nodes=removed_nodes)
            log.debug(f'izone={izone:d}: inode={inode:d} ielem={ielem}')
            # inode += 10
            # ielem += 10

    if removed_nodes:
        log.debug(f'calling remove_unused...')
        bdf_model = BDF(debug=debug)
        bdf_model.read_bdf(bdf_filename)
        remove_unused(bdf_model)
        bdf_model.write_bdf(bdf_filename)
        return bdf_model
    return None

def _write_loads(bdf_file: TextIO, zone: Zone, inode: int=1) -> None:
    """writes a tecplot Zone in Nastran format"""
    for ivar, var in enumerate(zone.variables):
        bdf_file.write(f'$ ivar={ivar+1}: {var!r}\n')

    bdf_file.write(f'$ nodal_results.shape={str(zone.nodal_results.shape)}\n')
    nvars = zone.nodal_results.shape[1]
    for ivar in range(nvars):
        bdf_file.write(f'$ ivar={ivar+1}\n')
        res = zone.nodal_results[:, ivar]
        for inodei, resi in enumerate(res):
            card = ['SLOAD', ivar+1, inode + inodei, resi]
            bdf_file.write(print_card_8(card))

def _write_grids(bdf_file: TextIO, zone: Zone, inode: int=1):
    """writes a tecplot Zone in Nastran format"""
    for inodei, node in enumerate(zone.xyz):
        card = ['GRID', inode + inodei, None,] + list(node)
        bdf_file.write(print_card_8(card))
    inode += inodei + 1
    return inode

def _write_shells(bdf_file: TextIO, zone: Zone, pid: int, log,
                  inode: int=0, ielem: int=1) -> int:
    """writes a tecplot Zone in Nastran format"""
    ielem0 = ielem
    itri = 0
    if len(zone.tri_elements):
        # tris only
        log.debug('tri')
        for itri, tri in enumerate(inode + zone.tri_elements):
            card = ['CTRIA3', ielem + itri, pid] + list(tri)
            bdf_file.write(print_card_8(card))
        ielem += itri + 1

    if len(zone.quad_elements):
        if len(zone.tri_elements) != 0:
            #log.debug('quad-1')
            # if there are tris, then we assume the quads are good
            for iquad, quad in enumerate(inode + zone.quad_elements):
                card = ['CQUAD4', ielem + iquad, pid] + list(quad)
                bdf_file.write(print_card_8(card))
        else:
            #log.debug('quad-2')
            # need to split out the CQUAD4 elements
            # ielem += itri + 1
            for iquad, quad in enumerate(inode + zone.quad_elements):
                if quad[2] == quad[3]:
                    # if it's a tri
                    card = ['CTRIA3', ielem + iquad, pid] + list(quad[:3])
                else:
                    card = ['CQUAD4', ielem + iquad, pid] + list(quad)
                bdf_file.write(print_card_8(card))
        ielem += iquad + 1
    assert ielem >= ielem0, 'ielem={ielem} ielem0={ielem0}'
    return ielem

def _write_solids(bdf_file: TextIO, zone: Zone, pid: int, log,
                  inode: int=0, ielem: int=1,
                  removed_nodes: bool=True) -> tuple[int, bool]:
    """writes a tecplot Zone in Nastran format"""
    if len(zone.tet_elements):
        #log.debug('tet')
        for itet, tet in enumerate(inode + zone.tet_elements):
            card = ['CTETRA', ielem + itet, pid] + list(tet)
            bdf_file.write(print_card_8(card))
        ielem += itet + 1

    if len(zone.hexa_elements):
        #log.debug('hexa')
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
                card = ['CTETRA', ielem + ihex, pid] + list(inode + nids)
                assert len(card) == 7, len(card)
            elif nnodes_unique == 5:
                card = ['CPYRAM', ielem + ihex, pid] + list(inode + nids)
                assert len(card) == 8, len(card)
            elif nnodes_unique == 6:
                card = ['CPENTA', ielem + ihex, pid] + list(inode + nids)
                assert len(card) == 9, len(card)
            elif nnodes_unique == 8:
                card = ['CHEXA', ielem + ihex, pid] + list(inode + hexa)
            bdf_file.write(print_card_8(card))
        ielem += ihex + 1
    return ielem, removed_nodes


def nastran_table_to_tecplot(bdf_model: BDF,
                             case: TableArray,
                             variables: list[str]) -> Tecplot:
    """assumes only triangles"""
    title = ('%s; %s' % (case.title, case.subtitle)).strip(' ;')
    ntimes, nnodes, nresults = case.data.shape

    xyz_list = []
    tris = []
    nid_map = {}
    #all_nids = list(bdf_model.nodes.keys())
    for inid, (nid, node) in enumerate(sorted(bdf_model.nodes.items())):
        xyz_list.append(node.get_position())
        nid_map[nid] = inid
    for eid, elem in sorted(bdf_model.elements.items()):
        if elem.type == 'CQUAD4':
            node_ids = elem.node_ids
            tri_nids = [node_ids[:3], node_ids[1:]]
        elif elem.type == 'CTRIA3':
            tri_nids = [elem.node_ids]
        else:
            raise NotImplementedError(elem)

        for nids in tri_nids:
            assert len(nids) == 3, tri_nids
            new_tri = []
            for nid in nids:
                new_tri.append(nid_map[nid])
            #except:
                #raise
            #new_tri = [nid_map[nid] for nid in nids]
            tris.append(new_tri)

    tecplot_model = Tecplot(log=bdf_model.log, debug=bdf_model.debug)
    zone = Zone(bdf_model.log)
    xyz = np.array(xyz_list, dtype='float64')

    results = np.full((nnodes, nresults), np.nan, dtype=xyz.dtype)
    zone.zone_data = np.hstack([xyz, results])
    zone.tri_elements = np.array(tris, dtype='int32') + 1

    tecplot_model.title = title
    zone.headers_dict['VARIABLES'] = variables
    zone.variables = variables
    tecplot_model.zones = [zone]
    return tecplot_model

def nastran_tables_to_tecplot_filenames(tecplot_filename_base: str,
                                        bdf_model: BDF,
                                        case: TableArray,
                                        variables: Optional[list[str]]=None,
                                        ivars: Optional[list[int]]=None) -> None:
    """case is intended as an op2 object like DisplacementArray"""
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
        zone.zone_data[:, 3:] = nodal_results[:, ivars]
        tecplot_model.write_tecplot(
            tecplot_filename, res_types=None, adjust_nids=False)

if __name__ == '__main__':  # pragma: no cover
    import sys
    tecplot_filename = sys.argv[1]
    bdf_filename = sys.argv[2]
    tecplot_to_nastran_filename(tecplot_filename, bdf_filename, log=None, debug=True)
