"""
Defines:
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)
 - tecplot_to_nastran(tecplot_filename, bdf_filename, debug=True)
"""
from __future__ import print_function
from six import string_types
from numpy import unique
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.converters.tecplot.tecplot import read_tecplot


def tecplot_to_nastran_filename(tecplot_filename, bdf_filename, log=None, debug=True):
    """Converts a Tecplot file to Nastran."""
    return tecplot_to_nastran(tecplot_filename, bdf_filename, log=log, debug=debug)


def tecplot_to_nastran(tecplot_filename, bdf_filename, log=None, debug=True):
    """Converts a Tecplot file to Nastran."""
    if isinstance(tecplot_filename, string_types):
        model = read_tecplot(tecplot_filename, log=log, debug=debug)
    else:
        model = tecplot_filename

    removed_nodes = False
    shell_pid = 1
    solid_pid = 2
    mid = 1
    istart = 1
    with open(bdf_filename, 'w') as bdf_file:
        bdf_file.write('$pyNastran : punch=True\n')
        for inode, node in enumerate(model.xyz):
            card = ['GRID', inode + 1, None,] + list(node)
            bdf_file.write(print_card_8(card))

        itri = 0
        if len(model.tri_elements):
            # tris only
            for itri, tri in enumerate(model.tri_elements):
                card = ['CTRIA3', itri + 1, shell_pid] + list(tri)
                bdf_file.write(print_card_8(card))
            #istart += bdf_model

        if len(model.quad_elements):
            if len(model.tri_elements) != 0:
                # if there are tris, then we assume the quads are good
                for iquad, quad in enumerate(model.quad_elements):
                    card = ['CQUAD4', iquad + 1, shell_pid] + list(quad)
                    bdf_file.write(print_card_8(card))
            else:
                # need to split out the CQUAD4 elements
                istart = itri + 1
                for iquad, quad in enumerate(model.quad_elements):
                    if quad[2] == quad[3]:
                        # if it's a tri
                        card = ['CTRIA3', istart + iquad, shell_pid] + list(quad[:3])
                    else:
                        card = ['CQUAD4', istart + iquad, shell_pid] + list(quad)
                    bdf_file.write(print_card_8(card))
            istart += iquad

        if len(model.tri_elements) + len(model.quad_elements):
            card = ['PSHELL', shell_pid, mid, 0.1]
            bdf_file.write(print_card_8(card))

        if len(model.tet_elements) + len(model.hexa_elements):
            card = ['PSOLID', solid_pid, mid]
            bdf_file.write(print_card_8(card))

        if len(model.tet_elements):
            for itet, tet in enumerate(model.tet_elements):
                card = ['CTETRA', istart + itet, solid_pid] + list(tet)
                bdf_file.write(print_card_8(card))

        if len(model.hexa_elements):
            # need to split out the CTETRA and CPENTA elements
            for ihex, hexa in enumerate(model.hexa_elements):
                uhexa = unique(hexa)
                nnodes_unique = len(uhexa)
                nids = hexa[:nnodes_unique]
                centroid_y = model.xyz[nids, 1].max()
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
