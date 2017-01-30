"""
Defines:
 - nodes, bars = export_mcids(bdf_filename, csv_filename_out=None)
"""
from __future__ import print_function
from six import string_types, iteritems
import numpy as np
from pyNastran.bdf.bdf import read_bdf


def export_mcids(bdf_filename, csv_filename_out=None):
    """
    Exports the element material coordinates systems, so you can
    load it into pyNastranGUI.

    Parameters
    ----------
    bdf_filename : str/BDF
        a bdf filename or BDF model
    csv_filename_out : str
        the path to the output csv

    TODO: add CTRIA3 support
    """
    if isinstance(bdf_filename, string_types):
        print('bdf_filename =', bdf_filename)
        model = read_bdf(bdf_filename, xref=False)
        print(model.get_bdf_stats())
        model.safe_cross_reference()
    else:
        model = bdf_filename

    skip_types = [
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CGAP',
        'CROD', 'CONROD',
        'CBEAM',
        'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
        'CRAC2D', 'CRAC3D',
    ]

    eid = 1
    nid = 1
    nodes = []
    bars = []
    for eidi, elem in sorted(iteritems(model.elements)):
        if elem.type in ['CQUAD4', 'CQUAD8']:
            xyz1 = elem.nodes_ref[0].get_position()
            xyz2 = elem.nodes_ref[1].get_position()
            xyz3 = elem.nodes_ref[2].get_position()
            xyz4 = elem.nodes_ref[3].get_position()

            dxyz21 = np.linalg.norm(xyz2 - xyz1)
            dxyz32 = np.linalg.norm(xyz3 - xyz2)
            dxyz43 = np.linalg.norm(xyz4 - xyz3)
            dxyz14 = np.linalg.norm(xyz1 - xyz4)
            dxyz = np.mean([dxyz21, dxyz32, dxyz43, dxyz14]) / 2.
            centroid, imat, jmat, _normal = elem.material_coordinate_system()

            iaxis = centroid + imat * dxyz
            jaxis = centroid + jmat * dxyz
            nodes.append((nid, centroid[0], centroid[1], centroid[2]))
            nodes.append((nid + 1, iaxis[0], iaxis[1], iaxis[2]))
            nodes.append((nid + 2, jaxis[0], jaxis[1], jaxis[2]))
            bars.append((eid, nid, nid + 1))
            bars.append((eid + 1, nid, nid + 2))
            nid += 3
            eid += 2
        elif elem.type in skip_types:
            pass
        else:
            raise NotImplementedError('element type=%r is not supported' % elem)

    if csv_filename_out:
        with open(csv_filename_out, 'w') as out_file:
            for node in nodes:
                out_file.write('GRID,%i,%s,%s,%s\n' % node)
            for bari in bars:
                out_file.write('BAR,%i,%i,%i\n' % bari)
    return nodes, bars
