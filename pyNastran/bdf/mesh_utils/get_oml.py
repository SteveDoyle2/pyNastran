"""
defines:
 - eids_oml = get_oml_eids(bdf_filename, eid_start, theta_tol=30.,
                           is_symmetric=True, consider_flippped_normals=True)

"""
from copy import deepcopy
import numpy as np

from pyNastran.bdf.bdf import read_bdf
#from pyNastran.bdf.bdf_interface.dev_utils import get_free_edges


def get_oml_eids(bdf_filename, eid_start, theta_tol=30.,
                 is_symmetric=True, consider_flippped_normals=True):
    """
    extracts the OML faces (outer mold line)

    Parameters
    ----------
    bdf_filename : str or BDF()
        the bdf filename
    eid_start : int
        the element to start from
    theta_tol : float; default=30.
        the angular tolerance in degrees
    is_symmetric : bool; default=True
        is the y=0 plane considered to be part of the OML
    consider_flippped_normals : bool; default=True
        if you extracted the free faces from tets, you can get flipped normals
        this considers a 180 degree error to be 0.0, which will cause other problems

    """
    #ninety = np.radians(90.)

    #2810 # start for bwb_saero.bdf
    #2811 # close
    #2819 # close
    #2818 # close

    #eids_oml = np.array([eid_start])
    eids_oml = set([eid_start])
    #---------------------------------
    theta_tol = np.radians(theta_tol)

    model = read_bdf(bdf_filename, xref=True)
    maps = model._get_maps(
        eids=None, map_names=None,
        consider_0d=False, consider_0d_rigid=False,
        consider_1d=False, consider_2d=True, consider_3d=False)
    edge_to_eid_map = maps['edge_to_eid_map']
    eid_to_edge_map = maps['eid_to_edge_map']
    nid_to_edge_map = maps['nid_to_edge_map']

    #free_edges = get_free_edges(model, maps=maps)
    #---------------------------------
    normals = {}
    etypes_skipped = set()
    for eid, elem in model.elements.items():
        if elem.type in ['CTRIA3', 'CQUAD4']:
            normals[eid] = elem.Normal()
        else:
            if elem.type in etypes_skipped:
                continue
            model.log.debug('elem.type=%r is not supported' % elem.type)
            etypes_skipped.add(elem.type)

    #eid_starts = eids_oml.tolist()
    eids_next = set([eid_start])
    while eids_next:
        eid_starts = deepcopy(eids_next)
        eids_oml_start = deepcopy(eids_oml)
        model.log.warning(len(eid_starts))
        while eid_starts:
            eid_start = eid_starts.pop()
            normal_start = normals[eid_start]

            # get the next set of edges
            edges = eid_to_edge_map[eid_start]


            #flattened = []
            #for row in matrix:
                #for n in row:
                    #flattened.append(n)
            # flattened = [n for row in matrix for n in row]
            #eids_to_consider = [edge_to_eid_map[edge] for edge in edges]
            list_eids_to_consider = []
            for edge in edges:
                eids_with_edge = edge_to_eid_map[edge]
                list_eids_to_consider += eids_with_edge
            #list_eids_to_consider = set([eid for eid in edge_to_eid_map[edge] for edge in edges])
            #print('list_eids_to_consider =', list_eids_to_consider)
            eids_to_consider = set(list_eids_to_consider)

            # don't do the same element twice; creates an infinite loop if you do
            #eids_to_check = np.setdiff1d(eids_to_consider, eids_oml)
            eids_to_check = eids_to_consider.difference(eids_oml)

            # don't check elements we're checking right now
            #eids_to_check = np.setdiff1d(eids_to_consider, eid_starts)
            eids_to_check = eids_to_consider.difference(eid_starts)

            #print('eids_to_check =', eids_to_check)
            for eid in eids_to_check:
                normal = normals[eid]
                # a o b = a * b * cos(theta)
                # cos(theta) = (a o b)/ (a b); where |a| = 1; |b| = 1
                cos_theta = normal @ normal_start
                theta = np.arccos(cos_theta)
                if theta < theta_tol:
                    eids_next.add(eid)
                    eids_oml.add(eid)
                elif consider_flippped_normals:
                    # handles flipped normals
                    cos_theta = normal @ -normal_start
                    theta = np.arccos(cos_theta)
                    if theta < theta_tol:
                        eids_next.add(eid)
                        eids_oml.add(eid)
            #print('eids_next =', eids_next)
        eids_next = eids_next.difference(eids_oml_start)
        #eids_next = eids_next.difference(eid_starts)
        #print('eids_next =', eids_next)
        #print('-------------------------------')
    model.log.warning('done...')

    with open('eids_oml.txt', 'w') as eids_file:
        eids_file.write('eids_oml = %s\n' % list(eids_oml))
    return eids_oml

def main():
    """runs the test problem"""
    bdf_filename = 'bwb_saero.bdf'
    eid_start = 2810
    eids_oml = get_oml_eids(bdf_filename, eid_start)

if __name__ == '__main__':  # pragma: no cover
    main()
