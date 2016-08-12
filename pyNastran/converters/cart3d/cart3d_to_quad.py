from __future__ import print_function
from collections import defaultdict
from numpy import cross, allclose
from numpy.linalg import norm
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.bdf.field_writer_8 import print_card_8

def main():
    """tests getting the normal groups"""
    cart3d = Cart3D(log=None, debug=False)
    result_names = []  # read the mesh only
    cart3d.read_cart3d('Cart3d_55000_2.4_10_0_0_0_0.i.triq', result_names=result_names)
    points = cart3d.points
    elements = cart3d.elements

    celements = elements.copy()
    normals, groups = get_normal_groups(points, celements)
    tris, quads = normal_groups_to_quads(celements, normals, groups)

    write_nastran_quads_tris(points, tris, quads, bdf_filename='tris_quads.bdf')

def normal_groups_to_quads(elements, normals, normal_groups):
    tris = []
    quads = []
    for group in normal_groups:
        if len(group) == 2:
            g1, g2 = group
            n1 = elements[g1] + 1
            for ni in elements[g2]:
                if ni not in n1:
                    break
            quad = [n1[0], n1[1], ni + 1, n1[2]]
            quads.append(quad)
        elif len(group) == 1:
            # null case, one element
            for groupi in group:
                tri = elements[groupi, :] + 1
                #print("t1", tri)
                tris.append(tri)
        else:
            # this is more complicated, so we'll skip it for now
            for groupi in group:
                tri = elements[groupi, :] + 1
                #print("t2", tri)
                tris.append(tri)

    #for quad in quads:
        #print(quad)
    #for tri in tris:
        #print("t", tri)
    return tris, quads

def write_nastran_quads_tris(nodes, tris, quads, bdf_filename='tris_quads.bdf'):
    with open(bdf_filename, 'wb') as bdf_file:
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')
        cp = 0
        for inode, node in enumerate(nodes):
            card = ['GRID', inode + 1, cp, ] + list(node)
            bdf_file.write(print_card_8(card))

        eid = 1
        pid = 1
        mid = 1
        for tri in tris:
            #print(tri)
            card = ['CTRIA3', eid, pid, ] + list(tri)
            #print(card)
            bdf_file.write(print_card_8(card))
            eid += 1

        for quad in quads:
            #print("quad[%s] = %s" % (eid, str(quad)))
            card = ['CQUAD4', eid, pid, ] + list(quad)
            bdf_file.write(print_card_8(card))
            eid += 1

        t = 0.1
        card = ['PSHELL', pid, mid, t]
        bdf_file.write(print_card_8(card))

        E = 1e7
        G = None
        nu = 0.3
        card = ['MAT1', mid, E, G, nu]
        bdf_file.write(print_card_8(card))

        bdf_file.write('ENDDATA\n')

def get_normal_groups(points, elements, rtol=1e-3, atol=1e-5):
    """
    Gets the normal groups for a cart3d model

    Parameters
    ----------
    points : ndarray
        the points from Cart3d()
    elements : ndarray
        the elements from Cart3d()

    Returns
    -------
    normal_groups : List[int, int, ...]
        the list of list of ints of the element IDs
        with the same normal that are touching

    .. note:: modifies elements
    """
    elements -= 1
    #npoints = points.shape[0]
    nelements = elements.shape[0]

    #print(points)
    #print(elements)

    n1 = elements[:, 0]

    #print(len(unique(n1)))
    #print(list(n1), max(n1), min(n1))
    p1 = points[n1, :]
    p2 = points[elements[:, 1]]
    p3 = points[elements[:, 2]]

    a = p2 - p1
    b = p3 - p1
    n = cross(a, b)
    ni = norm(n, axis=1)
    #print(n.shape)
    #print(ni.shape)

    normals = (n.T / ni).T
    #normals = n / norm(n, axis=1)  # this should work...but it doesn't

    #connectivity = zeros((n, 3), 'int32')
    cdict = defaultdict(list)

    # get edge mapping
    eid = 0
    for element in elements:
        c1 = tuple(sorted([element[0], element[1]]))
        c2 = tuple(sorted([element[1], element[2]]))
        c3 = tuple(sorted([element[0], element[2]]))
        #if eid == 0:
            #print(c1, c2, c3)
        cdict[c1].append(eid)
        cdict[c2].append(eid)
        cdict[c3].append(eid)
        eid += 1

    # sort it

    unique_normals = set([])
    for normal in normals:
        unique_normals.add(tuple(normal))

    #for normal
    # group elements by normal
    #normal_groups = defaultdict(list)
    #igroup = 0


    # group elements with same normal and connectivity
    #for normali in unique_normals:
        #i = where(normal == normali)[0]

    # ???
    groups = []
    set_elements = set([])
    for eid in range(nelements):
        #normali = normal[eid, :]
        if eid not in set_elements:
            same_normals = [eid]
            check_normals(eid, elements, normals, cdict, same_normals,
                          rtol=1e-4, atol=1e-3)  # same_normals is modified
            #print(same_normals, '\n')
            groups.append(same_normals)
            set_elements.update(same_normals)
        #break

    #for group in groups:
        #if len(group) > 1:
            #print(group)
    return normals, groups



def check_normals(eid, elements, normals, cdict, same_normals, rtol=1e-4, atol=1e-3):
    (n1, n2, n3) = elements[eid, :]
    c1 = tuple(sorted([n1, n2]))
    c2 = tuple(sorted([n2, n3]))
    c3 = tuple(sorted([n3, n1]))

    # remove the current element as it's been checked
    eids_alt = cdict[c1]
    #print("eids_alt", eids_alt)
    if eid in eids_alt:
        eids_alt.pop(eids_alt.index(eid))
        if eids_alt:
            eid_alt = eids_alt[0]
            if allclose(normals[eid_alt], normals[eid], rtol=rtol, atol=atol):
                # same normal, adjacent
                #print("  same", eid, eid_alt, list(normals[eid]), list(normals[eid_alt]))
                same_normals.append(eid_alt)
                check_normals(eid_alt, elements, normals, cdict, same_normals, rtol=rtol, atol=atol)

    eids_alt = cdict[c2]
    if eid in eids_alt:
        eids_alt.pop(eids_alt.index(eid))
        if eids_alt:
            eid_alt = eids_alt[0]
            if allclose(normals[eid_alt], normals[eid], rtol=rtol, atol=atol):
                # same normal, adjacent
                #print("  same", eid, eid_alt)
                same_normals.append(eid_alt)
                check_normals(eid_alt, elements, normals, cdict, same_normals, rtol=rtol, atol=atol)

    eids_alt = cdict[c3]
    if eid in eids_alt:
        eids_alt.pop(eids_alt.index(eid))
        if eids_alt:
            eid_alt = eids_alt[0]
            if allclose(normals[eid_alt], normals[eid], rtol=rtol, atol=atol):
                # same normal, adjacent
                #print("  same", eid, eid_alt)
                same_normals.append(eid_alt)
                check_normals(eid_alt, elements, normals, cdict, same_normals, rtol=rtol, atol=atol)


main()
