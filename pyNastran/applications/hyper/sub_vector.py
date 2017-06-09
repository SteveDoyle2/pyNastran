from six import iteritems
from numpy import hstack, repeat, sqrt, triu, arange

from pyNastran.dev.bdf_vectorized.bdf import BDF


def main():
    import os
    bdf_filename = 'threePlugs.bdf'
    if not os.path.exists(bdf_filename):
        import pyNastran
        from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename
        pkg_path = pyNastran.__path__[0]
        cart3d_filename = os.path.join(pkg_path, 'converters', 'cart3d', 'models', 'threePlugs.bin.tri')
        cart3d_to_nastran_filename(cart3d_filename, bdf_filename)
    run(bdf_filename)

def run(bdf_filename, debug=True):
    # http://www.3dpanelmethod.com/documents/Graduate%20Work.pdf
    model = BDF(debug=debug)
    model.read_bdf(bdf_filename)


    grids = model.grid
    #print(list(grids.node_id))
    xyz_global = grids.get_position_by_node_index()

    tris = model.elements.elements_shell.ctria3
    quads = model.elements.elements_shell.cquad4

    if tris.n:
        et = tris.element_id
        pt = tris.property_id
        At = tris.get_area_by_element_index()
        nt = tris.get_normal_by_element_index()
        ct = tris.get_centroid_by_element_index()

        #i = arange(tris.n)
        #n1, n2, n3 = tris._node_locations(xyz_cid0, i)
        # n3 = n4
        # is this right?
        #ut = (n1 + n2 - 2 * n3) / 2.
        #pt = (n2 - n1) / 2.
        #ot = cross(nt, ut)
    if quads.n:
        eq = quads.element_id
        pq = quads.property_id
        Aq = quads.get_area_by_element_index()
        nq = quads.get_normal_by_element_index()
        cq = quads.get_centroid_by_element_index()

        i = arange(quads.n)
        #n1, n2, n3, n4 = quads._node_locations(xyz_cid0, i)
        #uq = (n1 + n2 - n3 - n4) / 2.
        #pq = (n2 + n3 - n4 - n1) / 2.
        #oq = cross(nq, uq)

    if tris.n and quads.n:
        e = hstack([et, eq])
        p = hstack([pt, pq])
        A = hstack([At, Aq])
        n = hstack([nt, nq])
        c = hstack([ct, cq])
        o = hstack([ot, oq])
    elif tris.n:
        e = et
        p = pt
        A = At
        n = nt
        c = ct
        #o = ot
    elif quads.n:
        e = eq
        p = pq
        A = Aq
        n = nq
        c = cq
        #o = oq

    bcs = {}
    for key, set3i in iteritems(model.set3):
        if set3i.desc == 'ELEM':
            bcs[key] = set3i.IDs

    #A = array([1, 2, 3, 4, 5], dtype='float32')
    nelements = len(e)
    #print('c.shape', c.shape)

    # could we take advantage of upper triangular matrix?
    # it's a factor of 2 speedup, which is pretty minor relative to the
    # added code complexity
    #
    # we split the x, y, z components to make it easier to run our calculations
    csquarex = repeat(c[:, 0], nelements).reshape(nelements, nelements)
    csquarey = repeat(c[:, 1], nelements).reshape(nelements, nelements)
    csquarez = repeat(c[:, 2], nelements).reshape(nelements, nelements)

    csquarex -= c[:, 0]
    csquarey -= c[:, 1]
    csquarez -= c[:, 2]

    # 2-norm
    dist = sqrt(csquarex**2 + csquarey**2 + csquarez**2)

    print('dist', dist, dist.shape)
    #print(csquarex)
    #print("nelements", nelements)
    #print(csquarex.shape)

def run2(model):
    # quad - constant strength source
    k = 1. / (4. * pi)

    elements = model.elements
    n1 = elements.cquad4.nodes[:, 0]
    n2 = elements.cquad4.nodes[:, 1]
    n3 = elements.cquad4.nodes[:, 2]
    n4 = elements.cquad4.nodes[:, 3]

    x1 = xyz_cid0[n1, 0]
    x2 = xyz_cid0[n2, 0]
    x3 = xyz_cid0[n3, 0]
    x4 = xyz_cid0[n4, 0]

    y1 = xyz_cid0[n1, 1]
    y2 = xyz_cid0[n2, 1]
    y3 = xyz_cid0[n3, 1]
    y4 = xyz_cid0[n4, 1]

    d12 = sqrt((x2-x1)**2 + (y2-y1)**2)
    d23 = sqrt((x3-x2)**2 + (y3-y2)**2)
    d34 = sqrt((x4-x3)**2 + (y4-y3)**2)
    d41 = sqrt((x1-x4)**2 + (y1-y4)**2)

    m12 = (y2 - y1) / (x2 - x1)
    m23 = (y3 - y2) / (x3 - x2)
    m34 = (y4 - y3) / (x4 - x3)
    m41 = (y1 - y4) / (x1 - x4)


    r = sqrt( (x-x1)**2 + (y-y1)**2 + z**2)
    e1 = (x - x1)**2 + z**2
    e2 = (x - x2)**2 + z**2
    e3 = (x - x3)**2 + z**2
    e4 = (x - x4)**2 + z**2

    h1 = (x - x1) * (y - y1)
    h2 = (x - x2) * (y - y2)
    h3 = (x - x3) * (y - y3)
    h4 = (x - x4) * (y - y4)

    phi = k * (
      (
         (x-x1)*(y2-y1)-(y-y1)*(x2-x1)/d12 * ln((r1+r2+d12)/(r1+r2-d12))
        +(x-x2)*(y3-y2)-(y-y2)*(x3-x2)/d23 * ln((r2+r3+d23)/(r2+r3-d23))
        +(x-x3)*(y4-y3)-(y-y3)*(x4-x3)/d34 * ln((r3+r4+d34)/(r3+r4-d34))
        +(x-x4)*(y1-y4)-(y-y4)*(x1-x4)/d41 * ln((r4+r1+d41)/(r3+r4-d34))
      )
      + abs(z) * (
           atan2(m12*e1-h1, z*r1) - atan2(m12*e2-h2, z*r2)
          +atan2(m23*e2-h2, z*r2) - atan2(m23*e3-h3, z*r3)
          +atan2(m34*e3-h3, z*r3) - atan2(m34*e4-h4, z*r4)
          +atan2(m41*e4-h4, z*r4) - atan2(m41*e1-h1, z*r1)
        )
      )

    Cp = 1 - 1 / Vref**2 * (Q**2 - 2.0 * dphi_dt)

if __name__ == '__main__':
    main()
