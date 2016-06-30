from copy import deepcopy
from six.moves import range

from numpy import zeros

from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class Cart3d_Mesher(Cart3D):
    def __init__(self, log=None, debug=False):
        Cart3D.__init__(self, log=log, debug=debug)

    def read_cart3d(self, cart3d_filename):
        Cart3D.read_cart3d(self, cart3d_filename)
        #self.nodes = nodes
        #self.elements = elements - 1

    def _get_segment(self, a, eid, segments):
        if a in segments:
            a_elems = segments[a]
            i = a_elems.index(eid)
            a_elems.pop(i)
            eid_a = a_elems[0]
            return eid_a
        return None

    def _get_segments(self, nodes, elements):
        segments = {}  # key=eid,
        lengths = {}
        for eid, e in iteritems(elements):
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))

            if a in segments:
                segments[a].append(eid)
            else:
                segments[a] = [eid]
                lengths[a] = nodes[a[1]] - nodes[a[0]]

            if b in segments:
                segments[b].append(eid)
            else:
                segments[b] = [eid]
                lengths[b] = nodes[b[1]] - nodes[b[0]]

            if c in segments:
                segments[c].append(eid)
            else:
                segments[c] = [eid]
                lengths[c] = nodes[c[1]] - nodes[c[0]]
        return segments, lengths

    def make_quads(self, nodes, elements):
        raise NotImplementedError()
        segments, lengths = self._get_segments(nodes, elements)
        for eid, e in iteritems(elements):
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))
            #a.sort()
            #b.sort()
            #c.sort()
            print(eid, e)
            print(segments[a])
            print(lengths[a])
            print(len(segments[a]))

            eidA = self._get_segment(a, eid, segments)
            eidB = self._get_segment(b, eid, segments)
            eidC = self._get_segment(c, eid, segments)
            print("eidA=%s eidB=%s eidC=%s" % (eidA, eidB, eidC))
            if eidA:
                i = 0
                e2 = elements[eidA]
                self.check_quad(nodes, eid, eidA, e, e2, a, b, c, i)
                del segments[a]
            if eidB:
                i = 1
                e2 = elements[eidB]
                self.check_quad(nodes, eid, eidB, e, e2, a, b, c, i)
                del segments[b]
            if eidC:
                i = 2
                e2 = elements[eidC]
                self.check_quad(nodes, eid, eidC, e, e2, a, b, c, i)
                del segments[c]

            print("------")
            #break
        #for segment in segments:
        asdf

    def _check_quad(self, nodes, eid, eidA, e, e2, a, b, c, i):
        r"""
        ::
          A----B
          | \ e|
          |e2 \|
          C----D

        two tests
           1.  folding angle A-B x A-C
           2a. abs(A-C) - abs(B-D)  = 0  (abs to prevent 2L)
           2b. abs(A-B) - abs(C-D)  = 0
        """

        iplus1 = i + 1
        iplus2 = i + 2

        if iplus1 > 2:
            iplus1 -= 3
        if iplus2 > 2:
            iplus2 -= 3
        print(i, iplus1)
        print(iplus1, iplus2)
        print(iplus2, i)
        AD = nodes[e[i]] - nodes[e[iplus1]]
        AB = nodes[e[iplus1]] - nodes[e[iplus2]]
        BD = nodes[e[iplus2]] - nodes[e[i]]

        print(AD)
        print(AB)
        print(BD)
        print(e2)
        j = e2.index(e[i])

        jplus1 = j + 1
        jplus2 = j + 2
        if jplus1 > 2:
            jplus1 -= 3
        if jplus2 > 2:
            jplus2 -= 3

        print("DA = ", e[j], e[jplus1])
        DA = nodes[e[j]] - nodes[e[jplus1]]
        print(DA)

        asdf

    def project(self, bdf_filename, x0, growth_rate=1.3, nlayers=10):
        x = zeros(nlayers, dtype='float64')
        for i in range(nlayers):
            x[i] = x0 * growth_rate ** i

        print(self.nodes.shape)
        print(self.elements.shape)

        nnodes = self.nodes.shape[0]
        nelements = self.elements.shape[0]

        nnodes2 = nnodes * (nlayers + 1)
        npents = nelements * nlayers

        cnormals = self.get_normals(self.nodes, self.elements)
        nnormals = self.get_normals_at_nodes(self.nodes, self.elements, cnormals)

        nodes = zeros((nnodes2, 3), dtype='float64')
        pents = zeros((npents, 6), dtype='int32')
        ih1 = 0

        in1 = 0
        in2 = nnodes
        nodes[in1:in2, :] = self.nodes
        in1 += nnodes
        in2 += nnodes

        ih1 = 0
        ih2 = nelements
        print('x = %s' % x)
        for i in range(nlayers):
            nodes_old = self.nodes
            dx = nnormals * x[i]
            nodes_new = self.nodes + dx
            dn0 = nnodes * i
            dn1 = nnodes * (i + 1)
            elements_old = self.elements + nnodes * i
            elements_new = self.elements + nnodes * (i + 1)

            pents[ih1:ih2, 0:3] = deepcopy(elements_old)
            pents[ih1:ih2, 3:6] = deepcopy(elements_new)
            nodes[in1:in2, :] = deepcopy(nodes_new)

            in1 += nnodes
            in2 += nnodes
            ih1 += nelements
            ih2 += nelements

        with open(bdf_filename, 'wb') as f:
            f.write('CEND\n')
            f.write('BEGIN BULK\n')

            pents += 1
            cid = None
            for nid, grid in enumerate(nodes):
                if nid % 5000 == 0:
                    print('writing nid=%s' % (nid + 1))

                card = ['GRID', nid + 1, cid, ] + list(grid)
                f.write(print_card_16(card))

            pid = 0
            mid = 1
            for eid, penta in enumerate(pents):
                if (eid + 1) % nelements == 1:
                    pid += 1
                    card = ['PSOLID', pid, mid]
                    f.write(print_card_8(card))
                    print('bumping pid -> %s' % pid)
                if eid % 5000 == 0:
                    print('writing eid=%s' % (eid + 1))
                card = ['CPENTA', eid + 1, pid, ] + list(penta)
                f.write(print_card_8(card))

            card = ['MAT1', mid, 1.0e7, None, 0.3]
            f.write(print_card_8(card))
            f.write('ENDDATA\n')

def main():
    cart3d_filename = 'threePlugs_bin.tri'
    cart3d = Cart3d_Mesher(log=None, debug=False)
    cart3d.read_cart3d(cart3d_filename)

    x0 = 10.
    nlayers = 5
    bdf_filename = 'threePlugs_volume.bdf'
    cart3d.project(bdf_filename, x0, growth_rate=1.3, nlayers=nlayers)

    #cart3d.write_bdf(bdf_filename)



if __name__ == '__main__':
    main()

