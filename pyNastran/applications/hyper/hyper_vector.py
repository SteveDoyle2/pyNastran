from __future__ import print_function
from six import  iteritems
#from math import pi, degrees
from pyNastran.dev.bdf_vectorized.bdf import BDF, to_fields, BDFCard
from pyNastran.bdf.cards.utils import wipe_empty_fields
#from pyNastran.applications.hyper
from cards import FLOW, HYPER
from numpy import (array, cross, vstack, dot, arccos, pi, degrees, arange,
    hstack, append, unique, where, array_equal, argsort, searchsorted,
    zeros) #, radians
from numpy.linalg import norm


class Hypersonic(BDF):
    def __init__(self, log=None, debug=False):
        BDF.__init__(self, debug=debug, log=log)
        self.hyper = {}
        self.flow = {}
        self.cards_to_read.add('HYPER')
        self.cards_to_read.add('FLOW')

    def add_card(self, card_lines, card_name, comment='', is_list=True):
        card_name = card_name.upper()
        self._increase_card_count(card_name)
        if card_name in ['DEQATN']:
            card_obj = card_lines
            card = card_lines
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [self._parse_dynamic_syntax(field) if '%' in
                          field[0:1] else field for field in fields]

                card = wipe_empty_fields([interpret_value(field, fields)
                                          if field is not None
                                          else None for field in fields])
            else:  # leave everything as strings
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card)

        if card_name == 'HYPER':
            hyper = HYPER(card_obj, comment)
            self.hyper[hyper.pid] = hyper
            return
        elif card_name == 'FLOW':
            flow = FLOW(card_obj, comment)
            self.flow[flow.flow_id] = flow
            return
        BDF.add_card(self, card, card_name, comment=comment, is_list=True)

    def _write_common(self, size, card_writer):
        msg = ''
        for fid,flow in sorted(iteritems(self.flow)):
            msg += str(flow)
        for fid,hyper in sorted(iteritems(self.hyper)):
            msg += str(hyper)
        msg += BDF._write_common(self, size, card_writer)
        return msg

    def get_pressure(self, isubcase=1):
        """
        doesn't support solid elements
        """
        #print(self.flow)
        flow = self.flow[1]
        self.omega = flow.get_omega()
        #flow.get_xyzref()
        xyz_ref = array([0., 0., 0.])

        a = 1.0
        V, Vn = flow.get_V()
        mach = V / a
        mn = Vn / a
        q = 1.
        pinf = 0.

        grids = self.grid
        #print(list(grids.node_id))
        xyz_cid0 = grids.get_position_by_node_index()

        tris = self.elements.elements_shell.ctria3
        quads = self.elements.elements_shell.cquad4

        ntris = tris.n
        nquads = quads.n
        if tris.n:
            et = tris.element_id
            #pt = tris.property_id
            #At = tris.get_area_by_element_index()
            #nt = tris.get_normal_by_element_index()
            #ct = tris.get_centroid_by_element_index()

            i = arange(tris.n)
            n1, n2, n3 = tris._node_locations(xyz_cid0, i)
            #n3 = n4
            # is this right?
            #ut = (n1 + n2 - 2 * n3) / 2.
            #pt = (n2 - n1) / 2.
            #ot = cross(nt, ut)


            #i1, i2, i3 = tris.get_node_indices()
            n1, n2, n3 = tris._node_locations(xyz_cid0, i)

            # hacked an extra point in to make the equations work nice
            #xyz = vstack([n1, n2, n3, n3])
            xyz_avg_t = (n1 + n2 + n3) / 3.
            a = n2 - n1
            b = n3 - n1
            xyz_avgt = (n1 + n2 + n3) / 3.
            n = cross(a, b)
            norm_n = norm(n, axis=1)
            assert len(norm_n) == ntris, 'nnorm=%s nelem=%s' % (len(norm_n), ntris)
            #A = 0.5 * norm_n

            n /= norm_n.reshape(ntris, 1)

            t1 = a / norm(a, axis=1).reshape(ntris, 1)
            t2 = cross(n, t1)
            assert len(t1) == ntris, 'len(t1)=%s ntris=%s' %(len(t1), ntris)
            assert len(t2) == ntris, 'len(t2)=%s ntris=%s' %(len(t2), ntris)
            nt = n
            t1t = t1
            t2t = t2
            n1234 = zeros((ntris, 4), dtype='int32')
            n1234[:, :3] = tris.node_ids
            n1234[:, 3] = tris.node_ids[:, 0]
            n1234t = n1234

        if quads.n:
            eq = quads.element_id

            i = arange(quads.n)
            n1, n2, n3, n4 = quads._node_locations(xyz_cid0, i)
            #xyz = vstack([p1, p2, p3, p4])
            xyz_avg = (n1 + n2 + n3 + n4) / 4.
            a = n3 - n1
            b = n4 - n2

            #pq = quads.property_id
            #Aq = quads.get_area_by_element_index()
            #nq = quads.get_normal_by_element_index()
            #cq = quads.get_centroid_by_element_index()

            i = arange(quads.n)
            n1, n2, n3, n4 = quads._node_locations(xyz_cid0, i)
            uq = (n1 + n2 - n3 - n4) / 2.
            #pq = (n2 + n3 - n4 - n1) / 2.
            #oq = cross(nq, uq)
            xyz_avg_q = (n1 + n2 + n3 + n4) / 4.

            n = cross(a, b)
            norm_n = norm(n, axis=1)
            n /= norm_n.reshape(nquads, 1)
            t1 = a / norm(a, axis=1).reshape(nquads, 1)
            t2 = cross(n, t1)
            assert len(t1) == nquads, 'len(t1)=%s nquads=%s' %(len(t1), nquads)
            assert len(t2) == nquads, 'len(t2)=%s nquads=%s' %(len(t2), nquads)
            nq = n
            t1q = t1
            t2q = t2
            n1234 = quads.node_ids
            n1234q = n1234

        if tris.n and quads.n:
            eids = append(et, eq)
            #pids = append(pt, pq)
            #A = append(At, Aq)
            n = vstack([nt, nq])
            #c = append(ct, cq)
            #o = append(ot, oq)
            xyz_avg = vstack([xyz_avg_t, xyz_avg_q])
            t1 = vstack([t1t, t1q])
            t2 = vstack([t2t, t2q])
            n1234 = vstack([n1234t, n1234q])
        elif tris.n:
            eids = et
            #pids = pt
            #A = At
            n = nt
            #c = ct
            #o = ot
        elif quads.n:
            eids = eq
            #pids = pq
            #A = Aq
            n = nq
            #c = cq
            #o = oq
        nelem = tris.n + quads.n
        i = argsort(eids)
        eids = eids[i]
        #A = A[i]
        print(n)
        n = n[i, :]
        t1 = t1[i, :]
        t2 = t2[i, :]
        #o = o[i]
        #positions = {}
        #for nid, node in iteritems(self.nodes):
            #positions[nid] = node.get_position()

        #upids = unique(pids)
        print('looking for hyper cases=%s' % sorted(self.hyper.keys()))
        #n1234 =

        n1 = n1234[:, 0]
        n2 = n1234[:, 1]
        n3 = n1234[:, 2]
        n4 = n1234[:, 3]
        #n1234r = n1234.reshape(nquads * 4)
        print('nodes1234\n', n1234)
        n1234i = grids.get_node_index_by_node_id(n1234)
        print('nodes1234 index\n', n1234i)
        x = xyz_cid0[n1234i, 0]
        y = xyz_cid0[n1234i, 1]
        z = xyz_cid0[n1234i, 2]
        print('xyz\n', x)
        #y =
        #z =

        for key, hyper in sorted(iteritems(self.hyper)):
            for Type, set3_id in zip(hyper.Types, hyper.sets):
                print('Type=%s set3=%s' % (Type, set3_id))
                set3 = self.set3[set3_id]
                assert set3.desc == 'ELEM', set3.desc
                eids_set = set3.IDs
                print('eids = %s' % eids)
                print('eids_set = %s' % eids_set)
                j = searchsorted(eids, eids_set)
                assert array_equal(eids[j], eids_set)
                print('j = ', j)
                nj = len(j)


                # corner point projection distance
                #dka = n[j, 0] * (xyz_avg[j, 0] - x[j, 0])
                dk = zeros((nj, 4), dtype='float32')
                dk[:, 0] = n[j, 0] * (xyz_avg[j, 0] - x[j, 0])
                dk[:, 1] = n[j, 1] * (xyz_avg[j, 1] - x[j, 1])
                dk[:, 2] = n[j, 2] * (xyz_avg[j, 2] - x[j, 2])

                #dk = (n[j, 0] * (xyz_avg[j, 0] - x[:, 0]) +
                #      n[j, 1] * (xyz_avg[j, 1] - y[:, 1]) +
                #      n[j, 2] * (xyz_avg[j, 2] - z[:, 2]))
                assert dk.shape == (nelem, 4), 'dk.shape=%s nj=%s' % (str(dk.shape), nj)

                print('dk\n', dk)
                dk =  dk.reshape(nj, 1)
                x_prime = xyz_cid0[j, 0] + dot(n[j, 0], dk)
                y_prime = xyz_cid0[j, 1] + dot(n[j, 1], dk)
                z_prime = xyz_cid0[j, 2] + dot(n[j, 2], dk)

                print('nj.shape\n', n[j, ].shape)
                print('dk.shape\n', dk.shape)
                ndk = n[j, :] * dk
                xyz_prime = xyz_cid0[j, :] + ndk # n[j, :] * dk
                print('xyz_prime')
                print(xyz_prime)
                #x_prime = xyz_prime[:, 0]
                #y_prime = xyz_prime[:, 1]
                #z_prime = xyz_prime[:, 2]
                print('xyz_prime2')
                print(x_prime)
                print(y_prime)
                print(z_prime)

                # zeta_k_star
                #t1*dk - xyz_avg
                #zeta = t1 @ ((xyz + n @ dk) - xyz_avg))
                zeta = (t1[j, 0] * (x_prime - xyz_avg[j, 0]) +
                        t1[j, 1] * (y_prime - xyz_avg[j, 1]) +
                        t1[j, 2] * (z_prime - xyz_avg[j, 2]) )
                #eta = t2 @ ((xyz + n @ dk) - xyz_avg))
                eta  = (t2[j, 0] * (x_prime - xyz_avg[j, 0]) +
                        t2[j, 1] * (y_prime - xyz_avg[j, 1]) +
                        t2[j, 2] * (z_prime - xyz_avg[j, 2]) )

                if element.type in ['CTRIA3']:
                    # as eta4 becomes eta1
                    #     1
                    # -------- (e4*(n1-n2)+e2(n4-n1))
                    # 3(n2-n4)
                    inner = (zeta[3] * (eta[:, 0] - eta[:, 1]) +
                             zeta[1] * (eta[:, 3] - eta[:, 1]))
                    zeta_0 = 1./ (3. * (eta[:, 1] - eta[:, 3])) * inner
                else:
                    # p. 14
                    inner = (zeta[3] * (eta[:, 0] - eta[:, 1]) +
                             zeta[1] * (eta[:, 3] - eta[:, 1]))
                    zeta_0 = 1./ (3. * (eta[:, 1] - eta[:, 3])) * inner
                # p. 14
                eta_0 = -eta[:, 0] / 3.

                # p. 15
                zeta_k = zeta - zeta_0
                # p. 15
                eta_k = eta - eta_0

                # p. 15
                xyz_centroid = xyz_avg + t1 * zeta_0 + t2 * eta_0

                r = xyz_centroid - xyz_ref

                # p. 110
                Vbar = V - cross(self.omega, r)

                VdotN = dot(-V, n)
                if VdotN > 0.:
                    flow_type = 'leeward'
                else:
                    flow_type = 'windward'

                # p. 110
                Vlocal = norm(V)

                # p. 111
                delta_radians = pi / 2. - arccos(VdotN / Vlocal)
                cp = hyper.Cp(delta_radians)
                print('  eids=%i delta=%g flow=%s cp=%s' % (eids, degrees(delta_radians), flow_type, cp))
                #assert cp > 0
                assert q > 0, q
                #assert pinf > 0, pinf
                p = cp * q + pinf
                for pi, eidi in zip(p, eids):
                    card = ['PLOAD4', isubcase, eid, p]
                    self.add_card(card, 'PLOAD4', is_list=True)
                del cp


def main():
    h = Hypersonic(log=None, debug=True)
    h.read_bdf('hyper.bdf')

    isubcase = 1
    V = array([1., 0., 0.], dtype='float32')
    h.get_pressure(isubcase)
    h.write_bdf('hyper_pressure.bdf')


if __name__ == '__main__':  # pragma: no cover
    main()

