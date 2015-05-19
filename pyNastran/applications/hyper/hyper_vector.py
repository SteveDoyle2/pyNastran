from six import  iteritems
#from math import pi, degrees
from pyNastran.bdf.dev_vectorized.bdf import BDF, to_fields, BDFCard
from pyNastran.bdf.cards.utils import wipe_empty_fields
#from pyNastran.applications.hyper
from cards import FLOW, HYPER
from numpy import array, cross, vstack, dot, arccos, pi, degrees, arange, hstack #, radians
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
        xyz_global = grids.get_position_by_node_index()

        tris = self.elements.elements_shell.ctria3
        quads = self.elements.elements_shell.cquad4

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
            e = vstack([et, eq])
            pids = vstack([pt, pq])
            A = vstack([At, Aq])
            n = vstack([nt, nq])
            c = vstack([ct, cq])
            o = vstack([ot, oq])
        elif tris.n:
            e = et
            pids = pt
            A = At
            n = nt
            c = ct
            #o = ot
        elif quads.n:
            e = eq
            pids = pq
            A = Aq
            n = nq
            c = cq
            #o = oq

        positions = {}
        for nid, node in iteritems(self.nodes):
            positions[nid] = node.Position()

        upids = unique(pid)
        for upidi in upids:
            i = where(pids == upidi)[0]
            hyper = self.hyper[pid]
            gamma = hyper.gamma
            #_cp = hyper._cps

            #for iset, set_id in enumerate(hyper.sets):
                #SET = self.sets[set_id]
                #if eid not in SET.IDs:
                    #continue
                #hyper.Type = hyper.Types[iset]
                #hyper._cp = hyper._cps[iset]
                #hyper.set = hyper.sets[iset]
                #hyper.gamma = hyper.gammas[iset]


            for element in [ctria3, cquad4]:
                if element.type == 'CTRIA3':
                    i1, i2, i3 = ctria3.get_node_indices()
                    p1, p2, p3 = positions[i1], positions[i2], positions[i3]

                    # hacked an extra point in to make the equations work nice
                    xyz = vstack([p1, p2, p3, p3])
                    xyz_avg = (p1 + p2 + p3) / 3.
                    a = p2 - p1
                    b = p3 - p1
                if element.type == 'CQUAD4':
                    i1, i2, i3, i4 = cquad4.get_node_indices()
                    xyz = positions[i1, :], positions[i2, :], positions[i3, :], positions[i4, :]
                    p1, p2, p3, p4 = xyz[:, 0], xyz[:, 1], xyz[:, 2], xyz[:, 3]
                    #xyz = vstack([p1, p2, p3, p4])
                    xyz_avg = (p1 + p2 + p3 + p4) / 4.
                    a = p3 - p1
                    b = p4 - p2

                n = cross(a, b)
                norm_n = norm(n)
                #A = 0.5 * norm_n
                n /= norm_n

                t1 = a / norm(a)
                t2 = cross(n, t1)

                # corner point projection distance
                dk = (n[0] * (xyz_avg[0] - xyz[:, 0]) +
                      n[1] * (xyz_avg[1] - xyz[:, 1]) +
                      n[2] * (xyz_avg[2] - xyz[:, 2]))

                x_prime = xyz[:, 0] + n[0] * dk[0]
                y_prime = xyz[:, 1] + n[1] * dk[1]
                z_prime = xyz[:, 2] + n[2] * dk[2]

                # zeta_k_star
                #t1*dk - xyz_avg
                #zeta = t1 @ ((xyz + n @ dk) - xyz_avg))
                #zeta = (t1[:, 0] * (x_prime - xyz_avg[:, 0]) +
                #        t1[:, 1] * (y_prime - xyz_avg[:, 1]) +
                #        t1[:, 2] * (z_prime - xyz_avg[:, 2]) )
                #eta = t2 @ ((xyz + n @ dk) - xyz_avg))
                eta  = (t2[:, 0] * (x_prime - xyz_avg[:, 0]) +
                        t2[:, 1] * (y_prime - xyz_avg[:, 1]) +
                        t2[:, 2] * (z_prime - xyz_avg[:, 2]) )

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

