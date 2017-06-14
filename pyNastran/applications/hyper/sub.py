from __future__ import print_function
from math import pi, atan2, sqrt
from six import  iteritems
from numpy import array, zeros

from pyNastran.dev.bdf_vectorized.bdf import BDF, to_fields, BDFCard
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.applications.hyper.cards import FLOW, SUBSONIC

class Subsonic(BDF):
    def __init__(self, log=None, debug=False):
        BDF.__init__(self, debug=debug, log=log)
        self.subsonic = {}
        self.flow = {}
        self.cards_to_read.add('SUBSONIC')
        self.cards_to_read.add('FLOW')

    def add_card(self, card_lines, card_name, comment='', is_list=True, has_none=True):
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
            subsonic = SUBSONIC(card_obj, comment)
            self.subsonic[subsonic.pid] = subsonic
            return
        elif card_name == 'FLOW':
            flow = FLOW(card_obj, comment)
            self.flow[flow.flow_id] = flow
            return
        BDF.add_card(self, card, card_name, comment=comment, is_list=True)

    def _write_common(self, bdf_file, size=8, is_double=False):
        msg = ''
        for fid, flow in sorted(iteritems(self.flow)):
            msg += str(flow)
        for fid, subsonic in sorted(iteritems(self.subsonic)):
            msg += str(subsonic)
        bdf_file.write(msg)
        BDF._write_common(self, bdf_file, size=size, is_double=is_double)

    def get_pressure(self, isubcase=1):
        """
        doesn't support solid elements
        """
        #print(self.flow)
        #if 0:
            #flow = self.flow[1]
            #self.omega = flow.get_omega()
            ##flow.get_xyzref()
            #xyz_ref = array([0., 0., 0.])

            #a = 1.0
            #V, Vn = flow.get_V()
            ##mach = V / a
            ##mn = Vn / a
            #q = 1.
            #pinf = 0.

        xyz = self.grid.get_position_by_node_index()

        elements = self.elements.elements_shell
        nelements = elements.ctria3.n + elements.cquad4.n

        area = zeros(nelements, dtype='float32')
        centroid = zeros(nelements, dtype='float32')
        R = zeros((nelements, nelements), dtype='float32')
        print(self.elements.elements_shell.cquad4.__dict__.keys())
        nctria3 = elements.ctria3.n
        ncquad4 = elements.cquad4.n
        n0 = 0
        if nctria3 > 0:
            element_id = elements.ctria3.element_id
            area[n0:n0+nctri3] = elements.get_area_by_element_index(element_id, xyz_cid0=xyz)
            #area[n0:n0+nctri3] = elements.ctria3.get_area_by_element_index(xyz_cid0=xyz)
            #centroid[n0:n0+cntri3] = elements.ctria3.get_centroid_by_element_index(xyz_cid0=xyz)
        n0 += nctri3

        if ncquad4 > 0:
            area[n0:n0+ncquad4] = elements.get_area_by_element_index(element_id, xyz_cid0=xyz)
            #area[n0:n0+ncquad4] = elements.cquad4.get_area_by_element_index(xyz_cid0=xyz)
            centroid[n0:n0+ncquad4] = elements.cquad4.get_centroid_by_element_index(xyz_cid0=xyz)
        n0 += ncquad4

        for ne in range(nelements):
            Ri = centroid[ne:] - centroid[ne]
            R[ne, :] = Ri
            R[:, ne] = -Ri

        print(R)
        if 1:
            # quad - constant strength source
            k = 1/(4*pi)

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

            m12 = (y2-y1) / (x2-x1)
            m23 = (y3-y2) / (x3-x2)
            m34 = (y4-y3) / (x4-x3)
            m41 = (y1-y4) / (x1-x4)


            r = sqrt((x-x1)**2 + (y-y1)**2 + z**2)
            e1 = (x-x1)**2 + z**2
            e2 = (x-x2)**2 + z**2
            e3 = (x-x3)**2 + z**2
            e4 = (x-x4)**2 + z**2

            h1 = (x-x1) * (y-y1)
            h2 = (x-x2) * (y-y2)
            h3 = (x-x3) * (y-y3)
            h4 = (x-x4) * (y-y4)

            phi = k * ((
                (x-x1)*(y2-y1)-(y-y1)*(x2-x1)/d12 * ln((r1+r2+d12)/(r1+r2-d12)) +
                (x-x2)*(y3-y2)-(y-y2)*(x3-x2)/d23 * ln((r2+r3+d23)/(r2+r3-d23)) +
                (x-x3)*(y4-y3)-(y-y3)*(x4-x3)/d34 * ln((r3+r4+d34)/(r3+r4-d34)) +
                (x-x4)*(y1-y4)-(y-y4)*(x1-x4)/d41 * ln((r4+r1+d41)/(r3+r4-d34))
            ) + abs(z) * (
                atan2(m12*e1-h1, z*r1) - atan2(m12*e2-h2, z*r2) +
                atan2(m23*e2-h2, z*r2) - atan2(m23*e3-h3, z*r3) +
                atan2(m34*e3-h3, z*r3) - atan2(m34*e4-h4, z*r4) +
                atan2(m41*e4-h4, z*r4) - atan2(m41*e1-h1, z*r1)
            ))

        Cp = 1 - 1 / Vref**2 * (Q ** 2 - 2 * dphi_dt)
        #if 0:
            #print('  eid=%i delta=%g flow=%s cp=%s' % (eid, degrees(delta_radians), flow_type, cp))
            ##assert cp > 0
            #assert q > 0, q
            ##assert pinf > 0, pinf
            #p = cp * q + pinf
            #card = ['PLOAD4', isubcase, eid, p]
            #self.add_card(card, 'PLOAD4', is_list=True)
            #del cp


def main():
    model = Subsonic(log=None, debug=True)
    model.read_bdf('subsonic.bdf')

    isubcase = 1
    vel = array([1., 0., 0.], dtype='float32')
    model.get_pressure(isubcase)
    model.write_bdf('subsonic_pressure.bdf')


if __name__ == '__main__':  # pragma: no cover
    main()

