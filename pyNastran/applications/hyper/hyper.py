from math import pi
from pyNastran.bdf.bdf import BDF, to_fields, wipe_empty_fields, BDFCard
#from pyNastran.applications.hyper
from cards import FLOW, HYPER
from numpy import array, cross
from numpy.linalg import norm


class Hypersonic(BDF):
    def __init__(self, log=None, debug=False):
        BDF.__init__(self, debug=debug, log=log)
        self.hyper = {}
        self.flow = {}
        self.cardsToRead.add('HYPER')
        self.cardsToRead.add('FLOW')

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
        for fid,flow in sorted(self.flow.iteritems()):
            msg += str(flow)
        for fid,hyper in sorted(self.hyper.iteritems()):
            msg += str(hyper)
        msg += BDF._write_common(self, size, card_writer)
        return msg

    def get_pressure(self, isubcase=1):
        """
        doesn't support solid elements
        """
        #print self.flow
        flow = self.flow[1]
        a = 1.0
        V, Vn = flow.get_V()
        mach = V / a
        mn = Vn / a
        q = 1.
        pinf = 1.

        positions = {}
        for nid, node in self.nodes.iteritems():
            positions[nid] = node.Position()
        
        for eid, element in self.elements.iteritems():  
            pid = element.Pid()
            hyper = self.hyper[pid]

            for iset, set_id in enumerate(hyper.sets):
                SET = self.sets[set_id]
                if eid not in SET.IDs:
                    continue
                hyper.Type = hyper.Types[iset]
                hyper._cp = hyper._cps[iset]
                hyper.set = hyper.sets[iset]
                hyper.gamma = hyper.gammas[iset]
            if element.type in ['CTRIA3', 'CTRIA6']:
                n1, n2, n3 = element.nodeIDs()
                p1, p2, p3 = positions[n1], positions[n2], positions[n3]
                a = p2 - p1
                b = p3 - p1
                n = cross(a, b)
                n2 = norm(n)
                #sin2_theta = dot(Vinf, n/n2)**2 / uinf**2
                #if sin2_theta < 0:
                    #Cp = 0.
                #else:
                    #Cp = 2*sin2_theta
            elif element.type in ['CQUAD4', 'CQUAD8']:
                n1, n2, n3, n4 = element.nodeIDs()
                p1, p2, p3, p4 = positions[n1], positions[n2], positions[n3], positions[n4]
                xyz_avg = (p1 + p2 + p3 + p4) / 4.
                a = p2 - p1
                b = p3 - p1
                n = cross(a, b)
                n2 = norm(n)
                # a x b = |a|*|b| * sin(theta)
                # a o b = |a|*|b| * cos(theta)
                #sin_delta = n / (Vn * n2)
                #cos_delta = n / (Vn * n2)
            else:
                continue

            pi2 = pi / 2.
            theta = norm(n / (Vn * n2))
            delta_radians = pi2 - theta
            cp = hyper.Cp(delta_radians)
            print cp
            assert cp > 0
            assert q > 0
            assert pinf > 0
            p = cp * q + pinf
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


if __name__ == '__main__':
    main()

