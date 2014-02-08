from numpy import zeros, arange, where, searchsorted

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)


class Nodes(object):
    def __init__(self, model):
        self.model = model
        self.spoint = SPOINT(model)
        self.grdset = GRDSET(model)
        self.grid = GRID(model)
        self.point = POINT(model)

    def build(self):
        self.spoint.build()
        self.grid.build()
        self.point.build()

    def write_bdf(self, f, size=8, nids=None):
        f.write('$NODES\n')
        self.spoint.write_bdf(f, size, nids)
        self.grdset.write_bdf(f, size, nids)
        self.grid.write_bdf(f, size, nids)
        self.point.write_bdf(f, size, nids)

    def ndofs(self, sol):
        if self.model.sol in [101, 103, 144, 145]:
            ndofs = (6 * self.grid.n) + self.spoint.n
        elif self.model.sol in [159]:
            ndofs = self.grid.n + self.spoint.n
        else:
            raise NotImplementedError('sol=%r' % sol)
        return ndofs

    def get_stats(self):
        msg = []
        types = [self.spoint, self.grdset, self.grid, self.point]
        for node in types:
            if node.n:
                msg.append('  %-8s: %i' % (node.type, node.n))
        return msg


class GRDSET(object):
    type = 'GRDSET'
    def __init__(self, model):
        """
        Defines the GRID object.

        :param self: the GRID object
        :param model: the BDF object

        +--------+-----+----+----+----+----+----+----+------+
        |    1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
        +========+=====+====+====+====+====+====+====+======+
        | GRDSET |     | CP |    |    |    | CD | PS | SEID |
        +--------+-----+----+----+----+----+----+----+------+
        """
        self.model = model
        self._comment = ['']
        #: card count
        self.n = 0
        #: Grid point coordinate system
        self.cp = 0
        #: Analysis coordinate system
        self.cd = 0
        #: SPC constraint
        self.ps = -1
        #: Superelement ID
        self.seid = 0

    def add(self, card, comment):
        self._comment = comment
        self.n = 1
        self.cp = integer_or_blank(card, 2, 'cp', 0)
        self.cd = integer_or_blank(card, 6, 'cd', 0)
        self.ps = integer_or_blank(card, 7, 'ps', -1)
        self.seid = integer_or_blank(card, 8, 'seid', 0)

    def write_bdf(self, f, size=8):
        if self.n:
            card = ['GRDSET', None, self.cp, None, None, None, self.cd, self.seid]
            f.write(print_card(card, size))


class GRID(object):
    type = 'GRID'
    def __init__(self, model):
        """
        Defines the GRID object.

        :param self: the GRID object
        :param model: the BDF object

        +------+-----+----+----+----+----+----+----+------+
        |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
        +======+=====+====+====+====+====+====+====+======+
        | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
        +------+-----+----+----+----+----+----+----+------+
        """
        self.model = model
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            self.node_id = zeros(ncards, 'int32')
            self.xyz = zeros((ncards, 3), float_fmt)
            self.cp = zeros(ncards, 'int32')
            self.cd = zeros(ncards, 'int32')
            self.seid = zeros(ncards, 'int32')
            self.ps = zeros(ncards, 'int32')

            cp0 = self.model.grdset.cp
            cd0 = self.model.grdset.cd
            ps0 = self.model.grdset.ps
            seid0 = self.model.grdset.seid
            for i, card in enumerate(cards):
                #: Node ID
                self.node_id[i] = integer(card, 1, 'nid')

                #: Grid point coordinate system
                self.cp[i] = integer_or_blank(card, 2, 'cp', cp0)

                x = double_or_blank(card, 3, 'x1', 0.)
                y = double_or_blank(card, 4, 'x2', 0.)
                z = double_or_blank(card, 5, 'x3', 0.)
                #: node location in local frame
                self.xyz[i] = [x, y, z]

                #: Analysis coordinate system
                self.cd[i] = integer_or_blank(card, 6, 'cd', cd0)

                #: SPC constraint
                self.ps[i] = integer_or_blank(card, 7, 'ps', ps0)

                #: Superelement ID
                self.seid[i] = integer_or_blank(card, 8, 'seid', seid0)

    def position(self, node_ids=None):
        if node_ids is None:
            node_ids = self.node_id
            xyz = self.xyz.copy()
            # indexs
            n = arange(self.n)
        else:
            # indexs
            n = searchsorted(self.node_id, node_ids)
            assert len(node_ids) == len(n), 'n1=%s n2=%s'  %(len(node_ids), len(n))
            xyz = self.xyz[n, :].copy()
            print "n =", n

        cpn = self.cp[n]
        i = where(cpn != 0)[0]
        if i:
            n2 = n[i]
            cps = set(list(unique(cpn)))
            for cp in cps:
                i = where(cpn != 0)[0]
                T = self.model.coord.transform(cp)
                xyzi = xyz[n2[i], :]
                xyzi = dot(transpose(T), dot(xyzi, T))
        
        assert len(node_ids) == len(cpn), 'n1=%s n2=%s'  %(len(node_ids), len(cpn))
        return xyz

    def positions_wrt(self, node_ids=None, coord_ids=None):
        raise NotImplementedError()

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('GRID', self.n))
        return msg

    def write_bdf(self, f, size=8):
        if self.n:
            f.write('$GRID\n')
            cp0 = self.model.grdset.cp
            cd0 = self.model.grdset.cd
            ps0 = self.model.grdset.ps
            seid0 = self.model.grdset.seid

            Cp   = [cpi   if cpi   != cp0   else '' for cpi   in self.cp]
            Cd   = [cdi   if cdi   != cd0   else '' for cdi   in self.cd]
            Ps   = [psi   if psi   != ps0   else '' for psi   in self.ps]
            Seid = [seidi if seidi != seid0 else '' for seidi in self.seid]
            for (nid, cp, xyz, cd, ps, seid) in zip(self.node_id, Cp, self.xyz, Cd, Ps, Seid):
                card = ['GRID', nid, cp, xyz[0], xyz[1], xyz[2], cd, ps, seid]
                f.write(print_card(card, size))

    def __repr__(self):
        msg = "<GRID>\n"
        msg += '  nGRID = %i' % self.n