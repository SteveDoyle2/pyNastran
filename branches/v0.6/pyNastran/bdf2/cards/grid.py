from numpy import zeros

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)


class Nodes(object):
    def __init__(self, model):
        self.model = model
        self.spoint = SPOINT(model)
        self.grid = GRID(model)

    def build(self):
        self.spoint.build()
        self.grid.build()

    def write_bdf(self, f, size=8):
        f.write('$NODES\n')
        self.spoint.write_bdf(f)
        self.grid.write_bdf(f)

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
        types = [self.spoint, self.grid]
        for node in types:
            if node.n:
                msg.append('  %-8s: %i' % (node.type, node.n))
        return msg


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
            self.nid = zeros(ncards, 'int32')
            self.xyz = zeros((ncards, 3), float_fmt)
            self.cp = zeros(ncards, 'int32')
            self.cd = zeros(ncards, 'int32')
            self.seid = zeros(ncards, 'int32')
            self.ps = zeros(ncards, 'int32')

            for i, card in enumerate(cards):
                #: Node ID
                self.nid[i] = integer(card, 1, 'nid')

                #: Grid point coordinate system
                self.cp[i] = integer_or_blank(card, 2, 'cp', 0)

                x = double_or_blank(card, 3, 'x1', 0.)
                y = double_or_blank(card, 4, 'x2', 0.)
                z = double_or_blank(card, 5, 'x3', 0.)
                #: node location in local frame
                self.xyz[i] = [x, y, z]

                #: Analysis coordinate system
                self.cd[i] = integer_or_blank(card, 6, 'cd', 0)

                #: SPC constraint
                self.ps[i] = integer_or_blank(card, 7, 'ps', -1)

                #: Superelement ID
                self.seid[i] = integer_or_blank(card, 8, 'seid', 0)

    def positions(self, nids=None):
        if nids is None:
            nids = self.nids
        xyz = xyz.copy()

        n = arange(self.n)
        i = where(self.cid != 0)[0]
        if i:
            n = n[i]
            cids = set(list(unique(self.cid)))
            for cid in cids:
                i = where(self.cid != 0)[0]
                T = self.model.coord.transform(cid)
                xyzi = xyz[n[i], :]
                xyzi = dot(transpose(T), dot(xyzi, T))
        return xyz

    def positions_wrt(self, nids=None, cids=None):
        raise NotImplementedError()

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('GRID', self.n))
        return msg

    def write_bdf(self, f, size=8):
        if self.n:
            f.write('$GRID\n')
            cp = [cpi     if cpi   != 0 else '' for cpi in self.cp]
            cd = [cdi     if cdi   != 0 else '' for cdi in self.cd]
            seid = [seidi if seidi != 0 else '' for seidi in self.seid]
            for (nid, cp, xyz, cd, seid) in zip(self.nid, cp, self.xyz,
                    cd, seid):

                ps = None
                card = ['GRID', nid, cp, xyz[0], xyz[1], xyz[2], cd, seid]
                f.write(print_card(card, size))

    def __repr__(self):
        msg = "<GRID>\n"
        msg += '  nGRID = %i' % len(self.cp)