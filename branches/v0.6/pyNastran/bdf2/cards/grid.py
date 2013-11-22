from numpy import zeros

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)


class NODES(object):
    def __init__(self, model):
        self.model = model

    def build():
        self.spoint.build()
        self.grid.build()

    def write_bdf(self, f, size=8):
        f.write('$NODES\n')
        self.spoint.write_bdf(f)
        self.grid.write_bdf(f)

    def get_stats(self):
        msg = []
        types = [self.spoint, self.grid]
        for node in types:
            nnode = len(node.nid)
            if nnode:
                msg.append('  %-8s: %i' % (node.type, nnode))
        return msg


class SPOINT(object):
    type = 'SPOINT'
    def __init__(self, model):
        self.model = model
        self._spoint = []
        self._spoint_comments = []

    def add_spoint(self, card, comment):
        self._spoint.append(card)
        self._spoint_comment.append(comment)

    def build(self):
        self._spoint = []
        self._spoint_comment = []
        self.spoint = zeros(ncards, 'int32')

    def write_bdf(self, f, size=8):
        #..todo:: collapse the IDs
        card = ['SPOINT'] + list(self.spoint)
        f.write(print_card(card))

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
        self._grid = []
        self._grid_comment = []

    def add(self, card, comment):
        self._grid.append(card)
        self._grid_comment.append(comment)

    def build(self):
        cards = self._grid
        ncards = len(cards)

        self.nid = zeros(ncards, 'int32')
        self.xyz = zeros((ncards, 3), 'float64')
        self.cp = zeros(ncards, 'int32')
        self.cd = zeros(ncards, 'int32')
        self.seid = zeros(ncards, 'int32')
        #self.ps = zeros(ncards, 'int32')

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
            #self.ps[i] = str(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            self.seid[i] = integer_or_blank(card, 8, 'seid', 0)

    def get_stats(self):
        msg = []
        ngrid = len(self.nid)
        if ngrid:
            msg.append('  %-8s: %i' % ('GRID', ngrid))
        return msg

    def write_bdf(self, f, size=8):
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