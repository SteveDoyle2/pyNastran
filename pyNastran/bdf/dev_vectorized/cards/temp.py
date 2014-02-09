from numpy import zeros

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)


class TEMPD(object):
    type = 'TEMPD'
    def __init__(self, model):
        """
        Defines the TEMPD object.

        :param self: the TEMPD object
        :param model: the BDF object

        +-------+------+----+------+----+------+----+------+----+
        |    1  |  2   | 3  |  4   | 5  |  6   | 7  |  8   | 9  |
        +=======+======+====+======+====+======+====+======+====+
        | TEMPD | SID1 | T1 | SID2 | T2 | SID3 | T3 | SID4 | T4 |
        +-------+------+----+------+----+------+----+------+----+
        """
        self.model = model
        self._comment = ['']
        #: card count
        self.n = 0
        #: load case ID
        self.load_id = []
        #: Default temperature
        self.temperature_default = []

    def add(self, card, comment):
        self._comment = comment
        self.n += 1

        n = 1
        for i in range(1, len(card)):
            lid = integer(card, i, 'load_id%i' % n)
            t = double(card, i, 'tempd%i' % n)
            n += 1
            self.load_id.append(lid)
            self.temperature_default.append(t)

    def build(self):
        self.load_id = array(self.load_id)
        self.temperature_default = array(self.temperature_default)
        
    def write_bdf(self, f, size=8):
        if self.n:
            n = 0
            for lid, t in izip(self.load_id, self.temperature_default):
                card = ['TEMPD', lid, t]
                f.write(print_card(card, size))


class TEMP(object):
    type = 'TEMP'
    def __init__(self, model):
        """
        Defines the TEMP object.

        :param self: the TEMP object
        :param model: the BDF object

        +------+-----+----+----+----+----+----+----+------+
        |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
        +======+=====+====+====+====+====+====+====+======+
        | TEMP | SID | G1 | T1 | G2 | T2 | G3 | T3 |      |
        +------+-----+----+----+----+----+----+----+------+
        """
        self.model = model
        self._comments = []

    def add(self, card, comment):
        self.load_id = integer(card, 1, 'load_id')
        
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

            cp0 = self.model.grdset.cp
            cd0 = self.model.grdset.cd
            ps0 = self.model.grdset.ps
            seid0 = self.model.grdset.seid
            for i, card in enumerate(cards):
                #: Node ID
                self.nid[i] = integer(card, 1, 'nid')

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

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('TEMP', self.n))
        return msg

    def write_bdf(self, f, size=8):
        if self.n:
            t0 = self.model.tempd.cp
            seid0 = self.model.grdset.seid

            Cp   = [cpi   if cpi   != cp0   else '' for cpi   in self.cp]
            Cd   = [cdi   if cdi   != cd0   else '' for cdi   in self.cd]
            Ps   = [psi   if psi   != ps0   else '' for psi   in self.ps]
            Seid = [seidi if seidi != seid0 else '' for seidi in self.seid]
            for (nid, cp, xyz, cd, ps, seid) in zip(self.nid, Cp, self.xyz, Cd, Ps, Seid):
                card = ['GRID', nid, cp, xyz[0], xyz[1], xyz[2], cd, ps, seid]
                f.write(print_card(card, size))

    def __repr__(self):
        msg = "<GRID>\n"
        msg += '  nGRID = %i' % self.n