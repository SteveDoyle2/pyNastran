from numpy import zeros, where, unique, transpose, dot, array, arange

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)


class TEMPD:
    type = 'TEMPD'
    def __init__(self, model):
        """
        Defines the TEMPD object.

        Parameters
        ----------
        model : BDF
           the BDF object

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

    def add_card(self, card: BDFCard, comment: str=''):
        self.comment = comment
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

    def write_card(self, bdf_file, size=8):
        if self.n:
            #n = 0
            for lid, t in zip(self.load_id, self.temperature_default):
                card = ['TEMPD', lid, t]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))


class TEMP:
    type = 'TEMP'
    def __init__(self, model):
        """
        Defines the TEMP object.

        Parameters
        ----------
        model : BDF
           the BDF object

        +------+-----+----+----+----+----+----+----+------+
        |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
        +======+=====+====+====+====+====+====+====+======+
        | TEMP | SID | G1 | T1 | G2 | T2 | G3 | T3 |      |
        +------+-----+----+----+----+----+----+----+------+
        """
        self.model = model
        self._comments = []

    def add_card(self, card: BDFCard, comment: str=''):
        self.load_id = integer(card, 1, 'load_id')
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            float_fmt = self.model.float_fmt
            self.nid = zeros(ncards, 'int32')
            self.temp = zeros(ncards, float_fmt)

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
        xyz = self.xyz.copy()

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

    def write_card(self, bdf_file, size=8):
        if self.n:
            #t0 = self.model.tempd.cp
            seid0 = self.model.grdset.seid
            ps0 = self.model.grdset.ps
            cd0 = self.model.grdset.cd
            cp0 = self.model.grdset.cp
            Cp = [cpi if cpi != cp0 else '' for cpi in self.cp]
            Cd = [cdi if cdi != cd0 else '' for cdi in self.cd]
            Ps = [psi if psi != ps0 else '' for psi in self.ps]
            Seid = [seidi if seidi != seid0 else '' for seidi in self.seid]
            for (nid, cp, xyz, cd, ps, seid) in zip(self.nid, Cp, self.xyz, Cd, Ps, Seid):
                card = ['TEMP', nid, cp, xyz[0], xyz[1], xyz[2], cd, ps, seid]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def __repr__(self):
        msg = "<GRID>\n"
        msg += '  nGRID = %i' % self.n
