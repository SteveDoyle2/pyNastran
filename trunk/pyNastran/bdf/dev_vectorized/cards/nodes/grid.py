from numpy import zeros, arange, where, searchsorted, argsort, unique, asarray, array, dot, transpose

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter
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

    def allocate(self, card_count):
        ncards = card_count['GRID']
        float_fmt = self.model.float
        self.node_id = zeros(ncards, 'int32')
        self.xyz = zeros((ncards, 3), float_fmt)
        self.cp = zeros(ncards, 'int32')
        self.cd = zeros(ncards, 'int32')
        self.seid = zeros(ncards, 'int32')
        self.ps = zeros(ncards, 'int32')

    def build(self):
        print('--------building grid--------')
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
            i = argsort(self.node_id)
            self.cp = self.cp[i]
            self.xyz = self.xyz[i, :]
            self.cd = self.cd[i]
            self.ps = self.ps[i]
            self.seid = self.seid[i]

    def index_map(self, node_ids, msg=''):
        #return searchsorted(node_ids, self.node_id)
        #i_too_large = where(self.node_id[-1] < node_ids)[0]
        #if len(i_too_large):
            #raise RuntimeError('Cannot find GRID %s, %s' % (node_ids[i_too_large], msg))
        return searchsorted(self.node_id, node_ids)

    def get_positions(self, node_ids=None):
        """
        in the global frame
        """
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
            #print "n =", n

        cpn = self.cp[n]
        i = where(cpn != 0)[0]
        if len(i):
            n2 = n[i]
            cps = set(list(unique(cpn)))
            for cp in cps:
                #print self.model.coords
                T = self.model.coords.transform(cp)
                #print('T[%s] = \n%s\n' % (cp, T))
                j = where(self.cp[n] == cp)[0]
                #print('j = %s' % j)

                #if j.max() > len(n2):
                    #ii = where(i > len(n2))[0]
                    #i2 = i[ii]
                    ## save the bad data
                    #i = i2
                    #print('n2 = %s' % n2)
                    #print('i2 = %s' % i2)
                #j = n2[i]
                #print(j)
                xyzi = xyz[j, :]
                #xyzi = dot(transpose(T), dot(xyzi, T))
                xyz[j, :] = self.model.coords.get_global_position(xyzi, cp)


        assert len(node_ids) == len(cpn), 'n1=%s n2=%s'  %(len(node_ids), len(cpn))
        return xyz

    def get_positions_wrt(self, node_ids=None, coord_ids=None):
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

    def __getitem__(self, node_id):
        print('self.node_id = %s' % self.node_id)
        print('node_id = %s' % node_id)
        #node_id = slice_to_iter(node_id)
        i = where(self.node_id == node_id)[0]
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        #i = slice_to_iter(i)
        i = asarray(i)
        print('i = %s' % i, type(i))
        obj = GRID(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.node_id = self.node_id[i]
        obj.xyz = self.xyz[i, :]
        obj.cp = self.cp[i]
        obj.cd = self.cd[i]
        obj.ps = self.ps[i]
        obj.seid = self.seid[i]
        return obj
