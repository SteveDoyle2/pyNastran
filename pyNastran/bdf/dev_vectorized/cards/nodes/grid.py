from numpy import zeros, arange, where, searchsorted, argsort, unique, asarray, array, dot, transpose, append

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
        cards = self._cards
        ncards = len(cards)

        self.n = ncards
        if ncards:
            self.model.log.debug('--------building grid--------')
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

    def get_index_by_node_id(self, node_id=None):
        if node_ids is None:
            out_index = None
        else:
            out_index = searchsorted(self.node_id, node_ids)
            assert len(node_ids) == len(n), 'n1=%s n2=%s'  %(len(node_ids), len(n))
        return out_index

    def get_index_by_cp(self, cp=None, i=None):
        """Find all the j-indicies where cp=cpi for some given subset of i-indicies"""
        return self._get_index_by_param('cp', self.cp, cp, i)

    def get_index_by_cd(self, cd=None, i=None):
        """Find all the j-indicies where cd=cdi for some given subset of i-indicies"""
        return self._get_index_by_param('cd', self.cd, cd, i)

    def get_index_by_seid(self, seid=None, i=None):
        """Find all the j-indicies where seid=seidi for some given subset of i-indicies"""
        return self._get_index_by_param('seid', self.seid, seid, i)

    def _get_index_by_param(self, name, param_data, param, i):
        """
        You probably shouldn't be calling this method.
        It does the work associcated with get_index_by_cp / get_index_by_cd
        """
        if param is None:
            return i
        #param_all = unique(param_data)
        i, n = _index_to_nslice(i, self.n)
        out_index = array(n, dtype='int32')
        param_data_i = param_data[i]

        i0 = 0
        for parami in param:
            j = where(param_data_i == parami)[0]
            nj = len(j)
            out_index[i0:i0+nj] = j
            i0 += nj
        return out_index

    def get_index_by_cp2(self, cp=None, i=None):
        """Find all the j-indicies where cp=cpi for some given subset of i-indicies"""
        if cp is None:
            return i
        #cp_all = unique(self.cp)
        i, n = index_to_nslice(i, self.n)
        out_index = array(n, dtype='int32')
        Cp = self.cp[i]

        i0 = 0
        for cpi in cp:
            j = where(Cp == cpi)[0]
            nj = len(j)
            out_index[i0:i0+nj] = j
            i0 += nj
        return out_index

    def get_positions(self, node_id=None):
        i = self.get_index_by_node_id(node_id)
        return self.get_positions_by_index(i)

    def get_positions_by_index(self, i=None):
        """
        in the global frame
        """
        if i is None:
            xyz = self.xyz.copy()
            n = slice(None, None)
        else:
            i = n
            xyz = self.xyz[n, :].copy()

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

    def write_bdf(self, f, node_id=None, size=8, is_double=False):
        i = self.get_index_by_node_id(node_id)
        return self.write_bdf_by_index(f, i, size, is_double)

    def write_bdf_by_index(self, f, i=None, size=8, is_double=False):
        """
        Write the BDF cards

        :param f: a file object
        :param i: the indicies (default=None -> all)
        :param size: the field width (8/16)
        :param is_double: is this double precision (default=False)
        """
        if i is None:
            i = slice(None, None)
        if self.n:
            f.write('$GRID\n')
            # default to the GRDSET defaults
            #cp0 = self.model.grdset.cp
            #cd0 = self.model.grdset.cd
            #ps0 = self.model.grdset.ps
            #seid0 = self.model.grdset.seid

            # default to the GRID defaults
            cp0 = 0
            cd0 = 0
            ps0 = -1
            seid0 = 0
            blank = ' '*8 if size==8 else ' ' * 16
            Cp   = [cpi   if cpi   != cp0   else blank for cpi   in self.cp[i]]
            Cd   = [cdi   if cdi   != cd0   else blank for cdi   in self.cd[i]]
            Ps   = [psi   if psi   != ps0   else blank for psi   in self.ps[i]]
            Seid = [seidi if seidi != seid0 else blank for seidi in self.seid[i]]
            for (nid, cp, xyz, cd, ps, seid) in zip(self.node_id, Cp, self.xyz[i, :], Cd, Ps, Seid):
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

def _index_to_nslice(i, n):
    if i is None:
        i = slice(None, None)
        n = self.n
    else:
        n = len(i)
    return i, n

