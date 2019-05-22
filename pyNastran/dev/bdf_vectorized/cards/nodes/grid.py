from numpy import zeros, arange, where, argsort, unique, array

from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard



class Nodes:
    def __init__(self, model):
        self.model = model
        self.spoint = model.spoint
        self.epoint = model.epoint
        self.grdset = model.grdset
        self.grid = model.grid
        self.point = model.point

    def build(self):
        self.spoint.build()
        self.grid.build()
        self.point.build()

    def write_card(self, bdf_file, node_id=None, size=8, is_double=False, write_header=True):
        if write_header:
            bdf_file.write('$NODES\n')
        self.spoint.write_card(bdf_file, node_id, size, is_double, write_header)
        self.grdset.write_card(bdf_file, size, size, is_double, write_header)
        self.grid.write_card(bdf_file, size, size, is_double, write_header)
        self.point.write_card(bdf_file, size, size, is_double, write_header)

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


class GRDSET:
    type = 'GRDSET'
    def __init__(self, model):
        """
        Defines the GRID object.

        Parameters
        ----------
        model : BDF
           the BDF object

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

    def add_card(self, card, comment=''):
        self.comment = comment
        self.n = 1
        self.cp = integer_or_blank(card, 2, 'cp', 0)
        self.cd = integer_or_blank(card, 6, 'cd', 0)
        self.ps = integer_or_blank(card, 7, 'ps', -1)
        self.seid = integer_or_blank(card, 8, 'seid', 0)

    def write_card(self, bdf_file, size=8, is_double=False):
        if self.n:
            card = ['GRDSET', None, self.cp, None, None, None, self.cd, self.seid]
            if size == 8:
                bdf_file.write(print_card_8(card))
            else:
                bdf_file.write(print_card_16(card))


class GRID(VectorizedCard):
    type = 'GRID'
    def __init__(self, model):
        """
        Defines the GRID object.

        Parameters
        ----------
        model : BDF
           the BDF object

        +------+-----+----+----+----+----+----+----+------+
        |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
        +======+=====+====+====+====+====+====+====+======+
        | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
        +------+-----+----+----+----+----+----+----+------+
        """
        VectorizedCard.__init__(self, model)

    def shrink(self, refcheck=True):
        i = where(self.node_id == 0)[0]
        self.resize(i[0], refcheck=refcheck)

    def allocate(self, card_count):
        if 'GRID' in card_count:
            ncards = card_count['GRID']
            self.n = ncards
            #print('ngrid=%s' % self.n)
            float_fmt = self.model.float_fmt
            self.node_id = zeros(ncards, 'int32')
            self.xyz = zeros((ncards, 3), float_fmt)
            self.cp = zeros(ncards, 'int32')
            self.cd = zeros(ncards, 'int32')
            self.seid = zeros(ncards, 'int32')
            self.ps = zeros(ncards, 'int32')

    #def size_check(f, *args, **kwargs):
        #print('size_check', f.__name__)
        #return f(*args, **kwargs)

    #@size_check
    def add_card(self, card, comment=''):
        cp0 = self.model.grdset.cp
        cd0 = self.model.grdset.cd
        ps0 = self.model.grdset.ps
        seid0 = self.model.grdset.seid

        i = self.i
        nid = integer(card, 1, 'nid')
        if comment:
            self._comments[nid] = comment

        #: Node ID
        self.node_id[i] = nid

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
        self.i += 1

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'coord' : cid_map,
        }
        """
        nid_map = maps['node']
        cid_map = maps['coord']
        for i, (nid, cpi, cdi) in enumerate(zip(self.node_id, self.cp, self.cd)):
            self.node_id[i] = nid_map[nid]
            self.cp[i] = cid_map[cpi]
            self.cd[i] = cid_map[cdi]
        assert self.node_id.min() >= 0, self.node_id.min()

    def build(self):
        if self.n:
            i = argsort(self.node_id)
            self.node_id = self.node_id[i]
            self.cp = self.cp[i]
            self.xyz = self.xyz[i, :]
            self.cd = self.cd[i]
            self.ps = self.ps[i]
            self.seid = self.seid[i]

    #def get_node_index_by_node_id(self, node_id, msg=''):
        #return searchsorted(node_id, self.node_id)
        #i_too_large = where(self.node_id[-1] < node_id)[0]
        #if len(i_too_large):
            #raise RuntimeError('Cannot find GRID %s, %s' % (node_id[i_too_large], msg))
        #return self._get_sorted_index(self.node_id, node_id, 'node_id in GRID', check=True)

    def get_node_id_by_node_index(self, i=None, msg=''):
        return self.node_id[i]

    def get_node_index_by_node_id(self, node_id=None, msg=''):
        #assert msg != ''
        i = self._get_sorted_index(self.node_id, node_id, 'node_id',
                                   'node_id in GRID%s' % msg, check=True)
        return i

    def get_node_index_by_cp(self, cp=None, i=None):
        """Find all the j-indicies where cp=cpi for some given subset of i-indicies"""
        return self._get_index_by_param('cp', self.cp, cp, i)

    def get_node_index_by_cd(self, cd=None, i=None):
        """Find all the j-indicies where cd=cdi for some given subset of i-indicies"""
        return self._get_index_by_param('cd', self.cd, cd, i)

    def get_node_index_by_seid(self, seid=None, i=None):
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

    def get_position_by_node_id(self, node_id=None, msg=''):
        i = self.get_node_index_by_node_id(node_id, msg=msg)
        return self.get_position_by_node_index(i, msg=msg)

    def get_position_by_node_index(self, i=None, msg=''):
        """
        in the global frame
        """
        if i is None:
            xyz = self.xyz.copy()
            #n = slice(None, None)
            n = arange(self.n)
        else:
            n = i
            xyz = self.xyz[n, :].copy()

        cpn = self.cp[n]
        i = where(cpn != 0)[0]
        if len(i):
            n2 = n[i]
            cps = set(list(unique(cpn)))
            for cp in cps:
                #print(self.model.coords)
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
                xyz[j, :] = self.model.coords.get_global_position_by_xyz(xyzi, cp)

        #assert len(node_ids) == len(cpn), 'n1=%s n2=%s'  %(len(node_ids), len(cpn))
        return xyz

    def get_position_wrt_by_node_id(self, node_id=None, coord_id=0):
        msg = ', which is required by GRIDs'
        i = self.get_node_index_by_node_id(node_id, msg=msg)
        self.get_position_wrt_by_node_index(i, coord_id, msg)

    def get_position_wrt_by_node_index(self, i=None, coord_id=0, msg=''):
        if coord_id is None:
            return self.xyz[i, :]
        assert isinstance(coord_id, int), type(coord_id)
        raise NotImplementedError()

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('GRID', self.n))
        return msg

    def write_card(self, bdf_file, node_id=None, size=8, is_double=False, write_header=True):
        #self.model.log.debug('GRID node_id = %s' % node_id)
        #self.model.log.debug('GRID self.node_id = %s' % self.node_id)
        i = self.get_node_index_by_node_id(node_id, msg='GRID.write_card')
        return self.write_card_by_index(bdf_file, i, size, is_double, write_header)

    def write_card_by_index(self, bdf_file, i=None, size=8, is_double=False, write_header=True):
        """
        Write the BDF cards

        Parameters
        ----------
        bdf_file : file
            a file object
        i : List[int] (default=None -> all)
            the indicies
        size : int; default=8
            the field width (8/16)
        is_double: bool; default=False
            is this double precision
        write_header : bool; default=True
            should the card marker be written
        """
        if i is None:
            i = slice(None, None)
        if self.n:
            if write_header:
                bdf_file.write('$GRID\n')
            if max(self.node_id.max(), self.cd.max(), self.cp.max()) > self.model.max_int:
                size = 16
            #self.model.log.debug('GRID i = %s' % i)

            # default to the GRDSET defaults
            #cp0 = self.model.grdset.cp
            #cd0 = self.model.grdset.cd
            #ps0 = self.model.grdset.ps
            #seid0 = self.model.grdset.seid

            # default to the GRID defaults
            #cp0 = 0
            #cd0 = 0
            ps0 = -1
            seid0 = 0
            blank = ' ' * 8 if size == 8 else ' ' * 16
            Cp = [cpi if cpi != 0 else blank for cpi in self.cp[i]]
            Cd = [cdi if cdi != 0 else blank for cdi in self.cd[i]]
            Ps = [psi if psi != ps0 else blank for psi in self.ps[i]]
            Seid = [seidi if seidi != seid0 else blank for seidi in self.seid[i]]
            if size == 8:
                for (nid, cp, xyz, cd, ps, seid) in zip(
                        self.node_id, Cp, self.xyz[i, :], Cd, Ps, Seid):
                    msg = ('GRID    %8i%8s%s%s%s%8s%8s%s\n' % (
                        nid, cp,
                        print_float_8(xyz[0]),
                        print_float_8(xyz[1]),
                        print_float_8(xyz[2]),
                        cd, ps, seid)).rstrip() + '\n'

                    bdf_file.write(msg)
            else:
                if is_double:
                    for (nid, cp, xyz, cd, ps, seid) in zip(
                            self.node_id, Cp, self.xyz[i, :], Cd, Ps, Seid):
                        msg = (('GRID*   %16i%16s%16s%16s\n'
                                '*       %16s%16s%16s%16s\n' % (
                                    nid, cp,
                                    print_scientific_double(xyz[0]),
                                    print_scientific_double(xyz[1]),
                                    print_scientific_double(xyz[2]),
                                    cd, ps, seid))).rstrip() + '\n'
                        bdf_file.write(msg)
                else:
                    for (nid, cp, xyz, cd, ps, seid) in zip(
                            self.node_id, Cp, self.xyz[i, :], Cd, Ps, Seid):
                        msg = (('GRID*   %16i%16s%16s%16s\n'
                                '*       %-8s%16s%16s%16s\n' % (
                                    nid, cp,
                                    print_float_16(xyz[0]),
                                    print_float_16(xyz[1]),
                                    print_float_16(xyz[2]),
                                    cd, ps, seid))).rstrip() + '\n'
                        bdf_file.write(msg)
            #return self.comment + msg.rstrip() + '\n'


    def __repr__(self):
        msg = "<GRID>\n"
        msg += '  nGRID = %i' % self.n
        return msg

    def __getitem__(self, node_id):
        #raise NotImplementedError('Grid getitem')
        #return self.slice_by_index(i)
        return self.slice_by_node_id(node_id)

    def slice_by_node_id(self, node_id=None):
        #self.model.log.debug('self.node_id = %s' % self.node_id)
        #self.model.log.debug('node_id = %s' % node_id)
        #node_id = slice_to_iter(node_id)
        i = where(self.node_id == node_id)[0]
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        #self.model.log.debug('i = %s; type=%s' % (i, type(i)))
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
        #n = self.n
    else:
        n = len(i)
    return i, n

