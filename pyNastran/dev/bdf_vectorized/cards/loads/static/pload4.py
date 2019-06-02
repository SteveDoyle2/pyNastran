import numpy as np
from numpy import arange, array, zeros, searchsorted, unique, full, nan, where

from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_string_or_blank, string_or_blank)

from pyNastran.dev.bdf_vectorized.cards.loads.vectorized_load import VectorizedLoad

class PLOAD4(VectorizedLoad):
    type = 'PLOAD4'

    def __init__(self, model):
        """
        Defines the PLOAD4 object.

        Parameters
        ----------
        model : BDF
           the BDF object

        """
        VectorizedLoad.__init__(self, model)
        del self._comments
        #self.model = model
        #self.n = 0
        #self.i = 0
        #self._cards = []
        self._comments = []

        #self._load_id = []
        #self._element_ids = []
        #self._p = []

    #def __getitem__(self, i):
        #unique_lid = unique(self.load_id)
        ##print("force", unique_lid, i)
        ##if len(i):
        #f = PLOAD4(self.model)
        #f.load_id = self.load_id[i]
        #f.element_id = self.element_id[i]
        #f.p = self.p[i]
        #f.n = len(i)
        #return f
        ##raise RuntimeError('len(i) = 0')

    def __getitem__(self, load_id):
        self.model.log.debug('load_id = %s' % load_id)
        assert load_id > 0, load_id
        i = where(load_id == self.load_id)[0]
        return self.slice_by_index(i)

    def __mul__(self, value):
        raise NotImplementedError()
        #f = PLOAD4(self.model)
        #f.load_id = self.load_id
        #f.element_id = self.element_id
        #f.p = self.p[i]
        #f.n = self.n
        #return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card, comment=None):
        i = self.i
        self._comments.append(comment)
        #element_ids = {}
        #for i in range(ncards):
            #element_ids[i] = []


        self.load_id[i] = integer(card, 1, 'load_id')
        eid = integer(card, 2, 'element_id')
        self.element_ids[i] = eid
        p1 = double_or_blank(card, 3, 'p1', 0.0)
        p = [p1,
             double_or_blank(card, 4, 'p2', p1),
             double_or_blank(card, 5, 'p3', p1),
             double_or_blank(card, 6, 'p4', p1)]
        self.pressures[i, :] = p

        self.element_ids[i] = [eid]
        if(integer_string_or_blank(card, 7, 'g1/THRU') == 'THRU' and
           integer_or_blank(card, 8, 'eid2')):  # plates
            eid2 = integer(card, 8, 'eid2')
            if eid2:
                self.element_ids[i] = list(unique(expand_thru([eid, 'THRU', eid2],
                                           set_fields=False, sort_fields=False)))
            #self.g1 = None
            #self.g34 = None
        else:
            #: used for CPENTA, CHEXA
            self.element_ids[i] = [eid]
            #: used for solid element only
            self.g1[i] = integer_or_blank(card, 7, 'g1', -1)
            #: g3/g4 - different depending on CHEXA/CPENTA or CTETRA
            self.g34[i] = integer_or_blank(card, 8, 'g34', -1)

        #: Coordinate system identification number. See Remark 2.
        #: (Integer >= 0;Default=0)
        self.cid[i] = integer_or_blank(card, 9, 'cid', 0)
        self.nvector[i, :] = [
            double_or_blank(card, 10, 'N1'),
            double_or_blank(card, 11, 'N2'),
            double_or_blank(card, 12, 'N3'), ]
        self.surf_or_line[i] = string_or_blank(card, 13, 'sorl', 'SURF')
        self.line_load_dir[i] = string_or_blank(card, 14, 'ldir', 'NORM')
        assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1


    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.load_id = zeros(ncards, 'int32')
            #self.element_id = zeros(ncards, 'int32')
            self.pressures = zeros((ncards, 4), float_fmt)

            self.element_ids = {}
            for i in range(ncards):
                self.element_ids[i] = []

            self.g1 = full(ncards, nan, 'int32')
            self.g34 = full(ncards, nan, 'int32')
            self.line_load_dir = full(ncards, nan, '|U4')
            self.surf_or_line = full(ncards, nan, '|U4')
            self.cid = zeros(ncards, dtype='int32')
            self.nvector = zeros((ncards, 3), dtype=float_fmt)
        else:
            self.load_id = array([], dtype='int32')

    def build(self):
        if self.n:
            float_fmt = self.model.float_fmt

            #self.load_id = zeros(ncards, 'int32')
            ##: Element ID
            #self.element_id = zeros(ncards, 'int32')
            ##: Property ID
            #self.pressures = zeros((ncards, 4), 'int32')

            i = self.load_id.argsort()
            #self.element_id = self.element_id[i]
            self.pressures = self.pressures[i, :]
            #self.node_ids = self.node_ids[i, :]

            #element_ids = {}
            #for j in range(ncards):
                #element_ids[j] = element_ids[i[j]]

            self.g1 = self.g1[i]
            self.g34 = self.g34[i]
            self.line_load_dir = self.line_load_dir[i]
            self.surf_or_line = self.surf_or_line[i]
            self.cid = self.cid[i]
            self.nvector = self.nvector[i, :]
            #self.n += len(eids)

            self.load_cases = {}
            self.load_ids = unique(self.load_id)

    def get_load_ids(self):
        return self.load_ids

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOAD4', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            if load_id is None:
                i = arange(self.n)
            else:
                i = np.where(load_id == self.load_id)[0]

            #self.model.log.debug('i = %s' % i)
            #n = [None, None, None]
            #surf_or_line = None
            #line_load_dir = None
            cid = ['' if cidi == 0 else cidi for cidi in self.cid[i]]
            surf_or_line = ['' if surf_or_linei == 'SURF' else surf_or_linei for surf_or_linei in self.surf_or_line[i]]
            line_load_dir = ['' if line_load_diri == 'NORM' else line_load_diri for line_load_diri in self.line_load_dir[i]]
            for (ii, load_idi, p, n, cidi, surf_or_linei, line_load_diri) in zip(i,
                                                             self.load_id[i], self.pressures[i, :],
                                                             self.nvector[i, :], cid, surf_or_line, line_load_dir):
                #self.model.log.debug('ii = %s' % ii)
                element_id = self.element_ids[ii]
                #self.model.log.debug('element_id = %s' % element_id)
                #self.model.log.debug('p = %s\n' % p)
                #eidi = element_id[0]

                eidi = element_id[0]
                if len(element_id) == 1:
                    thru = [None, None]
                else:
                    thru = ['THRU', element_id[-1]]

                if p.max() == p.min():
                    card = ['PLOAD4', load_idi, eidi, p[0], None, None, None
                            ] + thru + [cidi] + list(n) + [surf_or_linei, line_load_diri]
                else:
                    card = ['PLOAD4', load_idi, eidi] + list(p) + thru + [
                        cidi] + list(n) + [surf_or_linei, line_load_diri]

                assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
                bdf_file.write(print_card_8(card))

    def slice_by_index(self, i):
        if i is None:
            return self
        i = np.asarray(i)
        #name = self.__class__.__name__
        #print('name = %r' % name)
        #obj = CQUAD4(self.model)
        #obj_class = type(name, (ShellElement, ), {})
        obj_class = self.__class__#.__class__
        obj = obj_class(self.model)
        #print(type(obj))

        obj.n = len(i)
        self.model.log.debug('i = %s' % i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.load_id = self.load_id[i]
        obj.pressures = self.pressures[i, :]
        obj.element_ids = {obj.element_ids[j] for j in range(i)}

        obj.g1 = self.g1[i]
        obj.g34 = self.g34[i]
        obj.line_load_dir = self.line_load_dir[i]
        obj.surf_or_line = self.surf_or_line[i]
        obj.cid = self.cid[i]
        obj.nvector = self.nvector[i, :]


        #self.load_id = zeros(ncards, 'int32')
        ##self.element_id = zeros(ncards, 'int32')
        return obj
