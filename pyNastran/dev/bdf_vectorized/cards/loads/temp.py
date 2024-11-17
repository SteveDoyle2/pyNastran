from collections import defaultdict
from numpy import zeros, array

from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.field_writer_8 import print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank)

from .loads import VectorizedCardDict

class TEMP:
    type = 'TEMP'
    def __init__(self, model):
        self.model = model
        self.temp_id = 0.
        self.n = 0
        self.default = 0.
        self.is_default = False
        self._i = 0

    def __len__(self):
        return self.n + 1 if self.is_default else 0

    def allocate(self, card_count):
        n = card_count['TEMP']
        float_fmt = self.model.float_fmt
        self.temp_id = 0
        #self.node_id = zeros(n, dtype='int32')
        self.node_id = self.model.grid.node_id
        self.temperature = zeros(n, dtype=float_fmt)

    def add(self, temp_id, card, comment):
        self.temp_id = temp_id
        assert len(card) % 2 == 0, card
        for i in range(1, len(card), 2):
            node_id = integer(card, i, 'node_id')
            temperature = double(card, i+1, 'temperature')
            #self.node_id[self.i] = node_id
            self.temperature[self.i] = temperature
            self.i += 1

    def add_default(self, temp_id, temperature, comment):
        self.temp_id = temp_id
        self.default = temperature
        self.is_default = True

    def write_card(self, bdf_file, size=8, is_double=False):
        #nids = ['%8i' % node_id for node_id in self.node_id]
        nleftover = self.n // 3 * 3
        if size == 8:
            temps = [print_float_8(tempi) % temp for temp in self.temperature]
            for i in range(0, self.n, 3):
                bdf_file.write('TEMP    %8i%8i%s %8i%s%8i%s' % (self.temp_id,
                                                                self.node_id[i], print_float_8(temps[i]),
                                                                self.node_id[i+1], print_float_8(temps[i+1]),
                                                                self.node_id[i+2], print_float_8(temps[i+2])))
            for i in range(nleftover, self.n - nleftover):
                bdf_file.write('TEMP    %8i%s' % (self.temp_id, print_float_8(temps[i])))

            if self.is_default:
                bdf_file.write('TEMPD   %8i%s\n' % (self.temp_id, print_float_8(self.default)))
        else:
            temps = [print_float_16(tempi) % temp for temp in self.temperature]
            for i in range(0, self.n, 3):
                bdf_file.write('TEMP*   %16i%16i%s\n' % (self.temp_id,
                                                         self.node_id[i], print_float_16(temps[i])))
                bdf_file.write('%8s%16i%s%16i%s' % ('*       ',
                                                    self.node_id[i+1], print_float_16(temps[i+1]),
                                                    self.node_id[i+2], print_float_16(temps[i+2])))
            for i in range(nleftover, self.n - nleftover):
                bdf_file.write('TEMP*   %16i%16i%s\n*\n' % (self.temp_id, self.node_id[i], print_float_16(temps[i])))
            if self.is_default:
                bdf_file.write('TEMPD*  %16i%s\n' % (self.temp_id, print_float_16(self.default)))

class TEMPP1:
    def __init__(self, model):
        """
        TEMPP1 SID  EID1 TBAR TPRIME T1   T2   NA    NA
               EID2 EID3 EID4 EID5   EID6 EID7  etc.

        TEMPP1 EID2 THRU EIDi EIDj THRU EIDk
        TEMPP1 1    THRU 10   30   THRU 61

        TEMPP1    500001    2639.107E+03-.12E+02
        """
        self.model = model
        self.n = 0
        self.i = 0

    def __len__(self):
        return self.n

    def allocate(self, card_count):
        if 'TEMPP1' in card_count:
            float_fmt = self.model.float_fmt
            n = card_count['TEMPP1']
            self.n = n
            #self.model.log.debug('TEMPP1 n=%s' % self.n)
            self.load_id = zeros(n, dtype='int32')
            self.element_id = zeros(n, dtype='int32')
            self.tbar = zeros(n, dtype=float_fmt)
            self.tprime = zeros(n, dtype=float_fmt)
            self.temp = zeros((n, 2), dtype=float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        assert self.n > 0, self.n
        i = self.i
        load_id = integer(card, 1, 'load_id')
        tbar = double(card, 3, 'Tbar')
        tprime = double(card, 4, 'Tprime')
        t1 = double_or_blank(card, 5, 'T1')
        t2 = double_or_blank(card, 6, 'T2')

        self.load_id[i] = load_id
        self.element_id[i] = integer(card, 2, 'element_id')
        self.tbar[i] = tbar
        self.tprime[i] = tprime
        self.temp[i, 0] = t1
        self.temp[i, 1] = t2
        self.i += 1

        if len(card) >= 7:
            # i must be < self.n
            eids = expand_thru(card[9:])
            for eid in eids:
                self.load_id[i] = load_id
                assert isinstance(eid, int), eid
                self.element_id[i] = eid
                self.tbar[i] = tbar
                self.tprime[i] = tprime
                self.temp[i, 0] = t1
                self.temp[i, 1] = t2
                self.i += 1
                assert self.i <= self.n

        assert len(card) <= 7, '%s; n=%s' % (card, len(card))
        #assert len(card) <= 7, len(card)
        self.eids = None

    def write_card(self, bdf_file, size=8, is_double=False):
        if self.n:
            if size == 8:
                for load_id, element_id, tbar, tprime in zip(self.load_id, self.element_id, self.tbar, self.tprime):
                    bdf_file.write('TEMPP1  %8i%8i%s%s\n' % (load_id, element_id, print_float_8(tbar), print_float_8(tprime)))
            else:
                for load_id, element_id, tbar, tprime in zip(self.load_id, self.element_id, self.tbar, self.tprime):
                    bdf_file.write('TEMPP1* %16i%16i%s%s\n*\n' % (load_id, element_id, print_float_16(tbar), print_float_16(tprime)))


class TEMPs(VectorizedCardDict):
    def __init__(self, model):
        VectorizedCardDict.__init__(self, model)
        self._objs = defaultdict(TEMP)
        self.tempp1 = model.tempp1

    def allocate(self, card_count):
        self.tempp1.allocate(card_count)

    def get_temperature_id(self):
        return array(self._objs.keys(), dtype='int32')

    def __len__(self):
        return len(self._tempp1) + len(self._objs)

    #def slice_by_temperature_id(self, load_id):
        #return

    def add_tempd(self, card, comment):
        assert (len(card) - 1) % 2 == 0, card
        for i in range(1, len(card), 2):
            temp_id = integer(card, i, 'temp_id')
            temperature = double(card, i+1, 'temperature')
            self._objs[temp_id].add_default(temp_id, temperature, comment)
        self.model.log.debug('TEMPs keys=%s' % self._objs.keys())

    def add_temp(self, card, comment):
        temp_id = integer(card, 1, 'temp_id')
        load = self._objs[temp_id].add(temp_id, card, comment)
        load.add_from_bdf(card, comment)
        self._objs[load.load_id].append(load)

    def add_tempp1(self, card, comment):
        self.tempp1.add(card, comment)

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None, sort_data=False):
        #self.model.log.debug('TEMPs keys=%s' % self._objs.keys())
        for load_id, load in sorted(self._objs.items()):
            #self.model.log.debug('TEMPs write_card load_id=%s' % load_id)
            load.write_card(bdf_file, size=size, is_double=is_double)
        self.tempp1.write_card(bdf_file, size=size, is_double=is_double)
