from numpy import array, zeros, argsort, concatenate, searchsorted, unique, where, nan, arange

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, double,
    double_or_blank, string_or_blank)


class PSHELL(object):
    type = 'PSHELL'
    def __init__(self, model):
        """
        Defines the PSHELL object.

        :param self: the PSHELL object
        :param model: the BDF object
        :param cards: the list of PSHELL cards
        """
        self.model = model
        self.n = 0
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
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.thickness = zeros(ncards, float_fmt)

            self.material_id2 = zeros(ncards, 'int32')
            self.twelveIt3 = zeros(ncards, float_fmt)
            self.material_id3 = zeros(ncards, 'int32')

            self.tst = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)
            self.z1 = zeros(ncards, float_fmt)
            self.z2 = zeros(ncards, float_fmt)
            self.material_id4 = zeros(ncards, 'int32')

            # ..todo:: incomplete
            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'property_id')
                self.material_id[i] = integer(card, 2, 'material_id')
                self.thickness[i] = double(card, 3, 'thickness')

                #: Material identification number for bending
                self.material_id2[i] = integer_or_blank(card, 4, 'material_id2', -1)
                #: Scales the moment of interia of the element based on the
                #: moment of interia for a plate
                #:
                #: ..math:: I = \frac{12I}{t^3} I_{plate}
                self.twelveIt3[i] = double_or_blank(card, 5, '12*I/t^3', 1.0)  # poor name
                self.material_id3[i] = integer_or_blank(card, 6, 'material_id3', -1)
                self.tst[i] = double_or_blank(card, 7, 'ts/t', 0.833333)
                #: Non-structural Mass
                self.nsm[i] = double_or_blank(card, 8, 'nsm', 0.0)

                tOver2 = self.thickness[i] / 2.
                self.z1[i] = double_or_blank(card, 9,  'z1', -tOver2)
                self.z2[i] = double_or_blank(card, 10, 'z2',  tOver2)
                self.material_id4[i] = integer_or_blank(card, 11, 'material_id4', -1)

                #if self.material_id2 is None:
                #    assert self.material_id3 is None
                #else: # material_id2 is defined
                #    #print "self.material_id2 = ",self.material_id2
                #    assert self.material_id2 >= -1
                #    #assert self.material_id3 >   0

                #if self.material_id is not None and self.material_id2 is not None:
                #    assert self.material_id4==None
                assert len(card) <= 12, 'len(PSHELL card) = %i' % len(card)

            # nan is float specific
            #self.material_id[where(self.material_id2 == -1)[0]] = nan

            # sort the NDARRAYs so we can use searchsorted
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.thickness = self.thickness[i]
            self.material_id2 = self.material_id2[i]
            self.twelveIt3 = self.twelveIt3[i]
            self.material_id3 = self.material_id3[i]
            self.tst = self.tst[i]
            self.nsm = self.nsm[i]
            self.z1 = self.z1[i]
            self.z2 = self.z2[i]
            self.material_id4 = self.material_id4[i]
        
            if len(unique(self.property_id)) != len(self.property_id):
                raise RuntimeError('There are duplicate PSHELL IDs...')
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, property_ids=None):
        """
        Writes the PSHELL properties.
        
        :param self:  the PSHELL object
        :param f:     file object
        :param size:  the bdf field size (8/16; default=8)
        :param property_ids:  the property_ids to write (default=None -> all)
        """
        if self.n:
            if property_ids is None:
                i = arange(len(self.property_id))
                #for (pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4) in zip(
                        #self.property_id, self.mid, self.thickness, self.mid2,
                        #self.twelveIt3, self.mid3, self.tst, self.nsm, self.z1,
                        #self.z2, self.mid4):
                    #card = ['PSHELL', pid, mid, t]
                    #f.write(print_card(card, size=size))
            else:
                i = searchsorted(self.property_id, property_ids)

            Mid2 = [midi if midi > 0 else '' for midi in self.material_id2[i]]
            Mid3 = [midi if midi > 0 else '' for midi in self.material_id3[i]]
            Mid4 = [midi if midi > 0 else '' for midi in self.material_id4[i]]
            Nsm       = ['' if nsmi == 0.0      else nsmi for nsmi in self.nsm[i]]
            Tst       = ['' if tsti == 0.833333 else tsti for tsti in self.tst[i]]
            TwelveIt3 = ['' if tw   == 1.0      else tw   for tw   in self.twelveIt3[i]]

            to2 = self.thickness[i] / 2
            Z1 = ['' if z1i == -to2[j] else z1i for j, z1i in enumerate(self.z1[i])]
            Z2 = ['' if z2i ==  to2[j] else z2i for j, z2i in enumerate(self.z2[i])]

            for (pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4) in zip(
                    self.property_id[i], self.material_id[i], self.thickness[i], Mid2,
                    TwelveIt3, Mid3, Tst, Nsm,  Z1, Z2, Mid4):
                card = ['PSHELL', pid, mid, t, mid2, twelveIt3, mid3,
                        tst, nsm, z1, z2, mid4]
                f.write(print_card(card, size=size))

    def get_thickness(self, property_ids=None):
        """
        Gets the thickness of the PHSELLs.
        
        :param self: the PSHELL object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        if property_ids is None:
            t = self.thickness
        else:
            t = self.thickness[searchsorted(property_ids, self.property_id)]
        return t

    def get_mid(self, property_ids=None):
        """
        Gets the material IDs of the PSHELLs.
        
        :param self: the PSHELL object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        if pids is None:
            mid = self.mid
        else:
            mid = self.mid[searchsorted(property_ids, self.property_id)]
        return t