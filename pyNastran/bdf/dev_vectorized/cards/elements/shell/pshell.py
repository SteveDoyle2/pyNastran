from six.moves import zip, StringIO
from numpy import array, zeros, searchsorted, unique, where

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank, double,
    double_or_blank)

from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

class PSHELL(Property):
    """
    +--------+-------+------+--------+------+----------+------+------+---------+
    | PSHELL |  PID  | MID1 |   T    | MID2 | 12I/T**3 | MID3 | TS/T |   NSM   |
    +--------+-------+------+--------+------+----------+------+------+---------+
    |        |  Z1   |  Z2  |  MID4  |
    +--------+-------+------+--------+------+----------+------+------+---------+
    | PSHELL | 41111 |  1   | 1.0000 |  1   |          |   1  |      | 0.02081 |
    +--------+-------+------+--------+------+----------+------+------+---------+
    """
    type = 'PSHELL'

    def __init__(self, model):
        """
        Defines the PSHELL object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Property.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.model.log.debug('allocate PSHELL')
            self.n = ncards
            float_fmt = self.model.float_fmt
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

    def add_card(self, card, comment=''):
        pid = integer(card, 1, 'property_id')
        if comment:
            self._comments[pid] = comment
        i = self.i
        #self.model.log.debug('adding %s; i=%s' %(card, i))
        self.property_id[i] = pid
        self.material_id[i] = integer(card, 2, 'material_id')
        self.thickness[i] = double(card, 3, 'thickness')

        #: Material identification number for bending
        self.material_id2[i] = integer_or_blank(card, 4, 'material_id2', -1)

        # .. todo:: poor name
        #: ..math:: I = \frac{12I}{t^3} I_{plate}
        #: Scales the moment of interia of the element based on the
        #: moment of interia for a plate
        self.twelveIt3[i] = double_or_blank(card, 5, '12*I/t^3', 1.0)

        self.material_id3[i] = integer_or_blank(card, 6, 'material_id3', -1)
        self.tst[i] = double_or_blank(card, 7, 'ts/t', 0.833333)

        #: Non-structural Mass
        self.nsm[i] = double_or_blank(card, 8, 'nsm', 0.0)

        tOver2 = self.thickness[i] / 2.
        self.z1[i] = double_or_blank(card, 9, 'z1', -tOver2)
        self.z2[i] = double_or_blank(card, 10, 'z2', tOver2)
        self.material_id4[i] = integer_or_blank(card, 11, 'material_id4', -1)

        #if self.material_id2 is None:
        #    assert self.material_id3 is None
        #else: # material_id2 is defined
        #    #print "self.material_id2 = ",self.material_id2
        #    assert self.material_id2 >= -1
        #    #assert self.material_id3 >   0

        #if self.material_id is not None and self.material_id2 is not None:
        #    assert self.material_id4==None
        assert len(card) <= 12, 'len(PSHELL card) = %i\ncard=%s' % (len(card), card)
        self.i += 1
        #self.model.log.debug('self.i = %s' % self.i)

    def build(self):
        if self.n:
            # sort the NDARRAYs so we can use searchsorted
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            #print('self.t1', self.thickness)
            self.thickness = self.thickness[i]
            #print('self.t2', self.thickness)
            self.material_id2 = self.material_id2[i]
            self.twelveIt3 = self.twelveIt3[i]
            self.material_id3 = self.material_id3[i]
            self.tst = self.tst[i]
            self.nsm = self.nsm[i]
            self.z1 = self.z1[i]
            self.z2 = self.z2[i]
            self.material_id4 = self.material_id4[i]

            if len(unique(self.property_id)) != len(self.property_id):
                pids = self.property_id
                raise RuntimeError('There are duplicate PSHELL IDs...pids=%s' % pids)
        else:
            self.property_id = array([], dtype='int32')
            self.material_id = array([], dtype='int32')

    def write_card(self, bdf_file, size=8, property_id=None):
        """
        Writes the PSHELL properties.

        Parameters
        ----------
        bdf_file: file
            a file object
        size : int; default=8
            the bdf field size (8/16)
        property_id : ???; default=None -> all
            the property_ids to write
        """
        if self.n:
            i = self.get_property_index_by_property_id(property_id)
            #print('i = %s' % i)
            #self.model.log.debug('material_id2 = %s' % self.material_id2[i])
            #self.model.log.debug('material_id3 = %s' % self.material_id3[i])
            #self.model.log.debug('material_id4 = %s' % self.material_id4[i])
            #self.model.log.debug('material_id2.shpae = %s' % str(self.material_id2[i].shape))

            Mid2 = [midi if midi > 0 else '' for midi in self.material_id2[i]]
            Mid3 = [midi if midi > 0 else '' for midi in self.material_id3[i]]
            Mid4 = [midi if midi > 0 else '' for midi in self.material_id4[i]]
            Nsm = ['' if nsmi == 0.0 else nsmi for nsmi in self.nsm[i]]
            Tst = ['' if tsti == 0.833333 else tsti for tsti in self.tst[i]]
            TwelveIt3 = ['' if tw == 1.0 else tw for tw in self.twelveIt3[i]]

            to2 = self.thickness[i] / 2
            Z1 = ['' if z1i == -to2[j] else z1i for j, z1i in enumerate(self.z1[i])]
            Z2 = ['' if z2i == to2[j] else z2i for j, z2i in enumerate(self.z2[i])]

            for (pid, mid, t, mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4) in zip(
                    self.property_id[i], self.material_id[i], self.thickness[i], Mid2,
                    TwelveIt3, Mid3, Tst, Nsm, Z1, Z2, Mid4):
                if pid in self._comments:
                    bdf_file.write(self._comments[pid])
                card = ['PSHELL', pid, mid, t, mid2, twelveIt3, mid3,
                        tst, nsm, z1, z2, mid4]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def get_nonstructural_mass_by_property_id(self, property_id=None):
        """
        Gets the nonstructural mass of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        #print('get_nonstructural_mass; pids = %s' % property_ids)
        int_flag = True if isinstance(property_id, int) else False
        if property_id is None:
            nsm = self.nsm
        else:
            i = self.get_property_index_by_property_id(property_id)
            #print('i = %s' % i)
            nsm = self.nsm[i]
        return nsm[0] if int_flag else nsm

    def get_thickness_by_property_id(self, property_id=None):
        """
        Gets the thickness of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        int_flag = True if isinstance(property_id, int) else False
        if property_id is None:
            t = self.thickness
        else:
            i = self.get_property_index_by_property_id(property_id)
            t = self.thickness[i]
        return t[0] if int_flag else t

    def get_mass_per_area_by_property_id(self, property_id=None):
        """
        Gets the mass per area of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        #mass_per_area = self.nsm + self.Rho() * self.t
        if property_id is None:
            t = self.thickness
            nsm = self.nsm
        else:
            i = self.get_property_index_by_property_id(property_id)
            t = self.thickness[i]
            nsm = self.nsm[i]

        density = self.get_density_by_property_id(property_id)
        return nsm + density * t

    def get_density_by_property_id(self, property_id=None):
        """
        Gets the density of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        int_flag = True if isinstance(property_id, int) else False
        material_id = self.material_id
        j = where(material_id == 0)[0]
        material_ids = material_id.copy()
        material_ids[j] = self.material_id2

        if property_id is not None:
            i = self.get_property_index_by_property_id(property_id)
            material_ids = material_ids[i]

        density = self.model.materials.get_density_by_material_id(material_ids)
        return density[0] if int_flag else density

    def get_material_id_by_property_id(self, property_id=None):
        """
        Gets the material IDs of the PSHELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        if property_id is None:
            mid = self.material_id
        else:
            i = self.get_property_index_by_property_id(property_id)
            mid = self.material_id[i]
        return mid

    def __getitem__(self, property_id):
        property_id, int_flag = slice_to_iter(property_id)
        #print('looking for %s property_ids' % str(property_id))
        i = searchsorted(self.property_id, property_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = self._validate_slice(i)

        obj = PSHELL(self.model)
        try:
            obj.n = len(i)
        except TypeError:
            msg = 'i=%s type(i)=%s' % (i, type(i))
            raise TypeError(msg)

        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.thickness = self.thickness[i]
        obj.material_id2 = self.material_id2[i]
        obj.twelveIt3 = self.twelveIt3[i]
        obj.material_id3 = self.material_id3[i]
        obj.tst = self.tst[i]
        obj.nsm = self.nsm[i]
        obj.z1 = self.z1[i]
        obj.z2 = self.z2[i]
        obj.material_id4 = self.material_id4[i]
        #print(obj)
        return obj

    def __repr__(self):
        f = StringIO()
        f.write('<PSHELL object> n=%s\n' % self.n)
        self.write_card(f)
        return f.getvalue()
