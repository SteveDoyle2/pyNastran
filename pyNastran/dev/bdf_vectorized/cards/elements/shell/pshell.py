from io import StringIO
import numpy as np
from numpy import array, zeros, searchsorted, unique, where

from pyNastran.dev.bdf_vectorized.utils import slice_to_iter
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property


class PSHELL(Property):
    """
    +--------+-------+------+--------+------+----------+------+------+---------+
    |   1    |   2   |   3  |    4   |  5   |    6     |   7  |  8   |    9    |
    +========+=======+======+========+======+==========+======+======+=========+
    | PSHELL |  PID  | MID1 |   T    | MID2 | 12I/T**3 | MID3 | TS/T |   NSM   |
    +--------+-------+------+--------+------+----------+------+------+---------+
    |        |  Z1   |  Z2  |  MID4  |      |          |      |      |         |
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

    @property
    def material_id(self):
        return self.material_ids[:, 0]

    @property
    def material_id2(self):
        return self.material_ids[:, 1]

    @property
    def material_id3(self):
        return self.material_ids[:, 2]

    @property
    def material_id4(self):
        return self.material_ids[:, 3]

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.model.log.debug('allocate PSHELL')
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.property_id = zeros(ncards, 'int32')

            #: Material identification number for bending
            #self.material_id = zeros(ncards, 'int32')
            #self.material_id2 = zeros(ncards, 'int32')
            #self.material_id3 = zeros(ncards, 'int32')
            #self.material_id4 = zeros(ncards, 'int32')
            self.material_ids = zeros((ncards, 4), 'int32')

            self.thickness = zeros(ncards, float_fmt)
            # .. todo:: poor name
            #: ..math:: I = \frac{12I}{t^3} I_{plate}
            #: Scales the moment of interia of the element based on the
            #: moment of interia for a plate
            self.twelveIt3 = zeros(ncards, float_fmt)

            self.tst = zeros(ncards, float_fmt)

            #: Non-structural Mass
            self.nsm = zeros(ncards, float_fmt)
            self.z1 = zeros(ncards, float_fmt)
            self.z2 = zeros(ncards, float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        pid = integer(card, 1, 'property_id')
        if comment:
            self.set_comment(pid, comment)
        i = self.i
        self.property_id[i] = pid
        self.thickness[i] = double(card, 3, 'thickness')

        self.twelveIt3[i] = double_or_blank(card, 5, '12*I/t^3', 1.0)
        self.tst[i] = double_or_blank(card, 7, 'ts/t', 0.833333)
        self.nsm[i] = double_or_blank(card, 8, 'nsm', 0.0)

        t_over_2 = self.thickness[i] / 2.
        self.z1[i] = double_or_blank(card, 9, 'z1', -t_over_2)
        self.z2[i] = double_or_blank(card, 10, 'z2', t_over_2)
        self.material_ids[i, :] = [
            integer(card, 2, 'material_id'),
            integer_or_blank(card, 4, 'material_id2', -1),
            integer_or_blank(card, 6, 'material_id3', -1),
            integer_or_blank(card, 11, 'material_id4', -1)
        ]
        #self.model.log.debug(self.material_ids[i, :])

        #ii = np.array([i])
        #file_obj = StringIO()
        #file_obj.write('<PSHELL object> n=%s\n' % self.n)
        #self.write_card_by_index(file_obj, i=ii)
        #print(file_obj.getvalue())


        assert len(card) <= 12, 'len(PSHELL card) = %i\ncard=%s' % (len(card), card)
        self.i += 1
        #self.model.log.debug('self.i = %s' % self.i)

    def build(self):
        if self.n:
            # sort the NDARRAYs so we can use searchsorted
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_ids = self.material_ids[i, :]
            self.thickness = self.thickness[i]
            self.twelveIt3 = self.twelveIt3[i]
            self.tst = self.tst[i]
            self.nsm = self.nsm[i]
            self.z1 = self.z1[i]
            self.z2 = self.z2[i]

            if len(unique(self.property_id)) != len(self.property_id):
                pids = self.property_id
                raise RuntimeError('There are duplicate PSHELL IDs...pids=%s' % pids)
        else:
            self.property_id = array([], dtype='int32')
            self.material_ids = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            pid_map = maps['property']
            mid_map = maps['material']
            for i, pid, mids in enumerate(zip(self.property_id, self.material_ids)):
                pid2 = pid_map[pid]
                mids2 = [mid_map[mid] for mid in mids]
                self.property_id[i] = pid2
                self.material_ids[i, :] = mids2

    def write_card_by_index(self, bdf_file, size=8, is_double=True, i=None):
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
        mid2 = [midi if midi > 0 else '' for midi in self.material_id2[i]]
        mid3 = [midi if midi > 0 else '' for midi in self.material_id3[i]]
        mid4 = [midi if midi > 0 else '' for midi in self.material_id4[i]]
        nsm = ['' if nsmi == 0.0 else nsmi for nsmi in self.nsm[i]]
        tst = ['' if tsti == 0.833333 else tsti for tsti in self.tst[i]]
        TwelveIt3 = ['' if tw == 1.0 else tw for tw in self.twelveIt3[i]]

        to2 = self.thickness[i] / 2
        z1 = ['' if z1i == -to2[j] else z1i for j, z1i in enumerate(self.z1[i])]
        z2 = ['' if z2i == to2[j] else z2i for j, z2i in enumerate(self.z2[i])]

        for (pid, mid1i, t, mid2i, twelveIt3, mid3i, tsti, nsmi, z1i, z2i, mid4i) in zip(
                self.property_id[i], self.material_id[i], self.thickness[i], mid2,
                TwelveIt3, mid3, tst, nsm, z1, z2, mid4):
            if pid in self._comments:
                bdf_file.write(self._comments[pid])

            #list_fields = ['PSHELL', self.pid, mid1, self.t, mid2,
                           #twelveIt3, mid3, tst, nsm, z1, z2, mid4]
            card = ['PSHELL', pid, mid1i, t, mid2i, twelveIt3, mid3i,
                    tsti, nsmi, z1i, z2i, mid4i]
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
            #self.model.log.debug('i = %s' % i)
            nsm = self.nsm[i]
        return nsm[0] if int_flag else nsm

    def get_thickness_by_property_id(self, property_id=None):
        """
        Gets the thickness of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        int_flag = True if isinstance(property_id, int) else False
        if property_id is None:
            thickness = self.thickness
        else:
            i = self.get_property_index_by_property_id(property_id)
            thickness = self.thickness[i]
        return thickness[0] if int_flag else thickness

    def get_mass_per_area_by_property_id(self, property_id=None):
        """
        Gets the mass per area of the PHSELLs.

        :param property_id: the property IDs to consider (default=None -> all)
        """
        #mass_per_area = self.nsm + self.Rho() * self.t
        if property_id is None:
            thickness = self.thickness
            nsm = self.nsm
        else:
            i = self.get_property_index_by_property_id(property_id)
            thickness = self.thickness[i]
            nsm = self.nsm[i]

        density = self.get_density_by_property_id(property_id)
        return nsm + density * thickness

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
        property_id = slice_to_iter(property_id)[0]
        #self.model.log.debug('looking for %s property_ids' % str(property_id))
        #self.model.log.debug(property_id)
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
        obj.material_ids = self.material_ids[i, :]
        obj.thickness = self.thickness[i]
        obj.twelveIt3 = self.twelveIt3[i]
        obj.tst = self.tst[i]
        obj.nsm = self.nsm[i]
        obj.z1 = self.z1[i]
        obj.z2 = self.z2[i]

        #print(obj)
        return obj

    def __repr__(self):
        file_obj = StringIO()
        file_obj.write('<PSHELL object> n=%s\n' % self.n)
        self.write_card(file_obj)
        return file_obj.getvalue()
