from itertools import count

from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_double, components_or_blank)

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.field_writer_8 import set_blank_if_default
#from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
#from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default

from numpy import array, zeros, unique, searchsorted, asarray, int64, where

#RigidElement
class RBE2:
    type = 'RBE2'

    def __init__(self, model):
        """
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |   1   |  2  |  3  |  4  |  5   |   6   |  7  |  8  |  9  |
        +=======+=====+=====+=====+======+=======+=====+=====+=====+
        |  RBE2 | EID | GN  | CM  | GM1  | GM2   | GM3 | GM4 | GM5 |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |       | GM6 | GM7 | GM8 | etc. | ALPHA |     |     |     |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        """
        self.model = model
        self.n = 0
        self.i = 0

    def shrink(self, refcheck=True):
        i = where(self.element_id == 0)[0]
        self.resize(i[0], refcheck=refcheck)

    def __iter__(self):
        pids = self.element_id
        for pid in pids:
            yield pid

    def values(self):
        pids = self.element_id
        for pid in pids:
            yield self.__getitem__(pid)

    def items(self):
        pids = self.element_id
        for pid in pids:
            yield pid, self.__getitem__(pid)

    def __getitem__(self, element_id):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.element_id, element_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        if isinstance(i, (int, int64)):
            i = [i]
        i = asarray(i)
        n = len(i)

        obj = RBE2(self.model)
        obj.element_id = self.element_id[i]
        obj.gn = self.gn[i]
        obj.cm = self.cm[i]
        obj.alpha = self.alpha[i]

        Gmi = {}
        j = 0
        for ii, gmi in self.gmi.items():
            #print(gmi)
            if ii in i:
                Gmi[j] = gmi
                j += 1

        obj.gmi = Gmi
        obj.n = n
        self.i = n



    def allocate(self, ncards):
        self.n = ncards
        #float_fmt = self.model.float_fmt

        #: Element identification number
        self.element_id = zeros(ncards, 'int32')

        #: Identification number of grid point to which all six independent
        #: degrees-of-freedom for the element are assigned. (Integer > 0)
        self.gn = zeros(ncards, 'int32')

        #: Component numbers of the dependent degrees-of-freedom in the
        #: global coordinate system at grid points GMi. (Integers 1 through
        #: 6 with no embedded blanks.)
        self.cm = zeros(ncards, '|U6')

        #: Grid point identification numbers at which dependent
        #: degrees-of-freedom are assigned. (Integer > 0)
        self.alpha = zeros(ncards, 'int32')
        self.gmi = {}


    def add_card(self, card: BDFCard, comment: str=''):
        #self.model.log.debug('RBE2.add')
        i = self.i
        #if comment:
            # self.comment = comment
        eid = integer(card, 1, 'element_id')
        gn = integer(card, 2, 'gn')
        cm = components_or_blank(card, 3, 'cm')
        #assert gn is not None, 'gn=%s' % self.gn
        #assert cm is not None, 'cm=%s' % self.cm

        self.element_id[i] = eid
        self.gn[i] = gn
        self.cm[i] = cm

        alpha = integer_or_double(card, len(card) - 1, 'alpha')
        if isinstance(alpha, float):
            # the last field is not part of Gmi
            self.alpha[i] = alpha
            n = 1
        else:
            # the last field is part of Gmi
            n = 0
            #self.alpha[i] = 0.0  # we don't need to set alpha

        j = 4
        gmis = []
        for k in range(len(card) - 4 - n):
            gmi = integer(card, j + k, 'Gm%i' % (k + 1))
            #print('GM%i = %s' % (k + 1, gmi))
            gmis.append(gmi)
        self.gmi[i] = gmis
        self.i += 1

        #self.gn = str(self.gn)
        #self.cm = str(self.cm)

    def add_op2(self, data):
        self.eid = data[0]
        self.gn = data[1]
        self.cm = data[2]
        self.gmi = data[3]
        self.alpha = data[4]
        self.model.log.debug("eid=%s gn=%s cm=%s Gmi=%s alpha=%s"
                             % (self.eid, self.gn, self.cm, self.gmi, self.alpha))
        raise NotImplementedError('RBE2 data...')

        assert self.gn is not None, 'gn=%s' % self.gn
        assert self.cm is not None, 'cm=%s' % self.cm
        #self.gn = str(self.gn)
        #self.cm = str(self.cm)

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.gn = self.gn[i]
            self.gn = self.gn[i]
            self.cm = self.cm[i]
            self.alpha = self.alpha[i]
            unique_eids = unique(self.element_id)

            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate RBE2 IDs...')
        else:
            self.element_id = array([], dtype='int32')

    #def convert_to_MPC():
        #pass

    #def raw_fields(self):
        #list_fields = ['RBE2', self.eid, self.gn, self.cm] + self.Gmi + [self.alpha]
        #return list_fields

    #def repr_fields(self):
        #alpha = set_blank_if_default(self.alpha, 0.)
        #list_fields = ['RBE2', self.eid, self.gn, self.cm] + self.Gmi + [alpha]
        #return list_fields

    def write_card(self, bdf_file, size, is_double):
        if self.n == 0:
            return
        if size == 8:
            for j, eid, gn, cm, alpha in zip(count(), self.element_id, self.gn, self.cm, self.alpha):
                gmi = self.gmi[j]
                salpha = set_blank_if_default(alpha, 0.)
                list_fields = ['RBE2', eid, gn, cm] + gmi + [salpha]
                bdf_file.write(print_card_8(list_fields))
        elif is_double:
            for j, eid, gn, cm, alpha in zip(count(), self.element_id, self.gn, self.cm, self.alpha):
                gmi = self.gmi[j]
                salpha = set_blank_if_default(alpha, 0.)
                list_fields = ['RBE2', eid, gn, cm] + gmi + [salpha]
                bdf_file.write(print_card_16(list_fields))
        else:
            for j, eid, gn, cm, alpha in zip(count(), self.element_id, self.gn, self.cm, self.alpha):
                gmi = self.gmi[j]
                salpha = set_blank_if_default(alpha, 0.)
                list_fields = ['RBE2', eid, gn, cm] + gmi + [salpha]
                bdf_file.write(print_card_double(list_fields))


        #card = self.repr_fields()
        #return self.comment + print_card_8(card)

    def __repr__(self):
        return '<%s object; n=%s>' % (self.type, self.n)
