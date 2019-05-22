from itertools import count
from numpy import array, zeros, unique, searchsorted, asarray, int64, where

from pyNastran.bdf.bdf_interface.assign_type import (double_or_blank,
    components, components_or_blank, integer_double_or_blank, blank, integer)

from pyNastran.bdf.field_writer_8 import (print_card_8, set_blank_if_default,
                                          #set_string8_blank_if_default
                                          )
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
#from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default


#RigidElement
class RBE3:
    type = 'RBE3'

    def __init__(self, model):
        """
        eid
        refgrid
        refc
        WtCG_groups = [wt,ci,Gij]
        Gmi
        Cmi
        alpha
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

        obj = RBE3(self.model)
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

        #: Grid point identification numbers at which dependent
        #: degrees-of-freedom are assigned. (Integer > 0)
        self.alpha = zeros(ncards, 'int32')
        self.gmi = {}


    def add_card(self, card, comment=''):
        #self.model.log.debug('RBE2.add')
        i = self.i
        #if comment:
            # self.comment = comment
        eid = integer(card, 1, 'element_id')
        #if comment:
            # self.comment = comment
        self.element_id[i] = integer(card, 1, 'eid')
        blank(card, 2, 'blank')
        self.refgrid[i] = integer(card, 3, 'refgrid')
        self.refc[i] = components_or_blank(card, 4, 'refc', 0)
        #iUM = fields.index('UM')

        fields = [field.upper() if isinstance(field, str) else field for field in card[5:]]
        iOffset = 5
        iWtMax = len(fields) + iOffset
        try:
            iAlpha = fields.index('ALPHA') + iOffset
            iWtMax = iAlpha  # the index to start parsing UM
            iUmStop = iAlpha  # the index to stop  parsing UM
        except ValueError:
            iAlpha = None
            iUmStop = iWtMax
        #print("iAlpha = %s" % iAlpha)
        try:
            iUm = fields.index('UM') + iOffset
            iWtMax = iUm
        except ValueError:
            iUm = None
        #print("iAlpha=%s iUm=%s" % (iAlpha, iUm))
        #print("iAlpha=%s iWtMax=%s" % (iAlpha, iWtMax))

        #print("iUM = %s" % iUM)
        WtCG_groups = []

        i = iOffset
        n = 1
        while i < iWtMax:
            Gij = []
            wtname = 'wt' + str(n)
            wt = double_or_blank(card, i, wtname)
            if wt is not None:
                cname = 'c'+str(n)
                ci = components_or_blank(card, i + 1, cname)

                #print("%s=%s %s=%s" % (wtname, wt, cname, ci))
                i += 2
                gij = 0

                j = 0
                while isinstance(gij, int) and i < iWtMax:
                    j += 1
                    gij_name = 'g%s,%s' % (n, j)
                    gij = integer_double_or_blank(card, i, gij_name)
                    if isinstance(gij, float):
                        break
                    #print("%s = %s" % (gij_name, gij))
                    if gij is not None:
                        Gij.append(gij)
                    i += 1
                wtCG_group = [wt, ci, Gij]
                WtCG_groups.append(wtCG_group)
                #print('----finished a group=%r----' % wtCG_group)
            else:
                i += 1
        self.WtCG_groups[i] = WtCG_groups
        Gmi = []
        Cmi = []
        #print("")
        if iUm:
            #print('UM = %s' % card.field(iUm))  # UM
            i = iUm + 1
            n = 1
            #print("i=%s iUmStop=%s" % (i, iUmStop))
            for j in range(i, iUmStop, 2):

                gm_name = 'gm' + str(n)
                cm_name = 'cm' + str(n)
                gmi = integer_or_blank(card, j, gm_name)
                if gmi is not None:
                    cmi = components(card, j + 1, cm_name)
                    #print("gmi=%s cmi=%s" % (gmi, cmi))
                    Gmi.append(gmi)
                    Cmi.append(cmi)
        self.Gmi[i] = Gmi
        self.Cmi[i] = Cmi
        if iAlpha:
            alpha = double_or_blank(card, iAlpha + 1, 'alpha', 0.0)
        else:
            alpha = 0.0
        self.alpha[i] = alpha

    #def add_op2(self, data):
        #self.eid = data[0]
        #self.gn = data[1]
        #self.cm = data[2]
        #self.Gmi = data[3]
        #self.alpha = data[4]
        #print("eid=%s gn=%s cm=%s Gmi=%s alpha=%s"
              #% (self.eid, self.gn, self.cm, self.Gmi, self.alpha))
        #raise NotImplementedError('RBE2 data...')

        #assert self.gn is not None, 'gn=%s' % self.gn
        #assert self.cm is not None, 'cm=%s' % self.cm
        #self.gn = str(self.gn)
        #self.cm = str(self.cm)

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            #self.gn = self.gn[i]
            #self.gn = self.gn[i]
            #self.cm = self.cm[i]
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

    def write_card(self, bdf_file, size=True, is_double=False):
        if self.n:
            if size == 8:
                for j, eid, alpha in zip(count(), self.element_id, self.alpha):
                    salpha = set_blank_if_default(alpha, 0.)
                    list_fields = ['RBE2', eid] + [salpha]
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
