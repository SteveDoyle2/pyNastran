from .conm1 import CONM1
from .conm2 import CONM2

#from .pmass import PMASS
#from .cmass1 import CMASS1
#from .cmass2 import CMASS2
#from .cmass3 import CMASS3
#from .cmass4 import CMASS4
#from .cmass5 import CMASS5

class Mass(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ElementsSolid object
        :param model: the BDF object
        """
        self.model = model
        self.conm1 = CONM1(model)
        self.conm2 = CONM2(model)
        #self.pmass = PMASS(model)
        #self.cmass1 = CMASS1(model)
        #self.cmass2 = CMASS2(model)
        #self.cmass3 = CMASS3(model)
        #self.cmass4 = CMASS4(model)
        #self.cmass5 = CMASS5(model)

    def build(self):
        types = self._get_types()
        for elems in types:
            elems.build()

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
        #    raise RuntimeError('There are duplicate CTRIA3/CQUAD4 IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_conm1(self, card, comment):
        self.conm1.add(card, comment)

    def add_conm2(self, card, comment):
        self.conm2.add(card, comment)

    def add_cmass1(self, card, comment):
        self.cmass1.add(card, comment)

    def add_cmass2(self, card, comment):
        self.cmass2.add(card, comment)

    def add_cmass3(self, card, comment):
        self.cmass3.add(card, comment)

    def add_cmass4(self, card, comment):
        self.cmass4.add(card, comment)

    def add_cmass5(self, card, comment):
        self.cmass5.add(card, comment)

    def add_pmass(self, card, comment):
        self.pmass.add(card, comment)

    #===========
    def write_bdf(self, f, size=8, element_ids=None):
        f.write('$ELEMENTS_MASS\n')
        types = self._get_types()
        for elems in types:
            print "MASS", elems.type
            elems.write_bdf(f, size=size, element_ids=element_ids)

    def _get_types(self):
        types = [self.conm1, self.conm2,
                # self.cmass1, self.cmass2, self.cmass3, self.cmass4, self.cmass5,
                # self.pmass,
                 ]
        return types

    def get_stats(self):
        msg = []
        types = self._get_types()
        for element in types:
            nele = element.n
            if nele:
                msg.append('  %-8s: %i' % (element.type, nele))
        return msg

    def _verify(self, xref=True):
        types = self._get_types()
        for elems in types:
            elems._verify(xref=xref)