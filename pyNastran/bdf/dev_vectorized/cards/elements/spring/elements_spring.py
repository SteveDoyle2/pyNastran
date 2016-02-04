from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas1 import CELAS1
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas2 import CELAS2

from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas3 import CELAS3
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas4 import CELAS4

class ElementsSpring(object):
    def __init__(self, model):
        """
        Defines the ElementsShell object.

        :param model: the BDF object
        """
        self.model = model

        self.n = 0
        self.celas1 = CELAS1(self.model)
        self.celas2 = CELAS2(self.model)
        self.celas3 = CELAS3(self.model)
        self.celas4 = CELAS4(self.model)

    def allocate(self, card_count):
        etypes = self._get_types(nlimit=False)
        for etype in etypes:
            if etype.type in card_count:
                print('etype.type =', etype.type, etype.n)
                etype.allocate(card_count[etype.type])
            #else:
                #assert hasattr(ptype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        types = self._get_types(nlimit=False)
        for elems in types:
            elems.build()
            self.n += elems.n

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
            #raise RuntimeError('There are duplicate CTRIA3/CQUAD4 IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_celas1(self, card, comment):
        self.celas1.add(card, comment)

    def add_celas2(self, card, comment):
        self.celas2.add(card, comment)

    def add_celas3(self, card, comment):
        self.celas3.add(card, comment)
        raise NotImplementedError()

    def add_celas4(self, card, comment):
        self.celas4.add(card, comment)
        raise NotImplementedError()

    def write_card(self, f, size=8, eids=None):
        f.write('$ELEMENTS\n')
        types = self._get_types()
        for element in types:
            element.write_card(f, size=size, eids=eids)

    def _get_types(self, nlimit=True):
        types = [self.celas1,
                 self.celas2,
                 self.celas3,
                 self.celas4,
                 ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    #print("etype.Type =", etype.Type)
                    types2.append(etype)
            types = types2
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
