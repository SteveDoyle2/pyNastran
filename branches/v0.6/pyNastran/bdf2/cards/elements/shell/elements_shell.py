from .ctria3 import CTRIA3
#from .ctria6 import CTRIA6
#from .cquad import CQUAD
#from .cquadx import CQUADX
from .cquad4 import CQUAD4
#from .cquad8 import CQUAD8
#from .cquad9 import CQUAD9


class ElementsShell(object):
    def __init__(self, model):
        """
        Defines the ElementsShell object.

        :param self: the ElementsShell object
        :param model: the BDF object
        """
        self.model = model

        self.ctria3 = CTRIA3(self.model)
        self.cquad4 = CQUAD4(self.model)
        #self.cquad8 = CQUAD8(self.model)

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

    def add_ctria3(self, card, comment):
        self.ctria3.add(card, comment)

    def add_ctria6(self, card, comment):
        self.ctria6.add(card, comment)

    def add_cquad(self, card, comment):
        self.cquad.add(card, comment)

    def add_cquad4(self, card, comment):
        self.cquad4.add(card, comment)

    def add_cquad8(self, card, comment):
        self.cquad8.add(card, comment)

    def add_cquad9(self, card, comment):
        self.cquad9.add(card, comment)

    def write_bdf(self, f, size=8, eids=None):
        f.write('$ELEMENTS\n')
        types = self._get_types()
        for element in types:
            element.write_bdf(f, size=size, eids=eids)

    def _get_types(self):
        types = [self.ctria3, self.cquad4] #, cquad8
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