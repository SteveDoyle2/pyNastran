from .ctetra import CTETRA4
from .cpenta6 import CPENTA6
from .chexa import CHEXA8


class ElementsSolid(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ElementsSolid object
        :param model: the BDF object
        """
        self.model = model

        self.ctetra4 = CTETRA4(self.model)
        self.cpenta6 = CPENTA6(self.model)
        self.chexa8  = CHEXA8(self.model)

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

    def add_ctetra4(self, card, comment):
        self.ctetra4.add(card, comment)

    def add_cpenta6(self, card, comment):
        self.cpenta6.add(card, comment)

    def add_chexa8(self, card, comment):
        self.chexa8.add(card, comment)

    #===========
    def add_ctetra10(self, card, comment):
        self.ctetra10.add(card, comment)

    def add_cpenta15(self, card, comment):
        self.cpenta15.add(card, comment)

    def add_chexa20(self, card, comment):
        self.chexa20.add(card, comment)

    #===========
    def write_bdf(self, f, size=8, eids=None):
        f.write('$ELEMENTS_SOLID\n')
        aaa
        self.ctetra4.write_bdf(f, size=size, eids=eids)
        #self.cpenta6.write_bdf(f, size=size, eids=eids)
        self.chexa8.write_bdf(f, size=size, eids=eids)

    def _get_types(self):
        types = [self.ctetra4, self.cpenta6, self.chexa8]
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