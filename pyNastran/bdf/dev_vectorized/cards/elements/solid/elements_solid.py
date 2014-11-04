from .ctetra4 import CTETRA4
from .cpenta6 import CPENTA6
from .chexa8 import CHEXA8

from .ctetra10 import CTETRA10
from .cpenta15 import CPENTA15
from .chexa20 import CHEXA20
from numpy import hstack

class ElementsSolid(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ElementsSolid object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0

        self.ctetra4 = CTETRA4(self.model)
        self.cpenta6 = CPENTA6(self.model)
        self.chexa8  = CHEXA8(self.model)

        self.ctetra10 = CTETRA10(self.model)
        self.cpenta15 = CPENTA15(self.model)
        self.chexa20  = CHEXA20(self.model)

    def allocate(self, card_count):
        etypes = self._get_types(nlimit=False)
        for etype in etypes:
            if etype.type in card_count:
                etype.allocate(card_count[etype.type])
            #else:
                #assert hasattr(ptype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        self.n = 0
        types = self._get_types(nlimit=False)
        for elems in types:
            if elems.n:
                self.model.log.debug('    building ES - %s' % elems.__class__.__name__)
            elems.build()
            self.n += elems.n

        self.element_id = hstack([
            self.ctetra4.element_id,
            self.cpenta6.element_id,
            self.chexa8.element_id,
            self.ctetra10.element_id,
            self.cpenta15.element_id,
            self.chexa20.element_id,
        ])
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

    #=========================================================================
    def get_mass(self, element_ids=None, total=False):
        types = self._get_types()
        massi = []
        for elems in types:
            if elems.n > 0:
                massi = elems.get_mass()

        if total:
            mass = massi.sum()
        else:
            mass = massi
        return mass

    #=========================================================================
    def write_bdf(self, f, size=8, element_ids=None):
        f.write('$ELEMENTS_SOLID\n')
        types = self._get_types(nlimit=True)
        for elems in types:
            #print(elems.type)
            elems.write_bdf(f, size=size, element_ids=element_ids)

    def _get_types(self, nlimit=True):
        types = [self.ctetra4, self.cpenta6, self.chexa8,
                 self.ctetra10, self.cpenta15, self.chexa20,
                 ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    #print("etype.Type =", etype.Type)
                    types2.append(etype)
            types = types2
        #print("solid nlimit=%s" % nlimit)
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