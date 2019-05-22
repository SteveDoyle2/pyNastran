from pyNastran.dev.bdf_vectorized.cards.aero.caero1 import CAERO1
#from .caero2 import CAERO2
#from .caero3 import CAERO3
#from .caero4 import CAERO4

class CAero:
    def __init__(self, model):
        """
        Defines the CAero object.

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.model = model

        self.caero1 = CAERO1(self.model)
        #self.caero2 = CAERO2(self.model)
        #self.caero3 = CAERO3(self.model)
        #self.caero4 = CAERO4(self.model)

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

    def add_caero1(self, card, comment):
        self.caero1.add(card, comment)

    def add_caero2(self, card, comment):
        self.caero2.add(card, comment)

    def add_caero3(self, card, comment):
        self.caero3.add(card, comment)

    def add_caero4(self, card, comment):
        self.caero4.add(card, comment)

    def add_caero5(self, card, comment):
        self.caero5.add(card, comment)

    #===========
    def write_card(self, bdf_file, size=8, element_id=None):
        bdf_file.write('$AERO\n')
        types = self._get_types()
        for elems in types:
            #print("AERO", elems.type)
            elems.write_card(bdf_file, size=size, element_id=element_id)

    def _get_types(self):
        types = [self.caero1, # self.caero2, self.caero3, self.caero4, self.caero5,
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
