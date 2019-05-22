from .paero1 import PAERO1
#from .paero2 import PAERO2
#from .paero3 import PAERO3
#from .paero4 import PAERO4
#from .paero5 import PAERO5

class PAero:
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model

        self.paero1 = PAERO1(self.model)
        #self.paero2 = PAERO2(self.model)
        #self.paero3 = PAERO3(self.model)
        #self.paero4 = PAERO4(self.model)
        #self.paero5 = PAERO5(self.model)

    def build(self):
        types = self._get_types()
        for elems in types:
            elems.build()

        #eid = concatenate(pshell.pid, pcomp.pid)
        #unique_eids = unique(eid)
        #if unique_eids != len(eid):
        #    raise RuntimeError('There are duplicate PAERO1/PAERO2/etc. IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_paero1(self, card, comment):
        self.paero1.add(card, comment)

    def add_paero2(self, card, comment):
        self.paero2.add(card, comment)

    def add_paero3(self, card, comment):
        self.paero3.add(card, comment)

    def add_paero4(self, card, comment):
        self.paero4.add(card, comment)

    def add_paero5(self, card, comment):
        self.paero5.add(card, comment)

    #===========
    def write_card(self, bdf_file, size=8, element_ids=None):
        bdf_file.write('$AERO\n')
        types = self._get_types()
        for elems in types:
            #print("AERO", elems.type)
            elems.write_card(bdf_file, size=size, element_ids=element_ids)

    def _get_types(self):
        types = [
            self.paero1,
            #self.paero2, self.paero3, self.paero4, self.paero5,
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
