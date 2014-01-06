from numpy import zeros, array, union1d, searchsorted, concatenate
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
        self.n = 0
        #self.pmass = PMASS(model)
        #self.cmass1 = CMASS1(model)
        #self.cmass2 = CMASS2(model)
        #self.cmass3 = CMASS3(model)
        #self.cmass4 = CMASS4(model)
        #self.cmass5 = CMASS5(model)

    def build(self):
        #self.n = 0
        types = self._get_types(nlimit=False)
        #print('bt', types)
        for elems in types:
            elems.build()
            self.n += elems.n
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

    #=========================================================================
    def get_indexs(self, element_ids=None):
        #mass_types = self._get_types()
        #for mass_type in mass_types:
            #element_ids.extend(mass_type.element_id)

        types = self._get_types()
        if types:
            _element_ids = concatenate([mtype.element_id for mtype in types])
            i = argsort(_element_ids)
            return _element_ids, i
        return None, None
    
    def get_mass(self, element_ids=None, total=False):
        assert element_ids is None
        mass_types = self._get_types()
        
        element_ids, i = self.get_indexs(element_ids)
        if element_ids is None:
            return 0.0

        n = len(element_ids)

        for mass_type in mass_types:
            element_ids2 = union1d(element_ids, mass_type.element_id)
            massi = mass_type.get_mass(element_ids2, total)

        if total:
            mass = massi.sum()
        else:
            mass = massi

        return mass

    #=========================================================================
    def write_bdf(self, f, size=8, element_ids=None):
        f.write('$ELEMENTS_MASS\n')
        types = self._get_types()
        
        print("types", types)
        for elems in types:
            print "MASS", elems.type
            elems.write_bdf(f, size=size, element_ids=element_ids)

    def _get_types(self, nlimit=True):
        mtypes = [self.conm1, self.conm2,
                # self.cmass1, self.cmass2, self.cmass3, self.cmass4, self.cmass5,
                # self.pmass,
        ]
        if nlimit:
            d = []
            for mtype in mtypes:
                    print('type=%s n=%s' % (mtype.type, mtype.n))
                #if mtype.n > 0:
                    d.append(mtype)
            #mtypes = d
            return d
        else:
            return mtypes
            #return [mtype if mtype.n > 0 for mtype in mtypes]
        #return mtypes

    def get_stats(self):
        msg = []
        types = self._get_types(nlimit=True)
        for element in types:
            nele = element.n
            if nele:
                msg.append('  %-8s: %i' % (element.type, nele))
        return msg

    def _verify(self, xref=True):
        types = self._get_types()
        for elems in types:
            elems._verify(xref=xref)

    def __repr__(self):
        return 'hi' + '\n'.join(self.get_stats())