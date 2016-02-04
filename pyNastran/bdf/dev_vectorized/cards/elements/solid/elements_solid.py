from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra4 import CTETRA4
from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta6 import CPENTA6
from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa8 import CHEXA8

from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra10 import CTETRA10
from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta15 import CPENTA15
from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa20 import CHEXA20
from numpy import hstack, zeros, where, searchsorted, argsort, full, nan, unique

class ElementsSolid(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param model: the BDF object
        """
        self.model = model
        self.n = 0

        self.ctetra4 = CTETRA4(self.model)
        self.cpenta6 = CPENTA6(self.model)
        self.chexa8 = CHEXA8(self.model)

        self.ctetra10 = CTETRA10(self.model)
        self.cpenta15 = CPENTA15(self.model)
        self.chexa20 = CHEXA20(self.model)

    def allocate(self, card_count):
        etypes = self._get_types(nlimit=False)
        for etype in etypes:
            if etype.type in card_count:
                self.model.log.debug('    allocate %s' % etype.type)
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
    def get_element_id_property_id(self):
        etypes = self._get_types()
        element_property = zeros((self.n, 2), dtype='int32')
        n0 = 0
        for etype in etypes:
            n = etype.n
            element_property[n0:n0+n, 0] = etype.element_id
            element_property[n0:n0+n, 1] = etype.property_id
            n0 += n
        return element_property

    def get_property_id_material_id(self):
        ptypes = {
            'PSOLID' : self.model.elements.element_solids.psolid.property_id,
            'PLSOLID' : self.model.elements.element_solids.plsolid.property_id,
        }
        #element_property = zeros((self.n, 2), dtype='int32')
        #n0 = n
        #for etype, elements in etypes.iteritems():
        #    n = element_property.n
        #    element_property[n0:n0+n, 0] = element_property.element_id
        #    element_property[n0:n0+n, 1] = element_property.property_id
        #    n0 += n
        #return element_property
        raise NotImplementedError()
        return None

    def get_material_id_by_element_id(self, element_id=None):
        if element_id is None:
            n = self.n
            element_id = self.element_id
        else:
            n = len(element_id)

        material_id = full(n, nan, dtype='int32')

        self.model.log.debug('----get_material_id_by_element_id----')
        self.model.log.debug('element_id = %s' % element_id)
        etypes = self._get_types()
        eid_pid = self.get_element_id_property_id()
        i = argsort(eid_pid[:, 0])
        eid_pid = eid_pid[i, :]

        j = searchsorted(eid_pid[:, 0], element_id)
        pids = eid_pid[j, 1]
        upids = unique(pids)
        ptypes = [
            self.model.elements.properties_solid.psolid,
            self.model.elements.properties_solid.plsolid,
        ]
        nupid = len(upids)
        upid_mids = zeros((nupid, 2), dtype='int32')
        for ptype in ptypes:
            for pid in upids:
                if pid in ptype.property_id:
                    i = where(pid == ptype.property_id)[0]
                    mid = ptype.material_id[i][0]
                    upid_mids[i] = mid
                    j = where(pid == pids)[0]
                    material_id[j] = mid
        self.model.log.debug('material_id = %s' % material_id)
        return material_id

    def get_mass_by_element_id(self, element_id=None, total=False):
        types = self._get_types()
        massi = []
        for elems in types:
            if elems.n > 0:
                massi = elems.get_mass_by_element_id()

        if total:
            mass = massi.sum()
        else:
            mass = massi
        return mass

    #=========================================================================
    def write_card(self, f, size=8, element_id=None):
        f.write('$ELEMENTS_SOLID\n')
        types = self._get_types(nlimit=True)
        for elems in types:
            self.model.log.debug(elems.type)
            elems.write_card(f, size=size, element_id=element_id)

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
