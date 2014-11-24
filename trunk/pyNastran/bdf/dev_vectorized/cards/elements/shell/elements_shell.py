from numpy import array, zeros, searchsorted, unique, concatenate, argsort, hstack

from .ctria3 import CTRIA3
from .ctria6 import CTRIA6
#from .cquad import CQUAD
#from .cquadx import CQUADX
from .cquad4 import CQUAD4
from .cquad8 import CQUAD8
#from .cquad9 import CQUAD9

#from .ctriax import CTRIAX
from .ctriax6 import CTRIAX6


class ElementsShell(object):
    def __init__(self, model):
        """
        Defines the ElementsShell object.

        :param self: the ElementsShell object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0

        self.ctria3 = CTRIA3(self.model)
        self.ctria6 = CTRIA6(self.model)
        #self.cquad = CQUAD(self.model)
        #self.cquadx = CQUADX(self.model)
        self.cquad4 = CQUAD4(self.model)
        self.cquad8 = CQUAD8(self.model)
        #self.cquad9 = CQUAD9(self.model)

        #self.ctriax = CTRIAX(self.model)
        self.ctriax6 = CTRIAX6(self.model)

    def allocate(self, card_count):
        etypes = self._get_types(nlimit=False)
        for etype in etypes:
            if etype.type in card_count:
                self.model.log.debug('    allocate %s' % etype.type)
                etype.allocate(card_count[etype.type])
            #else:
                #assert hasattr(etype, 'allocate'), '%s doesnt support allocate' % etype.type

    def build(self):
        #print('*build elements_shell')
        types = self._get_types(nlimit=False)
        self.n = 0
        for elems in types:
            elems.build()
            self.n += elems.n
        #print(' build elements_shell n=%s' % self.n)
        #eid = hstack(ctria3.element_id, cquad4.element_id)
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
        #self.cquad.add(card, comment)
        pass

    def add_cquadx(self, card, comment):
        #self.cquadx.add(card, comment)
        pass

    def add_cquad4(self, card, comment):
        self.cquad4.add(card, comment)

    def add_cquad8(self, card, comment):
        self.cquad8.add(card, comment)

    def add_cquad9(self, card, comment):
        #self.cquad9.add(card, comment)
        pass

    def add_ctriax(self, card, comment):
        self.ctriax.add(card, comment)

    def add_ctriax6(self, card, comment):
        self.ctriax6.add(card, comment)

    #=========================================================================
    # helper methods
    def _get_element_property_ids(self, element_ids):
        # figure out the element ids
        if element_ids is None:
            element_ids = []
            property_ids = []
            types = self._get_types()
            for element in types:
                if element.n:
                    element_ids.extend(element.element_id)
                    property_ids.extend(element.property_id)
            element_ids = array(element_ids)
            property_ids = array(property_ids)
        else:
            # figure out the property ids
            #raise NotImplementedError()
            element_ids = []
            property_ids = []
            types = self._get_types()
            for element in types:
                if element.n:
                    element_ids.extend(element.element_id)
                    property_ids.extend(element.property_id)
            element_ids2 = array(element_ids)
            property_ids2 = array(property_ids)
            i = searchsorted(element_ids, element_ids2)
            property_ids = property_ids2[i]

        return element_ids, property_ids

    def _unique_property_ids(self, element_ids=None):
        types = self._get_types()
        pids = []
        for element in types:
            if element.n:
                pids.extend(element.property_id)
        return unique(pids)
    #=========================================================================

    def get_area_by_element_id(self, element_id=None):
        pass

    def get_thickness_by_element_id(self, element_id=None):
        element_id, property_id = self._get_element_property_ids(element_id)
        assert len(element_id) == 1, element_id

        # lessen the work
        unique_pids = unique(property_id)
        if unique_pids is None:
            raise RuntimeError('there are no properties...')

        # makes searchsorted work
        unique_pids.sort()

        # now that we have the reduced seget the p
        pids = self._unique_property_ids(element_id)
        isort = argsort(pids)
        unique_thickness = self.model.properties_shell.get_thickness_by_property_id(unique_pids)

        n = len(element_id)
        thickness = zeros(n, 'float64')
        for i, pid in enumerate(property_id):
            #print("unique=%s pid=%s" % (unique_pids, pid))
            j = searchsorted(unique_pids, pid)
            thickness[i] = unique_thickness[j]
        return thickness

    def get_mass_by_element_id(self, element_id=None, total=False):
        types = self._get_types(nlimit=True)
        if element_id is None:
            element_id = []
            for etype in types:
                element_id.extend(etype.element_id)

        n = len(element_id)
        #print('element_ids =', element_id)
        #print('property_ids =', property_id)
        massi = zeros(n, dtype='float64')

        etypes = [etype.type for type in types]
        #print("etypes =", etypes)
        massi = zeros(n, dtype='float64')

        n0 = 0
        for elems in types:
            #if elems.n > 0:
            try:
                mass = elems.get_mass_by_element_id()
                massi[n0:n0 + elems.n] = elems.get_mass_by_element_id()
            except TypeError:  # ..todo:: remove this
                 massi[n0] = 1.
            except ValueError:  # ..todo:: remove this
                 massi[n0] = 1.
            n0 += elems.n
        assert massi.sum() > 0, elems.type
        #print("massii =", massi)
        if total:
            mass = massi.sum()
        else:
            mass = massi
        return mass

    def write_bdf(self, f, size=8, element_ids=None):
        f.write('$ELEMENTS_SHELL\n')
        types = self._get_types(nlimit=True)
        for element in types:
            if element.n:
                #print(element.type)
                element.write_bdf(f, size=size, element_ids=element_ids)

    def _get_types(self, nlimit=True):
        types = [self.ctria3, self.cquad4, self.ctria6, self.cquad8, self.ctriax6] #, cquad8
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
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
