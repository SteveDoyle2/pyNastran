from numpy import concatenate, argsort, searchsorted, ndarray

from .pbeam import PBEAM
from .pbeaml import PBEAML

class PropertiesBeam(object):
    def __init__(self, model):
        """
        Defines the PropertiesBar object.

        :param self: the PropertiesBar object
        :param model: the BDF object
        """
        self.model = model

        self.pbeam = PBEAM(self.model)
        self.pbeaml = PBEAML(self.model)
        self.n = 0

    def allocate(self, card_count):
        ptypes = self._get_types(nlimit=False)
        for ptype in ptypes:
            if ptype.type in card_count:
                ptype.allocate(card_count)
                del card_count[ptype.type]
            #else:
                #assert hasattr(etype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        #print('building beam')
        self.pbeam.build()
        self.pbeaml.build()

        npbeam = self.pbeam.n
        npbeaml = self.pbeaml.n
        self.n = npbeam + npbeaml

        #if npshell and npcomp and npcompg:
            #asdf
        #if npshell and npcomp:
            #pid = concatenate(self.pshell.property_id, self.pcomp.property_id)
            #unique_pids = unique(pid)
            #print unique_pids
            #if len(unique_pids) != len(pid):
                #raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        #else:
            #pass

    def rebuild(self):
        raise NotImplementedError()

    #=========================================================================
    def add_pbeam(self, card, comment):
        self.pbeam.add(card, comment)

    def add_pbeaml(self, card, comment):
        self.pbeaml.add(card, comment)

    #=========================================================================
    def get_area_by_element_id(self):
        return 0.0
    def get_mass_by_element_id(self):
        return 0.0

    #=========================================================================
    def _get_types(self, nlimit=True):
        types = [self.pbeam, self.pbeaml]
        if nlimit:
            types2 = []
            for ptype in types:
                if ptype.n:
                    types2.append(ptype)
            types = types2
        return types

    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_bdf(self, f, size=8, property_ids=None):
        f.write('$PROPERTIES_BEAM\n')
        types = self._get_types(nlimit=False)
        for prop in types:
            prop.write_bdf(f, size=size, property_ids=property_ids)