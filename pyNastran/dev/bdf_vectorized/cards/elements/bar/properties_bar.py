from pyNastran.dev.bdf_vectorized.cards.elements.bar.pbar import PBAR
from pyNastran.dev.bdf_vectorized.cards.elements.bar.pbarl import PBARL

class PropertiesBar:
    def __init__(self, model):
        """
        Defines the PropertiesBar object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model

        self.pbar = model.pbar
        self.pbarl = model.pbarl
        self.n = 0

    def allocate(self, card_count):
        ptypes = self._get_types(nlimit=False)
        for ptype in ptypes:
            if ptype.type in card_count:
                if hasattr(ptype, 'n'):
                    self.model.log.debug('    building %s.n = %s' % (ptype.__class__.__name__, ptype.n))
                else:
                    raise RuntimeError(ptype)
                ptype.allocate(card_count[ptype.type])
                del card_count[ptype.type]
            else:
                assert hasattr(ptype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        self.pbar.build()
        self.pbarl.build()

        npbar = self.pbar.n
        npbarl = self.pbarl.n
        self.n = npbar + npbarl
        #if npshell and npcomp and npcompg:
            #raise NotImplementedError()

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
    def add_pbar(self, card, comment):
        self.pbar.add(card, comment)

    def add_pbarl(self, card, comment):
        self.pbarl.add(card, comment)

    #=========================================================================
    def get_area_by_element_id(self):
        return 0.0
    def get_mass_by_element_id(self):
        return 0.0

    #=========================================================================
    def _get_types(self, nlimit=True):
        types = [self.pbar, self.pbarl]
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

    def write_card(self, bdf_file, size=8, property_id=None):
        bdf_file.write('$PROPERTIES_BAR\n')
        types = self._get_types(nlimit=False)
        for prop in types:
            #print('prop.type = %s' % prop.type)
            prop.write_card(bdf_file, size=size, property_id=property_id)

    def __repr__(self):
        msg = '<%s object; n=%s>\n' % (self.__class__.__name__, self.n)
        types = self._get_types()
        for prop in types:
            msg += '  <%s object; n=%s>\n' % (prop.type, prop.n)
        return msg
