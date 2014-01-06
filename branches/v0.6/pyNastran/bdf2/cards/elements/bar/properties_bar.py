from numpy import concatenate, argsort, searchsorted, ndarray

from .pbar import PBAR
from .pbarl import PBARL

class PropertiesBar(object):
    def __init__(self, model):
        """
        Defines the PropertiesBar object.

        :param self: the PropertiesBar object
        :param model: the BDF object
        """
        self.model = model

        self.pbar = PBAR(self.model)
        self.pbarl = PBARL(self.model)

    def build(self):
        self.pbar.build()
        self.pbarl.build()
        
        npbar = self.pbar.n
        npbarl = self.pbarl.n
        
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
    def add_pbar(self, card, comment):
        self.pbar.add(card, comment)

    def add_pbarl(self, card, comment):
        self.pbarl.add(card, comment)

    #=========================================================================
    def get_area(self):
        return 0.0
    def get_mass(self):
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

    def write_bdf(self, f, size=8, property_ids=None):
        f.write('$PROPERTIES_BAR\n')
        types = self._get_types()
        for prop in types:
            #print('*SHELL', prop.type)
            prop.write_bdf(f, size=size, property_ids=property_ids)