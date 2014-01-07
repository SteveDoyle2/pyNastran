from numpy import concatenate, argsort, searchsorted, ndarray

from .pshell import PSHELL
from .pcomp import PCOMP
from .pcompg import PCOMPG

class PropertiesShell(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        """
        self.model = model

        self.pshell = PSHELL(self.model)
        self.pcomp = PCOMP(self.model)
        self.pcompg = PCOMPG(self.model)

    def build(self):
        self.pshell.build()
        self.pcomp.build()
        self.pcompg.build()
        
        npshell = self.pshell.n
        npcomp  = self.pcomp.n
        npcompg = self.pcompg.n
        
        if npshell and npcomp and npcompg:
            asdf
        if npshell and npcomp:
            pid = concatenate(self.pshell.property_id, self.pcomp.property_id)
            unique_pids = unique(pid)
            print unique_pids
            if len(unique_pids) != len(pid):
                raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        else:
            pass

    def rebuild(self):
        raise NotImplementedError()

    #=========================================================================
    def add_pshell(self, card, comment):
        self.pshell.add(card, comment)

    def add_pcomp(self, card, comment):
        self.pcomp.add(card, comment)

    def add_pcompg(self, card, comment):
        raise NotImplementedError(card)
        self.pcompg.add(card, comment)

    #=========================================================================
    def get_thickness(self, property_ids=None):
        """
        Gets the thickness of the PSHELLs/PCOMPs.
        
        :param self: the ShellProperties object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        types = self._get_types(nlimit=True)
        _property_ids = concatenate([ptype.property_id for ptype in types])
        t = concatenate([ptype.get_thickness() for ptype in types] )
        if property_ids is None:
            return t
            #property_ids = _property_ids
        
        assert isinstance(property_ids, ndarray), type(property_ids)
        i = argsort(_property_ids)
        
        print(_property_ids[i])
        print(property_ids)
        j = searchsorted(property_ids, _property_ids[i])
        t2 = t[j]
        return t2

    #=========================================================================
    def _get_types(self, nlimit=True):
        types = [self.pshell, self.pcomp, self.pcompg]
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
        f.write('$PROPERTIES_SHELL\n')
        types = self._get_types()
        for prop in types:
            #print('*SHELL', prop.type)
            prop.write_bdf(f, size=size, property_ids=property_ids)