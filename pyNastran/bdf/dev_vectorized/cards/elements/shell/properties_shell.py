from numpy import concatenate, hstack, argsort, searchsorted, ndarray, unique, array

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
        self.n = 0

    def allocate(self, card_count):
        ptypes = self._get_types(nlimit=False)
        for ptype in ptypes:
            if ptype.type in card_count:
                ptype.allocate(card_count[ptype.type])
                del card_count[ptype.type]
            #else:
                #assert hasattr(etype, 'allocate'), '%s doesnt support allocate' % ptype.type

    def build(self):
        self.pshell.build()
        self.pcomp.build()
        self.pcompg.build()

        npshell = self.pshell.n
        npcomp  = self.pcomp.n
        npcompg = self.pcompg.n

        self.n = npshell + npcomp + npcompg
        pid = hstack([self.pshell.property_id, self.pcomp.property_id, self.pcompg.property_id])
        unique_pids = unique(pid)
        #print unique_pids
        if len(unique_pids) != len(pid):
            raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')

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
    def get_mid(self, property_ids):
        types = self._get_types(nlimit=True)
        _property_ids = concatenate([ptype.property_id for ptype in types])
        #print _property_ids
        return _property_ids

    def get_nonstructural_mass(self, property_ids=None):
        return self._getmethod(property_ids, 'get_nonstructural_mass')

    def get_mass_per_area(self, property_ids=None):
        return self._getmethod(property_ids, 'get_mass_per_area')

    def get_thickness(self, property_ids=None):
        """
        Gets the thickness of the PSHELLs/PCOMPs.

        :param self: the ShellProperties object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        return self._getmethod(property_ids, 'get_thickness')

    def get_density(self, property_ids=None):
        """
        Gets the density of the PSHELLs/PCOMPs.

        :param self: the ShellProperties object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        return self._getmethod(property_ids, 'get_density')

    def _getattr(self, property_ids, method):
        TypeMap = {
            'PSHELL': (self.pshell, self.pshell.property_id),
            'PCOMP' : (self.pcomp, self.pcomp.property_id),
            'PCOMPG': (self.pcompg, self.pcompg.property_id),
        }
        out = array([])
        for Type, (ptype, pids) in sorted(TypeMap.iteritems()):
            for pid in pids:
                if pid in property_ids:
                    value = getattr(ptype, method)
                    #assert len(value) == 1, value
                    #print('pid =', pid, value)
                    out = hstack([out, value])
        #print("out = %s" % out)
        return out

    def _getmethod(self, property_ids, method):
        #types = self._get_types(nlimit=True)
        #data = hstack([getattr(ptype, method)() for ptype in types] )
        #if property_ids is None:
            #return data

        TypeMap = {
            'PSHELL': (self.pshell, self.pshell.property_id),
            'PCOMP' : (self.pcomp, self.pcomp.property_id),
            'PCOMPG': (self.pcompg, self.pcompg.property_id),
        }
        #print('property_ids = %s' % property_ids)
        out = array([])
        for Type, (ptype, pids) in sorted(TypeMap.iteritems()):
            for pid in pids:
                if pid in property_ids:
                    value = getattr(ptype, method)([pid])
                    #assert len(value) == 1, value
                    #print('pid = %s %s' % (pid, value))
                    out = hstack([out, value])
        #print("out = %s" % out)
        return out
        #_property_ids = hstack([ptype.property_id for ptype in types])
        #assert isinstance(property_ids, ndarray), type(property_ids)
        #i = argsort(_property_ids)
        #j = searchsorted(property_ids, _property_ids[i])
        #data_short = data[j]
        #return data_short

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