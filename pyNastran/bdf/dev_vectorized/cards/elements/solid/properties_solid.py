from numpy import concatenate, argsort, searchsorted, ndarray

from .psolid import PSOLID
#from .plsolid import PLSOLID

class PropertiesSolid(object):
    def __init__(self, model):
        """
        Defines the PropertiesSolid object.

        :param self: the PropertiesSolid object
        :param model: the BDF object
        """
        self.model = model

        self.psolid = PSOLID(self.model)
        #self.plsolid = PLSOLID(self.model)
        self.n = 0

    def build(self):
        #print "building solid properties"
        types = self._get_types(nlimit=False)
        for prop in types:
            #print "prop =", prop
            prop.build()
            self.n += prop.n

        npsolid = self.psolid.n
        self.property_id = self.psolid.property_id
        #print('npsolid =', npsolid)
        #assert npsolid > 0
        #nplsolid = self.plsolid.n

        #if npshell and npcomp and pshear:
            #asf
        #elif npshell and npcomp:
            #pid = concatenate(self.pshell.pid, self.pcomp.pid)
            #unique_pids = unique(pid)
            #print unique_pids
            #if len(unique_pids) != len(pid):
                #raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        #else:
            #pass

    def rebuild(self):
        raise NotImplementedError()

    def add_psolid(self, card, comment):
        self.psolid.add(card, comment)

    def add_plsolid(self, card, comment):
        self.plsolid.add(card, comment)

    #=========================================================================
    def get_mid(self, property_ids):
        types = self._get_types(nlimit=True)
        _material_ids = concatenate([ptype.material_id for ptype in types])
        #print _material_ids
        return _material_ids

    #def get_mid(self, property_ids):
        #types = self._get_types(nlimit=True)
        #_material_ids = concatenate([ptype.material_id for ptype in types])
        #print _property_ids
        #return _property_ids

    def _get_types(self, nlimit=True):
        types = [self.psolid]
        if nlimit:
            types2 = []
            for ptype in types:
                if ptype.n:
                    types2.append(ptype)
            types = types2
        return types

    #=========================================================================
    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_bdf(self, f, size=8, property_ids=None):
        f.write('$PROPERTIES_SOLID\n')
        types = self._get_types()
        for prop in types:
            prop.write_bdf(f, size=size, property_ids=property_ids)