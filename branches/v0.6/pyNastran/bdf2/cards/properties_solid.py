from .psolid import PSOLID

class PropertiesSolid(object):
    def __init__(self, model):
        """
        Defines the PropertiesSolid object.

        :param self: the PropertiesSolid object
        :param model: the BDF object
        """
        self.model = model

        self.psolid = PSOLID(self.model)

    def build(self):
        self.psolid.build()
        
        npsolid = self.psolid.n
        #npshell = self.pshell.n
        #npshear = self.pshear.n
        
        #if npshell and npcomp and pshear:
        #    asf
        #elif npshell and npcomp:
        #    pid = concatenate(self.pshell.pid, self.pcomp.pid)
        #    unique_pids = unique(pid)
        #    print unique_pids
        #    if len(unique_pids) != len(pid):
        #        raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        #else:
        #    pass

    def rebuild(self):
        raise NotImplementedError()

    def add_psolid(self, card, comment):
        self.psolid.add(card, comment)

    def _get_types(self):
        return [self.psolid]

    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_bdf(self, f, size=8, pids=None):
        f.write('$PROPERTIES_SOLID\n')
        types = self._get_types()
        for prop in types:
            prop.write_bdf(f, size=size, pids=pids)