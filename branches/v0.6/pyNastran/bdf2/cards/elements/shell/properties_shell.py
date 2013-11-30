from .pshell import PSHELL
from .pcomp import PCOMP
from .pcomp import PCOMP as PCOMPG

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
            pid = concatenate(self.pshell.pid, self.pcomp.pid)
            unique_pids = unique(pid)
            print unique_pids
            if len(unique_pids) != len(pid):
                raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        else:
            pass

    def rebuild(self):
        raise NotImplementedError()

    def add_pshell(self, card, comment):
        self.pshell.add(card, comment)

    def add_pcomp(self, card, comment):
        self.pcomp.add(card, comment)

    def add_pcompg(self, card, comment):
        raise NotImplementedError(card)
        self.pcompg.add(card, comment)

    def _get_types(self):
        return [self.pshell, self.pcomp, self.pcompg]

    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def get_thickness(self, property_ids):
        """
        Gets the thickness of the PSHELLs/PCOMPs.
        
        :param self: the ShellProperties object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        _pids = concatenate(self.pshell.property_id,
                           self.pcomp.property_id)
        t = concatenate(self.pshell.get_thickness(),
                        self.pcomp.get_thickness() )
        i = argsort(_pids)
        
        t2 = t[searchsorted(pids, _pids)]
        return t2

    def write_bdf(self, f, size=8, property_ids=None):
        f.write('$PROPERTIES_SHELL\n')
        types = self._get_types()
        for prop in types:
            #print('*SHELL', prop.type)
            prop.write_bdf(f, size=size, property_ids=property_ids)