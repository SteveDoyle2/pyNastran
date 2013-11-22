from .pshell import PSHELL
from .pcomp import PCOMP
from .pshear import PSHEAR

class PropertiesShell(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        """
        self.model = model

        #: the list of PSHELL cards
        self._pshell = []

        #: the list of PCOMP cards
        self._pcomp = []

        #: the list of PSHEAR cards
        self._pshear = []

        self._pshell_comment = []
        self._pcomp_comment = []
        self._pshear_comment = []

        self.pshell = PSHELL(self.model)
        self.pcomp = PCOMP(self.model)
        self.pshear = PSHEAR(self.model)

    def build(self):
        self.pshell.build(self._pshell)
        self.pcomp.build(self._pcomp)
        self.pshear.build(self._pshear)
        
        self._pshell = []
        self._pcomp = []
        self._pshear = []

        self._pshell_comment = []
        self._pcomp_comment = []
        self._pshear_comment = []

        npcomp = self.pcomp.n
        npshell = self.pshell.n
        npshear = self.pshear.n
        
        if npshell and npcomp and pshear:
            asf
        elif npshell and npcomp:
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
        self._pshell.append(card)
        self._pshell_comment.append(comment)

    def add_pcomp(self, card, comment):
        self._pcomp.append(card)
        self._pcomp_comment.append(comment)

    def add_pshear(self, card, comment):
        self._pshear.append(card)
        self._pshear_comment.append(comment)

    def get_stats(self):
        msg = []
        types = [self.pshell, self.pcomp, self.pshear]
        for prop in types:
            nprop = prop.n
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def get_thickness(self, pids):
        """
        Gets the thickness of the PSHELLs/PCOMPs.
        
        :param self: the ShellProperties object
        :param pids: the property IDs to consider (default=None -> all)
        """
        _pids = concatenate(self.pshell.pid,
                           self.pcomp.pid)
        t = concatenate(self.pshell.get_thickness(),
                        self.pcomp.get_thickness() )
        i = argsort(_pids)
        
        t2 = t[searchsorted(pids, _pids)]
        return t2

    def write_bdf(self, f, size=8, pids=None):
        f.write('$PROPERTIES\n')
        self.pshell.write_bdf(f, size=size, pids=pids)
        self.pcomp.write_bdf(f, size=size, pids=pids)