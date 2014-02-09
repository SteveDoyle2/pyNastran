import StringIO

#from .pshell import PSHELL
#from .pcomp import PCOMP
#from .pcomp import PCOMP as PCOMPG
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG

class PropertiesShell(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        """
        self.model = model

        self.pshell = {}
        self.pcomp  = {}
        self.pcompg = {}

    def build(self):
        pass
        #self.pshell.build()
        #self.pcomp.build()
        #self.pcompg.build()
        
        #npshell = self.pshell.n
        #npcomp  = self.pcomp.n
        #npcompg = self.pcompg.n
        
        #if npshell and npcomp and npcompg:
            #asdf
        #if npshell and npcomp:
            #pid = concatenate(self.pshell.pid, self.pcomp.pid)
            #unique_pids = unique(pid)
            #print unique_pids
            #if len(unique_pids) != len(pid):
                #raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')
        #else:
            #pass

    def rebuild(self):
        return
        #raise NotImplementedError()

    def add_pshell(self, card, comment):
        card = PSHELL(card, comment=comment)
        self.pshell[card.property_id] = card

    def add_pcomp(self, card, comment):
        card = PCOMP(card, comment=comment)
        self.pcomp[card.property_id] = card

    def add_pcompg(self, card, comment):
        card = PCOMPG(card, comment=comment)
        self.pcompg[card.property_id] = card

    #=========================================================================
    def get_thickness(self, property_ids=None):
        """
        Gets the thickness of the PSHELLs/PCOMPs.
        
        :param self: the ShellProperties object
        :param property_ids: the property IDs to consider (default=None -> all)
        """
        _property_ids = concatenate(
            self.pshell.keys(),
           self.pcomp.keys(),
           self.pcompg.keys(), )
        t = concatenate(
            self.pshell.get_thickness(),
            self.pcomp.get_thickness(),
            self.pcompg.get_thickness(), )
        i = argsort(_pids)
        
        t2 = t[searchsorted(property_ids, _property_ids)]
        return t2

    #=========================================================================
    def _get_types(self):
        return [self.pshell, self.pcomp, self.pcompg]

    def get_stats(self):
        msg = []
        types = self._get_types()
        for prop in types:
            nprop = len(prop)
            if nprop:
                msg.append('  %-8s: %i' % (prop.type, nprop))
        return msg

    def write_bdf(self, f, size=8, property_ids=None):
        f.write('$PROPERTIES_SHELL\n')
        types = self._get_types()
        if property_ids is None:
            for props in types:
                for pid, prop in props:
                    prop.write_bdf(f, size=size)
        else:
            for props in types:
                for pid, prop in props:
                    if pid in property_ids:
                        prop.write_bdf(f, size=size)

    def __repr__(self):
        f = StringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue().rstrip()