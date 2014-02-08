from numpy import array, concatenate, searchsorted, unique

# ..todo:: incomplete

class Coord(object):
    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        :param pshells: the list of PSHELL cards
        :param pcomps: the list of PCOMP cards
        :param pshears: the list of PSHEAR cards
        """
        self.model = model
        
        self.cids = [0]
        #self.cord2r = CORD2R()
        #self.cord2c = CORD2C()
        #self.cord2s = CORD2S()

        #cid = concatenate(pshell.cid, pcomp.cid)
        #unique_cids = unique(cid)
        #if unique_cids != len(cid):
        #    raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')

    def build(self):
        pass
        #self.coord2r