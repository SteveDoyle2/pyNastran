class PCOMP(object):
    type = 'PCOMP'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
        self.model = model
        self.n = 0

    def build(self, cards):
        ncards = len(cards)
        self.n = ncards
        if ncards:
            #: Property ID
            self.pid = zeros(ncards, 'int32')

            ncards = len(cards)
            for card in cards:
                self.pid[i] = integer(card, 1, 'pid')
            i = self.pid.argsort()
            self.pid = self.pid[i]
            unique_pids = unique(self.pid)

            if len(unique_pids) != len(self.pid):
                raise RuntimeError('There are duplicate PSHELL IDs...')
        
    def write_bdf(self, f, size=8, pids=None):
        pass

