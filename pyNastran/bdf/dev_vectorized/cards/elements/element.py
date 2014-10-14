class Element(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % (self.type, self.n))
        return msg

    def __len__(self):
        return self.n

    def __iter__(self):
        eids = self.element_id
        for eid in eids:
            yield eid

    def values(self):
        eids = self.element_id
        for eid in eids:
            yield self.__getitem__(eid)

    def items(self):
        eids = self.element_id
        for eid in eids:
            yield eid, self.__getitem__(eid)
