import cStringIO


class SpringElement(object):

    def __len__(self):
        return self.n

    def __init__(self, model):
        """
        Defines the CELAS1 object.

        :param self: the CELAS1 object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % (self.type, self.n))
        return msg

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<%s object> n=%s\n' % (self.type, self.n))
        self.write_bdf(f)
        return f.getvalue()