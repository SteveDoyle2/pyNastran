import cStringIO

class RodElement(object):
    def __init__(self, model):
        """
        Defines the CONROD object.

        :param self: the CONROD object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<%s object> n=%s\n' % (self.type, self.n))
        self.write_bdf(f)
        return f.getvalue()