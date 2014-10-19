import cStringIO

class VectorizedCard(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []
        if self.type in model._element_name_to_element_type_mapper:
            self.op2_id = model._element_name_to_element_type_mapper[self.type]
        else:
            if self.type.startswith('C'):
                print('there is no op2_id to apply for element=%r' % self.type)

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % (self.type, self.n))
        return msg

    def __repr__(self):
        f = cStringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue().rstrip()

