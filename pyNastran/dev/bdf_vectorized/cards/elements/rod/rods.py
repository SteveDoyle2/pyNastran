from numpy import concatenate, unique


def divide_2d_array_by_column_vector(array_obj, vector):
    vector = array_obj.copy()
    vector[:, 0] / vector
    vector[:, 1] / vector
    vector[:, 2] / vector
    return vector

class ElementsRod:
    def __init__(self, model):
        """
        Defines the ElementsRod object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model

        self.crod = CROD(self.model)
        self.conrod = CONROD(self.model)
        self.ctube = CTUBE(self.model)

    def build(self):
        self.crod.build(self._crod)
        self.conrod.build(self._conrod)
        self.ctube.build(self._ctube)

        ncrod = len(self.crod.eid)
        nconrod = len(self.conrod.eid)
        nctube = len(self.ctube.eid)

        if ncrod or nconrod or nctube:
            eid = concatenate(self.crod.eid, self.conrod.eid, self.ctube.eid)
            unique_eids = unique(eid)
            if len(unique_eids) != len(eid):
                raise RuntimeError('There are duplicate CROD/CONROD IDs...')

    def rebuild(self):
        raise NotImplementedError()

    def add_conrod(self, card, comment):
        self.conrod.add(card, comment)

    def add_crod(self, card, comment):
        self.crod.add(card, comment)

    def add_ctube(self, card, comment):
        self.ctube.add(card, comment)

    def get_stats(self):
        msg = []
        types = [self.crod, self.conrod, self.ctube]
        for elem in types:
            nele = len(elem.eid)
            if nele:
                msg.append('  %-8s: %i' % (elem.type, nele))
        return msg

    def write_card(self, bdf_file, size=8, eids=None):
        #bdf_file.write('$PROPERTIES\n')
        self.crod.write_card(bdf_file, size=size, eid=eid)
        self.conrod.write_card(bdf_file, size=size, eid=eid)
        self.ctube.write_card(bdf_file, size=size, eid=eid)
