from numpy import array, zeros, argsort, concatenate, searchsorted, unique, where, nan, arange

from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, double,
    double_or_blank, string_or_blank)


def divide_2d_array_by_column_vector(array_obj, vector):
    vector = array_obj.copy()
    vector[:, 0] / vector
    vector[:, 1] / vector
    vector[:, 2] / vector
    return v

class ElementsRod(object):
    def __init__(self, model):
        """
        Defines the ElementsRod object.

        :param self: the ElementsRod object
        :param model: the BDF object
        """
        aaa
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

        if npshell and npcomp:
            eid = concatenate(self.pshell.eid, self.pcomp.eid)
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

    def write_bdf(self, f, size=8, eids=None):
        #f.write('$PROPERTIES\n')
        self.crod.write_bdf(f, size=size, eid=eid)
        self.conrod.write_bdf(f, size=size, eid=eid)