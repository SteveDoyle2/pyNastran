from six.moves import StringIO
from numpy import where, asarray
from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class RodElement(Element):
    def __init__(self, model):
        """
        Defines the CONROD object.

        :param self: the CONROD object
        :param model: the BDF object
        """
        Element.__init__(self, model)

    def __getitem__(self, element_id):
        #material_id = slice_to_iter(material_id)
        i = where(self.element_id == element_id)[0]
        return self.slice_by_index(i)

