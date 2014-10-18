import cStringIO
from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class RodElement(Element):
    def __init__(self, model):
        """
        Defines the CONROD object.

        :param self: the CONROD object
        :param model: the BDF object
        """
        Element.__init__(self, model)
