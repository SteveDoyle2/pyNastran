from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class RodElement(Element):
    def __init__(self, model):
        """
        Defines the CONROD object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Element.__init__(self, model)
