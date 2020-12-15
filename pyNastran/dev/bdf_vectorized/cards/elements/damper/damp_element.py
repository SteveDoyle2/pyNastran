from pyNastran.dev.bdf_vectorized.cards.elements.element import Element

class DamperElement(Element):
    def __init__(self, model):
        """
        Defines the DamperElement object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Element.__init__(self, model)
