from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class SpringElement(Element):
    def __init__(self, model):
        """
        Defines the SpringElement object.

        :param model: the BDF object
        """
        Element.__init__(self, model)
