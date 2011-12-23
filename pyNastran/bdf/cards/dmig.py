import os
import sys

from pyNastran.bdf.cards.baseCard import BaseCard

class DMIG(BaseCard):
    """
    Defines direct input matrices related to grid, extra, and/or scalar points. The matrix
    is defined by a single header entry and one or more column entries. A column entry
    is required for each column with nonzero elements.
    """
    type = 'DMIG'
    def __init__(self,card=None,data=None):
        self.name = card.field(1)
        #zero
        self.ifo   = card.field(3)
        self.tin   = card.field(4)
        self.tout  = card.field(5,0)
        self.polar = card.field(6)
        self.ncol  = card.field(8)

        #self.name = card.field(1)
        #self.Gj = []
        #self.Cj = []
        #self.Gi = []
        #self.Ci = []
        #self.Ai = []
        #self.Bi = []
        
