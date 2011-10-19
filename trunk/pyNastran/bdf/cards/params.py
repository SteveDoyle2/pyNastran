# my code
from baseCard import BaseCard

class PARAM(BaseCard):
    def __init__(self,card):
        self.key   = card.field(1)
        self.value = card.field(2)

    def __repr__(self):
        fields = ['PARAM',self.key,self.value]
        return self.printCard(fields)
