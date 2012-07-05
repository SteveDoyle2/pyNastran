from pyNastran.bdf.cards.baseCard import BaseCard

class PARAM(BaseCard):
    type = 'PARAM'
    def __init__(self, card):
        self.key   = card.field(1)
        self.value = card.field(2)

    #def isSameCard(self, param):
        #fields1 = [self.key, self.value ]
        #fields2 = [param.key, param.value]
        #for (field1, field2) in zip(fields1, fields2):
        #    if not self.isSame(field1, field2):
        #        return False
        #    ###
        ####
        #return True

    def rawFields(self):
        fields = ['PARAM', self.key, self.value]
        return fields
