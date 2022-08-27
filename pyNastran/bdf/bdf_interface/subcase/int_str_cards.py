from .subcase_base import IntStrCard

class ADACT(IntStrCard):
    type = 'ADACT'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super().__init__(value)

class AEROF(IntStrCard):
    type = 'AEROF'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class APRES(IntStrCard):
    type = 'APRES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPRSORT(IntStrCard):
    type = 'GPRSORT'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPSDCON(IntStrCard):
    type = 'GPSDCON'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class HARMONICS(IntStrCard):
    type = 'HARMONICS'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super().__init__(value)

class OFREQUENCY(IntStrCard):
    type = 'OFREQUENCY'
    alternate_names = {'OFREQ'}
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class OMODES(IntStrCard):
    type = 'OMODES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SUPER(IntStrCard):
    type = 'SUPER'
    allowed_strings = {'ALL'}
    #def __init__(self, value):
        #super().__init__(value)

#----------------------------

class GPSTRESS(IntStrCard):
    type = 'GPSTRESS'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SEALL(IntStrCard):
    type = 'SEALL'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SEDR(IntStrCard):
    type = 'SEDR'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPKE(IntStrCard):
    type = 'GPKE'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

INTSTR_CARDS = [
    ADACT, AEROF, APRES, GPRSORT, GPSDCON, HARMONICS, OFREQUENCY, OMODES,
    SUPER, SEALL, SEDR,
] + [GPSTRESS, GPKE, ]
INTSTR_CARD_DICT = {card.type : card for card in INTSTR_CARDS}
INTSTR_CARD_NAMES = tuple([card.type for card in INTSTR_CARDS])
