from pyNastran.bdf.bdf_interface.subcase_cards import StringCard

class AESYMXY(StringCard):
    type = 'AESYMXY'
    short_name = 'AESYMXY'
    allowed_values = ['SYMMETRIC', 'ANTISYMMETRIC', 'ASYMMTRIC']
    def __init__(self, value):
        super().__init__(value)

class AESYMXZ(StringCard):
    type = 'AESYMXZ'
    short_name = 'AESYMXZ'
    allowed_values = ['SYMMETRIC', 'ANTISYMMETRIC', 'ASYMMTRIC']
    def __init__(self, value):
        super().__init__(value)

class AXISYMMETRIC(StringCard):
    type = 'AXISYMMETRIC'
    short_name = type[:4]
    allowed_values = ['SINE', 'COSINE', 'FLUID']
    def __init__(self, value):
        super().__init__(value)

class AUTOSPC(StringCard):
    type = 'AUTOSPC'
    short_name = type[:4]
    allowed_values = ['YES', 'NO']
    def __init__(self, value):
        super().__init__(value)

class DSYM(StringCard):
    type = 'DSYM'
    short_name = type[:4]
    allowed_values = ['S', 'A', 'SS', 'SA', 'AS', 'AA']
    def __init__(self, value):
        super().__init__(value)

class SEQDEP(StringCard):
    type = 'SEQDEP'
    short_name = type[:4]
    allowed_values = ['YES', 'NO']
    def __init__(self, value):
        super().__init__(value)

#----------
# special

class K2PP(StringCard):
    type = 'K2PP'
    short_name = type[:4]
    def __init__(self, value):
        super().__init__(value, validate=False)

class ECHO(StringCard):
    """
    ECHO = NONE
    """
    type = 'ECHO'
    short_name = type[:4]
    allowed_values = ['BOTH', 'SORT', 'UNSORT', 'NONE', 'NOSORT', 'PUNCH']
    def __init__(self, value):
        super().__init__(value)

class ANALYSIS(StringCard):
    """ANALYSIS = HEAT"""
    type = 'ANALYSIS'
    short_name = type[:4]
    allowed_values = ['HEAT', 'STATICS', 'MODES', 'NLSTATIC', 'SAERO', 'DIVERG', 'FLUTTER',
                      'DCEIG', 'MCEIG', #  what?
                      'MFREQ', 'DFREQ', 'NLTRAN', 'MTRAN',
                      'STATIC', 'BUCK', 'FLUT', 'MTRANS'] # duplicates
    def __init__(self, value):
        super().__init__(value)

class THERMAL(StringCard): #  ???
    type = 'THERMAL'
    short_name = type[:4]
    allowed_values = ['ALL']
    def __init__(self, value):
        super().__init__(value)

class AECONFIG(StringCard): #  ???
    type = 'AECONFIG'
    short_name = type# [:4]
    allowed_values = ['FSWBASE', 'FSWHALF', 'FREEDOM4', 'RAERO', ]
    def __init__(self, value):
        super().__init__(value)

class RIGID(StringCard):
    type = 'RIGID'
    short_name = type
    allowed_values = ['LINEAR', 'LGELIM', 'LAGR', 'AUTO', 'LAGRANGE']
    def __init__(self, value):
        super().__init__(value)


special_cards = {'AESYMXY', 'AESYMXZ', 'AECONFIG'}
STR_CARDS = [AESYMXY, AESYMXZ, AECONFIG, AXISYMMETRIC, AUTOSPC, DSYM, SEQDEP] + [ECHO, ANALYSIS, K2PP, THERMAL, RIGID]

if 0:
    STR_CARD_DICT = {}
    STR_CARD_NAMES = []
    for card in STR_CARDS:
        if card.type in special_cards:
            name = card.type
        else:
            name = card.type[:4]
        STR_CARD_DICT[name] = card
        STR_CARD_NAMES.append(name)
#STR_CARD_DICT = {card.type[:4] : card if card.type not in special_cards else card.type : card
                 #for card in STR_CARDS}
#STR_CARD_NAMES = tuple([card.type for card in STR_CARDS])

STR_CARD_DICT = {card.short_name : card for card in STR_CARDS}
STR_CARD_NAMES = tuple([card.short_name for card in STR_CARDS])

STR_CARD_NAMES = tuple(STR_CARD_NAMES)
assert len(STR_CARD_DICT) == len(STR_CARD_DICT), f'ndict={len(STR_CARD_DICT)} STR_CARD_DICT.keys()={STR_CARD_DICT.keys()}'
