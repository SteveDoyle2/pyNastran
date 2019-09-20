from pyNastran.bdf.bdf_interface.subcase_cards import IntCard


class ADAPT(IntCard):
    type = 'ADAPT'
    def __init__(self, value):
        super().__init__(value)

class AUXMODEL(IntCard):
    type = 'AUXMODEL'
    def __init__(self, value):
        super().__init__(value)

class BC(IntCard):
    type = 'BC'
    def __init__(self, value):
        super().__init__(value)

class BCSET(IntCard):
    type = 'BCSET'
    def __init__(self, value):
        super().__init__(value)

class BGSET(IntCard):
    type = 'BGSET'
    def __init__(self, value):
        super().__init__(value)

class BOLTLD(IntCard):
    type = 'BOLTLD'
    def __init__(self, value):
        super().__init__(value)

class CLOAD(IntCard):
    type = 'CLOAD'
    def __init__(self, value):
        super().__init__(value)

class CMETHOD(IntCard):
    type = 'CMETHOD'
    def __init__(self, value):
        super().__init__(value)

class CSSCHD(IntCard):
    type = 'CSSCHD'
    def __init__(self, value):
        super().__init__(value)

class DEFORM(IntCard):
    type = 'DEFORM'
    def __init__(self, value):
        super().__init__(value)

class DESGLB(IntCard):
    type = 'DESGLB'
    def __init__(self, value):
        super().__init__(value)

class DESSUB(IntCard):
    type = 'DESSUB'
    def __init__(self, value):
        super().__init__(value)

class DIVERG(IntCard):
    type = 'DIVERG'
    def __init__(self, value):
        super().__init__(value)

class DLOAD(IntCard):
    type = 'DLOAD'
    def __init__(self, value):
        super().__init__(value)

class DRSPAN(IntCard):
    type = 'DRSPAN'
    def __init__(self, value):
        super().__init__(value)

class DTEMP(IntCard):
    type = 'DTEMP'
    def __init__(self, value):
        super().__init__(value)

class EBDSET(IntCard):
    type = 'EBDSET'
    def __init__(self, value):
        super().__init__(value)

class FMETHOD(IntCard):
    type = 'FMETHOD'
    def __init__(self, value):
        super().__init__(value)

class FREQUENCY(IntCard):
    type = 'FREQUENCY'
    def __init__(self, value):
        super().__init__(value)

class GUST(IntCard):
    type = 'GUST'
    def __init__(self, value):
        super().__init__(value)

class IC(IntCard):
    type = 'IC'
    def __init__(self, value):
        super().__init__(value)

class LINE(IntCard):
    type = 'LINE'
    def __init__(self, value):
        super().__init__(value)

class LOAD(IntCard):
    type = 'LOAD'
    def __init__(self, value):
        super().__init__(value)

class LOADSET(IntCard):
    type = 'LOADSET'
    def __init__(self, value):
        super().__init__(value)

class MAXLINES(IntCard):
    type = 'MAXLINES'
    def __init__(self, value):
        super().__init__(value)

class METHOD(IntCard):
    type = 'METHOD'
    def __init__(self, value):
        super().__init__(value)

class MFLUID(IntCard):
    type = 'MFLUID'
    def __init__(self, value):
        super().__init__(value)

class MODES(IntCard):
    type = 'MODES'
    def __init__(self, value):
        super().__init__(value)

class MODTRAK(IntCard):
    type = 'MODTRAK'
    def __init__(self, value):
        super().__init__(value)

class MPC(IntCard):
    type = 'MPC'
    def __init__(self, value):
        super().__init__(value)

class NLCNTL(IntCard):
    type = 'NLCNTL'
    def __init__(self, value):
        super().__init__(value)

class NLPARM(IntCard):
    type = 'NLPARM'
    alternate_names = {'NLPAR'}
    def __init__(self, value):
        super().__init__(value)

class NONLINEAR(IntCard):
    type = 'NONLINEAR'
    def __init__(self, value):
        super().__init__(value)

class NSM(IntCard):
    type = 'NSM'
    def __init__(self, value):
        super().__init__(value)

class OUTRCV(IntCard):
    type = 'OUTRCV'
    def __init__(self, value):
        super().__init__(value)

class PARTN(IntCard):
    type = 'PARTN'
    def __init__(self, value):
        super().__init__(value)

class REPCASE(IntCard):
    type = 'REPCASE'
    def __init__(self, value):
        super().__init__(value)

class RSMETHOD(IntCard):
    type = 'RSMETHOD'
    def __init__(self, value):
        super().__init__(value)

class SEFINAL(IntCard):
    type = 'SEFINAL'
    def __init__(self, value):
        super().__init__(value)

class SMETHOD(IntCard):
    type = 'SMETHOD'
    def __init__(self, value):
        super().__init__(value)

class SPC(IntCard):
    type = 'SPC'
    def __init__(self, value):
        super().__init__(value)

class STATSUB(IntCard):
    type = 'STATSUB'
    def __init__(self, value):
        super().__init__(value)

class SUPORT1(IntCard):
    type = 'SUPORT1'
    def __init__(self, value):
        super().__init__(value)

class SYM(IntCard):
    type = 'SYM'
    def __init__(self, value):
        super().__init__(value)

class SYMCOM(IntCard):
    type = 'SYMCOM'
    def __init__(self, value):
        super().__init__(value)

class TRIM(IntCard):
    type = 'TRIM'
    def __init__(self, value):
        super().__init__(value)

class TSTEP(IntCard):
    type = 'TSTEP'
    def __init__(self, value):
        super().__init__(value)

class TSTEPNL(IntCard):
    type = 'TSTEPNL'
    def __init__(self, value):
        super().__init__(value)

class TSTRU(IntCard):
    type = 'TSTRU'
    def __init__(self, value):
        super().__init__(value)

#------------
#not sure
class DESOBJ(IntCard):
    type = 'DESOBJ'
    def __init__(self, value):
        super().__init__(value)

class SDAMP(IntCard):
    type = 'SDAMP'
    def __init__(self, value):
        super().__init__(value)

class FREQ(IntCard):
    type = 'FREQ'
    def __init__(self, value):
        super().__init__(value)

class DYNRED(IntCard):
    type = 'DYNRED'
    def __init__(self, value):
        super().__init__(value)

class TFL(IntCard):
    type = 'TFL'
    def __init__(self, value):
        super().__init__(value)

class SUPORT(IntCard):
    type = 'SUPORT'
    def __init__(self, value):
        super().__init__(value)

INT_CARDS = [
    ADAPT, AUXMODEL, BC, BCSET, BGSET, BOLTLD, CLOAD, CMETHOD, CSSCHD,
    DEFORM, DESGLB, DESSUB, DIVERG, DLOAD, DRSPAN, DESSUB, DTEMP,
    EBDSET, FMETHOD, FREQUENCY, GUST, IC, LINE, LOAD, LOADSET, MAXLINES,
    METHOD, MFLUID, MODES, MODTRAK, MPC, NLCNTL, NLPARM, NONLINEAR, NSM,
    OUTRCV, PARTN, REPCASE, RSMETHOD, SEFINAL, SMETHOD, SPC,
    SUPORT1, SYM, SYMCOM, TRIM, TSTEP, TSTEPNL, TSTRU, TFL, SUPORT,

    # not sure
    DESOBJ, SDAMP, STATSUB, FREQ, DYNRED,
]
INT_CARD_DICT = {card.type : card for card in INT_CARDS}
for card in INT_CARDS:
    if hasattr(card, 'alternate_names'):
        for name in card.alternate_names:
            INT_CARD_DICT[name] = card


#INT_CARD_NAMES = tuple([card.type for card in INT_CARDS])
INT_CARD_NAMES = tuple(INT_CARD_DICT.keys())
