from pyNastran.bdf.bdf_interface.subcase_cards import IntCard


class ADACT(IntCard):
    type = 'ADACT'
    #def __init__(self, value):
        #super().__init__(value)

class ADAPT(IntCard):
    type = 'ADAPT'
    #def __init__(self, value):
        #super().__init__(value)

class AUXMODEL(IntCard):
    type = 'AUXMODEL'
    #def __init__(self, value):
        #super().__init__(value)

class BC(IntCard):
    type = 'BC'
    #def __init__(self, value):
        #super().__init__(value)

class BCONTACT(IntCard):
    type = 'BCONTACT'
    #def __init__(self, value):
        #super().__init__(value)

class BCSET(IntCard):
    type = 'BCSET'
    #def __init__(self, value):
        #super().__init__(value)

class BGSET(IntCard):
    type = 'BGSET'
    #def __init__(self, value):
        #super().__init__(value)

class BOLTLD(IntCard):
    type = 'BOLTLD'
    #def __init__(self, value):
        #super().__init__(value)

class CLOAD(IntCard):
    type = 'CLOAD'
    #def __init__(self, value):
        #super().__init__(value)

class CMETHOD(IntCard):
    type = 'CMETHOD'
    #def __init__(self, value):
        #super().__init__(value)

class CSSCHD(IntCard):
    type = 'CSSCHD'
    #def __init__(self, value):
        #super().__init__(value)

class DEFORM(IntCard):
    type = 'DEFORM'
    #def __init__(self, value):
        #super().__init__(value)

class DESGLB(IntCard):
    type = 'DESGLB'
    #def __init__(self, value):
        #super().__init__(value)

class DESSUB(IntCard):
    type = 'DESSUB'
    #def __init__(self, value):
        #super().__init__(value)

class DIVERG(IntCard):
    type = 'DIVERG'
    #def __init__(self, value):
        #super().__init__(value)

class DLOAD(IntCard):
    type = 'DLOAD'
    #def __init__(self, value):
        #super().__init__(value)

class DRSPAN(IntCard):
    type = 'DRSPAN'
    #def __init__(self, value):
        #super().__init__(value)

class DTEMP(IntCard):
    type = 'DTEMP'
    #def __init__(self, value):
        #super().__init__(value)

class EBDSET(IntCard):
    type = 'EBDSET'
    #def __init__(self, value):
        #super().__init__(value)

class ELAR(IntCard):
    type = 'ELAR'
    #def __init__(self, value):
        #super().__init__(value)

class FMETHOD(IntCard):
    type = 'FMETHOD'
    #def __init__(self, value):
        #super().__init__(value)

class FREQUENCY(IntCard):
    type = 'FREQUENCY'
    alternate_names = {'FREQ'}
    #def __init__(self, value):
        #super().__init__(value)

class GUST(IntCard):
    type = 'GUST'
    #def __init__(self, value):
        #super().__init__(value)

class IC(IntCard):
    type = 'IC'
    #def __init__(self, value):
        #super().__init__(value)

class LINE(IntCard):
    type = 'LINE'
    #def __init__(self, value):
        #super().__init__(value)

class LOAD(IntCard):
    type = 'LOAD'
    #def __init__(self, value):
        #super().__init__(value)

class LOADSET(IntCard):
    type = 'LOADSET'
    #def __init__(self, value):
        #super().__init__(value)

class MAXLINES(IntCard):
    type = 'MAXLINES'
    #def __init__(self, value):
        #super().__init__(value)

class METHOD(IntCard):
    type = 'METHOD'
    alternate_names = {'METH', 'METHO'}
    #def __init__(self, value):
        #super().__init__(value)

class MFLUID(IntCard):
    type = 'MFLUID'
    #def __init__(self, value):
        #super().__init__(value)

class MODES(IntCard):
    type = 'MODES'
    #def __init__(self, value):
        #super().__init__(value)

class MODTRAK(IntCard):
    type = 'MODTRAK'
    #def __init__(self, value):
        #super().__init__(value)

class MPC(IntCard):
    type = 'MPC'
    #def __init__(self, value):
        #super().__init__(value)

class NLCNTL(IntCard):
    type = 'NLCNTL'
    #def __init__(self, value):
        #super().__init__(value)

class NLPARM(IntCard):
    type = 'NLPARM'
    alternate_names = {'NLPAR'}
    #def __init__(self, value):
        #super().__init__(value)

class NLPCI(IntCard):
    type = 'NLPCI'
    #alternate_names = {'NLPAR'}
    #def __init__(self, value):
        #super().__init__(value)

class NONLINEAR(IntCard):
    type = 'NONLINEAR'
    #def __init__(self, value):
        #super().__init__(value)

class NSM(IntCard):
    type = 'NSM'
    #def __init__(self, value):
        #super().__init__(value)

class OUTRCV(IntCard):
    type = 'OUTRCV'
    #def __init__(self, value):
        #super().__init__(value)

class PARTN(IntCard):
    type = 'PARTN'
    #def __init__(self, value):
        #super().__init__(value)

class REPCASE(IntCard):
    type = 'REPCASE'
    #def __init__(self, value):
        #super().__init__(value)


class RMETHOD(IntCard):
    type = 'RMETHOD'
    #def __init__(self, value):
        #super().__init__(value)

class RSMETHOD(IntCard):
    type = 'RSMETHOD'
    #def __init__(self, value):
        #super().__init__(value)

class SEFINAL(IntCard):
    type = 'SEFINAL'
    #def __init__(self, value):
        #super().__init__(value)

class SMETHOD(IntCard):
    type = 'SMETHOD'
    #def __init__(self, value):
        #super().__init__(value)

class SPC(IntCard):
    type = 'SPC'
    #def __init__(self, value):
        #super().__init__(value)

class STATSUB(IntCard):
    type = 'STATSUB'
    #def __init__(self, value):
        #super().__init__(value)

class SUPORT1(IntCard):
    type = 'SUPORT1'
    alternate_names = {'SUPORT'}
    #def __init__(self, value):
        #super().__init__(value)

class SYM(IntCard):
    type = 'SYM'
    #def __init__(self, value):
        #super().__init__(value)

class SYMCOM(IntCard):
    type = 'SYMCOM'
    #def __init__(self, value):
        #super().__init__(value)

class TRIM(IntCard):
    type = 'TRIM'
    #def __init__(self, value):
        #super().__init__(value)

class TSTEP(IntCard):
    type = 'TSTEP'
    #def __init__(self, value):
        #super().__init__(value)

class TSTEPNL(IntCard):
    type = 'TSTEPNL'
    #def __init__(self, value):
        #super().__init__(value)

class TSTRU(IntCard):
    type = 'TSTRU'
    #def __init__(self, value):
        #super().__init__(value)

class RANDOM(IntCard):
    type = 'RANDOM'
    #def __init__(self, value):
        #super().__init__(value)

#------------
#not sure
class DESOBJ(IntCard):
    type = 'DESOBJ'
    #def __init__(self, value):
        #super().__init__(value)

class SDAMPING(IntCard):
    type = 'SDAMPING'
    alternate_names = {'SDAMP'}
    #def __init__(self, value):
        #super().__init__(value)

class DYNRED(IntCard):
    type = 'DYNRED'
    #def __init__(self, value):
        #super().__init__(value)

class TFL(IntCard):
    type = 'TFL'
    #def __init__(self, value):
        #super().__init__(value)

class MODESELECT(IntCard):
    type = 'MODESELECT'
    alternate_names = {'MODSEL'}
    #def __init__(self, value):
        #super().__init__(value)

class OTIME(IntCard):
    type = 'OTIME'
    #def __init__(self, value):
        #super().__init__(value)

#class OFREQUENCY(IntCard):
    #type = 'OFREQUENCY'
    #alternate_names = {'OFREQ'}
    #def __init__(self, value):
        #super().__init__(value)

class OMODES(IntCard):
    type = 'OMODES'
    #def __init__(self, value):
        #super().__init__(value)

class RGYRO(IntCard):
    type = 'RGYRO'
    #def __init__(self, value):
        #super().__init__(value)

class SEDV(IntCard):
    type = 'SEDV'


INT_CARDS = [
    ADACT, ADAPT, AUXMODEL, BCONTACT, BC, BCSET, BGSET, BOLTLD, CLOAD, CMETHOD, CSSCHD,
    DEFORM, DESGLB, DESSUB, DIVERG, DLOAD, DRSPAN, DESSUB, DTEMP,
    EBDSET, ELAR, FMETHOD, FREQUENCY,
    GUST, IC, LINE, LOAD, LOADSET,
    MAXLINES, METHOD, MFLUID, MODES, MODTRAK, MPC, MODESELECT,
    NLCNTL, NLPARM, NLPCI, NONLINEAR, NSM,
    OTIME, OMODES, OUTRCV, PARTN,
    RANDOM, RGYRO, REPCASE, RSMETHOD, RMETHOD,
    SUPORT1, SYM, SYMCOM, SEFINAL, SMETHOD, SPC, SEDV,
    TRIM, TSTEP, TSTEPNL, TSTRU, TFL,

    # not sure
    DESOBJ, SDAMPING, STATSUB, DYNRED,
]
INT_CARD_DICT = {card.type : card for card in INT_CARDS}
for card in INT_CARDS:
    if hasattr(card, 'alternate_names'):
        for name in card.alternate_names:
            INT_CARD_DICT[name] = card


#INT_CARD_NAMES = tuple([card.type for card in INT_CARDS])
INT_CARD_NAMES = tuple(INT_CARD_DICT.keys())
