class CardParseError(SyntaxError):
    pass

#-------------------------------------------------------
class ParamParseError(CardParseError):
    pass

#-----------------------------------------------------

class CardInstantiationError(RuntimeError):
    pass

