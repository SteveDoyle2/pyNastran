class CardParseError(SyntaxError):
    pass

#-------------------------------------------------------
class ParamParseError(CardParseError):
    pass

class FloatScientificCardParseError(CardParseError):
    pass

class ScientificCardParseError(CardParseError):
    pass

class WhitespaceCardParseError(CardParseError):
    pass

#-----------------------------------------------------

class CardInstantiationError(RuntimeError):
    pass

