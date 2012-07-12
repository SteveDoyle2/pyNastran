class BDFSyntaxError(SyntaxError):
    pass

class BDFRuntimeError(RuntimeError):
    pass

#-------------------------------------------------------
class InvalidSubcaseParseError(BDFSyntaxError):
    pass

class CardParseError(BDFSyntaxError):
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

class InvalidUnitVectorError(BDFRuntimeError):
    pass

class CardInstantiationError(BDFRuntimeError):
    pass

