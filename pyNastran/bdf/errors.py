class BDF_SyntaxError(SyntaxError):
    pass

class TabCharacterError(SyntaxError):
    pass

class TabCommaCharacterError(SyntaxError):
    pass

class ClosedBDFError(RuntimeError):
    pass

class MissingFileError(RuntimeError):
    pass

class ParamParseError(SyntaxError):
    pass

class InvalidSubcaseParseError(SyntaxError):
    pass

#-------------------------------------------------------
class CardParseError(SyntaxError):
    pass

class FloatScientificCardParseError(CardParseError):
    pass

class ScientificCardParseError(CardParseError):
    pass

class WhitespaceCardParseError(CardParseError):
    pass

#-----------------------------------------------------

class InvalidUnitVectorError(RuntimeError):
    pass

class CardInstantiationError(RuntimeError):
    pass

#-----------------------------------------------------
class StiffnessMatrixError(RuntimeError):
    pass

class InvalidRequestError(RuntimeError):
    pass

class InvalidFieldError(RuntimeError):
    pass

class InvalidResultCode(NotImplementedError):
    pass

class CoordTypeError(TypeError):
    pass

