class BDFSyntaxError(SyntaxError):
    pass

class BDFRuntimeError(RuntimeError):
    pass

class BDFIOError(IOError):
    pass


#-------------------------------------------------------
class BDF_SyntaxError(BDFSyntaxError):
    pass

class TabCharacterError(BDFSyntaxError):
    pass

class TabCommaCharacterError(BDFSyntaxError):
    pass

#-------------------------------------------------------
class ClosedBDFError(BDFIOError):
    pass

class MissingFileError(BDFIOError):
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

#-----------------------------------------------------
class StiffnessMatrixError(BDFRuntimeError):
    pass

class InvalidRequestError(BDFRuntimeError):
    pass

class InvalidFieldError(BDFRuntimeError):
    pass

#-------------------------------------------------------
class InvalidResultCode(NotImplementedError):
    pass

class CoordTypeError(TypeError):
    pass

