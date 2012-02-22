class BDF_SyntaxError(SyntaxError):
    pass

class TabCharacterError(SyntaxError):
    pass

class ClosedBDFError(RuntimeError):
    pass

class MissingFileError(RuntimeError):
    pass

class ParamParseError(SyntaxError):
    pass

class InvalidSubcaseParseError(SyntaxError):
    pass

class FloatScientificParseError(SyntaxError):
    pass

class ScientificParseError(SyntaxError):
    pass


class CardInstantiationError(RuntimeError):
    pass

class NotImplementedMethodError(NotImplementedError):
    pass

class StiffnessMatrixError(RuntimeError):
    pass

class InvalidFieldError(RuntimeError):
    pass

class InvalidResultCode(NotImplementedError):
    pass

class CoordTypeError(TypeError):
    pass

