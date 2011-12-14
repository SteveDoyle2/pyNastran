class TabCharacterError(SyntaxError):
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
