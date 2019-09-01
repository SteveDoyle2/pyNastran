class CrossReferenceError(RuntimeError):
    pass


class CardParseSyntaxError(SyntaxError):
    """
    Class that is used for testing.
    Users should just treat this as a SyntaxError.
    """
    pass

class DuplicateIDsError(RuntimeError):
    pass

class MissingDeckSections(RuntimeError):
    pass

class UnsupportedCard(NotImplementedError):
    pass

class DisabledCardError(RuntimeError):
    """lets bdf_test.py flag cards as auto-crashing and then skipping the deck (e.g., CGEN)"""
    pass

class SuperelementFlagError(SyntaxError):
    pass

class ReplicationError(SyntaxError):
    pass

class EnvironmentVariableError(SyntaxError):
    pass
