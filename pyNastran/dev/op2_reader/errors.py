from pyNastran.f06.errors import FatalError

class FortranMarkerError(Exception):
    pass

class SortCodeError(RuntimeError):
    pass

class DeviceCodeError(SyntaxError):
    pass

class MultipleSolutionNotImplementedError(NotImplementedError):
    pass
