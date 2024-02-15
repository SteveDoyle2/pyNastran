from io import StringIO
from typing import TextIO, Union

__all__ = ['TextIOLike']
TextIOLike = Union[TextIO, StringIO]
