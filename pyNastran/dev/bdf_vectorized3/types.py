from io import StringIO
from typing import TextIO

__all__ = ['TextIOLike']
TextIOLike = TextIO | StringIO
