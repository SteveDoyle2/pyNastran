from typing import Optional, Any
from .gui_result import GuiResult

Formi = tuple[str, Optional[int], Any]
Form = list[Formi]
FormDict = dict[tuple[Any, Any], Form]
HeaderDict = dict[tuple[Any, Any], str]
KeysMap = dict[str, Any]
Case = tuple[GuiResult, tuple[int, str]]
Cases = dict[int, Case]
