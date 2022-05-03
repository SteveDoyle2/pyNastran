import sys
from typing import List

def filter_no_args(msg: str, argv: List[str], quiet: bool=False):
    if len(argv) == 1:
        if quiet:
            sys.exit()
        sys.exit(msg)
