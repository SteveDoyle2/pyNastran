"""
defines:
 - cmd_line_plot_flutter()
 - cmd_line_plot_trim()
 - cmd_line_plot_optimization()

"""
from __future__ import annotations
import sys
import time
from typing import Optional, TYPE_CHECKING

from pyNastran.op2.tools.op2_merge import cmd_line_merge
from pyNastran.op2.tools.modal_combine import cmd_line_flutter_combine
# from pyNastran.op2.tools.envelope import cmd_line_envelope

if TYPE_CHECKING:
    from cpylog import SimpleLogger


def cmd_line(argv=None, log: Optional[SimpleLogger]=None) -> None:
    """the interface to ``f06`` on the command line"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        # USAGE_144 + USAGE_145 + USAGE_200 +
        # '\n'
        '  op2 modal_combine -h | --help\n'
        '  op2 envelope -h | --help\n'
        '  op2 merge -h | --help\n'
        '  op2 -v | --version\n'
        '\n'
    )
    if len(argv) == 1:
        sys.exit(msg)

    #assert sys.argv[0] != 'bdf', msg
    t0 = time.time()
    from cpylog import SimpleLogger
    log = SimpleLogger(level='debug')
    if argv[1] == 'merge':
        cmd_line_merge(argv=argv)
    elif argv[1] == 'modal_combine':
        cmd_line_flutter_combine(argv=argv, log=log)
    elif argv[1] == 'envelope':
        cmd_line_envelope(argv=argv, log=log)
    else:  # pragma: no cover
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % argv[1])
    print(f'dt = {time.time() - t0:.0f}')

if __name__ == '__main__':  # pragma: no cover
    cmd_line()
