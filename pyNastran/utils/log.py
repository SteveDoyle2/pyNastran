"""defines the make_log function"""
# coding: utf-8
import sys
import platform


def make_log():
    """
    Creates 'pyNastran.log' file with information about working environment,
    such as Python version, platform, architecture, etc. Useful for debugging.

    Returns
    -------
    msg : str
        the same string that goes to the log
    """
    smsg = [('sys.version', sys.version), ('sys.version_info', sys.version_info)]
    pmsg = [
        'machine', 'platform', 'processor', 'architecture', 'python_branch',
        'python_revision', 'win32_ver', 'version', 'uname', 'system',
        'python_build', 'python_compiler', 'python_implementation', 'system',
        'mac_ver', 'libc_ver', #'linux_distribution',
    ]

    fmt = '%-{0}s = %s\n'.format(max(map(len, pmsg + [j[0] for j in smsg])))
    msg = ''.join([fmt % (i, str(j).replace('\n', '; ')) for (i, j) in smsg])
    msg += ''.join([fmt % (i, str(getattr(platform, i)())) for i in pmsg])

    with open('pyNastran.log', 'w') as fil:
        fil.write(msg)
    return msg


if __name__ == '__main__':  # pragma: no cover
    make_log()
