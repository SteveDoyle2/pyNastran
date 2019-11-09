import os
import sys


# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
is_installed = False
if is_pynastrangui_exe or is_installed:
    # pyInstaller
    from pyNastran.version import __version__, __releaseDate__
else:
    import subprocess

    # this is still a requirement, but disabling it so readthedocs works
    if sys.version_info < (3, 7):  # pragma: no cover
        IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
        # makes sure we don't get the following bug:
        #   Issue #19099: The struct module now supports Unicode format strings.
        raise ImportError('Upgrade your Python to >= 3.7.0; version=(%s.%s.%s)' % (
            IMAJOR, MINOR1, MINOR2))

    def get_git_revision_short_hash():
        """determines the git revision; only works if the packages was checked
        out using git"""
        try:
            #ghash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

            # independent of pyNastran location as long as there is a git folder
            #   what about if you use setup_user.py install?
            #   what about if you don't have git?
            # can raise a subprocess.CalledProcessError, which means the return code != 0
            ghash = subprocess.check_output(['git', 'describe', '--always'],
                                            cwd=os.path.dirname(__file__))

            ghash = ghash.decode('utf-8').rstrip()
        except:
            # git isn't installed
            ghash = 'no.checksum.error'
        return 'dev.%s' % ghash

    revision = get_git_revision_short_hash()
    __version__ = '1.3.0+%s' % revision
    __releaseDate__ = '2019/6/xx'
    __releaseDate2__ = 'JUNE xx, 2019'

__author__ = 'Steven Doyle, Michael Redmond, Saullo Castro, Nikita Kalutsky'
__email__ = 'mesheb82@gmail.com'
__desc__ = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__longdesc__ = __desc__
__license__ = 'BSD-3'
__copyright__ = f'Copyright {__license__}; 2011-2019'
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__ = 'https://github.com/SteveDoyle2/pyNastran'
#__docs__ = 'http://pynastran.m4-engineering.com/master'  # still not setup...
if 'dev' in  __version__:
    __docs_rtd__ = 'https://pynastran-git.readthedocs.io/en/latest/quick_start/index.html'
    __docs__ = __docs_rtd__
else:
    # 1.2
    __docs_rtd__ = f'https://pynastran-git.readthedocs.io/en/{__version__[:3]}/quick_start/index.html'
    __docs__ = f'http://pynastran.m4-engineering.com/{__version__}'

__issue__ = 'https://github.com/SteveDoyle2/pyNastran/issues'
__discussion_forum__ = 'https://groups.google.com/forum/#!forum/pynastran-discuss'

is_release = True  ## True=turns on skipping of tables that aren't supported
