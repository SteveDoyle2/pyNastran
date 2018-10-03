import os
import sys

# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
if is_pynastrangui_exe:
    # pyInstaller
    from pyNastran.version import __version__, __releaseDate__
else:
    import subprocess

    def get_git_revision_short_hash():
        try:
            #ghash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

            # independent of pyNastran location as long as there is a git folder
            #   what about if you use setup_user.py install?
            #   what about if you don't have git?
            # can raise a subprocess.CalledProcessError, which means the return code != 0
            ghash = subprocess.check_output(['git', 'describe', '--always'], cwd=os.path.dirname(__file__))

            ghash = ghash.decode('utf-8').rstrip()
        except:
            # git isn't installed
            ghash = 'no.checksum.error'
        return 'dev.%s' % ghash

    revision = get_git_revision_short_hash()
    __version__ = '1.2.0+%s' % revision
    __releaseDate__ = '2018/6/xx'
    __releaseDate2__ = 'JUNE xx, 2018'

__author__  = 'Steven Doyle, Michael Redmond, Saullo Castro, hurlei, Paul Blelloch, Nikita Kalutsky'
__email__ = 'mesheb82@gmail.com'
__desc__ = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__longdesc__ = __desc__
__license__ = 'BSD-3'
__copyright__ = 'Copyright %s; 2011-2018' % __license__
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__ = 'https://github.com/SteveDoyle2/pyNastran'
#__docs__ = 'http://pynastran.m4-engineering.com/master'  # still not setup...
__docs__ = 'https://pynastran-git.readthedocs.io/en/latest/quick_start/index.html'
__issue__ = 'https://github.com/SteveDoyle2/pyNastran/issues'
__discussion_forum__ = 'https://groups.google.com/forum/#!forum/pynastran-discuss'

is_release = True  ## True=turns on skipping of tables that aren't supported
