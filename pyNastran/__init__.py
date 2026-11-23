# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
is_installed = False

import os
import sys

# Version scheme (PEP 440):
#   release : 1.5.0                       (IS_RELEASE = True)
#   dev     : 1.5.0.dev20261122+gad78709ca (.dev<commit date YYYYMMDD>+g<git hash>)
# A dev version sorts *below* the 1.5.0 release, so pip upgrades to the
# release once it exists.  pip skips dev versions unless asked, so depend
# on one with 'pyNastran>=1.5.0.dev0' (or use 'pip install --pre').
__version_release__ = '1.5.0'
IS_RELEASE = False
RELEASE_DATE = ('2026', '2', 'xx')  # (year, month, day); used for releases and without git

if is_pynastrangui_exe or is_installed:
    # pyInstaller / 'setup.py install' write a frozen version.py
    from pyNastran.version import __version__, __releaseDate__, __releaseDate2__
else:
    # this is still a requirement, but disabling it so readthedocs works
    if sys.version_info < (3, 10):  # pragma: no cover
        IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
        raise ImportError('Upgrade your Python to >= 3.10.0; version=(%s.%s.%s)' % (
            IMAJOR, MINOR1, MINOR2))

    def _git_hash_date() -> tuple[str, str]:
        """Returns the short git hash and commit date (YYYY-MM-DD) of the checkout.

        Returns ('', '') when pyNastran isn't in a git checkout (e.g., it was
        installed from a wheel) or git isn't installed.

        Notes
        -----
        One ``git log`` call gets both values (previously ``git describe`` +
        ``git show``), and git isn't run at all unless ``<repo>/.git`` exists
        (a folder, or a file for git worktrees).
        """
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(os.path.join(os.path.dirname(pkg_dir), '.git')):
            return '', ''
        import subprocess
        try:
            out = subprocess.check_output(
                ['git', 'log', '-1', '--format=%h %cs'],  # ad78709ca 2026-11-21
                cwd=pkg_dir, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            # git isn't installed or the repository is unreadable
            return '', ''
        ghash, date = out.decode('utf-8').split()
        return ghash, date

    def _installed_version() -> str:
        """Version from the installed package metadata ('' if not installed)."""
        from importlib.metadata import version, PackageNotFoundError
        try:
            return version('pyNastran')
        except PackageNotFoundError:
            return ''

    year, month, day = RELEASE_DATE
    if IS_RELEASE:
        __version__ = __version_release__
    else:
        ghash, date = _git_hash_date()
        if ghash:
            year, month, day = date.split('-')
            month, day = month.lstrip('0'), day.lstrip('0')
            __version__ = f'{__version_release__}.dev{date.replace("-", "")}+g{ghash}'
        else:
            # no git: report what pip installed (e.g., from a wheel)
            __version__ = _installed_version() or f'{__version_release__}.dev0'

    months = (
        'JANUARY', 'FEBRUARY', 'MARCH', 'APRIL', 'MAY', 'JUNE', 'JULY',
        'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER')
    __releaseDate__ = f'{year}/{month}/{day}'                    # 2024/3/11
    __releaseDate2__ = f'{months[int(month) - 1]} {day}, {year}'  # MARCH 11, 2024

__author__ = 'Steven Doyle'
__email__ = 'mesheb82@gmail.com'
__desc__ = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__ = 'BSD-3'
__copyright__ = f'Copyright {__license__}; 2011-2025'
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__ = 'https://github.com/SteveDoyle2/pyNastran'

DEV = 'dev' in __version__
IS_CI = (
    'TRAVIS' in os.environ or
    'TRAVIS_PYTHON_VERSION' in os.environ or
    'APPVEYOR' in os.environ or
    #'READTHEDOCS' in os.environ or
    'GITHUB_ACTOR' in os.environ
)
if DEV:
    __docs__ = 'https://pynastran-git.readthedocs.io/en/latest/quick_start/index.html'
else:
    # 1.3
    # we don't do separate doc releases for 1.3 vs 1.3.1
    __docs__ = f'https://pynastran-git.readthedocs.io/en/{__version__[:3]}/quick_start/index.html'

__issue__ = 'https://github.com/SteveDoyle2/pyNastran/issues'
__discussion_forum__ = 'https://groups.google.com/forum/#!forum/pynastran-discuss'

is_release = True  ## True=turns on skipping of tables that aren't supported
