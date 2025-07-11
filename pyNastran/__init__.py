# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
is_installed = False
if is_pynastrangui_exe or is_installed:
    # pyInstaller

    #from importlib.metadata import version, PackageNotFoundError
    #__version__ = version('pyNastran')
    #try:
        #__version__ = version('pyNastran')
    #except PackageNotFoundError:
        ## package is not installed
        #pass
    from pyNastran.version import __version__, __releaseDate__, __releaseDate2__
else:
    import os
    import sys
    import subprocess
    # this is still a requirement, but disabling it so readthedocs works
    if sys.version_info < (3, 10):  # pragma: no cover
        IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
        raise ImportError('Upgrade your Python to >= 3.10.0; version=(%s.%s.%s)' % (
            IMAJOR, MINOR1, MINOR2))

    def get_git_revision_short_hash_date() -> tuple[str, tuple[str, str, str]]:
        """determines the git revision; only works if the packages was checked
        out using git"""
        try:
            #ghash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

            # independent of pyNastran location as long as there is a git folder
            #   what about if you use setup_user.py install?
            #   what about if you don't have git?
            # can raise a subprocess.CalledProcessError, which means the return code != 0
            ghash_bytes = subprocess.check_output(['git', 'describe', '--always'],
                                                  cwd=os.path.dirname(__file__))
            ghash = ghash_bytes.decode('utf-8').rstrip()

            #git show -s --format=%cd --date=short
            #Mon Mar 11 23:00:39 2024 -0700
            #
            #git show --no-patch --format=%ai
            #2024-03-11 16:37:20 -0700
            date_bytes = subprocess.check_output(
                ['git', 'show', '-s', '--format=%cd', '--date=short'],
                cwd=os.path.dirname(__file__))
            date = date_bytes.decode('utf-8').rstrip()
            year, month, day = date.split('-')
            month = month.lstrip('0')
            day = day.lstrip('0')
            date_out = (year, month, day)
        except Exception:
            # git isn't installed
            ghash = 'no.git.checksum'
            date_out = ('2025', '5', 'xx')
        # 1.5.0+dev.0eccfa918
        return 'dev.%s' % ghash, date_out

    revision, (year, month, day) = get_git_revision_short_hash_date()
    __version_release__ = '1.5.0'

    # only for release; 1.4.0
    #__version__ = __version_release__
    # 1.4.0+dev.0eccfa918
    __version__ = f'{__version_release__}+{revision}'

    months = {
        '1': 'JANUARY', '2': 'FEBRUARY', '3': 'MARCH',
        '4': 'APRIL', '5': 'MAY', '6': 'JUNE',
        '7': 'JULY', '8': 'AUGUST', '9': 'SEPTEMBER',
        '10': 'OCTOBER', '11': 'NOVEMBER', '12': 'DECEMBER',
    }
    month_str = months[month]
    __releaseDate__ = f'{year}/{month}/{day}'        # 2024/3/11
    __releaseDate2__ = f'{month_str} {day}, {year}'  # MARCH 11, 2024

__author__ = 'Steven Doyle'
__email__ = 'mesheb82@gmail.com'
__desc__ = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__ = 'BSD-3'
__copyright__ = f'Copyright {__license__}; 2011-2025'
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__ = 'https://github.com/SteveDoyle2/pyNastran'

DEV = 'dev' in __version__
if DEV:
    __docs__ = 'https://pynastran-git.readthedocs.io/en/latest/quick_start/index.html'
else:
    # 1.3
    # we don't do separate doc releases for 1.3 vs 1.3.1
    __docs__ = f'https://pynastran-git.readthedocs.io/en/{__version__[:3]}/quick_start/index.html'

__issue__ = 'https://github.com/SteveDoyle2/pyNastran/issues'
__discussion_forum__ = 'https://groups.google.com/forum/#!forum/pynastran-discuss'

is_release = True  ## True=turns on skipping of tables that aren't supported

# separarate table_passer crasher; previously controlled by is_release
# TODO: currently no warn option
stop_on_op2_table_passer = False

# stop takes priority; previously controlled by is_release
# a missed table refers to the not_implemented_error_or_skip method
stop_on_op2_missed_table = False

## True=turns on skipping of tables that aren't supported
# a missed table refers to the not_implemented_error_or_skip method
warn_on_op2_missed_table = True
