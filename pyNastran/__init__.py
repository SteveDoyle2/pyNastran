import subprocess
def get_git_revision_short_hash():
    try:
        ghash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        ghash = ghash.decode('utf-8').rstrip()
    except:
        ghash = 'no_checksum_error'
    return 'dev-%s' % ghash

#revision = get_git_revision_short_hash()
__author__  = 'Steven Doyle, Al Danials, Marcin Gasiorek'
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__     = 'LGPLv3'
__copyright__   = 'Copyright %s; 2011-2015' % __license__
__pyqt_copyright__ = 'Copyright GPLv3; 2011-2015'
__releaseDate__ = '2015/4/xx'
__releaseDate2__ = 'APRIL xx, 2015'
#__version__     = '0.7.2_%s' % revision
__version__     = '0.7.2'
__website__     = 'https://github.com/SteveDoyle2/pyNastran'

is_release = True  ## turns on skipping of tables that aren't supported
