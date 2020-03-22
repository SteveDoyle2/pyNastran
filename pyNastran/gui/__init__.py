"""defines some common paths to the scripts and icons"""
import os
import sys
import pyNastran
from pyNastran import is_pynastrangui_exe

IS_DEV = (
    'TRAVIS' in os.environ or
    'APPVEYOR' in os.environ or
    'READTHEDOCS' in os.environ
)
#IS_TRAVIS = 'TRAVIS' in os.environ
#IS_RTD = 'READTHEDOCS' in os.environ

if is_pynastrangui_exe:  # pragma: no cover
    PKG_PATH = sys._MEIPASS #@UndefinedVariable
    SCRIPT_PATH = os.path.join(PKG_PATH, 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'icons')
else:
    PKG_PATH = pyNastran.__path__[0]
    SCRIPT_PATH = os.path.join(PKG_PATH, 'gui', 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'gui', 'icons')

IS_WINDOWS = 'nt' in os.name
IS_LINUX = 'posix' in os.name
IS_MAC = 'darwin' in os.name

font_dirname = ''
font_file = ''
bi_font_file = ''
if IS_WINDOWS:
    font_dirname = r'C:\Windows\Fonts'
    if os.path.exists(font_dirname):
        font_file = os.path.join(font_dirname, 'arial.ttf')
        assert os.path.exists(font_file), os.listdir(font_dirname)

        bi_font_file = os.path.join(font_dirname, 'arialbi.ttf')
        assert os.path.exists(bi_font_file), os.listdir(font_dirname)

elif IS_LINUX:
    # /etc/fonts/fonts.conf is an xml file that points to one of many font directories
    #
    #font_dirname = '/etc/fonts/fonts.conf'
    #/usr/local/share/fonts
    #/usr/share/fonts/

    #  we want a true-type font, so...
    #/usr/share/fonts/truetype/

    # only 1 font class for my system...
    #/usr/share/fonts/truetype/dejavu/
    font_dirname = '/usr/share/fonts/truetype/dejavu/'

    # let's go with a Sans-Serif font that's not bold
    # DejaVuSans.ttf            DejaVuSerif.ttf          DejaVuSansMono.ttf
    # DejaVuSans-Bold.ttf       DejaVuSerif-Bold.ttf     DejaVuSansMono-Bold.ttf
    #
    if os.path.exists(font_dirname):
        font_file = os.path.join(font_dirname, 'DejaVuSans.ttf')
        # TODO: where is the italics font?
        bi_font_file = os.path.join(font_dirname, 'DejaVuSans-Bold.ttf')

elif IS_MAC: # pragma: no cover
    # TODO: pulled off a screenshot...not tested...
    font_dirname = '/Library/Fonts'
    assert os.path.exists(font_dirname), font_dirname
    font_file = os.path.join(font_dirname, 'Arial.ttf')
    bi_font_file = os.path.join(font_dirname, 'Arialbi.ttf')
    assert os.path.exists(font_file), os.listdir(font_dirname)
    assert os.path.exists(bi_font_file), os.listdir(font_dirname)

if font_file and not os.path.exists(font_file):
    print('cant find %s' % font_file)
    #assert os.path.exists(font_file), os.listdir(font_dirname)
    font_file = ''

if bi_font_file and not os.path.exists(bi_font_file):
    print('cant find %s' % bi_font_file)
    #assert os.path.exists(font_file), os.listdir(font_dirname)
    bi_font_file = font_file

del font_dirname # , IS_WINDOWS, IS_LINUX, IS_MAC
