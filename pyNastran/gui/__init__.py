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

if is_pynastrangui_exe:
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
if IS_WINDOWS:
    font_dirname = r'C:\Windows\Fonts'
    if os.path.exists(font_dirname):
        font_file = os.path.join(font_dirname, 'arial.ttf')
        assert os.path.exists(font_file), os.listdir(font_dirname)
elif IS_LINUX:
    font_dirname = '/etc/fonts/fonts.conf'
elif IS_MAC:
    font_dirname = '/Library/Fonts'

font_file = ''
if font_dirname and os.path.exists(font_dirname):
    font_file = os.path.join(font_dirname, 'arial.ttf')
    assert os.path.exists(font_file), os.listdir(font_dirname)
del font_dirname, IS_WINDOWS, IS_LINUX, IS_MAC
