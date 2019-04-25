import os
import sys
import pyNastran

IS_DEV = (
    'TRAVIS' in os.environ or
    'APPVEYOR' in os.environ or
    'READTHEDOCS' in os.environ
)
#IS_TRAVIS = 'TRAVIS' in os.environ
#IS_RTD = 'READTHEDOCS' in os.environ

try:
    PKG_PATH = sys._MEIPASS #@UndefinedVariable
    SCRIPT_PATH = os.path.join(PKG_PATH, 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'icons')
except:
    PKG_PATH = pyNastran.__path__[0]
    SCRIPT_PATH = os.path.join(PKG_PATH, 'gui', 'scripts')
    ICON_PATH = os.path.join(PKG_PATH, 'gui', 'icons')
