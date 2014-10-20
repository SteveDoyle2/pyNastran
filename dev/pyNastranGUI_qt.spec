# -*- PyInstaller input file -*-
# -*- mode: python           -*-

#pyInstaller_path = 'C:\\Users\\steve\\Desktop\\pyInstaller'
import pyinstaller
pyInstaller_path = os.path.dirname(pyinstaller.__file__)
print "pyInstaller_path = %r" % pyInstaller_path


import os
import sys
import pyNastran

pkg_path = pyNastran.__path__[0]

#a1 = os.path.join(pkg_path, 'bdf','bdf.py')
#a2 = os.path.join(pkg_path, 'op2','op2.py')
#a3 = os.path.join(pkg_path, 'f06','f06.py')
#a4 = os.path.join(pkg_path, 'gui','gui_wx.py')
a4 = os.path.join(pkg_path, 'gui','gui_qt.py')

#analyze_files = [a1,a2,a3,a4]
analyze_files = [a4]

#print "pynastran = ",pyNastran.__file__
icon_path = os.path.join(pkg_path, 'gui', 'icons')

icon_names = os.listdir(icon_path)
icons = []
for icon_name in icon_names:
    i1 = os.path.join('icons', icon_name)
    i2 = os.path.join(icon_path, icon_name)
    icon = (i1, i2, 'DATA')
    icons.append(icon)
#i1  = ('icons/twireframe.png',os.path.join(icon_path,'twireframe.png'),'DATA')
#i2  = ('icons/tbdf.png',      os.path.join(icon_path,'tbdf.png'),      'DATA')
#i3  = ('icons/tcamera.png',   os.path.join(icon_path,'tcamera.png'),   'DATA')
#i4  = ('icons/texit.png',     os.path.join(icon_path,'texit.png'),     'DATA')
#i5  = ('icons/top2.png',      os.path.join(icon_path,'top2.png'),      'DATA')
#i6  = ('icons/topen.png',     os.path.join(icon_path,'topen.png'),     'DATA')
#i7  = ('icons/tredo.png',     os.path.join(icon_path,'tredo.png'),     'DATA')
#i8  = ('icons/tsave.png',     os.path.join(icon_path,'tsave.png'),     'DATA')
#i9  = ('icons/tsaveAs.png',   os.path.join(icon_path,'tsaveAs.png'),   'DATA')
#i10 = ('icons/tsolid.png',    os.path.join(iconPath,'tsolid.png'),    'DATA')
#i11 = ('icons/tundo.png',     os.path.join(iconPath,'tundo.png'),     'DATA')

#icon_main = os.path.join(iconPath,'guiDemo_128.ico')
icon_main = 'logo.ico' # pyInstaller doesnt like paths in the ico
                       # it also won't handle png files
assert os.path.exists(icon_main), '%s doesnt exist' % icon_main
#icons = [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]

#-------------------------------------------------------------------------
# main code

binaries = []
if sys.platform == 'win32':
    binaries = [('msvcp100.dll', 'C:\\Windows\\System32\\msvcp100.dll', 'BINARY'),
                ('msvcr100.dll', 'C:\\Windows\\System32\\msvcr100.dll', 'BINARY')]

python_path = 'C:\\Python27_x86'
a = Analysis(analyze_files,
             pathex=[pyInstaller_path,
                     python_path,
                     os.path.join(python_path, 'Lib'),
                     os.path.join(python_path, 'Lib', 'site-packages'),

                     pkg_path,
                     os.path.join(pkg_path, 'bdf'),
                     os.path.join(pkg_path, 'op2'),
                     os.path.join(pkg_path, 'f06'),
                     os.path.join(pkg_path, 'gui'),
                    ],
             hiddenimports=[],
             hookspath=None)
pyz = PYZ(a.pure)

#print "help(EXE) = \n",help(EXE)
exe = EXE(pyz,
          a.scripts,
          a.binaries + binaries + icons,
          a.zipfiles,
          a.datas,
          #exclude_binaries=True,
          name=os.path.join('build\\pyi.win32\\pyNastranGUI', 'pyNastranGUI.exe'),
          debug=True,
          strip=None,
          #upx=True,
          icon=icon_main,
          console=True )

#print '*'*80
#print "help(COLLECT) = \n",help(COLLECT)
#coll = COLLECT(exe,
#               a.binaries + binaries + icons,
#               a.zipfiles,
#               a.datas,
#               exclude_binaries=1,
#               #icon=icon_main,
#               strip=None,
#               upx=True,
#               name=os.path.join('dist', 'pyNastranGUI'))

