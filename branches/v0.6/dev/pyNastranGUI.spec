# -*- PyInstaller input file -*-
# -*- mode: python           -*-
import os
import sys
import pyNastran

pkgPath = pyNastran.__path__[0]

a1 = os.path.join(pkgPath,'bdf','bdf.py')
a2 = os.path.join(pkgPath,'op2','op2.py')
a3 = os.path.join(pkgPath,'f06','f06.py')
a4 = os.path.join(pkgPath,'gui','gui.py')

#analyzeFiles = [a1,a2,a3,a4]
analyzeFiles = [a4]

#print "pynastran = ",pyNastran.__file__
iconPath = os.path.join(pkgPath,'gui','icons')


i1  = ('icons/twireframe.png',os.path.join(iconPath,'twireframe.png'),'DATA')
i2  = ('icons/tbdf.png',      os.path.join(iconPath,'tbdf.png'),      'DATA')
i3  = ('icons/tcamera.png',   os.path.join(iconPath,'tcamera.png'),   'DATA')
i4  = ('icons/texit.png',     os.path.join(iconPath,'texit.png'),     'DATA')
i5  = ('icons/top2.png',      os.path.join(iconPath,'top2.png'),      'DATA')
i6  = ('icons/topen.png',     os.path.join(iconPath,'topen.png'),     'DATA')
i7  = ('icons/tredo.png',     os.path.join(iconPath,'tredo.png'),     'DATA')
i8  = ('icons/tsave.png',     os.path.join(iconPath,'tsave.png'),     'DATA')
i9  = ('icons/tsaveAs.png',   os.path.join(iconPath,'tsaveAs.png'),   'DATA')
i10 = ('icons/tsolid.png',    os.path.join(iconPath,'tsolid.png'),    'DATA')
i11 = ('icons/tundo.png',     os.path.join(iconPath,'tundo.png'),     'DATA')

#iconMain = os.path.join(iconPath,'guiDemo_128.ico')
iconMain = 'guiDemo.ico' # pyInstaller doesnt like paths in the ico

icons = [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]

#-------------------------------------------------------------------------
# main code

a = Analysis(analyzeFiles,
             pathex=['C:\\Users\\steve\\Desktop\\pyInstaller',
                     'C:\\Python27\\Lib\\site-packages',
                     'C:\\Python27\\Lib',
                     'C:\\Python27',
                     pkgPath,
                     os.path.join(pkgPath,'bdf'),
                     os.path.join(pkgPath,'op2'),
                     os.path.join(pkgPath,'f06'),
                     os.path.join(pkgPath,'gui'),
                    ],
             hiddenimports=[],
             hookspath=None)
pyz = PYZ(a.pure)

#print "help(EXE) = \n",help(EXE)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name=os.path.join('build\\pyi.win32\\pyNastranGUI', 'pyNastranGUI.exe'),
          debug=False,
          strip=None,
          upx=True,
          icon=iconMain,
          console=True )

#print '*'*80
#print "help(COLLECT) = \n",help(COLLECT)
coll = COLLECT(exe,
               a.binaries+icons,
               a.zipfiles,
               a.datas,
               exclude_binaries=1,
               icon=iconMain,
               strip=None,
               upx=True,
               name=os.path.join('dist', 'pyNastranGUI'))

