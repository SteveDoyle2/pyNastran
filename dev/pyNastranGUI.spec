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
a4 = os.path.join(pkg_path, 'gui','gui.py')

#analyze_files = [a1,a2,a3,a4]
analyze_files = [a4]

icon_path = os.path.join(pkg_path, 'gui', 'icons')

icon_names = os.listdir(icon_path)
icons = []
for icon_name in icon_names:
    i1 = os.path.join('icons', icon_name)   # where to put the icon in the exe
    i2 = os.path.join(icon_path, icon_name) # where the icon is now
    icon = (i1, i2, 'DATA')
    icons.append(icon)

#icon_main = os.path.join(iconPath,'guiDemo_128.ico')
icon_main = 'logo.ico' # pyInstaller doesn't like paths in the ico
                       # it also won't handle png files
assert os.path.exists(icon_main), '%s doesnt exist' % icon_main

#-------------------------------------------------------------------------
# main code

binaries = []
if sys.platform == 'win32':
    binaries = [('msvcp100.dll', 'C:\\Windows\\System32\\msvcp100.dll', 'BINARY'),
                ('msvcr100.dll', 'C:\\Windows\\System32\\msvcr100.dll', 'BINARY')]

#python_path = 'C:\\Python27_x86'
python_path = 'C:\\Anaconda'
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
             hiddenimports=[
                'vtk.vtkCommonPythonSIP', 'vtk.vtkFilteringPythonSIP',
                'PyQt4.QtOpenGL', 'vtk.vtkRenderingPythonSIP'],
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
          debug=False,
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

