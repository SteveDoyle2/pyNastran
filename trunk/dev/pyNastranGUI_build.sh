#!/bin/bash
#python pyinstaller-pyinstaller-67b940c/pyinstaller.py ../pyNastran/pyNastran/gui/gui.py
rm -rf build dist
python pyinstaller-pyinstaller-67b940c/pyinstaller.py gui.spec
#dist/pywin27/gui.exe
#dist/gui/gui.exe

