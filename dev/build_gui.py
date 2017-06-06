"""
don't even try to run this...

run:
 - python <path_to_pyinstaller>/pyinstaller.py pyNastranGUI.spec
"""
import os
import sys
from PyInstaller.__main__ import run, run_build
spec_file = 'pyNastranGUI.spec'
#sys.argv = [None, 'pyNastranGUI.spec']
#run(pyi_args=)

pyi_config = {
    'hasUPX' : True,
    'cachedir' : '',
}
opts = []
#run_build(opts, spec_file, pyi_config)
#run_build(opts, spec_file, noconfirm=True)
run_build(pyi_config, spec_file, noconfirm=True,
          distpath=os.path.join(os.curdir, 'dist'),
          workpath=os.path.join(os.curdir, 'build'),
          #clean_build=True,  # requires CACHE_DIR
          clean_build=False,
          #cachedir='',
)
