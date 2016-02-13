import sys
from PyInstaller.__main__ import run, run_build
spec_file = 'pyNastranGUI.spec'
#sys.argv = [None, 'pyNastranGUI.spec']
#run(pyi_args=)

pyi_config = None
opts = []
run_build(opts, spec_file, pyi_config)
