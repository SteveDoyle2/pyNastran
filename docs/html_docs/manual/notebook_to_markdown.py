from __future__ import print_function
import os
import shutil

import pyNastran
from pyNastran.utils import print_bad_path
pkg_path = pyNastran.__path__[0]


def create_rst_from_ipython_notebooks():
    #curdir = os.getcwd()
    notebook_dir = os.path.join(pkg_path, '..', 'docs', 'quick_start', 'demo')
    pydocs_dir = os.path.join(pkg_path, '..', 'docs', 'html_docs', 'quick_start')

    assert os.path.exists(pydocs_dir), print_bad_path(quick_start_pydocs_dir)
    assert os.path.exists(notebook_dir), print_bad_path(notebook_dir)
    os.chdir(notebook_dir)
    for fname in os.listdir(notebook_dir):
        fnamei = os.path.basename(fname)
        base = os.path.splitext(fnamei)[0]
        fname2 = base + '.rst'

        if fnamei.startswith('.'):
            continue
        if not fnamei.endswith('.ipynb'):
            continue
        os.system('ipython nbconvert --to rst %s' % fname)
        if not os.path.exists(fname2):
            print('%s was not made...' % fname2)
            continue

        moved_fname2 = os.path.join(pydocs_dir, fname2)
        try:
            if os.path.exists(moved_fname2):
                os.remove(moved_fname2)
            os.rename(fname2, moved_fname2)
        except:
            pass

