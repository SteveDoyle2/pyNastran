"""
First you must:
 - pip install nbconvert
 - install pandoc: https://pandoc.org/installing.html
 - make sure 'which pandoc' finds it; Powershell doesn't work for me; cmd does
"""
from __future__ import print_function
import os

import pyNastran
from pyNastran.utils import print_bad_path
pkg_path = pyNastran.__path__[0]


def create_rst_from_ipython_notebooks():
    """converts a Jupyter notebook"""
    #curdir = os.getcwd()
    notebook_dir = os.path.join(pkg_path, '..', 'docs', 'quick_start', 'demo')
    pydocs_dir = os.path.join(pkg_path, '..', 'docs', 'html_docs', 'quick_start')

    assert os.path.exists(pydocs_dir), print_bad_path(pydocs_dir)
    assert os.path.exists(notebook_dir), print_bad_path(notebook_dir)
    os.chdir(notebook_dir)
    for fname in os.listdir(notebook_dir):
        fnamei = os.path.basename(fname)
        base = os.path.splitext(fnamei)[0]
        rst_name = base + '.rst'

        if fnamei.startswith('.'):
            continue
        if not fnamei.endswith('.ipynb'):
            continue
        assert  os.path.exists(fnamei)
        cmd = 'jupyter nbconvert --to rst %s' % fnamei
        os.system(cmd)
        if not os.path.exists(rst_name):
            print('%s was not made...' % rst_name)
            continue

        moved_rst_name = os.path.join(pydocs_dir, rst_name)
        try:
            if os.path.exists(moved_rst_name):
                os.remove(moved_rst_name)
            os.rename(rst_name, moved_rst_name)
        except:
            pass

if __name__ == '__main__':  # pragma: no cover
    create_rst_from_ipython_notebooks()
