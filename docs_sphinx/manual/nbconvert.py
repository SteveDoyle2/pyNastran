from __future__ import print_function
import os
import shutil

import pyNastran
from pyNastran.utils import print_bad_path
pkg_path = pyNastran.__path__[0]

def create_rst_from_python_files():
    curdir = os.getcwd()
    notebook_dir = os.path.join(pkg_path, '..', 'quick_start', 'demo')
    pydocs_dir = os.path.join(pkg_path, '..', 'docs_sphinx', 'manual', 'py_docs')

    assert os.path.exists(pydocs_dir), print_bad_path(quick_start_pydocs_dir)
    assert os.path.exists(notebook_dir), print_bad_path(notebook_dir)
    os.chdir(notebook_dir)
    for fname in os.listdir(notebook_dir):
        #fname1 = os.path.join(notebook_dir, fname)
        fname2 = os.path.join(pydocs_dir, fname)
        if '.ipynb' in fname:
            os.system('ipython nbconvert --to rst %s' % fname)
        assert os.path.exists(fname2), fname2
        #print(fname1)
        #print(fname2)

        shutil.copyfile(fname1, fname2)
    return
    fnames = [
        os.path.join('py_docs', fname) for fname in os.listdir('py_docs')
        if '__init__' not in fname
        and fname.endswith('.py')
        and not fname.endswith('.pyc')
        and not fname.endswith('.html')
        and not fname.endswith('.txt')
        #and not fname.endswith('.rst')
    ]

    print(fnames)
    for py_filename in fnames:
        root = os.path.splitext(py_filename)[0]
        base = os.path.splitext(os.path.basename(py_filename))[0]

        rst_filename = base + '.rst'
        html_filename = root + '.html'
        print(rst_filename)
        #print(py_filename)
        create_rst_from_python(py_filename, rst_filename)
        create_html_from_rst(rst_filename, html_filename)

    print(os.path.abspath(py_filename))
    for py_filename in fnames:
        root = os.path.splitext(py_filename)[0]
        base = os.path.splitext(os.path.basename(py_filename))[0]
        rst_filename = base + '.rst'
        print('ran %s; saved %s' % (py_filename, rst_filename))
