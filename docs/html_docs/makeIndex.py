"""
regenerates the *.rst files based on the code
"""
from __future__ import print_function
import os
from six import iteritems
import pyNastran
from pyNastran.utils.log import get_logger2
log = get_logger2(log=None, debug=True, encoding='utf-8')


pkg_path = pyNastran.__path__[0]

ignore_dirs = ['src', 'dmap', 'solver', '__pycache__',
               'op4_old', 'calculix', 'bars', 'case_control',
               'pch', 'old', 'solver', 'test']
mods_skip = ['spike', 'shell_backup']

make_files = True

def get_folders_files(dirname, skip_file_suffix=('.pyc', '.pyx', '.bdf', '.op2', '.f06', '.op4'),
                      skip_directories=('.svn', '.idea', '.settings', '.git')):
    """
    Return list of directories and files in a given tree path. By default discards:

    * directories ".svn", ".idea", ".settings"

    * files that ends with ".pyc", .pyx", ".bdf"
    """
    dirname = os.path.join(dirname)
    files = []
    folders = []
    for root, dirs, fil in os.walk(dirname):
        folders.append(root)
        files += [os.path.join(root, i) for i in fil
                  if not i.endswith(skip_file_suffix)]
        dirs[:] = [d for d in dirs if not d in skip_directories]
    return (folders, files)

def get_classes(filenames):
    d = {}
    for filename in filenames:
        #log.info("filename =", filename)
        with open(filename, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if 'class' in line:
                log.info(line.rstrip())
                sline = line.split('(')
                if len(sline) == 2:
                    if 'raise' in sline[0]:
                        continue
                    pre, post = sline
                else:
                    log.info("*******")
                    continue

                pre, post = pre.split('class')
                name = post.strip()
                log.info("name = %s" % name)
                dirname = os.path.dirname(filename)
                d[dirname] = name

def get_doc(filenames):
    """regenerate the *.rst files"""
    dirs = {}
    for filename in filenames:
        log.debug('trying %s' % filename)
        if '.py' not in filename[-3:]:
            continue

        dirname = os.path.dirname(filename)
        dirname = os.path.join(dirname)
        if dirname == r'../pyNastran':
            continue

        if dirname not in dirs:
            dirs[dirname] = []
        dirs[dirname].append(filename)

    i = 0
    for dirname, files in sorted(iteritems(dirs)):
        if len(files) == 0:
            continue
        log.info("******************************************************************************")
        basepath = os.path.basename(dirname)
        if basepath in ignore_dirs:
            log.debug('skipping ignore_dirs %r' % dirname)
            continue
        if 'dev_' in basepath:
            log.info('skipping dev_* %r' % dirname)
            continue

        #absname = os.path.join(dirname, filename)
        dirnamebase = dirname.strip('./')
        rst_dirname = os.path.dirname(dirnamebase)

        dot_dirname = dirnamebase.replace('\\', '.')

        # ../../pyNastran\applications\aero_panel_buckling
        log.info("dirname = %s" % dirname)

        # pyNastran\applications\aero_panel_buckling
        log.info("dirnamebase = %s" % dirnamebase)

        # pyNastran.applications.aero_panel_buckling
        log.info('dot_dirname = %s' % dot_dirname)

        # pyNastran\applications
        #log.info('rst_dirname = %s' % rst_dirname)

        rst_sline = rst_dirname.split('\\', 2)

        #log.info('rst_sline = ', rst_sline)
        if len(rst_sline) == 1 and rst_sline[0] in ('', 'pyNastran'):
            # what's the quote directory?
            # pyNastran
            rst_dirname2 = ''
        elif len(rst_sline) == 2:
            # pyNastran\applications
            # pyNastran\bdf -> pyNastran
            rst_dirname2 = rst_sline[1]
        elif len(rst_sline) == 2 or len(rst_sline) == 3:
            # pyNastran\bdf\cards -> bdf
            # pyNastran\bdf\cards\elements -> bdf
            rst_dirname2 = rst_sline[1]
            #rst_sline2 =  = rst_sline[0].split('\\', 1)
            #raise NotImplementedError()
        else:
            raise NotImplementedError()
            #rst_dirname2 = ''
        #if rst_dirname2:
        log.debug('rst_dirname2 = %r' % rst_dirname2)

        new_dirname = os.path.join('index', rst_dirname2)
        if not os.path.exists(new_dirname):
            os.makedirs(new_dirname)

        # TODO: update this...
        rst_dirname2 = ''
        rst_path = os.path.join('index', rst_dirname2, '%s.rst' % dot_dirname)

        msg = '%s Package\n' % basepath  # loads Package
        msg += '=' * (len(msg)-1) + '\n'
        msg += '\n'

        log.info("----------------------")
        #log.info('absname2 = ' % absname2)
        log.debug('rst_path = %s' % rst_path)
        log.debug("basepath = %s" % basepath)

        if len(files) == 0:
            continue

        if not make_files:
            continue

        sub_dirs = get_sub_dirs(dirname)
        write_rst_file(rst_path, files, dirname, sub_dirs, msg)

        #if i == 5:
            #break
        i += 1

def write_rst_file(rst_filename, files, dirname, sub_dirs, msg):
    log.info('working on rst_filename = %s' % rst_filename)
    with open(rst_filename, 'wb') as rst_file:
        rst_file.write(msg)
        for fname in files:
            if os.path.getsize(fname) == 0:
                continue
            assert os.path.exists(fname), fname

            py_filename_base = os.path.basename(fname)
            modname = py_filename_base[:-3]

            # path to the module relative to the root of the pyNastran directory
            # pyNastran\utils\numpy_utils
            pynastran_path = fname[:-3].strip('/.')

            # corresponding "import path"
            # pyNastran.utils.numpy_utils
            automodname = pynastran_path.replace('\\', '.')

            #log.info("fname = ", fname)
            #log.info("py_filename_base = %s" % py_filename_base)
            log.debug("  pynastran_path = %s.py" % pynastran_path)
            #log.info("automodname    = %s" % automodname)

            if modname in mods_skip:
                continue

            msg = ':mod:`%s` Module\n' % modname  # `loads\  # fname2
            msg += '-' * (len(msg)-1) + '\n'
            msg += '\n'
            msg += '.. inheritance-diagram:: %s\n' % automodname
            msg += '\n'
            # pyNastran.bdf.bdfInterface.bdf_card, pyNastran.bdf.cards.loads.static_loads
            msg += '.. automodule:: %s\n' % automodname
            msg += '    :members:\n'
            #msg += '    :private-members:\n'
            msg += '    :undoc-members:\n'
            msg += '    :show-inheritance:\n'

            #msg += '    :autofunction:\n'
            #msg += '    :autodata:\n'
            #msg += '    :automethod:\n'
            #msg += '    :autoattribute:\n'

            msg += '\n'
            rst_file.write(msg)

        if sub_dirs:
            log.debug('  sub_dirs = %s' % sub_dirs)
            msg += 'Subpackages\n'
            msg += '-----------\n'
            msg += '\n'
            msg += '.. toctree::\n'
            msg += '\n'
            for diri in sub_dirs:
                sub_rst_filename = diri.strip('/.').replace('\\', '.')
                log.debug("    sub_rst_filename =", sub_rst_filename)
                msg += '    %s\n' % sub_rst_filename
            msg += '\n'
            rst_file.write(msg)

def get_sub_dirs(dirname):
    maybe_dirs = os.listdir(dirname)
    sub_dirs = []
    #print("dirname_check = %s" % dirname)

    for idir in maybe_dirs:
        if idir in ignore_dirs:
            continue

        dir2 = os.path.join(dirname, idir)
        if os.path.isdir(dir2):
            init_file = os.path.join(idir, '__init__.py')
            log.warning(init_file)
            if idir not in ignore_dirs and os.path.exists(init_file):
                #print("adding %s" % dir2)
                sub_dirs.append(dir2)
    return sub_dirs

def run():
    """regenerates the *.rst files"""
    folders, files = get_folders_files('../../pyNastran')
    #c = get_classes(files)
    get_doc(files)


if __name__ == '__main__':  # pragma: no cover
    run()
