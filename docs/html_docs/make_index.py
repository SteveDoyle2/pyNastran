"""
Regenerates the *.rst files based on the code.
There is some manual work that has to happen afterwards.
"""
from __future__ import print_function, unicode_literals
import os
from io import open
from cpylog import get_logger

log_debug = get_logger(log=None, debug=True, encoding='utf-8')
log_warning = get_logger(log=None, debug=None, encoding='utf-8')


IGNORE_DIRS = ['src', 'dmap', 'solver', '__pycache__',
               'op4_old', 'calculix', 'bars', 'case_control',
               'pch', 'old', 'solver', 'test', 'dev', 'bkp', 'bdf_vectorized']
MODS_SKIP = ['spike', 'shell_backup']

SKIP_DIRECTORIES = ['.svn', '.idea', '.settings', '.git', 'test', 'bkp', '__pycache__', 'dev',
                    'htmlcov', 'vtk_examples', 'SnakeRiverCanyon', 'M100', 'SWB']
SKIP_FILE_SUFFIX = [
    '.pyc', '.pyx', # python
    '.bdf', '.op2', '.f06', '.op4', '.dat', '.inp', # nastran
    '.err', '.log', '.rej', '.db', '.db.jou', '.ses', '.ses.01', # patran
    '.pptx',
    '.png', '.gif', # pictures
    '.txt', '.csv', '.out', '.coverage', '.whl', # generic
    '.mapbc', '.front', '.flo', 'cogsg', '.bc', '.d3m', '.inpt', '.nml', # usm3d/fun3d
    '.ele', '.node', '.smesh', '.off',
    '.mk5', '.wgs', '.stl', '.fgrid', '.su2', '.obj', # other formats
    '.tri', '.cntl', '.c3d',  # cart3d
    '.surf', '.tags', '.ugrid', '.bedge', # aflr
    '.plt', # tecplot
    '.p3d',
    '.tex', '.bib', # latex
    ]
MAKE_FILES = True

def get_folders_files(dirname, skip_file_suffix=None, skip_directories=None):
    """
    Return list of directories and files in a given tree path. By default discards:

    * directories ".svn", ".idea", ".settings"

    * files that ends with ".pyc", .pyx", ".bdf"
    """
    if skip_directories is None:
        skip_directories = SKIP_DIRECTORIES
    if skip_file_suffix is None:
        skip_file_suffix = tuple(SKIP_FILE_SUFFIX)

    dirname = os.path.join(dirname)
    files = []
    folders = []
    for root, dirs, fil in os.walk(dirname):
        folders.append(root)
        files += [os.path.join(root, i) for i in fil
                  if not i.endswith(skip_file_suffix)]
        dirs[:] = [d for d in dirs if not d in skip_directories]
        #if len(dirs):
            #print('root = %s' % root)
            #print(dirs)
            #print('------------------')
    return folders, files

def get_classes(py_filenames, log):
    d = {}
    for py_filename in py_filenames:
        #log.info("py_filename =", py_filename)
        with open(py_filename, 'r') as py_file:
            lines = py_file.readlines()

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
                dirname = os.path.dirname(py_filename)
                d[dirname] = name

def get_doc(filenames, log):
    """regenerate the *.rst files"""
    dirs = {}
    for filename in filenames:
        if '.py' not in filename[-3:]:
            log.info('skipping %s' % filename)
            continue
        log.debug('trying %s' % filename)

        dirname = os.path.dirname(filename)
        dirname = os.path.join(dirname)
        if dirname == r'../pyNastran':
            continue

        if dirname not in dirs:
            dirs[dirname] = []
        dirs[dirname].append(filename)

    i = 0
    for dirname, files in sorted(dirs.items()):
        # temporarily suppressing log messages for classes that aren't the bdf
        #if 'cards' not in dirname: # 'bdf' not in dirname and
        if 'cards' not in dirname:
            log = log_warning
        else:
            log = log_debug

        if len(files) == 0:
            continue
        log.info("******************************************************************************")
        basepath = os.path.basename(dirname)
        if basepath in IGNORE_DIRS:
            log.debug('skipping ignore_dirs %r' % dirname)
            continue
        if 'dev_' in basepath:
            log.info('skipping dev_* %r' % dirname)
            continue

        #absname = os.path.join(dirname, filename)
        dirname_base = dirname.strip('./')
        rst_dirname = os.path.dirname(dirname_base)

        dot_dirname = dirname_base.replace('\\', '.')

        # ../../pyNastran\applications\aero_panel_buckling
        log.info("dirname = %s" % dirname)

        # pyNastran\applications\aero_panel_buckling
        log.info("dirname_base = %s" % dirname_base)

        # pyNastran.applications.aero_panel_buckling
        log.info('dot_dirname = %s' % dot_dirname)

        # pyNastran\applications
        #log.info('rst_dirname = %s' % rst_dirname)

        rst_sline = rst_dirname.split('\\', 2)

        log.info('rst_sline = %s' % rst_sline)
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
            rst_dirname2 = os.path.join(rst_sline[1])
            #rst_dirname2 = os.path.join(*rst_sline[1:])
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
        #print(msg)

        if len(files) == 0:
            continue

        if not MAKE_FILES:
            continue

        sub_dirs = get_sub_dirs(dirname, log)
        write_rst_file(rst_path, files, dirname, sub_dirs, msg, log)

        #if i == 5:
            #break
        i += 1

def write_rst_file(rst_filename, files, unused_dirname, sub_dirs, msg, log):
    log.info('working on rst_filename = %s' % rst_filename)
    with open(rst_filename, 'w') as rst_file:
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

            if modname in MODS_SKIP:
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
                log.debug("    sub_rst_filename = %s" % sub_rst_filename)
                msg += '    %s\n' % sub_rst_filename
            msg += '\n'
            rst_file.write(msg)

def get_sub_dirs(dirname, log):
    maybe_dirs = os.listdir(dirname)
    sub_dirs = []
    #print("dirname_check = %s" % dirname)

    for idir in maybe_dirs:
        if idir in IGNORE_DIRS:
            continue

        dir2 = os.path.join(dirname, idir)
        if os.path.isdir(dir2):
            init_file = os.path.join(idir, '__init__.py')
            log.warning(init_file)
            if idir not in IGNORE_DIRS and os.path.exists(init_file):
                #print("adding %s" % dir2)
                sub_dirs.append(dir2)
    return sub_dirs

def run():
    """regenerates the *.rst files"""
    log = get_logger(log=None, debug=True, encoding='utf-8')
    unused_folders, files = get_folders_files('../../pyNastran')
    #c = get_classes(files)
    get_doc(files, log)


if __name__ == '__main__':  # pragma: no cover
    run()
