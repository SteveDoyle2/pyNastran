from __future__ import print_function, unicode_literals
import os
from io import open
from cpylog import get_logger
import shutil

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
    for root, dirs, filenames in os.walk(dirname):
        folders.append(root)
        for filename in filenames:
            if filename.endswith(skip_file_suffix):
                continue
            if 'test_' in os.path.basename(filename):
                continue
            files.append(os.path.join(root, filename))
        #files += [os.path.join(root, filename) for filename in filenames
                  #if not filename.endswith(skip_file_suffix)]
        dirs[:] = [d for d in dirs if not d in skip_directories]
        #if len(dirs):
            #print('root = %s' % root)
            #print(dirs)
            #print('------------------')
    return folders, files


def get_classes_functions_in_file(py_filename):
    with open(py_filename, 'r', encoding='utf8') as f:
        lines = f.readlines()

    function_list = []
    class_list = []
    for line in lines:
        line = line.split('#')[0].rstrip()
        if line.startswith('class '):
            # class ASDF(object):
            class_name = line.split('(')[0].split(' ')[1]
            is_object = False
            if '(object):' in line:
                is_object = True
            class_list.append((class_name, is_object))
        elif line.startswith('def '):
            function_name = line.split('(')[0].split(' ')[1]
            if function_name.startswith('_'):
                continue
            function_list.append(function_name)


    #for class_name in class_list:
        #print(class_name)
    return class_list, function_list

def get_pyfilenames():
    folders, filenames = get_folders_files('../../pyNastran')

    filenames_classes = []
    for py_filename in filenames:
        py_filename2, dot_path = get_location_filename_for_pyfilename(py_filename)

        class_names, function_names = get_classes_functions_in_file(py_filename)
        #for class_name, is_object in class_names:
            #print('  %s (class)' % class_name)
        #for function_name in function_names:
            #print('  %s (function)' % function_name)
        filenames_classes.append((py_filename, py_filename2, dot_path, class_names))
    return filenames_classes

def get_location_filename_for_pyfilename(py_filename):
    """../../pyNastran/utils/nastran_utils.py -> pyNastran/utils/nastran_utils.py"""
    path = py_filename.lstrip('../\\')
    no_py = os.path.splitext(path)[0]
    dot_path = no_py.replace('\\', '.').replace('/', '.')

    #print(dot_path)
    return path, dot_path

def filenames_to_rsts(filenames_classes, make_rsts=False):
    for py_filename, py_filename2, dot_path, class_names in filenames_classes:
        if not class_names:
            continue

        base_folder = os.path.dirname(py_filename2)
        #print('%-20s %s %s' % (base_folder[:20], py_filename2, dot_path))
        folder = os.path.join('rsts', base_folder)
        if 'cards' in folder:
            while not folder.endswith('cards'):
                folder = os.path.dirname(folder)

        if not os.path.exists(folder):
            os.makedirs(folder)

        rst_filename = os.path.join(folder, 'index.rst')
        mode = 'w'
        rst_lines = '.. toctree::\n\n'
        if os.path.exists(rst_filename):
            rst_lines = ''
            mode = 'a'

        for class_name, is_object in class_names:
            create_rst_file_for_class(folder, dot_path, class_name, is_object)
            print('  %s' % str(class_name))

            #pyNastran.bdf.cards.aset
            rst_lines += '    %s.%s\n' % (dot_path, class_name)

        #print(rst_lines)
        with open(rst_filename, mode) as rst_file:
            rst_file.write(rst_lines)

def create_rst_file_for_class(folder, dot_path, class_name, is_object):
    split_path = dot_path.split('.')
    split_path[-1] += '.rst'
    #rst_filename = os.path.join(*split_path)
    dot_class_path = '%s.%s.rst' % (dot_path, class_name)
    rst_filename = os.path.join(folder, dot_class_path)
    #dirname = os.path.dirname(rst_filename)
    #if not os.path.exists(dirname):
        #os.makedirs(dirname)

    lines = ''
    if is_object:
        lines = '%s\n' % class_name
        lines += '%s\n' % (len(class_name) * '-')
        lines += '.. autoclass:: %s.%s\n' % (dot_path, class_name)
        lines += '    :inherited-members:\n'
        lines += '    :members:\n'
        #lines += '    :private-members:\n'
    else:
        lines = '%s\n' % class_name
        lines += '%s\n' % (len(class_name) * '-')
        lines += '.. autoclass:: %s.%s\n' % (dot_path, class_name)
        lines += '    :show-inheritance:\n'
        lines += '    :inherited-members:\n'
        lines += '    :members:\n'
        #lines += '    :private-members:\n'

    #ASET
    #----
    #.. autoclass:: pyNastran.bdf.cards.bdf_sets.ASET
        #:show-inheritance:
        #:inherited-members:
        #:members:
        #:private-members:

    #print(rst_filename)
    if lines:
        with open(rst_filename, 'w') as rst_file:
            rst_file.write(lines)

def main():
    if os.path.exists('rsts'):
        shutil.rmtree('rsts')
    filenames_classes = get_pyfilenames()
    filenames_to_rsts(filenames_classes, make_rsts=False)
    #py_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\bdf\cards\bdf_sets.py'
    #get_classes_in_file(py_filename)

if __name__ == '__main__':
    main()
