import os
import pyNastran
pkg_name = pyNastran.__path__[0]
#from shutil import rmtree

def get_base_import_names():
    publics = []
    for module_name in ['__future__', 'os', 'os.path', 'six', 'six.moves',
                        'numpy', 'numpy.linalg', 'scipy', 'vtk', 'docopt', 'wx',
                        'types', 'math', 'h5py', 'sympy', 'tables', 'struct']:
        module = importlib.import_module(module_name)
        public = [key for key in sorted(module.__dict__)
                  if not key.startswith('__') and not key.endswith('__')]
        publics += public

    publics += ['scipy', 'main', 'wx', 'range', 'iteritems', 'zip', 'unittest',
        'traceback', 'logging', 'io', 'run', 'string', 'mp', 'quad', 'splev',
        'splrep', 'pkg_path', 'is_wx', 'pyNastran', 'testpath', 'rootpath',
        'log_formater', 'log_handler_file', 'logger', 'get_logger', 'get_logger2',
        'bdf', 'comment_bad', 'comment_good', 'model', 'root_path',
        'unpack', 'pack',
        'make_log', 'platform', 'stderr_logging', 'utils', 'date',
        #'coo_matrix', 'dim', '_TemporarilyImmutableSet', 'BaseSet', 'ImmutableSet',
        #'ifilterfalse',  'Array', 'Atom', 'BoolAtom', 'BoolCol', 'CArray',
        #'ClosedFileError', 'ClosedNodeError', 'Complex128Col', 'Complex32Atom',
        #'Complex32Col', 'Complex64Atom', 'Complex64Col', 'ComplexAtom',

        ]
    return publics

def get_folders_files(dirname, skip_file_sufix = ('.pyc', '.pyx', '.bdf', '.cpp'),
                            skip_directories = ('.svn', '.idea', '.settings',
                            'applications', # bad
                            'dev_vectorized', # done
                            'converters', # bad
                            'nastranwrapper-using_pynastran', # bad
                            #'bdf', # done
                            #'op2', # done
                            'op4', # cop4
                            'pch', # bad
                            )):
    """
    Return list of directories and files in a given tree path. By default discards:

    * directories ".svn", ".idea", ".settings"

    * files that ends with ".pyc", .pyx", ".bdf"
    """
    files = []
    folders = []
    for root, dirs, fil in os.walk(dirname):
        folders.append(root)
        files += [os.path.join(root, i) for i in fil
                  if not i.endswith(skip_file_sufix)]
        dirs[:] = [d for d in dirs if not d in skip_directories]

    return (folders, files)

def cleanup():
    rmtree(pkg_name, 'nastranwrapper-using_pynastran')
    rmtree(pkg_name, 'pch')
    rmtree(pkg_name, 'latex')
    rmtree(pkg_name, 'converters', 'calculix')
    rmtree(pkg_name, 'converters', 'cgns')
    rmtree(pkg_name, 'converters', 'iges')
    rmtree(pkg_name, 'converters', 'plot3d')
    rmtree(pkg_name, 'converters', 'tecplot')
    os.remove(pkg_name, 'converters', 'toCodeAster.py')
    os.remove(pkg_name, 'converters', 'toSfepy.py')
    os.remove(pkg_name, 'converters', 'nastran', 'bdf_to_p3d.py')

    rmtree(pkg_name, 'applications')
    #rmtree(pkg_name, 'applications', 'solver')
    #rmtree(pkg_name, 'applications', 'hyper')
    #rmtree(pkg_name, 'applications', 'contact')
    #rmtree(pkg_name, 'applications', 'aero_tet_mesh')


def get_class_function_list(filename):
    f = open(filename, 'r')
    lines = f.read().split('\n')
    class_names = []
    function_names = []
    for line in lines:
        if line.startswith('class '):
            class_name = line.split('(')[0].split(' ')[1]
            class_names.append(class_name)
            #print('  class ' + class_name)
        if line.startswith('def '):
            function_name = line.split('(')[0].split(' ')[1]
            function_names.append(function_name)
            #print('  func ' + function_name)
    return class_names, function_names

import sys
import importlib
from pyNastran.utils import object_methods
import inspect

def get_function_args(function, function_name):
    f = inspect.getargspec(function_name)
    print(f)
    msg = '%s('
    for argi in f.args:
        msg += argi + ', '
    msg += ')'
    print(msg)
    asf

def main():
    base_imports = get_base_import_names()

    #rmtree(pkg_name, 'nastranwrapper-using_pynastran')
    #cleanup()
    folders, fnames = get_folders_files(pkg_name)
    n = len(sys.path)
    for fname in fnames:
        class_names, function_names = get_class_function_list(fname)
        if len(class_names + function_names):

            dirname = os.path.dirname(fname)
            basename = os.path.basename(fname).split('.')[0]
            basename2 = fname.split('.')[0]
            #print('basename =', basename)
            if 'gui' in fname or 'test' in fname:
                continue
            sys.path.insert(n, dirname)
            module = importlib.import_module(basename)

            # remove __methods__
            public = [key for key in sorted(module.__dict__)
                      if not key.startswith('__') and not key.endswith('__')
                      and key not in base_imports and 'test' not in key.lower()]
            if public:
                print(fname)
                print(public)

                for class_name in class_names:
                    if class_name not in public:
                        continue
                    try:
                        class_obj = getattr(module, class_name)
                    except AttributeError:
                        print('bad_module=%s' % fname)
                        raise
                    object_methods(class_obj)

                if function_names:
                    for function_name in function_names:
                        func = getattr(module, function_name)
                        fargs = get_function_args(func, function_name)

                    #print(function_names)
                print('')

                #for function_name in function_names:
                    #print(object_methods(class_obj))
            #sys.path.pop()
            #del module
            del sys.modules[basename]

if __name__ == '__main__':
    main()