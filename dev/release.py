import os
import pyNastran
pkg_name = pyNastran.__path__[0]
#from shutil import rmtree


def get_folders_files(dirname, skip_file_sufix = ('.pyc', '.pyx', '.bdf', '.cpp'),
                            skip_directories = ('.svn', '.idea', '.settings',
                            'applications', # bad
                            'dev_vectorized', # done
                            'converters', # bad
                            'nastranwrapper-using_pynastran', # bad
                            #'bdf', # done
                            'op2', # done
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

def main():
    #rmtree(pkg_name, 'nastranwrapper-using_pynastran')
    #cleanup()
    folders, fnames = get_folders_files(pkg_name)
    for fname in fnames:
        class_names, function_names = get_class_function_list(fname)
        if len(class_names + function_names):
            #print(fname)

            dirname = os.path.dirname(fname)
            basename = os.path.basename(fname).split('.')[0]
            basename2 = fname.split('.')[0]
            #print('basename =', basename)
            sys.path.append(dirname)
            print(importlib.import_module(basename))

            for class_name in class_names:
                pass
            for function_name in function_names :
                pass

if __name__ == '__main__':
    main()