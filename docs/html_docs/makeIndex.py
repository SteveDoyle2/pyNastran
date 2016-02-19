from __future__ import print_function
from six import iteritems
import pyNastran
import os

pkg_path = pyNastran.__path__[0]

ignore_dirs = ['src', 'dmap', 'solver', '__pycache__', 'converters',
               'op4_old', 'cart3d', 'panair', 'calculix', 'bars', 'case_control',
               'gui', 'gui_qt', 'pch', 'mesh', 'LaWGS', 'old', 'solver', 'test']
mods_skip = ['spike', 'shell_backup']

def get_folders_files(dirname, skip_file_sufix=('.pyc', '.pyx', '.bdf'),
                      skip_directories=('.svn', '.idea', '.settings')):
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

def get_classes(filenames):
    d = {}
    for filename in filenames:
        #print "filename =", filename
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if 'class' in line:
                print(line.rstrip())
                sline = line.split('(')
                if len(sline) == 2:
                    if 'raise' in sline[0]:
                        continue
                    pre, post = sline
                else:
                    print("*******")
                    continue

                pre, post = pre.split('class')
                name = post.strip()
                print("name =", name)
                dirname = os.path.dirname(filename)
                d[dirname] = name

def get_doc(filenames):
    dirs = {}
    for filename in filenames:
        if '.py' not in filename[-3:]:
            continue

        dirname = os.path.dirname(filename)
        dirname = os.path.join(dirname)
        if dirname == r'../pyNastran':
            continue

        if dirname not in dirs:
            dirs[dirname] = []
        dirs[dirname].append(filename)

    for dirname, files in sorted(iteritems(dirs)):
        print("***************************************")
        basepath = os.path.basename(dirname)
        if basepath in ignore_dirs:
            continue
        if 'dev_' in basepath:
            continue

        print("dirname = ", dirname)
        #absname = os.path.join(dirname, filename)
        dirnamebase = dirname
        dirname2 = dirnamebase.strip('./')
        dirname2 = dirname2.replace('\\', '.')

        print('dirname2 = ', dirname2)
        fpath = os.path.join('index', '%s.rst' % dirname2) # basepath, dirname2

        msg = '%s Package\n' % basepath  # loads Package
        msg += '=' * (len(msg)-1) + '\n'
        msg += '\n'

        print("----------------------")
        #print 'absname2 = ', absname2)
        print('fpath = ', fpath)
        print("basepath = ", basepath)

        if len(files) == 0:
            continue

        f = open(fpath, 'wb')
        f.write(msg)
        for fname in files:
            if os.path.getsize(fname) == 0:
                continue
            assert os.path.exists(fname), fname

            fnamebase = os.path.basename(fname)
            fnamebase2 = fnamebase[:-3]

            fname2 = fname[:-3]
            fname2 = fname2.strip('/.')

            modname = fnamebase2
            automodname = fname2.replace('\\', '.')

            #print "fname = ", fname
            print("fname2 = ", fname2)

            if modname in mods_skip:
                continue

            msg = ':mod:`%s` Module\n' % modname  # `loads\  # fname2
            msg += '-' * (len(msg)-1) + '\n'
            msg += '\n'
            msg += '.. inheritance-diagram:: %s\n' % automodname
            msg += '\n'
            msg += '.. automodule:: %s\n' % automodname # pyNastran.bdf.bdfInterface.bdf_card, pyNastran.bdf.cards.loads.static_loads
            msg += '    :members:\n'
            msg += '    :private-members:\n'
            msg += '    :undoc-members:\n'
            msg += '    :show-inheritance:\n'

            #msg += '    :autofunction:\n'
            #msg += '    :autodata:\n'
            #msg += '    :automethod:\n'
            #msg += '    :autoattribute:\n'

            msg += '\n'
            f.write(msg)

        dirs = get_sub_dirs(dirname)
        if dirs:
            msg += 'Subpackages\n'
            msg += '-----------\n'
            msg += '\n'
            msg += '.. toctree::\n'
            msg += '\n'
            for diri in dirs:
                diri2 = diri.strip('/.')
                diri2 = diri2.replace('\\', '.')
                print("diri2 =", diri2)
                msg += '    %s\n' % diri2
            msg += '\n'
            f.write(msg)
        f.close()

def get_sub_dirs(dirname):
    maybe_dirs = os.listdir(dirname)
    dirs = []
    print("dirname_check =", dirname)

    for idir in maybe_dirs:
        dir2 = os.path.join(dirname, idir)
        if os.path.isdir(dir2):
            if idir not in ignore_dirs:
                print("adding %s" % dir2)
                dirs.append(dir2)
    return dirs

def run():
    folders, files = get_folders_files('../pyNastran')
    #c = get_classes(files)
    get_doc(files)


if __name__ == '__main__':  # pragma: no cover
    run()
