from __future__ import print_function
import os
import os.path
import re

def get_folders_files(dirname, skip_file_sufix = ('.pyc', '.pyx', '.bdf'),
                            skip_directories = ('.svn', '.idea', '.settings')):
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


def print_tree(files, padding = 4):
    """
    Prints given list of files (usually returned by get_folders_files) in
    a graphical form of a tree.
    @param files list of files to be printed
    @padding number of spaces to use when padding tree nodes, default is 4
    """
    def divide_path(path): # divide full path into parts
        part = path 
        res = []
        while part: 
            mark = os.sep if os.path.isdir(path) else ""
            path, part = os.path.split(path)
            # if part is empty take the root (or drive letter)
            res.append(part + mark if part else path )
        return list(reversed(res))
    
    paths =  sorted([divide_path(i) for i in files])
    
    npaths = [(0, paths[0])]
    for i in range(1,len(paths)): #max common prefix
        l, max_len = 0, min(len(paths[i-1]), len(paths[i]))
        while l < max_len and paths[i-1][l] == paths[i][l]:
            l += 1
        npaths.append( (l, paths[i][l:]) ) # padding size and new names

    #print the tree
    for padd, path in npaths:
        for i, part in enumerate(path):
            print(' ' * padding * (padd + i) + part)
    

def update_copyright():
    """
    Update copytight info in all .py pyNastran files.
    """
    import pyNastran
    dirname = os.path.relpath(os.path.dirname(pyNastran.__file__))
    print("dirname = %s" %(dirname))
    files = [i for i in get_folders_files(dirname)[1] if i.endswith('.py')]
    
    with open(os.path.join(os.path.dirname(dirname),'licenses','header.txt'),'r') as header:
        header_txt = header.read()
        
    for fname in files:
        print("adding fname=%s" % (fname))
        with open(fname,'r') as fin:
            lines = fin.read()
        #remove copyright lines 
        for regx in ("^## GNU .*\n(##.*\n)*?##.*<http://www.gnu.org/licenses/>.\n([# ]*\n)*",
                     "^# GNU .*\n(#.*\n)*?#.*file will be lost.\n([# ]*\n)*"):
            lines = re.sub(regx, '', lines)
            
        with open(fname,'wb') as fout:
            fout.write(header_txt + lines)
      
    print("---FILES---")
    print_tree(files)

if __name__ == '__main__':
    pass
    #update_copyright()
    