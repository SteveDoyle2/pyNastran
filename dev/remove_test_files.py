import os
import shutil
import pyNastran

def main():
    PKG_PATH = pyNastran.__path__[0]
    BASE = os.path.abspath(os.path.join(PKG_PATH, '..'))
    test_filename = os.path.join(BASE, 'dev', 'test_files.txt')
    assert os.path.exists(test_filename), test_filename

    with open(test_filename, 'r') as test_file:
        lines = test_file.readlines()

    files = []
    folders = []
    files_folders = [line.rstrip() for line in lines if line.rstrip()]
    for folder_fname0 in files_folders:
        if os.path.isabs(folder_fname0):
            print(f'*skipped absolute path: {folder_fname0}')
            continue
        folder_fname = os.path.abspath(os.path.join(BASE, folder_fname0))
        if os.path.commonpath([BASE, folder_fname]) != BASE:
            print(f'*skipped path outside base: {folder_fname0}')
            continue
        if not os.path.exists(folder_fname):
            print(f'*skipped missing path: {folder_fname0}')
            continue
        assert os.path.exists(folder_fname), folder_fname0
        if os.path.isdir(folder_fname):
            folders.append(folder_fname0)
        else:
            files.append(folder_fname0)

    delete = '--delete' in os.sys.argv
    if not delete:
        print('dry run mode; pass --delete to actually remove files/folders')
    for folder0 in folders:
        folder = os.path.abspath(os.path.join(BASE, folder0))
        if os.path.commonpath([BASE, folder]) != BASE:
            print(f'*skipped folder due to incorrect path: {folder0}')
        elif os.path.exists(folder) and delete:
            print(f'deleting folder: {folder0}')
            shutil.rmtree(folder)
            if os.path.exists(folder):
                print(f'*skipped folder due to incorrect path: {folder0}')
        elif os.path.exists(folder):
            print(f'dry run folder: {folder0}')
        else:
            print(f'*skipped folder due to incorrect path: {folder0}')
    for fname0 in files:
        fname = os.path.abspath(os.path.join(BASE, fname0))
        if os.path.commonpath([BASE, fname]) != BASE:
            print(f'*skipped file due to incorrect path: {fname0}')
        elif os.path.exists(fname) and delete:
            print(f'deleting file: {fname0}')
            os.remove(fname)
            if os.path.exists(fname):
                print(f'*skipped file due to incorrect path: {fname0}')
        elif os.path.exists(fname):
            print(f'dry run file: {fname0}')
        else:
            print(f'*skipped file due to incorrect path: {fname0}')
    x = 1

if __name__ == '__main__':
    main()