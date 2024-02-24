import os
import shutil
import pyNastran

def main():
    PKG_PATH = pyNastran.__path__[0]
    BASE = os.path.join(PKG_PATH, '..')
    test_filename = os.path.join(BASE, 'dev', 'test_files.txt')
    assert os.path.exists(test_filename), test_filename

    with open(test_filename, 'r') as test_file:
        lines = test_file.readlines()

    files = []
    folders = []
    files_folders = [line.rstrip() for line in lines if line.rstrip()]
    for folder_fname0 in files_folders:
        folder_fname = os.path.abspath(os.path.join(BASE, folder_fname0))
        if not os.path.exists(folder_fname):
            continue
        assert os.path.exists(folder_fname), folder_fname0
        if os.path.isdir(folder_fname):
            folders.append(folder_fname0)
        else:
            files.append(folder_fname0)

    delete = False
    for folder0 in folders:
        folder = os.path.abspath(os.path.join(BASE, folder0))
        if os.path.exists(folder) and not delete:
            print(f'deleting folder: {folder0}')
            shutil.rmtree(folder)
        else:
            print(f'*failed to delete folder: {folder0}')
    for fname0 in files:
        fname = os.path.abspath(os.path.join(BASE, fname0))
        if os.path.exists(folder) and not delete:
            print(f'deleting file: {fname0}')
            os.remove(fname)
        else:
            print(f'*failed to delete file: {fname0}')
    x = 1

if __name__ == '__main__':
    main()
