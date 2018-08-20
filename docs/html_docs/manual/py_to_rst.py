from __future__ import print_function
import os
import shutil
from importlib import import_module
from docutils import core

import pyNastran
from pyNastran.utils import check_path
pkg_path = pyNastran.__path__[0]

class Section(object):
    def __init__(self, func_name):
        if 'def ' in func_name:
            func_name = func_name.split('def ')[1]
        self.func_name = func_name
        self.name = None
        self.level = None
        self.blocks = []
        self.code_block = ''
        self.note = ''
        self.iblock = 0
        self.new_block()

    def save_block(self):
        #print('saving %s' % self.iblock)
        self.iblock += 1
        self.blocks.append([self.note, self.code_block])
        self.code_block = ''
        self.note = ''

    def new_block(self):
        if self.code_block:
            self.blocks.append([self.note, self.code_block])
            #print('saving %s' % self.iblock)
            self.iblock += 1
            self.code_block = ''
            self.note = ''

    def write(self):
        msg = ''

        n = len(self.name)
        if self.level == 1:
            msg += '=' * n + '\n'
            msg += '%s\n' % self.name
            msg += '=' * n + '\n'
        elif self.level == 2:
            msg += '-' * n + '\n'
            msg += '%s\n' % self.name
            msg += '-' * n + '\n'
        elif self.level == 2:
            msg += '%s\n' % self.name
            msg += '=' * n + '\n'
        elif self.level == 3:
            msg += '%s\n' % self.name
            msg += '-' * n + '\n'
        else:
            raise RuntimeError(self.func_name, self.level, msg)

        #print('blocks =', self.blocks)
        for note, block in self.blocks:
            msg += note
            if note:
                msg += '\n'
            if block and block[0] == ' ':
                msg += '.. code-block:: python\n'
            if block:
                msg += '\n%s' % block + '\n'
        return msg.rstrip('> \n') + '\n\n'

    def __repr__(self):
        return self.write()

def run_import(fname):
    #p, m = fname.rsplit('.', 1)
    print('fname = %r' % fname)
    p = fname[:-3].replace(os.sep, '.')
    print('p = %r' % p)
    m = 'main'
    print(p)
    mod = import_module(p)
    main = getattr(mod, m)
    main()

def create_html_from_rst(rst_filename, html_filename):
    #os.system('')
    with open (rst_filename, 'r') as source, open(html_filename, 'w') as destination:
        core.publish_file(source=source, destination=destination, writer_name='html')


def create_rst_from_python_files():
    quick_start_pydocs_dir = os.path.join(pkg_path, '..', 'docs', 'quick_start', 'py_docs')
    pydocs_dir = os.path.join(pkg_path, '..', 'docs', 'html_docs', 'quck_start')
    if not os.path.exists(pydocs_dir):
        os.makedirs(pydocs_dir)

    check_path(quick_start_pydocs_dir, 'quick_start_pydocs_dir')
    for fname in os.listdir(quick_start_pydocs_dir):
        fname1 = os.path.join(quick_start_pydocs_dir, fname)
        fname2 = os.path.join(pydocs_dir, fname)
        print(fname1)
        print(fname2)
        shutil.copyfile(fname1, fname2)
    return
    #fnames = [
        #os.path.join('py_docs', fname) for fname in os.listdir('py_docs')
        #if '__init__' not in fname
        #and fname.endswith('.py')
        #and not fname.endswith('.pyc')
        #and not fname.endswith('.html')
        #and not fname.endswith('.txt')
        #and not fname.endswith('.rst')
    #]

    #print(fnames)
    #for py_filename in fnames:
        #root = os.path.splitext(py_filename)[0]
        #base = os.path.splitext(os.path.basename(py_filename))[0]

        #rst_filename = base + '.rst'
        #html_filename = root + '.html'
        #print(rst_filename)
        #print(py_filename)
        #create_rst_from_python(py_filename, rst_filename)
        #create_html_from_rst(rst_filename, html_filename)

    #print(os.path.abspath(py_filename))
    #for py_filename in fnames:
        #root = os.path.splitext(py_filename)[0]
        #base = os.path.splitext(os.path.basename(py_filename))[0]
        #rst_filename = base + '.rst'
        #print('ran %s; saved %s' % (py_filename, rst_filename))


def create_rst_from_python(py_filename, rst_filename, debug=False):
    sections = []
    run_import(py_filename)
    with open(py_filename, 'r') as py_file:
        is_function_active = False
        is_header_active = False

        i = 0
        while 1:
            line = py_file.readline()
            if 'def main():' in line:
                section.save_block()
                break

            if not line:
                if debug:
                    print('%r' % line)
                i += 1
            if i == 2:
                break

            if debug:
                print(line)
            if is_function_active and is_header_active:
                #print('a')
                if 'Level' in line:
                    name, level = line.strip().split(' - Level ')
                    #print('name =', name)
                    #print(section)
                    section.level = int(level)
                    section.name = name
                elif '"""' in line:
                    is_header_active = False
                    nspaces = line.index('"')
                else:
                    raise RuntimeError('a', line)

            elif is_function_active and not is_header_active:
                #print('b')
                line2 = line[nspaces:]
                if debug:
                    print('b; line2=%r' % line2)

                if line2 == 'pass\n':
                    section.save_block()
                    continue

                if line2.startswith('#'):
                    if line2 == '##\n':
                        if debug:
                            print('b0')
                        section.code_block += '    \n'
                    elif line2.startswith('###'):
                        if debug:
                            print('b1')
                        new_line = line2[4:]
                        if len(new_line) == 0:
                            new_line = '\n'
                        if debug:
                            print('&&&& %r' % new_line)
                        section.code_block += new_line
                    elif line2.startswith('##'):
                        if debug:
                            print('b1')
                            print('&&&& %r' % line2[3:])
                        section.code_block += '    ' + line2[3:]
                        #section.save_block()
                    elif line2.startswith('#'):
                        if debug:
                            print('b2')
                        section.save_block()
                        section.note += line2[2:]
                    else:
                        raise RuntimeError('b', line)

                elif 'def ' in line:
                    section.save_block()

                    func_name = line.split('(')[0]
                    if debug:
                        print('func_name = %r' % func_name)
                    if func_name == 'main':
                        break
                    section = Section(func_name)
                    #print('sec**', section)
                    sections.append(section)
                    is_function_active = False
                    is_header_active = False
                else:
                    if debug:
                        print('b3 = %r' % line)
                    if not line2.endswith('\n'):
                        #line2 += '\n'
                        continue
                    section.code_block += '    >>> ' + line2
                    if debug:
                        print('line2 = %r' % line2)
                    #print(section.code_block)

            elif not is_function_active:
                if debug:
                    print('c')
                if 'import ' in line:
                    pass
                elif line.startswith('def '):
                    func_name = line.split('(')[0]
                    section = Section(func_name)
                    #print('sec**', section)
                    sections.append(section)
                elif line.strip().startswith('"""'):
                    is_function_active = True
                    is_header_active = True
                elif line.strip() == '':
                    pass
                else:
                    print(line)
                    #raise RuntimeError(line)
            else:
                raise RuntimeError('c', line)

    print('-----------------------------------')
    rst = open(rst_filename, 'w')
    fname = os.path.basename(rst_filename).split('.rst')[0][:-4] # remove the _doc
    n = len(fname)
    print(fname)
    msg = '=' * n + '\n'
    msg += '%s\n' % fname
    msg += '=' * n + '\n\n'
    rst.write(msg)

    for section in sections:
        #print(section)
        rst.write(str(section))

if __name__ == '__main__':
    create_rst_from_python_files()

