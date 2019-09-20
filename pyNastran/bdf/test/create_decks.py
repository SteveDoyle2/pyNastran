import  os
import sys
from pyNastran.op2.test.op2_test import get_files_of_type # get_all_files,
from pyNastran.bdf.bdf_interface.utils import to_fields
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.nastran_utils import run_nastran


def main():
    nastran_path = r'C:\Program Files\Siemens\SimcenterNastran_2019.2\bin\nastran.exe'
    print(sys.argv)
    export_dir = sys.argv[1]
    print('export_dir = %r' % export_dir)

    #folders_file = os.path.abspath('.')
    #files2 = get_all_files(folders_file, '.bdf')
    #files2 += get_all_files(folders_file, '.nas')
    #files2 += get_all_files(folders_file, '.dat')
    #files2 = list(set(files2))
    #files2.sort()

    max_size = 4.2 # MB
    dirname = '.'
    files2 = get_files_of_type(dirname, '.bdf', max_size=max_size)
    files2 += get_files_of_type(dirname, '.nas', max_size=max_size)
    files2 += get_files_of_type(dirname, '.dat', max_size=max_size)
    files2 += get_files_of_type(dirname, '.ecd', max_size=max_size)
    files2 = list(set(files2))
    files2.sort()

    skip_files = []
    files = [fname for fname in files2
             if not os.path.basename(fname).startswith('out_')
             and '.test_op2.' not in fname # removing test output files
             and '.test_bdf.' not in fname
             and '.test_bdfv.' not in fname
             and 'tecplot' not in fname
             and os.path.basename(fname) not in skip_files]
    for fname in files:
        #print(fname)
        fname2 = update_with_post(fname, export_dir)
        keywords = 'scr=yes old=no news=no'
        run_nastran(fname2, nastran_cmd=nastran_path, keywords=keywords, run=True, run_in_bdf_dir=True)

def update_with_post(fname, dirname):
    with open(fname, 'r') as bdf_file:
        lines = bdf_file.readlines()

    post = -2
    found_post = False
    lines2 = None
    for i, line in enumerate(lines):
        line_upper = line.upper().split('$')[0].rstrip()
        if line_upper.startswith(('LABEL', 'TITLE', 'SUBTITLE')):
            pass
        elif 'PRGPOST' in line_upper or 'POSTEXT' in line_upper:
            pass
        elif 'PARAM' in line_upper and 'POST' in line_upper:
            nlines_extra = 0
            if '*' in line_upper:
                #print(fname)
                #print(line_upper)
                nlines_extra = 1

            card_lines = [line_upper]
            param_post_fields = [field.strip()
                                 for field in to_fields(card_lines, 'PARAM')[:3]]
            if ' ' in param_post_fields[1].strip():  # PARAM, POST 0
                field2 = param_post_fields[1].split()[1]
                try:
                    unused_found_post_value = int(field2)
                except ValueError:
                    print(fname)
                    print(line_upper)
                    print(param_post_fields)
                    raise

            else:
                try:
                    unused_found_post_value = int(param_post_fields[2])
                except ValueError:
                    print(fname)
                    print(line_upper)
                    print(param_post_fields)
                    raise
            #if post_value == 0:

            param_post_fields = ['PARAM', 'POST', post]
            #print(fname)
            #print(card_lines)
            line = print_card_8(param_post_fields)
            #print(f'post {nlines_extra}!')
            lines2 = lines[:i] + [line] +  lines[i+nlines_extra+1:]
            found_post = True

            #if nlines_extra:
                #print(fname)
                #print(lines[i-2:i+3])
                #print(lines2[i-2:i+3])

        elif 'ENDDATA' in line_upper:
            #print(fname)
            if not found_post:
                #print('enddata not found post')
                lines2 = lines[i:] + [f'PARAM,POST,{post}\n', 'ENDDATA\n'] + lines[i:]
            break
        elif 'ENDATA' in line_upper:
            print(fname)
            print('*en_data')
            ddd

    if lines2 is None:
        lines2 = lines[i:] + [f'PARAM,POST,{post}\n'] + lines[i:]

    assert lines2 is not None, fname

    basename = os.path.basename(fname)
    fname2 = os.path.join(dirname, basename)
    with  open(fname2, 'w') as bdf_file:
        bdf_file.writelines(lines2)
    return fname2

if __name__ == '__main__':
    main()
