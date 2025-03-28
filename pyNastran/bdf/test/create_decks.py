"""simple script for mass running of Nastran jobs to create OP2 files"""
import  os
import sys
from pyNastran.op2.test.op2_test import get_files_of_type # get_all_files,
from pyNastran.bdf.bdf_interface.utils import to_fields
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.nastran_utils import run_nastran


# no files are written if RUN_NASTRAN is False
RUN_NASTRAN = True

# nx 2019
#NASTRAN_EXE = r'C:\Program Files\Siemens\SimcenterNastran_2019.2\bin\nastran.exe'

# msc 2020
#NASTRAN_EXE = r'C:\MSC.Software\MSC_Nastran_2020\bin\nastran.exe' #bad
NASTRAN_EXE = r'C:\MSC.Software\MSC_Nastran_2020\bin\nast20200.exe'



solutions = (
    'SOL 1',
    'SOL 3', # modal
    'SOL 5',
    'SOL 7',
    'SOL 8',
    'SOL 21',
    'SOL 24',
    'SOL 25',
    'SOL 26',
    'SOL 27',
    'SOL 28',
    'SOL 29',
    'SOL 47',
    'SOL 48',
    'SOL 60',
    'SOL 61',
    'SOL 62',
    'SOL 63',
    'SOL 64',
    'SOL 66',
    'SOL 67',
    'SOL 68',
    'SOL 74',
    'SOL 75',
    'SOL 76',
    'SOL 81',
    'SOL 82',
    'SOL 83',
    'SOL 88',
    'SOL 89',
    'SOL 91',
    'SOL 99',
    #'SOL 3',
    'SOL 100',
    'SOL USERDMAP',
    'SOL MAIN',
    'SOL XXX',
    'SOL 101', 'SOL STATIC', 'SOL SESTATIC',
    'SOL 103', 'SOL SEMODES',
    'SOL 105', 'SOL BUCKLING',
    'SOL 106', 'SOL NLSTATIC', # nonlinear static
    'SOL 107', 'SOL SEMFREQ', # direct complex frequency response
    'SOL 108', 'SOL SEDFREQ', # direct frequency response
    'SOL 109', 'SOL SEDTRAN', # time linear?
    'SOL 110', 'SOL SEDCEIG', # modal complex eigenvalue
    'SOL 111', 'SOL SEMFREQ', # modal frequency response
    'SOL 112', 'SOL SEMTRAN', # modal transient response
    'SOL 114',  #
    'SOL 115',  #
    'SOL 118',  #
    'SOL 126',
    'SOL 129', 'SOL NLTRAN', # time nonlinear
    'SOL 144', # static aero
    'SOL 145', # flutter
    'SOL 146', # gust
    'SOL 153', 'SOL NLSCSH', # nonlinear thermal
    'SOL 159', # nonlinear transient thermal
    'SOL 190', 'SOL DBTRANS',
    'SOL 200',  # optimization
)
def main():
    print(sys.argv)
    #export_dir = sys.argv[1]
    export_dir = r'C:\MSC.Software\msc_nastran_runs'
    print('export_dir = %r' % export_dir)

    #folders_file = os.path.abspath('.')
    #files2 = get_all_files(folders_file, '.bdf')
    #files2 += get_all_files(folders_file, '.nas')
    #files2 += get_all_files(folders_file, '.dat')
    #files2 = list(set(files2))
    #files2.sort()

    max_size = 4.2 # MB
    dirname = r'C:\MSC.Software\msc_nastran_docs_2020\tpl'
    files2 = get_files_of_type(dirname, '.bdf', max_size_mb=max_size)
    files2 += get_files_of_type(dirname, '.nas', max_size_mb=max_size)
    files2 += get_files_of_type(dirname, '.dat', max_size_mb=max_size)
    files2 += get_files_of_type(dirname, '.ecd', max_size_mb=max_size)
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

    keywords = 'scr=yes old=no news=no notify=no'
    for fname in files:
        #print(fname)
        try:
            fname2 = update_with_post(fname, export_dir)
        except Exception as e:
            print(str(e))
            continue
        if fname2 is None:
            continue
        if not RUN_NASTRAN:
            continue
        run_nastran(fname2, nastran_cmd=NASTRAN_EXE, keywords=keywords, run=True, run_in_bdf_dir=True)

def load_lines(fname):
    encodings = ['latin1', 'cp1252']
    for encoding in encodings:
        try:
            with open(fname, 'r', encoding=encoding) as bdf_file:
                lines = bdf_file.readlines()
            break
        except UnicodeDecodeError:
            continue
    return lines

def update_with_post(fname: str, dirname: str) -> str:
    #print('-'*40)
    print('fname =', fname)
    basename = os.path.basename(fname)
    base = os.path.splitext(basename)[0]
    bdf_name2 = os.path.join(dirname, basename)
    f04_name2 = os.path.join(dirname, base + '.f04')
    f06_name2 = os.path.join(dirname, base + '.f06')
    op2_name2 = os.path.join(dirname, base + '.op2')
    log_name2 = os.path.join(dirname, base + '.log')

    lines = load_lines(fname)
    if len(lines) == 0:
        return None

    post = -1
    found_post = False
    lines2 = None
    for i, line in enumerate(lines):
        #if i < 50:
            #print(line.rstrip())
        line_upper = line.upper().split('$')[0].rstrip()

        line_upper = line_upper.replace('\t', ' ').strip()
        if line_upper.startswith('SOL '):
            while '  ' in line_upper:
                line_upper = line_upper.replace('  ', ' ')

            if any(sol in line_upper for sol in solutions):
                continue

            print('%r' % line_upper)
            if '601' in line_upper:
                # SOL 601 doesn't have a valid license
                #SOL 601,129
                if os.path.exists(bdf_name2):
                    print(f'removing {bdf_name2}')
                    os.remove(bdf_name2)
                if os.path.exists(f04_name2):
                    os.remove(f04_name2)
                if os.path.exists(f06_name2):
                    os.remove(f06_name2)
                if os.path.exists(op2_name2):
                    os.remove(op2_name2)
                if os.path.exists(log_name2):
                    os.remove(log_name2)
                return None
        if line_upper.startswith(('LABEL', 'TITLE', 'SUBTITLE', 'SUBTI')):
            pass
        elif 'PRGPOST' in line_upper or 'POSTEXT' in line_upper or 'AKMPOST' in line_upper or 'POSTADF' in line_upper:
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
                lines2 = lines[:i] + [f'PARAM,POST,{post}\n', 'ENDDATA\n'] + lines[i:]
            break
        elif 'ENDATA' in line_upper:
            print(fname)
            print('*en_data')
            raise RuntimeError('ENDATA...')
    #print('------------------')
    #print('i = ', i)
    if lines2 is None:
        lines2 = lines[:i] + [f'PARAM,POST,{post}\n'] + lines[i:]

    assert lines2 is not None, fname
    #print("bdf_name2 =", bdf_name2)
    #for line in lines2:
        #print(line)
    with open(bdf_name2, 'w') as bdf_file:
        try:
            bdf_file.writelines(lines2)
        except UnicodeEncodeError:
            # just run the original file
            return fname
    #sss
    if not RUN_NASTRAN:
        return None
    return bdf_name2

if __name__ == '__main__':   # pragma: no cover
    main()
