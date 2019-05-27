def read_dmap(lines):
    lines2 = []
    append_flag = False
    for line in lines:
        line = line.strip()
        if line == '$':
            continue
        elif len(line) == 0:
            continue


        if line.endswith(','):
            lines2.append(line)
            append_flag = True
        elif append_flag:
            lines2[-1] += line
            append_flag = False
        else:
            lines2.append(line)
            append_flag = False

    lines = [line.upper() for line in lines2]
    del lines2

    #for line in lines:
        #print(line)

    sol = None
    files = {}
    lines2 = []
    code = ""
    echo = True
    level = 0
    spaces = ''
    for line in lines:
        line = ' '.join(line.split())

        comment = ''
        if '$' in line:
            line, comment = line.split('$', 1)
            comment = comment.strip()
        line = line.strip()
        if len(line) == 0:
            continue

        if echo:
            print(comment.rstrip())
            print('line=%r' % line.rstrip())

        if line.startswith('ASSIGN'):
            line = line.replace(' =', '=').replace('= ', '=')
            sline = line.split()
            if len(sline) == 3:
                unused_assign, output, unit = sline
                word, unit_num = unit.split('=')
                output, fname = output.split('=')
                assert word == 'UNIT', word
                assert output == 'OUTPUT2', output
                unit_num = int(unit_num)
                files[unit_num] = [output, fname, unit]
                code = 'from pyNastran.op2.op2 import OP2\n'
                code += 'op2 = OP2()\n'
                #code += ''
            else:
                raise SyntaxError(sline)
        elif line.startswith('IF'):
            assert line.endswith(') THEN'), line
            line = line.replace(') THEN', '):')
            if '=' in line:
                code += line.replace('=', '==').lower() + '\n'
            else:
                raise RuntimeError(line)
            level += 1

        elif line.startswith('TYPE PARM,'):
            line = line[5:]
            #print('line = %r' % line)
            sline = line.split(',')

            # broken for:
            #   TYPE PARAM,,CS,Y,ALPHA=(1.,1.)
            #
            # works for
            #    TYPE PARAM,,CS,Y,ALPAH=1.
            unused_base, unused_NDDL, unused_Type, unused_y, word_val = sline
            param_name, default_value = word_val.split('=')
            code += 'model.params[%r].set_value(%s)\n' % (param_name, default_value)

        elif line.startswith('LAMX '):
            sline = line[5:].split('/')
            print('sline =', sline)
            if len(sline) == 4:
                pre, unused_b, unused_c, resflag = sline
                unused_resflg = int(resflag)
                if pre.startswith('FREQMASS,'):
                    raise NotImplementedError('FREQMASS')
                elif pre.startswith(',,LAMA'):
                    code += spaces + 'LAMX = LAMA\n'
                    if resflag == 1:
                        code += spaces + 'LAMX.msg = "Before Augmentation of Residual Vectors"\n'
                    elif resflag == 2:
                        code += spaces + 'LAMX.msg = "After Augmentation of Residual Vectors"\n'
                    else:
                        code += spaces + 'LAMX.msg = ""\n'
                else:
                    raise RuntimeError(line)
        elif line.startswith('ECHOOFF'):
            #echo = False
            pass
        elif line.startswith('MATMOD'):  # super hard
            pass
        elif line.startswith('MATGEN'):
            line = line[7:].lstrip(',/')
            sline = line.split('/')
            print(sline)
            name = sline[0]
            option = int(sline[1])
            if option == 6:
                sline = sline[2:]
                code += spaces + '%s = \n'
            else:
                raise NotImplementedError(option)

        elif line.startswith('SOL '):
            word, sol = line.split()
            sol = int(sol)
        elif line.startswith('DIAG '):
            # request DIAGnostic output
            pass
        elif line.startswith('COMPILE USERDMAP'):
            pass
        elif line.startswith('COMPILE SEDRCVR'):
            pass
        elif line.startswith('ALTER'):
            pass
        elif line.startswith('DMIIN DMI,DMINDX/'):
            sline = line.split('/')
            pre, data, unused_post = sline
            matrices = data.replace(' ', '').rstrip(',')
            matrices_to_load = matrices.split(',')
            print('matrices_to_load =', matrices_to_load)

            code += 'from pyNastran.bdf.bdf import BDF\n'
            code += 'model = BDF()\n'
            code += 'model.read_bdf(%s, xref=True)\n\n' % 'bdf_filename'
            for matrix_name in matrices_to_load:
                code += "%s = model.dmis[%r].get_matrix()\n" % (matrix_name, matrix_name)
            del sline
            code += '\n'
        elif line.startswith('MPYAD '):
            sline = line[6:].replace(' ', '').split(',')
            a_name, b_name, other = sline
            sline2 = other.strip('/').split('/')
            if len(sline2) == 2:
                out_name, transpose_method = sline2
                transpose_method = int(transpose_method)
                if transpose_method == 1:
                    code += '%s = %s.T @ %s\n' % (out_name, a_name, b_name)
                elif transpose_method == 0:
                    code += '%s = %s @ %s\n' % (out_name, a_name, b_name)
                else:
                    raise RuntimeError(transpose_method)
            else:
                raise RuntimeError(sline2)

            print('line=%r' % line)
        elif line.startswith('OUTPUT2 '):
            sline = line[8:].replace(' ', '').split('/')
            if len(sline) == 4:
                output_matrices, unused_dunno, unused_method, unit_num = sline
                matrix_names = output_matrices.split(',')
                code += '\n'
                for name in matrix_names:
                    code += 'op2.write(%s)\n' % name
            else:
                raise RuntimeError(sline)

            #print('line=%r' % line2)
            #adsf
        elif line.startswith('CEND'):
            pass
        else:
            msg = 'line=%r\n\n' % line
            msg += code
            raise SyntaxError(msg)

    print('-----------')
    print(code)

def main():
    lines = [
        "ASSIGN OUTPUT2='myMatrix.op2',",
        '    UNIT=15',
        '$',
        'SOL 100',
        'DIAG 8,44',
        'COMPILE USERDMAP',
        'ALTER 2',
        'DMIIN DMI,DMINDX/A,B,MYDOF,,,,,,,/ $',
        'MPYAD A,B,/ATB/1///$',
        'MPYAD B,A,/BTA/1///$',
        'OUTPUT2 A,B,ATB,BTA,MYDOF//0/15$',
        'CEND',
    ]
    #with open('isat.bdf', 'r') as f:
        #lines = f.readlines()
    read_dmap(lines)

if __name__ == '__main__':  # pragma: no cover
    main()
