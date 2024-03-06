def read_ansys(ansys_filename, log=None, debug=False):
    model = Ansys()
    model.read_ansys(ansys_filename)
    return model

class Ansys:
    def __init__(self, log=None, debug=False):
        pass

    def read_ansys(self, ansys_filename):
        with open(ansys_filename, 'r') as ansys_file:
            lines = ansys_file.readlines()

        #nodes = []
        elements = {}

        i = 0
        nlines = len(lines)
        while i < nlines:
            line = lines[i].strip()
            if line.startswith(r'/nolist'):
                print('line = %s' % line)
                i += 4
                # line = (1i9,3e20.9e3)
                snodes = []
                i += 1
                line = lines[i]
                nnodes = 0
                while not line.startswith('-1'):
                    #print('a =', line)
                    #snode = [float(val) for val in line.strip().split()[1:]]
                    snode = line.strip().split()[1:]
                    if len(snode) != 3:
                        print(snode)
                        print(line)
                        print(lines[i])
                        print(lines[i-1])
                        print(lines[i-2])
                        asdf1
                    snodes.append(snode)
                    line = lines[i]
                    #print(line)
                    i += 1
                    nnodes += 1
                #print(snodes[:5])
                #nodes = array(snodes, dtype='float32')
                print('****%r' % line)
                # nnodes = 793310
                #asdf2
                #line = lines[i]
                #print(line)
                i -= 1
                #asdf
            elif line.startswith('/wb,elem,start'):
                print('line = %s' % line)
                i += 1
                line = lines[i]
                while line.startswith('/com'):
                    i += 1

                    et_line = lines[i].strip()
                    fmt_line = lines[i+2].strip()
                    i += 3
                    line = lines[i]
                    if fmt_line == '(19i9)':
                        # eblock,19,solid,,71892
                        while not line.startswith('-1'):
                            #        1        1        1        1        0        0        0        0       10        0   697401  1297419  1304724  1297455  1302783  2097856  2097997  2097853  2097855
                            #  2109421  2097995
                            #       27       27       27       27        0        0        0        0       10        0   387759   631841   659167   639072   631842   675592   723723   675588   675585
                            #   675599   675595
                            line = lines[i].strip() + lines[i+1].strip()
                            i += 2
                            print(line)

                            sline = line.split()
                            a = sline[0]
                            b = sline[1]
                            c = sline[2]
                            d = sline[3]
                            assert a == b, 'a=%r b=%r c=%r d=%r' % (a, b, c, d)
                            assert a == c, 'a=%r b=%r c=%r d=%r' % (a, b, c, d)
                            assert a == d, 'a=%r b=%r c=%r d=%r' % (a, b, c, d)

                            e = sline[3]
                            f = sline[4]
                            g = sline[5]
                            h = sline[6]
                            assert e == f, 'e=%r f=%r g=%r h=%r' % (e, f, g, h)

                            #asdf
                    else:  # pragma: no cover
                        raise NotImplementedError(fmt_line)
                    print(line)
                    asdf
            else:
                if line.startswith('/'):
                    print(line)
            i += 1


def main():
    model = Ansys()
    ansys_filename = 'ds.dat'
    model.read_ansys(ansys_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
