#from numpy import array #zeros,

class Ansys(object):
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
                    else:
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

"""
/com,*********** Create Remote Point "Internal Remote Point 39" ***********
! -------- Remote Point Used by "Fixed - Line Body To EndCap 14054021-1 d" --------
*set,_npilot,803315
_npilot474=_npilot
et,332,170
type,332
real,332
mat,332
keyo,332,2,1              ! don't fix pilot node
keyo,332,4,0              ! MPC for all DOF's
tshape,pilo
en,501901,803315        ! create pilot node for rigid link
tshape
en,501902,803315,127827
/com,*********** Create Remote Point "Internal Remote Point 40" ***********
! -------- Remote Point Used by "Fixed - Line Body To EndCap 14054021-1 d" --------
*set,tid,334
*set,cid,333
et,cid,175
et,tid,170
keyo,tid,2,1               ! Don't fix the pilot node
keyo,tid,4,111111
keyo,cid,12,5              ! Bonded Contact
keyo,cid,4,0               ! Rigid CERIG style load
keyo,cid,2,2               ! MPC style contact
mat,333
real,333
type,333
en,501903,418114
en,501904,418115
en,501905,418116
en,501906,418117
en,501907,418118
en,501908,418119
en,501909,418120
en,501910,418121
en,501911,418122
en,501912,418123
en,501913,418124
en,501914,427511
en,501915,427512
en,501916,427518
en,501917,427524
en,501918,427528
en,501919,427533
en,501920,427539
en,501921,427544
en,501922,427551
en,501923,427562
en,501924,427569
*set,_npilot,803316
_npilot475=_npilot
type,tid
mat ,cid
real,cid
tshape,pilo
en,501925,_npilot
tshape



et,2,187
et,27,187   # element, group 27, element_type=187 -> tet10
et,30,188

etype   nastran_name
187     tet10
186     hexa20
188     beam

eblock,19,solid,,213
eblock,19,solid,,8
#----------------------------------------------------------------
et,_jid,184
et,tid,170
et,cid,174


keyo,tid,2,1               ! Don't fix the pilot node
keyo,tid,4,111111
keyo,cid,12,5              ! Bonded Contact
keyo,cid,4,2               ! Rigid CERIG style load
keyo,cid,2,2               ! MPC style contact
eblock,10,,,16

"""
