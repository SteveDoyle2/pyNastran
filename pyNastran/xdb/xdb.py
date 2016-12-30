import os
import struct
from pyNastran.op2.fortran_format import FortranFormat


def read_xdb(xdb_filename, etype, nsubcases=1, npload4s=1, debug=False, log=None):
    xdb = XDB(debug=debug, log=log)
    xdb.read_xdb(xdb_filename, etype, nsubcases, npload4s)
    return xdb

class XDB(FortranFormat):
    def __init__(self, debug=False, log=None):
        FortranFormat.__init__(self)
        self.n = 0
        self._endian = '<'

    def read_xdb(self, xdb_filename, etype, nsubcases, npload4s):
        self.nbytes = os.path.getsize(xdb_filename)
        with open(xdb_filename, mode='rb') as self.f:

            print('(4100)')
            #self.show(80)
            # 6 tri
            #    1024, 16777215, 16777097, 118, 0, 2, 39, 39, 3, 1, 1, 4, 40, 20010, 72501883, 5, 0, 0, 0, 0
            # 1 tri
            #    1024, 16777215, 16777127, 88,  0, 2, 29, 29, 3, 1, 1, 4, 40, 20010, 72501883, 5, 0, 0, 0, 0
            # 2-quad
            #    1024, 16777215, 16777121, 94,  0, 2, 31, 31, 3, 1, 1, 4, 40, 20010, 72501883, 5, 0, 0, 0, 0
            self.f.read(4100)
            self.n += 4100

            # CTR3-----
            table_name = self.read_table_name()
            print('table_name = %r (4100)' % table_name)
            self.f.read(4100)
            self.n += 4100

            # CTR3/CQD4-----
            # DDLFORDB-----
            # DISPR-----
            # EQEXINE-----
            # EQEXING-----
            # GRIDX-----
            # LIMITS-----
            # MAT1-----
            # PATHINT-----
            # PATHLINK-----
            # PATHQUAL-----
            # PLOAD4-----
            # PRODUCT-----
            # SPCFR-----
            # PROJECT-----
            # PSHELL-----
            # SID-----
            # SOLVE-----
            # SPC1-----
            # SPCFR-----
            # SUBCASE-----
            # SUBCASES-----
            # SUBTITL-----
            # SUBGRID-----
            for i in range(32 + npload4s + nsubcases):  # 26 + npload4s + nsubcases?
                table_name = self.read_table_name()
                self.read_table_header(table_name, etype, npload4s)

            # SUPERS-----
            #self.show(1000, types='s')
            table_name = self.read_table_name()
            self.read_table(table_name, etype, nsubcases)

            # CTR3-----
            table_name = self.read_table_name()
            dn = 4092
            data = self.f.read(dn)
            self.n += dn
            #self.read_table(table_name, etype, nsubcases)

            # DDLFORDB-----
            table_name = self.read_table_name()
            print('table_name = %r (24616)' % table_name)

            data = self.f.read(24616)
            self.n += 24616

            # NASTRAN 8-----
            data = self.f.read(12)
            word = struct.unpack('12s', data)
            self.n += 12
            print('word = %r (16)' % word)

            data = self.f.read(16)
            self.n += 16

            # INTEL64 FAMILY 6 MOD                    WINDOWS 7           5-----
            data = self.f.read(64)
            word = struct.unpack('64s', data)
            self.n += 64
            print('word = %r (3944)' % word)

            data = self.f.read(3944)
            self.n += 3944

            # PROJECT-----
            table_name = self.read_table_name()
            print('table_name = %r (4108)' % table_name)

            data = self.f.read(4108)
            self.n += 4108

            # spaces-----
            data = self.f.read(56)
            word = struct.unpack('56s', data)
            self.n += 56
            print('word = %r (8)' % word)

            data = self.f.read(8)
            self.n += 8

            # spaces-----
            data = self.f.read(56)
            word = struct.unpack('56s', data)
            self.n += 56
            print('word = %r (8052)' % word)

            data = self.f.read(8052)
            self.n += 8052

            # PATHQUAL-----
            table_name = self.read_table_name()
            print('table_name = %r (12280)' % table_name)


            data = self.f.read(12280)
            self.n += 12280

            # PATHLINK-----
            table_name = self.read_table_name()
            print('table_name = %r (4108)' % table_name)

            data = self.f.read(4108)
            self.n += 4108

            # MODEL-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # SOLID-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # HIGHQUAL-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # AUXMID-----
            table_name = self.read_table_name()
            print('table_name = %r (8080)' % table_name)

            data = self.f.read(8080)
            self.n += 8080

            # PATHINT-----
            table_name = self.read_table_name()
            print('table_name = %r (12292)' % table_name)

            data = self.f.read(12292)
            self.n += 12292

            # SUPERS-----
            table_name = self.read_table_name()
            print('table_name = %r (4076)' % table_name)

            data = self.f.read(4076)
            self.n += 4076

            # PSHELL-----
            # PSHELL-----
            # GRIDX-----
            # CTR3-----
            # LIMITS-----
            # EXEXING-----
            # EXEXINE-----
            # SID-----
            # PLOAD4-----
            for isubcase in range(3 + nsubcases + npload4s):
                table_name = self.read_table_name()
                self.read_table(table_name, etype, nsubcases)

            # SPC1-----
            # SUBCASES-----
            table_name = self.read_table_name()
            print('table_name = %r (12292)' % table_name)

            data = self.f.read(12292)
            self.n += 12292

            # SUBCASE-----
            table_name = self.read_table_name()
            print('table_name = %r (4076)' % table_name)

            data = self.f.read(4076)
            self.n += 4076

            # SUBGRID-----
            table_name = self.read_table_name()
            print('table_name = %r' % table_name)
            self.read_table(table_name, etype, nsubcases)

            # MODEL-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # SOLID-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # HIGHQUAL-----
            table_name = self.read_table_name()
            print('table_name = %r (20)' % table_name)

            data = self.f.read(20)
            self.n += 20

            # AUXMID-----
            table_name = self.read_table_name()
            print('table_name = %r (3094)' % table_name)

            data = self.f.read(3984)
            self.n += 3984

            # PATHLINK-----
            table_name = self.read_table_name()
            print('table_name = %r (944)' % table_name)

            data = self.f.read(944)
            self.n += 944

            # SUBTITL/SUBCTITL-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (4104)' % table_name)
            data = self.f.read(4104)
            self.n += 4104

            for isubcase in range(nsubcases):
                # SUBCASE 1-----
                data = self.f.read(384)
                table_name, = struct.unpack('384s', data)
                self.n += 384
                print('  table_name = %r (384)' % table_name.strip())

                data = self.f.read(4)
                self.n += 4

            if etype == 'tri':
                if nsubcases == 6:
                    dn = 4908
                elif nsubcases == 1:
                    dn = 6848
                else:
                    raise NotImplementedError(nsubcases)
            elif etype == 'quad':
                if nsubcases == 2:
                    dn = 6460
                else:
                    raise NotImplementedError(nsubcases)
                #self.show(6468, types='s')
                #aaa
            else:
                raise NotImplementedError(etype)
            #self.show(dn, types='s')
            print('(%s)' % dn)
            data = self.f.read(dn)
            self.n += dn
            #self.show(100, types='s')
            #aa

            # SUPERS-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (932)' % table_name)

            data = self.f.read(932)
            self.n += 932

            # DISPR-----
            table_name = self.read_table_name()
            print('table_name = %r (3736)' % table_name)

            data = self.f.read(3736)
            self.n += 3736

            # GRIDX-----
            table_name = self.read_table_name()
            print('table_name = %r (8536)' % table_name)

            data = self.f.read(8536)
            self.n += 8536

            # DISPR-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (3148)' % table_name)

            data = self.f.read(3148)
            self.n += 3148

            # SUBCASE-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (9124)' % table_name)

            data = self.f.read(9124)
            self.n += 9124

            # DISPR-----
            table_name = self.read_table_name()
            print('table_name = %r (3736)' % table_name)

            data = self.f.read(3736)
            self.n += 3736

            # LIMITS-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (8536)' % table_name)

            if etype == 'tri':
                if nsubcases == 6:
                    dn = 8536
                elif nsubcases == 1:
                    dn = 8532
                else:
                    raise NotImplementedError(nsubcases)
            elif etype == 'quad':
                if nsubcases == 2:
                    dn = 8536
                else:
                    raise NotImplementedError(nsubcases)
            else:
                raise NotImplementedError(etype)
            #self.show(dn, types='s')
            data = self.f.read(dn)
            self.n += dn
            assert self.n == self.f.tell()

            #if etype == 'tri':
                #if nsubcases == 1:
                    #self.show(1500, types='s')
                    #return
            #elif etype == 'quad':
                #if nsubcases == 2:
                    #dn = (self.nbytes - self.n) // 5
                    #dn = 3744
                    #self.show(dn, types='s')
                    #print('dn = ', dn)
                    #return

            # DISPR-----
            table_name = self.read_table_name()
            #print('602')
            if table_name is None:
                return
            print('table_name = %r (3736)' % table_name)
            data = self.f.read(3736)
            self.n += 3736

            # EQEXING-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            self.read_table(table_name, etype, nsubcases)

            # EQEXINE-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            self.read_table(table_name, etype, nsubcases)

            # SID-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            self.read_table(table_name, etype, nsubcases)

            for isubcase in range(nsubcases):
                # PLOAD4-----
                table_name = self.read_table_name()
                self.read_table(table_name, etype, nsubcases)

            # SPC1-----
            table_name = self.read_table_name()
            if table_name is None:
                return
            print('table_name = %r (8532)' % table_name)
            data = self.f.read(8532)
            self.n += 8532

            self.show(5000, types='s')
        print('done!')

    def read_table_header(self, table_name, etype, nsubcases):
        if table_name in [b'CQD4', b'DDLFORDB', b'DISPR', b'GRIDX', b'LIMITS',
                          b'MAT1', b'PATHINT', b'PATHLINK', b'PATHQUAL', b'PLOAD4',
                          b'SPC1', b'PRODUCT', b'PROJECT', b'PSHELL', b'SID',
                          b'SOLVE', b'SPCFR', b'SUBCASE', b'SUBCASES', b'SUBCTITL',
                          b'SUBGRID', b'EQEXINE', b'EQEXING', b'SPCFR', b'CTR3',
                          b'CBAR', b'CCON', b'CELAS2', b'CONM2' ,b'CSTM', b'EBARR',
                          b'ECONR', b'EELSR', b'EQD4R', b'ETR3R', b'FBARR', b'FCONR',
                          b'MPCFR', b'PATHBCD', b'PATHINT',
                          ]:
            dn = 88
            self.show(dn, types='i')
            self.f.read(dn)
            self.n += dn
        else:
            raise NotImplementedError('table_name=%r' % table_name)
        print('read_table_header table_name=%r (%s)' % (table_name, dn))

    def read_table(self, table_name, etype, nsubcases):
        #print('read_table...%r' % table_name)
        if table_name == b'SUBCTITL':
            dn = 88
            #self.show(dn + 8, types='s')
        elif table_name == b'SUPERS':
            print('etype=%r, nsubcases=%s' % (etype, nsubcases))
            if etype == 'tri':
                if nsubcases == 6:
                    dn = 436
                else:
                    dn = 1396
            elif etype == 'quad':
                if nsubcases == 2:
                    dn = 1204
                elif nsubcases == 1:
                    dn = 88
                else:
                    raise NotImplementedError(nsubcases)
            else:
                raise NotImplementedError(etype)
            #self.show(dn + 8, types='s')
            #print('etype=%r, nsubcases=%s' % (etype, nsubcases))

        elif table_name == b'SUBGRID':
            if etype == 'tri':
                dn = 7752
                #print('table_name = %r (%s)' % (table_name, dn))
            elif etype == 'quad':
                dn = 7252
                #self.show(dn + 8, types='s')
            else:
                raise NotImplementedError(etype)
            #self.show(dn + 8, types='s')

        elif table_name == b'CQD4':
            dn = 4092
            #self.show(dn + 8, types='s')
        elif table_name == b'CTR3':
            #dn = 4092
            dn = 12280
            #self.show(dn2 + 8, types='s')
            #aaa
        elif table_name in [b'SID', b'PLOAD4', b'EQEXING', b'EQEXINE', b'LIMITS',
                            b'CTR3', b'GRIDX', b'MAT1', b'PSHELL', b'SPC1',
                            ]:
            dn = 12280
        else:
            raise NotImplementedError('table_name=%r' % table_name)
        self.f.read(dn)
        self.n += dn
        print('table_name = %r (%s)' % (table_name, dn))

    def read_table_name(self, rewind=False):
        if self.n == self.nbytes:
            return None
        data = self.f.read(8)
        if len(data) == 0:
            print(self.n, self.nbytes)
            self.nbytes = self.f.tell()
            return None
        if rewind:
            self.f.seek(self.n)
        else:
            self.n += 8
        return struct.unpack(self._endian + '8s', data)[0].rstrip()

def test_ctria3():
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '..', 'models')

    #xdb_filename = os.path.join(model_path, 'pload4', 'ctria3_6subcases.xdb')
    #model = read_xdb(xdb_filename, 'tri', nsubcases=6, npload4s=6)

    #xdb_filename = os.path.join(model_path, 'pload4', 'ctria3.xdb')
    #model = read_xdb(xdb_filename, 'tri', nsubcases=1, npload4s=6)

    #xdb_filename = r'C:\Users\nikita.kalutskiy\Desktop\A318_FEM\Static\w1000bostat.xdb'
    #model = read_xdb(xdb_filename, 'quad', nsubcases=3, npload4s=6)

    xdb_filename = os.path.join(model_path, 'pload4', 'cquad4_1subcase.xdb')
    model = read_xdb(xdb_filename, 'quad', nsubcases=1, npload4s=6)

if __name__ == '__main__':
    test_ctria3()
