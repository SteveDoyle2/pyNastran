"""Main XDB class"""
import os
import struct

from pyNastran.dev.xdb.xdb_object import XDB_obj
from pyNastran.dev.xdb.debug_output import debug_output


def read_xdb(xdb_filename, etype, nsubcases=1, npload4s=1, debug=False, log=None):
    """function interface to the XDB class"""
    xdb = XDB(debug=debug, log=log)
    xdb.read_xdb(xdb_filename, etype, nsubcases, npload4s)
    return xdb

class XDB:
    def __init__(self, debug=False, log=None):
        self.n = 0
        self._endian = '<'
        self.debug = debug

    def read_xdb(self, xdb_filename, etype, nsubcases, npload4s):
        """reads an *.xdb file"""
        self.nbytes = os.path.getsize(xdb_filename)
        xdb_objects = []
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
            for i in range(34 + npload4s + nsubcases):  # 26 + npload4s + nsubcases?
                table_name = self.read_table_name()

                xdb_obj = self.read_table_header(table_name, etype, npload4s)
                xdb_objects.append(xdb_obj)

            if self.debug:
                debug_output(xdb_objects)

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
            else:
                raise NotImplementedError(etype)
            #self.show(dn, types='s')
            print('(%s)' % dn)
            data = self.f.read(dn)
            self.n += dn
            #self.show(100, types='s')

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
        """Reading control information"""

        tables_a = [
            b'CQD4', #Connectivity Data
            b'CTR3', #Connectivity Data
            b'CBAR', b'CCON', b'CELAS2', b'CONM2', b'CSTM',
            b'RBE3', b'RBEPOOL',

            b'DDLFORDB',
            b'DISPR', # Displacements real
            b'GRIDX', # Node locations in both reference and analysis CSs
            b'LOADR',

            b'LIMITS',
            b'MAT1', # Material data
            b'PATHINT', # Path Attribute Field (integer)
            b'PATHLINK', # Qualifier data type, the number of values associated with the qualifier, and a pointer to the appropriate path value object
            b'PATHQUAL', # The most global data base object (utilize no attributes). Keyed object is used to obtain the valid list of qualifiers for the database.

            b'PLOAD4',
            b'SPC1',
            b'PRODUCT', # MSC product definition, i.e. NASTRAN
            b'PROJECT', # Project description. The most global data base object (utilize no attributes).
            b'PSHELL', # Information from the Bulk Data user input
            b'SID',
            b'SOLVE',
            b'SPCFR', #SPC forces real
            b'SUBCASE', #This field corresponds to the MSC.Nastran Case Control Section definition of SUBCASE structure
            b'SUBCASES',
            b'SUBCTITL', #Contain the information from the TITLE, SUBTITLE and LABEL statements found in the Case Control Section

            b'EQEXINE', b'EQEXING',

            #Strain Recovery Data (Real)
            b'EBARR',
            b'ECONR',
            b'EELSR',
            b'EQD4R',
            b'ETR3R',

            #Stress Recovery Data (Real)
            b'SQD4R',
            b'STR3R',

            #Force Recovery Data
            b'FBARR',
            b'FCONR',

            b'MPCFR', b'PATHBCD', b'PATHINT', b'FORCE', b'PBAR', b'SBARR',
            b'SUBELEM',
        ]
        if table_name in tables_a:
            dn = 88

        elif table_name in [b'SUPERS']:
            dn = 1684

        #Grid object presence indicator for output data recovery:
        elif table_name in [b'SUBGRID']:
            dn = 148
        else:
            raise NotImplementedError('table_name=%r' % table_name)

        self.f.read(dn)
        self.n += dn

        strings, ints, floats = self.show(dn, types='i')


        # XDB Object Header Parsing
        x_obj = XDB_obj(table_name, ints)

        print('read_table_header table_name=%r (%s)' % (table_name, dn))

        return x_obj


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
        elif table_name in [b'SID', b'PLOAD4', b'EQEXING', b'EQEXINE', b'LIMITS',
                            b'CTR3', b'GRIDX', b'MAT1', b'PSHELL', b'SPC1']:
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
    xdb_filename = r'D:\!Work\access\a101x.xdb'
    #r'D:\!Work\bar1.xdb'
    #os.path.join(model_path, 'support_structure', 'w1000bostat.xdb')
    model = read_xdb(xdb_filename, 'quad', nsubcases=3, npload4s=6, debug=True)


    #xdb_filename = os.path.join(model_path, 'pload4', 'cquad4_1subcase.xdb')
    #model = read_xdb(xdb_filename, 'quad', nsubcases=1, npload4s=6)

if __name__ == '__main__':
    test_ctria3()
