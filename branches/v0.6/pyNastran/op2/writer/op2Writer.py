from struct import pack
from pyNastran.op2.writer.oesWriter import Oes1Writer
from pyNastran.op2.writer.ougWriter import Ougv1Writer


class Op2Writer(Ougv1Writer, Oes1Writer):
    hollerith = pack('i', 584)  # assumes a 128 character string

    def pack(self, format, *vals):
        return [str(val) for val in vals]

    def writeStart(self):
        self.write_markers([3, 7])
        self.write_string_block('NASTRAN FORT TAPE ID CODE - ', 28)
        self.write_markers([2, -1])

    def writeOp2(self, op2Name):
        op2 = open(op2Name, 'wb')
        self.writeStart()
        #self.writeGEOM1()
        #self.writeGEOM2()
        #self.writeGEOM3()
        #self.writeGEOM4()

        #self.writeEPT()
        #self.writeMPTS()

        #self.writeOQG1() # spc forces

        #self.writeOGP()
        self.writeOUGV1()  # displacements/temps/heatFlux
        #self.writeOGF1()
        self.writeOES1()  # stress

    def print_header(self, word, nChars):
        self.device_code = 1  # print the OP2...

        self.write_string_block(word, nChars)

        msg += write_markers([-1, 7])
        out = [101, 0, 4136, 0, 0, 0, 1]  # TODO what this is - DMAP -> "no def or month,year,one,one"...huh???
        msg += pack('iiiiiii', *out)

        msg += write_markers([-2, 1, 0])

        # approach_code=1, table_code=1
        self.iTable = -3

    def write_string_block(self, word, nChars):
        """
        'OUG'      - 1 word  = 4 characters
        'OUGV'     - 1 word  = 4 characters
        'OUGV1'    - 2 words = 8 characters
        '12345678' - 2 words = 8 characters
        nWords = round(ceil(len(word)/4.))
        nChars = nWords*4 != len(word)
        just set nChars and dont overthink it too much
        """
        m = self.write_markers([1])
        #word = self.read_string_block()
        ## creating a %Xs - where x is the number of
        ## words*4,left justifying it, and setting the value
        value = '%%-%ss' % (nChars) % (word)

        out = pack('c' * nChars, value)
        out = ''
        return m + out + m

    def packTitle(isubcase):
        Title = self.Title
        subtitle, label = self.iSubcaseNameMap[isubcase]
        titleSubtitleLabel = '%128s%128s%128s' % (Title, subtitle, label)
        msg = pack('c' * 384, list(titleSubtitleLabel)) + self.hollerith
        return msg

    def write_markers(markers):
        """
        takes -5,1,0  -> [4,5,4,  4,1,4,  4,0,4]
        and puts it into binary
        """
        out = []
        for marker in markers:
            out += [4, marker, 4]
        n = len(out)
        return pack('i' * n, *out)

    def combineApproachDeviceCodes(self, approach_code):
        aCode = approach_code * 10 + self.device_code
        return aCode

    def combineTableSortCodes(self, table_code, sort_code):
        tCode = sort_code * 1000 + table_code
        return tCode

    def aCode_tCode(self, approach_code, table_code, sort_code):
        aCode = self.combineApproachDeviceCodes(approach_code)
        tCode = self.combineTableDeviceCodes(table_code, sort_code)
        return (aCode, tCode)
