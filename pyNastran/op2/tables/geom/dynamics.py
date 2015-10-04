from struct import unpack

from pyNastran.bdf.cards.loads.loads import DAREA


class DYNAMICS(object):
    def _read_dynamics_4(self, data):
        return self._read_geom_4(self._dynamics_map, data)

    def __init__(self):
        self.card_count = {}
        self._dynamics_map = {
            (5307,53,379) : ['ACSRCE', self._readFake], # 1
            (27, 17, 182): ['DAREA', self.readDArea],  # 2
            (37, 18, 183): ['DELAY', self.readDelay],  # 3
            (57, 5, 123): ['DLOAD', self.readDLoad],  # 4
            (77, 19, 184): ['DPHASE', self.readDPhase],  # 5
            (107, 1, 86): ['EIGB', self.readEigb],   # 7
            (207, 2, 87): ['EIGC', self.readEigc],   # 8
            (257, 4, 158): ['EIGP', self.readEigp],   # 9
            (307, 3, 85): ['EIGR', self.readEigr],   # 10

            (308, 8, 348): ['EIGRL', self.readEigrl],  # 11
            (707, 7, 124): ['EPOINT', self.readEPoint],  # 12
            (1307, 13, 126): ['FREQ', self.readFreq],   # 13
            (1007, 10, 125): ['FREQ1', self.readFreq1],  # 14
            (1107, 11, 166): ['FREQ2', self.readFreq2],  # 15
            (1407, 14, 39): ['FREQ3', self.readFreq3],  # 16
            (1507, 15, 40): ['FREQ4', self.readFreq4],  # 17
            (1607, 16, 41): ['FREQ5', self.readFreq5],  # 18

            (3707, 37, 556) : ['NLRGAP', self._readFake], # 19

            (3107, 31, 127): ['NONLIN1', self._readFake], # 20
            (3207, 32, 128): ['NONLIN2', self._readFake], # 21
            (3207, 33, 129): ['NONLIN3', self._readFake], # 22
            (3207, 34, 130): ['NONLIN4', self._readFake], # 23
            (2107, 21, 195): ['RANDPS', self._readFake], # 24
            (2207, 22, 196): ['RANDT1', self._readFake], # 25
            (5107, 51, 131): ['RLOAD1', self.readRLoad1],  # 26
            (5207, 52, 132): ['RLOAD2', self.readRLoad2],  # 27
            (8910, 89, 606): ['ROTORB', self._readFake],  # 28
            (8210, 82, 599): ['ROTORD', self._readFake],  # 29
            (8410, 84, 600): ['ROTORG', self._readFake],  # 30
            (5707, 57, 135): ['SEQEP', self._readFake],  # 31
            (6207, 62, 136): ['TF', self._readFake],  # 32
            (6607, 66, 137): ['TIC', self._readFake],  # 33
            (7107, 71, 138): ['TLOAD1', self.readTLoad1],  # 37
            (7207, 72, 139): ['TLOAD2', self.readTLoad2],  # 38
            (8307, 83, 142): ['TSTEP', self.readTStep],  # 39

            (10701, 107, 117) : ['', self._readFake],
            (10801, 108, 242) : ['', self._readFake],
            (3807, 38, 505) : ['', self._readFake],
            (4807, 48, 306) : ['', self._readFake],
            (11001, 110, 310) : ['', self._readFake],
            (10901, 109, 260) : ['', self._readFake],
            (3307, 33, 129) : ['', self._readFake],
            (11101, 111, 368) : ['', self._readFake],
        }

#ACSRCE (5307,53,379)

    def readDArea(self, data, n):
        """DAREA(27,17,182) - the marker for Record 2"""
        #print("reading DAREA")
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DAREA', nentries)
        for i in range(nentries):
            eData = data[n:n+ntotal]
            out = unpack('iiff', eData)
            #(sid,p,c,a) = out
            darea = DAREA(data=out)
            self.add_DAREA(darea)
            n += ntotal
        return n

    def readDelay(self, data, n):
        """DELAY(37,18,183) - Record 3"""
        self.binary_debug.write('skipping DELAY in DYNAMICS\n')
        return len(data)

    def readDLoad(self, data, n):
        """DLOAD(57,5,123) - Record 4"""
        self.binary_debug.write('skipping DLOAD in DYNAMICS\n')
        return len(data)

    def readDPhase(self, data, n):
        """DPHASE(77,19,184) - Record 5"""
        self.binary_debug.write('skipping DPHASE in DYNAMICS\n')
        return len(data)

#DYNRED(4807,48,306)

    def readEigb(self, data, n):
        """EIGB(107,1,86) - Record 7"""
        self.binary_debug.write('skipping EIGB in DYNAMICS\n')
        return len(data)

    def readEigc(self, data, n):
        """EIGC(207,2,87) - Record 8"""
        self.binary_debug.write('skipping EIGC in DYNAMICS\n')
        return len(data)

    def readEigp(self, data, n):
        """EIGP(257,4,158) - Record 9"""
        self.binary_debug.write('skipping EIGP in DYNAMICS\n')
        return len(data)

    def readEigr(self, data, n):
        """EIGR(307,3,85) - Record 10"""
        self.binary_debug.write('skipping EIGR in DYNAMICS\n')
        return len(data)

    def readEigrl(self, data, n):
        """EIGRL(308,8,348) - Record 11"""
        self.binary_debug.write('skipping EIGRL in DYNAMICS\n')
        return len(data)

    def readEPoint(self, data, n):
        """EPOINT(707,7,124) - Record 12"""
        self.binary_debug.write('skipping EPOINT in DYNAMICS\n')
        return len(data)

    def readFreq(self, data, n):
        """FREQ(1307,13,126) - Record 13"""
        self.binary_debug.write('skipping FREQ in DYNAMICS\n')
        return len(data)

    def readFreq1(self, data, n):
        """FREQ1(1007,10,125) - Record 14"""
        self.binary_debug.write('skipping FREQ1 in DYNAMICS\n')
        return len(data)

    def readFreq2(self, data, n):
        """FREQ2(1107,11,166) - Record 15"""
        self.binary_debug.write('skipping FREQ2 in DYNAMICS\n')
        return len(data)

    def readFreq3(self, data, n):
        """FREQ3(1407,14,39) - Record 16"""
        self.binary_debug.write('skipping FREQ3 in DYNAMICS\n')
        return len(data)

    def readFreq4(self, data, n):
        """FREQ4(1507,15,40) - Record 17"""
        self.binary_debug.write('skipping FREQ4 in DYNAMICS\n')
        return len(data)

    def readFreq5(self, data, n):
        """FREQ5(1607,16,41) - Record 18"""
        self.binary_debug.write('skipping FREQ5 in DYNAMICS\n')
        return len(data)

#NLRSFD
#NOLIN1
#NOLIN2
#NOLIN3
#NOLIN4
#RANDPS
#RANDT1

    def readRLoad1(self, data, n):
        """RLOAD1(5107,51,131) - Record 26"""
        self.binary_debug.write('skipping RLOAD1 in DYNAMICS\n')
        return len(data)

    def readRLoad2(self, data, n):
        """RLOAD2(5107,51,131) - Record 27"""
        self.binary_debug.write('skipping RLOAD2 in DYNAMICS\n')
        return len(data)

#
#RLOAD2(5207,52,132)
#RGYRO
#ROTORG
#RSPINR
#RSPINT
#SEQEP(5707,57,135)
#TF
#TIC
#TIC
#TIC3

    def readTLoad1(self, data, n):
        """TLOAD1(7107,71,138) - Record 37"""
        self.binary_debug.write('skipping TLOAD1 in DYNAMICS\n')
        return len(data)

    def readTLoad2(self, data, n):
        """TLOAD2(7207,72,139) - Record 37"""
        self.binary_debug.write('skipping TLOAD2 in DYNAMICS\n')
        return len(data)

    def readTStep(self, data, n):
        """TSTEP(8307,83,142) - Record 38"""
        self.binary_debug.write('skipping TSTEP in DYNAMICS\n')
        return len(data)

#UNBALNC
