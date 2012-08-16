import sys
from struct import unpack

from pyNastran.bdf.cards.loads.loads import DAREA


class DYNAMICS(object):
    def readTable_DYNAMICS(self):
        self.iTableMap = {
            (27, 17, 182): self.readDArea,  # 2
            (37, 18, 183): self.readDelay,  # 3
            (57, 5, 123): self.readDLoad,  # 4
            (77, 19, 184): self.readDPhase,  # 5
            (107, 1, 86): self.readEigb,   # 7
            (207, 2, 87): self.readEigc,   # 8
            (257, 4, 158): self.readEigp,   # 9
            (307, 3, 85): self.readEigr,   # 10
            (308, 8, 348): self.readEigrl,  # 11
            (707, 7, 124): self.readEPoint,  # 12


            (1307, 13, 126): self.readFreq,   # 13
            (1007, 10, 125): self.readFreq1,  # 14
            (1107, 11, 166): self.readFreq2,  # 15
            (1407, 14, 39): self.readFreq3,  # 16
            (1507, 15, 40): self.readFreq4,  # 17
            (1607, 16, 41): self.readFreq5,  # 18

            (3107, 31, 127): self.readFake,
            (5107, 51, 131): self.readRLoad1,  # 26
            (5207, 52, 132): self.readRLoad2,  # 27
            (6207, 62, 136): self.readFake,
            (6607, 66, 137): self.readFake,
            (7107, 71, 138): self.readTLoad1,  # 37
            (7207, 72, 139): self.readTLoad2,  # 38
            (8307, 83, 142): self.readTStep,  # 39
            (2107, 21, 195): self.readFake,
            (2207, 22, 196): self.readFake,
        }
        self.readRecordTable('DYNAMICS')

#ACSRCE (5307,53,379)

    def readDArea(self, data):
        """DAREA(27,17,182) - the marker for Record 2"""
        #print "reading DAREA"
        while len(data) >= 16:  # 4*4
            eData = data[:16]
            data = data[16:]
            out = unpack('iiff', eData)
            #(sid,p,c,a) = out
            darea = DAREA(data=out)
            self.addDArea(darea)
        ###

    def readDelay(self, data):
        """DELAY(37,18,183) - Record 3"""
        self.skippedCardsFile.write('skipping DELAY in DYNAMICS\n')

    def readDLoad(self, data):
        """DLOAD(57,5,123) - Record 4"""
        self.skippedCardsFile.write('skipping DLOAD in DYNAMICS\n')

    def readDPhase(self, data):
        """DPHASE(77,19,184) - Record 5"""
        self.skippedCardsFile.write('skipping DPHASE in DYNAMICS\n')

#DYNRED(4807,48,306)

    def readEigb(self, data):
        """EIGB(107,1,86) - Record 7"""
        self.skippedCardsFile.write('skipping EIGB in DYNAMICS\n')

    def readEigc(self, data):
        """EIGC(207,2,87) - Record 8"""
        self.skippedCardsFile.write('skipping EIGC in DYNAMICS\n')

    def readEigp(self, data):
        """EIGP(257,4,158) - Record 9"""
        self.skippedCardsFile.write('skipping EIGP in DYNAMICS\n')

    def readEigr(self, data):
        """EIGR(307,3,85) - Record 10"""
        self.skippedCardsFile.write('skipping EIGR in DYNAMICS\n')

    def readEigrl(self, data):
        """EIGRL(308,8,348) - Record 11"""
        self.skippedCardsFile.write('skipping EIGRL in DYNAMICS\n')

    def readEPoint(self, data):
        """EPOINT(707,7,124) - Record 12"""
        self.skippedCardsFile.write('skipping EPOINT in DYNAMICS\n')

    def readFreq(self, data):
        """FREQ(1307,13,126) - Record 13"""
        self.skippedCardsFile.write('skipping FREQ in DYNAMICS\n')

    def readFreq1(self, data):
        """FREQ1(1007,10,125) - Record 14"""
        self.skippedCardsFile.write('skipping FREQ1 in DYNAMICS\n')

    def readFreq2(self, data):
        """FREQ2(1107,11,166) - Record 15"""
        self.skippedCardsFile.write('skipping FREQ2 in DYNAMICS\n')

    def readFreq3(self, data):
        """FREQ3(1407,14,39) - Record 16"""
        self.skippedCardsFile.write('skipping FREQ3 in DYNAMICS\n')

    def readFreq4(self, data):
        """FREQ4(1507,15,40) - Record 17"""
        self.skippedCardsFile.write('skipping FREQ4 in DYNAMICS\n')

    def readFreq5(self, data):
        """FREQ5(1607,16,41) - Record 18"""
        self.skippedCardsFile.write('skipping FREQ5 in DYNAMICS\n')

#NLRSFD
#NOLIN1
#NOLIN2
#NOLIN3
#NOLIN4
#RANDPS
#RANDT1

    def readRLoad1(self, data):
        """RLOAD1(5107,51,131) - Record 26"""
        self.skippedCardsFile.write('skipping RLOAD1 in DYNAMICS\n')

    def readRLoad2(self, data):
        """RLOAD2(5107,51,131) - Record 27"""
        self.skippedCardsFile.write('skipping RLOAD2 in DYNAMICS\n')

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

    def readTLoad1(self, data):
        """TLOAD1(7107,71,138) - Record 37"""
        self.skippedCardsFile.write('skipping TLOAD1 in DYNAMICS\n')

    def readTLoad2(self, data):
        """TLOAD2(7207,72,139) - Record 37"""
        self.skippedCardsFile.write('skipping TLOAD2 in DYNAMICS\n')

    def readTStep(self, data):
        """TSTEP(8307,83,142) - Record 38"""
        self.skippedCardsFile.write('skipping TSTEP in DYNAMICS\n')

#UNBALNC
