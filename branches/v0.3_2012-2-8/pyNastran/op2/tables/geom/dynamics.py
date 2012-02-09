import os
import sys
import struct
from struct import unpack

from pyNastran.op2.op2Errors import *

class DYNAMICS(object):
    def readTable_DYNAMICS(self):
        self.iTableMap = {
                            (37, 18, 183): self.readFake,
                            (57,   5,123): self.readFake,
                            (107,  1, 86): self.readFake,
                            (207,  2, 87): self.readFake,
                            (307,  3, 85): self.readFake,
                            (308,  8,348): self.readFake,
                            (707,  7,124): self.readFake,
                            (1007,10,125): self.readFake,
                            (1307,13,126): self.readFake,
                            (3107,31,127): self.readFake,
                            (5107,51,131): self.readFake,
                            (5207,52,132): self.readFake,
                            (6207,62,136): self.readFake,
                            (6607,66,137): self.readFake,
                            (7107,71,138): self.readFake,
                            (7207,72,139): self.readFake,
                            (8307,83,142): self.readFake,
                            (2107,21,195): self.readFake,
                            (2207,22,196): self.readFake,
                            
                            
                         }
        self.readRecordTable('DYNAMICS')
