#pylint: disable=C0103,C0111,C0301,W0612,W0613,R0914,C0326,R0201
from six import b
from six.moves import range
from struct import unpack, Struct

from pyNastran.bdf.cards.loads.staticLoads import (FORCE, FORCE1, FORCE2, GRAV,
                                                   MOMENT, MOMENT1, MOMENT2,
                                                   LOAD, PLOAD1, PLOAD2,  #PLOAD3,
                                                   PLOAD4)  # PLOAD3,
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, TEMP, TEMPD

class GEOM3(object):
    def add_thermal_load(self, load):
        raise RuntimeError('this should be overwritten')
    def add_load(self, load):
        raise RuntimeError('this should be overwritten')

    def _read_geom3_4(self, data):
        return self._read_geom_4(self._geom3_map, data)

    def __init__(self):
        self.card_count = {}
        self._geom3_map = {
            (4201, 42,  18): ['FORCE', self._readFORCE],   # record 3
            (4001, 40,  20): ['FORCE1', self._readFORCE1],  # record 4
            (4101, 41,  22): ['FORCE2', self._readFORCE2],  # record 5
            (4401, 44,  26): ['GRAV', self._readGRAV],    # record 7 - buggy
            (4551, 61,  84): ['LOAD', self._readLOAD],        # record 8
            (3709, 37, 331): ['LOADCYH', self._readLOADCYH], # record 9 - not done
            (3609, 36, 188): ['LSEQ', self._readLSEQ],       # record 12 - not done
            (4801, 48,  19): ['MOMENT', self._readMOMENT],    # record 13 - not tested
            (4601, 46,  21): ['MOMENT1', self._readMOMENT1],  # record 14 - not tested
            (4701, 47,  23): ['MOMENT2', self._readMOMENT2],  # record 15 - not tested
            (5101, 51,  24): ['PLOAD', self._readPLOAD],      # record 16 - not done
            (6909, 69, 198): ['PLOAD1', self._readPLOAD1],     # record 17 - buggy
            (6802, 68, 199): ['PLOAD2', self._readPLOAD2],     # record 18 - buggy
            (7109, 81, 255): ['PLOAD3', self._readPLOAD3],   # record 19 - not done
            (7209, 72, 299): ['PLOAD4', self._readPLOAD4],     # record 20 - buggy - g1/g3/g4
            (7309, 73, 351): ['PLOADX1', self._readPLOADX1], # record 22
            (4509, 45, 239): ['QBDY1', self._readQBDY1],     # record 24
            (4909, 49, 240): ['QBDY2', self._readQBDY2],     # record 25
            (2109, 21, 414): ['QBDY3', self._readQBDY3],     # record 26
            (5509, 55, 190): ['RFORCE', self._readRFORCE],   # record 30 - not done
            (5401, 54,  25): ['SLOAD', self._readSLOAD],      # record 31 - not done
            (5701, 57,  27): ['TEMP', self._readTEMP],        # record 32
            (5641, 65,  98): ['TEMPD', self._readTEMPD],      # record 33
            (8409, 84, 204): ['TEMPRB', self._readTEMPRB],   # record 40 - not done
            (8109, 81, 201): ['TEMPP1', self._readTEMPP1],   # record 37 - not done
            (8209, 82, 202): ['TEMPP2', self._readTEMPP2],   # record 38 - not done
            (8309, 83, 203): ['TEMPP3', self._readTEMPP3],   # record 39 - not done
            (8409, 84, 204): ['TEMP4', self._readTEMPP4],     # record 40 - not done
            (2309, 23, 416): ['', self._read_fake],
            (4309, 43, 233): ['', self._read_fake],
            (6609, 66, 9031): ['', self._read_fake],
            (8100, 81, 381): ['', self._read_fake],
            (11302, 113, 600): ['', self._read_fake],
            (11402, 114, 601): ['', self._read_fake],
            (2209, 22, 241): ['', self._read_fake],
            (6409, 64, 9032): ['', self._read_fake],
            (3809, 38, 332): ['', self._read_fake],    # record
            (6209, 62, 390): ['', self._read_fake],    # record
            (10901, 109, 427): ['', self._read_fake],  # record
            (10801, 108, 428): ['', self._read_fake],  # record
            (11329, 113, 9602): ['', self._read_fake],  # record
            (11429, 114, 9603): ['', self._read_fake],  # record
            (11529, 115, 9604): ['', self._read_fake],  # record
            (7002, 70, 254) : ['', self._read_fake],  # record
            (7601, 76, 608) : ['', self._read_fake],  # record
        }

# ACCEL
# ACCEL1

    def _readFORCE(self, data, n):
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        #print("reading FORCE")
        ntotal = 28  # 7*4
        nEntries = (len(data) - n) // ntotal
        s = Struct(b(self._endian + 'iiiffff'))
        for i in range(nEntries):
            out = s.unpack(data[n:n + 28])
            (sid, g, cid, f, n1, n2, n3) = out
            self.binary_debug.write('  FORCE=%s\n' % str(out))
            load = FORCE(None, [sid, g, cid, f, n1, n2, n3])
            self.add_load(load)
            n += 28
        return n

    def _readFORCE1(self, data, n):
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        #print("reading FORCE1")
        ntotal = 20  # 5*4
        s = Struct(b(self._endian + 'iifii'))
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 20]
            out = s.unpack(eData)
            (sid, g, f, n1, n2) = out
            self.binary_debug.write('  FORCE1=%s\n' % str(out))
            load = FORCE1(None, [sid, g, f, n1, n2])
            self.add_load(load)
            n += 20
        self.card_count['FORCE1'] = nEntries
        return n

    def _readFORCE2(self, data, n):
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        #print("reading FORCE2")
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + 'iif4i'))
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            out = s.unpack(data[n:n + 28])
            (sid, g, f, n1, n2, n3, n4) = out
            self.binary_debug.write('  FORCE2=%s\n' % str(out))
            load = FORCE2(None, [sid, g, f, n1, n2, n3, n4])
            self.add_load(load)
            n += 28
        self.card_count['FORCE2'] = nEntries
        return n

# GMLOAD

    def _readGRAV(self, data, n):
        """
        GRAV(4401,44,26) - the marker for Record 7
        """
        #print("reading GRAV")
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + 'ii4fi'))
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 28]
            out = s.unpack(eData)
            (sid, cid, a, n1, n2, n3, mb) = out
            grav = GRAV(None, out)
            self.add_load(grav)
            n += 28
        self.card_count['GRAV'] = nEntries
        return n

    def _readLOAD(self, data, n):
        """
        (4551, 61, 84) - the marker for Record 8
        .. todo:: add object
        """
        #print("reading LOAD")
        ntotal = 16  # 4*4
        nEntries = (len(data) - n) // ntotal
        count = 0
        while (len(data) - n) >= 16:
            eData = data[n:n+16]
            n += 16
            out = unpack('iffi', eData)
            (sid, s, si, l1) = out
            self.binary_debug.write('  LOAD=%s\n' % str(out))
            Si = [si]
            L1 = [l1]
            #print(Si, L1)
            while 1:
                eData = data[n:n+8]
                n += 8
                (si, l1) = unpack('fi', eData)
                siTest, = self.struct_i.unpack(eData[0:4])
                #print(si,siTest, l1)
                #print(type(si))

                if [siTest, l1] == [-1, -1]:
                    break
                Si.append(si)
                L1.append(l1)
                self.binary_debug.write('       [%s,%s]\n' % (si, l1))
                #print(Si, L1)

            data_in = [sid, s, Si, L1]
            load = LOAD(None, data_in)
            self.add_load(load)
            count += 1
            if count > 1000:
                raise RuntimeError('Iteration limit...probably have a bug.')
        self.card_count['LOAD'] = nEntries
        return n

    def _readLOADCYH(self, data, n):
        self.binary_debug.write('skipping LOADCYG in GEOM3\n')
        return n

# LOADCYN
# LOADCYT

    def _readLSEQ(self, data, n):
        self.binary_debug.write('skipping LSEQ in GEOM3\n')
        return n

    def _readMOMENT(self, data, n):
        """
        MOMENT(4801,48,19) - the marker for Record 13
        """
        #print("reading MOMENT")
        ntotal = 28
        s = Struct(b(self._endian + '3i4f'))
        nEntries = (len(data) - n) // 28  # 7*4
        for i in range(nEntries):
            eData = data[n:n + 28]
            out = s.unpack(eData)
            self.binary_debug.write('  MOMENT=%s\n' % str(out))
            (sid, g, cid, m, n1, n2, n3) = out
            load = MOMENT(None, out)
            self.add_load(load)
            n += 28
        self.card_count['MOMENT'] = nEntries
        return n

    def _readMOMENT1(self, data, n):
        """
        MOMENT1(4601,46,21) - the marker for Record 14
        """
        #print("reading MOMENT1")
        ntotal = 20  # 5*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 20]
            out = unpack(b(self._endian + 'iifii'), eData)
            self.binary_debug.write('  MOMENT1=%s\n' % str(out))
            (sid, g, m, n1, n2) = out
            load = MOMENT1(None, out)
            self.add_load(load)
            n += 20
        self.card_count['MOMENT1'] = nEntries
        return n

    def _readMOMENT2(self, data, n):
        """
        MOMENT2(4701,47,23) - the marker for Record 15
        """
        #print("reading MOMENT2")
        ntotal = 28  # 7*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 28]
            out = unpack(b(self._endian + 'iif4i'), eData)
            self.binary_debug.write('  MOMENT2=%s\n' % str(out))
            (sid, g, m, n1, n2, n3, n4) = out

            load = MOMENT2(None, out)
            self.add_load(load)
            n += 28
        self.card_count['MOMENT2'] = nEntries
        return n

    def _readPLOAD(self, data, n):
        return n

    def _readPLOAD1(self, data, n):
        """
        PLOAD2(6802,68,199) - the marker for Record 17
        """
        #print("reading PLOAD1")
        ntotal = 32  # 8*4
        s = Struct(b(self._endian + '4i4f'))
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 32]
            out = s.unpack(eData)
            self.binary_debug.write('  PLOAD1=%s\n' % str(out))
            (sid, eid, Type, scale, x1, p1, x2, p2) = out
            #print("PLOAD1 = ", out)
            load = PLOAD1(None, out)
            self.add_load(load)
            n += 32
        self.card_count['PLOAD1'] = nEntries
        return n

    def _readPLOAD2(self, data, n):
        """
        PLOAD2(6802,68,199) - the marker for Record 18
        """
        #print("reading PLOAD2")
        ntotal = 12  # 3*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 12]
            out = unpack('ifi', eData)
            self.binary_debug.write('  PLOAD2=%s\n' % str(out))
            (sid, p, eid) = out
            load = PLOAD2(None, out)
            self.add_load(load)
            n += 12
        self.card_count['PLOAD2'] = nEntries
        return n

    def _readPLOAD3(self, data, n):
        """
        PLOAD3(7109,71,255) - the marker for Record 19
        """
        #print("reading PLOAD3")
        ntotal = 20  # 5*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 20]
            out = unpack('if3i', eData)
            self.binary_debug.write('  PLOAD3=%s\n' % str(out))
            (sid, p, eid, n1, n2) = out
            load = PLOAD3(None, out)  # undefined
            self.add_load(load)
            n += 20
        self.card_count['PLOAD3'] = nEntries
        return n

    def _readPLOAD4(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20
        """
        #print("reading PLOAD4")
        ntotal = 48  # 13*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 48]
                         #iiffffiiifffi   ssssssssssssssss
            out = unpack('2i4f3i3f', eData)
            self.binary_debug.write('  PLOAD4=%s\n' % str(out))
            (sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3) = out
            #s1,s2,s3,s4,s5,s6,s7,s8,L1,L2,L3,L4,L5,L6,L7,L8
            #sdrlA = s1+s2+s3+s4
            #sdrlB = s5+s6+s7+s8
            #ldirA = L1+L2+L3+L4
            #ldirB = L5+L6+L7+L8
            sdrlA = None
            sdrlB = None
            ldirA = None
            ldirB = None
            load = PLOAD4(None, [sid, eid, [p1, p2, p3, p4], g1, g34,
                                 cid, [n1, n2, n3], sdrlA, sdrlB, ldirA, ldirB])
            self.add_load(load)
            n += 48
        self.card_count['PLOAD4'] = nEntries
        return n

# PLOADX - obsolete
    def _readPLOADX1(self, data, n):
        self.binary_debug.write('skipping PLOADX1 in GEOM3\n')
        return n

# PRESAX

    def _readQBDY1(self, data, n):
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        #print("reading QBDY1")
        ntotal = 12  # 3*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 12]
            out = unpack('ifi', eData)
            self.binary_debug.write('  QBDY1=%s\n' % str(out))
            (sid, q0, eid) = out
            load = QBDY1(None, out)
            self.add_thermal_load(load)
            n += 12
        self.card_count['QBDY1'] = nEntries
        return n

    def _readQBDY2(self, data, n):
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        #print("reading QBDY2")
        ntotal = 40  # 10*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 40]
            out = unpack('ii8f', eData)
            self.binary_debug.write('  QBDY2=%s\n' % str(out))
            (sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
            load = QBDY2(None, out)
            self.add_thermal_load(load)
            n += 40
        self.card_count['QBDY2'] = nEntries
        return n

    def _readQBDY3(self, data, n):
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        #print("reading QBDY3")
        ntotal = 16  # 4*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 16]
            out = unpack('ifii', eData)
            (sid, q0, cntrlnd, eid) = out
            load = QBDY3(None, out)
            self.add_thermal_load(load)
            n += 16
        self.card_count['QBDY3'] = nEntries
        return n

    def _readTEMP(self, data, n):
        """
        TEMP(5701,57,27) - the marker for Record 32
        .. warning:: buggy
        """
        #print("reading TEMP")
        ntotal = 12  # 3*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 12]
            out = unpack('iif', eData)
            self.binary_debug.write('  TEMP=%s\n' % str(out))
            (sid, g, T) = out
            if g < 10000000:
                load = TEMP(None, out)
                self.add_thermal_load(load)
            else:
                self.log.debug('TEMP = %s' % (out))
            n += 12
        self.card_count['TEMP'] = nEntries
        return n

    def _readTEMPD(self, data, n):
        """
        TEMPD(5641,65,98) - the marker for Record 33
        .. todo:: add object
        """
        #print("reading TEMPD")
        ntotal = 8  # 2*4
        nEntries = (len(data) - n) // ntotal
        for i in range(nEntries):
            eData = data[n:n + 8]
            out = unpack('if', eData)
            self.binary_debug.write('  TEMPD=%s\n' % str(out))
            (sid, T) = out
            load = TEMPD(None, 0, out)
            #self.add_thermal_load(load)
            n += 8
        self.card_count['TEMPD'] = nEntries
        return n

# QHBDY
# QVECT
# QVOL

    def _readRFORCE(self, data, n):
        self.binary_debug.write('skipping RFORCE in GEOM3\n')
        return n

    def _readSLOAD(self, data, n):
        self.binary_debug.write('skipping SLOAD in GEOM3\n')
        return n

# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
# TEMPEST
# TEMPF
# TEMP1C

    def _readTEMPP1(self, data, n):
        self.binary_debug.write('skipping TEMPP1 in GEOM3\n')
        return n

    def _readTEMPP2(self, data, n):
        self.binary_debug.write('skipping TEMPP2 in GEOM3\n')
        return n

    def _readTEMPP3(self, data, n):
        self.binary_debug.write('skipping TEMPP3 in GEOM3\n')
        return n

    def _readTEMPP4(self, data, n):
        """
        TEMPP4(4201,42,18) - the marker for Record 40
        """
        return n

    def _readTEMPRB(self, data, n):
        self.binary_debug.write('skipping TEMPRB in GEOM3\n')
        return n

# PFACE
# PEDGE
