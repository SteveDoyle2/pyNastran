#pylint: disable=C0103,C0111,C0301,W0612,W0613,R0914,C0326
from struct import unpack, Struct
from six import b
from six.moves import range
import numpy as np

from pyNastran.bdf.cards.loads.static_loads import (
    FORCE, FORCE1, FORCE2, GRAV,
    MOMENT, MOMENT1, MOMENT2,
    LOAD, PLOAD, PLOAD1, PLOAD2,  #PLOAD3,
    PLOAD4)  # PLOAD3,
from pyNastran.bdf.cards.loads.loads import LSEQ, SLOAD #, DAREA, RANDPS, RFORCE, RFORCE1, LOADCYN
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, TEMP, TEMPD
from pyNastran.op2.tables.geom.geom_common import GeomCommon

class GEOM3(GeomCommon):
    """defines methods for reading op2 loads"""

    def _read_geom3_4(self, data, ndata):
        return self._read_geom_4(self._geom3_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom3_map = {
            (11302, 113, 600): ['ACCEL', self._read_fake], # record 2 - not done
            (11402, 114, 601): ['ACCEL1', self._read_fake], # record 3 - not done
            (4201, 42,  18): ['FORCE', self._read_force],    # record 3
            (4001, 40,  20): ['FORCE1', self._read_force1],  # record 4
            (4101, 41,  22): ['FORCE2', self._read_force2],  # record 5
            (4401, 44,  26): ['GRAV', self._read_grav],      # record 7 - buggy
            (4551, 61,  84): ['LOAD', self._read_load],        # record 8
            (3709, 37, 331): ['LOADCYH', self._read_loadcyh],  # record 9 - not done
            (3609, 36, 188): ['LSEQ', self._read_lseq],        # record 12 - not done
            (4801, 48,  19): ['MOMENT', self._read_moment],    # record 13 - not tested
            (4601, 46,  21): ['MOMENT1', self._read_moment1],  # record 14 - not tested
            (4701, 47,  23): ['MOMENT2', self._read_moment2],  # record 15 - not tested
            (5101, 51,  24): ['PLOAD', self._read_pload],      # record 16 - not done
            (6909, 69, 198): ['PLOAD1', self._read_pload1],    # record 17 - buggy
            (6802, 68, 199): ['PLOAD2', self._read_pload2],    # record 18 - buggy
            (7109, 81, 255): ['PLOAD3', self._read_pload3],    # record 19 - not done
            (7209, 72, 299): ['PLOAD4', self._read_pload4],    # record 20 - buggy - g1/g3/g4
            (7309, 73, 351): ['PLOADX1', self._read_ploadx1],  # record 22
            (4509, 45, 239): ['QBDY1', self._read_qbdy1],      # record 24
            (4909, 49, 240): ['QBDY2', self._read_qbdy2],      # record 25
            (2109, 21, 414): ['QBDY3', self._read_qbdy3],      # record 26
            (5509, 55, 190): ['RFORCE', self._read_rforce],    # record 30 - not done
            (5401, 54,  25): ['SLOAD', self._read_sload],      # record 31 - not done
            (5701, 57,  27): ['TEMP', self._read_temp],        # record 32
            (5641, 65,  98): ['TEMPD', self._read_tempd],      # record 33
            (8409, 84, 204): ['TEMPRB', self._read_temprb],    # record 40 - not done
            (8109, 81, 201): ['TEMPP1', self._read_tempp1],    # record 37 - not done
            (8209, 82, 202): ['TEMPP2', self._read_tempp2],    # record 38 - not done
            (8309, 83, 203): ['TEMPP3', self._read_tempp3],    # record 39 - not done
            (8409, 84, 204): ['TEMP4', self._read_tempp4],     # record 40 - not done
            (2309, 23, 416): ['QVOL', self._read_fake],
            (4309, 43, 233): ['QHBDY', self._read_fake],
            (6609, 66, 9031): ['PEDGE', self._read_fake],
            (8100, 81, 381): ['CHACAB', self._read_fake],
            (2209, 22, 241): ['QVECT', self._read_fake],
            (6409, 64, 9032): ['PFACE', self._read_fake],
            (3809, 38, 332): ['LOADCYN', self._read_fake],    # record
            (6209, 62, 390): ['TEMPF', self._read_fake],    # record
            (10901, 109, 427): ['', self._read_fake],  # record
            (10801, 108, 428): ['', self._read_fake],  # record
            (11329, 113, 9602): ['', self._read_fake],  # record
            (11429, 114, 9603): ['', self._read_fake],  # record
            (11529, 115, 9604): ['', self._read_fake],  # record
            (7002, 70, 254) : ['BOLTFOR', self._read_fake],  # record
            (7601, 76, 608) : ['BOLTLD', self._read_fake],  # record
        }

# ACCEL
# ACCEL1

    def _read_force(self, data, n):
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        ntotal = 28  # 7*4
        nentries = (len(data) - n) // ntotal
        s = Struct(b(self._endian + 'iiiffff'))
        for i in range(nentries):
            out = s.unpack(data[n:n + 28])
            (sid, g, cid, f, n1, n2, n3) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE=%s\n' % str(out))
            load = FORCE(sid, g, cid, f, np.array([n1, n2, n3]))
            self.add_load(load)
            n += 28
        self.card_count['FORCE'] = nentries
        return n

    def _read_force1(self, data, n):
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        ntotal = 20  # 5*4
        s = Struct(b(self._endian + 'iifii'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            (sid, g, f, n1, n2) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE1=%s\n' % str(out))
            load = FORCE1.add_op2_data([sid, g, f, n1, n2])
            self.add_load(load)
            n += 20
        self.card_count['FORCE1'] = nentries
        return n

    def _read_force2(self, data, n):
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + 'iif4i'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = s.unpack(data[n:n + 28])
            (sid, g, f, n1, n2, n3, n4) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE2=%s\n' % str(out))
            load = FORCE2.add_op2_data([sid, g, f, n1, n2, n3, n4])
            self.add_load(load)
            n += 28
        self.card_count['FORCE2'] = nentries
        return n

# GMLOAD

    def _read_grav(self, data, n):
        """
        GRAV(4401,44,26) - the marker for Record 7
        """
        ntotal = 28  # 7*4
        s = Struct(b(self._endian + 'ii4fi'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 28]
            out = s.unpack(edata)
            (sid, cid, a, n1, n2, n3, mb) = out
            grav = GRAV.add_op2_data(out)
            self.add_load(grav)
            n += 28
        self.card_count['GRAV'] = nentries
        return n

    def _read_load(self, data, n):
        """
        (4551, 61, 84) - the marker for Record 8
        .. todo:: add object
        """
        ntotal = 16  # 4*4
        nentries = (len(data) - n) // ntotal
        count = 0
        while (len(data) - n) >= 16:
            edata = data[n:n+16]
            n += 16
            out = unpack('iffi', edata)
            (sid, s, si, l1) = out
            if self.is_debug_file:
                self.binary_debug.write('  LOAD=%s\n' % str(out))
            Si = [si]
            L1 = [l1]
            #print(Si, L1)
            while 1:
                edata = data[n:n+8]
                n += 8
                (si, l1) = unpack('fi', edata)
                siTest, = self.struct_i.unpack(edata[0:4])
                #print(si,siTest, l1)
                #print(type(si))

                if [siTest, l1] == [-1, -1]:
                    break
                Si.append(si)
                L1.append(l1)
                if self.is_debug_file:
                    self.binary_debug.write('       [%s,%s]\n' % (si, l1))
                #print(Si, L1)

            data_in = [sid, s, Si, L1]
            load = LOAD(sid, s, Si, L1)
            self.add_load(load)
            count += 1
            if count > 1000:
                raise RuntimeError('Iteration limit...probably have a bug.')
        self.card_count['LOAD'] = nentries
        return n

    def _read_loadcyh(self, data, n):
        self.log.debug('skipping LOADCYH in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping LOADCYH in GEOM3\n')
        return len(data)

# LOADCYN
# LOADCYT

    def _read_lseq(self, data, n):
        ntotal = 20  # 5*4
        s = Struct(b(self._endian + '5i'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = s.unpack(data[n:n + ntotal])
            (sid, darea, load_id, temperature_id, undef) = out
            if self.is_debug_file:
                self.binary_debug.write('  LSEQ=%s\n' % str(out))
            load = LSEQ.add_op2_data(out)
            self.add_LSEQ(load)
            n += ntotal
        self.card_count['LSEQ'] = nentries
        return n

    def _read_moment(self, data, n):
        """
        MOMENT(4801,48,19) - the marker for Record 13
        """
        ntotal = 28
        s = Struct(b(self._endian + '3i4f'))
        nentries = (len(data) - n) // 28  # 7*4
        for i in range(nentries):
            edata = data[n:n + 28]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT=%s\n' % str(out))
            (sid, g, cid, m, n1, n2, n3) = out
            load = MOMENT.add_op2_data(out)
            self.add_load(load)
            n += 28
        self.card_count['MOMENT'] = nentries
        return n

    def _read_moment1(self, data, n):
        """
        MOMENT1(4601,46,21) - the marker for Record 14
        """
        ntotal = 20  # 5*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 20]
            out = unpack(b(self._endian + 'iifii'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT1=%s\n' % str(out))
            (sid, g, m, n1, n2) = out
            load = MOMENT1.add_op2_data(out)
            self.add_load(load)
            n += 20
        self.card_count['MOMENT1'] = nentries
        return n

    def _read_moment2(self, data, n):
        """
        MOMENT2(4701,47,23) - the marker for Record 15
        """
        ntotal = 28  # 7*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 28]
            out = unpack(b(self._endian + 'iif4i'), edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT2=%s\n' % str(out))
            (sid, g, m, n1, n2, n3, n4) = out

            load = MOMENT2.add_op2_data(out)
            self.add_load(load)
            n += 28
        self.card_count['MOMENT2'] = nentries
        return n

    def _read_pload(self, data, n):
        """
        PLOAD(5101,51,24)
        """
        ntotal = 24  # 6*4
        s = Struct(b(self._endian + 'i f 4i'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD=%s\n' % str(out))
            (sid, pressure, n1, n2, n3, n4) = out
            load = PLOAD.add_op2_data(out)
            self.add_load(load)
            n += 24
        self.card_count['PLOAD1'] = nentries
        return n

    def _read_pload1(self, data, n):
        """
        PLOAD1(????) - the marker for Record 17
        """
        ntotal = 32  # 8*4
        s = Struct(b(self._endian + '4i4f'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD1=%s\n' % str(out))
            (sid, eid, Type, scale, x1, p1, x2, p2) = out
            #print("PLOAD1 = ", out)
            load = PLOAD1.add_op2_data(out)
            self.add_load(load)
            n += 32
        self.card_count['PLOAD1'] = nentries
        return n

    def _read_pload2(self, data, n):
        """
        PLOAD2(6802,68,199) - the marker for Record 18
        """
        ntotal = 12  # 3*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 12]
            out = unpack('ifi', edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD2=%s\n' % str(out))
            (sid, p, eid) = out
            load = PLOAD2.add_op2_data(out)
            self.add_load(load)
            n += 12
        self.card_count['PLOAD2'] = nentries
        return n

    def _read_pload3(self, data, n):
        """
        PLOAD3(7109,71,255) - the marker for Record 19
        """
        ntotal = 20  # 5*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 20]
            out = unpack('if3i', edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD3=%s\n' % str(out))
            (sid, p, eid, n1, n2) = out
            load = PLOAD3.add_op2_data(out)  # undefined
            self.add_load(load)
            n += 20
        self.card_count['PLOAD3'] = nentries
        return n

    def _read_pload4(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20
        """
        ntotal = 48  # 13*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 48]
                         #iiffffiiifffi   ssssssssssssssss
            out = unpack('2i 4f 3i 3f', edata)
            if self.is_debug_file:
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
            load = PLOAD4.add_op2_data(
                [sid, eid, [p1, p2, p3, p4], g1, g34,
                 cid, [n1, n2, n3], sdrlA, sdrlB, ldirA, ldirB])
            self.add_load(load)
            n += 48
        self.card_count['PLOAD4'] = nentries
        return n

# PLOADX - obsolete
    def _read_ploadx1(self, data, n):
        self.log.debug('skipping PLOADX1 in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping PLOADX1 in GEOM3\n')
        return len(data)

# PRESAX

    def _read_qbdy1(self, data, n):
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        ntotal = 12  # 3*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 12]
            out = unpack('ifi', edata)
            if self.is_debug_file:
                self.binary_debug.write('  QBDY1=%s\n' % str(out))
            (sid, q0, eid) = out
            load = QBDY1.add_op2_data(out)
            self.add_thermal_load(load)
            n += 12
        self.card_count['QBDY1'] = nentries
        return n

    def _read_qbdy2(self, data, n):
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        ntotal = 40  # 10*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 40]
            out = unpack('ii8f', edata)
            if self.is_debug_file:
                self.binary_debug.write('  QBDY2=%s\n' % str(out))
            (sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
            load = QBDY2.add_op2_data(out)
            self.add_thermal_load(load)
            n += 40
        self.card_count['QBDY2'] = nentries
        return n

    def _read_qbdy3(self, data, n):
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        ntotal = 16  # 4*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 16]
            out = unpack('ifii', edata)
            (sid, q0, cntrlnd, eid) = out
            load = QBDY3.add_op2_data(out)
            self.add_thermal_load(load)
            n += 16
        self.card_count['QBDY3'] = nentries
        return n

    def _read_temp(self, data, n):
        """
        TEMP(5701,57,27) - the marker for Record 32
        .. warning:: buggy
        """
        ntotal = 12  # 3*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 12]
            out = unpack('iif', edata)
            if self.is_debug_file:
                self.binary_debug.write('  TEMP=%s\n' % str(out))
            (sid, g, T) = out
            if g < 10000000:
                load = TEMP.add_op2_data(out)
                self.add_thermal_load(load)
            else:
                self.log.debug('TEMP = %s' % (out))
            n += 12
        self.card_count['TEMP'] = nentries
        return n

    def _read_tempd(self, data, n):
        """
        TEMPD(5641,65,98) - the marker for Record 33
        .. todo:: add object
        """
        ntotal = 8  # 2*4
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = unpack('if', edata)
            if self.is_debug_file:
                self.binary_debug.write('  TEMPD=%s\n' % str(out))
            (sid, T) = out
            load = TEMPD(None, 0, out)
            #self.add_thermal_load(load)
            n += ntotal
        self.card_count['TEMPD'] = nentries
        return n

# QHBDY
# QVECT
# QVOL

    def _read_rforce(self, data, n):
        self.log.debug('skipping RFORCE in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping RFORCE in GEOM3\n')
        return len(data)

    def _read_sload(self, data, n):
        ntotal =  12  # 3*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(b(self._endian + '2i f'))
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  SLOAD=%s\n' % str(out))
            #(sid, nid, scale_factor) = out
            load = SLOAD.add_op2_data(out)
            self.add_load(load)
            n += ntotal
        self.card_count['SLOAD'] = nentries
        return n

# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
# TEMPEST
# TEMPF
# TEMP1C

    def _read_tempp1(self, data, n):
        self.log.debug('skipping TEMPP1 in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPP1 in GEOM3\n')
        return len(data)

    def _read_tempp2(self, data, n):
        self.log.debug('skipping TEMPP2 in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPP2 in GEOM3\n')
        return len(data)

    def _read_tempp3(self, data, n):
        self.log.debug('skipping TEMPP3 in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPP3 in GEOM3\n')
        return len(data)

    def _read_tempp4(self, data, n):
        """
        TEMPP4(4201,42,18) - the marker for Record 40
        """
        self.log.debug('skipping TEMPP4 in GEOM3\n')
        return len(data)

    def _read_temprb(self, data, n):
        self.log.debug('skipping TEMPRB in GEOM3\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPRB in GEOM3\n')
        return len(data)

# PFACE
# PEDGE
