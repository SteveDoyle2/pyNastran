"""
defines readers for BDF objects in the OP2 GEOM3/GEOM3S table
"""
#pylint: disable=C0103,C0111,C0301,W0612,W0613,R0914,C0326
from __future__ import print_function
from struct import Struct
import numpy as np

from pyNastran.bdf.cards.loads.static_loads import (
    FORCE, FORCE1, FORCE2, GRAV,
    MOMENT, MOMENT1, MOMENT2,
    LOAD, PLOAD, PLOAD1, PLOAD2,  #PLOAD3,
    PLOAD4, PLOADX1)  # PLOAD3,
from pyNastran.bdf.cards.loads.loads import LSEQ, SLOAD, RFORCE #, DAREA, RANDPS, RFORCE1, LOADCYN
from pyNastran.bdf.cards.thermal.loads import (
    QBDY1, QBDY2, QBDY3, TEMP, TEMPD, TEMPP1, QVOL)
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class GEOM3(GeomCommon):
    """defines methods for reading op2 loads"""

    def _read_geom3_4(self, data, ndata):
        return self._read_geom_4(self._geom3_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._geom3_map = {
            (11302, 113, 600): ['ACCEL', self._read_accel],    # record 2 - not done
            (11402, 114, 601): ['ACCEL1', self._read_accel1],  # record 3 - not done
            (4201, 42,  18): ['FORCE', self._read_force],      # record 3
            (4001, 40,  20): ['FORCE1', self._read_force1],    # record 4
            (4101, 41,  22): ['FORCE2', self._read_force2],    # record 5
            (4401, 44,  26): ['GRAV', self._read_grav],        # record 7 - buggy
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
            (7001, 70, 278) : ['PLOADX', self._read_ploadx],
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
            (2309, 23, 416): ['QVOL', self._read_qvol],
            (4309, 43, 233): ['QHBDY', self._read_qhbdy],
            (6609, 66, 9031): ['PEDGE', self._read_pedge],
            (8100, 81, 381): ['CHACAB', self._read_fake],
            (2209, 22, 241): ['QVECT', self._read_qvect],
            (6409, 64, 9032): ['PFACE', self._read_pface],
            (3809, 38, 332): ['LOADCYN', self._read_loadcyn],    # record
            (6209, 62, 390): ['TEMPF', self._read_tempf],    # record
            (10901, 109, 427): ['', self._read_fake],  # record
            (10801, 108, 428): ['', self._read_fake],  # record
            (11329, 113, 9602): ['', self._read_fake],  # record
            (11429, 114, 9603): ['', self._read_fake],  # record
            (11529, 115, 9604): ['', self._read_fake],  # record
            (7002, 70, 254) : ['BOLTFOR', self._read_boltfor],  # record
            (7601, 76, 608) : ['BOLTLD', self._read_boltld],  # record

            # ???
            (6701,67,978): ['PLOADE1', self._read_fake],  # record

            # nx-specific
            (3909, 39, 333): ['LOADCYT', self._read_fake],  # record
        }

    def _read_accel(self, data, n):
        """ACCEL"""
        self.log.info('skipping ACCEL in GEOM3')
        return len(data)

    def _read_accel1(self, data, n):
        """
        ACCEL1(7401,74,601)

        1 SID    I Load set identification number
        2 CID    I Coordinate system identification number
        3 A     RS Acceleration vector scale factor
        4 N(3)  RS Components of a vector coordinate system defined by CID
        7 GRIDID I Grid ID or THRU or BY code
        Words 7 repeats until (-1) occurs.
        NX/MSC
        """
        ntotal = 28  # 7*4
        ints = np.frombuffer(data[n:], dtype='int32').copy()
        floats = np.frombuffer(data[n:], dtype='float32').copy()
        strings = np.frombuffer(data[n:], dtype='|S4')
        i_minus_1s = np.where(ints == -1)[0]

        i0 = 0
        #self.show_data(data[n:])
        for i_minus_1 in i_minus_1s:
            sid = ints[i0]
            cid = ints[i0 + 1]
            scale = floats[i0 + 2]
            n1 = floats[i0 + 3]
            n2 = floats[i0 + 4]
            n3 = floats[i0 + 5]
            nids = [
                strings[i].decode('utf8').strip() if strings[i] in [b'THRU', b'BY  '] else ints[i]
                for i in range(i0+6, i_minus_1)]
            assert nids[-1] > 0
            if self.is_debug_file:
                self.binary_debug.write('  ACCEL1=%s\n' % str([sid, cid, scale, n1, n2, n3, nids]))
            accel = self.add_accel1(sid, scale, [n1, n2, n3], nids, cid=cid)
            accel.validate()
            i0 = i_minus_1 + 1
        return len(data)

    def _read_force(self, data, n):
        """
        FORCE(4201,42,18) - the marker for Record 3
        """
        ntotal = 28  # 7*4
        nentries = (len(data) - n) // ntotal
        s = Struct(self._endian + b'iiiffff')
        for i in range(nentries):
            out = s.unpack(data[n:n + 28])
            (sid, g, cid, f, n1, n2, n3) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE=%s\n' % str(out))
            force = FORCE(sid, g, f, cid=cid, xyz=np.array([n1, n2, n3]))
            self._add_load_object(force)
            n += 28
        self.card_count['FORCE'] = nentries
        return n

    def _read_force1(self, data, n):
        """
        FORCE1(4001,40,20) - the marker for Record 4
        """
        ntotal = 20  # 5*4
        s = Struct(self._endian + b'iifii')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            (sid, g, f, n1, n2) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE1=%s\n' % str(out))
            load = FORCE1.add_op2_data([sid, g, f, n1, n2])
            self._add_load_object(load)
            n += 20
        self.card_count['FORCE1'] = nentries
        return n

    def _read_force2(self, data, n):
        """
        FORCE2(4101,41,22) - the marker for Record 5
        """
        ntotal = 28  # 7*4
        s = Struct(self._endian + b'iif4i')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = s.unpack(data[n:n + 28])
            (sid, g, f, n1, n2, n3, n4) = out
            if self.is_debug_file:
                self.binary_debug.write('  FORCE2=%s\n' % str(out))
            load = FORCE2.add_op2_data([sid, g, f, n1, n2, n3, n4])
            self._add_load_object(load)
            n += 28
        self.card_count['FORCE2'] = nentries
        return n

    def _read_gmload(self, data, n):
        """GMLOAD"""
        self.log.info('skipping GMLOAD in GEOM3')
        return len(data)

    def _read_grav(self, data, n):
        """
        GRAV(4401,44,26) - the marker for Record 7

        Word Name Type Description
        1 SID I Load set identification number
        2 CID I Coordinate system identification number
        3 A RS Acceleration vector scale factor
        4 N(3) RS Components of a vector coordinate system defined by CID
        7 MB I Bulk Data Section with CID definition: -1=main, 0=partitioned
        """
        ntotal = 28  # 7*4
        s = Struct(self._endian + b'ii4fi')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 28]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  GRAV=%s\n' % str(out))
            (sid, cid, a, n1, n2, n3, mb) = out
            grav = GRAV.add_op2_data(out)
            self._add_load_object(grav)
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
        struct_i2fi = Struct('iffi')
        struct_fi = Struct('fi')
        while (len(data) - n) >= 16:
            edata = data[n:n+16]
            n += 16
            out = struct_i2fi.unpack(edata)
            (sid, s, si, l1) = out
            if self.is_debug_file:
                self.binary_debug.write('  LOAD=%s\n' % str(out))
            Si = [si]
            L1 = [l1]
            #print(Si, L1)
            while 1:
                edata = data[n:n+8]
                n += 8
                (si, l1) = struct_fi.unpack(edata)
                siTest, = self.struct_i.unpack(edata[0:4])
                #print(si, siTest, l1)

                if [siTest, l1] == [-1, -1]:
                    break
                Si.append(si)
                L1.append(l1)
                if self.is_debug_file:
                    self.binary_debug.write('       [%s,%s]\n' % (si, l1))
                #print(Si, L1)

            data_in = [sid, s, Si, L1]
            load = LOAD(sid, s, Si, L1)
            self._add_load_combination_object(load)
            count += 1
            if count > 1000:
                raise RuntimeError('Iteration limit...probably have a bug.')
        self.card_count['LOAD'] = nentries
        return n

    def _read_loadcyh(self, data, n):
        """LOADCYH"""
        self.log.info('skipping LOADCYH in GEOM3')
        if self.is_debug_file:
            self.binary_debug.write('skipping LOADCYH in GEOM3\n')
        return len(data)

    def _read_loadcyn(self, data, n):
        """LOADCYN"""
        self.log.info('skipping LOADCYN in GEOM3')
        return len(data)

    def _read_loadcyt(self, data, n):
        """LOADCYT"""
        self.log.info('skipping LOADCYT in GEOM3')
        return len(data)

    def _read_lseq(self, data, n):
        ntotal = 20  # 5*4
        struct_5i = Struct(self._endian + b'5i')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            out = struct_5i.unpack(data[n:n + ntotal])
            (sid, darea, load_id, temperature_id, undef) = out
            if self.is_debug_file:
                self.binary_debug.write('  LSEQ=%s\n' % str(out))
            load = LSEQ.add_op2_data(out)
            n += ntotal
        self.card_count['LSEQ'] = nentries
        return n

    def _read_moment(self, data, n):
        """
        MOMENT(4801,48,19) - the marker for Record 13
        """
        ntotal = 28
        s = Struct(self._endian + b'3i4f')
        nentries = (len(data) - n) // 28  # 7*4
        for i in range(nentries):
            edata = data[n:n + 28]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT=%s\n' % str(out))
            (sid, g, cid, m, n1, n2, n3) = out
            load = MOMENT.add_op2_data(out)
            self._add_load_object(load)
            n += 28
        self.card_count['MOMENT'] = nentries
        return n

    def _read_moment1(self, data, n):
        """
        MOMENT1(4601,46,21) - the marker for Record 14
        """
        ntotal = 20  # 5*4
        nentries = (len(data) - n) // ntotal
        s = Struct(self._endian + b'iifii')
        for i in range(nentries):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT1=%s\n' % str(out))
            (sid, g, m, n1, n2) = out
            load = MOMENT1.add_op2_data(out)
            self._add_load_object(load)
            n += 20
        self.card_count['MOMENT1'] = nentries
        return n

    def _read_moment2(self, data, n):
        """
        MOMENT2(4701,47,23) - the marker for Record 15
        """
        ntotal = 28  # 7*4
        nentries = (len(data) - n) // ntotal
        s = Struct(self._endian + b'iif4i')
        for i in range(nentries):
            edata = data[n:n + 28]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  MOMENT2=%s\n' % str(out))
            (sid, g, m, n1, n2, n3, n4) = out

            load = MOMENT2.add_op2_data(out)
            self._add_load_object(load)
            n += 28
        self.card_count['MOMENT2'] = nentries
        return n

    def _read_pload(self, data, n):
        """
        PLOAD(5101,51,24)
        """
        ntotal = 24  # 6*4
        s = Struct(self._endian + b'i f 4i')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 24]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD=%s\n' % str(out))
            (sid, pressure, n1, n2, n3, n4) = out
            load = PLOAD.add_op2_data(out)
            self._add_load_object(load)
            n += 24
        self.card_count['PLOAD1'] = nentries
        return n

    def _read_pload1(self, data, n):
        """
        PLOAD1(????) - the marker for Record 17
        """
        ntotal = 32  # 8*4
        s = Struct(self._endian + b'4i4f')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n + 32]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD1=%s\n' % str(out))
            (sid, eid, load_type, scale, x1, p1, x2, p2) = out
            load = PLOAD1.add_op2_data(out)
            self._add_load_object(load)
            n += 32
        self.card_count['PLOAD1'] = nentries
        return n

    def _read_pload2(self, data, n):
        """
        PLOAD2(6802,68,199) - the marker for Record 18
        """
        ntotal = 12  # 3*4
        nentries = (len(data) - n) // ntotal
        struct_ifi = Struct('ifi')
        for i in range(nentries):
            edata = data[n:n + 12]
            out = struct_ifi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD2=%s\n' % str(out))
            (sid, p, eid) = out
            load = PLOAD2.add_op2_data(out)
            self._add_load_object(load)
            n += 12
        self.card_count['PLOAD2'] = nentries
        return n

    def _read_pload3(self, data, n):
        """
        PLOAD3(7109,71,255) - the marker for Record 19
        """
        ntotal = 20  # 5*4
        nentries = (len(data) - n) // ntotal
        s = Struct('if3i')
        for i in range(nentries):
            edata = data[n:n + 20]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD3=%s\n' % str(out))
            (sid, p, eid, n1, n2) = out
            load = PLOAD3.add_op2_data(out)  # undefined
            self._add_load_object(load)
            n += 20
        self.card_count['PLOAD3'] = nentries
        return n

    def _read_rbar(self, data, n):
        """RBAR(6601,66,292) - Record 22"""
        n = self._read_dual_card(data, n, self._read_rbar_nx, self._read_rbar_msc,
                                 'RBAR', self._add_rigid_element_object)
        return n

    def _read_pload4(self, data, n):
        """PLOAD4(7209,72,299) - the marker for Record 20"""
        n = self._read_dual_card(data, n, self._read_pload4_nx, self._read_pload4_msc,
                                 'PLOAD4', self._add_load_object)
        return n

    def _read_pload4_msc(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20

        Word Name Type Description
        1 SID          I Load set identification number
        2 EID          I Element identification number
        3 P(4)        RS Pressures
        7 G1           I Grid point identification number at a corner of the face
        8 G34          I Grid point identification number at a diagonal from G1 or CTETRA corner
        9  CID         I Coordinate system identification number
        10 N(3)       RS Components of a vector coordinate system defined by CID
        13 SDRL(2) CHAR4 Load set on element SURF or LINE
        15 LDIR(2) CHAR4 Load direction
        """
        ntotal = 64  # 16*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        loads = []
        s = Struct('2i 4f 3i 3f 8s 8s')
        for i in range(nentries):
            edata = data[n:n + 64]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD4=%s\n' % str(out))
            (sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3, surf_or_line, line_load_dir) = out

            surf_or_line = surf_or_line.rstrip().decode('latin1')
            line_load_dir = line_load_dir.rstrip().decode('latin1')
            if line_load_dir == '':
                # TODO: not 100%
                line_load_dir = 'NORM'

            # forces NX pload4 function to get called if it should be
            assert surf_or_line in ['SURF', 'LINE']
            assert line_load_dir in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM'], 'line_load_dir=%r' % line_load_dir

            load = PLOAD4.add_op2_data(
                [sid, eid, [p1, p2, p3, p4], g1, g34,
                 cid, [n1, n2, n3], surf_or_line, line_load_dir])
            load.validate()
            loads.append(load)
            n += ntotal
        self.card_count['PLOAD4'] = nentries
        return n, loads

    def _read_pload4_nx(self, data, n):  ## inconsistent with DMAP
        """
        PLOAD4(7209,72,299) - the marker for Record 20

        Word Name Type Description
        1 SID          I Load set identification number
        2 EID          I Element identification number
        3 P(4)        RS Pressures
        7 G1           I Grid point identification number at a corner of the face
        8 G34          I Grid point identification number at a diagonal from G1 or CTETRA corner
        9  CID         I Coordinate system identification number
        10 N(3)       RS Components of a vector coordinate system defined by CID
        """
        ntotal = 48  # 12*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        loads = []
        s = Struct('2i 4f 3i 3f')
        for i in range(nentries):
            edata = data[n:n + 48]
            out = s.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOAD4=%s\n' % str(out))
            (sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3) = out

            surf_or_line = None
            line_load_dir = None
            load = PLOAD4.add_op2_data(
                [sid, eid, [p1, p2, p3, p4], g1, g34,
                 cid, [n1, n2, n3], surf_or_line, line_load_dir])
            load.validate()
            loads.append(load)
            n += 48
        self.card_count['PLOAD4'] = nentries
        return n, loads

    def _read_ploadx(self, data, n):
        self.log.info('skipping PLOADX in GEOM3')
        return len(data)

    def _read_ploadx1(self, data, n):
        ntotal = 28  # 7*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'2i2f iif')
        for i in range(nentries):
            edata = data[n:n + 28]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  PLOADX1=%s\n' % str(out))
            load = PLOADX1.add_op2_data(out)
            self._add_load_object(load)
            n += 28
        self.card_count['PLOADX1'] = nentries
        return n

# PRESAX

    def _read_qbdy1(self, data, n):
        """
        QBDY1(4509,45,239) - the marker for Record 24
        """
        ntotal = 12  # 3*4
        nentries = (len(data) - n) // ntotal
        struct_ifi = Struct('ifi')
        for i in range(nentries):
            edata = data[n:n + 12]
            out = struct_ifi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  QBDY1=%s\n' % str(out))
            (sid, q0, eid) = out
            load = QBDY1.add_op2_data(out)
            self._add_thermal_load_object(load)
            n += 12
        self.card_count['QBDY1'] = nentries
        return n

    def _read_qbdy2(self, data, n):
        """
        QBDY2(4909,49,240) - the marker for Record 25
        """
        ntotal = 40  # 10*4
        nentries = (len(data) - n) // ntotal
        struct_2i8f = Struct('ii8f')
        for i in range(nentries):
            edata = data[n:n + 40]
            out = struct_2i8f.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  QBDY2=%s\n' % str(out))
            (sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
            load = QBDY2.add_op2_data(out)
            self._add_thermal_load_object(load)
            n += 40
        self.card_count['QBDY2'] = nentries
        return n

    def _read_qbdy3(self, data, n):
        """
        QBDY3(2109,21,414) - the marker for Record 26
        """
        ntotal = 16  # 4*4
        nentries = (len(data) - n) // ntotal
        struct_if2i = Struct('ifii')
        for i in range(nentries):
            edata = data[n:n + 16]
            out = struct_if2i.unpack(edata)
            (sid, q0, cntrlnd, eid) = out
            load = QBDY3.add_op2_data(out)
            self._add_thermal_load_object(load)
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
        struct_2if = Struct('iif')
        for i in range(nentries):
            edata = data[n:n + 12]
            out = struct_2if.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  TEMP=%s\n' % str(out))
            (sid, g, T) = out
            if g < 10000000:
                load = TEMP.add_op2_data(out)
                self._add_thermal_load_object(load)
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
        struct_if = Struct('if')
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struct_if.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  TEMPD=%s\n' % str(out))
            (sid, T) = out
            load = TEMPD.add_op2_data(out)
            #self.add_thermal_load(load)
            n += ntotal
        self.card_count['TEMPD'] = nentries
        return n

    def _read_qhbdy(self, data, n):
        self.log.info('skipping QHBDY in GEOM3')
        return len(data)

    def _read_qvect(self, data, n):
        self.log.info('skipping QVECT in GEOM3')
        return len(data)

    def _read_qvol(self, data, n):
        """
        Record 30 -- QVOL(2309,23,416)

        1 SID      I Load set identification number
        2 QVOL    RS Power input per unit volume produced by a
                     conduction element
        3 CNTRLND  I Control point used for controlling heat
                     generation
        4 EID      I Element identification number
        """
        ntotal =  16  # 4*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'if2i')
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  QVOL=%s\n' % str(out))
            #(sid, qvol, cntrlnd, eid) = out
            load = QVOL.add_op2_data(out)
            self._add_load_object(load)
            n += ntotal
        self.card_count['TEMPP1'] = nentries
        return n

    def _read_rforce(self, data, n):
        ntotal =  40  # 10*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'3i 4f ifi')
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  RFORCE=%s\n' % str(out))
            #(sid, nid, scale_factor) = out
            load = RFORCE.add_op2_data(out)
            self._add_load_object(load)
            n += ntotal
        self.card_count['RFORCE'] = nentries
        return n

    def _read_sload(self, data, n):
        ntotal =  12  # 3*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'2i f')
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  SLOAD=%s\n' % str(out))
            #(sid, nid, scale_factor) = out
            load = SLOAD.add_op2_data(out)
            self._add_load_object(load)
            n += ntotal
        self.card_count['SLOAD'] = nentries
        return n

# TEMP(5701,57,27) # 32
# TEMPD(5641,65,98) # 33
# TEMPEST
    def _read_tempf(self, data, n):
        self.log.info('skipping TEMPF in GEOM3')
        return len(data)
# TEMP1C

    def _read_tempp1(self, data, n):
        ntotal =  24  # 6*4
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'2i 4f')
        for i in range(nentries):
            edata = data[n:n + ntotal]
            out = struc.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  TEMPP1=%s\n' % str(out))
            #(sid, nid, scale_factor) = out
            load = TEMPP1.add_op2_data(out)
            self._add_load_object(load)
            n += ntotal
        self.card_count['TEMPP1'] = nentries
        return n

    def _read_tempp2(self, data, n):
        self.log.info('skipping TEMPP2 in GEOM3')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPP2 in GEOM3\n')
        return len(data)

    def _read_tempp3(self, data, n):
        self.log.info('skipping TEMPP3 in GEOM3')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPP3 in GEOM3\n')
        return len(data)

    def _read_tempp4(self, data, n):
        """
        TEMPP4(4201,42,18) - the marker for Record 40
        """
        self.log.info('skipping TEMPP4 in GEOM3')
        return len(data)

    def _read_temprb(self, data, n):
        self.log.info('skipping TEMPRB in GEOM3')
        if self.is_debug_file:
            self.binary_debug.write('skipping TEMPRB in GEOM3\n')
        return len(data)

    def _read_pface(self, data, n):
        self.log.info('skipping PFACE in GEOM3')
        return len(data)

    def _read_pedge(self, data, n):
        self.log.info('skipping PEDGE in GEOM3')
        return len(data)

    def _read_boltfor(self, data, n):
        self.log.info('skipping BOLTFOR in GEOM3')
        return len(data)

    def _read_boltld(self, data, n):
        self.log.info('skipping BOLTLD in GEOM3')
        return len(data)
