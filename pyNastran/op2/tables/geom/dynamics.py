from __future__ import print_function
from struct import unpack, Struct
from six import b
import numpy as np

from pyNastran.bdf.cards.nodes import EPOINTs
from pyNastran.bdf.cards.loads.loads import DAREA
from pyNastran.bdf.cards.methods import EIGB
from pyNastran.bdf.cards.dynamic import FREQ1, TF
from pyNastran.bdf.cards.loads.dloads import RLOAD2
from pyNastran.op2.tables.geom.geom_common import GeomCommon

class DYNAMICS(GeomCommon):
    """defines methods for reading op2 dynamics loads/methods"""
    def _read_dynamics_4(self, data, ndata):
        return self._read_geom_4(self._dynamics_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._dynamics_map = {
            (5307, 53, 379) : ['ACSRCE', self._read_fake], # 1
            (27, 17, 182): ['DAREA', self._read_fake],  # 2

            (37, 18, 183): ['DELAY', self._read_delay],  # 3
            (57, 5, 123): ['DLOAD', self._read_dload],  # 4
            (77, 19, 184): ['DPHASE', self._read_dphase],  # 5
            (107, 1, 86): ['EIGB', self._read_eigb],   # 7
            (207, 2, 87): ['EIGC', self._read_eigc],   # 8
            (257, 4, 158): ['EIGP', self._read_eigp],   # 9
            (307, 3, 85): ['EIGR', self._read_eigr],   # 10

            (308, 8, 348): ['EIGRL', self._read_eigrl],  # 11
            (707, 7, 124): ['EPOINT', self._read_epoint],  # 12
            (1307, 13, 126): ['FREQ', self._read_freq],   # 13
            (1007, 10, 125): ['FREQ1', self._read_freq1],  # 14
            (1107, 11, 166): ['FREQ2', self._read_freq2],  # 15
            (1407, 14, 39): ['FREQ3', self._read_freq3],  # 16
            (1507, 15, 40): ['FREQ4', self._read_freq4],  # 17
            (1607, 16, 41): ['FREQ5', self._read_freq5],  # 18

            (3707, 37, 556) : ['NLRGAP', self._read_fake], # 19

            (3107, 31, 127): ['NONLIN1', self._read_fake], # 20
            (3207, 32, 128): ['NONLIN2', self._read_fake], # 21
            (3207, 33, 129): ['NONLIN3', self._read_fake], # 22
            (3207, 34, 130): ['NONLIN4', self._read_fake], # 23
            (2107, 21, 195): ['RANDPS', self._read_fake], # 24
            (2207, 22, 196): ['RANDT1', self._read_fake], # 25
            (5107, 51, 131): ['RLOAD1', self._read_rload1],  # 26
            (5207, 52, 132): ['RLOAD2', self._read_rload2],  # 27
            (8910, 89, 606): ['ROTORB', self._read_fake],  # 28
            (8210, 82, 599): ['ROTORD', self._read_fake],  # 29
            (8410, 84, 600): ['ROTORG', self._read_fake],  # 30
            (5707, 57, 135): ['SEQEP', self._read_fake],  # 31
            (6207, 62, 136): ['TF', self._read_tf],  # 32
            (6607, 66, 137): ['TIC', self._read_fake],  # 33
            (7107, 71, 138): ['TLOAD1', self._read_tload1],  # 37
            (7207, 72, 139): ['TLOAD2', self._read_tload2],  # 38
            (8307, 83, 142): ['TSTEP', self._read_tstep],  # 39

            (10701, 107, 117) : ['RGYRO', self._read_fake],
            (10801, 108, 242) : ['ROTORG', self._read_fake],
            (3807, 38, 505) : ['NLRSFD', self._read_fake],
            (4807, 48, 306) : ['DYNRED', self._read_fake],
            (11001, 110, 310) : ['RSPINT', self._read_fake],
            (10901, 109, 260) : ['RSPINR', self._read_fake],
            (3307, 33, 129) : ['NONLIN3', self._read_fake],
            (11101, 111, 368) : ['UNBALNC', self._read_fake],
        }

#ACSRCE (5307,53,379)

    def _read_area(self, data, n):
        """DAREA(27,17,182) - the marker for Record 2"""
        #print("reading DAREA")
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DAREA', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iiff', edata)
            #(sid,p,c,a) = out
            darea = DAREA.add_op2_data(data=out)
            self._add_darea_object(darea)
            n += ntotal
        return n

    def _read_delay(self, data, n):
        """DELAY(37,18,183) - Record 3"""
        if self.is_debug_file:
            self.binary_debug.write('skipping DELAY in DYNAMICS\n')
        return len(data)

    def _read_dload(self, data, n):
        """
        DLOAD(57,5,123) - Record 4

        1 SID  I Load set identification number
        2  S   RS Overall scale factor
        3  SI  RS Scale factor i
        4  LI  I Load set identification number i
        Words 3 through 4 repeat until (-1,-1) occurs
        """
        ndata = len(data)
        nfields = (ndata - n) // 4

        datan = data[n:]
        ints = np.fromstring(data[n:], self.idtype)
        floats = np.fromstring(data[n:], self.fdtype)
        iminus1 = np.where(ints == -1)[0]
        delta = iminus1[1:] - iminus1[:-1]
        idelta = np.where(delta == 1)[0]
        iminus1_delta = iminus1[delta] - 1
        istart = 0
        nentries = 0
        for iend in iminus1_delta:
            datai = data[n+istart*4 : n+iend*4]
            sid = ints[istart]
            global_scale = floats[istart + 1]
            #print('  sid=%s global_scale=%s' % (sid, global_scale))
            deltai = iend - istart - 2 # subtract 2 for sid, global scale
            assert deltai % 2 == 0, (self.show_data(data[n+istart*4 : n+iend*4], 'if'))

            scales = []
            load_ids = []
            for iscale in range(deltai // 2):
                scale = floats[istart + 2 + 2*iscale]
                load_id = ints[istart + 3 + 2*iscale]
                scales.append(scale)
                load_ids.append(load_id)
            dload = self.add_dload(sid, global_scale, scales, load_ids)
            istart = iend + 1
            nentries += 1
        self._increase_card_count('DLOAD', nentries)
        return n

    def _read_dphase(self, data, n):
        """DPHASE(77,19,184) - Record 5"""
        if self.is_debug_file:
            self.binary_debug.write('skipping DPHASE in DYNAMICS\n')
        return len(data)

#DYNRED(4807,48,306)

    def _read_eigb(self, data, n):
        """EIGB(107,1,86) - Record 7"""
        #if self.is_debug_file:
            #self.binary_debug.write('skipping EIGB in DYNAMICS\n')
        #return len(data)
        ntotal = 60
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('EIGB', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            #self.show_data(edata[44:])
            # int, 8s, 2f, 3i, i, 8s, 4i
            out = unpack('i8s ff 3i i 8s 4i', edata)
            sid, method, L1, L2, nep, ndp, ndn, dunno, norm, g, c, dunnoA, dunnoB = out
            if self.is_debug_file:
                self.binary_debug.write('EIGB=%s\n' % str(out))
            #print('out = %s' % str(out))
            method = method.strip().decode('latin1')
            norm = norm.strip().decode('latin1')
            eigb = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, g, c)
            self._add_method_object(eigb)
            n += ntotal
        return n

    def _read_eigc(self, data, n):
        """EIGC(207,2,87) - Record 8"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGC in DYNAMICS\n')
        return len(data)

    def _read_eigp(self, data, n):
        """EIGP(257,4,158) - Record 9"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGP in DYNAMICS\n')
        return len(data)

    def _read_eigr(self, data, n):
        """EIGR(307,3,85) - Record 10"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGR in DYNAMICS\n')
        return len(data)

    def _read_eigrl(self, data, n):
        """EIGRL(308,8,348) - Record 11"""
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGRL in DYNAMICS\n')
        return len(data)

    def _read_epoint(self, data, n):
        """EPOINT(707,7,124) - Record 12"""
        npoints = (len(data) - n) // 4
        fmt = b(self._endian + '%ii' % npoints)
        nids = unpack(fmt, data[n:])
        if self.is_debug_file:
            self.binary_debug.write('EPOINT=%s\n' % str(nids))
        epoint = EPOINTs.add_op2_data(list(nids))
        self._add_epoint_object(epoint)
        self.card_count['EPOINT'] = npoints
        self._increase_card_count('EPOINT', count_num=npoints)
        return n

    def _read_freq(self, data, n):
        """FREQ(1307,13,126) - Record 13"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ in DYNAMICS\n')
        return len(data)

    def _read_freq1(self, data, n):
        """FREQ1(1007,10,125) - Record 14

        1 SID I  Set identification number
        2 F1  RS First frequency
        3 DF  RS Frequency increment
        4 NDF I  Number of frequency increments
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iffi', edata)
            sid, f1, df, ndf = out
            if self.is_debug_file:
                self.binary_debug.write('FREQ1=%s\n' % str(out))
            #print('out = %s' % str(out))
            freq = FREQ1(sid, f1, df, ndf=ndf)
            self._add_freq_object(freq)
            n += ntotal
        self._increase_card_count('FREQ1', nentries)
        return n

    def _read_freq2(self, data, n):
        """FREQ2(1107,11,166) - Record 15"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ2 in DYNAMICS\n')
        return len(data)

    def _read_freq3(self, data, n):
        """FREQ3(1407,14,39) - Record 16"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ3 in DYNAMICS\n')
        return len(data)

    def _read_freq4(self, data, n):
        """FREQ4(1507,15,40) - Record 17"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ4 in DYNAMICS\n')
        return len(data)

    def _read_freq5(self, data, n):
        """FREQ5(1607,16,41) - Record 18"""
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ5 in DYNAMICS\n')
        return len(data)

#NLRSFD
#NOLIN1
#NOLIN2
#NOLIN3
#NOLIN4
#RANDPS
#RANDT1

    def _read_rload1(self, data, n):
        """RLOAD1(5107,51,131) - Record 26"""
        if self.is_debug_file:
            self.binary_debug.write('skipping RLOAD1 in DYNAMICS\n')
        return len(data)

    def _read_rload2(self, data, n):
        """
        RLOAD2(5107,51,131) - Record 27

        1 SID     I  Load set identification number
        2 DAREA   I  DAREA Bulk Data entry identification number
        3 DELAYI  I  DELAY Bulk Data entry identification number
        4 DPHASEI I  DPHASE Bulk Data entry identification number
        5 TBI     I  TABLEDi Bulk Data entry identification number for B(f)
        6 TPI     I  TABLEDi Bulk Data entry identification number for Phi(f)
        7 TYPE    I  Nature of the dynamic excitation
        8 DELAYR  RS If DELAYI = 0, constant value for delay
        9 DPHASER RS If DPHASEI = 0, constant value for phase
        10 TBR    RS If TBI = 0, constant value for B(f)
        11 TPR    RS If TPI = 0, constant value for PHI(f)
        """
        ntotal = 44
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('7i 4f', edata)
            sid, darea, delayi, dphasei, tbi, tpi, Type, delayr, dphaser, tbr, tpr = out
            if self.is_debug_file:
                self.binary_debug.write('RLOAD2=%s\n' % str(out))

            tb = tbi
            tp = tpi
            delay = delayi
            dphase = dphasei
            if tbi == 0:
                tb = tbr
            if dphasei == 0:
                dphase = dphaser
            if tpi == 0:
                tp = tpr
            dload = RLOAD2(sid, darea, delay=delay, dphase=dphase, tb=tb, tp=tp,
                           Type=Type, comment='')
            self._add_dload_entry(dload)
            n += ntotal
        self._increase_card_count('RLOAD2', nentries)
        return n

#RLOAD2(5207,52,132)
#RGYRO
#ROTORG
#RSPINR
#RSPINT
#SEQEP(5707,57,135)

    def _read_tf(self, data, n):
        nfields = (len(data) - n) // 4

        # subtract of the header (sid, nid, component, b0, b1, b2)
        # divide by 5 (nid1, component1, a0, a1, a2)
        #nrows = (nfields - 6) // 5
        #print('n=%s nrows=%s' % (n, nrows))
        #print(self.show_data(data))
        #nid1, component1, a0, a1, a2

        ndata = len(data)
        struct1 = Struct(b'3i3f')
        struct2 = Struct(b'2i3f')
        while n < ndata:
            n2 = n + 24 # 20=4*6
            sid, nid, component, b0, b1, b2 = struct1.unpack(data[n:n2])
            if self.is_debug_file:
                self.binary_debug.write('TF header -> %s\n' % ([sid, nid, component, b0, b1, b2]))

            nids = []
            components = []
            a = []
            irow = 0
            while 1:
                n3 = n2 + 20 # 20=4*5
                nid1, component1, a0, a1, a2 = struct2.unpack(data[n2:n3])

                if self.is_debug_file:
                    self.binary_debug.write('  i=%s     -> %s\n' % (
                        irow, [nid1, component1, a0, a1, a2]))

                if nid1 == -1 and component1 == -1:
                    break
                assert nid1 > -1
                assert component1 > -1
                nids.append(nid1)
                components.append(component1)
                a.append([a0, a1, a2])
                n2 = n3
                irow += 1

            tf = TF(sid, nid, component, b0, b1, b2, nids, components, a)
            #if self.is_debug_file:
                #self.binary_debug.write('%s\n' % str(tf))
            self._add_tf_object(tf)
            self._increase_card_count('TF')
            n = n3
        return n

#TIC
#TIC3

    def _read_tload1(self, data, n):
        """TLOAD1(7107,71,138) - Record 37"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TLOAD1 in DYNAMICS\n')
        return len(data)

    def _read_tload2(self, data, n):
        """TLOAD2(7207,72,139) - Record 37"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TLOAD2 in DYNAMICS\n')
        return len(data)

    def _read_tstep(self, data, n):
        """TSTEP(8307,83,142) - Record 38"""
        if self.is_debug_file:
            self.binary_debug.write('skipping TSTEP in DYNAMICS\n')
        return len(data)

#UNBALNC
