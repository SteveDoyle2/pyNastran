"""
defines readers for BDF objects in the OP2 DYNAMIC/DYNAMICS table
"""
from __future__ import print_function
from struct import unpack, Struct
from six import b
import numpy as np

from pyNastran.bdf.cards.nodes import EPOINTs
from pyNastran.bdf.cards.loads.loads import DAREA
from pyNastran.bdf.cards.methods import EIGB
from pyNastran.bdf.cards.dynamic import FREQ1, TF, DELAY, DPHASE
from pyNastran.bdf.cards.loads.dloads import TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.tables.geom.dit import get_iend_from_ints

class DYNAMICS(GeomCommon):
    """defines methods for reading op2 dynamics loads/methods"""
    def _read_dynamics_4(self, data, ndata):
        return self._read_geom_4(self._dynamics_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._dynamics_map = {
            (5307, 53, 379) : ['ACSRCE', self._read_acsrce], # 1
            (27, 17, 182): ['DAREA', self._read_darea],  # 2

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
            (2107, 21, 195): ['RANDPS', self._read_randps], # 24
            (2207, 22, 196): ['RANDT1', self._read_fake], # 25
            (5107, 51, 131): ['RLOAD1', self._read_rload1],  # 26
            (5207, 52, 132): ['RLOAD2', self._read_rload2],  # 27
            (8910, 89, 606): ['ROTORB', self._read_fake],  # 28
            (8210, 82, 599): ['ROTORD', self._read_fake],  # 29
            (8410, 84, 600): ['ROTORG', self._read_fake],  # 30
            (5707, 57, 135): ['SEQEP', self._read_fake],  # 31
            (6207, 62, 136): ['TF', self._read_tf],  # 32
            (6607, 66, 137): ['TIC', self._read_tic],  # 33
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

    def _read_acsrce(self, data, n):
        """ACSRCE(5307,53,379)"""
        self.log.info('skipping ACSRCE in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping ACSRCE in DYNAMICS\n')
        return len(data)

    def _read_darea(self, data, n):
        """
        DAREA(27,17,182) - the marker for Record 2

        1 SID I  Load set identification number
        2 P   I  Grid, scalar, or extra point identification number
        3 C   I  Component number
        4 A   RS Scale factor
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DAREA', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iiif', edata)
            #(sid,p,c,a) = out
            darea = DAREA.add_op2_data(data=out)
            self._add_darea_object(darea)
            n += ntotal
        return n

    def _read_delay(self, data, n):
        """
        DELAY(37,18,183) - Record 3

        1 SID I  Set identification number
        2 P   I  Grid, scalar, or extra point identification number
        3 C   I  Component number
        4 T   RS Time delay
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DELAY', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iiif', edata)
            sid, nodes, components, delays = out
            if self.is_debug_file:
                self.binary_debug.write('  DELAY=%s\n' % str(out))
            delay = DELAY(sid, nodes, components, delays)
            self._add_delay_object(delay)
            n += ntotal
        return n

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
        istart = 0
        iminus1_delta = get_iend_from_ints(ints)
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
            #print(dload)
            istart = iend + 2
            nentries += 1
        self._increase_card_count('DLOAD', nentries)
        return n

    def _read_dphase(self, data, n):
        """
        DPHASE(77,19,184) - Record 5

        1 SID I Load set identification number
        2 P   I Grid, scalar, or extra point identification number
        3 C   I Component number
        4 TH  RS Phase lead
        """
        #self.log.info('skipping DPHASE in DYNAMICS\n')
        #if self.is_debug_file:
            #self.binary_debug.write('skipping DPHASE in DYNAMICS\n')
        #return len(data)
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('DPHASE', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iiif', edata)
            sid, nodes, components, delays = out
            if self.is_debug_file:
                self.binary_debug.write('  DPHASE=%s\n' % str(out))
            delay = DPHASE(sid, nodes, components, delays)
            self._add_dphase_object(delay)
            n += ntotal
        return n

#DYNRED(4807,48,306)

    def _read_eigb(self, data, n):
        """EIGB(107,1,86) - Record 7"""
        ntotal = 60
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('EIGB', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            #self.show_data(edata[44:])
            # int, 8s, 2f, 3i, i, 8s, 4i
            out = unpack('i8s ff 3i i 8s 4i', edata)
            sid, method, L1, L2, nep, ndp, ndn, dunno, norm, g, c, dunno_a, dunno_b = out
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
        self.log.info('skipping EIGC in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGC in DYNAMICS\n')
        return len(data)

    def _read_eigp(self, data, n):
        """
        EIGP(257,4,158) - Record 9

        1 SID    I Load set identification number
        2 ALPHA RS Location of pole on real axis
        3 OMEGA RS Location of pole on imaginary axis
        4 M      I Multiplicity of complex root at pole
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self._increase_card_count('EIGP', nentries)
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('i2fi', edata)
            sid, alpha, omega, m = out
            if self.is_debug_file:
                self.binary_debug.write('EIGP=%s\n' % str(out))
            #print('out = %s' % str(out))
            eigp = self.add_eigp(sid, alpha, omega, m, alpha2=None, omega2=None, m2=None)
            n += ntotal
        return n

    def _read_eigr(self, data, n):
        """
        EIGR(307,3,85) - Record 10

        1 SID I Set identification number
        2 METHOD(2) CHAR4 Method of eigenvalue extraction
        4 F1 RS Lower bound of frequency range of interest
        5 F2 RS Upper bound of frequency range of interest
        6 NE  I Number of estimated roots
        7 ND  I Number of desired roots
        8 UNDEF(2 ) None
        10 NORM(2) CHAR4 Method for normalizing eigenvectors
        12 G  I Grid or scalar point identification number
        13 C  I Component number
        14 UNDEF(5 ) None
        """
        #self.log.info('skipping EIGR in DYNAMICS\n')
        #if self.is_debug_file:
            #self.binary_debug.write('skipping EIGR in DYNAMICS\n')
        #return len(data)
        ntotal = 72
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('i 8s 2f 4i 8s 7i', edata)
            (sid, method, f1, f2, ne, nd, null_a, null_b, norm, g, c,
             null_c, null_d, null_e, null_f, null_g) = out
            if self.is_debug_file:
                self.binary_debug.write('EIGR=%s\n' % str(out))
            method = method.strip().decode('latin1')
            norm = norm.strip().decode('latin1')
            eigr = self.add_eigr(sid, method=method, f1=f1, f2=f2, ne=ne, nd=nd,
                                 norm=norm, G=g, C=c, comment='')
            eigr.validate()
            n += ntotal
        self._increase_card_count('EIGR', nentries)
        return n

    def _read_eigrl(self, data, n):
        """
        EIGRL(308,8,348) - Record 11

        1 SID     I Set identification number
        2 V1     RS Lower bound of frequency range of interest
        3 V2     RS Upper bound of frequency range of interest
        4 ND      I Number of desired eigenvectors
        5 MSGLVL  I Diagnostic level
        6 MAXSET  I Number of vectors in block or set
        7 SHFSCL RS Estimate of first flexible mode
        8 FLAG1   I V1 specification flag - set to 1 if V1 is specified
        9 FLAG2   I V2 specification flag - set to 1 if V2 is specified
        10 NORM(2) CHAR4 Method for normalizing eigenvectors
        12 ALPH  RS Constant for quadratic frequency segment distribution
        13 NUMS   I Number of frequency segments
        14 FI    RS Frequency at the upper boundary of the i-th segment
        Word 14 repeats NUMS times
        """
        self.log.info('skipping EIGRL in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping EIGRL in DYNAMICS\n')
        return len(data)

        #self.show_data(data[n+36:n+100], 'ifs')
        #print(len(data[n:]) / 4.)
        #ndata = len(data)
        #while n < ndata:
            #edata = data[n:n+52] # 13*52
            #out = unpack('i 2f 3i f 2i 8s f i', edata)
            #sid, v1, v2, nd, msglvl, maxset, shfscl, flag1, flag2, norm, alpha, nums = out
            #norm = norm.strip().decode('latin1')
            #print("norm = %r" % norm)
            #print("nums = ", nums)
            #assert nums < 10000, nums
            #edata2 = data[n+52:n+52+nums*4]
            #self.show_data(edata2, 'if')
            #fi = unpack('%if' % nums, edata2)
            #print(out, fi)

            #if self.is_debug_file:
                #self.binary_debug.write('  EGIRL=%s\n' % str(out))
            #eigrl = self.add_eigrl(sid, v1=None, v2=None, nd=None, msglvl=0,
                                   #maxset=None, shfscl=None,
                                   #norm=None, options=None, values=None)
            #n += ntotal
        #self._increase_card_count('EGIRL', nentries)
        #return n

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
        ints = np.fromstring(data[n:], self.idtype)
        floats = np.fromstring(data[n:], self.fdtype)
        iminus1 = np.where(ints == -1)[0]
        istart = 0
        for iend in iminus1:
            sid = ints[istart]
            freqs = floats[istart + 1:iend]
            istart = iend + 1
            self.add_freq(sid, freqs)
        self._increase_card_count('FREQ', count_num=len(iminus1))
        return n

    def _read_freq1(self, data, n):
        """
        FREQ1(1007,10,125) - Record 14

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
        """
        FREQ2(1107,11,166) - Record 15

        1 SID  I Set identification number
        2 F1  RS First frequency
        3 F2  RS Last frequency
        4 NF   I Number of logarithmic intervals
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('iffi', edata)
            sid, f1, f2, nf = out
            if self.is_debug_file:
                self.binary_debug.write('  FREQ2=%s\n' % str(out))
            #print('out = %s' % str(out))
            freq = self.add_freq2(sid, f1, f2, nf)
            n += ntotal
        self._increase_card_count('FREQ2', nentries)
        return n

    def _read_freq3(self, data, n):
        """FREQ3(1407,14,39) - Record 16"""
        self.log.info('skipping FREQ3 in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ3 in DYNAMICS\n')
        return len(data)

    def _read_freq4(self, data, n):
        """
        FREQ4(1507,15,40) - Record 17

        1 SID   I Set identification number
        2 F1   RS Lower bound of modal frequency range
        3 F2   RS Upper bound of modal frequency range
        4 FSPD RS Frequency spread
        5 NFM   I Number of evenly spaced frequencies per spread
        """
        ntotal = 20 # 4*5
        nentries = (len(data) - n) // ntotal
        struc = Struct(b('i 3f i'))
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, fspread, nfm = out
            if self.is_debug_file:
                self.binary_debug.write('  FREQ4=%s\n' % str(out))
            freq = self.add_freq4(sid, f1, f2, fspread=fspread, nfm=nfm)
            n += ntotal
        self._increase_card_count('FREQ4', nentries)
        return n

    def _read_freq5(self, data, n):
        """FREQ5(1607,16,41) - Record 18"""
        self.log.info('skipping FREQ5 in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping FREQ5 in DYNAMICS\n')
        return len(data)

#NLRSFD
#NOLIN1
#NOLIN2
#NOLIN3
#NOLIN4

    def _read_randps(self, data, n):
        """common method for reading NX/MSC RLOAD1"""
        n = self._read_dual_card(data, n, self._read_randps_nx, self._read_randps_msc,
                                 'RLOAD1', self._add_dload_entry)
        return n

    def _read_randps_nx(self, data, n):
        """
        RANDPS(2107,21,195)

        NX
        1 SID   I Set identification number
        2 J     I Subcase identification number of the excited set
        3 K     I Subcase identification number of the applied load set
        4 X    RS X component
        5 Y    RS Y component
        6 TIDI  I Identification number of a TABRNDi entry that defines G(f)
        7 TIDR RS If TIDI = 0, constant value for G(f)
        """
        ntotal = 28
        nentries = (len(data) - n) // ntotal
        struc = Struct(b('3i 2f if'))
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, j, k, x, y, tidi, tidf = out
            tid = tidi
            if tidi == 0:
                tid = tidf
            if self.is_debug_file:
                self.binary_debug.write('  RANDPS=%s\n' % str(out))
            #self.log.debug('  RANDPS=%s\n' % str(out))
            self.add_randps(sid, j, k, x=x, y=y, tid=tid)
            n += ntotal
        self._increase_card_count('RANDPS', nentries)
        return n, []

    def _read_randps_msc(self, data, n):
        """
        RANDPS(2107,21,195)

        MSC
        1 SID I  Set identification number
        2 J   I  Subcase identification number of the excited set
        3 K   I  Subcase identification number of the applied load set
        4 X   RS X component
        5 Y   RS Y component
        6 TID I  Identification number of a TABRNDi entry that defines G(F)

        NX
        1 SID I Set identification number
        2 J I Subcase identification number of the excited set
        3 K I Subcase identification number of the applied load set
        4 X RS X component
        5 Y RS Y component
        6 TIDI I Identification number of a TABRNDi entry that defines
        G(f)
        7 TIDR RS If TIDI = 0, constant value for G(f)
        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        struc = Struct(b('3i2fi'))
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, j, k, x, y, tid = out
            if self.is_debug_file:
                self.binary_debug.write('  RANDPS=%s\n' % str(out))
            #self.log.debug('  RANDPS=%s\n' % str(out))
            self.add_randps(sid, j, k, x=x, y=y, tid=tid)
            n += ntotal
        self._increase_card_count('RANDPS', nentries)
        return n, []

#RANDT1

    def _read_rload1(self, data, n):
        """common method for reading NX/MSC RLOAD1"""
        n = self._read_dual_card(data, n, self._read_rload1_nx, self._read_rload1_msc,
                                 'RLOAD1', self._add_dload_entry)
        return n

    def _read_rload1_nx(self, data, n):
        """
        RLOAD1(5107,51,131) - Record 26

        NX
        1 SID      I Load set identification number
        2 DAREA    I DAREA Bulk Data entry identification number
        3 DELAYI   I DELAY Bulk Data entry identification number
        4 DPHASEI  I DPHASE Bulk Data entry identification number
        5 TCI      I TABLEDi Bulk Data entry identification number for C(f)
        6 TDI      I TABLEDi Bulk Data entry identification number for D(f)
        7 TYPE     I Nature of the dynamic excitation
        8 DELAYR  RS If DELAYI = 0, constant value for delay
        9 DPHASER RS If DPHASEI = 0, constant value for phase
        10 TCR    RS If TCI = 0, constant value for C(f)
        11 TDR    RS If TDI = 0, constant value for D(f)
        """
        dloads = []
        ntotal = 44
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('7i 4f', edata)
            sid, darea, delayi, dphasei, tci, tdi, Type, delayr, dphaser, tcr, tdr = out # 44

            if self.is_debug_file:
                self.binary_debug.write('  RLOAD1=%s\n' % str(out))

            tc = tci
            td = tdi
            delay = delayi
            dphase = dphasei
            if delayi == 0:
                delay = delayr
            if dphasei == 0:
                dphase = dphaser
            if tci == 0:
                tc = tcr
            if tdi == 0:
                td = tdr
            dload = RLOAD1(sid, darea, delay=delay, dphase=dphase, tc=tc, td=td,
                           Type=Type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_rload1_msc(self, data, n):
        """
        RLOAD1(5107,51,131) - Record 26

        MSC
        1 SID    I  Load set identification number
        2 DAREA  I  DAREA Bulk Data entry identification number
        3 DPHASE RS DPHASE Bulk Data entry identification number
        4 DELAY  RS DELAY Bulk Data entry identification number
        5 TC     I  TABLEDi Bulk Data entry identification number for C(f)
        6 TD     I  TABLEDi Bulk Data entry identification number for D(f)
        7 TYPE   I  Nature of the dynamic excitation
        8 T      RS Time delay
        9 PH     RS Phase lead
        """
        dloads = []
        ntotal = 36
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('2i 2f 3i 2f', edata)
            sid, darea, dphaser, delayr, tci, tdi, Type, tau, phi = out # 36

            if self.is_debug_file:
                self.binary_debug.write('  RLOAD1=%s\n' % str(out))

            tc = tci
            td = tdi
            delay = delayr
            dphase = dphaser
            dload = RLOAD1(sid, darea, delay=delay, dphase=dphase, tc=tc, td=td,
                           Type=Type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_rload2(self, data, n):
        """common method for reading NX/MSC RLOAD2"""
        n = self._read_dual_card(data, n, self._read_rload2_nx, self._read_rload2_msc,
                                 'RLOAD2', self._add_dload_entry)
        return n

    def _read_rload2_nx(self, data, n):
        """
        RLOAD2(5207,52,132) - Record 27
        NX
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
        dloads = []
        ntotal = 44
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('7i 4f', edata)
            sid, darea, delayi, dphasei, tbi, tpi, Type, delayr, dphaser, tbr, tpr = out
            if self.is_debug_file:
                self.binary_debug.write('  RLOAD2=%s\n' % str(out))

            tb = tbi
            tp = tpi
            delay = delayi
            dphase = dphasei
            if tbi == 0:
                tb = tbr
            if tpi == 0:
                tp = tpr
            if delayi == 0:
                delay = delayr
            if dphasei == 0:
                dphase = dphaser
            dload = RLOAD2(sid, darea, delay=delay, dphase=dphase, tb=tb, tp=tp,
                           Type=Type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_rload2_msc(self, data, n):
        """
        RLOAD2(5207,52,132) - Record 27
        MSC
        1 SID    I  Load set identification number
        2 DAREA  I  DAREA Bulk Data entry identification number
        3 DPHASE I  DPHASE Bulk Data entry identification number
        4 DELAY  I  DELAY Bulk Data entry identification number
        5 TB     I  TABLEDi Bulk Data entry identification number for B(f)
        6 TP     I  TABLEDi Bulk Data entry identification number for Phi(f)
        7 TYPE   I  Nature of the dynamic excitation
        8 T      RS Time delay
        9 PH     RS Phase lead
        """
        dloads = []
        ntotal = 36
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('7i 2f', edata)
            sid, darea, dphasei, delayi, tbi, tpi, Type, tau, phase = out
            if self.is_debug_file:
                self.binary_debug.write('  RLOAD2=%s\n' % str(out))

            dload = RLOAD2(sid, darea, delay=delayi, dphase=dphasei, tb=tbi, tp=tpi,
                           Type=Type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

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

    def _read_tic(self, data, n):
        """TIC"""
        self.log.info('skipping TIC in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TIC in DYNAMICS\n')
        return len(data)

#TIC3

    def _read_tload1(self, data, n):
        """
        TLOAD1(7107,71,138) - Record 37

        1 SID    I  Load set identification number
        2 DAREA  I  DAREA Bulk Data entry identification number
        3 DELAYI I  DELAY Bulk Data entry identification number
        4 TYPE   I  Nature of the dynamic excitation
        5 TID    I  Identification number of TABLEDi entry that gives F(t)
        6 DELAYR RS If DELAYI = 0, constant value for delay
        6 U0     RS Initial displacement factor for enforced motion (MSC; NX undocumented)
        7 V0     RS Initial velocity factor for enforced motion (MSC; NX undocumented)
        8 T      RS Time delay (MSC)
        """
        ntotal = 8*4
        #self.show_data(data[n:], 'if')
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('5i 3f', edata)
            sid, darea, delayi, Type, tid, delayr, us0, vs0 = out
            if self.is_debug_file:
                self.binary_debug.write('TLOAD1=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD1(sid, darea, tid, delay=delay, Type=Type,
                           us0=us0, vs0=vs0)
            self._add_dload_entry(dload)
            n += ntotal
        self._increase_card_count('TLOAD1', nentries)
        return n

    def _read_tload2(self, data, n):
        return self._read_tload2_nx(data, n)

    def _read_tload2_nx(self, data, n):
        """
        TLOAD2(7207,72,139) - Record 37

        NX
        1 SID      I  Load set identification number
        2 DAREA    I DAREA Bulk Data entry identification number
        3 DELAYI   I DELAY Bulk Data entry identification number
        4 TYPE     I Nature of the dynamic excitation
        5 T1      RS Time constant 1
        6 T2      RS Time constant 2
        7 F       RS Frequency
        8 P       RS Phase angle
        9 C       RS Exponential coefficient
        10 B      RS Growth coefficient
        11 DELAYR RS If DELAYI = 0, constant value for delay
        12 US0    RS not documented in NX
        13 VS0    RS not documented in NX
        """
        ntotal = 52
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = unpack('4i 7f 2f', edata)
            sid, darea, delayi, Type, t1, t2, freq, p, c, growth, delayr, us0, vs0 = out
            if self.is_debug_file:
                self.binary_debug.write('  TLOAD2=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD2(sid, darea, delay=delay, Type=Type, T1=t1,
                           T2=t2, frequency=freq,
                           phase=p, c=c, b=growth,
                           us0=us0, vs0=vs0)
            self._add_dload_entry(dload)
            n += ntotal
        self._increase_card_count('TLOAD2', nentries)
        return n

    def _read_tstep(self, data, n):
        """TSTEP(8307,83,142) - Record 38"""
        self.log.info('skipping TSTEP in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TSTEP in DYNAMICS\n')
        return len(data)

#UNBALNC
