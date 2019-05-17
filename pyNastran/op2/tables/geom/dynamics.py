"""defines readers for BDF objects in the OP2 DYNAMIC/DYNAMICS table"""
from __future__ import print_function
from struct import unpack, Struct
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
            #(3207, 33, 129): ['NONLIN3', self._read_fake], # 22
            (3307, 33, 129) : ['NONLIN3', self._read_fake],
            (3407, 34, 130): ['NONLIN4', self._read_fake], # 23
            (2107, 21, 195): ['RANDPS', self._read_randps], # 24
            (2207, 22, 196): ['RANDT1', self._read_randt1], # 25
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

            (10701, 107, 117) : ['RGYRO', self._read_rgyro],
            (10801, 108, 242) : ['ROTORG', self._read_fake],
            (3807, 38, 505) : ['NLRSFD', self._read_fake],
            (4807, 48, 306) : ['DYNRED', self._read_fake],
            (11001, 110, 310) : ['RSPINT', self._read_rspint],
            (10901, 109, 260) : ['RSPINR', self._read_fake],
            (11101, 111, 368) : ['UNBALNC', self._read_fake],

            # F:\work\pyNastran\examples\Dropbox\move_tpl\nlttlhxb.op2
            (7507, 75, 626) : ['TEMPD/TTEMP/TMPSET', self._read_fake],

            #F:\work\pyNastran\examples\Dropbox\move_tpl\rcross01.op2
            (3201, 24, 54) : ['RCROSS', self._read_fake],
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
        self.increase_card_count('DAREA', nentries)
        struc = Struct(self._endian + b'3if')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
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
        self.increase_card_count('DELAY', nentries)
        struc = Struct(self._endian + b'3if')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
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
        #ndata = len(data)
        #nfields = (ndata - n) // 4

        datan = data[n:]
        ints = np.frombuffer(datan, self.idtype).copy()
        floats = np.frombuffer(datan, self.fdtype).copy()
        istart = 0
        iminus1_delta = get_iend_from_ints(ints)
        nentries = 0
        for iend in iminus1_delta:
            #datai = data[n+istart*4 : n+iend*4]
            sid = ints[istart]
            global_scale = floats[istart + 1]
            #print('  sid=%s global_scale=%s' % (sid, global_scale))
            deltai = iend - istart - 2 # subtract 2 for sid, global scale
            assert deltai % 2 == 0, (self.show_data(data[n+istart*4 : n+iend*4], 'if'))
            out = [sid, global_scale]
            scales = []
            load_ids = []
            for iscale in range(deltai // 2):
                scale = floats[istart + 2 + 2*iscale]
                load_id = ints[istart + 3 + 2*iscale]
                scales.append(scale)
                load_ids.append(load_id)
                out.append(scale)
                out.append(load_id)
            if self.is_debug_file:
                self.binary_debug.write('  DLOAD=%s\n' % str(out))

            dload = self.add_dload(sid, global_scale, scales, load_ids)
            istart = iend + 2
            nentries += 1
        self.increase_card_count('DLOAD', nentries)
        return n

    def _read_dphase(self, data, n):
        """
        DPHASE(77,19,184) - Record 5

        1 SID I Load set identification number
        2 P   I Grid, scalar, or extra point identification number
        3 C   I Component number
        4 TH  RS Phase lead

        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        self.increase_card_count('DPHASE', nentries)
        struc = Struct(self._endian + b'3if')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
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
        self.increase_card_count('EIGB', nentries)
        struc = Struct(self._endian + b'i8s ff 3i i 8s 4i')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            #self.show_data(edata[44:])
            # int, 8s, 2f, 3i, i, 8s, 4i
            out = struc.unpack(edata)
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
        """
        EIGC(207,2,87) - Record 8

        Word Name    Type Description
        1 SID           I Load set identification number
        2 METHOD(2) CHAR4 Method of eigenvalue extraction
        4 NORM(2)   CHAR4 Method for normalizing eigenvectors
        6 G             I Grid or scalar point identification number
        7 C             I Component number
        8 E            RS Convergence criterion
        9 ND1           I Number of desired eigenvectors
        10 CONTFLG      I Continuation flag
        CONTFLG=0 With continuation
        11 AAJ         RS Location of A on real axis
        12 WAJ         RS Location of A on imaginary axis
        13 ABJ         RS Location of B on real axis
        14 WBJ         RS Location of B on imaginary axis
        15 LJ          RS Width of search region
        16 NEJ          I Number of estimated roots
        17 NDJ          I Number of desired eigenvectors
        Words 11 through 17 repeat until (-1,-1,-1,-1,-1,-1,-1) occ
        CONTFLG =-1 Without continuation
        End CONTFLG

        """
        if 0:
            ntotal = 60
            nentries = (len(data) - n) // ntotal
            self.increase_card_count('EIGB', nentries)
            struct1 = Struct(self._endian + b'i 4s 4s 2i f 2i')
            #struct2 = Struct(self._endian + b'5f2i')
            for i in range(nentries):
                edata = data[n:n+ntotal]
                #self.show_data(edata[44:])
                # int, 8s, 2f, 3i, i, 8s, 4i
                out = struct1.unpack(edata)
                sid, method, L1, L2, nep, ndp, ndn, dunno, norm, g, c, dunno_a, dunno_b = out
                if self.is_debug_file:
                    self.binary_debug.write('EIGC=%s\n' % str(out))
                #print('out = %s' % str(out))
                method = method.strip().decode('latin1')
                norm = norm.strip().decode('latin1')
                eigb = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, g, c)
                self._add_method_object(eigb)
                n += ntotal
            #return n

        #-------------------------------------------------
        ndata = len(data)
        nfields = (ndata - n) // 4
        datan = data[n:]
        ints = unpack(self._endian + b'%ii' % nfields, datan)
        floats = unpack(self._endian + b'%if' % nfields, datan)
        strings = unpack(self._endian + b'4s'* nfields, datan)

        i = 0
        nentries = 0
        while i < nfields:
            sid = ints[i+0]
            method = (strings[i+1] + strings[i+2]).strip().decode('latin1')
            norm = (strings[i+3] + strings[i+4]).strip().decode('latin1')
            grid = ints[i+5]
            component = ints[i+6]
            epsilon = floats[i+7]
            neigenvalues = ints[i+8]
            control_flag = ints[i+9]

            #print('sid=%s method=%r norm=%r grid=%s component=%s epsilon=%s, '
            #      'neigenvalues=%s ctlflag=%s' % (
            #          sid, method, norm, grid, component, epsilon, neigenvalues, control_flag))

            alphaAjs = []
            omegaAjs = []
            alphaBjs = []
            omegaBjs = []
            LJs = []
            NEJs = []
            NDJs = []

            # dummy
            mblkszs = []
            iblkszs = []
            ksteps = []
            NJIs = []
            if control_flag == -1:
                datai = [sid, method, norm, grid, component, epsilon, neigenvalues, -1]
            elif control_flag == 0:
                intsi = ints[i+10:i+17]
                assert len(intsi) == 7, 'len=%s intsi=%s' % (len(intsi), intsi)
                datai = [sid, method, norm, grid, component, epsilon, neigenvalues, 0]
                while intsi != (-1, -1, -1, -1, -1, -1, -1):
                    aaj, waj, abj, wbj, lj = floats[i+10:i+15]
                    #print('aaj=%s waj=%s abj=%s wbj=%s lj=%s'  % (
                        #aaj, waj, abj, wbj, lj))
                    nej, ndj = ints[i+15:i+17]
                    datai.extend([(aaj, waj, abj, wbj, lj, nej)])

                    alphaAjs.append(aaj)
                    omegaAjs.append(waj)
                    alphaBjs.append(abj)
                    omegaBjs.append(wbj)
                    if norm == 'MAX':
                        # faked
                        mblkszs.append(0.)
                        iblkszs.append(0)
                        ksteps.append(0)
                        NJIs.append(0)

                    LJs.append(lj)
                    NEJs.append(nej)
                    NDJs.append(ndj)
                    #print('aaj=%s waj=%s abj=%s wbj=%s lj=%s nej=%s ndj=%s'  % (
                        #aaj, waj, abj, wbj, lj, nej, ndj
                    #))
                    i += 7
                    intsi = ints[i+10:i+17]
                    #print('intsi = ', intsi)
                    assert len(intsi) == 7, 'len=%s intsi=%s' % (len(intsi), intsi)
                    #print("intsi = ", intsi)
                    #print()
                    #11 AAJ RS Location of A on real axis
                    #12 WAJ RS Location of A on imaginary axis
                    #13 ABJ RS Location of B on real axis
                    #14 WBJ RS Location of B on imaginary axis
                    #15 LJ RS Width of search region
                    #16 NEJ I Number of estimated roots
                    #17 NDJ I Number of desired eigenvectors
                assert len(intsi) == 7, intsi
                #print('intsi = ', intsi)
                #raise NotImplementedError('EIGC control_flag=%s' % control_flag)
            else:
                raise NotImplementedError('EIGC control_flag=%s' % control_flag)
            datai.extend([-1, -1, -1, -1, -1, -1, -1])  # creates a +7

            if self.is_debug_file:
                self.binary_debug.write('  EIGC=%s\n' % str(datai))

            if grid == 0:
                grid = None
                assert component == 0, component
                component = None

            eigc = self.add_eigc(sid, method, grid, component, epsilon, neigenvalues, norm=norm,
                                 mblkszs=mblkszs, iblkszs=iblkszs, ksteps=ksteps, NJIs=NJIs,
                                 alphaAjs=alphaAjs, omegaAjs=omegaAjs,
                                 alphaBjs=alphaBjs, omegaBjs=omegaBjs,
                                 LJs=LJs, NEJs=NEJs, NDJs=NDJs,
                                 shift_r1=None, shift_i1=None, isrr_flag=None,
                                 nd1=None, comment='')
            eigc.validate()

            #while intsi != (-1, -1, -1):
                #gridi, compi, coeffi = ints[i+4], ints[i+5], floats[i+6]
                #mpc_data.extend([gridi, compi, coeffi])
                #nodes.append(gridi)
                #components.append(compi)
                #coefficients.append(coeffi)
                #i += 3

            # +10 is for the prefix; +7 is for the -1s
            i += 10 + 7 # 3 + 4 from (-1,-1,-1) and (sid,grid,comp,coeff)
            nentries += 1
            #print('--------------')
        self.increase_card_count('EIGC', nentries)
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
        self.increase_card_count('EIGP', nentries)
        struct1 = Struct('i2fi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            sid, alpha, omega, m = out
            if self.is_debug_file:
                self.binary_debug.write('EIGP=%s\n' % str(out))
            #print('out = %s' % str(out))
            self.add_eigp(sid, alpha, omega, m, alpha2=None, omega2=None, m2=None)
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
        ntotal = 72
        nentries = (len(data) - n) // ntotal
        struct1 = Struct('i 8s 2f 4i 8s 7i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
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
        self.increase_card_count('EIGR', nentries)
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
        ndata = len(data)
        s = Struct('i 2f 3i f 2i 8s f i')
        nbytes = 52
        nentries = 0
        while n < ndata:
            #edata = data[n:n+46] # 11*4+2 = 46
            edata = data[n:n+nbytes] # 13*4 = 52
            out = s.unpack(edata)
            sid, v1, v2, nd, msglvl, maxset, shfscl, flag1, flag2, norm, alpha, nums = out
            norm = norm.decode('latin1').rstrip('\x00 ')
            if nums != 538976288:
                assert nums < 10000, nums
                #edata2 = data[n+46:n+46+nums*4]
                edata2 = data[n+nbytes:n+nbytes+nums*4]
                self.show_data(edata2, 'if')
                fi = unpack('%if' % nums, edata2)
                print(out, fi)
                raise NotImplementedError(nums)

            if self.is_debug_file:
                self.binary_debug.write('  EIGRL=%s\n' % str(out))
            options = []
            values = []
            eigrl = self.add_eigrl(sid, v1=v1, v2=v2, nd=nd, msglvl=msglvl,
                                   maxset=maxset, shfscl=shfscl,
                                   norm=norm, options=options, values=values)
            #print(eigrl)
            if nums == 538976288:
                n = len(data)
            else:
                n += nbytes + nums * 4
            nentries += 1
        self.increase_card_count('EIGRL', nentries)
        return n

    def _read_epoint(self, data, n):
        """EPOINT(707,7,124) - Record 12"""
        npoints = (len(data) - n) // 4
        fmt = self._endian + b'%ii' % npoints
        nids = unpack(fmt, data[n:])
        if self.is_debug_file:
            self.binary_debug.write('EPOINT=%s\n' % str(nids))
        epoint = EPOINTs.add_op2_data(list(nids))
        self._add_epoint_object(epoint)
        self.card_count['EPOINT'] = npoints
        self.increase_card_count('EPOINT', count_num=npoints)
        return n

    def _read_freq(self, data, n):
        """FREQ(1307,13,126) - Record 13"""
        ints = np.frombuffer(data[n:], self.idtype).copy()
        floats = np.frombuffer(data[n:], self.fdtype).copy()
        iminus1 = np.where(ints == -1)[0]
        istart = 0
        for iend in iminus1:
            sid = ints[istart]
            freqs = floats[istart + 1:iend]
            istart = iend + 1
            self.add_freq(sid, freqs)
        self.increase_card_count('FREQ', count_num=len(iminus1))
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
        struc = Struct(self._endian + b'iffi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, df, ndf = out
            if self.is_debug_file:
                self.binary_debug.write('FREQ1=%s\n' % str(out))
            #print('out = %s' % str(out))
            freq = FREQ1(sid, f1, df, ndf=ndf)
            self._add_freq_object(freq)
            n += ntotal
        self.increase_card_count('FREQ1', nentries)
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
        struc = Struct(self._endian + b'iffi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, nf = out
            if self.is_debug_file:
                self.binary_debug.write('  FREQ2=%s\n' % str(out))
            #print('out = %s' % str(out))
            self.add_freq2(sid, f1, f2, nf)
            n += ntotal
        self.increase_card_count('FREQ2', nentries)
        return n

    def _read_freq3(self, data, n):
        """
        FREQ3(1407,14,39) - Record 16

        1 SID      I Set identification number
        2 F1      RS Lower bound of modal frequency range
        3 F2      RS Upper bound of modal frequency range
        4 TYPE CHAR4 Type of interpolation: LINE or LOG
        5 NEF      I Number of frequencies
        6 BIAS    RS Clustering bias parameter

        """
        ntotal = 24 # 4*6
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'i 2f 4s if')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, freq_type, nef, bias = out
            freq_type = freq_type.strip().decode('latin1')
            if freq_type == 'LINE':
                freq_type = 'LINEAR'

            if self.is_debug_file:
                self.binary_debug.write('  FREQ3=%s\n' % str(out))

            self.add_freq3(sid, f1, f2=f2, Type=freq_type, nef=nef, cluster=bias)
            n += ntotal
        self.increase_card_count('FREQ3', nentries)
        return n

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
        struc = Struct(self._endian + b'i 3f i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, fspread, nfm = out
            if self.is_debug_file:
                self.binary_debug.write('  FREQ4=%s\n' % str(out))
            self.add_freq4(sid, f1, f2, fspread=fspread, nfm=nfm)
            n += ntotal
        self.increase_card_count('FREQ4', nentries)
        return n

    def _read_freq5(self, data, n):
        """
        FREQ5(1607,16,41) - Record 18

        1 SID  I Load set identification number
        2 F1  RS Lower bound of modal frequency range
        3 F2  RS Upper bound of modal frequency range
        4 FRI RS Fractions of natural frequencies

        """
        ints = np.frombuffer(data, dtype='int32').copy()
        floats = np.frombuffer(data, dtype='float32').copy()
        i_minus_1s = np.where(ints == -1)[0]
        nentries = len(i_minus_1s)

        i0 = 0
        for i_minus_1 in i_minus_1s:
            sid = ints[i0]
            floatsi = floats[i0 + 1:i_minus_1]
            if self.is_debug_file:
                self.binary_debug.write('  FREQ5=(%s, %s)\n' % (sid, floatsi.tolist()))
            f1 = floatsi[0]
            f2 = floatsi[1]
            fractions = floatsi[2:]
            self.add_freq5(sid, fractions, f1=f1, f2=f2)
            #print('freq =', freq)
            i0 = i_minus_1 + 1
        self.increase_card_count('FREQ5', nentries)
        return len(data)

        #ntotal = 20 # 4*5
        #nentries = (len(data) - n) // ntotal
        #struc = Struct(self._endian + b'i 3f')
        #for i in range(nentries):
            #edata = data[n:n+ntotal]
            #out = struc.unpack(edata)
            #sid, f1, f2, fspread, nfm = out
            #if self.is_debug_file:
                #self.binary_debug.write('  FREQ5=%s\n' % str(out))
            #freq = self.add_freq4(sid, f1, f2, fspread=fspread, nfm=nfm)
            #n += ntotal
        #self.increase_card_count('FREQ5', nentries)


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
        struc = Struct(self._endian + b'3i 2f if')
        for unused_i in range(nentries):
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
        self.increase_card_count('RANDPS', nentries)
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
        6 TIDI I Identification number of a TABRNDi entry that defines G(f)
        7 TIDR RS If TIDI = 0, constant value for G(f)

        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'3i2fi')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, j, k, x, y, tid = out
            if self.is_debug_file:
                self.binary_debug.write('  RANDPS=%s\n' % str(out))
            #self.log.debug('  RANDPS=%s\n' % str(out))
            self.add_randps(sid, j, k, x=x, y=y, tid=tid)
            n += ntotal
        self.increase_card_count('RANDPS', nentries)
        return n, []

    def _read_randt1(self, data, n):
        """
        RANDT1(2207,22,196)

        Word Name Type Description
        1 SID  I  Set identification number
        2 N    I  Number of time lag intervals
        3 TO   RS Starting time lag
        4 TMAX RS Maximum time lag

        """
        ntotal = 16  # 4*4
        struct1 = Struct(self._endian + b'2i 2f')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            if self.is_debug_file:
                self.binary_debug.write('  RANDT1=%s\n' % str(out))
            sid, nlags, to, tmax = out
            self.add_randt1(sid, nlags, to, tmax)
            n += ntotal
        self.card_count['RANDT1'] = nentries
        return n

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
        struc = Struct(self._endian + b'7i 4f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, dphasei, tci, tdi, load_type, delayr, dphaser, tcr, tdr = out # 44

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
                           Type=load_type, comment='')
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
        struc = Struct(self._endian + b'2i 2f 3i 2f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, dphaser, delayr, tci, tdi, load_type, tau, phi = out # 36

            if self.is_debug_file:
                self.binary_debug.write('  RLOAD1=%s\n' % str(out))

            tc = tci
            td = tdi
            delay = delayr
            dphase = dphaser
            dload = RLOAD1(sid, darea, delay=delay, dphase=dphase, tc=tc, td=td,
                           Type=load_type, comment='')
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
        struc = Struct(self._endian + b'7i 4f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, dphasei, tbi, tpi, load_type, delayr, dphaser, tbr, tpr = out
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
                           Type=load_type, comment='')
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
        struc = Struct(self._endian + b'7i 2f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, dphasei, delayi, tbi, tpi, load_type, tau, phase = out
            if self.is_debug_file:
                self.binary_debug.write('  RLOAD2=%s\n' % str(out))

            dload = RLOAD2(sid, darea, delay=delayi, dphase=dphasei, tb=tbi, tp=tpi,
                           Type=load_type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_rgyro(self, data, n):
        """
        FREQ4(1507,15,40) - Record 17

        1 SID          I RGYRO identification number
        2 ATYPE(2) CHAR4 ASYNC/SYNC flag
        4 REFROT       I Reference rotor identification number
        5 UNIT(2)  CHAR4 RPM/FREQ flag for speed input
        7 SPDLOW      RS Lower limit of speed range
        8 SPDHGH      RS Upper limit of speed range
        9 SPEED       RS Specific speed

        """
        ntotal = 36 # 4*9
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'i 8s i 8s 3f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, asynci, refrot, unit, speed_low, speed_high, speed = out
            asynci = asynci.decode('latin1')
            unit = unit.decode('latin1')
            if self.is_debug_file:
                self.binary_debug.write('  RGYRO=%s\n' % str(out))
            self.add_rgyro(sid, asynci, refrot, unit, speed_low, speed_high, speed)
            n += ntotal
        self.increase_card_count('RGYRO', nentries)
        return n

#ROTORG
#RSPINR

    def _read_rspint(self, data, n):
        """
        RSPINT(11001,110,310) - Record 31

        1 RID         I Rotor identification number
        2 GRIDA       I Grid A for rotation direction vector
        3 GRIDB       I Grid B for rotation direction vector
        4 GR         RS Rotor damping coefficient
        5 UNIT(2) CHAR4 RPM/FREQ flag for speed input
        7 TABLEID     I Table identification number for speed history

        """
        #self.show_data(data[12:], 'ifs')
        ntotal = 28 # 4*7
        nentries = (len(data) - n) // ntotal
        struc = Struct(self._endian + b'3if 8s i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            rid, grida, gridb, gr, unit, table_id = out
            unit = unit.decode('latin1')
            if self.is_debug_file:
                self.binary_debug.write('  RSPINT=%s\n' % str(out))
            self.add_rspint(rid, grida, gridb, gr, unit, table_id)
            n += ntotal
        self.increase_card_count('RSPINT', nentries)
        return n

#SEQEP(5707,57,135)

    def _read_tf(self, data, n):
        """TF"""
        # subtract of the header (sid, nid, component, b0, b1, b2)
        # divide by 5 (nid1, component1, a0, a1, a2)
        #nrows = (nfields - 6) // 5
        #print('n=%s nrows=%s' % (n, nrows))
        #print(self.show_data(data))
        #nid1, component1, a0, a1, a2

        ndata = len(data)
        struct1 = Struct(self._endian + b'3i3f')
        struct2 = Struct(self._endian + b'2i3f')
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
            self.increase_card_count('TF')
            n = n3
        return n

    def _read_tic(self, data, n):
        """
        TIC(6607,66,137)

        1 SID I Load set identification number
        2 G   I Grid, scalar, or extra point identification number
        3 C   I Component number for point GD
        4 U0 RS Initial displacement
        5 V0 RS Initial velocity

        """
        ntotal = 20  # 5*4
        struct1 = Struct(self._endian + b'3i 2f')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            if self.is_debug_file:
                self.binary_debug.write('  TIC=%s\n' % str(out))
            sid, nid, comp, u0, v0 = out
            self.add_tic(sid, [nid], [comp], u0=u0, v0=v0)
            n += ntotal
        self.card_count['TIC'] = nentries
        return n

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
        struc = Struct(self._endian + b'5i 3f')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, tid, delayr, us0, vs0 = out
            if self.is_debug_file:
                self.binary_debug.write('TLOAD1=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD1(sid, darea, tid, delay=delay, Type=load_type,
                           us0=us0, vs0=vs0)
            self._add_dload_entry(dload)
            n += ntotal
        self.increase_card_count('TLOAD1', nentries)
        return n

    def _read_tload2(self, data, n):
        """TLOAD2"""
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
        struc = Struct(self._endian + b'4i 7f 2f')
        for i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, t1, t2, freq, p, c, growth, delayr, us0, vs0 = out
            if self.is_debug_file:
                self.binary_debug.write('  TLOAD2=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD2(sid, darea, delay=delay, Type=load_type, T1=t1,
                           T2=t2, frequency=freq,
                           phase=p, c=c, b=growth,
                           us0=us0, vs0=vs0)
            self._add_dload_entry(dload)
            n += ntotal
        self.increase_card_count('TLOAD2', nentries)
        return n

    def _read_tstep(self, data, n):
        """TSTEP(8307,83,142) - Record 38"""
        self.log.info('skipping TSTEP in DYNAMICS\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping TSTEP in DYNAMICS\n')
        return len(data)

#UNBALNC
