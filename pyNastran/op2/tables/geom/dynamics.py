"""defines readers for BDF objects in the OP2 DYNAMIC/DYNAMICS table"""
from __future__ import annotations
from struct import unpack, Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.nodes import EPOINTs
from pyNastran.bdf.cards.loads.loads import DAREA
from pyNastran.bdf.cards.methods import EIGB
from pyNastran.bdf.cards.dynamic import FREQ1, TF, DELAY, DPHASE
from pyNastran.bdf.cards.loads.dloads import TLOAD1, TLOAD2, RLOAD1, RLOAD2, ACSRCE
from pyNastran.op2.op2_interface.op2_reader import mapfmt # , reshape_bytes_block
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.tables.geom.dit import get_iend_from_ints

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom

class DYNAMICS(GeomCommon):
    """defines methods for reading op2 dynamics loads/methods"""
    def read_dynamics_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.dynamics_map, data, ndata)

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def read_fake(self, data: bytes, n: int) -> int:
        return self.op2._read_fake(data, n)

    def read_stop(self, data: bytes, n: int) -> int:
        return self.op2.reader_geom1.read_stop(data, n)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        self.dynamics_map = {
            (5307, 53, 379) : ['ACSRCE', self.read_acsrce], # 1
            (27, 17, 182): ['DAREA', self.read_darea],  # 2

            (37, 18, 183): ['DELAY', self.read_delay],  # 3
            (57, 5, 123): ['DLOAD', self.read_dload],  # 4
            (77, 19, 184): ['DPHASE', self.read_dphase],  # 5
            (107, 1, 86): ['EIGB', self.read_eigb],   # 7
            (207, 2, 87): ['EIGC', self.read_eigc],   # 8
            (257, 4, 158): ['EIGP', self.read_eigp],   # 9
            (307, 3, 85): ['EIGR', self.read_eigr],   # 10

            (308, 8, 348): ['EIGRL', self.read_eigrl],  # 11
            (707, 7, 124): ['EPOINT', self.read_epoint],  # 12
            (1307, 13, 126): ['FREQ', self.read_freq],   # 13
            (1007, 10, 125): ['FREQ1', self.read_freq1],  # 14
            (1107, 11, 166): ['FREQ2', self.read_freq2],  # 15
            (1407, 14, 39): ['FREQ3', self.read_freq3],  # 16
            (1507, 15, 40): ['FREQ4', self.read_freq4],  # 17
            (1607, 16, 41): ['FREQ5', self.read_freq5],  # 18

            (3707, 37, 556) : ['NLRGAP', self.read_nrlgap], # 19

            (3107, 31, 127): ['NOLIN1', self.read_nolin1], # 20
            (3207, 32, 128): ['NOLIN2', self.read_nolin2], # 21
            #(3207, 33, 129): ['NOLIN3', self.read_fake], # 22
            (3307, 33, 129) : ['NOLIN3', self.read_nolin3],
            (3407, 34, 130): ['NOLIN4', self.read_nolin4], # 23
            (2107, 21, 195): ['RANDPS', self.read_randps], # 24
            (2207, 22, 196): ['RANDT1', self.read_randt1], # 25
            (5107, 51, 131): ['RLOAD1', self.read_rload1],  # 26
            (5207, 52, 132): ['RLOAD2', self.read_rload2],  # 27
            (8910, 89, 606): ['ROTORB', self.read_fake],  # 28
            (8210, 82, 599): ['ROTORD', self.read_rotord],  # 29
            (8410, 84, 600): ['ROTORG', self.read_rotorg],  # 30
            (5707, 57, 135): ['SEQEP', self.read_fake],  # 31
            (6207, 62, 136): ['TF', self.read_tf],  # 32
            (6607, 66, 137): ['TIC', self.read_tic],  # 33
            (7107, 71, 138): ['TLOAD1', self.read_tload1],  # 37
            (7207, 72, 139): ['TLOAD2', self.read_tload2],  # 38
            (8307, 83, 142): ['TSTEP', self.read_tstep],  # 39
            (17500, 175, 618): ['TSTEP', self.read_tstep1],

            (10701, 107, 117) : ['RGYRO', self.read_rgyro],
            (10801, 108, 242) : ['ROTORG', self.read_fake],
            (3807, 38, 505) : ['NLRSFD', self.read_nlrsfd],
            (4807, 48, 306) : ['DYNRED', self.read_dynred],
            (11001, 110, 310) : ['RSPINT', self.read_rspint],
            (10901, 109, 260) : ['RSPINR', self.read_fake],
            (11101, 111, 368) : ['UNBALNC', self.read_fake],

            # F:\work\pyNastran\examples\Dropbox\move_tpl\nlttlhxb.op2
            (7507, 75, 626) : ['TEMPD/TTEMP/TMPSET', self.read_fake],

            #F:\work\pyNastran\examples\Dropbox\move_tpl\rcross01.op2
            (3201, 24, 54) : ['RCROSS', self.read_rcross],

            (9010, 90, 569): ['CBEAR', self.read_fake],

            (9010, 90, 569): ['CBEAR', self.read_fake],
            (9110, 91, 570): ['PBEAR', self.read_fake],
            (9310, 93, 633): ['ROTPARM', self.read_fake],
            (9010, 90, 569): ['CBEAR', self.read_fake],

            (9407, 94, 659): ['ACADAPT', self.read_fake],
            (11701, 117, 656): ['CAMPBLL', self.read_campbll],

            (2601, 26, 58): ['FRFFLEX', self.read_fake],
            (2807, 28, 79): ['FRFOMAP', self.read_fake],
            (2701, 27, 62): ['FRFOTM', self.read_fake],
            (3501, 35, 56): ['RCROSSC', self.read_fake],
            (5407, 54, 649): ['ALOAD', self.read_fake],
            (5807, 59, 653): ['ACPLNW', self.read_acplnw],
            (7307, 73, 647): ['TLOAD3', self.read_fake],
            (9607, 96, 660): ['ACORDER', self.read_fake],
            (11801, 118, 657): ['NLHARM', self.read_fake],
            (12001, 120, 661): ['NLFREQ1', self.read_nlfreq1],
            (12201, 122, 667): ['JCONADD', self.read_fake],
            (13207, 132, 697): ['DAMPING', self.read_fake],
            (13007, 130, 676): ['FRFXIT1', self.read_fake],

            (12607, 126, 672): ['???-1', self.read_fake],
            (14207, 142, 1001): ['NLFREQ1?', self.read_fake],

            (12707, 127, 673): ['SPC?', self.read_fake],
            (12907, 129, 675): ['FRFXIT', self.read_fake],
            #(12607, 126, 672): ['???', self.read_stop],
        }

    def read_nlfreq1(self, data: bytes, n: int) -> int:
        """
        NLFREQ1, 20, 0.0, 0.01, 150
        data = (20, 0.0, 0.01, 150)
        """
        op2: OP2Geom = self.op2
        op2.log.warning('geom skipping NLFREQ1')
        return len(data)

    def read_campbll(self, data: bytes, n: int) -> int:
        """
        CAMPBLL
        CAMPBLL CID    VPARM DDVALID TYPE
                MODTRK CORU  SWITR   NUMMOD PRTCOR
        CAMPBLL,15,SPEED,22,RPM
        ints    = (15, 1, 22, RPM, 0, 0, 0, 0, 0)
        floats  = (15, 1, 22, RPM, 0.0, 0.0, 0.0, 0.0, 0.0)

                  (15, 1, 22, 'FREQ    ', 1, 0.85, 0, 0, 1)
        """
        op2: OP2Geom = self.op2
        assert self.size == 4, f'CAMPBLL size={self.size}'
        ntotal = 40 * self.factor
        nentries = (len(data) - n) // ntotal
        op2.increase_card_count('CAMPBLL', nentries)
        struc = Struct(op2._endian + b'3i 8s if 3i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            #(sid,p,c,a) = out
            #op2.show_data(edata, types='ifs')
            campbell_id, variable_param, ddval, cambell_type, modtrak, corr_threshold, switch, num_modes, print_flag = out
            #print(campbell_id, variable_param, ddval, cambell_type, modtrak, corr_threshold, switch, num_modes)
            assert isinstance(corr_threshold, float), corr_threshold
            check = (switch, num_modes)
            assert check == (0, 0), check
            assert variable_param == 1, 'variable_param=1 (speed?)'
            assert modtrak in {0, 1}, modtrak
            assert print_flag in {0, 1}, print_flag

            # darea = DAREA.add_op2_data(data=out)
            # op2._add_methods._add_darea_object(darea)
            n += ntotal
        op2.log.warning('geom skipping CAMPBLL')
        return n

    def _read_aload(self, data: bytes, n: int) -> int:
        """
        Record - ALOAD(5407,54,649)
        Combination of acoustic loads
        Word Name Type Description
        1 SID         I Load set identification number
        2 UNDEF(3) None N/A
        5 LI          I Load set identification number i
        Word 5 repeats until -1 occurs
        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data, dtype=op2.idtype8)
        asdf
        return len(data)

    def read_acplnw(self, data: bytes, n: int) -> int:
        """
        Record - ACPLNW(5807,59,653)

        Acoustic plane wave source
        Word      Name Type Description
        1 SID      I Load set identification number
        2 FORM     I Complex format: 0 for real/imaginary, 1 for magnitude/phase
        3 A1      RS Scale factor 1
        4 TRTM    RS Real part or magnitude of plane wave
        5 TIDTRTM  I TABLEDi bulk entry identification number for real part or magnitude of plane wave
        6 TITP    RS Imaginary part or phase (in degrees) of plane wave
        7 TIDTITP  I TABLEDi bulk entry identification number for imaginary part or phase (in degrees) of plane wave
        8 CID1     I Coordinate system identification number for location of plane wave source
        9  X      RS X-coordinate of plane wave source
        10 Y      RS Y-coordinate of plane wave source
        11 Z      RS Z-coordinate of plane wave source
        12 CID2    I Coordinate system identification number for direction of plane wave
        13 NX     RS Direction cosine between plane wave direction and X-axis
        14 NY     RS Direction cosine between plane wave direction and Y-axis
        15 NZ     RS Direction cosine between plane wave direction and Z-axis
        16 UNDEF(3) None
        """
        op2: OP2Geom = self.op2

        ntotal = 18 * self.size
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0, 'ACPLNW'

        loads = []
        struc = Struct(mapfmt(op2._endian + b'2i 2f if 2i 3f i 3f 3i', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, form, scale, real_float, tid_real, imag_float, tid_imag,
             cid1, x, y, z, cid2, nx, ny, nz, dummy1, dummy2, dummy3) = out
            assert sid > 0, sid
            assert form in {0, 1}, form

            form_str = 'REAL' if form == 0 else 'PHASE'

            real = real_float if tid_real == 0 else tid_real
            imag = imag_float if tid_imag == 0 else tid_imag

            assert cid1 >= 0, cid1
            assert cid2 >= 0, cid2
            assert cid2 >= 0, cid2
            assert (dummy1, dummy2, dummy3) == (0, 0, 0), (dummy1, dummy2, dummy3)
            xyz = [x, y, z]
            nxyz = [nx, ny, nz]


            # ACPLNW SID  FORM A  TR/TM TI/TP
            #        CID1 X1   X2 X3    CID2 N1 N2 N3
            op2.add_acplnw(sid, form_str, scale, real, imag, cid1, xyz, cid2, nxyz)
            #acsrce = ACSRCE(sid, excite_id, rho, b, delay=delay, dphase=dphase, power=0)
            n += ntotal
            #loads.append(acsrce)
        return n # , loads

    def read_acsrce(self, data: bytes, n: int) -> int:
        """common method for reading NX/MSC ACSRCE"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(data, n, self._read_acsrce_nx, self._read_acsrce_msc,
                                             'ACSRCE', op2._add_methods._add_dload_entry)
        return n

    def _read_acsrce_msc(self, data, n):
        """
        ACSRCE(5307,53,379)

        MSC 2005r2 - 2016.1
        -------------------
        Word Name Type Description
        1 SID    I Load set identification number
        2 DAREA  I DAREA Bulk Data entry identification number
        3 DPHASE I DPHASE Bulk Data entry identification number
        4 DELAY  I DELAY Bulk Data entry identification number
        5 TC     I TABLEDi Bulk Data entry identification number for C(f)
        6 RHO   RS Density of the fluid
        7 B     RS Bulk modulus of the fluid
        8 T     RS Time delay
        9 PH    RS Phase lead
        """
        op2: OP2Geom = self.op2
        ntotal = 36 * self.factor
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0, 'ACSRCE'

        loads = []
        struc = Struct(mapfmt(op2._endian + b'5i4f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, excite_id, dphase, delay, tc, rho, b, delay, dphase_float) = out
            assert sid > 0, sid
            assert dphase_float >= 0, dphase
            assert delay >= 0, delay
            assert tc >= 0, tc
            assert rho > 0, rho
            assert b > 0, b
            if dphase == 0:
                dphase = dphase_float
            acsrce = ACSRCE(sid, excite_id, rho, b, delay=delay, dphase=dphase, power=0)
            n += ntotal
            loads.append(acsrce)
        return n, loads

    def _read_acsrce_nx(self, data, n):
        """
        NX 12 - 2019.2
        --------------
        Word Name Type Description
        1 SID    I Load set identification number
        2 DAREA  I DAREA Bulk Data entry identification number
        3 DPHASE I DPHASE Bulk Data entry identification number
        4 DELAY  I DELAY Bulk Data entry identification number
        5 TC     I TABLEDi Bulk Data entry identification number for C(f)
        6 RHO     RS Density of the fluid
        7 B       RS Bulk modulus of the fluid
        8 DELAYR  RS If DELAYI  = 0, constant value for delay
        9 DPHASER RS If DPHASEI = 0, constant value for phase
        10 TCR    RS If TCI     = 0, constant value for C(f)

        """
        op2: OP2Geom = self.op2
        ntotal = 40 * self.factor
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0, 'ACSRCE'

        loads = []
        struc = Struct(mapfmt(op2._endian + b'5i5f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, excite_id, dphase, delay, tc, rho, b, delay_float, dphase_float, tc_float) = out
            assert sid > 0, sid
            assert rho > 0, rho
            assert b > 0, b
            if dphase == 0:
                dphase = dphase_float
            if delay == 0:
                delay = delay_float
            if tc == 0:
                tc = tc_float
            acsrce = ACSRCE(sid, excite_id, rho, b, delay=delay, dphase=dphase, power=0)
            n += ntotal
            loads.append(acsrce)
        return n, loads

    def read_darea(self, data: bytes, n: int) -> int:
        """
        DAREA(27,17,182) - the marker for Record 2

        1 SID I  Load set identification number
        2 P   I  Grid, scalar, or extra point identification number
        3 C   I  Component number
        4 A   RS Scale factor

        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor
        nentries = (len(data) - n) // ntotal
        op2.increase_card_count('DAREA', nentries)
        struc = Struct(mapfmt(op2._endian + b'3if', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            #(sid,p,c,a) = out
            darea = DAREA.add_op2_data(data=out)
            op2._add_methods._add_darea_object(darea)
            n += ntotal
        return n

    def read_delay(self, data: bytes, n: int) -> int:
        """
        DELAY(37,18,183) - Record 3

        1 SID I  Set identification number
        2 P   I  Grid, scalar, or extra point identification number
        3 C   I  Component number
        4 T   RS Time delay

        """
        op2: OP2Geom = self.op2
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        op2.increase_card_count('DELAY', nentries)
        struc = Struct(op2._endian + b'3if')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, nodes, components, delays = out
            if op2.is_debug_file:
                op2.binary_debug.write('  DELAY=%s\n' % str(out))
            delay = DELAY(sid, nodes, components, delays)
            op2._add_methods._add_delay_object(delay)
            n += ntotal
        return n

    def read_dload(self, data: bytes, n: int) -> int:
        """
        DLOAD(57,5,123) - Record 4

        1 SID  I Load set identification number
        2  S   RS Overall scale factor
        3  SI  RS Scale factor i
        4  LI  I Load set identification number i
        Words 3 through 4 repeat until (-1,-1) occurs

        """
        op2: OP2Geom = self.op2
        #ndata = len(data)
        #nfields = (ndata - n) // 4

        datan = data[n:]
        ints = np.frombuffer(datan, op2.idtype8).copy()
        floats = np.frombuffer(datan, op2.fdtype8).copy()
        istart = 0
        iminus1_delta = get_iend_from_ints(ints)
        nentries = 0
        for iend in iminus1_delta:
            #datai = data[n+istart*4 : n+iend*4]
            sid = ints[istart]
            global_scale = floats[istart + 1]
            #print('  sid=%s global_scale=%s' % (sid, global_scale))
            deltai = iend - istart - 2 # subtract 2 for sid, global scale
            assert deltai % 2 == 0, (op2.show_data(data[n+istart*4 : n+iend*4], 'if'))
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
            if op2.is_debug_file:
                op2.binary_debug.write('  DLOAD=%s\n' % str(out))

            op2.add_dload(sid, global_scale, scales, load_ids)
            istart = iend + 2
            nentries += 1
        op2.increase_card_count('DLOAD', nentries)
        return len(data)

    def read_dphase(self, data: bytes, n: int) -> int:
        """
        DPHASE(77,19,184) - Record 5

        1 SID I Load set identification number
        2 P   I Grid, scalar, or extra point identification number
        3 C   I Component number
        4 TH  RS Phase lead

        """
        op2: OP2Geom = self.op2
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        op2.increase_card_count('DPHASE', nentries)
        struc = Struct(op2._endian + b'3if')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, nodes, components, delays = out
            if op2.is_debug_file:
                op2.binary_debug.write('  DPHASE=%s\n' % str(out))
            delay = DPHASE(sid, nodes, components, delays)
            op2._add_methods._add_dphase_object(delay)
            n += ntotal
        return n

    def read_dynred(self, data: bytes, n: int) -> int:
        """
        DYNRED(4807,48,306)

        Word Name Type Description
        1 SID   I Load set identification number
        2 FMAX RS Highest frequency of interest
        3 NIRV  I Number of initial random vectors
        4 NIT   I Number of iterations
        5 IDIR  I Starting point to generate initial random vectors
        6 NQDES I Number of generalized degrees-of-freedom
        7 UNDEF(2 ) none
        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\n2h39.op2
        #DYNRED  202     1000.

        #C:\NASA\m4\formats\git\examples\move_tpl\gdr20w.op2
        #DYNRED  24      160.            10                              YES
        #DYNRED  23      160.            10              75
        #(23, 160.0, 6, 10, 0, 75, 75, 0)
        #(24, 160.0, 6, 10, 0,  0,  0, YES)
        ntotal = 32 # 4*8
        #op2.show_data(data[n:])
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        #struc = Struct(op2._endian + b'i f 5i 4s')
        struc = Struct(op2._endian + b'i f 5i i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            #print(out)
            (sid, fmax, nirv, nit, idir, nqdes, unused_zero1, zero_yes) = out
            if zero_yes == 542328153:  # integer version of yes
                zero_yes = 'YES'
            #print(sid, fmax, nirv, nit, idir, nqdes, unused_zero1, zero_yes)
            #assert unused_zero1 == 0, out
            #assert unused_yes == b'YES ', out
            if op2.is_debug_file:
                op2.binary_debug.write('  DYNRED=%s\n' % str(out))
            op2.add_dynred(sid, fmax, nirv, nit, idir, nqdes)
            n += ntotal
        op2.increase_card_count('DYNRED', nentries)
        return n


    def read_eigb(self, data: bytes, n: int) -> int:
        """EIGB(107,1,86) - Record 7

        NX
        Word Name Type Description
        1 SID           I Set identification number
        2 METHOD(2) CHAR4 Method of eigenvalue extraction
        4 L1           RS Lower bound of eigenvalue range of interest
        5 L2           RS Upper bound of eigenvalue range of interest
        6 NEP           I Estimate of number of roots in positive range
        7 NDP           I Desired number of positive roots
        8 NDN           I Desired number of negative roots
        9 UNDEF None
        10 NORM(2)  CHAR4 Method for normalizing eigenvectors
        12 G            I Grid or scalar point identification number
        13 C            I Component number
        14 UNDEF(5) None
        """
        op2: OP2Geom = self.op2
        #ntotal = 60 * self.factor  # 4*15
        ntotal = 72 * self.factor  # 4*18
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries
        op2.increase_card_count('EIGB', nentries)
        if self.size == 4:
            struc = Struct(op2._endian + b'i 8s ff 3i i 8s 7i')  # 4i -> 7i
        else:
            struc = Struct(op2._endian + b'q 16s dd 3q q 16s 7q')  # 4q -> 7q
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            #op2.show_data(edata[44:])
            # int, 8s, 2f, 3i, i, 8s, 4i
            out = struc.unpack(edata)
            sid, method, L1, L2, nep, ndp, ndn, dunno, norm, g, c, dunno_a, dunno_b, dunno_c, dunno_d, dunno_e = out
            if op2.is_debug_file:
                op2.binary_debug.write('EIGB=%s\n' % str(out))
            #print('out = %s' % str(out))
            method = method.strip().decode('latin1')
            norm = norm.strip().decode('latin1')
            eigb = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, g, c)
            op2._add_methods._add_method_object(eigb)
            n += ntotal
        return n

    def read_eigc(self, data: bytes, n: int) -> int:
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

        data  = (2, CLAN, '', MAX, '', 0,     0, 1e-08, 32, -1,
                 3, HESS, '', MAX, '', 0,     0, 1e-08, 32, -1)
        """
        op2: OP2Geom = self.op2
        if 0:
            ntotal = 60
            nentries = (len(data) - n) // ntotal
            op2.increase_card_count('EIGB', nentries)
            struct1 = Struct(op2._endian + b'i 4s 4s 2i f 2i')
            #struct2 = Struct(op2._endian + b'5f2i')
            for i in range(nentries):
                edata = data[n:n+ntotal]
                #op2.show_data(edata[44:])
                # int, 8s, 2f, 3i, i, 8s, 4i
                out = struct1.unpack(edata)
                sid, method, L1, L2, nep, ndp, ndn, dunno, norm, g, c, dunno_a, dunno_b = out
                if op2.is_debug_file:
                    op2.binary_debug.write('EIGC=%s\n' % str(out))
                #print('out = %s' % str(out))
                method = method.strip().decode('latin1')
                norm = norm.strip().decode('latin1')
                eigb = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, g, c)
                op2._add_methods._add_method_object(eigb)
                n += ntotal
            #return n

        #-------------------------------------------------
        ndata = len(data)
        nfields = (ndata - n) // self.size
        datan = data[n:]
        if self.size == 4:
            ints = unpack(op2._endian + b'%ii' % nfields, datan)
            floats = unpack(op2._endian + b'%if' % nfields, datan)
            strings = unpack(op2._endian + b'4s'* nfields, datan)
        else:
            ints = unpack(op2._endian + b'%iq' % nfields, datan)
            floats = unpack(op2._endian + b'%id' % nfields, datan)
            strings = unpack(op2._endian + b'8s'* nfields, datan)

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

            #print('sid=%s method=%r norm=%r grid=%s component=%s epsilon=%g, '
                  #'neigenvalues=%s ctlflag=%s' % (
                      #sid, method, norm, grid, component, epsilon, neigenvalues, control_flag))

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
                i += 10
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

                # +10 is for the prefix; +7 is for the -1s
                i += 10 + 7 # 3 + 4 from (-1,-1,-1) and (sid,grid,comp,coeff)
            else:
                raise NotImplementedError('EIGC control_flag=%s' % control_flag)
            datai.extend([-1, -1, -1, -1, -1, -1, -1])  # creates a +7

            if op2.is_debug_file:
                op2.binary_debug.write('  EIGC=%s\n' % str(datai))

            if grid == 0:
                grid = None
                assert component == 0, component
                component = None

            eigc = op2.add_eigc(
                sid, method, grid, component, epsilon, neigenvalues, norm=norm,
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

            nentries += 1
            #print('--------------')
        op2.increase_card_count('EIGC', nentries)
        return len(data)

    def read_eigp(self, data: bytes, n: int) -> int:
        """
        EIGP(257,4,158) - Record 9

        1 SID    I Load set identification number
        2 ALPHA RS Location of pole on real axis
        3 OMEGA RS Location of pole on imaginary axis
        4 M      I Multiplicity of complex root at pole

        """
        op2: OP2Geom = self.op2
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        op2.increase_card_count('EIGP', nentries)
        struct1 = Struct('i2fi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            sid, alpha, omega, m = out
            if op2.is_debug_file:
                op2.binary_debug.write('EIGP=%s\n' % str(out))
            #print('out = %s' % str(out))
            op2.add_eigp(sid, alpha, omega, m, alpha2=None, omega2=None, m2=None)
            n += ntotal
        return n

    def read_eigr(self, data: bytes, n: int) -> int:
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
        op2: OP2Geom = self.op2
        ntotal = 72 * self.factor
        nentries = (len(data) - n) // ntotal
        if self.size == 4:
            struct1 = Struct('i 8s 2f 4i 8s 7i')
        else:
            struct1 = Struct('q 16s 2d 4q 16s 7q')

        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (sid, method, f1, f2, ne, nd, null_a, null_b, norm, g, c,
             null_c, null_d, null_e, null_f, null_g) = out
            if op2.is_debug_file:
                op2.binary_debug.write('EIGR=%s\n' % str(out))
            method = method.strip().decode('latin1')
            norm = norm.strip().decode('latin1')
            eigr = op2.add_eigr(sid, method=method, f1=f1, f2=f2, ne=ne, nd=nd,
                                 norm=norm, G=g, C=c, comment='')
            eigr.validate()
            n += ntotal
        op2.increase_card_count('EIGR', nentries)
        return n

    def read_eigrl(self, data: bytes, n: int) -> int:
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

        I think the EIGRL is just a bizarre card that adds 0s when there are
        blank lines after the card...
        # optistruct
        'EIGRL          5                      10                            MASS'

          strings   = (b'5 0   0 \n\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00MASS
                                                                                 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0',)
            ints    = (5, 0,   0,   10, 0,   0,   0,   0,   0,   1397965133,     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
            floats  = (5, 0.0, 0.0, 10, 0.0, 0.0, 0.0, 0.0, 0.0, 907333664768.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        """
        op2: OP2Geom = self.op2
        #op2.show_data(data[n:], types='ifs')
        ndata = len(data)
        size = self.size
        if size == 4:
            s = Struct('i 2f 3i f 2i 8s f i')
            nbytes = 52
            types = 'if'
        else:
            s = Struct('q 2d 3q d 2q 16s d i') #  TODO: why is there an i?
            nbytes = 100
            types = 'qds'

        nentries = 0
        is_none = False
        while n < ndata:
            #edata = data[n:n+46] # 11*4+2 = 46
            edata = data[n:n+nbytes] # 13*4 = 52
            out = s.unpack(edata)
            sid, v1, v2, nd, msglvl, maxset, shfscl, flag1, flag2, norm, alpha, nums = out

            #large: sid=12 v1=0.0 v2=0.0 nd=1 msglvl=0 maxset=0 shfscl=0.0 flag1=0 flag2=-1 norm=b'MASS\x00\x00\x00\x00' alpha=0.0 nums=0
            #small: sid=12 v1=0.0 v2=0.0 nd=1 msglvl=0 maxset=0 shfscl=0.0 flag1=0 flag2=-1 norm=b'MASS\x00\x00\x00\x00' alpha=0.0 nums=538976288

            #op2.log.warning(f'sid={sid} v1={v1} v2={v2} nd={nd} msglvl={msglvl} maxset={maxset} '
                            #f'shfscl={shfscl} flag1={flag1} flag2={flag2} norm={norm} alpha={alpha} nums={nums}')

            norm = norm.decode('latin1').rstrip('\x00 ')
            #print('self._nastran_format =', self._nastran_format)
            if nums == 0: # and self._nastran_format == 'nx':
                #op2.log.warning(f'sid={sid} v1={v1} v2={v2} nd={nd} msglvl={msglvl} maxset={maxset} '
                                #f'shfscl={shfscl} flag1={flag1} flag2={flag2} norm={norm} alpha={alpha} nums={nums}')
                #pass
                nums = 14
                is_none = True

            #elif nums is None:
                #nums = 14
                #nums_total = nums * 4
                #assert len(data) == nbytes, f'ndata={len(data)} nbytes={nbytes} size={size}'
                #edata2 = data[n+nbytes:n+nbytes+nums_total]
                #op2.show_data(edata2, types=types)
                #fi = unpack(mapfmt('%if' % nums, size), edata2)
                #print(out, fi)
                #nbytes += nums_total

            elif nums != 538976288:
                op2.log.warning(f'sid={sid} v1={v1} v2={v2} nd={nd} msglvl={msglvl} maxset={maxset} '
                                 f'shfscl={shfscl} flag1={flag1} flag2={flag2} norm={norm} alpha={alpha} nums={nums}')
                assert nums < 10000, nums
                #edata2 = data[n+46:n+46+nums*4]
                edata2 = data[n+nbytes:n+nbytes+nums*4]
                op2.show_data(edata2, types=types)
                fi = unpack(mapfmt('%if' % nums, size), edata2)
                op2.log.warning(f'  frequencies={fi}')
                raise NotImplementedError(nums)

            if op2.is_debug_file:
                op2.binary_debug.write('  EIGRL=%s\n' % str(out))
            options = []
            values = []

            assert flag1 in {-1, 0}, (flag1, flag2)
            assert flag2 in {-1, 0}, (flag1, flag2)
            v1 = v1 if flag1 == 0 else None
            v2 = v2 if flag2 == 0 else None

            op2.add_eigrl(sid, v1=v1, v2=v2, nd=nd, msglvl=msglvl,
                           maxset=maxset, shfscl=shfscl,
                           norm=norm, options=options, values=values)
            #print(eigrl)
            if nums == 538976288:
                n = len(data)
            else:
                n += nbytes + nums * 4
            nentries += 1
        op2.increase_card_count('EIGRL', nentries)
        if is_none:
            assert n == len(data), f'n={n} ndata={len(data)}'
        return n

    def read_epoint(self, data: bytes, n: int) -> int:
        """EPOINT(707,7,124) - Record 12"""
        op2: OP2Geom = self.op2
        npoints = (len(data) - n) // 4
        fmt = op2._endian + b'%ii' % npoints
        nids = unpack(fmt, data[n:])
        if op2.is_debug_file:
            op2.binary_debug.write('EPOINT=%s\n' % str(nids))
        epoint = EPOINTs.add_op2_data(list(nids))
        op2._add_methods._add_epoint_object(epoint)
        op2.increase_card_count('EPOINT', count_num=npoints)
        return len(data)

    def read_freq(self, data: bytes, n: int) -> int:
        """FREQ(1307,13,126) - Record 13"""
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data[n:], op2.idtype8).copy()
        floats = np.frombuffer(data[n:], op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]
        istart = 0
        for iend in iminus1:
            sid = ints[istart]
            freqs = floats[istart + 1:iend]
            istart = iend + 1
            op2.add_freq(sid, freqs)
        op2.increase_card_count('FREQ', count_num=len(iminus1))
        return len(data)

    def read_freq1(self, data: bytes, n: int) -> int:
        """
        FREQ1(1007,10,125) - Record 14

        1 SID I  Set identification number
        2 F1  RS First frequency
        3 DF  RS Frequency increment
        4 NDF I  Number of frequency increments

        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor
        nentries = (len(data) - n) // ntotal
        struc = Struct(mapfmt(op2._endian + b'iffi', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, df, ndf = out
            if op2.is_debug_file:
                op2.binary_debug.write('FREQ1=%s\n' % str(out))
            #print('out = %s' % str(out))
            freq = FREQ1(sid, f1, df, ndf=ndf)
            op2._add_methods._add_freq_object(freq)
            n += ntotal
        op2.increase_card_count('FREQ1', nentries)
        return n

    def read_freq2(self, data: bytes, n: int) -> int:
        """
        FREQ2(1107,11,166) - Record 15

        1 SID  I Set identification number
        2 F1  RS First frequency
        3 F2  RS Last frequency
        4 NF   I Number of logarithmic intervals
        """
        op2: OP2Geom = self.op2
        ntotal = 16 * self.factor
        nentries = (len(data) - n) // ntotal
        struc = Struct(mapfmt(op2._endian + b'iffi', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, nf = out
            if op2.is_debug_file:
                op2.binary_debug.write('  FREQ2=%s\n' % str(out))
            #print('out = %s' % str(out))
            op2.add_freq2(sid, f1, f2, nf)
            n += ntotal
        op2.increase_card_count('FREQ2', nentries)
        return n

    def read_freq3(self, data: bytes, n: int) -> int:
        """
        FREQ3(1407,14,39) - Record 16

        1 SID      I Set identification number
        2 F1      RS Lower bound of modal frequency range
        3 F2      RS Upper bound of modal frequency range
        4 TYPE CHAR4 Type of interpolation: LINE or LOG
        5 NEF      I Number of frequencies
        6 BIAS    RS Clustering bias parameter

        """
        op2: OP2Geom = self.op2
        ntotal = 24 # 4*6
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'i 2f 4s if')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, freq_type, nef, bias = out
            freq_type = freq_type.strip().decode('latin1')
            if freq_type == 'LINE':
                freq_type = 'LINEAR'

            if op2.is_debug_file:
                op2.binary_debug.write('  FREQ3=%s\n' % str(out))

            op2.add_freq3(sid, f1, f2=f2, Type=freq_type, nef=nef, cluster=bias)
            n += ntotal
        op2.increase_card_count('FREQ3', nentries)
        return n

    def read_freq4(self, data: bytes, n: int) -> int:
        """
        FREQ4(1507,15,40) - Record 17

        1 SID   I Set identification number
        2 F1   RS Lower bound of modal frequency range
        3 F2   RS Upper bound of modal frequency range
        4 FSPD RS Frequency spread
        5 NFM   I Number of evenly spaced frequencies per spread

        """
        op2: OP2Geom = self.op2
        ntotal = 20 # 4*5
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'i 3f i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, f1, f2, fspread, nfm = out
            if op2.is_debug_file:
                op2.binary_debug.write('  FREQ4=%s\n' % str(out))
            op2.add_freq4(sid, f1, f2, fspread=fspread, nfm=nfm)
            n += ntotal
        op2.increase_card_count('FREQ4', nentries)
        return n

    def read_freq5(self, data: bytes, n: int) -> int:
        """
        FREQ5(1607,16,41) - Record 18

        1 SID  I Load set identification number
        2 F1  RS Lower bound of modal frequency range
        3 F2  RS Upper bound of modal frequency range
        4 FRI RS Fractions of natural frequencies

        """
        op2: OP2Geom = self.op2
        ints = np.frombuffer(data, dtype='int32').copy()
        floats = np.frombuffer(data, dtype='float32').copy()
        i_minus_1s = np.where(ints == -1)[0]
        nentries = len(i_minus_1s)

        i0 = 0
        for i_minus_1 in i_minus_1s:
            sid = ints[i0]
            floatsi = floats[i0 + 1:i_minus_1]
            if op2.is_debug_file:
                op2.binary_debug.write('  FREQ5=(%s, %s)\n' % (sid, floatsi.tolist()))
            f1 = floatsi[0]
            f2 = floatsi[1]
            fractions = floatsi[2:]
            op2.add_freq5(sid, fractions, f1=f1, f2=f2)
            #print('freq =', freq)
            i0 = i_minus_1 + 1
        op2.increase_card_count('FREQ5', nentries)
        return len(data)

        #ntotal = 20 # 4*5
        #nentries = (len(data) - n) // ntotal
        #struc = Struct(op2._endian + b'i 3f')
        #for i in range(nentries):
            #edata = data[n:n+ntotal]
            #out = struc.unpack(edata)
            #sid, f1, f2, fspread, nfm = out
            #if op2.is_debug_file:
                #op2.binary_debug.write('  FREQ5=%s\n' % str(out))
            #freq = op2.add_freq4(sid, f1, f2, fspread=fspread, nfm=nfm)
            #n += ntotal
        #op2.increase_card_count('FREQ5', nentries)

    def read_nrlgap(self, data: bytes, n: int) -> int:
        r"""
        C:\NASA\m4\formats\git\examples\move_tpl\nlrgap2.op2
        NLRGAP  SID     GA      GB      PLANE   TABK     TABG    TABU RADIUS
        nlrgap  200     9       10      XY      -961     965     967     5.0
        (200, 9, 10, 101, 961, 965, 967, 5.0)

        (400, 1, 2, 200, 401, 420, 0, 0,
         500, 1, 2, 200, 501, 520, 0, 0)

        """
        op2: OP2Geom = self.op2
        ntotal = 32 # 4*8
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0
        assert ndatai % ntotal == 0
        struc = Struct(op2._endian + b'7i f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            #sid, ga, gb, plane, table_k, table_g, table_u, radius = out
            n += ntotal
        return len(data)

    def read_nlrsfd(self, data: bytes, n: int) -> int:
        """this card is faked..."""
        card_name = 'NLRSFD'
        card_obj = None
        methods = {
            80 : self._read_nlrsfd_80,
            128 : self._read_nlrsfd_128,
        }
        #try:
        #n = self._read_rgyro_52(RGYRO, data, n)
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_double_card(
            card_name, card_obj,
            self.op2.reader_geom3._add_op2_rigid_element,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n

    def _read_nlrsfd_80(self, card_obj, data: bytes, n: int) -> tuple[int, list[Any]]:
        """
        NLRSFD(3807,38,505)

        Word Name Type Description
        1 SID          I Load set identification number
        2 GA           I Inner grid id
        3 GB           I Outer grid id
        4 PLANE(2) CHAR4 Radial gap orientation plan
        6 BDIA        RS Inner journal diameter
        7 BLEN        RS Damper length
        8 BCLR        RS Damper radial clearance
        9 SOLN(2)  CHAR4 Solution option
        11 VISCO      RS Lubricant viscosity
        12 PVAPCO     RS Lubricant vapor pressure
        13 NPORT       I Number of lubication ports
        14 PRES1      RS Boundary pressure for port 1
        15 THETA1     RS Angular position for port 1
        16 PRES2      RS Boundary pressure for port 2
        17 THETA2     RS Angular position for port 2
        18 NPNT        I Number of finite diff points
        19 OFFSET1    RS Offset in the SFD direction 1
        20 OFFSET2    RS Offset in the SFD direction 2

        data = (1, 101, 103, 'XY      ', 6.44, 0.727, 0.003, 'SHORT   ', 7.0e-7, 0.0, 1, 0.0, 270.0, 0.0, 0.0, 31, 0.0, 0.0)
        """
        assert self.factor == 1, self.factor
        op2: OP2Geom = self.op2

        ntotal = 80 * self.factor # 4*20
        struc = Struct(op2._endian + b'3i 8s 3f 8s 2f i 4f i 2f')

        op2.show_data(data[12:])
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, f'ndatai={ndatai}; nenteries={nentries}; leftover={ndatai % ntotal}'
        assert ndatai % ntotal == 0, f'ndatai={ndatai}; nenteries={nentries}; leftover={ndatai % ntotal}'
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, ga, gb, plane, bdia, blen, bclr, soln,
             visco, pvapco, nport,
             pres1, theta1, pres2, theat2, npnt,
             offset1, offset2) = out
            #print(out)

            plane = plane.rstrip().decode('latin1')
            soln = soln.rstrip().decode('latin1')
            #NLRSFD SID     GA     GB    PLANE BDIA   BLEN  BCLR   SOLN
            #       VISCO   PVAPCO NPORT PRES1 THETA1 PRES2 THETA2 NPNT
            #       OFFSET1 OFFSET2
            if op2.is_debug_file:
                op2.binary_debug.write('  NLRSFD=%s\n' % str(out))
            op2.add_nlrsfd(sid, ga, gb, plane, bdia, blen, bclr, soln,
                            visco, pvapco, nport,
                            pres1, theta1, pres2, theat2, npnt,
                            offset1, offset2)
            n += ntotal
        op2.increase_card_count('NLRSFD', nentries)
        return n, []

    def _read_nlrsfd_128(self, card_obj, data: bytes, n: int) -> tuple[int, list[Any]]:
        """
        NLRSFD(3807,38,505)

        Word Name Type Description
        1 SID          I Load set identification number
        2 GA           I Inner grid id
        3 GB           I Outer grid id
        4 PLANE(2) CHAR4 Radial gap orientation plan
        6 BDIA        RS Inner journal diameter
        7 BLEN        RS Damper length
        8 BCLR        RS Damper radial clearance
        9 SOLN(2)  CHAR4 Solution option
        11 VISCO      RS Lubricant viscosity
        12 PVAPCO     RS Lubricant vapor pressure
        13 NPORT       I Number of lubication ports
        14 PRES1      RS Boundary pressure for port 1
        15 THETA1     RS Angular position for port 1
        16 PRES2      RS Boundary pressure for port 2
        17 THETA2     RS Angular position for port 2
        18 NPNT        I Number of finite diff points
        19 OFFSET1    RS Offset in the SFD direction 1
        20 OFFSET2    RS Offset in the SFD direction 2

        ndata = 256:
          strings = (1, 901, 1600, 'XY      ', 7.6, 0.4, 0.003, 'SHORT   '\xa2\xe7;5\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x87C\x00\x00\x00\x00\x00\x00\x00\x00\x1f\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\xe9\x03\x00\x00\xdc\x05\x00\x00XY      \xcd\xcc\xcc@\xac\x1c:?\xa6\x9bD;SHORT   \xa2\xe7;5\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x87C\x00\x00\x00\x00\x00\x00\x00\x00\x1f\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',)
          ints    = (1, 901, 1600, 'XY      ', 7.6, 0.4, 0.003, 'SHORT   ', 7.e-7, 0,   1, 0,   270.0, 0,   0,   31, 0,   0,   538976288, 538976288, 538976288, 538976288, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1001, 1500, 538990936, 538976288, 1087163597, 1060773036, 994352038, 1380927571, 538976340, 893118370, 0, 1, 0, 1132920832, 0, 0, 31, 0, 0, 538976288, 538976288, 538976288, 538976288, 0, 0, 0, 0, 0, 0, 0, 0)
          floats  = (1, 901, 1600, 'XY      ', 7.6, 0.4, 0.003, 'SHORT   ', 7.e-7, 0.0, 1, 0.0, 270.0, 0.0, 0.0, 31, 0.0, 0.0, 1.3563156426940112e-19, 1.3563156426940112e-19, 1.3563156426940112e-19, 1.3563156426940112e-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.401298464324817e-45, 1.4026997627891419e-42, 2.1019476964872256e-42, 1.358208852320992e-19, 1.3563156426940112e-19, 6.400000095367432, 0.7269999980926514, 0.003000000026077032, 222567907328.0, 1.3563223635364882e-19, 6.999999868639861e-07, 0.0, 1.401298464324817e-45, 0.0, 270.0, 0.0, 0.0, 4.344025239406933e-44, 0.0, 0.0, 1.3563156426940112e-19, 1.3563156426940112e-19, 1.3563156426940112e-19, 1.3563156426940112e-19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        assert self.factor == 1, self.factor
        op2: OP2Geom = self.op2

        ntotal = 128 * self.factor # 4*32
        struc = Struct(op2._endian + b'3i 8s 3f 8s 2f i 4f i 2f 8s8s 8f')

        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0, f'ndatai={ndatai}; nenteries={nentries}; leftover={ndatai % ntotal}'
        assert ndatai % ntotal == 0, f'ndatai={ndatai}; nenteries={nentries}; leftover={ndatai % ntotal}'
        for unused_i in range(nentries):
            edata1 = data[n:n+ntotal]
            out = struc.unpack(edata1)
            (sid, ga, gb, plane, bdia, blen, bclr, soln,
             visco, pvapco, nport,
             pres1, theta1, pres2, theat2, npnt,
             offset1, offset2, word1_bytes, word2_bytes, *floats) = out
            word1 = word1_bytes.decode('latin1').rstrip()
            word2 = word2_bytes.decode('latin1').rstrip()
            #op2.show_data(edata1, types='ifs')
            #print(out)
            assert word1 in {'', 'TESTNLR'}, word1
            assert word2 in {'', 'NLRSFDA', 'NLRSFD'}, word2
            assert min(floats) >= 0, floats
            #print('blank', blank)
            #print('*other =', other)
            plane = plane.rstrip().decode('latin1')
            soln = soln.rstrip().decode('latin1')
            #NLRSFD SID     GA     GB    PLANE BDIA   BLEN  BCLR   SOLN
            #       VISCO   PVAPCO NPORT PRES1 THETA1 PRES2 THETA2 NPNT
            #       OFFSET1 OFFSET2
            if op2.is_debug_file:
                op2.binary_debug.write('  NLRSFD=%s\n' % str(out))
            op2.add_nlrsfd(sid, ga, gb, plane, bdia, blen, bclr, soln,
                            visco, pvapco, nport,
                            pres1, theta1, pres2, theat2, npnt,
                            offset1, offset2)
            n += ntotal
        op2.increase_card_count('NLRSFD', nentries)
        return n, []


    def read_nolin1(self, data: bytes, n: int) -> int:
        """
        NOLIN1(3107,31,127)

        Word Name Type Description
        1 SID I Load set identification number
        2 GI  I Grid, scalar, or extra point identification number of I
        3 CI  I Component number for GI.
        4 S  RS Scale factor
        5 GJ  I Grid, scalar, or extra point identification number of J
        6 CJ  I Component number for GJ
        7 T   I Identification number of a TABLEDi Bulk Data entry.
        8 UNDEF none

        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\d10903.op2
        #NOLIN1 SID GI CI S GJ CJ TID
        #NOLIN1  115     2               -1.0    2               2
        #NOLIN1  115     13              -1.0    13              3

        ntotal = 32 # 4*8
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'3if4i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, gi, ci, s, gj, cj, t, unused_zero) = out
            assert unused_zero == 0, out
            if op2.is_debug_file:
                op2.binary_debug.write('  NOLIN1=%s\n' % str(out))
            op2.add_nolin1(sid, gi, ci, s, gj, cj, t)
            n += ntotal
        op2.increase_card_count('NOLIN1', nentries)
        return n

    def read_nolin2(self, data: bytes, n: int) -> int:
        """
        NOLIN2(3207,32,128)

        Word Name Type Description
        1 SID I Load set identification number
        2 GI  I Grid, scalar, or extra point identification number of I
        3 CI  I Component number for GI.
        4 S  RS Scale factor
        5 GJ  I Grid, scalar, or extra point identification number of J
        6 CJ  I Component number for GJ
        7 GK  I Grid, scalar, or extra point identification number of K
        8 CK  I Component number for GK

        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\n12905b.op2
        #NOLIN2  SID     GI      CI         S    GJ      CJ      GK      CK
        #NOLIN2  400     2       1       -.01    3       0       2       1
        #NOLIN2  400     2       2       -.01    3       0       2       2
        #NOLIN2  400     3       0       -.005   2       1       2       1
        #NOLIN2  400     3       0       -.005   2       2       2       2

        ntotal = 32 # 4*8
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'3if4i')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, g, ci, s, gj, cj, gk, ck) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  NOLIN2=%s\n' % str(out))
            op2.add_nolin2(sid, g, ci, s, gj, cj, gk, ck)
            n += ntotal
        op2.increase_card_count('NOLIN2', nentries)
        return n

    def read_nolin3(self, data: bytes, n: int) -> int:
        """
        NOLIN3(3307,33,129)

        Word Name Type Description
        1 SID I Load set identification number
        2 GI  I Grid, scalar, or extra point identification number of I
        3 CI  I Component number for GI.
        4 S  RS Scale factor
        5 GJ  I Grid, scalar, or extra point identification number of J
        6 CJ  I Component number for GJ
        7 A  RS Exponent of the forcing function
        8 UNDEF none
        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\d10903.op2
        #NOLIN3  115     4               -3.5+3  4       10      2.
        #NOLIN3  115     11              -3.5+3  11      10      2.

        ntotal = 32 # 4*8
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'3if2ifi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, gi, ci, s, gj, cj, a, zero) = out
            assert zero == 0, out
            if op2.is_debug_file:
                op2.binary_debug.write('  NOLIN3=%s\n' % str(out))
            op2.add_nolin3(sid, gi, ci, s, gj, cj, a)
            n += ntotal
        op2.increase_card_count('NOLIN3', nentries)
        return n

    def read_nolin4(self, data: bytes, n: int) -> int:
        """
        NOLIN4(3407,34,130)

        Word Name Type Description
        1 SID I Load set identification number
        2 GI  I Grid, scalar, or extra point identification number of I
        3 CI  I Component number for GI.
        4 S  RS Scale factor
        5 GJ  I Grid, scalar, or extra point identification number of J
        6 CJ  I Component number for GJ
        7 A  RS Exponent of the forcing function
        8 UNDEF none
        """
        op2: OP2Geom = self.op2
        #C:\NASA\m4\formats\git\examples\move_tpl\d10903.op2
        #NOLIN4  115     4               -3.5+3  4       10      2.
        #NOLIN4  115     11              -3.5+3  11      10      2.

        ntotal = 32 # 4*8
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'3if2ifi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, gi, ci, s, gj, cj, a, zero) = out
            assert zero == 0, out
            if op2.is_debug_file:
                op2.binary_debug.write('  NOLIN4=%s\n' % str(out))
            op2.add_nolin4(sid, gi, ci, s, gj, cj, a)
            n += ntotal
        op2.increase_card_count('NOLIN4', nentries)
        return n

    def read_randps(self, data: bytes, n: int) -> int:
        """common method for reading NX/MSC RLOAD1"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(data, n, self._read_randps_nx, self._read_randps_msc,
                                             'RLOAD1', op2._add_methods._add_dload_entry)
        return n

    def _read_randps_nx(self, data: bytes, n: int) -> tuple[int, list]:
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
        op2: OP2Geom = self.op2
        ntotal = 28
        ndatai = (len(data) - n)
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0
        struc = Struct(op2._endian + b'3i 2f if')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, j, k, x, y, tidi, tidf = out
            tid = tidi
            if tidi == 0:
                tid = tidf
            if op2.is_debug_file:
                op2.binary_debug.write('  RANDPS=%s\n' % str(out))
            #op2.log.debug('  RANDPS=%s\n' % str(out))
            op2.add_randps(sid, j, k, x=x, y=y, tid=tid)
            n += ntotal
        op2.increase_card_count('RANDPS', nentries)
        return n, []

    def _read_randps_msc(self, data: bytes, n: int) -> tuple[int, list[None]]:
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
        op2: OP2Geom = self.op2
        ntotal = 24
        ndatai = (len(data) - n)
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0
        struc = Struct(op2._endian + b'3i2fi')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, j, k, x, y, tid = out
            if op2.is_debug_file:
                op2.binary_debug.write('  RANDPS=%s\n' % str(out))
            #op2.log.debug('  RANDPS=%s\n' % str(out))
            randps = op2.add_randps(sid, j, k, x=x, y=y, tid=tid)
            #cards.append(randps)
            n += ntotal
        op2.increase_card_count('RANDPS', nentries)
        return n, []

    def read_randt1(self, data: bytes, n: int) -> int:
        """
        RANDT1(2207,22,196)

        Word Name Type Description
        1 SID  I  Set identification number
        2 N    I  Number of time lag intervals
        3 TO   RS Starting time lag
        4 TMAX RS Maximum time lag

        """
        op2: OP2Geom = self.op2
        ntotal = 16  # 4*4
        struct1 = Struct(op2._endian + b'2i 2f')
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            if op2.is_debug_file:
                op2.binary_debug.write('  RANDT1=%s\n' % str(out))
            sid, nlags, to, tmax = out
            op2.add_randt1(sid, nlags, to, tmax)
            n += ntotal
        op2.card_count['RANDT1'] = nentries
        return n

    def read_rload1(self, data: bytes, n: int) -> int:
        """common method for reading NX/MSC RLOAD1"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(data, n, self._read_rload1_nx, self._read_rload1_msc,
                                             'RLOAD1', op2._add_methods._add_dload_entry)
        return n

    def _read_rload1_nx(self, data: bytes, n: int) -> tuple[int, list[RLOAD1]]:
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
        op2: OP2Geom = self.op2
        dloads = []
        ntotal = 44 * self.factor
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        struc = Struct(mapfmt(op2._endian + b'7i 4f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, dphasei, tci, tdi, load_type, delayr, dphaser, tcr, tdr = out # 44

            if op2.is_debug_file:
                op2.binary_debug.write('  RLOAD1=%s\n' % str(out))

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
        op2: OP2Geom = self.op2
        dloads = []
        ntotal = 36 * self.factor
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        struc = Struct(mapfmt(op2._endian + b'2i 2f 3i 2f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, dphaser, delayr, tci, tdi, load_type, tau, phi = out # 36

            if op2.is_debug_file:
                op2.binary_debug.write('  RLOAD1=%s\n' % str(out))

            tc = tci
            td = tdi
            delay = delayr
            dphase = dphaser
            dload = RLOAD1(sid, darea, delay=delay, dphase=dphase, tc=tc, td=td,
                           Type=load_type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def read_rload2(self, data: bytes, n: int) -> int:
        """common method for reading NX/MSC RLOAD2"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_rload2_nx, self._read_rload2_msc,
            'RLOAD2', op2._add_methods._add_dload_entry)
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
        op2: OP2Geom = self.op2
        dloads = []
        ntotal = 44 * self.factor
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0
        struc = Struct(mapfmt(op2._endian + b'7i 4f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, dphasei, tbi, tpi, load_type, delayr, dphaser, tbr, tpr = out
            if op2.is_debug_file:
                tau = 0.
                op2.binary_debug.write('  RLOAD2=%s\n' % str(out))
                assert sid > 0 and darea > 0, (f'RLOAD2; sid={sid} darea={darea} dphasei={dphasei} delayi={delayi} '
                                               f'tbi={tbi} tpi={tpi} load_type={load_type} tau={tau} dphasei={dphasei}')

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
        op2: OP2Geom = self.op2
        dloads = []
        ntotal = 36
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0
        struc = Struct(op2._endian + b'7i 2f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, dphasei, delayi, tbi, tpi, load_type, tau, phase = out
            if op2.is_debug_file:
                op2.binary_debug.write('  RLOAD2=%s\n' % str(out))
            assert sid > 0 and darea > 0, (f'RLOAD2; sid={sid} darea={darea} dphasei={dphasei} delayi={delayi} '
                                           f'tbi={tbi} tpi={tpi} load_type={load_type} tau={tau} phase={phase}')
            dload = RLOAD2(sid, darea, delay=delayi, dphase=dphasei, tb=tbi, tp=tpi,
                           Type=load_type, comment='')
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def read_rgyro(self, data: bytes, n: int) -> int:
        """this card is faked..."""
        card_name = 'RGYRO'
        card_obj = None
        methods = {
            36 : self._read_rgyro_36,
            52 : self._read_rgyro_52,
        }
        #try:
        #n = self._read_rgyro_52(RGYRO, data, n)
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_double_card(
            card_name, card_obj,
            self.op2.reader_geom3._add_op2_rigid_element,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n

    def _read_rgyro_36(self, card_obj, data: bytes, n: int) -> tuple[int, list]:
        """
        RGYRO(1507,15,40) - Record 17

        1 SID          I RGYRO identification number
        2 ATYPE(2) CHAR4 ASYNC/SYNC flag
        4 REFROT       I Reference rotor identification number
        5 UNIT(2)  CHAR4 RPM/FREQ flag for speed input
        7 SPDLOW      RS Lower limit of speed range
        8 SPDHGH      RS Upper limit of speed range
        9 SPEED       RS Specific speed

        RGYRO RID SYNCFLG REFROTR SPDUNIT SPDLOW SPDHIGH SPEED ROTRSEID
              WR3WRL WR4WRL WRHWRL
        RGYRO    1       SYNC      4      RPM                   1000.0
        """
        op2: OP2Geom = self.op2
        ntotal = 36 # 4*9
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'i 8s i 8s 3f')
        cards = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, asynci, refrot, unit, speed_low, speed_high, speed = out
            asynci = asynci.decode('latin1')
            unit = unit.decode('latin1')
            if op2.is_debug_file:
                op2.binary_debug.write('  RGYRO=%s\n' % str(out))
            op2.add_rgyro(sid, asynci, refrot, unit, speed_low, speed_high, speed)
            #RGYRO(sid, asynci, refrot, unit, speed_low, speed_high, speed)
            n += ntotal
        op2.increase_card_count('RGYRO', nentries)
        return n, cards

    def _read_rgyro_52(self, card_obj, data: bytes, n: int) -> tuple[int, list]:
        """
        RGYRO(10701, 107, 117) - Record 17

        1 SID            I RGYRO identification number
        2 ATYPE(2)   CHAR4 ASYNC/SYNC flag
        4 REFROT         I Reference rotor identification number
        5 UNIT(2)    CHAR4 RPM/FREQ flag for speed input
        7 SPDLOW        RS Lower limit of speed range
        8 SPDHGH        RS Upper limit of speed range
        9 SPEED         RS Specific speed
        10 rotor_seid    I ?
        11 whirl3       RS ?
        12 whirl4       RS ?
        13 whirl_hybrid RS ?

        RGYRO RID SYNCFLG REFROTR SPDUNIT SPDLOW SPDHIGH SPEED ROTRSEID
              WR3WRL WR4WRL WRHWRL
        RGYRO    1       SYNC      4      RPM                   1000.0

        ints    = (1, 206.349, 1.3563e-19  4, 541937746,  1.3563e-19, 0,   0,   1000.0, 0, 0, 0, 0)
        floats  = (1, 206.349, 1.3563e-19, 4, 1.7390e-19, 1.3563e-19, 0.0, 0.0, 1000.0, 0.0, 0.0, 0.0, 0.0)
        """
        cards = []
        op2: OP2Geom = self.op2
        ntotal = 52 # 4*13
        nentries = (len(data) - n) // ntotal
        struc = Struct(op2._endian + b'i 8s i 8s 3f i3f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, asynci, refrot, unit, speed_low, speed_high, speed,
             rotor_seid, whirl3, whirl4, whirl_hybrid) = out
            asynci = asynci.decode('latin1')
            unit = unit.decode('latin1')
            if op2.is_debug_file:
                op2.binary_debug.write('  RGYRO=%s\n' % str(out))
            op2.add_rgyro(sid, asynci, refrot, unit, speed_low, speed_high, speed)
            n += ntotal
        op2.increase_card_count('RGYRO', nentries)
        return n, cards

    def read_rotord(self, data: bytes, n: int) -> int:
        """
        ROTORD(8210, 82, 599)

        1 SID               I Set identification number
        2 NUMROT            I Number of rotors (i=1, NUMROT below)
        3 RSTART           RS Starting rotor speed (in RPM)
        4 RSTEP            RS Rotor speed step size (in RPM)
        5 NUMSTEP           I Number of steps for rotor speed
        6,7 REFSYS      CHAR8 Rotational reference system (YHROT,YH or YHFIX,YH
        8 CMOUT            RS Complex mode output request
        9,10 RUNIT      CHAR8 Revolution input/output units for rotor (YHRPM or YHRAD or YHCPS or YHHZ)
        11,12 FUNIT     CHAR8 Frequency output units (YHRPM or YHRAD or YHCPS or YHHZ)
        13,14 ZSTEIN    CHAR8 Flag to incorpoate Steiner inertia terms (YHYES,YH or YHNO,YH)
        15 ORBEPS          RS Threshold value for detection of whirl direction
        16 ROTPRT           I Optional printout flag
        17 SYNC             I Synchronous analysis selection flag
        18 ETYPE            I Excitation type
        19 EORDER          RS Excitation order
        20-21                 (not used)
        22+8*(i-1) RIDi     I Rotor ID for ith rotor
        23+8*(i-1) RSETi    I Set number for rotor speed for multiple rotors
        24+8*(i-1) RSPEEDi RS Relative rotor speed for multiple rotors
        25+8*(i-1) RCOORDI  I Coordinate system ID number specifying rotation axis as Z
        26+8*(i-1) W3i     RS Damping coefficient like PARAM,W3
        27+8*(i-1) W4i     RS Damping coefficient like PARAM,W4
        28+8*(i-1) RFORCEi I RFORCE bulk data ID number specifying rotational force to be applied
        29+8*(i-1) BRGSETi I CBEAR bearing set GROUP ID number

        verified in NX 2019
        """
        op2: OP2Geom = self.op2
        nentries = 0

        # up to EORDER
        ntotal = 84 # 19*4
        struct1 = Struct(op2._endian + b'2i 2f i 8s f 8s8s8s f 3i f 2i')

        ntotal2 = 32 # 8 * 4
        struct2 = Struct(op2._endian + b'2i f i 2f 2i')
        while n < len(data):
            #if op2.is_debug_file:
                #op2.binary_debug.write('  ROTORD=%s\n' % str(out))
            edata = data[n:n+ntotal]
            (sid, numrot, rstart, rstep, numstep, refsys, cmout, runit, funit,
             zstein, orbeps, rotprt, sync, etype, eorder, unused_a, unused_b) = struct1.unpack(edata)
            runit = runit.decode('latin1').strip()
            funit = funit.decode('latin1').strip()
            refsys = refsys.decode('latin1').strip()
            zstein = zstein.decode('latin1').strip()
            #print(sid, numrot, rstart, rstep, numstep, refsys, cmout,
                  #runit, funit, zstein, orbeps, rotprt, sync, etype, eorder)
            n += ntotal

            rids = []
            rsets = []
            rspeeds = []
            rcords = []
            w3s = []
            w4s = []
            rforces = []
            brgsets = []
            while n < len(data):
                edata = data[n:n+ntotal2]
                (rid, rset, rspeed, rcoord, w3, w4, rforce, brgset) = struct2.unpack(edata)
                #print((rid, rset, rspeed, rcoord, w3, w4, rforce, brgset))
                rids.append(rid)
                rsets.append(rset)
                rspeeds.append(rspeed)
                rcords.append(rcoord)
                w3s.append(w3)
                w4s.append(w4)
                rforces.append(rforce)
                brgsets.append(brgset)
                n += ntotal2

            op2.add_rotord(sid, rstart, rstep, numstep,
                            rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
                            refsys=refsys, cmout=cmout, runit=runit, funit=funit,
                            zstein=zstein, orbeps=orbeps, roprt=rotprt, sync=sync, etype=etype,
                            eorder=eorder, threshold=0.02, maxiter=10, comment='')
            nentries += 1
        op2.increase_card_count('ROTORD', nentries)
        return n

    def read_rotorg(self, data: bytes, n: int) -> int:
        """
        ROTORG(8410, 84, 600)

        1 SID I Set identification number of rotor
        2 ID  I Grid or scalar point ID's which define a rotor

        verified in NX 2019
        """
        op2: OP2Geom = self.op2
        nentries = 0
        ints = np.frombuffer(data[12:], dtype=op2.idtype)
        izero = np.where(ints==-1)[0]
        istart = 0
        for iend in izero:
            sid = ints[istart]
            grids = ints[istart+1:iend]
            assert -1 not in grids, grids.tolist()
            assert sid > 0, sid
            rotorg = op2.add_rotorg(sid, grids.tolist())
            str(rotorg)
            istart = iend + 1
            nentries += 1
        op2.increase_card_count('ROTORG', nentries)
        return len(data)

#RSPINR

    def read_rspint(self, data: bytes, n: int) -> int:
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_rspint_32, self._read_rspint_56,
            'RSPINT', self._add_rspint_obj)
        return n

    def _add_rspint_obj(self, rspints: list[int]) -> None:
        """TODO: remove this..."""
        pass

    def _read_rspint_32(self, data: bytes, n: int) -> tuple[int, list[None]]:
        r"""
        RSPINT(11001,110,310) - Record 31

        1 RID         I Rotor identification number
        2 GRIDA       I Grid A for rotation direction vector
        3 GRIDB       I Grid B for rotation direction vector
        4 GR         RS Rotor damping coefficient
        5 UNIT(2) CHAR4 RPM/FREQ flag for speed input
        7 TABLEID     I Table identification number for speed history
        8 ZERO ???

        C:\MSC.Software\msc_nastran_runs\r1d101m1.op2

        RSPINT   4       9       10     FREQ     20
         0.
        ints    = (4, 9, 10, 'FREQ    ', 20,
                   may be floats
                   0, 0, 0, 0, 0, 0, 0, 0)
        """
        op2: OP2Geom = self.op2
        #strings = (b'\xf9*\x00\x00n\x00\x00\x006\x01\x00\x00\x05\x00\x00\x00y\xe6\x00\x00~\xe6\x00\x00\x00\x00\x00\x00RPM     \xa0\x0f\x00\x00\x00\x00\x00\x00\x06\x00\x00\x00\x89\r\x01\x00\x97\r\x01\x00\x00\x00\x00\x00RPM     \xf8*\x00\x00\x00\x00\x00\x00',)
        #data = (11001, 110, 310,
                #5, 59001, 59006, 0, 541937746, 538976288, 4000, 0,
                #6, 69001, 69015, 0, 541937746, 538976288, 11000, 0)
        #op2.show_data(data)

        # per DMAP ???
        #ntotal = 28 * self.factor # 4*7
        #struc = Struct(op2._endian + b'3i f 8s i')  # per QRG

        # MSC
        # C:\NASA\m4\formats\git\examples\move_tpl\nltrot99.op2
        ntotal = 32 * self.factor # 4*7
        struc = Struct(op2._endian + b'3i f 8s 2i')  # per QRG
        assert self.size == 4, f'RSPINT size={self.size}'
        #op2.show_data(data[n:], types='ifs')
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0

        rspints = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            n += ntotal
            out = struc.unpack(edata)
            # rid, grida, gridb, gr, unit, table_id = out
            rid, grida, gridb, gr, unit, table_id, zero = out
            unit = unit.decode('latin1').rstrip()
            op2.log.debug(f'RSPINT: rid={rid} nids=[{grida}, {gridb}] gr={gr} unit={unit!r} table_id={table_id} zero={zero}')
            assert zero == 0
            if op2.is_debug_file:
                op2.binary_debug.write('  RSPINT=%s\n' % str(out))
            op2.add_rspint(rid, grida, gridb, gr, unit, table_id)
            rspint = None
            rspints.append(rspint)
        return n, rspints

    def _read_rspint_56(self, data: bytes, n: int) -> tuple[int, list[None]]:
        r"""
        RSPINT(11001,110,310) - Record 31

        1 RID         I Rotor identification number
        2 GRIDA       I Grid A for rotation direction vector
        3 GRIDB       I Grid B for rotation direction vector
        4 GR         RS Rotor damping coefficient
        5 UNIT(2) CHAR4 RPM/FREQ flag for speed input
        7 TABLEID     I Table identification number for speed history
        8 ZERO ???

        C:\MSC.Software\msc_nastran_runs\r1d101m1.op2

        RSPINT ROTORID GRIDA   GRIDB   SPDUNT SPTID SPDOUT ROTSEID
              GR       ALPHAR1 ALPHAR2 WR3R   WR4R  WRHR
        RSPINT   4       9       10     FREQ     20
         0.
        ints    = (4, 9, 10, 'FREQ    ', 20,
                   may be floats
                   0, 0, 0, 0, 0, 0, 0, 0)

        ROTORID : int
            Identification number of rotor. (Integer > 0; Required).
        GRIDA/GRIDB : int
            Positive rotor spin direction is defined from GRIDA to GRIDB.
        SPDUNIT : str
            Specifies whether the spin rates are given in terms of
             - RPM (revolutions/minute)
             - FREQ: frequency (revolutions(cycles)/sec).
        SPTID : int
            Rotor spin rate.
        SPDOUT : int
            EPOINT to output the rotor speed vs. time. Output will be in SPDUNITs
        GR : float; default=0.0
            Rotor structural damping factor.
        ALPHAR1 : float; default=0.0
            Scale factor applied to the rotor mass matrix for Rayleigh damping.
        ALPHAR2 : float; default=0.0
            Scale factor applied to the rotor stiffness matrix for Rayleigh damping.
        ROTRSEID : int; default=0
            Identification number of the superelement in which the rotor specified in the
            ROTORID field is defined.
        WR3R : float; default=0.0
            Specifies 'average' excitation frequency for calculation of rotor damping and
            circulation terms for rotor structural damping specified through GR field.
        RSPINT ROTORID GRIDA   GRIDB   SPDUNT SPTID SPDOUT ROTSEID
               GR      ALPHAR1 ALPHAR2 WR3R   WR4R  WRHR
        RSPINT 4       9       10      RPM    2
               0.1     1.0E-2  1.0E-6
        #RSPINT,4,9,10,RPM,2,
        #,0.1,1.0E-2,1.0E-6
        floats=['0.1', '0.01', '1e-06']

        """
        op2: OP2Geom = self.op2
        #strings = (b'\xf9*\x00\x00n\x00\x00\x006\x01\x00\x00\x05\x00\x00\x00y\xe6\x00\x00~\xe6\x00\x00\x00\x00\x00\x00RPM     \xa0\x0f\x00\x00\x00\x00\x00\x00\x06\x00\x00\x00\x89\r\x01\x00\x97\r\x01\x00\x00\x00\x00\x00RPM     \xf8*\x00\x00\x00\x00\x00\x00',)
        #data = (11001, 110, 310,
                #5, 59001, 59006, 0, 541937746, 538976288, 4000, 0,
                #6, 69001, 69015, 0, 541937746, 538976288, 11000, 0)
        #op2.show_data(data)

        # per DMAP ???
        #ntotal = 28 * self.factor # 4*7
        #struc = Struct(op2._endian + b'3i f 8s i')  # per QRG

        # MSC
        # C:\NASA\m4\formats\git\examples\move_tpl\nltrot99.op2
        ntotal = 56 * self.factor # 4*14
        struc = Struct(op2._endian + b'3i 8s i   2i 3f 3i')  # per QRG
        assert self.size == 4, f'RSPINT size={self.size}'
        #op2.show_data(data[n:], types='ifs')
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0

        rspints = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            n += ntotal
            out = struc.unpack(edata)
            #op2.show_data(edata[24:], types='ifs')
            (rid, grida, gridb, speed_unit_bytes, table_id, spd_out, zero,
             gr, alpha_r1, alpha_r2, wr3r, wr4r, wrhr) = out
            assert zero == 0, zero #  (SPDOUT, ROTSEID)
            speed_unit = speed_unit_bytes.decode('latin1').rstrip()
            #print('RSPINT', rid, grida, gridb, speed_unit, table_id, zero)
            op2.log.debug(f'RSPINT: rid={rid} nids=[{grida}, {gridb}] gr={gr:g} speed_unit={speed_unit!r} '
                          f'table_id={table_id} spd_out={spd_out} alpha_r=[{alpha_r1:g}, {alpha_r2:g}], '
                          f'wr=[{wr3r:g}, {wr4r:g}, {wrhr:g}]')
            if op2.is_debug_file:
                op2.binary_debug.write('  RSPINT=%s\n' % str(out))
            op2.add_rspint(rid, grida, gridb, gr, speed_unit, table_id)
            rspint = None
            rspints.append(rspint)
        return n, rspints

#SEQEP(5707,57,135)

    def read_tf(self, data: bytes, n: int) -> int:
        """TF"""
        op2: OP2Geom = self.op2
        # subtract of the header (sid, nid, component, b0, b1, b2)
        # divide by 5 (nid1, component1, a0, a1, a2)
        #nrows = (nfields - 6) // 5
        #print('n=%s nrows=%s' % (n, nrows))
        #print(op2.show_data(data))
        #nid1, component1, a0, a1, a2

        ndata = len(data)
        struct1 = Struct(op2._endian + b'3i3f')
        struct2 = Struct(op2._endian + b'2i3f')
        while n < ndata:
            n2 = n + 24 # 20=4*6
            sid, nid, component, b0, b1, b2 = struct1.unpack(data[n:n2])
            if op2.is_debug_file:
                op2.binary_debug.write('TF header -> %s\n' % ([sid, nid, component, b0, b1, b2]))

            nids = []
            components = []
            a = []
            irow = 0
            while 1:
                n3 = n2 + 20 # 20=4*5
                nid1, component1, a0, a1, a2 = struct2.unpack(data[n2:n3])

                if op2.is_debug_file:
                    op2.binary_debug.write('  i=%s     -> %s\n' % (
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
            #if op2.is_debug_file:
                #op2.binary_debug.write('%s\n' % str(tf))
            op2._add_methods._add_tf_object(tf)
            op2.increase_card_count('TF')
            n = n3
        return n

    def read_tic(self, data: bytes, n: int) -> int:
        """
        TIC(6607,66,137)

        1 SID I Load set identification number
        2 G   I Grid, scalar, or extra point identification number
        3 C   I Component number for point GD
        4 U0 RS Initial displacement
        5 V0 RS Initial velocity

        """
        op2: OP2Geom = self.op2
        ntotal = 20 *self.factor  # 5*4
        struct1 = Struct(mapfmt(op2._endian + b'3i 2f', self.size))
        nentries = (len(data) - n) // ntotal
        for unused_i in range(nentries):
            out = struct1.unpack(data[n:n+ntotal])
            if op2.is_debug_file:
                op2.binary_debug.write('  TIC=%s\n' % str(out))
            sid, nid, comp, u0, v0 = out
            op2.add_tic(sid, [nid], [comp], u0=u0, v0=v0)
            n += ntotal
        op2.card_count['TIC'] = nentries
        return n

#TIC3
    def read_tload1(self, data: bytes, n: int) -> int:
        """
        common method for reading NX/MSC TLOAD1

        (7107, 71, 138,
        100000001, 2, 0, 0, 100000000, 0, 0, 0, 0,  # 9
        100000002, 4, 0, 0, 100000000, 0, 0, 0, 0,
        100000004, 6, 0, 0, 100000001, 0, 0, 0, 0,
        100000005, 8, 0, 0, 100000001, 0, 0, 0, 0,
        100000007, 10, 0, 0, 100000002, 0, 0, 0, 0,
        100000008, 12, 0, 0, 100000002, 0, 0, 0, 0,
        100000010, 14, 0, 0, 100000003, 0, 0, 0, 0,
        100000011, 16, 0, 0, 100000003, 0, 0, 0, 0)
        """
        op2: OP2Geom = self.op2
        # NX - 24
        # MSC - 32
        #$       sid excite  delay    type    tid/f us vs
        #TLOAD1  5   2                         2    0. 0.
                   #sid excite ?    ?    tid ?    ?    ?    ?
        #ints    = (5,  2,     0,   0,   2,  0,   0,   0,   0)
        #floats  = (5,  2,     0.0, 0.0, 2,  0.0, 0.0, 0.0, 0.0)

        op2: OP2Geom = self.op2
        card_name = 'TLOAD1'
        card_obj = TLOAD1
        methods = {
            24 : self._read_tload1_nx_24,
            32 : self._read_tload1_msc_32,
            36 : self._read_tload1_36,
        }
        #try:
        n = op2.reader_geom2._read_double_card(
            card_name, card_obj, op2._add_methods._add_dload_entry,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n

        # ndatai = (len(data) - n) // self.factor
        # if ndatai % 24 == 0 and ndatai % 32 and ndatai % 36:
        #     n, dloads = self._read_tload1_nx_24(data, n)
        # elif ndatai % 32 == 0 and ndatai % 24 and ndatai % 36:
        #     n, dloads = self._read_tload1_msc_32(data, n)
        # elif ndatai % 36 == 0 and ndatai % 24 and ndatai % 32:
        #     n, dloads = self._read_tload1_36(data, n)
        # else:
        #     n = op2.reader_geom2._read_dual_card(
        #         data, n,
        #         self._read_tload1_nx_24, self._read_tload1_msc_32,
        #         'TLOAD1', op2._add_methods._add_dload_entry)
        #     return n
        # for dload in dloads:
        #     op2._add_methods._add_dload_entry(dload)
        # op2.card_count['TLOAD1'] = len(dloads)
        return n

    def _read_tload1_nx_24(self, card_obj, data: bytes, n: int) -> tuple[int, list[TLOAD1]]:
        """
        Record - TLOAD1(7107,71,138)

        Word Name Type Description
        1 SID     I Load set identification number
        2 DAREA   I DAREA Bulk Data entry identification number
        3 DELAYI  I DELAY Bulk Data entry identification number
        4 TYPE    I Nature of the dynamic excitation
        5 TID     I Identification number of TABLEDi entry that gives F(t)
        6 DELAYR RS If DELAYI = 0, constant value for delay

        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor # 6*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries

        dloads = []
        struc = Struct(mapfmt(op2._endian + b'5i f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, tid, delayr = out
            assert sid > 0, sid
            if op2.is_debug_file:
                op2.binary_debug.write('TLOAD1=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD1(sid, darea, tid, delay=delay, Type=load_type)
            #dload.validate()
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_tload1_msc_32(self, card_obj, data: bytes, n: int) -> tuple[int, list[TLOAD1]]:
        r"""
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

        # C:\MSC.Software\msc_nastran_runs\pbxsfsd.op2
        # why are there 9 fields?
        $       sid excite  delay    type    tid/f us vs
        TLOAD1  5   2                         2    0. 0.
                   sid excite ?    ?    tid ?    ?    ?    ?
        ints    = (5,  2,     0,   0,   2,  0,   0,   0,   0)
        floats  = (5,  2,     0.0, 0.0, 2,  0.0, 0.0, 0.0, 0.0)
        """
        op2: OP2Geom = self.op2
        ntotal = 32 * self.factor # 8*4
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries

        dloads = []
        struc = Struct(mapfmt(op2._endian + b'5i 3f', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, tid, delayr, us0, vs0 = out
            assert sid > 0, sid
            if op2.is_debug_file:
                op2.binary_debug.write('TLOAD1=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD1(sid, darea, tid, delay=delay, Type=load_type,
                           us0=us0, vs0=vs0)
            #dload.validate()
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_tload1_36(self, card_obj, data: bytes, n: int) -> tuple[int, list[TLOAD1]]:
        r"""
        TLOAD1(7107,71,138) - Record 37

        1 SID    I  Load set identification number
        2 DAREA  I  DAREA Bulk Data entry identification number
        3 DELAYI I  DELAY Bulk Data entry identification number
        4 TYPE   I  Nature of the dynamic excitation
        5 TID    I  Identification number of TABLEDi entry that gives F(t)
        6 DELAYR RS If DELAYI = 0, constant value for delay
        6 U0     RS Initial displacement factor for enforced motion (MSC; NX undocumented)
        7 V0     RS Initial velocity factor for enforced motion (MSC; NX undocumented)
        8 F      RS Value of F(t) for all times

        # C:\MSC.Software\msc_nastran_runs\pbxsfsd.op2
        # why are there 9 fields?
        $       sid excite  delay    type    tid/f us vs
        TLOAD1  5   2                         2    0. 0.
                   sid excite ?    ?    tid ?    ?    ?    ?
        ints    = (5,  2,     0,   0,   2,  0,   0,   0,   0)
        floats  = (5,  2,     0.0, 0.0, 2,  0.0, 0.0, 0.0, 0.0)

        C:\MSC.Software\msc_nastran_runs\mmc_112_rbe2.op2
        $	sid	excite	delay	type	tid/f	us0	vs0
        TLOAD1         2       2             45.       2
        data = (2, 2, 0, 0, 0, 2.0, 0.0, 0.0, f=45.0)
        """
        op2: OP2Geom = self.op2
        ntotal = 36 * self.factor # 8*4
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nentries > 0, nentries

        dloads = []
        struc = Struct(mapfmt(op2._endian + b'5i 3ff', self.size))
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, tid, delayr, us0, vs0, value = out
            assert sid > 0, sid
            #print(f'sid={sid} darea={darea} delayi={delayi} load_type={load_type} tid={tid} delayr={delayr} us0={us0} vs0={vs0} value={value}')
            if op2.is_debug_file:
                op2.binary_debug.write('TLOAD1=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            tid_f = tid
            if tid == 0:
                tid_f = value
            dload = TLOAD1(sid, darea, tid_f, delay=delay, Type=load_type,
                           us0=us0, vs0=vs0)
            #str(dload)
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def read_tload2(self, data: bytes, n: int) -> int:
        """common method for reading NX/MSC TLOAD2"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n,
            self._read_tload2_nx, self._read_tload2_msc,
            'TLOAD2', op2._add_methods._add_dload_entry)
        return n

    def _read_tload2_nx(self, data: bytes, n: int) -> tuple[int, list[TLOAD2]]:
        """
        TLOAD2(7207,72,139) - Record 37

        NX 2019.2
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

        My guess...
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
        op2: OP2Geom = self.op2
        ntotal = 44
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries

        dloads = []
        struc = Struct(op2._endian + b'4i 7f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, t1, t2, freq, p, c, growth, delayr = out
            if op2.is_debug_file:
                op2.binary_debug.write('  TLOAD2=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD2(sid, darea, delay=delay, Type=load_type, T1=t1,
                           T2=t2, frequency=freq,
                           phase=p, c=c, b=growth)
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def _read_tload2_msc(self, data: bytes, n: int) -> tuple[int, list[TLOAD2]]:
        """
        TLOAD2(7207,72,139) - Record 37

        MSC
        1 SID    I Load set identification number
        2 DAREA  I DAREA Bulk Data entry identification number
        3 DELAY  I DELAY Bulk Data entry identification number
        4 TYPE   I Nature of the dynamic excitation
        5 T1    RS Time constant 1
        6 T2    RS Time constant 2
        7 F     RS Frequency
        8 P     RS Phase angle
        9 C     RS Exponential coefficient
        10 B    RS Growth coefficient
        11 U0   RS Initial displacement factor for enforced motion
        12 V0   RS Initial velocity factor for enforced motion
        13 T    RS Time delay
        """
        op2: OP2Geom = self.op2
        ntotal = 52
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries

        dloads = []
        struc = Struct(op2._endian + b'4i 7f 2f')
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            sid, darea, delayi, load_type, t1, t2, freq, p, c, growth, delayr, us0, vs0 = out
            if op2.is_debug_file:
                op2.binary_debug.write('  TLOAD2=%s\n' % str(out))
            delay = delayi
            if delayi == 0:
                delay = delayr
            dload = TLOAD2(sid, darea, delay=delay, Type=load_type, T1=t1,
                           T2=t2, frequency=freq,
                           phase=p, c=c, b=growth,
                           us0=us0, vs0=vs0)
            dloads.append(dload)
            n += ntotal
        return n, dloads

    def read_tstep(self, data: bytes, n: int) -> int:
        """TSTEP(8307,83,142) - Record 38

        Word Name Type Description
        1 SID  I Set identification number
        2 N    I Number of time steps of value DTi
        3 DT  RS Time increment
        4 NO   I Skip factor for output
        Words 2 through 4 repeat until (-1,-1,-1) occurs
        """
        op2: OP2Geom = self.op2
        #op2.show_data(data)
        ndata = len(data)
        nfields = (ndata - n) // self.size
        datan = data[n:]
        #ints = np.frombuffer(datan, op2.idtype8).copy()
        #floats = np.frombuffer(datan, op2.fdtype8).copy()
        ints = unpack(mapfmt(op2._endian + b'%ii' % nfields, self.size), datan)
        floats = unpack(mapfmt(op2._endian + b'%if' % nfields, self.size), datan)
        # strings = unpack(op2._endian + b'4s'* nfields, datan)

        i = 0
        #nentries = 0
        ntimes = []
        dt = []
        no = []
        while i < nfields:
            sid = ints[i]
            ntimesi = ints[i + 1]
            dti = ints[i + 1]
            noi = ints[i + 3]
            dtf = floats[i + 2]
            i += 4
            n += 16
            while (ntimesi, dti, noi) != (-1, -1, -1):
                ntimes.append(ntimesi)
                dt.append(dtf)
                no.append(noi)
                #print('sid=%s ntimes=%s dt=%s no=%s' % (sid, ntimes, dt, no))
                ntimesi = ints[i]
                dti = ints[i + 1]
                noi = ints[i + 2]
                dtf = floats[i + 1]
                i += 3
                n += 12
            #print('sid=%s ntimes=%s dt=%s no=%s' % (sid, ntimes, dt, no))
            op2.add_tstep(sid, ntimes, dt, no)
        return len(data)

    def read_tstep1(self, data: bytes, n: int) -> int:
        """
        Record - TSTEP1(17500,175,618)

        Word Name Type Description
        1 SID   I Set identification number
        2 TEND RS End time
        3 NINC  I Number of load increments
        4 NOUT  I >0, output at each NOUT; =-1, YES; =-2, END;
        =-3, ALL
        Words 2 through 4 repeat until (-1,-1,-1) occurs
        """
        op2: OP2Geom = self.op2
        op2.log.info('skipping TSTEP1 in DYNAMICS')
        if op2.is_debug_file:
            op2.binary_debug.write('skipping TSTEP1 in DYNAMICS\n')
        return len(data)

#UNBALNC

    def read_rcross(self, data: bytes, n: int) -> int:
        """
        """
        #C:\NASA\m4\formats\git\examples\move_tpl\rcross01.op2
        #RCROSS  SID    RTYPE1   ID1     COMP1   RTYPE2  ID2     COMP2   CURID
        #RCROSS  100     FORCE   202     167     FORCE   202     167     1202001
        #RCROSS  100     FORCE   3302    6               3302    6       1330201
        #RCROSS  100     FORCE   202     167     FORCE   3302    6       1330202
        #RCROSS  100     FORCE   3302    6       FORCE   202     167     1330204
        #RCROSS  100     FORCE   6402    6               6402    6       1640201
        #RCROSS  100     FORCE   6412    6               6412    6       1641201
        #RCROSS  100     FORCE   6402    6               6412    6       1641202
        #RCROSS  100     FORCE   6412    6               6402    6       1641204
        #RCROSS  100     FORCE   6901    17              6901    17      1690101
        #RCROSS  100     STRESS  201     16      STRESS  201     16      2201001
        #RCROSS  100     STRESS  3301    7       STRESS  3301    7       2330101
        #RCROSS  100     STRESS  201     16      STRESS  3301    7       2330102
        #RCROSS  100     STRESS  3301    7       STRESS  201     16      2330104
        #RCROSS  100     STRESS  6401    7               6401    7       2640101
        #RCROSS  100     STRESS  6411    7               6411    7       2641101
        #RCROSS  100     STRESS  6701    86              6701    86      2670101
        #RCROSS  100     STRESS  7001    12              7001    12      2700101
        #RCROSS  100     STRESS  7401    12              7401    12      2740101
        #RCROSS  100     STRESS  7501    12              7501    12      2750101
        #RCROSS  100     STRESS  7511    12              7511    12      2751101
        #RCROSS  100     STRESS  9901    7               9901    7       2990101
        #RCROSS  100     STRESS  7511    12              9901    7       2990102
        #RCROSS  100     STRESS  9901    7               7511    12      2990104
        #RCROSS  100     FORCE   9902    6       STRESS  9901    7       2990106
        #RCROSS  100     STRESS  9901    7       FORCE   9902    6       2990108
        #RCROSS  100     DISP    3306    3       DISP    3306    3       3330601
        #RCROSS  100     DISP    3306    5       DISP    3306    5       3330602
        #RCROSS  100     DISP    6423    1               6423    1       3642301
        #RCROSS  100     DISP    6423    4               6423    4       3642302
        #RCROSS  100     DISP    3306    3       DISP    3306    5       3330603
        #RCROSS  100     DISP    3306    5       DISP    3306    3       3330604
        #RCROSS  100     DISP    3306    3               6423    4       3642303
        #RCROSS  100     DISP    6423    4               3306    3       3642304
        #RCROSS  100     SPCF    3305    3       SPCF    3305    3       4330501
        #RCROSS  100     SPCF    3305    5       SPCF    3305    5       4330502
        #RCROSS  100     SPCF    3305    3               3305    5       4330503
        #RCROSS  100     SPCF    3305    5               3305    3       4330504
        #RCROSS  100     SPCF    6413    1       SPCF    6413    1       4641301
        #RCROSS  100     SPCF    6413    4       SPCF    6413    4       4641302
        #RCROSS  100     SPCF    3305    3               6413    4       4641303
        #RCROSS  100     SPCF    6413    4               3305    3       4641304
        #RCROSS  100     SPCF    3305    3       DISP    3306    3       4641306
        #RCROSS  100     DISP    3306    3       SPCF    3305    3       4641308
        #RCROSS  200     SPCF    4444    2       DISP    9999    1       9999999
        op2: OP2Geom = self.op2

        #op2.show_data(data[n:])
        ntotal = 32 # 4*8
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        struc = Struct(op2._endian + b'i 4s 2i 4s 2i i')

        allowed = {'DISP', 'VELO', 'ACCEL', 'OLOAD', 'SPCF', 'MPCF', 'STRESS', 'STRAIN', 'FORCE', 'PRESS'}
        #DISP Displacement Vector
        #VELO Velocity Vector
        #ACCEL Acceleration Vector
        #OLOAD Applied Load Vector
        #SPCF Single-point Constraint Force Vector
        #MPCF Multi-point Constraint Force Vector
        #STRESS Element Stress
        #STRAIN Element Strain
        #FORCE Element Force

        map_type = {
            'DISP' : 'DISP',
            'VELO' : 'VELO',
            'SPCF' : 'SPCF',
            'MPCF' : 'MPCF',

            'ACCE' : 'ACCEL',
            'FORC' : 'FORCE',
            'STRE' : 'STRESS',
            'PRES' : 'PRESS',
            'STRA' : 'STRAIN',
            'OLOA' : 'OLOAD',
        }
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = struc.unpack(edata)
            (sid, rtype1, id1, comp1, rtype2, id2, comp2, curid) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  RCROSS=%s\n' % str(out))

            rtype1 = map_type[rtype1.decode('latin1').rstrip()]
            rtype2 = map_type[rtype2.decode('latin1').rstrip()]
            assert rtype1 in allowed, rtype1
            assert rtype2 in allowed, rtype2
            op2.add_rcross(sid, rtype1, id1, comp1, rtype2, id2, comp2, curid)
            n += ntotal
        op2.increase_card_count('RCROSS', nentries)
        return n
