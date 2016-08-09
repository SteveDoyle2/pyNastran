#pylint: disable=C0111,C0103,C0301,W0612,W0613,R0914,R0201
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack, Struct
from six import b
from six.moves import range

from pyNastran.bdf.cards.materials import (CREEP, MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11, MATHP)
from pyNastran.bdf.cards.material_deps import MATS1 # MATT1
from pyNastran.bdf.cards.dynamic import NLPARM, TSTEPNL # TSTEP
from pyNastran.op2.tables.geom.geom_common import GeomCommon
#from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 #PHBDY, CONV, CONVM, RADBC)
from pyNastran.bdf.cards.thermal.thermal import RADM


class MPT(GeomCommon):
    """defines methods for reading op2 materials & time-stepping methods"""

    def _read_mpt_4(self, data, ndata):
        return self._read_geom_4(self._mpt_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self.big_materials = {}
        self._mpt_map = {
            (1003, 10, 245): ['CREEP', self._read_creep],  # record 1
            (103, 1, 77): ['MAT1', self._read_mat1],       # record 3-msc-dmap2014
            (203, 2, 78): ['MAT2', self._read_mat2],       # record 3
            (1403, 14, 122): ['MAT3', self._read_mat3],    # record 4
            (2103, 21, 234): ['MAT4', self._read_mat4],    # record 5
            (2203, 22, 235): ['MAT5', self._read_mat5],    # record 6
            (2503, 25, 288): ['MAT8', self._read_mat8],    # record 7
            (2603, 26, 300): ['MAT9', self._read_mat9],    # record 8 - buggy
            (2801, 28, 365): ['MAT10', self._read_mat10],  # record 9
            (2903, 29, 371) : ['MAT11', self._read_mat11],  # record ??? - NX specific - buggy?

            #(2603, 26, 300): ['MAT9', self._read_fake],    # record 8 - buggy
            #(2903, 29, 371) : ['MAT11', self._read_fake],  # record ??? - NX specific - buggy?

            (4506, 45, 374): ['MATHP', self._read_mathp],   # record 11
            (503, 5, 90): ['MATS1', self._read_mats1],      # record 12
            (703, 7, 91): ['MATT1', self._read_matt1],      # record 13 - not done
            (803, 8, 102): ['MATT2', self._read_matt2],     # record 14 - not done
            (1503, 14, 189): ['MATT3', self._read_matt3],   # record 15 - not done
            (2303, 23, 237): ['MATT4', self._read_matt4],   # record 16 - not done
            (2403, 24, 238): ['MATT5', self._read_matt5],   # record 17 - not done
            (2703, 27, 301): ['MATT9', self._read_matt9],   # record 19 - not done
            (8802, 88, 413): ['RADM', self._read_radm],     # record 25 - not done
            # record 26
            (3003, 30, 286): ['NLPARM', self._read_nlparm],   # record 27
            (3104, 32, 350): ['NLPCI', self._read_nlpci],     # record 28
            (3103, 31, 337): ['TSTEPNL', self._read_tstepnl], # record 29
            (3303, 33, 988) : ['MATT11', self._read_fake],

            (903, 9, 336) : ['MATT8', self._read_fake],
            (8902, 89, 423) : ['RADMT', self._read_fake],
            (9002, 90, 410) : ['RADBND', self._read_fake],
        }

    def add_op2_material(self, mat):
        if mat.mid > 100000000:
            raise RuntimeError('bad parsing...')
        self.add_structural_material(mat, allow_overwrites=True)
        #print(str(mat)[:-1])

    def _read_creep(self, data, n):
        """
        CREEP(1003,10,245) - record 1
        """
        nmaterials = (len(data) - n) // 64
        s = Struct(b(self._endian + 'i2f4ifi7f'))
        for i in range(nmaterials):
            edata = data[n:n+64]
            out = s.unpack(edata)
            (mid, T0, exp, form, tidkp, tidcp, tidcs, thresh,
             Type, ag1, ag2, ag3, ag4, ag5, ag6, ag7) = out
            mat = CREEP.add_op2_data(out)
            self.add_creep_material(mat, allow_overwrites=True)
            n += 64
        self.card_count['CREEP'] = nmaterials
        return n

    def _read_mat1(self, data, n):
        """
        MAT1(103,1,77) - record 2
        """
        ntotal = 48  # 12*4
        s = Struct(b(self._endian + 'i10fi'))
        nmaterials = (len(data) - n) // ntotal
        for i in range(nmaterials):
            edata = data[n:n+48]
            out = s.unpack(edata)
            (mid, E, G, nu, rho, A, TRef, ge, St, Sc, Ss, mcsid) = out
            mat = MAT1.add_op2_data(out)
            self.add_op2_material(mat)
            n += ntotal
        self.card_count['MAT1'] = nmaterials
        return n

    def _read_mat2(self, data, n):
        """
        MAT2(203,2,78) - record 3
        """
        ntotal = 68  # 17*4
        s = Struct(b(self._endian + 'i15fi'))
        nmaterials = (len(data) - n) // ntotal
        for i in range(nmaterials):
            edata = data[n:n+68]
            out = s.unpack(edata)
            (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
             TRef, ge, St, Sc, Ss, mcsid) = out
            #print("MAT2 = ",out)
            mat = MAT2.add_op2_data(out)

            if mid > 1e8 or mid < 0:  # just a checker for out of range materials
                self.big_materials[mid] = mat
            else:
                self.add_op2_material(mat)
            n += ntotal
        self.card_count['MAT2'] = nmaterials
        return n

    def _read_mat3(self, data, n):
        """
        MAT3(1403,14,122) - record 4
        """
        s = Struct(b(self._endian + 'i8fi5fi'))
        nmaterials = (len(data) - n) // 64
        for i in range(nmaterials):
            out = s.unpack(data[n:n+64])
            (mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
             blank, ax, ath, az, TRef, ge, blank) = out
            mat = MAT3.add_op2_data([mid, ex, eth, ez, nuxth, nuthz,
                                     nuzx, rho, gzx, ax, ath, az, TRef, ge])
            self.add_op2_material(mat)
            n += 64
        self.card_count['MAT3'] = nmaterials
        return n

    def _read_mat4(self, data, n):
        """
        MAT4(2103,21,234) - record 5
        """
        s = Struct(b(self._endian + 'i10f'))
        nmaterials = (len(data) - n) // 44
        for i in range(nmaterials):
            out = s.unpack(data[n:n+44])
            (mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
            mat = MAT4.add_op2_data(out)
            self.add_thermal_material(mat, allow_overwrites=True)
            n += 44
        self.card_count['MAT4'] = nmaterials
        return n

    def _read_mat5(self, data, n):
        """
        MAT5(2203,22,235) - record 6
        """
        s = Struct(b(self._endian + 'i9f'))
        nmaterials = (len(data) - n) // 40
        for i in range(nmaterials):
            out = s.unpack(data[n:n+40])
            (mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
            mat = MAT5.add_op2_data(out)
            self.add_thermal_material(mat, allow_overwrites=True)
            n += 40
        self.card_count['MAT5'] = nmaterials
        return n

    def _read_mat8(self, data, n):
        """
        MAT8(2503,25,288) - record 7
        """
        s = Struct(b(self._endian + 'i18f'))
        nmaterials = (len(data) - n) // 76
        for i in range(nmaterials):
            out = s.unpack(data[n:n+76])
            (mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
             TRef, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
            mat = MAT8.add_op2_data(out)
            self.add_op2_material(mat)
            n += 76
        self.card_count['MAT8'] = nmaterials
        return n

    def _read_mat9(self, data, n):
        """
        MAT9(2603,26,300) - record 9
        """
        ntotal = 140
        s = Struct(b(self._endian + 'i 30f iiii'))
        nmaterials = (len(data) - n) // ntotal
        for i in range(nmaterials):
            out = s.unpack(data[n:n+ntotal])
            assert len(out) == 35, out
            (mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
             rho, a1, a2, a3, a4, a5, a6, TRef, ge,
             blank1, blank2, blank3, blank4) = out
            assert blank1 == 0, blank1
            data_in = [mid, [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
                             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21],
                       rho, [a1, a2, a3, a4, a5, a6],
                       TRef, ge]
            mat = MAT9.add_op2_data(data_in)
            self.add_op2_material(mat)
            n += ntotal
        self.card_count['MAT9'] = nmaterials
        return n

    def _read_mat10(self, data, n):
        """
        MAT10(2801,28,365) - record 9
        """
        #self.log.debug("reading MAT10")
        ntotal = 20  # 5*4
        s = Struct(b(self._endian + 'i4f'))
        nmaterials = (len(data) - n) // ntotal
        assert nmaterials > 0, nmaterials
        for i in range(nmaterials):
            edata = data[n:n+20]
            out = s.unpack(edata)
            (mid, bulk, rho, c, ge) = out
            assert mid > 0, out
            mat = MAT10.add_op2_data(out)
            assert mat.mid > 0, mat
            self.add_op2_material(mat)
            n += 20
        self.card_count['MAT10'] = nmaterials
        return n

    def _read_mat11(self, data, n):
        """
        MAT11(2903,29,371)
        """
        ntotal = 128  # 23*4
        struc = Struct(b(self._endian + 'i 15f 16i'))
        nmaterials = (len(data) - n) // ntotal
        assert nmaterials > 0, nmaterials
        for i in range(nmaterials):
            edata = data[n:n+ntotal]

            out = struc.unpack(edata)
            (mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
             rho, a1, a2, a3, tref, ge) = out[:16]
            mat = MAT11.add_op2_data(out)
            self.add_op2_material(mat)
            n += ntotal
        self.card_count['MAT11'] = nmaterials
        return n

    def _read_mat11_old(self, data, n):
        """
        MAT11(2903,29,371)
        """
        ntotal = 80  # 20*4
        s = Struct(b(self._endian + 'i 15f 4s 4s 4s 4s'))
        nmaterials = (len(data) - n) // ntotal
        assert nmaterials > 0, nmaterials
        for i in range(nmaterials):
            edata = data[n:n+80]
            out = s.unpack(edata)
            (mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
             rho, a1, a2, a3, tref, ge,
             blank1, blank2, blank3, blank4) = out
            mat = MAT11.add_op2_data(out)
            assert mid > 0, mat
            self.add_op2_material(mat)
            n += 80
        self.card_count['MAT11'] = nmaterials
        return n

    def _read_mathp(self, data, n):
        """MATHP(4506,45,374) - Record 11"""
        nmaterials = 0
        s1 = Struct(b(self._endian + 'i7f3i23fi'))
        s2 = Struct(b(self._endian + '8i'))
        n2 = n
        while n2 < n:
            edata = data[n:n+140]
            n += 140
            out1 = s1.unpack(edata)
            (mid, a10, a01, d1, rho, alpha, tref, ge, sf, na, nd, kp,
             a20, a11, a02, d2,
             a30, a21, a12, a03, d3,
             a40, a31, a22, a13, a04, d4,
             a50, a41, a32, a23, a14, a05, d5,
             continue_flag) = out1
            data_in = [out1]

            if continue_flag:
                edata = data[n:n+32]  # 7*4
                n += 32
                out2 = s2.unpack(edata)
                (tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = out2
                data_in.append(out2)
                mat = MATHP.add_op2_data(data_in)
            self.add_op2_material(mat)
            nmaterials += 1
        self.card_count['MATHP'] = nmaterials
        return n

    def _read_mats1(self, data, n):
        """
        MATS1(503,5,90) - record 12
        """
        ntotal = 44  # 11*4
        s = Struct(b(self._endian + '3ifiiff3i'))
        nmaterials = (len(data) - n) // ntotal
        for i in range(nmaterials):
            edata = data[n:n+44]
            out = s.unpack(edata)
            (mid, tid, Type, h, yf, hr, limit1, limit2, a, bmat, c) = out
            data_in = [mid, tid, Type, h, yf, hr, limit1, limit2]
            mat = MATS1.add_op2_data(data_in)
            self.add_material_dependence(mat, allow_overwrites=True)
        self.card_count['MATS1'] = nmaterials
        return n

    def _read_matt1(self, data, n):
        self.log.debug('skipping MATT1 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT1 in MPT\n')
        return len(data)

    def _read_matt2(self, data, n):
        self.log.debug('skipping MATT2 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT2 in MPT\n')
        return len(data)

    def _read_matt3(self, data, n):
        self.log.debug('skipping MATT3 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT3 in MPT\n')
        return len(data)

    def _read_matt4(self, data, n):
        self.log.debug('skipping MATT4 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT4 in MPT\n')
        return len(data)

    def _read_matt5(self, data, n):
        self.log.debug('skipping MATT5 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT5 in MPT\n')
        return len(data)

# MATT8 - unused
    def _read_matt9(self, data, n):
        self.log.debug('skipping MATT9 in MPT\n')
        if self.is_debug_file:
            self.binary_debug.write('skipping MATT9 in MPT\n')
        return len(data)

# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO
# RADBND

    def _read_radm(self, data, n):
        """
        RADM(8802,88,413) - record 25
        .. todo:: add object
        """
        struct_i = self.struct_i
        nmaterials = 0
        ndata = len(data)
        while n < ndata:  # 1*4
            packs = []
            edata = data[n:n+4]
            number, = struct_i.unpack(edata)
            n += 4

            iformat = 'i %if' % (number)
            struct_i_nf = Struct(b(self._endian + iformat))
            #mid, absorb, emiss1, emiss2, ...
            ndata_per_pack = 1 + number
            nstr_per_pack = ndata_per_pack * 4

            nfields = (ndata - n) // 4
            npacks = nfields // ndata_per_pack
            for ipack in range(npacks):
                edata = data[n:n+nstr_per_pack]
                pack = list(struct_i_nf.unpack(edata))
                packs.append(pack)
                n += nstr_per_pack

                mat = RADM.add_op2_data(pack)
                self.add_thermal_BC(mat, mat.radmid)
                nmaterials += 1

        self.card_count['RADM'] = nmaterials
        return n

# RADMT

    def _read_nlparm(self, data, n):
        """
        NLPARM(3003,30,286) - record 27
        """
        #print("reading NLPARM")
        ntotal = 76  # 19*4
        s = Struct(b(self._endian + 'iif5i3f3iffiff'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+76]
            out = s.unpack(edata)
            #(sid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,maxR,rTolB) = out
            self.add_NLPARM(NLPARM.add_op2_data(out))
            n += ntotal
        self.card_count['NLPARM'] = nentries
        return n

    def _read_nlpci(self, data, n):
        if self.is_debug_file:
            self.binary_debug.write('skipping NLPCI in MPT\n')
        return len(data)

    def _read_tstepnl(self, data, n):
        """
        TSTEPNL(3103,31,337) - record 29
        """
        #print("reading TSTEPNL")
        ntotal = 88  # 19*4
        s = Struct(b(self._endian + 'iif5i3f3if3i4f'))
        nentries = (len(data) - n) // ntotal
        for i in range(nentries):
            edata = data[n:n+88]
            out = s.unpack(edata)
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out
            self.add_TSTEPNL(TSTEPNL.add_op2_data(out))
            n += ntotal
        self.card_count['TSTEPNL'] = nentries
        return n
