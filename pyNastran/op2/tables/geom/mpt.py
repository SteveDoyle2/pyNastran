"""defines readers for BDF objects in the OP2 MPT/MPTS table"""
#pylint: disable=C0111,C0103,C0301,W0612,R0914,R0201
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.materials import (CREEP, MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11, MATHP)
from pyNastran.bdf.cards.material_deps import (
    MATS1, MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, MATT11)
from pyNastran.bdf.cards.dynamic import (TSTEPNL,
                                         NLPARM, NLPARM_CONV_MAP, NLPARM_KMETHOD_MAP) # TSTEP
#from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 #PHBDY, CONV, CONVM, RADBC)
from pyNastran.bdf.cards.thermal.radiation import RADM
from pyNastran.op2.op2_interface.op2_reader import mapfmt # , reshape_bytes_block
from .geom2 import DoubleCardError
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom


class MPT:
    """defines methods for reading op2 materials & time-stepping methods"""
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

    def read_mpt_4(self, data: bytes, ndata: int):
        return self.op2._read_geom_4(self.mpt_map, data, ndata)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2

        #F:\work\pyNastran\examples\Dropbox\move_tpl\chkout01.op2
        self.mpt_map = {
            (1003, 10, 245) : ['CREEP', self.read_creep],  # record 1
            (103, 1, 77) : ['MAT1', self.read_mat1],       # record 3-msc-dmap2014
            (203, 2, 78) : ['MAT2', self.read_mat2],       # record 3
            (1403, 14, 122) : ['MAT3', self.read_mat3],    # record 4
            (2103, 21, 234) : ['MAT4', self.read_mat4],    # record 5
            (2203, 22, 235) : ['MAT5', self.read_mat5],    # record 6
            (2503, 25, 288) : ['MAT8', self.read_mat8],    # record 7
            (2603, 26, 300) : ['MAT9', self.read_mat9],    # record 8 - buggy
            (2801, 28, 365) : ['MAT10', self.read_mat10],  # record 9
            (2903, 29, 371) : ['MAT11', self.read_mat11],  # record ??? - NX specific - buggy?

            (4506, 45, 374) : ['MATHP', self.read_mathp],   # record 11
            (503, 5, 90) : ['MATS1', self.read_mats1],      # record 12
            (703, 7, 91) : ['MATT1', self.read_matt1],      # record 13 - not done
            (803, 8, 102) : ['MATT2', self.read_matt2],     # record 14
            #(1503, 14, 189) : ['MATT3', self.read_matt3],   # record 15 - not done
            (1503, 15, 189)  : ['MATT3', self.read_matt3],
            (2303, 23, 237) : ['MATT4', self.read_matt4],   # record 16 - not done
            (2403, 24, 238) : ['MATT5', self.read_matt5],   # record 17 - not done
            (2703, 27, 301) : ['MATT9', self.read_matt9],   # record 19 - not done
            (8802, 88, 413) : ['RADM', self.read_radm],     # record 25 - not done
            # record 26
            (3003, 30, 286) : ['NLPARM', self.read_nlparm],   # record 27
            (3104, 32, 350) : ['NLPCI', self.read_nlpci],     # record 28
            (3103, 31, 337) : ['TSTEPNL', self.read_tstepnl], # record 29
            (3303, 33, 988) : ['MATT11', self.read_matt11],

            (903, 9, 336) : ['MATT8', self.read_matt8],
            (8902, 89, 423) : ['RADMT', self.read_radmt],
            (9002, 90, 410) : ['RADBND', self.read_radbnd],
            (4801, 48, 961): ['MATPOR', self.read_matpor],
            (5101, 51, 642): ['MATDMG', self.read_fake],
            (14403, 144, 840): ['NLSTEP', self.read_fake],
            (4603, 46, 623): ['MATCRP', self.read_fake],
            (4701, 50, 965): ['MAT10C', self.read_fake],
            (3403, 34, 902): ['MATFT', self.read_fake],
            (2008, 20, 249): ['MATTC', self.read_fake],
            (4201, 42, 966): ['MATSR', self.read_fake],
            (8310, 83, 403): ['MATG', self.read_fake],

            (5303, 53, 906): ['MATCZ', self.read_fake],
            (10900, 109, 596): ['MATORT', self.read_matort],
            (10600, 106, 593): ['MATF', self.read_fake],
            (11000, 110, 597): ['MATTORT', self.read_fake],
            (13103, 131, 695): ['NLMOPTS', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            #(8310, 83, 403): ['???', self.read_fake],
            (15703, 157, 975): ['MAT1F', self.read_fake],
            (15903, 159, 977): ['MAT2F', self.read_fake],
            (16303, 163, 981): ['MAT2SP', self.read_fake],
            (12000, 120, 108): ['CONTRLT', self.read_fake],
            (15603, 156, 967): ['MONCARL', self.read_fake],
            (16803, 168, 1019): ['MAT8', self.read_fake],
            (16503, 165, 983): ['MAT9', self.read_mat9],
        }

    def add_op2_material(self, mat):
        #if mat.mid > 100000000:
            #raise RuntimeError('bad parsing...')
        self.op2._add_methods.add_structural_material_object(mat, allow_overwrites=False)
        #print(str(mat)[:-1])

    def read_matort(self, data: bytes, n: int) -> int:
        """
        MATORT(10900,109,596) - SOL 600 only

        Word Name Type Description
        1 MID      I   Material identification number
        2 E1      RS   modulus of elasticity in longitudinal direction
        3 E2      RS   modulus of elasticity in lateral direction
        4 E3      RS   modulus of elasticity in thickness direction
        5 NU12    RS   Poisson's ratio
        6 NU23    RS   Poisson's ratio
        7 NU31    RS   Poisson's ratio
        8 RHO     RS   Mass density
        9 G12     RS   shear modulus in plane 1-2
        10 G23    RS   shear modulus in plane 2-3
        11 G31    RS   shear modulus in plane 3-1
        12 A1     RS   coefficient of thermal expansion in direction 1
        13 A2     RS   coefficient of thermal expansion in direction 2
        14 A3     RS   coefficient of thermal expansion in direction 3
        15 TREF   RS   reference temperature
        16 GE     RS   structure damping coefficient
        17 IYLD    I   Yield criteria
        18 IHARD   I   Work harding rules
        19 SY     RS   Equivalent (von Mises) tensile yield stress
        20 SORNL  RS   For ORNL only, 10th cycle equivalent yield stress
        21 Y1     RS   Hill's yield stress in direction 1
        22 Y2     RS   Hill's yield stress in direction 2
        23 Y3     RS   Hill's yield stress in direction 3
        24 UNDEF  none
        25 YSHR1  RS   Hill's yield shear stress in direction 1
        26 YSHR2  RS   Hill's yield shear stress in direction 2
        27 YSHR3  RS   Hill's yield shear stress in direction 3
        """
        op2: OP2Geom = self.op2
        ntotal = 27 * self.size
        s = Struct(mapfmt(op2._endian + b'i 15f 2i 5f i 3f', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, e1, e2, e2, nu12, nu23, nu31, rho,
             g12, g23, g31, a1, a2, a3, tref, ge,
             iyield, ihard, sy, sornl, y1, y2, y3,
             undef, yshear1, yshear2, yshear3) = out
            assert undef == 0, undef
            #mat = op2.add_matort
            #self.add_op2_material(mat)
            n += ntotal
        op2.log.warning('geom skipping MATORT in MPT')
        #op2.card_count['MATORT'] = nmaterials
        return n

    def read_creep(self, data: bytes, n: int) -> int:
        """
        CREEP(1003,10,245) - record 1
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor
        nmaterials = (len(data) - n) // ntotal
        s = Struct(mapfmt(op2._endian + b'i2f4ifi7f', self.size))
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #(mid, T0, exp, form, tidkp, tidcp, tidcs, thresh,
             #Type, ag1, ag2, ag3, ag4, ag5, ag6, ag7) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  CREEP=%s\n' % str(out))
            mat = CREEP.add_op2_data(out)
            op2._add_methods.add_creep_material_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['CREEP'] = nmaterials
        return n

    def read_mat1(self, data: bytes, n: int) -> int:
        """
        MAT1(103,1,77) - record 2
        """
        op2: OP2Geom = self.op2
        ntotal = 48 * self.factor  # 12*4
        s = Struct(mapfmt(op2._endian + b'i10fi', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #(mid, E, G, nu, rho, A, tref, ge, St, Sc, Ss, mcsid) = out
            mat = MAT1.add_op2_data(out)
            self.add_op2_material(mat)
            n += ntotal
        op2.card_count['MAT1'] = nmaterials
        return n

    def read_mat2(self, data: bytes, n: int) -> int:
        """
        MAT2(203,2,78) - record 3

        ints    = (100000001, 10101010.0, 1010101.0, 0,   10101010.0, 0,   4545454.5, 0.05, 0.001, 0.001, 0, 0, 0, 0, 0, 0, 0,               -200000001, 0,   0,   0,   0,   0,   0,   0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        floats  = (100000001, 10101010.0, 1010101.0, 0.0, 10101010.0, 0.0, 4545454.5, 0.05, 0.001, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -200000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        """
        op2: OP2Geom = self.op2
        card_name = 'MAT2'
        card_obj = MAT2
        methods = {
            68 : self._read_mat2_68,
            92 : self._read_mat2_92,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, self.add_op2_material,
                methods, data, n)
        except DoubleCardError:
            raise
        return n

        #op2: OP2Geom = self.op2
        #ndatai = len(data) - n
        #if ndatai % 68 == 0:
        #    ntotal = 68  # 17*4
        #    s = Struct(op2._endian + b'i15fi')
        #else:
        #    ntotal = (17 + 6) * 4
        #    nleftover = ndatai % ntotal
        #    s = Struct(op2._endian + b'i15fi 6i')
        #    op2.log.warning(f'unexpected MAT2 format; ndatai={ndatai} ntotal={ntotal} nmaterials={ndatai // ntotal} '
        #                     f'leftover={ndatai % ntotal}')
        #    assert nleftover == 0, nleftover
        #nmaterials = ndatai // ntotal
        #
        #nbig_materials = 0
        #for unused_i in range(nmaterials):
        #    edata = data[n:n+ntotal]
        #    out = s.unpack(edata)
        #    if op2.is_debug_file:
        #        op2.binary_debug.write('  MAT2=%s\n' % str(out))
        #    if ntotal == 68:
        #        (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
        #         tref, ge, St, Sc, Ss, mcsid) = out
        #        mat = MAT2.add_op2_data(out)
        #    else:
        #        (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
        #         tref, ge, St, Sc, Ss, mcsid, *blanks) = out
        #        mat = MAT2.add_op2_data(out)
        #        op2.log.debug(f'\n{mat}')
        #    #print("MAT2 = ",out)
        #    if mid < 0:
        #        ndata = 4692
        #          #ints    = (100000001, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000002, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000003, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000004, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000005, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000006, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000007, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000008, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000009, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000010, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000011, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000012, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000013, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000014, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000015, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000016, 1260995102, 1246259552, 866102869, 1260995102, 904798606, 1248743307, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           100000017, 1258615469, 1244073420, 858401944, 1258615469, 900961598, 1245077109, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           200000001, 1267294939, 1241962828, 1238754656, 1254165792, 1238754656, 1244446583, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           200000002, 1267294939, 1241962828, 1238754656, 1254165792, 1238754656, 1244446583, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           200000003, 1267294939, 1241962828, 1238754656, 1254165792, 1238754656, 1244446583, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           200000004, 1267294939, 1241962828, 1238754656, 1254165792, 1238754656, 1244446583, 1030590824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        #          #           200000005, 1267294939, 1241962828, 1238754656, 1254165792, 1238754656, 1244446583, 1030590824, 0, 0, 0, 0, 0, 0)
        #        op2.show_data(data[12:], types='i', force=True)
        #        #op2.show_data(data[n:n+176], types='i')
        #        raise RuntimeError(mat)
        #        #, f'\n{mat}'
        #    if 0 < mid <= 1e8:  # just a checker for out of range materials
        #        self.add_op2_material(mat)
        #    else:
        #        nbig_materials += 1
        #        op2.big_materials[mid] = mat
        #    n += ntotal
        #
        #ncards = nmaterials - nbig_materials
        #if ncards:
        #    op2.card_count['MAT2'] = ncards
        #return n

    def _read_mat2_68(self, material: MAT2, data: bytes, n: int) -> tuple[int, list[MAT2]]:
        op2: OP2Geom = self.op2
        ntotal = 68 * self.factor  # 17*4
        s = Struct(op2._endian + mapfmt(b'i15fi', self.size))
        ndatai = len(data) - n
        assert ndatai % ntotal == 0
        nmaterials = ndatai // ntotal

        mats = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT2=%s\n' % str(out))

            #(mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
             #tref, ge, St, Sc, Ss, mcsid) = out
            mid = out[0]
            #print(mid)
            if mid < 0:
                n += ntotal
                continue
            assert mid > 0, mid
            mat = MAT2.add_op2_data(out)
            mats.append(mat)
            n += ntotal
        return n, mats


    def _read_mat2_92(self, material: MAT2, data: bytes, n: int) -> tuple[int, list[MAT2]]:
        """
        MAT2 MID   G11  G12  G13  G22  G23  G33  RHO
             A1    A2   A3   TREF GE   ST   SC   SS
             MCSID GE11 GE12 GE13 GE22 GE23 GE33

        MAT2    1       2.7866+38.3447+21.0139+23.3020+31.0139+21.1069+31.-9
                                                0.
                        .15     .15     .15     .15     .15     .15
        PDISTB   1      PCOMP   1        T       2      .2
                 1      1.
        PDISTB   2      PCOMP   1        T       5      .2
                 1      1.
        PCOMP    1                                              1.+5
                 1      .1      90.      YES     1      .8      45.      YES
                 1      .2      -45.     YES     1      .2      -45.     YES
                 1      .9      45.      YES     1      .1      90.      YES

        #ints    = (1, 2786.6, 834.47, 101.39, 3302.0, 101.39, 1106.9, 1e-9, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           100000001, 2953.90625, 832.65, -81.64, 3138.31, -81.64, 1105.08, 1e-9, 0, 0, 0, 0,         100000.0, 0, 0, 0, 0,         1041865114, 0.15, 0.15, 0.15, 0.15, 0.15, 200000001, 1161345214, 1146109282, -1024369403, 1162084660, -1024369403, 1149906036, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           400000001, 1069826929, -1472888832, 1064921892, -1077656719, 1064921892, -1491992576, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
        #floats  = (1, 2786.6, 834.47, 101.39, 3302.0, 101.39, 1106.9, 1e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           100000001, 2953.90625, 832.65, -81.64, 3138.31, -81.64, 1105.08, 1e-9, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 9.081061202086803e-32, 2955.54638671875, 832.9591064453125, -120.68167877197266, 3136.0751953125, -120.68167877197266, 1105.38916015625, 1e-9, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           400000001, 1.5333081483840942, -1.0075273948473296e-14, 0.974, -1.5333081483840942, 0.974, -2.0261570199409107e-15, 1e-9, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
        #
        #ints    = (1, 1160653210, 1146134036, 1120585646, 1162764288, 1120585646, 1149918413, 814313567,
        #             0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           100000001, 1161338496, 1146104342, -1029486047, 1162093848, -1029486047, 1149903566, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 200000001, 1161345214, 1146109282, -1024369403, 1162084660, -1024369403, 1149906036, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           400000001, 1069826929, -1472888832, 1064921892, -1077656719, 1064921892, -1491992576, 814313567,
        #             0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
        #floats  = (1, 2786.60009765625, 834.469970703125, 101.38999938964844, 3302.0, 101.38999938964844, 1106.9000244140625, 1e-9,
        #           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           100000001, 2953.90625, 832.6575927734375, -81.64478302001953, 3138.318359375, -81.64478302001953, 1105.087646484375, 1e-9, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 9.081061202086803e-32, 2955.54638671875, 832.9591064453125, -120.68167877197266, 3136.0751953125, -120.68167877197266, 1105.38916015625, 1e-9,
        #           0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
        #           400000001, 1.5333081483840942, -1.0075273948473296e-14, 0.974, -1.5333081483840942, 0.974, -2.0261570199409107e-15, 1e-9, 0.0, 0.0, 0.0, 0.0, 100000.0,
        #           0.0, 0.0, 0.0, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)

        ints    = (1, 1160653210, 1146134036, 1120585646, 1162764288, 1120585646, 1149918413, 814313567, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 100000001, 1161338496, 1146104342, -1029486047, 1162093848, -1029486047, 1149903566, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 200000001, 1161345214, 1146109282, -1024369403, 1162084660, -1024369403, 1149906036, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 400000001, 1069826929, -1472888832, 1064921892, -1077656719, 1064921892, -1491992576, 814313567, 0, 0, 0, 0, 1203982336, 0, 0, 0, 0, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114, 1041865114)
        floats  = (1.401298464324817e-45, 2786.60009765625, 834.469970703125, 101.38999938964844, 3302.0, 101.38999938964844, 1106.9000244140625, 9.999999717180685e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 2.3122342657588655e-35, 2953.90625, 832.6575927734375, -81.64478302001953, 3138.318359375, -81.64478302001953, 1105.087646484375, 9.999999717180685e-10, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 9.081061202086803e-32, 2955.54638671875, 832.9591064453125, -120.68167877197266, 3136.0751953125, -120.68167877197266, 1105.38916015625, 9.999999717180685e-10, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 1.3927371822189304e-24, 1.5333081483840942, -1.0075273948473296e-14, 0.9742910861968994, -1.5333081483840942, 0.9742910861968994, -2.0261570199409107e-15, 9.999999717180685e-10, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448, 0.15000000596046448)
        """
        op2: OP2Geom = self.op2
        ntotal = 92 * self.factor  # 23*4
        s = Struct(op2._endian + mapfmt(b'i15fi 6f', self.size))
        ndatai = len(data) - n
        assert ndatai % ntotal == 0
        nmaterials = ndatai // ntotal

        mats = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT2=%s\n' % str(out))

            #print(out)
            #(100000001, 2953.90625, 832.6575927734375, -81.64478302001953, 3138.318359375, -81.64478302001953, 1105.087646484375, 1e-9,
            #   0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0,
            #   0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
            #MAT2 MID   G11  G12  G13  G22  G23  G33  RHO
            #     A1    A2   A3   TREF GE   ST   SC   SS
            #     MCSID GE11 GE12 GE13 GE22 GE23 GE33
            #MAT2        10000001  2.9539E+03  8.3266E+02 -8.1645E+01  3.1383E+03 -8.1645E+01  1.1051E+03  1.0000E-09
            #+         0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.0000E+05  0.0000E+00  0.0000E+00  0.0000E+00
            #+                  0  1.5000E-01  1.5000E-01  1.5000E-01  1.5000E-01  1.5000E-01  1.5000E-01

            (mid, g1, g2, g3, g4, g5, g6, rho, aj1, aj2, aj3,
             tref, ge, St, Sc, Ss, mcsid,
             ge1, ge2, ge3, ge4, ge5, ge6) = out
            #ge_list = (ge1, ge2, ge3, ge4, ge5, ge6)
            n += ntotal
            if mid == -1:
                continue
            assert mid > 0, mid
            #assert ge == 0.0, f'mid={mid} ge={ge} ge_list={ge_list}'
            #assert max(blanks) == min(blanks) == 0, (mid, ge, blanks)
            mat = MAT2.add_op2_data(out)
            mats.append(mat)
        return n, mats

    def read_mat3(self, data: bytes, n: int) -> int:
        """
        MAT3(1403,14,122) - record 4
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor
        s = Struct(mapfmt(op2._endian + b'i8fi5fi', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            out = s.unpack(data[n:n+ntotal])
            (mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
             blank, ax, ath, az, tref, ge, blank) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT3=%s\n' % str(out))
            mat = MAT3.add_op2_data([mid, ex, eth, ez, nuxth, nuthz,
                                     nuzx, rho, gzx, ax, ath, az, tref, ge])
            self.add_op2_material(mat)
            n += ntotal
        op2.card_count['MAT3'] = nmaterials
        return n

    def read_mat4(self, data: bytes, n: int) -> int:
        """
        MAT4(2103,21,234) - record 5
        """
        op2: OP2Geom = self.op2
        ntotal = 44 * self.factor
        s = Struct(mapfmt(op2._endian + b'i10f', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            out = s.unpack(data[n:n+ntotal])
            #(mid, k, cp, rho, h, mu, hgen, refenth, tch, tdelta, qlat) = out
            mat = MAT4.add_op2_data(out)
            op2._add_methods.add_thermal_material_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['MAT4'] = nmaterials
        return n

    def read_mat5(self, data: bytes, n: int) -> int:
        """
        MAT5(2203,22,235) - record 6
        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'i9f')
        nmaterials = (len(data) - n) // 40
        for unused_i in range(nmaterials):
            out = s.unpack(data[n:n+40])
            #(mid, k1, k2, k3, k4, k5, k6, cp, rho, hgen) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT5=%s\n' % str(out))
            mat = MAT5.add_op2_data(out)
            op2._add_methods.add_thermal_material_object(mat, allow_overwrites=False)
            n += 40
        op2.card_count['MAT5'] = nmaterials
        return n

    def read_mat8(self, data: bytes, n: int) -> int:
        """
        MAT8(2503,25,288) - record 7
        """
        op2: OP2Geom = self.op2
        ntotal = 76 * self.factor
        s = Struct(mapfmt(op2._endian + b'i18f', self.size))
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert ndatai % ntotal == 0
        for unused_i in range(nentries):
            out = s.unpack(data[n:n+ntotal])
            #(mid, E1, E2, nu12, G12, G1z, G2z, rho, a1, a2,
            # tref, Xt, Xc, Yt, Yc, S, ge, f12, strn) = out
            mat = MAT8.add_op2_data(out)
            self.add_op2_material(mat)
            n += ntotal
        op2.card_count['MAT8'] = nentries
        return n

    def read_mat9(self, data: bytes, n: int) -> int:
        """
        MAT9(2603,26,300) - record 9
        """
        card_name = 'MATT9'
        card_obj = MATT9
        methods = {
            #140 : self._read_mat9_140,
            140 : self._read_mat9_140,
            224 : self._read_mat9_224,
        }
        #try:
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_double_card(
            card_name, card_obj,
            self.add_op2_material,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n

    def _read_mat9_140(self, card_obj, data: bytes, n: int) -> tuple[int, list[MAT9]]:
        op2: OP2Geom = self.op2
        #op2.log.info('geom skipping MAT9')
        #return len(data)
        materials = []
        ndatai = len(data) - n
        size = self.size
        s2 = Struct(mapfmt(op2._endian + b'i 30f iifi', size))
        ntotal = 35 * size
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} nmaterials={nmaterials} leftover={ndatai % ntotal}'
        if op2.is_debug_file:
            op2.binary_debug.write(
                '  MAT9=(mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, '
                'g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21, '
                'rho, a1, a2, a3, a4, a5, a6, tref, ge, '
                'blank1, blank2, blank3, blank4)\n')
        for unused_i in range(nmaterials):
            out = s2.unpack(data[n:n+ntotal])
            if op2.is_debug_file:
                op2.binary_debug.write('    MAT9=%s\n' % str(out))
            #print(out)
            (mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
             rho, a1, a2, a3, a4, a5, a6, tref, ge,
             blank1, blank2, blank3, blank4) = out
            assert mid > 0, mid
            assert blank1 == 0, blank1
            assert blank2 == 0, blank2
            assert blank3 == 0, blank3
            assert blank4 == 0, blank4
            g = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
                 g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21]
            a = [a1, a2, a3, a4, a5, a6]
            data_in = [mid, g, rho, a, tref, ge]
            #print(f'mid={mid} g={g} rho={rho} a={a} tref={tref} ge={ge}')

            mat = MAT9.add_op2_data(data_in)
            #print(mat)
            materials.append(mat)
            n += ntotal
        return n, materials

    def _read_mat9_224(self, card_obj, data: bytes, n: int) -> tuple[int, list[MAT9]]:
        op2: OP2Geom = self.op2
        #op2.log.info('geom skipping MAT9')
        #return len(data)
        materials = []
        ndatai = len(data) - n
        size = self.size
        ntotal = (35 + 21) * size # 56*4
        s2 = Struct(mapfmt(op2._endian + b'i 30f iiii 21f', size))
        #s2 = Struct(op2._endian + b'i 30f iiii 3f 3i 2f 3i f 3i f 2i f if')
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0, f'ndatai={ndatai} ntotal={ntotal} nmaterials={nmaterials} leftover={ndatai % ntotal}'
        # $       2       3       4       5       6       7       8       9       0
        # $       ID      G11     G12     G13     G14     G15     G16     G22
        # MAT9    2500003 1.0769+74.6154+64.6154+60.      0.      0.      1.0769+7
        # $       G23     G24     G25     G26     G33     G34     G35     G36
        #         4.6154+60.      0.      0.      1.0769+70.      0.      0.
        # $       G44     G45     G46     G55     G56     G66     RHO     A1
        #         3.0769+60.      0.      3.0769+60.      3.0769+63.7E-8
        # $       A2      A3      A4      A5      A6      T       GE
        #                                                         .013
        # $       GE11    GE12    GE13    GE14    GE15    GE16    GE22    GE23
        #         .017    .046    .046                            .017    .046
        # $       GE24    GE25    GE26    GE33    GE34    GE35    GE36    GE44
        #         .017                            .03
        # $       GE45    GE46    GE55    GE56    GE66
        #                         .03             .03
        if op2.is_debug_file:
            op2.binary_debug.write(
                '  MAT9=(mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, '
                'g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21, '
                'rho, a1, a2, a3, a4, a5, a6, tref, ge, '
                'ge1, e2 ge3, ge4, ge5, g6, ge7, ge8, ge9, ge10, ge11, ge12, '
                'ge13, ge14, ge15, ge16, ge17, ge18, ge19, ge20, ge21)\n')
        for unused_i in range(nmaterials):
            out = s2.unpack(data[n:n+ntotal])
            if op2.is_debug_file:
                op2.binary_debug.write('    MAT9=%s\n' % str(out))
            (mid, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21,
             rho, a1, a2, a3, a4, a5, a6, tref, ge,
             blank1, blank2, blank3, blank4, *ge_list) = out
            assert mid > 0, mid
            assert blank1 == 0, blank1
            assert blank2 == 0, blank2
            assert blank3 == 0, blank3
            assert blank4 == 0, blank4
            data_in = [mid, [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10,
                             g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, g21],
                       rho, [a1, a2, a3, a4, a5, a6],
                       tref, ge]
            mat = MAT9.add_op2_data(data_in)
            materials.append(mat)
            n += ntotal
        return n, materials

    def read_mat10(self, data: bytes, n: int) -> int:
        """
        MAT10(2801,28,365) - record 9

        Word Name Type Description
        1 MID   I Material identification number
        2 BULK RS Bulk modulus
        3 RHO  RS Mass density
        4 C    RS Speed of sound
        5 GE   RS Structural damping coefficient

        """
        op2: OP2Geom = self.op2
        card_name = 'MAT10'
        card_obj = MAT10
        methods = {
            20 : self._read_mat10_20,
            24 : self._read_mat10_24,
            #56 : self._read_ctria6_current_56,
        }
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, self.add_op2_material,
                methods, data, n)
        except DoubleCardError:
            raise

        #n = self._read_split_card(data, n,
                                  #self._read_ctria6_current, self._read_ctria6_v2001,
                                  #'CTRIA6', CTRIA6, self.add_op2_element)
        return n

    def _read_mat10_20(self, material: MAT10, data: bytes, n: int) -> tuple[int, list[MAT10]]:
        op2: OP2Geom = self.op2
        ntotal = 20 * self.factor # 5*4
        s = Struct(mapfmt(op2._endian + b'i4f', self.size))
        ndatai = (len(data) - n)
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nmaterials > 0, nmaterials
        materials = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            n += ntotal

            (mid, bulk, rho, c, ge) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT10=%s\n' % str(out))
            if mid == 0 and bulk == 0. and rho == 0. and c == 0. and ge == 0.:
                op2.log.debug('  skipping empty MAT10...')
                continue
            mat = MAT10.add_op2_data(out)
            assert mat.mid > 0, mat
            materials.append(mat)
        return n, materials

    def _read_mat10_24(self, material: MAT10, data: bytes, n: int) -> tuple[int, list[MAT10]]:
        """
        1 MID   I  Material identification number
        2 BULK  RS Bulk modulus
        3 RHO   RS Mass density
        4 C     RS Speed of sound
        5 GE    RS Structural damping coefficient
        6 ALPHA RS
        data = (25, 1.0, 0.1, 3.16, 0.02, 0)
             = (25, 1.0, 0.1, 3.16, 0.02, 0.0)
        """
        op2: OP2Geom = self.op2
        ntotal = 24 * self.factor # 6*4
        s = Struct(mapfmt(op2._endian + b'i5f', self.size))
        ndatai = (len(data) - n)
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nmaterials > 0, nmaterials
        materials = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            n += ntotal

            (mid, bulk, rho, c, ge, alpha) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT10=%s\n' % str(out))
            if mid == 0 and bulk == 0. and rho == 0. and c == 0. and ge == 0. and alpha == 0.0:
                op2.log.debug('  skipping empty MAT10...')
                continue
            mat = MAT10.add_op2_data(out)
            assert mat.mid > 0, mat
            materials.append(mat)
        return n, materials

    def read_mat11(self, data: bytes, n: int) -> int:
        """
        MAT11(2903,29,371)
        """
        op2: OP2Geom = self.op2
        ntotal = 128 * self.factor  # 23*4
        struc = Struct(mapfmt(op2._endian + b'i 15f 16i', self.size))
        nmaterials = (len(data) - n) // ntotal
        assert nmaterials > 0, nmaterials
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]

            out = struc.unpack(edata)
            #(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
             #rho, a1, a2, a3, tref, ge) = out[:16]
            if op2.is_debug_file:
                op2.binary_debug.write('  MA11=%s\n' % str(out))
            mat = MAT11.add_op2_data(out)
            self.add_op2_material(mat)
            n += ntotal
        op2.card_count['MAT11'] = nmaterials
        return n

    def _read_mat11_old(self, data: bytes, n: int) -> int:
        """
        MAT11(2903,29,371)
        """
        op2: OP2Geom = self.op2
        ntotal = 80  # 20*4
        s = Struct(op2._endian + b'i 15f 4s 4s 4s 4s')
        nmaterials = (len(data) - n) // ntotal
        assert nmaterials > 0, nmaterials
        for unused_i in range(nmaterials):
            edata = data[n:n+80]
            out = s.unpack(edata)
            (mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
             rho, a1, a2, a3, tref, ge,
             blank1, blank2, blank3, blank4) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MAT11-old=%s\n' % str(out))
            mat = MAT11.add_op2_data(out)
            assert mid > 0, mat
            self.add_op2_material(mat)
            n += 80
        op2.card_count['MAT11'] = nmaterials
        return n

    def read_mathp(self, data: bytes, n: int) -> int:
        """
        MATHP(4506,45,374) - Record 11

        NX/MSC
        1 MID       I Material identification number
        2 A10      RS Material constant related to distortional deformation
        3 A01      RS Material constant related to distortional deformation
        4 D1       RS Material constant related to volumetric deformation
        5 RHO      RS Mass density
        6 ALPHA    RS Coefficient of volumetric thermal expansion
        7 TREF     RS Reference temperature
        8 GE       RS Structural damping element coefficient
        9 SF        I ???
        10 NA       I Order of the distortional strain energy polynomial function
        11 ND       I Order of the volumetric strain energy polynomial function
        12 KP      RS ???
        13 A20     RS Material constant related to distortional deformation
        14 A11     RS Material constant related to distortional deformation
        15 A02     RS Material constant related to distortional deformation
        16 D2      RS Material constant related to volumetric deformation
        17 A30     RS Material constant related to distortional deformation
        18 A21     RS Material constant related to distortional deformation
        19 A12     RS Material constant related to distortional deformation
        20 A03     RS Material constant related to distortional deformation
        21 D3      RS Material constant related to volumetric deformation
        22 A40     RS Material constant related to distortional deformation
        23 A31     RS Material constant related to distortional deformation
        24 A22     RS Material constant related to distortional deformation
        25 A13     RS Material constant related to distortional deformation
        26 A04     RS Material constant related to distortional deformation
        27 D4      RS Material constant related to volumetric deformation
        28 A50     RS Material constant related to distortional deformation
        29 A41     RS Material constant related to distortional deformation
        30 A32     RS Material constant related to distortional deformation
        31 A23     RS Material constant related to distortional deformation
        32 A14     RS Material constant related to distortional deformation
        33 A05     RS Material constant related to distortional deformation
        34 D5      RS Material constant related to volumetric deformation
        35 CONTFLG  I Continuation flag
        CONTFLG =1 With continuation
        36 TAB1 I TABLES1 identification number which defines tension/compression
        37 TAB2 I TABLES1 identification number which defines equibiaxial tension
        38 TAB3 I TABLES1 identification number which defines simple shear
        39 TAB4 I TABLES1 identification number which defines pure shear
        40 UNDEF(3) None
        43 TAB5 I TABLES1 identification number which defines volumetric compression
        CONTFLG =0 Without continuation
        End CONTFLG
        """
        op2: OP2Geom = self.op2
        nmaterials = 0
        s1 = Struct(mapfmt(op2._endian + b'i7f3i23fi', self.size))
        s2 = Struct(mapfmt(op2._endian + b'8i', self.size))
        n2 = len(data)
        ntotal1 = 140 * self.factor
        ntotal2 = 32 * self.factor  # 7*4
        while n < n2:
            edata = data[n:n+ntotal1]
            out1 = s1.unpack(edata)
            n += ntotal1
            (mid, a10, a01, d1, rho, alpha, tref, ge, sf, na, nd, kp,
             a20, a11, a02, d2,
             a30, a21, a12, a03, d3,
             a40, a31, a22, a13, a04, d4,
             a50, a41, a32, a23, a14, a05, d5,
             continue_flag) = out1

            if n == n2:
                # we have to hack the continue_flag because it's wrong...
                # C:\Users\sdoyle\Dropbox\move_tpl\ehq45.op2
                continue_flag = 0
                out1 = (mid, a10, a01, d1, rho, alpha, tref, ge, sf, na, nd, kp,
                        a20, a11, a02, d2,
                        a30, a21, a12, a03, d3,
                        a40, a31, a22, a13, a04, d4,
                        a50, a41, a32, a23, a14, a05, d5,
                        continue_flag)
            data_in = [out1]

            if continue_flag:
                edata = data[n:n+ntotal2]
                out2 = s2.unpack(edata)
                n += ntotal2
                #(tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = out2
                data_in.append(out2)
            mat = MATHP.add_op2_data(data_in)

            if op2.is_debug_file:
                op2.binary_debug.write('  MATHP=%s\n' % str(out1))
            op2._add_methods.add_hyperelastic_material_object(mat)
            nmaterials += 1
        assert nmaterials > 0, 'MATP nmaterials=%s' % nmaterials
        op2.card_count['MATHP'] = nmaterials
        return n

    def read_mats1(self, data: bytes, n: int) -> int:
        """
        MATS1(503,5,90) - record 12
        """
        op2: OP2Geom = self.op2
        ntotal = 48 * self.factor  # 12*4
        s = Struct(mapfmt(op2._endian + b'3ifiiff4i', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, tid, Type, h, yf, hr, limit1, limit2, strmeas, a, bmat, c) = out
            assert a == 0, a
            assert bmat == 0, bmat
            assert c == 0, c
            data_in = [mid, tid, Type, h, yf, hr, limit1, limit2, strmeas]
            if op2.is_debug_file:
                op2.binary_debug.write('  MATS1=%s\n' % str(out))
            mat = MATS1.add_op2_data(data_in)
            op2._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['MATS1'] = nmaterials
        return n

    def read_matt1(self, data: bytes, n: int) -> int:
        """
        MATT1(703,7,91)
        checked NX-10.1, MSC-2016
        """
        op2: OP2Geom = self.op2
        s = Struct(mapfmt(op2._endian + b'12i', self.size))
        ntotal = 48 *  self.factor # 12*4
        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT1=%s\n' % str(out))
            #(mid, tableid, ...., None) = out
            mat = MATT1.add_op2_data(out)
            op2._add_methods.add_material_dependence_object(mat)
            n += ntotal
        op2.increase_card_count('MATT1', ncards)
        return n

    def read_matt2(self, data: bytes, n: int) -> int:
        card_name = 'MATT2'
        card_obj = MATT2
        methods = {
            68 : self._read_matt2_68,
            92 : self._read_matt2_92,
        }
        #try:
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_double_card(
            card_name, card_obj,
            op2._add_methods.add_material_dependence_object,
            methods, data, n)
        #except DoubleCardError:
            #raise
        return n

    def _read_matt2_68(self, card_obj, data: bytes, n: int) -> tuple[int, list[MATT2]]:
        """
        1 MID         I Material identification number
        2 TID(15)     I TABLEMi entry identification numbers
        17        UNDEF none Not used
        """
        op2: OP2Geom = self.op2
        ntotal = 68 * self.factor # 17*4
        s = Struct(mapfmt(op2._endian + b'17i', self.size))
        nmaterials = (len(data) - n) // ntotal
        cards = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, g11_table, g12_table, g13_table, g22_table,
             g23_table, g33_table, rho_table,
             a1_table, a2_table, a3_table, unused_zeroa, ge_table,
             st_table, sc_table, ss_table, unused_zerob) = out
            assert mid > 0, mid
            assert unused_zeroa == 0, f'unused_zeroa={unused_zeroa} out={out}'
            assert unused_zerob == 0, f'unused_zerob={unused_zerob} out={out}'
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT2=%s\n' % str(out))

            mat = MATT2(mid, g11_table, g12_table, g13_table, g22_table,
                        g23_table, g33_table, rho_table,
                        a1_table, a2_table, a3_table, ge_table,
                        st_table, sc_table, ss_table, comment='')
            cards.append(mat)
            #op2._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        return n, cards

    def _read_matt2_92(self, card_obj, data: bytes, n: int) -> tuple[int, list[MATT2]]:
        """
        1 MID         I Material identification number
        2 TID(15)     I TABLEMi entry identification numbers
        17        UNDEF none Not used
        """
        op2: OP2Geom = self.op2
        ntotal = 92 * self.factor # 23*4
        s = Struct(mapfmt(op2._endian + b'17i 6i', self.size))
        nmaterials = (len(data) - n) // ntotal
        cards = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, g11_table, g12_table, g13_table, g22_table,
             g23_table, g33_table, rho_table,
             a1_table, a2_table, a3_table, unused_zeroa, ge_table,
             st_table, sc_table, ss_table, unused_zerob, *unused_zeroc) = out
            assert mid > 0, mid
            assert unused_zeroa == 0, f'unused_zeroa={unused_zeroa} out={out}'
            assert unused_zerob == 0, f'unused_zerob={unused_zerob} out={out}'
            assert max(unused_zeroc) == 0
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT2=%s\n' % str(out))

            mat = MATT2(mid, g11_table, g12_table, g13_table, g22_table,
                        g23_table, g33_table, rho_table,
                        a1_table, a2_table, a3_table, ge_table,
                        st_table, sc_table, ss_table, comment='')
            cards.append(mat)
            #op2._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        return n, cards

    def read_matt3(self, data: bytes, n: int) -> int:
        r"""
        Word Name Type Description
        1 MID     I Material identification number
        2 TID(15) I entry identification numbers

        test_op2 -g C:\MSC.Software\msc_nastran_runs\varmat4c.op2
        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\m402mat3_matt3_ex_ey_nuxy_gxy.op2
        """
        op2: OP2Geom = self.op2
        ntotal = 64 * self.factor # 16*4
        s = Struct(mapfmt(op2._endian + b'16i', self.size))
        ndatai = len(data) - n
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0
        assert nmaterials > 0
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)

            (mid, *tables, a, b, c, d) = out
            #mid, _ex, _eth, _ez, _nuxth, _nuthz, _nuzx, rho, _gxz, ax, _ath, az, _ge, *other = out
            if op2._nastran_format == 'msc':
                mid, a, b, c, d, e, f, rho, g, h, i, ax, j, az, k, m = out
                assert sum([a, b, c, d, e, f, g, h, i, j, k, m]) == 0, out
                #assert rho == 92, rho
                #assert ax == 93, out
                #assert az == 93, az
                #assert rho == 92, rho
            elif op2._nastran_format == 'nx':
                # $ NX
                # $ MID	T(EX)	T(EY)	T(EZ)	T(NUXY)	T(NUYZ)	T(NUZX)	T(RHO)
                # $ T(GXY)	T(GZX)	T(AX)	T(AY)	T(AZ)	T(GE)
                # MATT3          2       1               2       3                        +
                # +              4
                # $            mid      ex              ez    nuxy
                # MATT3          2       1               2       3                        +
                # $            gxy
                # +              4
                #$            mid      ex      ey      ez    nuxy    nuyz    nuzx      rho
                #$            gxy     gzx     ax      ay      az       ge

                mid, ex, ey, ez, nuxy, nuyz, nuzx, rho, gxy, gzx, ax, ay, ax, gea, geb, gec = out
                #assert ex == 1, f'ex={ex} out={out}'
                #assert ey == 0, f'ey={ey} out={out}'
                #assert ez == 2, out
                #assert nuxy == 3, out
                #assert gxy == 4, out
                #assert ge == 0, ge
                assert sum([nuyz, nuzx, rho, gzx, ax, ay, ax, gea, geb, gec]) == 0, [mid, nuyz, nuzx, rho, gzx, ax, ay, ax, gea, geb, gec]
            else:  # pragma: no cover
                raise NotImplementedError(op2._nastran_format)
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT3=%s\n' % str(out))
            #mat = MATT3(mid, ex_table=None, eth_table=None, ez_table=None, nuth_table=None,
                     #nuxz_table=None, rho_table=None, gzx_table=None,
                     #ax_table=None, ath_table=None, az_table=None, ge_table=None,)
            mat = MATT3(mid, *tables, comment='')
            op2._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['MATT3'] = nmaterials
        return n

    def read_matt4(self, data: bytes, n: int) -> int:
        """
        MATT4(2303,23,237)
        checked NX-10.1, MSC-2016
        """
        op2: OP2Geom = self.op2
        struct_7i = Struct(mapfmt(op2._endian + b'7i', self.size))
        ntotal = 28 * self.factor # 7*4
        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = struct_7i.unpack(edata)

            if op2.is_debug_file:
                op2.binary_debug.write('  MATT4=%s\n' % str(out))
            #(mid, tk, tcp, null, th, tmu, thgen) = out
            mat = MATT4.add_op2_data(out)
            op2._add_methods.add_material_dependence_object(mat)
            n += ntotal
        op2.increase_card_count('MATT4', ncards)
        return n

    def read_matt5(self, data: bytes, n: int) -> int:
        """
        MATT5(2403,24,238)
        checked NX-10.1, MSC-2016

        """
        op2: OP2Geom = self.op2
        s = Struct(op2._endian + b'10i')
        ntotal = 40 # 10*4
        ncards = (len(data) - n) // ntotal
        for unused_i in range(ncards):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT4=%s\n' % str(out))
            #(mid, kxx_table, kxy_table, kxz_table, kyy_table, kyz_table, kzz_table,
            # cp_table, null, hgen_table) = out
            mat = MATT5.add_op2_data(out)
            op2._add_methods.add_material_dependence_object(mat)
            n += ntotal
        op2.increase_card_count('MATT5', ncards)
        return n

    def read_matt8(self, data: bytes, n: int) -> int:
        """common method to read MSC/NX MATT8s"""
        op2: OP2Geom = self.op2
        n = op2.reader_geom2._read_dual_card(
            data, n, self._read_matt8_18, self._read_matt8_19,
            'MATT8', op2._add_methods.add_material_dependence_object)
        return n

    def _read_matt8_19(self, data: bytes, n: int) -> int:
        """
        MATT8 (903, 9, 336)
        (903, 9, 336,
        2, 1, 2, 0, 3, 4, 5, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, xxx)

        Word Name Type Description
        1 MID I
        2 TID(9)  I TABLEMi entry identification numbers
        11 UNDEF None
        12 TID(7) I TABLEMi entry identification numbers
        19 UNDEF None
        """
        op2: OP2Geom = self.op2
        ntotal = 76 * self.factor  # 35*4
        s = Struct(mapfmt(op2._endian + b'i18i', self.size))
        ndatai = len(data) - n
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0

        matt8s = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, e1_table, e2_table, nu12_table, g12_table, g1z_table, g2z_table, trho,
             ta1, ta2, blank,
             xt_table, xc_table, yt_table, yc_table,
             s_table, ge_table, f12_table, final) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT8=%s\n' % str(out))
            mat = MATT8(mid,
                        e1_table=e1_table, e2_table=e2_table,
                        nu12_table=nu12_table, g12_table=g12_table,
                        g1z_table=g1z_table, g2z_table=g2z_table, rho_table=trho,
                        a1_table=ta1, a2_table=ta2,
                        xt_table=xt_table, xc_table=xc_table,
                        yt_table=yt_table, yc_table=yc_table,
                        s_table=s_table, ge_table=ge_table, f12_table=f12_table)
            assert blank == 0, f'blank={blank} out={out}'
            assert final == 0, f'final={final} out={out}'
            #assert xc_table == 0, f'xc_table={xc_table} out={out}'
            #assert yt_table == 0, f'yt_table={yt_table} out={out}'
            #assert yc_table == 0, f'yc_table={yc_table} out={out}'
            #assert s_table == 0, f's_table={s_table} out={out}'
            #assert ge_table == 0, f'ge_table={ge_table} out={out}'
            #assert f12_table == 0, f'f12_table={f12_table} out={out}'
            str(mat)
            matt8s.append(mat)
            n += ntotal

        return n, matt8s

    def _read_matt8_18(self, data: bytes, n: int) -> int:
        """
        MATT8 (903, 9, 336)
        (903, 9, 336,
        2, 1, 2, 0, 3, 4, 5, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0)

        Word Name Type Description
        1 MID I
        2 TID(9)  I TABLEMi entry identification numbers
        11 UNDEF None
        12 TID(7) I TABLEMi entry identification numbers
        19 UNDEF None
        """
        op2: OP2Geom = self.op2
        ntotal = 72 * self.factor  # 35*4
        s = Struct(mapfmt(op2._endian + b'18i', self.size))
        ndatai = len(data) - n
        nmaterials = ndatai // ntotal
        assert ndatai % ntotal == 0
        matt8s = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, e1_table, e2_table, nu12_table, g12_table, g1z_table, g2z_table, trho,
             ta1, ta2, blank,
             xt_table, xc_table, yt_table, yc_table,
             s_table, ge_table, f12_table) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT8=%s\n' % str(out))
            mat = MATT8(mid,
                        e1_table=e1_table, e2_table=e2_table,
                        nu12_table=nu12_table, g12_table=g12_table,
                        g1z_table=g1z_table, g2z_table=g2z_table, rho_table=trho,
                        a1_table=ta1, a2_table=ta2,
                        xt_table=xt_table, xc_table=xc_table,
                        yt_table=yt_table, yc_table=yc_table,
                        s_table=s_table, ge_table=ge_table, f12_table=f12_table)
            assert blank == 0, f'blank={blank} out={out}'
            #assert xc_table == 0, f'xc_table={xc_table} out={out}'
            #assert yt_table == 0, f'yt_table={yt_table} out={out}'
            #assert yc_table == 0, f'yc_table={yc_table} out={out}'
            #assert s_table == 0, f's_table={s_table} out={out}'
            #assert ge_table == 0, f'ge_table={ge_table} out={out}'
            #assert f12_table == 0, f'f12_table={f12_table} out={out}'
            str(mat)
            matt8s.append(mat)
            n += ntotal

        return n, matt8s

    def read_matt9(self, data: bytes, n: int) -> int:
        """common method for reading MATT9s"""
        op2: OP2Geom = self.op2
        card_name = 'MATT9'
        card_obj = MATT9
        methods = {
            140 : self._read_matt9_140,
            224 : self._read_matt9_224,
        }
        add_method = op2._add_methods.add_material_dependence_object
        #self._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
        try:
            n = op2.reader_geom2._read_double_card(
                card_name, card_obj, add_method,
                methods, data, n)
        except DoubleCardError:
            raise
            #op2.log.warning(f'try-except {card_name}')
            #n = self._read_split_card(data, n,
                                      #self._read_cquad8_current, self._read_cquad8_v2001,
                                      #card_name, self.add_op2_element)

        return n

    def _read_matt9_224(self, card_obj, data: bytes, n: int) -> tuple[int, list[MATT9]]:
        r"""
        Word Name Type Description
        1 MID    I Material identification number
        2 TC(21) I TABLEMi identification numbers for material property matrix
        23 TRHO  I TABLEMi identification number for mass density
        24 TA(6) I TABLEMi identification numbers for thermal expansion coefficients
        30 UNDEF None
        31 TGE   I TABLEMi identification number for structural damping coefficient
        32 UNDEF(4) None
        ????

        # C:\MSC.Software\msc_nastran_runs\freefld.op2
        MATT9,1101,2 ,3 ,4 ,,,,8 ,+P101
        +P101,9 ,,,,13
        $ mid, g11, g12, g13, g14, g15, g16, g22
        MATT9,1251,2 ,3 ,4 ,5 ,6 ,7 ,8 ,+P251
        $ g23, g24, g25, g26, g33, g34, g35, g36
        +P251,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,+P252
        $     44  45 46, 55  56 66
        +P252,17 ,18 ,,  20 , , 22
           g36,

        ints    = (
            2703, 27, 301,
            mid   11 12  13 14 15 16 22 23 24 25 26 33  ?
            1101, 2,  3, 4, 0, 0, 0, 8, 9, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1102, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1103, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1151, 2,  3, 4, 0, 0, 0, 8, 9, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1152, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1153, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1201, 2,  3, 4, 0, 0, 0, 8, 9, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1203, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1204, 2,  3, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...)
            ...
           (1251, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 0, 20, 0, 22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),

        """
        op2: OP2Geom = self.op2
        #self.show_data(data, types='if')
        ntotal = 224 * self.factor  # 56*4
        s = Struct(mapfmt(op2._endian + b'56i', self.size))
        nmaterials = (len(data) - n) // ntotal
        materials = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print(out)
            (mid, g11, g12, g13, g14, g15, g16, g22,
             g23, g24, g25, g26, g33, g34, g35, g36,
             g44, g45, g46, g55, g56, g66, *other) = out

            tc_tables = (g11, g12, g13, g14, g15, g16, g22,
                         g23, g24, g25, g26, g33, g34, g35, g36,
                         g44, g45, g46, g55, g56, g66)
            #(mid, tc_tables, *other) = out
            #print(mid, tc_tables, *other)
            if sum(other) != 0:
                op2.log.warning(f'mATT9 mid={mid} other={other} flags are dropped...')
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT9=%s\n' % str(out))

            #MATT9(mid, g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                                 #g15_table=None, g16_table=None, g22_table=None, g23_table=None,
                                 #g24_table=None, g25_table=None, g26_table=None, g33_table=None,
                                 #g34_table=None, g35_table=None, g36_table=None, g44_table=None,
                                 #g45_table=None, g46_table=None, g55_table=None, g56_table=None,
                                 #g66_table=None, rho_table=None,
                                 #a1_table=None, a2_table=None, a3_table=None,
                                 #a4_table=None, a5_table=None, a6_table=None, ge_table=None, comment='')
            tc_tables = [g11, g12, g13, g14, g15, g16, g22, g23, g24, g25, g26, g33]
            assert mid > 0, (mid, tc_tables, *other)
            trho = 0
            ta1 = 0
            ta2 = 0
            ta3 = 0
            ta4 = 0
            ta5 = 0
            ta6 = 0
            tge = 0
            mat = MATT9(mid, *tc_tables, trho, ta1, ta2, ta3, ta4, ta5, ta6, tge, comment='')
            materials.append(mat)
            n += ntotal
        return n, materials

    def _read_matt9_140(self, card_obj, data: bytes, n: int) -> tuple[int, list[MATT9]]:
        """
        Word Name Type Description
        1 MID    I Material identification number
        2 TC(21) I TABLEMi identification numbers for material property matrix
        23 TRHO  I TABLEMi identification number for mass density
        24 TA(6) I TABLEMi identification numbers for thermal expansion coefficients
        30 UNDEF None
        31 TGE   I TABLEMi identification number for structural damping coefficient
        32 UNDEF(4) None

        """
        op2: OP2Geom = self.op2
        ntotal = 140 * self.factor  # 35*4
        s = Struct(mapfmt(op2._endian + b'35i', self.size))
        nmaterials = (len(data) - n) // ntotal
        materials = []
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            (mid, *tc_tables, trho, ta1, ta2, ta3, ta4, ta5, ta6, a, tge, b, c, d, e) = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT9=%s\n' % str(out))
            assert a == 0, out
            assert b == 0, out
            assert c == 0, out
            assert d == 0, out
            assert e == 0, out
            #MATT9(mid, g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                                 #g15_table=None, g16_table=None, g22_table=None, g23_table=None,
                                 #g24_table=None, g25_table=None, g26_table=None, g33_table=None,
                                 #g34_table=None, g35_table=None, g36_table=None, g44_table=None,
                                 #g45_table=None, g46_table=None, g55_table=None, g56_table=None,
                                 #g66_table=None, rho_table=None,
                                 #a1_table=None, a2_table=None, a3_table=None,
                                 #a4_table=None, a5_table=None, a6_table=None, ge_table=None, comment='')
            mat = MATT9(mid, *tc_tables, trho, ta1, ta2, ta3, ta4, ta5, ta6, tge, comment='')
            materials.append(mat)
            n += ntotal
        return n, materials

    def junk_read_mat11(self, data: bytes, n: int) -> int:
        """
        Solid element orthotropic material property definition.
        Defines the material properties for a 3-D orthotropic material for
        isoparametric solid elements.

        Word Name Type Description

        1 MID I Material identification number
        2 E1 RS Modulus of elasticity in the longitudinal direction or 1-direction
        3 E2 RS Modulus of elasticity in the lateral direction or 2-direction
        4 E3 RS Modulus of elasticity in the thickness direction or 3-direction
        5 NU12 RS Poisson's ratio (ε2/ε1 for uniaxial loading in the 1-direction)
        6 NU13 RS Poisson's ratio (ε3/ε1 for uniaxial loading in the 1-direction)
        7 NU23 RS Poisson's ratio (ε3/ε2 for uniaxial loading in the 2-direction)
        8 G12 RS In-plane shear modulus
        9 G13 RS Transverse shear modulus for shear in the 1–3 plane
        10 G23 RS Transverse shear modulus for shear in the 2–3 plane
        11 RHO RS Mass density
        12 A1 RS Thermal expansion coefficient in the longitudinal direction
        13 A2 RS Thermal expansion coefficient in the lateral direction
        14 A3 RS Thermal expansion coefficient in the thickness direction
        15 TREF RS Reference temperature for calculation of thermal loads
        16 GE RS Structural damping coefficient
        17 UNDEF(16) None

        """
        op2: OP2Geom = self.op2
        ntotal = 128 # 32*4
        struct1 = Struct(mapfmt(op2._endian + b'i 15f 4f 12i', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho, a1, a2, a3, tref, ge, *other = out
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT11=%s\n' % str(out))
            print(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho, a1, a2, a3, tref, ge)
            #assert a == 0, out
            #assert b == 0, out
            #assert c == 0, out
            #assert d == 0, out
            #assert e == 0, out
            #MATT9(mid, g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                                 #g15_table=None, g16_table=None, g22_table=None, g23_table=None,
                                 #g24_table=None, g25_table=None, g26_table=None, g33_table=None,
                                 #g34_table=None, g35_table=None, g36_table=None, g44_table=None,
                                 #g45_table=None, g46_table=None, g55_table=None, g56_table=None,
                                 #g66_table=None, rho_table=None,
                                 #a1_table=None, a2_table=None, a3_table=None,
                                 #a4_table=None, a5_table=None, a6_table=None, ge_table=None, comment='')
            #mat = MATT11(mid, *tc_tables, trho, ta1, ta2, ta3, ta4, ta5, ta6, tge, comment='')
            #self._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['MAT11'] = nmaterials
        op2.log.warning('geom skipping MAT11 in MPT')
        return n

    def read_matt11(self, data: bytes, n: int) -> int:
        """
        Record – MATT11(3303,33,988)

        Solid orthotropic material temperature dependence.
        Defines the temperature dependent material property for a
        3D orthotropic material for isoparametric solid elements.

        Word Name Type Description

        1 MID   I Material identification number
        2 TE1   I TABLEMi ID for modulus of elasticity in the 1-direction
        3 TE2   I TABLEMi ID for modulus of elasticity in the 2-direction
        4 TE3   I TABLEMi ID for modulus of elasticity in the 3-direction
        5 TNU12 I TABLEMi ID for Poisson's ratio (ε2/ε1 for uniaxial loading in the 1-direction)
        6 TNU13 I TABLEMi ID for Poisson's ratio (ε3/ε1 for uniaxial loading in the 1-direction)
        7 TNU23 I TABLEMi ID for Poisson's ratio (ε3/ε2 for uniaxial loading in the 2-direction)
        8 TRHO  I TABLEMi ID for mass density
        9 TG12  I TABLEMi ID for shear modulus in 1–2 plane
        10 TG13 I TABLEMi ID for shear modulus in 1–3 plane
        11 TG23 I TABLEMi ID for shear modulus in 2–3 plane
        12 TA1  I TABLEMi ID for thermal expansion coefficient in the 1-direction
        13 TA2  I TABLEMi ID for thermal expansion coefficient in the 2-direction
        14 TA3  I TABLEMi ID for thermal expansion coefficient in the 3-direction
        15 UNDEF None
        16 TGE RS TABLEMi ID for structural damping coefficient
        17 UNDEF(16) None
        ints = (1, 10, 20, 20, 30, 30, 30, 40, 40, 50, 60, 70, 70, 70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        op2: OP2Geom = self.op2
        ntotal = 128 * self.factor # 32*4
        #self.show_data(data[n:], types='if')
        struct1 = Struct(mapfmt(op2._endian + b'32i', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (mid, e1_table, e2_table, e3_table, nu12_table, nu13_table, nu23_table,
             rho_table, g12_table, g13_table, g23_table,
             a1_table, a2_table, a3_table, undef, ge_table, *other) = out
            assert undef == 0, undef
            if op2.is_debug_file:
                op2.binary_debug.write('  MATT11=%s\n' % str(out))
            assert min(other) == 0, other
            assert max(other) == 0, other
            assert sum(other) == 0, other
            mat = MATT11(mid, e1_table=e1_table, e2_table=e2_table, e3_table=e3_table,
                         nu12_table=nu12_table, nu13_table=nu13_table, nu23_table=nu23_table,
                         g12_table=g12_table, g13_table=g13_table, g23_table=g23_table,
                         a1_table=a1_table, a2_table=a2_table, a3_table=a3_table,
                         rho_table=rho_table, ge_table=ge_table)
            op2._add_methods.add_material_dependence_object(mat, allow_overwrites=False)
            n += ntotal
        op2.card_count['MATT11'] = nmaterials
        return n

    def read_matpor(self, data: bytes, n: int) -> int:
        """
        Record – MATPOR(4801,48,961)

        Material properties for porous materials used as acoustic absorbers
        Word Name Type Description
        1 MID     I Material identification number
        2 MODEL   I Absorber model
        MODEL = 1 Craggs model
        3 RHO    RS Mass density
        4 C      RS Speed of sound
        5 UNDEF(4) None
        9 RES    RS Flow resistivity
        10 POR   RS Porosity
        11 TORT  RS Tortuosity
        12 UNDEF(5) None
        MODEL = 2 Delaney-Miki model
        3 RHO       RS Mass density
        4 C         RS Speed of sound
        5 FRAME      I FRAME = 1 for rigid; FRAME = 2 for limp
        6 UNDEF(3) None
        9 RES       RS Speed of sound
        10 POR      RS Porosity
        11 UNDEF    RS
        12 DENS     RS Frame density (for limp frames only)
        13 UNDEF(4) RS
        MODEL = 3 Johnson-Champoux-Allard model
        3 RHO       RS Mass density
        4 C         RS Speed of sound
        5 FRAME      I FRAME = 1 for rigid; FRAME = 2 for limp
        6 GAMMA     RS Ratio of constant pressure specific heat to constant volume specific heat
        7 PR        RS Prandtl number
        8 MU        RS Dynamic viscosity
        9 RES       RS Speed of sound
        10 POR      RS Porosity
        11 TORT     RS Tortuosity
        12 DENS     RS Frame density (for limp frames only)
        13 L1       RS Characteristic viscous length
        14 L2       RS Characteristic thermal length
        15 UNDEF(2) None
        """
        op2: OP2Geom = self.op2
        ntotal0 = 2 * self.size # 32*4
        ntotal = 16 * self.size # 32*4
        struct0 = Struct(mapfmt(op2._endian + b'2i', self.size))
        struct1 = Struct(mapfmt(op2._endian + b'2i 2f 4i 3f 5i', self.size)) # cragg
        struct2 = Struct(mapfmt(op2._endian + b'2i 2fi 3i 2fif 4f', self.size))
        struct3 = Struct(mapfmt(op2._endian + b'2i 2fi 9f 2i', self.size))
        nmaterials = (len(data) - n) // ntotal
        for unused_i in range(nmaterials):
            edata0 = data[n:n+ntotal0]
            edata = data[n:n+ntotal]
            mid, model = struct0.unpack(edata0)
            if model == 1:
                out = struct1.unpack(edata)
                mid, model, rho, c, dummy1, dummy2, dummy3, dummy4, resistivity, porosity, tortuosity, *dummy = out # cragg
                assert (dummy1, dummy2, dummy3, dummy4) == (0, 0, 0, 0), (dummy1, dummy2, dummy3, dummy4)
                assert dummy == [0, 0, 0, 0, 0], dummy
                op2.add_matpor_craggs(mid, rho, c, resistivity, porosity, tortuosity)
            elif model == 2:
                out = struct2.unpack(edata)
                (mid, model, rho, c, frame_int, dummy1, dummy2, dummy3,
                 resistivity, porosity, dummy4, density, *dummy) = out # Delaney-Miki
                if frame_int == 1:
                    frame = 'RIGID'
                elif frame_int == 2:
                    frame = 'LIMP'
                assert (dummy1, dummy2, dummy3, dummy4) == (0, 0, 0, 0), (dummy1, dummy2, dummy3, dummy4)
                assert dummy == [0.0, 0.0, 0.0, 0.0], dummy
                #assert (dummy1, dummy2, dummy3, dummy4) == (0, 0, 0, 0), (dummy1, dummy2, dummy3, dummy4)
                #assert dummy == (0, 0, 0), dummy
                op2.add_matpor_delmiki(mid, rho, c, resistivity, porosity, frame, density=density, comment='')
            elif model == 3:
                out = struct3.unpack(edata)
                #MATPOR         2     JCA   1.205   343.0    LIMP     1.4    0.71 1.84E-5+
                #+        80000.0    0.93     2.5    30.0  1.0E-5  1.0E-4
                (mid, model, rho, c, frame_int, gamma, prandtl_number,
                 dynanmic_viscosity, resistivity, porosity, tortuosity, density,
                 L1, L2, *dummy) = out # johnson
                assert dummy == [0, 0], dummy

                # FRAME = 1 for rigid
                # FRAME = 2 for limp
                if frame_int == 1:
                    frame = 'RIGID'
                elif frame_int == 2:
                    frame = 'LIMP'
                op2.add_matpor_jca(mid, rho, c, resistivity, porosity, tortuosity, frame,
                                   gamma, prandtl_number, dynanmic_viscosity,
                                   L1, L2, density=density)
            else:  # pragma: no cover
                raise RuntimeError(model)

            if op2.is_debug_file:
                op2.binary_debug.write('  MATPOR=%s\n' % str(out))
            #print(mid, te1, te2, te3, tnu12, tnu13, tnu23, trho, tg12, tg13, tg23, ta1, ta2, ta3, blank, tge)

        op2.card_count['MATPOR'] = nmaterials
        return len(data)

# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO

    def read_radbnd(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RADBND in MPT')
        return len(data)


    def read_radm(self, data: bytes, n: int) -> int:
        """
        RADM(8802,88,413) - record 25
        .. todo:: add object
        """
        op2: OP2Geom = self.op2
        struct_i = op2.struct_i
        nmaterials = 0
        ndata = len(data)
        while n < ndata:  # 1*4
            packs = []
            edata = data[n:n+4]
            number, = struct_i.unpack(edata)
            n += 4

            iformat = b'i %if' % number
            struct_i_nf = Struct(op2._endian + iformat)
            #mid, absorb, emiss1, emiss2, ...
            ndata_per_pack = 1 + number
            nstr_per_pack = ndata_per_pack * 4

            nfields = (ndata - n) // 4
            npacks = nfields // ndata_per_pack
            for unused_ipack in range(npacks):
                edata = data[n:n+nstr_per_pack]
                pack = list(struct_i_nf.unpack(edata))
                packs.append(pack)
                n += nstr_per_pack

                mat = RADM.add_op2_data(pack)
                op2._add_methods.add_thermal_bc_object(mat, mat.radmid)
                nmaterials += 1
        op2.card_count['RADM'] = nmaterials
        return n

    def read_radmt(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping RADMT in MPT')
        return len(data)

    def read_nlparm(self, data: bytes, n: int) -> int:
        r"""
        NLPARM(3003,30,286) - record 27

        NX 2019.2
        Word Name Type Description
        1 SID       I Set identification number
        2 NINC      I Number of increments
        3 DT       RS Incremental time interval for creep analysis
        4 KMETHOD   I Method for controlling stiffness updates
        5 KSTEP     I Number of iterations before the stiffness update
        6 MAXITER   I Limit on number of iterations for each load increment
        7 CONV      I Flags to select convergence criteria
        8 INTOUT    I Intermediate output flag
        9 EPSU     RS Error tolerance for displacement U criterion
        10 EPSP    RS Error tolerance for displacement P criterion
        11 EPSW    RS Error tolerance for displacement W criterion
        12 MAXDIV   I Limit on probable divergence conditions
        13 MAXQN    I Maximum number of quasi-Newton correction vectors
        14 MAXLS    I Maximum number of line searches
        15 FSTRESS RS Fraction of effective stress
        16 LSTOL   RS Line search tolerance
        17 MAXBIS   I Maximum number of bisections
        18 MAXR    RS Maximum ratio for the adjusted arc-length increment
        19 RTOLB   RS Maximum value of incremental rotation

        ndata = 80:
                  sid nic dt   km ks max con int  epu   epp   epw   mx mx  mx fstr  lso  mx mx    rtolb
        ints    = (1, 10, 0,   1, 5, 25, -1, 0,   0.01, 0.01, 0.01, 3, 25, 4, 0.20, 0.5, 5, 20.0, 20.0, 0)
        floats  = (1, 10, 0.0, 1, 5, 25, -1, 0.0, 0.01, 0.01, 0.01, 3, 25, 4, 0.20, 0.5, 5, 20.0, 20.0, 0.0)

        # C:\MSC.Software\msc_nastran_runs\lcdf07a.op2
        ints    = (10000001, 1, 0,   1, 500, 25, 14, 0,   0.01, 0.01, 0.01, 5, 25, 0,   0.2, 0.5, 5, 20.0, 20.0, 0)
        floats  = (10000001, 1, 0.0, 1, 500, 25, 14, 0.0, 0.01, 0.01, 0.01, 5, 25, 0.0, 0.2, 0.5, 5, 20.0, 20.0, 0.0)

        """
        op2: OP2Geom = self.op2
        ndatai = (len(data) - n) // self.factor
        ndata_80 = ndatai % 80
        ndata_76 = ndatai % 76
        if ndata_80 == 0 and ndata_76:
            n, nlparms = self._read_nlparm_80(data, n)
        elif ndata_76 == 0 and ndata_80:
            n, nlparms = self._read_nlparm_76(data, n)
        elif ndata_76 == 0 and ndata_80 == 0:
            n = op2.reader_geom2._read_dual_card(
                data, n,
                self._read_nlparm_76, self._read_nlparm_80,
                'NLPARM', op2._add_methods.add_nlparm_object)
            #n = self._read_nlparm_76(data, n)
            return n
        else:
            raise NotImplementedError(f'ndatai={ndatai} ndata_76={ndata_76} ndata_80={ndata_80}')

        assert isinstance(n, int), n
        nentries = len(nlparms)
        if nentries > 0:
            op2.card_count['NLPARM'] = nentries
        return n

    def _read_nlparm_76(self, data: bytes, n: int) -> tuple[int, list[NLPARM]]:
        """
        Word Name Type Description
        1 SID       I Set identification number
        2 NINC      I Number of increments
        3 DT       RS Incremental time interval for creep analysis
        4 KMETHOD   I Method for controlling stiffness updates
        5 KSTEP     I Number of iterations before the stiffness update
        6 MAXITER   I Limit on number of iterations for each load increment
        7 CONV      I Flags to select convergence criteria
        8 INTOUT    I Intermediate output flag
        9 EPSU     RS Error tolerance for displacement U criterion
        10 EPSP    RS Error tolerance for displacement P criterion
        11 EPSW    RS Error tolerance for displacement W criterion
        12 MAXDIV   I Limit on probable divergence conditions
        13 MAXQN    I Maximum number of quasi-Newton correction vectors
        14 MAXLS    I Maximum number of line searches
        15 FSTRESS RS Fraction of effective stress
        16 LSTOL   RS Line search tolerance
        17 MAXBIS   I Maximum number of bisections
        18 MAXR    RS Maximum ratio for the adjusted arc-length increment
        19 RTOLB   RS Maximum value of incremental rotation

        """
        op2: OP2Geom = self.op2
        ntotal = 76 * self.factor  # 19*4
        s = Struct(mapfmt(op2._endian + b'iif5i3f3iffiff', self.size))
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0
        #assert ndatai % ntotal == 0
        nlparms = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            n += ntotal

            out = s.unpack(edata)
            #(sid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,maxR,rTolB) = out
            conv = out[6]
            if op2.is_debug_file:
                op2.binary_debug.write('  NLPARM=%s\n' % str(out))
            if conv in [10, 14]:
                nentries -= 1
                op2.log.warning('  skipping NLPARM=%s\n' % str(out))
                continue

            nlparm = NLPARM.add_op2_data(out)
            nlparms.append(nlparm)
        return n, nlparms

    def _read_nlparm_80(self, data: bytes, n: int) -> tuple[int, list[NLPARM]]:
        """
        Word Name Type Description
        1 SID       I Set identification number
        2 NINC      I Number of increments
        3 DT       RS Incremental time interval for creep analysis
        4 KMETHOD   I Method for controlling stiffness updates
        5 KSTEP     I Number of iterations before the stiffness update
        6 MAXITER   I Limit on number of iterations for each load increment
        7 CONV      I Flags to select convergence criteria
        8 INTOUT    I Intermediate output flag
        9 EPSU     RS Error tolerance for displacement U criterion
        10 EPSP    RS Error tolerance for displacement P criterion
        11 EPSW    RS Error tolerance for displacement W criterion
        12 MAXDIV   I Limit on probable divergence conditions
        13 MAXQN    I Maximum number of quasi-Newton correction vectors
        14 MAXLS    I Maximum number of line searches
        15 FSTRESS RS Fraction of effective stress
        16 LSTOL   RS Line search tolerance
        17 MAXBIS   I Maximum number of bisections
        18 MAXR    RS Maximum ratio for the adjusted arc-length increment
        19 RTOLB   RS Maximum value of incremental rotation
        20 ZERO  RS/I Dummy field?
        """
        op2: OP2Geom = self.op2
        ntotal = 80 * self.factor  # 20*4
        s = Struct(mapfmt(op2._endian + b'iif5i3f3iffiff i', self.size))
        ndatai = len(data) - n
        nentries = ndatai // ntotal
        assert nentries > 0
        #assert ndatai % ntotal == 0
        nlparms = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            n += ntotal

            out = s.unpack(edata)
            #(sid,ninc,dt,kMethod,kStep,maxIter,conv,intOut,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,maxR,rTolB) = out
            #conv = out[6]
            if op2.is_debug_file:
                op2.binary_debug.write('  NLPARM=%s\n' % str(out))
            #if conv in [10, 14]:
                #nentries -= 1
                #op2.log.warning('  skipping NLPARM=%s\n' % str(out))
                #continue

            nlparm = NLPARM.add_op2_data(out)
            nlparms.append(nlparm)
        return n, nlparms

    def read_nlpci(self, data: bytes, n: int) -> int:
        self.op2.log.info('geom skipping NLPCI in MPT')
        return len(data)

    def read_tstepnl(self, data: bytes, n: int) -> int:
        """Common method to read MSC/NX TSTEPNLs"""
        ndatai = (len(data) - n) * self.factor
        n108 = ndatai % 108 # nx
        n88 = ndatai % 88 # msc
        n92 = ndatai % 92 # msc-2
        if n108 == 0 and n88 != 0:
            n, tstepnls = self._read_tstepnl_nx(data, n)
        elif n108 != 0 and n88 == 0:
            n, tstepnls = self._read_tstepnl_msc(data, n)
        elif n108 != 0 and n88 != 0 and n92 == 0:
            n, tstepnls = self._read_tstepnl_msc_92(data, n)
            #self.show_data(data[n:])
            #n, tstepnls = self._read_tstepnl_msc(data, n)
        else:
            raise RuntimeError(f'ndatai={ndatai} n108={n108} n88={n88}')

        #n = self._read_dual_card(data, n, self._read_tstepnl_nx, self._read_tstepnl_msc,
                                 #'TSTEPNL', self._add_methods.add_tstepnl_object)
        return n

    def _read_tstepnl_nx(self, data: bytes, n: int) -> tuple[int, list[TSTEPNL]]:
        """
        TSTEPNL(3103,31,337) - record 29

        NX 2019.2
        23 KDAMP    I Flags to include differential stiffness to form structural damping
        24 KUPDATE  I Method for dynamic matrix update
        25 KUSTEP   I Number of iterations before the stiffness update
        26 TINTOPT  I Time integration method
        27 GAMMA   RS Amplitude decay factor for 2nd order transient integration

        """
        op2: OP2Geom = self.op2
        ntotal = 108 * self.factor  # 27*4
        s = Struct(mapfmt(op2._endian + b'iif5i3f3if3i4f 4if', self.size))
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries
        tstepnls = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out
            method = out[4]
            if method in [4]:
                op2.log.warning('method=4; skipping TSTEPNL=%r' % str(out))
            else:
                tstepnl = TSTEPNL.add_op2_data(out)
                tstepnls.append(tstepnl)
            n += ntotal
        return n, tstepnls

    def _read_tstepnl_msc(self, data: bytes, n: int) -> tuple[int, list[TSTEPNL]]:
        """
        TSTEPNL(3103,31,337) - record 29

        MSC 2005.2
        1 SID       I Set identification number
        2 NDT       I Number of time steps of value DT
        3 DT       RS Time increment
        4 NO        I Time step interval for output
        5 METHOD    I Method for dynamic matrix update
        6 KSTEP     I Time step interval or number of converged bisections
        7 MAXITER   I Limit on number of iterations
        8 CONV      I Flags to select convergence criteria
        9 EPSU     RS Error tolerance for displacement U criterion
        10 EPSP    RS Error tolerance for displacement P criterion
        11 EPSW    RS Error tolerance for displacement W criterion
        12 MAXDIV   I Limit on probable divergence conditions
        13 MAXQN    I Maximum number of quasi-Newton correction vectors
        14 MAXLS    I Maximum number of line searches
        15 FSTRESS RS Fraction of effective stress
        16 MAXBIS   I Maximum number of bisections
        17 ADJUST   I Time step skip factor for automatic time step adjustment
        18 MSTEP    I Number of steps to obtain the dominant period response
        19 RB      RS Define bounds for maintaining the same time step
        20 MAXR    RS Maximum ratio for the adjusted arc-length increment
        21 UTOL    RS Tolerance on displacement or temperature increment
        22 RTOLB   RS Maximum value of incremental rotation

        """
        op2: OP2Geom = self.op2
        ntotal = 88 * self.factor  # 22*4
        s = Struct(mapfmt(op2._endian + b'iif5i3f3if3i4f', self.size))
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries
        tstepnls = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out
            method = out[4]
            if method in [4]:
                op2.log.warning('method=4; skipping TSTEPNL=%r' % str(out))
            else:
                tstepnl = TSTEPNL.add_op2_data(out)
                tstepnls.append(tstepnl)
            n += ntotal
        return n, tstepnls

    def _read_tstepnl_msc_92(self, data: bytes, n: int) -> tuple[int, list[TSTEPNL]]:
        r"""
        TSTEPNL(3103,31,337) - record 29

        MSC 2005.2
        1 SID       I Set identification number
        2 NDT       I Number of time steps of value DT
        3 DT       RS Time increment
        4 NO        I Time step interval for output
        5 METHOD    I Method for dynamic matrix update
        6 KSTEP     I Time step interval or number of converged bisections
        7 MAXITER   I Limit on number of iterations
        8 CONV      I Flags to select convergence criteria
        9 EPSU     RS Error tolerance for displacement U criterion
        10 EPSP    RS Error tolerance for displacement P criterion
        11 EPSW    RS Error tolerance for displacement W criterion
        12 MAXDIV   I Limit on probable divergence conditions
        13 MAXQN    I Maximum number of quasi-Newton correction vectors
        14 MAXLS    I Maximum number of line searches
        15 FSTRESS RS Fraction of effective stress
        16 MAXBIS   I Maximum number of bisections
        17 ADJUST   I Time step skip factor for automatic time step adjustment
        18 MSTEP    I Number of steps to obtain the dominant period response
        19 RB      RS Define bounds for maintaining the same time step
        20 MAXR    RS Maximum ratio for the adjusted arc-length increment
        21 UTOL    RS Tolerance on displacement or temperature increment
        22 RTOLB   RS Maximum value of incremental rotation


        TSTEPNL  1       10      .01     1       ADAPT   2       10     U
        tstepnl,2,300,1.0e-6,1,adapt,,,u
        TSTEPNL  3  10  .001        ADAPT  2  10  Upw
        TSTEPNL  4  10  .0001       ADAPT  2  10  Upw
        +									+
        +		9990
        TSTEPNL  5       10      .2      1       ADAPT   2       10     U
        TSTEPNL ID     NDT    DT    NO     METHOD KSTEP MAXITER  CONV
                EPSU   EPSP   EPSW  MAXDIV MAXQN  MAXLS FSTRESS
                MAXBIS ADJUST MSTEP RB     MAXR   UTOL  RTOLB    MINITER

                   2i     f        5i              3f                   3i        f    3i            4f                   i
        ints    = (1, 10, 0.01,    1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0,   0.6, 32.0, 0.1, 20.0, 0,
                   2, 300, 1.0e-8, 1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0,   0.6, 32.0, 0.1, 20.0, 0,
                   3, 10, 0.001,   1, 3, 2, 10, 7, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0,   0.6, 32.0, 0.1, 20.0, 0,
                   4, 10, 1.e-06,  1, 3, 2, 10, 7, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 9990, 0,   0.6, 32.0, 0.1, 20.0, 0,
                   5, 10, 0.2,     1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0,   0.6, 32.0, 0.1, 20.0, 0)
        floats  = (1, 10,  0.01,   1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0.0, 0.6, 32.0, 0.1, 20.0, 0.0,
                   2, 300, 1.0e-8, 1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0.0, 0.6, 32.0, 0.1, 20.0, 0.0,
                   3, 10, 0.001,   1, 3, 2, 10, 7, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0.0, 0.6, 32.0, 0.1, 20.0, 0.0,
                   4, 10, 1.e-06,  1, 3, 2, 10, 7, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 9990, 0.0, 0.6, 32.0, 0.1, 20.0, 0.0,
                   5, 10, 0.2,     1, 3, 2, 10, 4, 0.01, 0.001, 1.0e-8, 2, 10, 2, 0.2, 5, 5,    0.0, 0.6, 32.0, 0.1, 20.0, 0.0)

        C:\MSC.Software\msc_nastran_runs\lcdf03p.4.op2
        ints    = (100000000, 1, 1.0, 1, 1, 500, 25, 10, 0.001, 0.001, 1.0e-7, 3, 25, 0,   0.2, 5, 5, 20, 0.75, 20.0, 0.1, 20.0, 0)
        floats  = (100000000, 1, 1.0, 1, 1, 500, 25, 10, 0.001, 0.001, 1.0e-7, 3, 25, 0.0, 0.2, 5, 5, 20, 0.75, 20.0, 0.1, 20.0, 0.0)

        """
        op2: OP2Geom = self.op2
        ntotal = 92 * self.factor  # 23*4
        #s = Struct(mapfmt(op2._endian + b'iif5i3f3if3i4f', self.size))
        s = Struct(mapfmt(op2._endian + b'2i f 5i 3f 3i f 3i 4f i', self.size))

        #self.show_data(data, types='ifs')
        nentries = (len(data) - n) // ntotal
        assert (len(data) - n) % ntotal == 0
        assert nentries > 0, nentries
        tstepnls = []
        for unused_i in range(nentries):
            edata = data[n:n+ntotal]
            out = s.unpack(edata)
            #print(out)
            (tstep_id, ndt, dt, no, method_int, kstep, max_iter, conv_int,
             eps_u, eps_p, eps_w, max_div, max_qn, max_ls,
             fstress, max_bisect, adjust, mstep, rb, max_r,
             utol, rtol_b, min_iter) = out
            #(sid,ndt,dt,no,kMethod,kStep,maxIter,conv,epsU,epsP,epsW,
            # maxDiv,maxQn,maxLs,fStress,lsTol,maxBisect,adjust,mStep,rb,maxR,uTol,rTolB) = out

            try:
                method = NLPARM_KMETHOD_MAP[method_int]
            except KeyError:
                raise NotImplementedError('tstepnl=%s method_int=%r' % (tstep_id, method_int))

            try:
                conv = NLPARM_CONV_MAP[conv_int]
            except KeyError:
                raise NotImplementedError('tstepnl=%s conv_int=%r' % (tstep_id, conv_int))

            tstepnl = TSTEPNL(
                tstep_id, ndt, dt, no, method=method, kstep=kstep, max_iter=max_iter, conv=conv,
                eps_u=eps_u, eps_p=eps_p, eps_w=eps_w, max_div=max_div, max_qn=max_qn, max_ls=max_ls,
                fstress=fstress, max_bisect=max_bisect, adjust=adjust, mstep=mstep, rb=rb, max_r=max_r,
                utol=utol, rtol_b=rtol_b, min_iter=min_iter, comment='')
            tstepnl.validate()
            #method = out[4]
            #if method in [4]:
                #op2.log.warning('method=4; skipping TSTEPNL=%r' % str(out))
            #else:
                #tstepnl = TSTEPNL.add_op2_data(out)
            tstepnls.append(tstepnl)
            n += ntotal
        return n, tstepnls
