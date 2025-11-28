from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import (
        ComplexPlateStressArray, ComplexPlateStrainArray)

    from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import (
        ComplexPlateVMStressArray, ComplexPlateVMStrainArray)

    from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
    from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, \
        RandomPlateVMStrainArray
    from pyNastran.op2.op2 import OP2


def oes_cquad4_33_complex_vm_17(op2: OP2, data: bytes,
                                obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                                nelements: int, ntotal: int,
                                is_magnitude_phase: bool) -> int:
    """
    OESVM1/2 - Table of element stresses or strains with von Mises
    OSTRVM1/2  for frequency response results.

    """
    #self.log.info(f'nelements = {nelements}')
    n = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i, von_mises1,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i, von_mises2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, von_mises1,
                   fd2, sx2, sy2, txy2, von_mises2)
        n += ntotal
    return n


def oes_cquad4_33_complex_15(op2: OP2,
                             data: bytes,
                             obj: ComplexPlateStressArray | ComplexPlateStrainArray,
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool) -> int:
    """
    NX/MSC 2018.2
    2 FD1    RS Z1 = Fibre distance
    3 EX1R   RS Normal in x at Z1
    4 EX1I   RS Normal in x at Z1
    5 EY1R   RS Normal in y at Z1
    6 EY1I   RS Normal in y at Z1
    7 EXY1R  RS Shear in xy at Z1
    8 EXY1I  RS Shear in xy at Z1
    9 FD2    RS Z2 = Fibre distance
    10 EX2R  RS Normal in x at Z2
    11 EX2I  RS Normal in x at Z2
    12 EY2R  RS Normal in y at Z2
    13 EY2I  RS Normal in y at Z2
    14 EXY2R RS Shear in xy at Z2
    15 EXY2I RS Shear in xy at Z2
    """
    n = 0
    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'14f')
    cen = 0  # 'CEN/4'
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*15=60
        n += ntotal
        out = s1.unpack(edata)  # 15
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=%s\n' % (eid, str(out)))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_random_vm_57(op2: OP2,
                            data: bytes,
                            obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                            nelements: int, ntotal: int, nnodes: int,
                            dt: Any) -> int:
    """
    name    ntotal numwide
    CQUAD8  228    57
    CTRIA6  184
    """

    #ntotal = 228 # 57*4
    n = 0
    cen = 0
    if op2.size == 4:
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
    else:
        print('ntotal =', ntotal, 228*2)
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
        raise NotImplementedError('64-bit')

    nnodes_corners = nnodes - 1
    #print('random-57: nelements =', nelements)
    factor = op2.factor

    ntotal1 = 52 * factor
    ntotal2 = 44 * factor
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        ni = 0
        edata = data[n:n+ntotal]
        #self.show_data(edata)
        #out = struct1.unpack(edata)  # len=17*4
        (eid_device, grid1, four,
         fd1, ox1, oy1, txy1, vm1,
         fd2, ox2, oy2, txy2, vm2) = struct1.unpack(edata[ni:ni+ntotal1])  # 13*4=52
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, cen)
        add_sort_x(dt, eid, cen,
                   fd1, ox2, oy2, txy2, vm1,
                   fd2, ox2, oy2, txy2, vm2)
        ni += ntotal1
        #print(eid, grid1,
              #fd1, ox1, oy1, txy1, vm1,
              #fd2, ox2, oy2, txy2, vm2)

        for unused_jnode in range(nnodes_corners): # 4*44=176
            (grid,
             fd1, ox1, oy1, txy1, vm1,
             fd2, ox2, oy2, txy2, vm2) = struct2.unpack(edata[ni:ni+ntotal2])  # 11*4=44
            #print(eid, grid)
            ni += ntotal2
            #print(grid,
                  #fd1, ox1, oy1, txy1, vm1,
                  #fd2, ox2, oy2, txy2, vm2)
            add_sort_x(dt, eid, grid,
                       fd1, ox2, oy2, txy2, vm1,
                       fd2, ox2, oy2, txy2, vm2)
        assert ni == ntotal, f'ni={ni}'
        n += ntotal
        #print()
    #assert n == ndata
    #n = self._not_implemented_or_skip(data, ndata, msg)
    #nelements = None
    #ntotal = None
    return n


def oes_cquad4_144_complex_77(op2: OP2,
                              data: bytes,
                              obj: ComplexPlateStressArray | ComplexPlateStrainArray,
                              nelements: int, nnodes: int,
                              dt: Any,
                              is_magnitude_phase: bool) -> int:
    #print('nelements =', nelements)
    n = 0
    grid_center = 0
    if op2.size == 4:
        s1 = op2.struct_2i  # 2
    else:
        s1 = op2.struct_2q  # 2
    s2 = Struct(op2._endian + mapfmt(b'i14f', op2.size)) # 15

    ntotal1 = 8 * op2.factor
    ntotal2 = 60 * op2.factor  # 4*15

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        (eid_device, _) = s1.unpack(data[n:n+ntotal1])
        n += ntotal1

        edata = data[n:n+ntotal2]
        n += ntotal2
        out = s2.unpack(edata)  # len=15*4

        eid = op2._check_id(eid_device, 'element', out)
        if op2.is_debug_file:
            op2.binary_debug.write('%s\n' % (str(out)))
        (grid,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
        #grid_center = 'CEN/%i' % grid   # this is correct, but fails

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sy1 = complex(sy1r, sy1i)
            sx2 = complex(sx2r, sx2i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)

        for unused_node_id in range(nnodes):  # nodes pts
            edata = data[n:n + ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('%s\n' % (str(out)))
            (grid,
             fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if is_magnitude_phase:
                sx1 = polar_to_real_imag(sx1r, sx1i)
                sx2 = polar_to_real_imag(sx2r, sx2i)
                sy1 = polar_to_real_imag(sy1r, sy1i)
                sy2 = polar_to_real_imag(sy2r, sy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sx2 = complex(sx2r, sx2i)
                sy1 = complex(sy1r, sy1i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_complex_vm_87(op2: OP2, data: bytes,
                             obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                             nelements: int, nnodes: int,
                             is_magnitude_phase: bool) -> int:
    """
    etype name      numwide nodes
    33  CQUAD4-33 -> 17      (not in this function)
    82  CQUADR    -> 87      (1+4 nodes)
    144 CQUAD144  -> 87      (1+4 nodes)
    ??? CTRIA6    -> 70?     (1+3 nodes)

    2 TERM    CHAR4
    3 GRID    I
    4 FD1     RS Fiber distance at z1
    5 EX1R    RS Normal in x at z1
    6 EX1I    RS Normal in x at z1
    7 EY1R    RS Normal in y at z1
    8 EY1I    RS Normal in y at z1
    9 ETXY1R  RS Shear in xy at z1
    10 ETXY1I RS Shear in xy at z1
    11 VM1    RS von Mises at z1
    12 FD2    RS Fiber distance at z2
    13 EX2R   RS Normal in x at z2
    14 EX2I   RS Normal in x at z2
    15 EY2R   RS Normal in y at z2
    16 EY2I   RS Normal in y at z2
    17 ETXY2R RS Shear in xy at z2
    18 ETXY2I RS Shear in xy at z2
    19 VM2    RS von Mises at z2

    """
    #strings = (b'20.0 CEN/ 4 \xcd\xccL\xbd\x91@\x1aD\xebX\xac\xc1b\x80FB\x96\xc9\xdd\xbf\xe7\xbe\n,\xaah\xc5\xa9b\x87\x14D\xcd\xccL=\x91@\x1a\xc4\xebX\xacAb\x80F\xc2\x96\xc9\xdd?\xe7\xbe\n\xac\xaah\xc5)b\x87\x14D\x01\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D\x02\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x03\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x04\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D',)
    # ints    = (20.0, CEN/, 4,
    #           -0.05, 617.0, -21.5, 49.6, -1.7, 738901735,              -1446680406,            1142196066,
    #           0.05, -617.0,  21.5, -49.6, 1.7, -1408581913, 700803242, 1142196066, 1, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618, 2, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 3, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 4, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618)
    #floats  = (20.0, CEN/, 4,
    #           -0.05,  617.0,   -21.5,   49.6,  -1.7,  1.9716951595721843e-12, -8.766713754677219e-14,  594.1,
    #            0.05, -617.0,    21.5,  -49.6,   1.7, -1.9716951595721843e-12,  8.766713754677219e-14,  594.1,
    #        1.4012984e-45,
    #           -0.05,  1278.68, -44.6,  198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #            0.05, -1278.68,  44.6, -198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6,
    #        2.802596928649634e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        4.203895392974451e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        5.605193857299268e-45,
    #           -0.05, 1278.6,    -44.6, 198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #           0.05, -1278.6,     44.6,-198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6)
    # -------------------------
    #ints    = (0,    'CEN/',
    #           5,     -0.025, 1131762990, -1056779073, 1137493831, -1051039598, 1133740074, -1055043205, 1142680183,
    #                   0.025, -1016740689, 1089802185, -1009242998, 1097240578, -1013972178, 1092196624, 1142689762,
    #           80000, -0.025, -1039337422, 1066818106, -1022560206, 1083595322, 1118381010, -1070307877, 1128400256,
    #                   0.025, 1111837962, -1076728084, 1128615178, -1059950868, -1029871408, 1076355832, 1130378350,
    #           80002, -0.025, 1144778127, -1043828805, 1161555343, -1027051589, 1134339017, -1054404396, 1160567412,
    #                   0.025, -1002241635, 1104149606, -985464419, 1120926822, -1014451576, 1091685317, 1160927605,
    #           80202, -0.025, -1013123432, 1093101862, -996473336, 1109743497, 1143660503, -1045020817, 1154190671,
    #                   0.025, 1133381240, -1055425923, 1150079339, -1038733090, -1002862704, 1103487198, 1154274046,
    #           80200, -0.025, 1141944942, -1046850567, 1095277553, -1093274063, 1120227570, -1068338410, 1142258188,
    #                   0.025, -1007929795, 1098641187, 1116060214, -1072783145, -1033402643, 1073147904, 1140969037,
    #           2.5,  'CEN/',
    #           5,     -0.025, 1131832225, -1056704978, 1137555116, -1050908425, 1133787322, -1054942075, 1142729036,
    #                   0.025, -1016671483, 1089950313, -1009181671, 1097371840, -1013924936, 1092297741, 1142738537,
    #           80000, -0.025, -1039289241, 1066921230, -1022512025, 1083698446, 1118448079, -1070164312, 1128472712,
    #                   0.025, 1111886248, -1076624736, 1128663464, -1059847520, -1029804359, 1076499355, 1130448349,
    #           80002, -0.025, 1144846119, -1043683269, 1161623335, -1026906053, 1134386659, -1054302422, 1160629636,
    #                   0.025, -1002173624, 1104295182, -985396408, 1121072398, -1014403979, 1091787195, 1160989873,
    #           80202, -0.025, -1013031023, 1093299683, -996403997, 1109891922, 1143714425, -1044905405, 1154266576,
    #                   0.025, 1133473626, -1055228152, 1150148641, -1038584742, -1002808755, 1103602670, 1154348580,
    #           80200, -0.025, 1141995401, -1046742558, 1095362185, -1093092785, 1120294524, -1068195089, 1142310151,
    #                   0.025, -1007828983, 1098856980, 1116049937, -1072805159, -1033268970, 1073434045, 1141020308,
    #           10.0, 'CEN/',
    #           5,     -0.025, 1132706983, -1055459734, 1138546383, -1048703317, 1134552107, -1053240134, 1143519784,
    #                   0.025, -1015552254, 1091479459, -1008189778, 1099242966, -1013160238, 1093999495, 1143528113,
    #           80000, -0.025, -1038510246, 1068653763, -1021733030, 1085430979, 1119535604, -1067741982, 1129646715,
    #                   0.025, 1112666820, -1074888834, 1129444036, -1058111618, -1028717126, 1078921059, 1131585244,
    #           80002, -0.025,  1145947464, -1041231474, 1162724680, -1024454258, 1135157762, -1052586493, 1161637710,
    #                   0.025, -1001071999, 1106747574, -984294783, 1123524790, -1013633545, 1093501695, 1161998610,
    #           80202, -0.025, -1011529318,  1096648297, -995279005, 1112398467, 1144586718, -1042964811, 1155498288,
    #                   0.025, 1134974980, -1051880286, 1151273087, -1036079365, -1001936046, 1105544152, 1155559868,
    #           80200, -0.025,  1142813256, -1044921274, 1096762940, -1090229644, 1121380799, -1065774838, 1143152359,
    #                   0.025, -1006413905, 1100701913, 1115879378, -1073189541, -1031449359, 1076004431, 1141852280)
    #floats  = (0.0,  'CEN/',
    #           5,       -0.025,  245.33273315429688, -8.176939964294434,   409.5568542480469, -13.650529861450195,  295.00128173828125,  -9.832392692565918,
    #           623.8668, 0.025, -229.76829528808594,  7.658176898956299, -432.34796142578125,  14.410158157348633, -288.02484130859375,   9.599868774414062, 624.2481689453125,
    #           1.12-40, -0.025, -35.24237823486328,   1.174628496170044, -140.96951293945312,   4.698513984680176,   84.56996154785156, -2.8187167644500732,
    #           194.4375, 0.025,  49.325233459472656, -1.644010066986084,  197.30093383789062,  -6.576040267944336,   -78.7047119140625,   2.623228073120117, 224.20480346679688,
    #           1.12106, -0.025,  751.7118530273438, -25.054555892944336,   3006.847412109375, -100.21822357177734,   313.2795715332031, -10.441608428955078,
    #           2765.65,  0.025, -780.0252075195312, 25.998241424560547, -3120.100830078125, 103.99296569824219, -273.394775390625, 9.112248420715332, 2853.591064453125,
    #           1.1238,  -0.025, -313.926513671875, 10.463171005249023, -1240.1884765625, 41.33548355102539, 683.4974975585938, -22.780973434448242,
    #           1628.4,   0.025, 284.050537109375, -9.46740436553955, 1126.5443115234375, -37.54772186279297, -742.1181640625, 24.734798431396484, 1638.593505859375,
    #           1.12,    -0.025, 578.7879638671875, -19.29100227355957, 12.538071632385254, -0.4178939163684845, 98.65809631347656, -3.2882742881774902,
    #           597.90,   0.025, -472.4237365722656, 15.745882987976074, 66.86369323730469, -2.228566884994507, -57.88176345825195, 1.92919921875, 519.2234497070312,
    #           2.5,  'CEN/',
    #           5,          -0.025, 246.38917541503906, -8.247602462768555, 411.4271240234375, -13.775626182556152, 296.44317626953125, -9.928837776184082,
    #           626.671875,  0.025, -230.8242950439453, 7.728809833526611, -434.2195129394531, 14.53533935546875, -289.466552734375, 9.696301460266113, 627.2251586914062,
    #           1.1217e-40, -0.025, -35.42617416381836, 1.1869218349456787, -141.70469665527344, 4.747687339782715, 85.08165740966797, -2.852945327758789,
    #           195.153125,  0.025, 49.509429931640625, -1.6563301086425781, 198.0377197265625, -6.6253204345703125, -79.21625518798828, 2.6574466228485107, 225.2729034423828,
    #           1.1212e-40, -0.025, 755.8617553710938, -25.332143783569336, 3023.447021484375, -101.32857513427734, 314.7334899902344, -10.538858413696289,
    #           2780.65625,  0.025, -784.17626953125, 26.27590560913086, -3136.705078125, 105.10362243652344, -274.8473205566406, 9.209406852722168, 2868.793212890625,
    #           1.1358e-40, -0.025, -316.7466125488281, 10.651827812194824, -1248.6527099609375, 41.90167999267578, 686.7886352539062, -23.0011043548584,
    #           1637.40625,  0.025, 286.86993408203125, -9.656013488769531, 1135.0040283203125, -38.113624572753906, -745.4109497070312, 24.95504379272461, 1647.69189453125,
    #           1.1383e-40, -0.025, 581.8677368164062, -19.497013092041016, 12.618782997131348, -0.42329642176628113, 99.16891479492188, -3.3224446773529053,
    #           601.471875,  0.025, -475.5002746582031, 15.951679229736328, 66.78528594970703, -2.223318338394165, -58.391685485839844, 1.96330988407135, 522.352783203125,
    #           10.0, 'CEN/',
    #           5,          -0.025, 263.4738464355469, -9.435159683227539, 441.6781921386719, -15.878581047058105, 319.7825622558594, -11.551935195922852,
    #           674.140625,  0.025, -247.90237426757812, 8.915926933288574, -464.48968505859375, 16.639568328857422, -312.80328369140625, 11.319220542907715, 675.4170532226562,
    #           1.1337e-40, -0.025, -38.397804260253906, 1.3934558629989624, -153.9121704101562, 5.57382345199585, 93.37881469726562, -3.430473804473877,
    #           213.310547,  0.025, 52.48707580566406, -1.8632657527923584, 209.94830322265625, -7.453063011169434, -87.51118469238281, 3.234825849533081, 242.62054443359375,
    #           1.1202e-40, -0.025, 823.08251953125, -30.00857162475586, 3292.330078125, -120.03428649902344, 338.26568603515625, -12.17529582977295,
    #           3026296875,  0.025, -851.4141235351562, 30.953472137451172, -3405.656494140625, 123.81388854980469, -298.3591003417969, 10.844481468200684, 3115.06689453125,
    #           1.1298e-40, -0.025, -362.57501220703125, 13.845314979553223, -1385.9808349609375, 51.4633903503418, 740.0291748046875, -26.70249366760254,
    #           1788109375,  0.025, 332.6876220703125, -12.848787307739258, 1272.2655029296875, -47.67087936401367, -798.6768798828125, 28.658126831054688, 1795.55419921875,
    #           1.1293e-40, -0.025, 631.78564453125, -22.970836639404297, 13.954647064208984, -0.517249345779419, 107.45653533935547, -3.899477481842041,
    #           652.410938,  0.025, -525.3700561523438, 19.42228889465332, 65.48402404785156, -2.1316745281219482, -66.66590118408203, 2.539447546005249, 573.13232421875)
    n = 0
    ntotal1 = 8 * op2.factor  # 2*4
    ntotal2 = 68 * op2.factor # 17*4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ndata', len(data))
    #print('ntotal', ntotal1+nnodes*ntotal2)
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(f'{i}/{nelements}')
        edata1 = data[n:n+ntotal1]
        #self.show_data(edata1)
        n += ntotal1

        fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'4s', op2.size)
        #fmt2 = mapfmt(op2._endian + b'i 5f2if 5f2if', op2.size)
        fmt2 = mapfmt(op2._endian + b'i 16f', op2.size)
        s1 = Struct(fmt1)
        s2 = Struct(fmt2)

        eid_device, cen = s1.unpack(edata1)
        #self.show_data(data, types='ifs')
        assert cen == b'CEN/', (eid_device, cen)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        for inode in range(nnodes):
            edata2 = data[n:n+ntotal2]

            #DISTANCE             NORMAL-X                     NORMAL-Y                    SHEAR-XY             VON MISES
            out = s2.unpack(edata2)
            (nid,
             fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
             fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2) = out
            fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
            if inode == 0:
                nid = 0

            if is_magnitude_phase:
                oxx1 = polar_to_real_imag(oxx1r, oxx1i)
                oxx2 = polar_to_real_imag(oxx2r, oxx2i)
                oyy1 = polar_to_real_imag(oyy1r, oyy1i)
                oyy2 = polar_to_real_imag(oyy2r, oyy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                oxx1 = complex(oxx1r, oxx1i)
                oxx2 = complex(oxx2r, oxx2i)
                oyy1 = complex(oyy1r, oyy1i)
                oyy2 = complex(oyy2r, oyy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            #print((eid, dt), nid, f'{fd1:g}, {fd2:g}')
            #print('real', f'{oxx1r:g}, {oxx1i:g}, {oyy1r:g}, {oyy1i:g}, {txy1r:g}, {txy1i:g}, {ovm1:g}')
            #print('imag', oxx2i, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2)
            add_sort_x(dt, eid, nid,
                       fd1, oxx1, oyy1, txy1, ovm1,
                       fd2, oxx2, oyy2, txy2, ovm2)
            #print(out)
            #assert nid in (1, 2, 3, 4, 17), ((eid, dt), nid, fd1, fd2)
            n += ntotal2
        #print('-----------')
    return n


def oes_cquad4_33_random_9(op2: OP2, data: bytes,
                           obj: RandomPlateStressArray | RandomPlateStrainArray,
                           nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('cquad33_9 - SORT1')
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    #cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 no VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, 0,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_ctria3_real_17(op2: OP2, data: bytes,
                       obj: RealPlateStressArray | RealPlateStrainArray,
                       ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0  # 'CEN/3'
    #assert op2.sort_method == 1, op2.code_information()

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    #print(add_new_eid_sort_x)
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'16f', op2.size))
    #print('nelements =', nelements, nelements*2, obj.ntotal)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        #print('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
            #eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(dt, eid, cen,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += ntotal
    return n


def oes_quad4_33_real_17(op2: OP2, data: bytes,
                         obj: RealPlateStressArray | RealPlateStrainArray,
                         ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0 # CEN/4
    #assert op2.sort_method == 1, op2.code_information()
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size)
    struct1 = Struct(fmt)
    #if op2.sort_method == 1:
        #print('**nelements =', nelements, nelements*2, obj.ntotal)
    #else:
        #print('**ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(
            dt, eid, cen,
            fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
            fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2)
        n += ntotal
    return n

def oes_ctria3_complex_vm_17(op2: OP2,
                             data: bytes,
                             obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                             nelements: int, ntotal: int,
                             dt,
                             is_magnitude_phase: bool) -> int:
    """
    #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
    #                                                          (REAL/IMAGINARY)
    #
    #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
    #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
    #0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
    #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00

    """
    n = 0
    cen = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
                     fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2, ) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            oxx1 = polar_to_real_imag(oxx1r, oxx1i)
            oxx2 = polar_to_real_imag(oxx2r, oxx2i)
            oyy1 = polar_to_real_imag(oyy1r, oyy1i)
            oyy2 = polar_to_real_imag(oyy2r, oyy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            oxx1 = complex(oxx1r, oxx1i)
            oxx2 = complex(oxx2r, oxx2i)
            oyy1 = complex(oyy1r, oyy1i)
            oyy2 = complex(oyy2r, oyy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, oxx1, oyy1, txy1, ovm1,
                   fd2, oxx2, oyy2, txy2, ovm2)
        n += ntotal

    #msg = op2.code_information()
    #self.log.warning(f'  skipping {op2.table_name_str} {op2.element_name}')
    #msg = '%s-%s' % (op2.table_name_str, op2.element_name)
    return n


def oes_ctria3_random_9(op2: OP2, data: bytes,
                        obj: RandomPlateStressArray | RandomPlateStrainArray,
                        nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ctria3_9 - SORT1')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_cquad4_33_random_vm_11(op2: OP2, data: bytes,
                               obj: RandomPlateVMStressArray | RandomPlateVMStrainArray,
                               nelements: int, ntotal: int) -> int:
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')
    cen = 0 # CEN/4
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, ovm1,
                   fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_ctria3_random_vm_11(op2: OP2, data: bytes,
                            obj: RandomPlateStressArray | RandomPlateStrainArray,
                            nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out
        #print('CTRIA3', op2.nonlinear_factor, out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                  fd1, sx1, sy1, txy1, ovm1,
                  fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_cquad4_144_real(op2: OP2, data: bytes, ndata: int,
                        obj: RealPlateStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    if op2.size == 4:
        center_format = op2._endian + op2._analysis_code_fmt + b'4si16f'
        node_format = op2._endian + b'i16f'
    else:
        center_format = op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'8sq16d'
        node_format = op2._endian + b'q16d'
    cs = Struct(center_format)
    ns = Struct(node_format)

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    n76 = 76 * op2.factor
    n68 = 68 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+n76]

        out = cs.unpack(edata)  # len=17*4
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(f'tri eid={eid} dt={dt} nid={grid_center} unused_j={unused_j} grid={grid}')
        #print(out[:3])
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
        add_new_eid_sort_x(dt, eid, grid_center,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += n76
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + n68])
            (grid,
             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            #print(f' {inode}: eid={eid} dt={dt} nid={grid}')
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                       fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
            n += n68
    return n


def oes_cquad4_144_random(op2: OP2, data: bytes,
                          obj: RandomPlateStressArray | RandomPlateStrainArray,
                          nelements: int, nnodes: int, ndata: int) -> int:
    n = 0
    center_format = op2._endian + op2._analysis_code_fmt + b'4s i8f'
    node_format = op2._endian + b'i8f'
    cs = Struct(center_format)
    ns = Struct(node_format)

    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1,\n'
            '  #                                fd2, sx2, sy2, txy2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1,\n'
            '  #                fd2, sx2, sy2, txy2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    ntotal1 = 44 * op2.factor
    ntotal2 = 36 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        #self.show_data(edata)
        out = cs.unpack(edata)  # len=17*4
        #print(out)
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        eid_device = out[0]
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        #print(out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1, fd2, sx2, sy2, txy2)
        n += ntotal1
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + ntotal2])
            #print(out)
            (grid,
             fd1, sx1, sy1, txy1,
             fd2, sx2, sy2, txy2,) = out
            #print(out)
            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1,
                           fd2, sx2, sy2, txy2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            # leaving off grid
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
            n += ntotal2
    return n


def oes_ctria3_complex_15(op2: OP2, data: bytes,
                          obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                          nelements: int, ntotal: int,
                          is_magnitude_phase: bool,
                          element_name_type: str) -> int:
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'14f', op2.size))
    cen = 0  # CEN/3
    n = 0
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  OESC %s - eid=%i; C=[%s]\n' % (
                element_name_type, eid,
                ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sy1 = complex(sy1r, sy1i)
            sx2 = complex(sx2r, sx2i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        obj.add_sort1(dt, eid, cen,
                      fd1, sx1, sy1, txy1,
                      fd2, sx2, sy2, txy2)
        n += ntotal
    return n
