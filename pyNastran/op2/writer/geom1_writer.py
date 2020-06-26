from collections import defaultdict
from struct import pack, Struct
from pyNastran.op2.errors import SixtyFourBitError

def write_geom1(op2, op2_ascii, obj, endian=b'<'):
    #if not hasattr(obj, 'nodes'):
        #return
    if not hasattr(obj, 'nodes'):
        return
    nnodes = len(obj.nodes)
    ncoords = len(obj.coords)
    ngeom1 = nnodes or ncoords
    if not ngeom1:
        return
    write_geom_header(b'GEOM1', op2, op2_ascii)
    itable = -3

    if nnodes:
        max_nid = max(obj.nodes)
        if max_nid > 99999999:  #  is the max 2147483647?  2^31-1
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max GRID nid={max_nid}')

        #nvalues = nnodes * 8
        #nbytes = nvalues * 4
        #assert nnodes == 72, nnodes
        nfields = 8 # nid, cp, x, y, z, cd, ps, seid
        nvalues = nfields * nnodes + 3 # 3 comes from the keys
        #assert nbytes == 2316, nbytes
        #op2.write(pack('6i', *[4, 0, 4, 4, 1, 4]))

        key = (4501, 45, 1)
        nbytes = write_block(op2, op2_ascii, nvalues, key)

        spack = Struct('ii 3f 3i')
        for unused_nid, node in sorted(obj.nodes.items()):
            xyz = node.xyz
            ps = node.ps
            if ps == '':
                psi = 0
            else:
                psi = int(ps)

            seid = node.seid
            if seid == '':
                seidi = 0
            else:
                seidi = int(seid)
            data = [node.nid, node.Cp(), xyz[0], xyz[1], xyz[2], node.Cd(), psi, seidi]
            op2.write(spack.pack(*data))
            op2_ascii.write('  nid=%s cp=%s xyz=(%s, %s, %s) cd=%s ps=%s seid=%s\n' % tuple(data))
        op2.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')
        #-------------------------------------

    if ncoords:
        out = defaultdict(list)
        for cid, coord in obj.coords.items():
            if coord.type == 'GMCORD':
                obj.log.warning(f'skipping {coord.type}')
                continue
            out[coord.type].append(cid)

        coord_type_key_map = {
            'CORD1C' : (1701, 17, 6),
            'CORD1R' : (1801, 18, 5),
            'CORD1S' : (1901, 19, 7),
            'CORD2C' : (2001, 20, 9),
            'CORD2R' : (2101, 21, 8),
            'CORD2S' : (2201, 22, 10),
            'CORD3G' : (14301, 143, 651),
        }
        for coord_type, cids in sorted(out.items()):
            max_cid = max(cids)
            if max_cid > 99999999:  #  is the max 2147483647?  2^31-1
                raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max {coord_type}={max_cid}')

            key = coord_type_key_map[coord_type]
            ncards = len(cids)
            if '2' in coord_type:
                coord_int = 2
            elif '1' in coord_type:
                coord_int = 1
            else:  # pragma: no cover
                raise NotImplementedError(coord_type)

            if coord_type[-1] == 'R':
                coord_rcs_int = 1
            elif coord_type[-1] == 'C':
                coord_rcs_int = 2
            elif coord_type[-1] == 'S':
                coord_rcs_int = 3
            else:  # pragma: no cover
                raise NotImplementedError(coord_type)

            if coord_type in ['CORD2R', 'CORD2C', 'CORD2S']:
                nvalues = 13 * ncards + 3
                spack = Struct(b'4i 9f')
                nbytes = write_block(op2, op2_ascii, nvalues, key)

                for cid in sorted(cids):
                    coord = obj.coords[cid]
                    data = ([cid, coord_rcs_int, coord_int, coord.Rid(), ] +
                            list(coord.e1) + list(coord.e2) + list(coord.e3))
                    op2.write(spack.pack(*data))
                    op2_ascii.write(' cid=%s data=%s' % (cid, str(data[1:])))
            elif coord_type in ['CORD1R', 'CORD1C', 'CORD1S']:
                nvalues = 6 * ncards + 3
                spack = Struct(b'6i')
                nbytes = write_block(op2, op2_ascii, nvalues, key)
                nids = []
                for cid in cids:
                    coord = obj.coords[cid]
                    nids.extend([coord.G1(), coord.G2(), coord.G3()])
                max_nid = max(nids)
                if max_nid > 99999999:
                    raise SixtyFourBitError(f'64-bit OP2 writing is not supported; {coord_type}: max nid={max_nid}')
                del nids

                for cid in sorted(cids):
                    coord = obj.coords[cid]
                    data = [cid, coord_rcs_int, coord_int, coord.G1(), coord.G2(), coord.G3()]
                    op2.write(spack.pack(*data))
                    op2_ascii.write(' cid=%s data=%s' % (cid, str(data[1:])))
            else:
                raise NotImplementedError(coord_type)
            op2.write(pack('i', nbytes))
            itable -= 1
            data = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4]
            op2.write(pack('9i', *data))
            op2_ascii.write(str(data) + '\n')

    #_write_markers(op2, op2_ascii, [2, 4])
    #-------------------------------------
    close_geom_table(op2, op2_ascii, itable)

def write_block(op2, op2_ascii, nvalues, key):
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))
    op2.write(pack('3i', *key))
    op2_ascii.write(str(key) + '\n')
    return nbytes

def init_table(table_name):
    data = [
        4, 2, 4,
        #4, 2,4,
        8, b'%8s' % table_name, 8,
        4, -1, 4,
        #4, 1, 4,
        #4, 0, 4,
    ]
    return data

def write_geom_header(table_name, op2, op2_ascii, endian=b'<'):
    op2_ascii.write('----------\n')
    data = init_table(table_name)
    op2.write(pack('4i 8s i 3i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        4, 7, 4,
        28, 1, 2, 3, 4, 5, 6, 7, 28,
    ]
    #struct_3i = Struct(endian + b'3i')
    op2.write(pack('3i 9i', *data))
    op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    data = [
        4, -2, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        #4, 0, 4,
        4, 2, 4,
        8, 1, 2, 8,
    ]
    op2.write(pack('3i 4i', *data))
    op2_ascii.write(str(data) + '\n')
    #data = [8, 1, 2, 8]
    #op2.write(pack('4i', *data))
    #-------------------------------------

    data = [
        4, -3, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')


def close_geom_table(op2, op2_ascii, itable):
    data = [
        4, 3, 4,
        12, 1, 2, 3, 12]
    op2.write(pack('3i 5i', *data))
    op2_ascii.write(str(data) + '\n')
    itable -= 1
    #-------------------------------------

    data = [
        4, itable, 4,
        4, 1, 4,
        4, 0, 4]
    op2.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        4, 0, 4,
        #4, 2, 4
    ]
    op2.write(pack('3i', *data))
    op2_ascii.write(str(data) + '\n')
    itable -= 1

def fill_defaultdict(typed_dict, direct_dict=None):  # pragma: no cover
    if not isinstance(typed_dict, tuple):
        typed_dict = (typed_dict, )

    out = defaultdict(list)
    for dicti in typed_dict:
        for idi, obj in dicti.items():
            out[obj.type].append(idi)
    if direct_dict:
        if not isinstance(direct_dict, tuple):
            direct_dict = (direct_dict, )
        for dicti in direct_dict:
            values = list(dicti.keys())
            value0 = values[0]
            obj0 = dicti[value0]
            out[obj0.type] = values
    return out

