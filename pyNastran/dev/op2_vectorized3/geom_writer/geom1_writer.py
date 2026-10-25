from __future__ import annotations
from struct import pack, Struct
from typing import TYPE_CHECKING
from pyNastran.bdf import MAX_INT
from pyNastran.op2.errors import SixtyFourBitError
from pyNastran.op2.writer.geom1 import (
    write_block, init_table, write_geom_header,
    close_geom_table, fill_defaultdict)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def write_geom1v(op2_file, op2_ascii, obj: OP2, endian: bytes=b'<',
                 nastran_format: str='nx'):
    #if not hasattr(obj, 'nodes'):
        #return
    if not hasattr(obj, 'grid'):
        return

    grid = obj.grid
    nnode = len(grid)
    ncoord = 0  # len(obj.coords)
    ngeom1 = nnode or ncoord
    if not ngeom1:
        return
    cards_written = {}
    write_geom_header(b'GEOM1', op2_file, op2_ascii)
    itable = -3

    if nnode:
        max_nid = max(grid.node_id)
        if max_nid > MAX_INT:  #  is the max 2147483647?  2^31-1
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max GRID nid={max_nid}')

        nfields = 8 # nid, cp, x, y, z, cd, ps, seid
        nvalues = nfields * nnode + 3 # 3 comes from the keys
        cards_written['GRID'] = nnode

        key = (4501, 45, 1)
        nbytes = write_block(op2_file, op2_ascii, nvalues, key)

        data = columnstack_astype([
            grid.node_id, grid.cp, grid.xyz,
            grid.cd, grid.ps, grid.seid])

        op2_file.write(np.writebuffer(data.ravel())
        op2_file.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2_file.write(pack('9i', *data))
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
            cards_written[coord_type] = ncards
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
                nbytes = write_block(op2_file, op2_ascii, nvalues, key)

                for cid in sorted(cids):
                    coord = obj.coords[cid]
                    data = ([cid, coord_rcs_int, coord_int, coord.Rid(), ] +
                            list(coord.e1) + list(coord.e2) + list(coord.e3))
                    op2_file.write(spack.pack(*data))
                    op2_ascii.write(' cid=%s data=%s' % (cid, str(data[1:])))
            elif coord_type in ['CORD1R', 'CORD1C', 'CORD1S']:
                nvalues = 6 * ncards + 3
                spack = Struct(b'6i')
                nbytes = write_block(op2_file, op2_ascii, nvalues, key)
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
                    op2_file.write(spack.pack(*data))
                    op2_ascii.write(' cid=%s data=%s' % (cid, str(data[1:])))
            else:
                raise NotImplementedError(coord_type)
            op2_file.write(pack('i', nbytes))
            itable -= 1
            data = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4]
            op2_file.write(pack('9i', *data))
            op2_ascii.write(str(data) + '\n')

    #_write_markers(op2_file, op2_ascii, [2, 4])
    #-------------------------------------
    close_geom_table(op2_file, op2_ascii, itable)
    obj.log.debug(str(cards_written))


def columnstack_astype(myarrays: list[np.ndarray],
                       dtype: str='int32'):
    """
    # nid, cp, x, y, z, cd, ps, seid
    data = np.zeros((nnode, 8), dtype='int32')
    xyz = grid.xyz
    data[:, 0] = grid.node_id
    data[:, 1] = grid.cp
    data[:, 2] = xyz[:, 0]
    data[:, 3] = xyz[:, 1]
    data[:, 4] = xyz[:, 2]
    data[:, 5] = grid.cd
    data[:, 6] = grid.ps
    data[:, 7] = grid.seid
    data = columnstack_astype([
        grid.node_id, grid.cp, grid.xyz,
        grid.cd, grid.ps, grid.seid])
    """
    nrow = len(myarrays[0])
    ncols = []
    for myarray in myarrays:
        if myarray.ndim == 1:
            ncol = 1
        else:
            assert myarray.ndim == 2, myarray.shape
            ncol = myarray.shape[1]
        ncols.append(ncol)
    ncol_total = sum(ncols)
    assert ncol_total > 0, ncol_total

    data = np.zeros((nrow, ncol_total, dtype=dtype)
    dtype_name = data.dtype.name
    assert dtype_name in ['int32', 'int64', 'float32', 'float64'], dtype_name

    i = 0
    for ncol, myarray in zip(ncols, myarrays):
        mydtype_name = myarray.dtype.name
        if dtype_name != mydtype_name:
            myarray = myarray.viewvalues(dtype_name)
        
        if ncol == 1:
            data[i, :] = myarray
        else:
            data[i:i+ncol] = myarray
    return data
