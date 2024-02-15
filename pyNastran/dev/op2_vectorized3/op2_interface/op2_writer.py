asdf
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.errors import SixtyFourBitError
#from .geom1_writer import write_geom_header, close_geom_table
integer_types = int

from pyNastran.op2.writer.op2_writer import _write_result_tables
from pyNastran.op2.writer.geom1_writer import (
    write_geom_header, close_geom_table, write_block, # write_op2_header,
    MAX_INT)
from pyNastran.op2.op2_interface.write_utils import set_table3_field, view_dtype, view_idtype_as_fdtype
from pyNastran.op2.writer.op2_writer import TrashWriter, write_op2_header, _set_skips

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran2.dev.op2_vectorized3.op2_geom import OP2


def op2_write_op2(self, op2_out_filename: str,
                  post: int=-1,
                  endian: bytes=b'<',
                  includes: Optional[list[str]]=None,
                  skips: Optional[list[str]]=None,
                  nastran_format: Optional[str]=None) -> int:
    """
    Writes an OP2 file based on the data we have stored in the object

    Parameters
    ----------
    op2_out_filename : str
        the name of the F06 file to write
    post : int; default=-1
        the PARAM,POST flag
    endian : bytes; default='<'
        little endian is strongly recommended
    includes : list[str]; default=None
        list of results to include; exclusive with skips; default for both includes/skips=None -> all included
    skips : list[str]; default=None
        list of results to skip; exclusive with includes
    nastran_format : str; default=None -> 'msc'
        supported formats: ['msc', 'nx', 'optistruct']
    #is_mag_phase : bool; default=False
        #should complex data be written using Magnitude/Phase
        #instead of Real/Imaginary (default=False; Real/Imag)
        #Real objects don't use this parameter.

    """
    if nastran_format is None:
        nastran_format = self._nastran_format
    assert nastran_format in {'msc', 'nx', 'optistruct'}, nastran_format
    skips = _set_skips(self, includes, skips)

    #print('writing %s' % op2_outname)

    if isinstance(op2_out_filename, str):
        op2_file = open(op2_out_filename, 'wb')
        #fop2_ascii = open(op2_outname + '.txt', 'w')
        fop2_ascii = TrashWriter()
        #print('op2 out = %r' % op2_outname)
        close = True
    else:
        assert isinstance(op2_out_filename, file), f'type(op2_out_filename) = {op2_out_filename}'
        op2_file = op2_out_filename
        op2_outname = op2_out_filename.name
        close = False
        #print('op2_outname =', op2_outname)

    try:
        total_case_count, table_names = _write_op2(
            op2_file, fop2_ascii, self,
            skips,
            post=post, endian=endian,
            nastran_format=nastran_format)
    except Exception:  # NotImplementedError
        if close:
            op2_file.close()
            fop2_ascii.close()
        raise
    return total_case_count, table_names

def _write_op2(op2_file, fop2_ascii, obj: OP2,
               skips: Set[str],
               post: int=-1, endian: bytes=b'<',
               nastran_format: str='nx') -> tuple[int, list[str]]:
    """actually writes the op2"""
    date = obj.date
    #op2_ascii.write('writing [3, 7, 0] header\n')

    struct_3i = Struct(endian + b'3i')
    write_op2_header(obj, op2_file, fop2_ascii, struct_3i, post=post, endian=endian)
    #if 'CASECC' not in skips:
        #write_casecc(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    obj.log.debug(f'nastran_format={nastran_format}')
    if 'GEOM1' not in skips:  # nodes
        write_geom1(op2_file, fop2_ascii, obj, endian=endian)
    #if 'GEOM2' not in skips:  # elements
        #write_geom2(op2_file, fop2_ascii, obj, endian=endian)
    #if 'GEOM3' not in skips:  # constraints
        #write_geom3(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    #if 'GEOM4' not in skips:  # loads
        #write_geom4(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    #if 'EPT' not in skips:    # properties
        #write_ept(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    #if 'MPT' not in skips:    # materials
        #write_mpt(op2_file, fop2_ascii, obj, endian=endian)

    #if 'EDT' not in skips:  # aero
        #write_edt(op2_file, fop2_ascii, obj, endian=endian, nastran_format=nastran_format)
    #if 'EDOM' not in skips:  # optimization
        #write_edom(op2_file, fop2_ascii, obj, endian=endian)
    #if 'DIT' not in skips:  # tables
        #write_dit(op2_file, fop2_ascii, obj, endian=endian)
    #if 'DYNAMIC' not in skips:
        #write_dynamic(op2_file, fop2_ascii, obj)
    if 'grid_point_weight' not in skips:
        for key, weight in obj.grid_point_weight.items():
            weight.write_op2(op2_file, fop2_ascii, date, endian=endian)

    #is_mag_phase = False
    # we write all the other tables
    # nastran puts the tables in order of the Case Control deck,
    # but we're lazy so we just hardcode the order

    case_count, table_names = _write_result_tables(obj, op2_file, fop2_ascii, struct_3i, endian, skips)
    return case_count, table_names


def write_geom1(op2_file, op2_ascii, model: BDF, endian=b'<'):
    if not hasattr(model, 'nodes'):
        return
    nnodes = len(model.grid.node_id)
    ncoords = len(model.coord.coord_id)
    ngeom1 = nnodes or ncoords
    if not ngeom1:
        return
    write_geom_header(b'GEOM1', op2_file, op2_ascii)
    itable = -3

    if nnodes:
        grid = model.grid
        max_nid = max(grid.node_id)
        if max_nid > MAX_INT:  #  is the max 2147483647?  2^31-1
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max GRID nid={max_nid}')

        #nvalues = nnodes * 8
        #nbytes = nvalues * 4
        #assert nnodes == 72, nnodes
        nfields = 8 # nid, cp, x, y, z, cd, ps, seid
        nvalues = nfields * nnodes + 3 # 3 comes from the keys
        #assert nbytes == 2316, nbytes
        #op2_file.write(pack('6i', *[4, 0, 4, 4, 1, 4]))

        key = (4501, 45, 1)
        nbytes = write_block(op2_file, op2_ascii, nvalues, key)

        #spack = Struct('ii 3f 3i')
        idtype = 'int32'
        data = np.zeros((nnodes, 7), dtype=idtype)
        data[:, 0] = grid.node_id
        data[:, 1] = grid.cp
        data[:, [2, 3, 4]] = view_dtype(grid.xyz, idtype)
        data[:, 5] = grid.cd
        data[:, 6] = grid.ps
        data[:, 7] = grid.seid
        op2_file.write(data)
        #for unused_nid, node in sorted(obj.nodes.items()):
            #xyz = node.xyz
            #ps = node.ps
            #if ps == '':
                #psi = 0
            #else:
                #psi = int(ps)

            #seid = node.seid
            #if seid == '':
                #seidi = 0
            #else:
                #seidi = int(seid)
            #data = [node.nid, node.Cp(), xyz[0], xyz[1], xyz[2], node.Cd(), psi, seidi]
            #op2_file.write(spack.pack(*data))
            #op2_ascii.write('  nid=%s cp=%s xyz=(%s, %s, %s) cd=%s ps=%s seid=%s\n' % tuple(data))
        op2_file.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2_file.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')
        #-------------------------------------

    if ncoords:
        asdf
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

