"""
Defines various tables that don't fit in other sections:
  - OP2Reader
    - read_cmodeext(self)
    - read_cmodeext_helper(self)
    - read_aemonpt(self)
    - read_monitor(self)
    - read_r1tabrg(self)
    - read_hisadd(self)

    - read_cstm(self)
    - read_dit(self)
    - read_extdb(self)
    - read_fol(self)
    - read_frl(self)
    - read_gpl(self)
    - read_ibulk(self)
    - read_intmod(self)
    - read_meff(self)
    - read_omm2(self)
    - read_sdf(self)
    - read_tol(self)
    - _skip_pcompts(self)
    - _read_pcompts(self)

  - Matrix
    - _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout)
    - _skip_matrix_mat(self)
    - read_matrix(self, table_name)
    - _read_matpool_matrix(self)
    - _read_matrix_mat(self)
    - grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                                make_matrix_symmetric)

  - Others
    - _get_marker_n(self, nmarkers)
    - read_markers(self)
    - _skip_subtables(self)
    - _skip_table_helper(self)
    - _print_month(self, month, day, year, zero, one)
    - read_results_table(self)

"""
from __future__ import annotations
import os
import sys
from copy import deepcopy
from itertools import count
from struct import unpack, Struct # , error as struct_error
from typing import Optional, Callable, TYPE_CHECKING

import numpy as np
import scipy  # type: ignore

from cpylog import SimpleLogger
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import FortranMarkerError, SortCodeError, EmptyRecordError
from pyNastran.op2.result_objects.gpdt import GPDT, BGPDT
from pyNastran.op2.result_objects.eqexin import EQEXIN
from pyNastran.op2.result_objects.matrix import Matrix
from pyNastran.op2.result_objects.matrix_dict import MatrixDict
from pyNastran.op2.result_objects.design_response import DSCMCOL
from pyNastran.op2.result_objects.qualinfo import QualInfo
from pyNastran.op2.op2_interface.msc_tables import MSC_GEOM_TABLES
from pyNastran.op2.op2_interface.nx_tables import NX_GEOM_TABLES

from pyNastran.op2.op2_interface.utils import (
    mapfmt, reshape_bytes_block,
    reshape_bytes_block_size)

from .version import parse_nastran_version
from .gpdt import (_get_gpdt_nnodes2,
                   _read_gpdt_8_14, _read_gpdt_4_7, _read_gpdt_4_10)
from .extdb import _read_extdb_extdb, _read_extdb_geomx, _read_extdb_phip

from pyNastran.op2.op2_interface.utils_matpool import (
    read_matpool_dmig, read_matpool_dmig_4, read_matpool_dmig_8)
from pyNastran.op2.result_objects.monpnt import MONPNT1, MONPNT3


from pyNastran.op2.result_objects.design_response import (
    WeightResponse, DisplacementResponse, StressResponse, StrainResponse, ForceResponse,
    FlutterResponse, FractionalMassResponse, Convergence, Desvars, DSCMCOL)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

#class MinorTables:
    #def __init__(self, op2_reader):
        #self.op2_reader = op2_reader

class SubTableReadError(Exception):
    pass

GEOM_TABLES = MSC_GEOM_TABLES + NX_GEOM_TABLES

DENSE_MATRICES = [
    b'KELM', b'MELM', b'BELM',
    b'KELMP', b'MELMP',

    b'EFMASSS', b'EFMFACS', b'EFMFSMS',
    b'MEFMASS', b'MEFWTS', b'MPFACS', b'RBMASSS',
]
# https://pyyeti.readthedocs.io/en/latest/modules/nastran/generated/pyyeti.nastran.bulk.wtextseout.html
EXTSEOUT = [
    'KAA', 'MAA', 'BAA', 'K4XX', 'PA', 'GPXX', 'GDXX', 'RVAX',
    'VA', 'MUG1', 'MUG1O', 'MES1', 'MES1O', 'MEE1', 'MEE1O', 'MGPF',
    'MGPFO', 'MEF1', 'MEF1O', 'MQG1', 'MQG1O', 'MQMG1', 'MQMG1O',
]


class OP2Reader:
    """Stores methods that aren't useful to an end user"""
    def __init__(self, op2: OP2):
        #: should an h5_file be created
        self.load_as_h5 = False
        #: the h5 file object used to reduce memory usage
        self.h5_file = None
        self.size = 4

        # Hack to dump the IBULK/CASECC decks in reverse order
        # It's in reverse because that's how Nastran writes it.
        #
        # We could write it correctly, but there's a chance the
        # deck could crash.
        self._dump_deck = False

        self.op2 = op2  # type: OP2

        self.mapped_tables = {
            b'RST': (self.read_rst, 'restart file?'),
            b'GPL' : (self.read_gpl, 'grid point list'),
            b'GPLS' : (self.read_gpls, 'grid point list (superelement)'),

            # GPDT  - Grid point definition table
            b'GPDT' : (self.read_gpdt, 'grid point locations'),
            b'GPDTS' : (self.read_gpdt, 'grid point locations (superelement)'),

            # BGPDT - Basic grid point definition table.
            b'BGPDT' : (self.read_bgpdt, 'grid points in cid=0 frame'),
            b'BGPDTS' : (self.read_bgpdt, 'grid points in cid=0 (superelement)'),
            b'BGPDTOLD' : (self.read_bgpdt, 'grid points in cid=0 frame'),
            b'BGPDTVU' : (self.read_bgpdt, 'VU grid points in cid=0 frame'),

            # optimization
            b'DESCYC' : (self.read_descyc, 'design iteration'),
            b'DBCOPT' : (self.read_dbcopt, 'design variable history table'),
            b'DSCMCOL' : (self.read_dscmcol, 'creates op2_results.responses.dscmcol'),
            b'DESTAB' :  (self._read_destab, 'creates op2_results.responses.desvars'),

            #b'MEFF' : self.read_meff,
            b'INTMOD' : (self.read_intmod, '???'),
            b'HISADD' : (self.read_hisadd, 'optimization history; op2_results.responses.convergence_data'),
            b'EXTDB' : (self.read_extdb, 'external superlelements'),
            b'OMM2' : (self.read_omm2, 'max/min table'),
            b'STDISP' : (self.read_stdisp, 'aero-structural displacement?'),
            b'TOL' : (self.read_tol, 'time output list?'),
            b'PCOMPT' : (self._read_pcompts, 'NX: LAM option input from the PCOMP bulk entry'),
            b'PCOMPTS' : (self._read_pcompts, 'NX: LAM option input from the PCOMP bulk entry (superelement)'),
            b'MONITOR' : (self.read_monitor, 'MONITOR point output'),
            b'AEMONPT' : (self.read_aemonpt, 'aero matrix'),
            b'FOL' : (self.read_fol, 'frequency output list'),
            b'FRL' : (self.read_frl, 'frequency response list'),
            b'SDF' : (self.read_sdf, 'aero-structural displacement?'),
            b'IBULK' : (self.read_ibulk, 'explicit bulk data'),
            b'ICASE' : (self.read_icase, 'explicit case control'),
            b'CASECC': (self.read_casecc, 'case control'),
            b'XCASECC': (self.read_xcasecc, 'case control'),

            b'CDDATA' : (self.read_cddata, 'Cambell diagram'),
            b'CMODEXT' : (self._read_cmodext, 'Component mode synthesis - external'),

            #MSC
            #msc / units_mass_spring_damper
            b'UNITS' : (self._read_units, 'units'),
            #b'CPHSF': self._read_cphsf,

            # element matrices
            #b'KELM' : self._read_element_matrix,
            #b'MELM' : self._read_element_matrix,
            #b'BELM' : self._read_element_matrix,
            #b'KELMP' : self._read_element_matrix,
            #b'MELMP' : self._read_element_matrix,

            # element dictionaries
            b'KDICT' : (self._read_dict, 'matrix'),
            b'MDICT' : (self._read_dict, 'matrix'),
            b'BDICT' : (self._read_dict, 'matrix'),
            b'KDICTP' : (self._read_dict, 'matrix'),
            b'MDICTP' : (self._read_dict, 'matrix'),
            b'KDICTDS' : (self._read_dict, 'matrix'),
            b'KDICTX' : (self._read_dict, 'matrix'),
            b'XDICT' : (self._read_dict, 'matrix'),
            b'XDICTB' : (self._read_dict, 'matrix'),
            b'XDICTDS' : (self._read_dict, 'matrix'),
            b'XDICTX' : (self._read_dict, 'matrix'),

            # coordinate system transformation matrices
            b'CSTM' : (self.read_cstm, 'coordinate transforms'),
            b'CSTMS' : (self.read_cstm, 'coordinate transforms (superelement)'),
            b'TRMBD': (self.read_trmbd, 'euler angles for transforming from material to (deformed) basic csys'),
            b'TRMBU': (self.read_trmbu, 'euler angles for transforming from material to (undeformed) basic csys'),

            b'R1TABRG': (self.read_r1tabrg, 'DRESP1 optimization table'),
            # Qualifier info table???
            b'QUALINFO' : (self.read_qualinfo, 'Qualifier info table'),

            # Equivalence between external and internal grid/scalar numbers
            b'EQEXIN' : (self.read_eqexin, 'internal/external ids'),
            b'EQEXINS' : (self.read_eqexin, 'internal/external ids (superelement)'),

            b'XSOP2DIR' : (self.read_xsop2dir, 'list of external superelement matrices?'),
            #b'TEF1': (self._read_mef1, 'superelement matrix'),
            #b'MEF1': (self._read_mef1, 'superelement matrix'),
            #b'MUG1B': (self._read_mef1, 'superelement matrix'),
            b'MATPOOL' : (self._read_matrix_matpool, 'matrices'),

            b'OBC1': (self.read_obc1, 'Contact pressures and tractions at grid points'),
            b'OBG1': (self.read_obc1, 'Glue normal and tangential tractions at grid point in cid=0 frame'),
            b'PTMIC' : (self._read_ptmic, 'property of VATV microphone points'),

            # OVG: Table of aeroelastic x-y plot data for V-g or V-f curves
            b'MKLIST': (self._read_mklist, 'M/K aero pairs'),
        }
        desc_map = {
            b'PERF': 'aero matrix',
            b'META': 'string matrix',
            b'TUG1' : 'external superelement',
            b'TQG1' : 'external superelement',
            b'MKQG1' : 'external superelement',
            b'MATRV': 'external superelement',
            b'MUG1B': 'external superelement',
        }
        desc_map.update({key: values[1] for key, values in self.mapped_tables.items()})
        #desc_map = {key: values[1] for key, values in self.mapped_tables.items()}
        self.desc_map = desc_map

    def read_nastran_version(self, mode: str) -> None:
        """
        reads the version header
          ints = (3, 4, 12, 1, 28, 12, 12, 4, 7, 4,
          28, 1414742350, 541999442, 1414680390, 1346458656, 1145643077, 1146045216, 539828293, 28,
          4, 2, 4, 8, 1482184792, 1482184792, 8, 4, -1, 4, 4, 0, 4, 4, 2, 4, 8, 1297040711, 538976305, 8, 4, -1
          """
        #try:
        op2 = self.op2
        markers = self.get_nmarkers(1, rewind=True)
        #except Exception:
            #self._goto(0)
            #try:
                #self.f.read(4)
            #except Exception:
                #raise FatalError("The OP2 is empty.")
            #raise
        if self.is_debug_file:
            if self.read_mode == 1:
                self.binary_debug.write('read_mode = %s (vectorized; 1st pass)\n' % self.read_mode)
            elif self.read_mode == 2:
                self.binary_debug.write('read_mode = %s (vectorized; 2nd pass)\n' % self.read_mode)

        if markers == [3,]:  # PARAM, POST, -1
            if self.is_debug_file:
                self.binary_debug.write('marker = 3 -> PARAM,POST,-1?\n')
            if op2.post is None:
                op2.post = -1
            self.read_markers([3])
            data = self._read_block()   # TODO: is this the date...pretty sure at least for MSC
            ndata = len(data)
            if ndata == 4:
                one = Struct(self._endian + b'i').unpack(data)[0]
                assert one == 1, one
            elif ndata == 12:
                date = self.op2.struct_3i.unpack(data)
                op2.log.debug(f'date = {date}')
            elif ndata == 24:
                date = self.op2.struct_3q.unpack(data)
                op2.log.debug(f'date = {date}')
            else:
                #self.show_data(data, types='ifsqd', endian=None, force=False)
                assert ndata in [4, 12, 24], f'ndata={ndata} data={data}'

            self.read_markers([7])
            data = self.read_string_block()  # 'NASTRAN FORT TAPE ID CODE - '
            if data == b'NASTRAN FORT TAPE ID CODE - ':
                macro_version = 'nastran'
            elif b'IMAT v' in data:
                imat_version = data[6:11].encode('utf8')
                macro_version = 'IMAT %s' % imat_version
                mode = 'msc'
            else:
                version_ints = Struct(self._endian + b'7i').unpack(data)
                if version_ints == (1, 2, 3, 4, 5, 6, 7):
                    macro_version = 'MSFC'
                    mode = 'msc'
                else:
                    self.show_data(data)
                    raise NotImplementedError(data)

            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)
            #print('macro_version = %r' % macro_version)

            data = self._read_record()
            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)
            version = data.strip()
            #version_str = version.decode(self._encoding)
            #print('version = %r' % version_str)

            if macro_version == 'nastran':
                mode = parse_nastran_version(
                    data, version, self._encoding, self.op2.log)
                # don't uncomment this...it breaks tests
                #op2._nastran_format = mode
            elif macro_version.startswith('IMAT'):
                assert version.startswith(b'ATA'), version
                op2._nastran_format = macro_version
            elif macro_version == 'MSFC':
                op2._nastran_format = macro_version
                #self.show_data(version)
                #assert version.startswith(b'ATA'), version

            if self.is_debug_file:
                self.binary_debug.write(data.decode(self._encoding) + '\n')
            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -2
            if self.is_debug_file:
                self.binary_debug.write('marker = 2 -> PARAM,POST,-2?\n')
            if op2.post is None:
                op2.post = -2
        else:
            raise NotImplementedError(markers)

        #print(f'mode = {mode!r} fmt={op2._nastran_format!r}')
        if mode == 'autodesk' or op2._nastran_format == 'autodesk':
            op2.post = -4
            mode = 'autodesk'
        elif mode == 'nasa95' or op2._nastran_format == 'nasa95':
            op2.post = -4
            mode = 'nasa95'
        elif mode == 'optistruct' or op2._nastran_format == 'optistruct':
            op2.post = -4
            #pass # mode = 'optistruct'
        elif isinstance(op2._nastran_format, str):
            if op2._nastran_format not in {'msc', 'nx', 'optistruct'}:
                raise RuntimeError(f'nastran_format={op2._nastran_format} mode={mode} and must be "msc", "nx", "optistruct", or "autodesk"')
            mode = op2._nastran_format
        elif mode is None:
            self.log.warning("No mode was set, assuming 'msc'")
            mode = 'msc'
        if self.read_mode == 1:
            self.log.debug(f'mode = {mode!r}')
        self.op2.set_mode(mode)
        self.op2.set_table_type()


    def _read_mklist(self):
        """
        reads the MKLIST table and puts it in:
           - op2.op2_results.mklist
        """
        op2 = self.op2
        size = self.size
        unused_table_name = self._read_table_name(rewind=False)
        assert size == 4, size
        self.read_markers([-1])

        #(101, 82, 0, 0, 0, 0, 0)
        data = self._read_record()

        itable = -2
        #self.log.warning(f'MKLIST itable={itable}')
        self.read_3_markers([itable, 1, 0])
        #(b'MKLIST  ',)
        data = self._read_record()
        itable -= 1

        #self.log.warning(f'MKLIST itable={itable}')
        self.read_3_markers([itable, 1, 0])
        #(0.45, 0.00283, ...
        # 0.45, 0.02833)
        data = self._read_record()
        mk_data = np.frombuffer(data, dtype=op2.fdtype)
        nmk = len(mk_data)
        assert nmk > 0, nmk
        assert nmk % 2 == 0, nmk
        mk = mk_data.reshape(nmk // 2, 2)
        if not hasattr(op2.op2_results, 'mklist'):
            op2.op2_results.mklist = []
        op2.op2_results.mklist.append(mk)
        itable -= 1

        self.read_markers([itable, 1, 0, 0])

    def _read_ptmic(self):
        """
        5.73 PTMIC
        Table of the properties of microphone grids for ATV analysis.
        Record 0 – HEADER
        Word Name Type Description
        1 NAME(2) CHAR4 Data block name

        Record 1 – MICROPHONE POINT
        Word Name Type Description
        1 EXTID I Microphone point external ID
        2 PID   I Property ID of microphone point
        Words 1 and 2 repeat until End of Record

        Record 2 – PROPERTY LIST BY FREQUENCY
        Word Name Type Description
        1 EXTID I  Property ID of microphone point
        2 RHOR  RS Real part of density
        3 RHOI  RS Imaginary part of density
        4 CR    RS Real part of velocity
        5 CI    RS Imaginary part of velocity
        Words 1 through 5 repeat until End of Record
        """
        op2 = self.op2
        if self.read_mode == 1:
            self.read_geom_table()
            return

        size = self.size
        unused_table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._read_record()
        itable = -1
        struct1 = Struct(mapfmt(b'7i', size))

        # (104, 90, 1, 6, 0, 0, 0)
        a, nrows, b, nrecord2, c, d, e = struct1.unpack(data)
        assert a == 104, a
        assert b == 1, b
        assert c == 0, c
        assert d == 0, d
        assert e == 0, e
        self.log.warning(f'ptmic itable={itable} - (104, nrows={nrows}, 1, nrecord2={nrecord2}, 0, 0, 0)')

        # (104, 90, 1, 6, 0, 0, 0)  # q
        #op2.show_data(data, 'iq')

        itable = -2
        self.log.warning(f'ptmic itable={itable} - PTMIC')
        self.read_3_markers([itable, 1, 0])
        data = self._read_record()
        # PTMIC
        #op2.show_data(data, 's')

        itable = -3
        self.log.warning(f'ptmic itable={itable} (406, 1, 420, 1, ...)')
        self.read_3_markers([itable, 1, 0])
        data = self._read_record()
        extid_pid = np.frombuffer(data, dtype=op2.idtype8)
        #print(extid_pid, len(extid_pid))
        assert len(extid_pid) % 2 == 0, len(extid_pid)
        nrows_test = len(extid_pid) // 2
        assert nrows == nrows_test

        extid_pid = extid_pid.reshape(nrows, 2)
        #print(extid_pid, extid_pid.shape)

        # Record 2 repeats for each material type region in which a microphone point is defined.

        itable = -4

        ntotal = 5 * size
        marker = self.get_nmarkers(1, rewind=True)[0]
        while marker != 0:
            #self.log.warning(f'ptmic itable={itable}')
            self.read_3_markers([itable, 1, 0])
            ni = op2.n
            #self.show(20, types='ifsdq')
            try:
                marker_test = self.get_nmarkers(1, rewind=True)
                #print('marker_testA =', marker_test)
            except Exception:
                print('ni =', ni, op2.n)
                raise
            if marker_test == [0]:
                break

            data = self._read_record()
            #self.show_data(data, types='ifsdq')
            nvalues = len(data) // size
            ncomplex = nvalues // 5
            assert nvalues % 5 == 0
            n = 0
            struct_i4f = Struct(mapfmt(b'i 4f', size))
            for unused_i in range(ncomplex):
                edata = data[n:n+ntotal]
                out = struct_i4f.unpack(edata)
                str(out)
                #print(out)
                n += ntotal
            assert ncomplex == 1
            itable -= 1

            markers = self.get_marker1(1, True)
            if markers == -10:
                markers = self.get_marker1(1, False)

        markers_check = self.get_nmarkers(1, rewind=False)[0]

        #self.log.warning(f'ptmic itable={itable}')
        #op2.show_data(data, 'idsq')
        #op2.show(100, 'ifsqd')
        #aa

    def read_xsop2dir(self):
        """
        Matrix datablocks are named with the matrix name, and they contain
        only the matrices with no row and column DOF information.

        External superelements are written to the OP2 in the XSOP2DIR datablock
        which contains the directory of matrices, and multiple EXTDB datablocks
        that contain these matrices.  The information in XSOP2DIR can be used
        to rename the EXTDB datablocks to their corresponding matrix name.

        The MATPOOL datablock contains matrices of several different formats.
        MATPOOL matrices are DMIG-formatted matrices. Matrices from all of
        these sources can be stored as structured array.

        """
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\extse04c_cnv1_0.op2
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\extse04c_cnv1_0.op2
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\atv005mat.op2
        op2 = self.op2
        self.xsop2dir_names = []
        if self.read_mode == 1 or op2.make_geom is False:
            table_name = self._read_table_name(rewind=True)
            self._skip_table(table_name, warn=False)
            return

        #op2 = self.op2
        unused_table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        # (101, 14, 0, 0, 0, 0, 0)
        # (101, 17, 0, 0, 0, 0, 0) # modes_elements_dmig.op2
        data = self._read_record()
        values = np.frombuffer(data, dtype=op2.idtype8)
        ntables = values[1]
        #assert ntables == 7, values

        #self.read_3_markers([-2, 1, 0])
        #data = self._read_record()

        itable = -2
        self.read_3_markers([itable, 1, 0])
        #marker = self.get_marker1(rewind=True, macro_rewind=False)

        skip_tables = {'PVT', 'PVT0'}
        if self.size == 4:
            #self.show(1000, types='ifs', endian=None, force=False)
            struct_8s = Struct(self._endian + b'8s')
            struct_16s = Struct(self._endian + b'16s')
            for unused_i in range(ntables + 1):
            #while marker != 0:
                itable -= 1
                data, ndata = self._read_record_ndata()
                if ndata == 8:
                    name = struct_8s.unpack(data)[0].decode('latin1').rstrip()
                else:
                    assert ndata == 16, self.show_data(data, types='ifs', endian=None, force=False)
                    name = struct_16s.unpack(data)[0].decode('latin1').rstrip()
                self.read_3_markers4([itable, 1, 0])
                #marker = self.get_marker1_4(rewind=True, macro_rewind=False)
                if name in skip_tables:
                    self.log.warning(f'skipping {name}')
                    continue
                self.log.warning(f'{len(self.xsop2dir_names)} {name}')
                self.xsop2dir_names.append(name)
        else:
            struct_16s = Struct(self._endian + b'16s')
            for unused_i in range(ntables + 1):
            #while marker != 0:
                itable -= 1
                data = self._read_record()
                name = struct_16s.unpack(data)[0]
                name = reshape_bytes_block(name).decode('latin1').rstrip()
                self.read_3_markers([itable, 1, 0])
                #marker = self.get_marker1_8(rewind=True, macro_rewind=False)
                if name in skip_tables:
                    self.log.warning(f'skipping {name}')
                    continue
                self.log.warning(f'{len(self.xsop2dir_names)} {name}')
                self.xsop2dir_names.append(name)
        self.read_markers([0])
        #self.log.warning(f'XSOP2DIR itable={itable}')
        # (101, 14, 0, 0, 0, 0, 0)
        #b'XSOP2DIR',
        #b'PVT0    '
        #b'GEOM1EX '
        # b'GEOM2EX '
        # b'GEOM4EX '
        # b'GEOM1EXA'
        #b'MATK    '
        #b'MATM    '
        #b'MATV    '
        #b'TUG1    '
        #b'MUG1    '
        #b'TES1    '
        #b'MES1    '

    def read_eqexin(self):
        """isat_random.op2"""
        op2 = self.op2
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        fmt1 = mapfmt(self._endian + b'7i', self.size)

        idata = unpack(fmt1, data)
        #if idata[0] == 101 and op2._nastran_format == 'nasa95':
        assert idata[0] == 101, idata
        unused_nnodes = idata[1]
        assert idata[2] == 0, idata
        assert idata[3] == 0, idata
        assert idata[4] == 0, idata
        assert idata[5] == 0, idata
        assert idata[6] == 0, idata
        #else:
        op2.log.debug(f'eqexin idata={idata}')
        #print('----------------------')

        self.read_3_markers([-2, 1, 0])
        self.read_table_name(['EQEXIN', 'EQEXINS', 'EQEXNOUT', 'SAEQEXIN'])
        #print('----------------------')
        # ints - sort order
        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        eqexin1 = np.frombuffer(data, dtype=op2.idtype)

        # ints - nid, dof_type
        self.read_3_markers([-4, 1, 0])
        data = self._read_record()
        eqexin2 = np.frombuffer(data, dtype=op2.idtype)

        self.read_markers([-5, 1, 0, 0])

        result_name = 'eqexin'
        if op2._results.is_saved(result_name):
            op2._results._found_result(result_name)

            nid, dof, doftype = eqexin_to_nid_dof_doftype(eqexin1, eqexin2)
            op2.op2_results.eqexin = EQEXIN(nid, dof, doftype)
        #print('nid = %s' % nid.tolist()) [1,2,3,...]
        #print('dof = %s' % dof.tolist()) [1,7,13,...]
        #print('doftype = %s' % doftype.tolist())

    def read_aemonpt(self):
        r"""
        reads the AEMONPT table
        D:\NASA\git\examples\backup\aeroelasticity\loadf.op2
        """
        #self.log.debug("table_name = %r" % op2.table_name)
        unused_table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(self._endian + b'7i', data)
            assert a == 101, a
            assert bi in [0, 1, 3], bi
            assert c == 27, c
            assert d == 1, d
            assert e in [6, 11], e
            assert f == 0, f
            assert g == 0, g
        #print('-----------------------')
        #print('record 2')
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()

        word, = unpack(self._endian + b'8s', data)
        assert word in [b'AECFMON ', b'AEMON   '], word
        #self.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_3_markers([-3, 1, 0])
        #self.show_data(data)

        structi = Struct(self._endian + b'8s 56s 5i 4s 8s 3i')
        itable = -4
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            data = self._read_record()
            ndata = len(data)
            if ndata != 108:
                self.show_data(data, types='ifs', endian=None, force=False)
                assert ndata == 108, ndata
            else:
                n = 8 + 56 + 20 + 12 + 12
                out = structi.unpack(data[:n])
                (aero_bytes, name_bytes, comps, cp, bi, c, d, coeff, word_bytes, e, f, g) = out
                aero = aero_bytes.decode('latin1').rstrip()
                name = name_bytes.decode('latin1').rstrip()
                word = word_bytes.decode('latin1').rstrip()
                print('aero=%r' % aero)
                print('  name=%r' % name)
                print('  comps=%r cp=%s b,c,d=(%s, %s, %s)' % (comps, cp, bi, c, d))
                print('  coeff=%r' % coeff)
                print('  word=%r (e, f, g)=(%s, %s, %s)' % (word, e, f, g)) # (1, 2, 0)
                assert cp in [0, 2], cp
                assert bi == 0, bi
                assert c == 0, c
                assert d == 0, d
                assert e in [1, 7, 9], e
                assert f in [0, 2], f
                assert g == 0, g

            self.read_3_markers([itable, 1, 0])

            markers = self.get_nmarkers(1, rewind=True)
            itable -= 1
        #self.show(100, types='ifs', endian=None, force=False)
        #markers = self.get_nmarkers(1, rewind=False)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)

    def read_monitor(self):
        r"""
        reads the MONITOR table; new version
        D:\NASA\git\examples\backup\aeroelasticity\loadf.op2"""
        op2 = self.op2
        self.log.debug("table_name = %r" % op2.table_name)
        unused_table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        #(101, 2, 27, 0, 9, 0, 0)
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None, force=False)
        self.read_3_markers([-2, 1, 0])

        #b'STMON   '
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None, force=False)

        markers = self.get_nmarkers(1, rewind=True)
        self.read_3_markers([-3, 1, 0])

        itable = -4
        markers = self.get_nmarkers(1, rewind=True)
        if self.read_mode == 2:
            op2.monitor_data = []
        while markers[0] != 0:
            if self.read_mode == 1:
                self._skip_record()
            else:
                assert self.read_mode == 2, self.read_mode
                data = self._read_record()
                #self.show_data(data, types='ifs', endian=None, force=False)
                structi = Struct(self._endian + b'8s 56s 2i 3f 4s 8s 3i')
                ndata = len(data)
                assert ndata == 108, ndata
                (aero_bytes, name_bytes, comps, cp, x, y, z, coeff, word_bytes, column, cd,
                 ind_dof) = structi.unpack(data)
                aero = aero_bytes.decode('latin1').rstrip()
                name = name_bytes.decode('latin1').rstrip()
                word = word_bytes.decode('latin1').rstrip()
                print('aero=%r' % aero)
                print('  name=%r' % name)
                print('  comps=%s cp=%s (x, y, z)=(%s, %s, %s)' % (comps, cp, x, y, z))
                print('  coeff=%r' % coeff)
                print('  word=%r (column, cd, ind_dof)=(%s, %s, %s)' % (word, column, cd, ind_dof))
                assert cp in [0, 2], cp
                assert x == 0.0, x
                assert y == 0.0, y
                assert z == 0.0, z
                assert column in [1, 7], column
                assert cd in [0, 2], cd
                assert ind_dof == 0, ind_dof
                monitor = {
                    'name' : name,
                    'cp' : cp,
                    'cd' : cd,
                    'xyz' : [x, y, z],
                    'comps' : comps,
                }
                op2.monitor_data.append(monitor)
            self.read_3_markers([itable, 1, 0])

            markers = self.get_nmarkers(1, rewind=True)
            itable -= 1
        markers = self.get_nmarkers(1, rewind=False)

    def _create_objects_from_matrices(self) -> None:
        """
        creates the following objects:
          - monitor3 : MONPNT3 object from the MP3F matrix
          - monitor1 : MONPNT1 object from the PMRF, PERF, PFRF, AGRF, PGRF, AFRF matrices

        """
        op2 = self.op2
        op2_results = self.op2.op2_results
        #assert len(self._frequencies) > 0, self._frequencies
        if 'MP3F' in op2.matrices:
            op2_results.monitor3 = MONPNT3(op2._frequencies, op2.matrices['MP3F'])
        # these are totally wrong...it doesn't go by component;
        # it goes by inertial, external, flexibility, etc.
        if 'PERF' in op2.matrices and 'PGRF' in op2.matrices:
            #op2.monitor1 = MONPNT1(
                #op2._frequencies, op2.matrices, [
                # :)       ?       :)      :)      ?       ?
                #'PMRF', 'AFRF', 'PFRF', 'PGRF', 'AGRF', 'PERF', ])

            op2_results.monitor1 = MONPNT1(
                op2._frequencies, op2.matrices,
                #  :)       ?       :)      :)2     ?       ?
                ['PMRF', 'PERF', 'PFRF', 'AGRF', 'PGRF', 'AFRF', ])


    def _read_dict(self):
        """testing the KDICT"""
        op2 = self.op2
        if self.read_mode == 1:
            table_name = self._read_table_name(rewind=True)
            self._skip_table(table_name, warn=False)
            return
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)

        unused_ints = unpack(self._endian + b'7i', data)
        #print('date?  = (?, month, day, ?, ?, ?) =', ints)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        name, = unpack(self._endian + b'8s', data)
        name = name.decode('ascii').strip()
        #print('name = %r' % name)

        self.read_3_markers([-3, 1, 0])
        matdict = MatrixDict(name)
        itable = -4
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            #-------------------------------------------------------------------
            data = self._read_record()
            n0 = 0
            ndata = len(data)

            eltype, numwids, numgrid, dof_per_grid, form = unpack(self._endian + b'5i', data[:20])
            n0 += 20
            #print('etype', eltype)
            #print('numwids', numwids)
            #print('numgrid', numgrid)
            #print('dof_per_grid', dof_per_grid)
            #print('form', form)
            eids = []
            ge = []
            address = []
            xforms = []
            sils = []
            while n0 < ndata:
                #print("--------")
                #print('n0=%s ndata=%s' % (n0, ndata))
                eid, unused_nactive_grid, gei, address1, address2 = unpack(self._endian + b'2i f 2i', data[n0:n0+20])
                eids.append(eid)
                ge.append(gei)
                address.append([address1, address2])
                n0 += 20
                #print('eid', eid)
                #print('nactive_grid', nactive_grid)
                #print('ge', ge)
                #print('address1', address1)
                #print('address2', address2)
                if form in [3, 4, 5, 6]:
                    nbytes = numgrid * 4
                    sil = unpack('%ii' % numgrid, data[n0:n0+nbytes])
                    n0 += nbytes
                    #print('sil =', sil)
                    assert min(sil) >= 0, sil
                    assert max(sil) < 100000000, sil
                    sils.append(sil)
                else:
                    raise NotImplementedError(form)
                if form in [4]:
                    xform = unpack(b'9d', data[n0:n0+36*2])
                    xform_array = np.array(xform, dtype='float64').reshape(3, 3)
                    #print(xform_array)
                    xforms.append(xform_array)
                    n0 += 36*2
                    #print('xform =', xform, len(xform))
                #self.show_data(data[:n0+80], types='if')
            if len(xforms) == 0:
                xforms = None
            eids = np.array(eids, dtype='int32')
            ge = np.array(ge, dtype='float32')
            address = np.array(address, dtype='int32')
            sils = np.array(sils, dtype='int32')
            #print('eids, address:\n', eids, '\n', address)
            #print('sils (etype=%s):\n%s' % (eltype, sils))
            matdict.add(eltype, numwids, numgrid, dof_per_grid, form,
                        eids, ge, address, sils, xform=xforms)
            #-------------------------------------------------------------------
            self.read_3_markers([itable, 1, 0])
            itable -= 1
        self.read_markers([0])
        self.op2.matdicts[name] = matdict

    def _read_destab(self) -> None:
        """reads the DESTAB table"""
        op2 = self.op2

        result_name = 'responses.desvars'
        if op2._results.is_not_saved(result_name):
            data = self._skip_table('DESTAB', warn=False)
            return

        #if self.read_mode == 1:
            #return ndata
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_destab - %s\n' % op2.table_name)

        self.read_markers([-1])
        data = self._read_record()

        # (101, 3, 3, 0, 3, 0, 0)

        itable = -2
        markers = self.read_3_markers([itable, 1, 0])
        data = self._read_record()
        if self.size == 4:
            destab = op2.struct_8s.unpack(data)[0].rstrip()
            structi = Struct('2i 8s 4f')
        else:
            destab = op2.struct_16s.unpack(data)[0]
            destab = reshape_bytes_block(destab).rstrip()
            structi = Struct('2q 16s 4d')
        assert destab == b'DESTAB', destab

        itable -= 1
        markers = self.read_3_markers([itable, 1, 0])

        desvars = []
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break

            data = self._read_record()

            #self.show(100, types='ifs', endian=None)
            #self.show_data(data[:8], types='ifs', endian=None)

            #1 IDVID I Internal design variable identification number
            #2 DVID I External design variable identification number
            #3 LABEL1 CHAR4 First part of design Variable
            #4 LABEL2 CHAR4 Second part of design Variable
            #5 VMIN RS Lower bound
            #6 VMAX RS Upper bound
            #7 DELX RS Move limit for a design cycle
            # 8 ???

            #C:\NASA\m4\formats\git\examples\move_tpl\accopt3.op2
            #---------------------------------------------------------------------------------------------------------
                #INTERNAL       DESVAR                         LOWER                               UPPER
                    #ID            ID          LABEL            BOUND             VALUE             BOUND
            #---------------------------------------------------------------------------------------------------------
                #1             1      THICK           2.0000E-02        5.0000E-02        8.0000E-02
                #2             2      SPRING          1.0000E-02        6.2500E-02        7.5000E-02
                #3             3      SPRING          5.0000E-02        1.2500E-01        1.5000E-01
            # internal, desvar, label, lower,   upper,  ???,  ???
            # (1, 1, b'THICK   ', 0.019999, 0.07999999, -0.5, 0.0)
            # (2, 2, b'SPRING  ', 0.009999, 0.07500000, -0.5, 0.0)
            # (3, 3, b'SPRING  ', 0.050000, 0.15000000, -0.5, 0.0)

            ##        id       label   xinit   xlb   xub delxv
            #DESVAR  1       THICK   .05     0.02    .08
            #DESVAR  2       SPRING  .0625   .01     .075
            #DESVAR  3       SPRING  .125    .05     .15

            # C:\NASA\m4\formats\git\examples\move_tpl\betaadj.op2
            #---------------------------------------------------------------------------------------------------------
                #INTERNAL       DESVAR                         LOWER                               UPPER
                    #ID            ID          LABEL            BOUND             VALUE             BOUND
            #---------------------------------------------------------------------------------------------------------
                #11            20      BETA            1.0000E-03        8.0000E-01        1.0000E+20
            # int, des, label,     lower,     upper,     ???,  ???
            #(11, 20, b'BETA    ', 0.0010000, 1.000e+20, 0.200, 0.0)
            #        id       label   xinit   xlb   xub delxv
            #desvar  20       beta    0.8     0.001	xub	0.20
            desvar = structi.unpack(data)
            #internal_id, desvar_id, label, lower, upper, delxv, dunno = desvar
            #print(desvar)
            #assert np.allclose(desvar[5], -0.5), desvar  # -0.5 is the default
            assert np.allclose(desvar[6], 0.0), desvar
            desvars.append(desvar)
            itable -= 1
            markers = self.read_3_markers([itable, 1, 0])

        self.op2.op2_results.responses.desvars = Desvars(desvars)
        #if self.read_mode == 2:
            #self.log.warning('DESTAB results were read, but not saved')
        markers = self.read_markers([0])


    def _read_cmodext(self):
        r"""
        fails if a streaming block???:
         - nx_spike\mnf16_0.op2

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)

        size = self.size
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()
        out = unpack(mapfmt(b'7i', size), data)
        self.log.debug(str(out))
        # (101, 500, 37232, 2, 1, 18355, 158)
        #self.show_data(data, types='sqd')

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        marker = -2
        markers = self.read_markers([marker, 1, 0])

        data = self._read_record()
        if self.size == 4:
            unused_table_name, oneseventy_a, oneseventy_b = unpack('8sii', data)
        else:
            table_name, oneseventy_a, oneseventy_b = unpack('16sqq', data)
            unused_table_name = reshape_bytes_block(table_name)
        assert oneseventy_a == 170, oneseventy_a
        assert oneseventy_b == 170, oneseventy_b
        #print('170*4 =', 170*4)
        #self.show_data(data)
        marker -= 1
        n = op2.n
        while marker < 0:
            #marker = self._read_cmodext_helper(marker) # -3
            marker, is_done = self._read_cmodext_helper(marker)
            #print(f'--end marker={marker}')
            #print(f'--end marker2={marker2}')
            if is_done:
                break
            marker2 = self.get_nmarkers(1, rewind=True)[0]
            n = op2.n
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        #marker = self._read_cmodext_helper(marker); print(marker)
        print(f'table end; {marker2-1}')
        #marker = self._read_cmodext_helper(marker, debug=True)
        #op2.show_ndata(200)
        #sss

    def _read_cmodext_helper(self, marker_orig, debug=False) -> tuple[int, bool]:
        r"""

        64-bit:
          C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\mbdrecvr_c_0.op2
          C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cntlmtl05_0.op2
        """
        op2 = self.op2
        marker = marker_orig
        #markers = self.read_nmarkers([marker, 1, 1]) # -3
        debug = False
        if debug:
            op2.show_ndata(100)
        markers = self.get_nmarkers(3, rewind=False)
        assert markers == [marker_orig, 1, 1], markers
        #print('markers =', markers)

        #marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        val_old = 0
        if debug:
            print('-----------------------------')
        #i = 0
        #icheck = 7
        #factor = self.factor
        size = self.size
        if size == 4:
            expected_marker = 3
            sint = 'i'
            ifs = 'ifs'
            #dq = 'if'
            structi = Struct(b'i')
            struct_qd = Struct(b'i f')
            struct_q2d = Struct(b'i 2f')
            struct_q3d = Struct(b'i 3f')
            struct_d4q = Struct(b'i 4f')
        else:
            expected_marker = 3
            sint = 'q'
            ifs = 'qds'
            #dq = 'dq'
            structi = Struct(b'q')
            struct_qd = Struct(b'q d')
            #struct_dq = Struct(b'd q')
            struct_q2d = Struct(b'q 2d')
            struct_q3d = Struct(b'q 3d')
            struct_d4q = Struct(b'd 4q')

        marker = self.get_nmarkers(1, rewind=True, macro_rewind=True)[0]
        #print('AAAAA', marker)
        if marker is None:
            asdf
        elif marker == 3:
            pass
        elif marker < 0:
            self.log.warning('finished CMODEXT :)!')
            done = True
            return None, done
        else:
            raise RuntimeError(marker)

        i = 0
        while 1:
            #print('------------------------------')
            #print('i = %i' % i)
            marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
            if marker != expected_marker:
                self.log.info(f'CMODEXT i={i} marker={marker}')

            if self.size == 4:
                assert marker in [expected_marker], marker
            else:
                #if marker == -4:
                    #self.show(200, types='isq', endian=None, force=False)
                assert marker in [1, 2, 3, 10, 13, 16], marker
            data = self._read_block()
            ndata = len(data)
            nvalues = ndata // size
            if nvalues == 2:
                #self.show_data(data, types=dq)
                #print(ndata)
                #self.show_data(data, types=ifs)
                a, b = struct_qd.unpack(data)
                print(f' 2 {i}: {a:-4d} {b:.6g}')
            elif nvalues == 3:
                #self.show_data(data, types=dq)
                #print(ndata)
                #self.show_data(data, types=ifs)
                a, b, c = struct_q2d.unpack(data)
                print(f' 3 {i}: {a:-4d} {b:.6g}')
            elif nvalues == 4:
                #self.show_data(data, types=dq)
                #print(ndata)
                #self.show_data(data, types=ifs)
                a, b, c, d = struct_q3d.unpack(data)
                #print(f'4 {i}: {a:-4d} {b:10g} {c:10g} {d:10g}')
            elif nvalues == 5:
                #self.show_data(data, types=dq)
                #print(ndata)
                out = struct_d4q.unpack(data)
                print('5 {i}:', out)
            #elif ndata == 24 * factor:
                #self.show_data(data, types=dq)
                #print(ndata)
                #self.show_data(data, types=ifs)
                #a, b, c, d = struct_q2d.unpack(data)
            elif ndata == 11:
                # ints    = (697, 0, 0, 1072693248, -427255041, -1128310260, -1065349346, 1019092170, -458662240, 1015082339, -1959380180, 1018327314, 863132413, -1129941844, 477962639, 992946870, -935437031, 992276971, -1286826876, -1158808745, 1765409376, 986144275)
                # floats  = (9.767050296343975e-43, 0.0, 0.0, 1.875, -3.2255050747205135e+23, -0.023358367383480072, -4.001845359802246, 0.023207087069749832, -2.4994545586342025e+22, 0.015738194808363914, -3.509644095891333e-32, 0.02178243175148964, 5.642776912395675e-08, -0.02031930536031723, 8.375405145620162e-22, 0.0026728338561952114, -194932.390625, 0.002516860840842128, -4.76325254794574e-08, -0.0018156570149585605, 1.4054463629294865e+25, 0.001521053141914308)
                # doubles (float64) = (3.444e-321, 1.0, -4.350930636102116e-16, 4.1789445167220847e-16, 2.936419423485894e-17, 2.559298496168014e-16, -1.5581808467238972e-16, 1.2890303437976399e-23, 8.662650619686277e-24, -7.750115449925671e-25, 1.5100882399748836e-25)
                # long long (int64) = (697, 4607182418800017408, -4846055662573544705, 4376967544989290270, 4359745452588490400, 4373682512589110060, -4853063265498801411, 4264674333793526159, 4261797142378470681, -4977045659085663100, 4235457412028039776)

                #ints    = (697, 0, 0, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z)
                #floats  = (697, 0.0,
                           #0.0, 1.875,
                           #-3.2255050747205135e+23, -0.023358367383480072, -4.001845359802246, 0.023207087069749832, -2.4994545586342025e+22, 0.015738194808363914, -3.509644095891333e-32, 0.02178243175148964, 5.642776912395675e-08, -0.02031930536031723, 8.375405145620162e-22, 0.0026728338561952114, -194932.390625, 0.002516860840842128, -4.76325254794574e-08, -0.0018156570149585605, 1.4054463629294865e+25, 0.001521053141914308)
                #doubles (float64) = (697, 1.0, y, y, y, y, y, y, y, y, y)
                #long long (int64) = (697, 1.0, x, x, x, x, x, x, x, x, x)
                self.show_data(data, types='ifsqd')
                #self.show_data(data[4:], types='qd')
                raise RuntimeError(f'marker={marker} ndata={ndata}; nvalues={nvalues}')
            elif nvalues == 17:
                inti = structi.unpack(data[:size])
                floats = np.frombuffer(data, dtype=op2.fdtype8)[1:]
                print(inti, floats)
            else:
                self.show_data(data, types=ifs)
                raise RuntimeError(f'marker={marker} ndata={ndata}; nvalues={nvalues}')
                #sdf
                #ints = [
                    #1387, 0, -1574726656, -1076024976, 12958534, -1079706775, -1204798216, -1076795484, -674762419, -1125074255, -1234714250, 1025640630, 367990681, 1024314941, -1752243687, 1021642278,
                    #-3, 1017741311,
                    #-3, 1019878559,
                    #-2, 1019293439,
                    #-5, -1134034945,
                    #-3, 1017733119,
                    #-3, 1017241599]
                #floats = [
                    #1387, 0.0, -2.2168969812013176e-18, -1.7278270721435547, 1.815877379410093e-38, -1.2889224290847778, -4.201845149509609e-05, -1.6359753608703613, -439678385782784.0, -0.029385896399617195, -3.45342914442881e-06, 0.03955908864736557, 9.656856178702571e-26, 0.034620512276887894, -9.233609993079022e-25, 0.02795703336596489,
                    #nan, 0.02069091610610485,
                    #nan, 0.024671850726008415,
                    #nan, 0.023581979796290398,
                    #nan, -0.014160155318677425,
                    #nan, 0.02067565731704235,
                    #nan, 0.01976012997329235]
                #self.show_data(data, types='ifsqd')
                #ints = np.frombuffer(data, dtype='int32').tolist()
                #floats = np.frombuffer(data, dtype='float32').tolist()
                #print(ints)
                #print(floats)
                #np.frombuffer(data, dtype='ifsdq')
                print('------------------------------')
                self.show_data(data, types=ifs)
            #else:
                #self.show_data(data, types=ifs)

            val = unpack(sint, data[:size])[0]
            if debug:
                print('val=%s delta=%s' % (val, val - val_old))
                #self.show_data(data, types=ifs)
                #self.show_data(data[4:], types=ifs)
            assert len(data) > size
            #print('i=%s val=%s delta=%s' % (i, val, val - val_old))
            val_old = val

            marker2 = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
            #print(marker2)
            if marker2 < 0:
                self.log.warning(f'breaking marker2={marker2}')
                self.show(300, types=ifs)
                break
            #if marker2 == 696:
                #self.log.warning('breaking 696')
                #break
            i += 1
        if debug:
            print('----------------------------------------')

        marker = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
        #if debug:
        print('****marker  = %s' % marker)
        #print('****marker2 = %s' % marker2)
        #assert marker == 696, marker
        #data = self._read_block()
        #self.show_data(data)

        #marker = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
        assert marker == (marker_orig - 1), f'marker={marker} marker_orig={marker_orig}'

        if debug:
            op2.show_ndata(200)
        done = False
        return marker, done

        #data = self._read_record()
        #marker -= 1
        #op2.show_ndata(100)

        ##marker -= 1
        ##marker_end = op2.get_marker1(rewind=False)

    def read_cstm(self):
        """
        Reads the CSTM table, which defines the transform from global to basic.

        Returns 14-column matrix 2-d array of the CSTM data:
        ::
          [
           [ id1 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           [ id2 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           ...
          ]

        T is transformation from local to basic for the coordinate system.

        """
        op2 = self.op2
        is_geometry = op2.is_geometry
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        #print(self.show_data(data, types='ifsqd'))
        # 101, 466286, 15,      1, 1, 180, 0
        # 101, 466286, ncoords, 1, 1, 180, 0
        factor = self.factor
        #if self.size == 4:
            #idtype = 'int32'
            #fdtype = 'float32'
        #else:
            #idtype = 'int64'
            #fdtype = 'float64'
        assert len(data) == 28 * factor, len(data)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record() # CSTM
        #print(self.show_data(data, types='s'))
        assert len(data) == 8 * factor, len(data)

        self.read_3_markers([-3, 1, 0])

        coord_type_map = {
            1 : 'CORD2R',
            2 : '???',
            3 : 'CORD2S',
            5 : 'GMSURF',
            7 : '???',
            8 : '???',
        }
        #1. Coordinate system type:
        #- 0 = unknown (seriously?)
        #- 1 = rectangular
        #- 2 = cylindrical
        #- 3 = spherical
        #- 4 = convective coordinate system defined on a GMCURV+GMSURF pair
        #- 5 = convective coordinate system defined on a GMSURF
        #- 6 = convective coordinate system defined on a FEEDGE+FEFACE pair
        #- 7 = convective coordinate system defined on a FEFACE
        i = 0
        itable = -4

        blocks = []
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            data = self._read_record()
            blocks.append(data)
            self.read_markers([itable, 1, 0])
            itable -= 1
        markers = self.get_nmarkers(1, rewind=False)

        if not self.read_mode == 1: # or b'GEOM1' in op2.table_names:

        #if not is_geometry: # or self.read_mode == 1 or b'GEOM1' in op2.table_names:
            return
        nblocks = len(blocks)
        if nblocks == 1:
            # vectorized
            ints = np.frombuffer(blocks[0], dtype=op2.idtype8)
            floats = np.frombuffer(blocks[0], dtype=op2.fdtype8)
            #doubles = np.frombuffer(blocks[1], dtype='float64')
            nints = len(ints)

            assert nints % 14 == 0, 'nints=%s' % (nints)
            ncstm = get_table_size_from_ncolumns('CSTM', nints, 14)
            ints = ints.reshape(ncstm, 14)[:, :2]
            floats = floats.reshape(ncstm, 14)[:, 2:]
            #assert ncstm == 1, 'ncoords = %s' % ncstm
            #print(self.coords)
            #for i, unused_coord in enumerate(ints):
            cid = ints[i, 0]
            coord_type_int = ints[i, 1]
            if coord_type_int in coord_type_map:
                unused_coord_type = coord_type_map[coord_type_int]
            else:  # pragma: no cover
                msg = 'cid=%s coord_type_int=%s is not supported\n' % (cid, coord_type_int)
                if hasattr(self, 'coords'):
                    print(op2.coords)
                raise RuntimeError(msg)

            if not is_geometry:
                return
            for intsi, valuesi in zip(ints, floats):
                cid = intsi[0]
                coord_type_int = intsi[1]
                assert len(valuesi) == 12, valuesi
                #if coord_type_int in coord_type_map:
                    #unused_coord_type = coord_type_map[coord_type_int]
                #else:  # pragma: no cover
                    #msg = 'cid=%s coord_type_int=%s is not supported\n' % (cid, coord_type_int)
                    #if hasattr(self, 'coords'):
                        #print(op2.coords)
                    #raise RuntimeError(msg)

                origin = valuesi[:3]
                i = valuesi[3:6]
                unused_j = valuesi[6:9]
                k = valuesi[9:12]
                zaxis = origin + k
                xzplane = origin + i
                assert len(origin) == 3, origin
                assert len(zaxis) == 3, zaxis
                assert len(xzplane) == 3, xzplane
                if coord_type_int == 1:
                    coord = op2.add_cord2r(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                elif coord_type_int == 2:
                    coord = op2.add_cord2c(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                elif coord_type_int == 3:
                    coord = op2.add_cord2s(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                elif coord_type_int == 5:
                    #- 7 = convective coordinate system defined on a FEFACE
                    print('COORD_GMSURF', cid, origin, zaxis, xzplane)
                    coord = None
                elif coord_type_int == 7:
                    #- 7 = convective coordinate system defined on a FEFACE
                    print('COORD_FEFACE', cid, origin, zaxis, xzplane)
                    coord = None
                elif coord_type_int == 8:
                    #- 7 = convective coordinate system defined on a FEFACE
                    print('COORD_???', cid, origin, zaxis, xzplane)
                    coord = None
                else:  # pragma: no cover
                    raise NotImplementedError(f'coord_type_int={coord_type_int}')
                str(coord)
        elif nblocks == 2:
            # cstm style
            #block4 - 4 values  - cid, type, int_index, double_index
            #block5 - 12 values - ox, oy, oz, T11, T12, T13, T21, T22, T23, T31, T32, T33
            ints = np.frombuffer(blocks[0], dtype=self.op2.idtype8)
            doubles = np.frombuffer(blocks[1], dtype='float64')
            nints = len(ints)
            ndoubles = len(doubles)
            ncoords = nints // 4
            ints = ints.reshape(ncoords, 4)
            #doubles2 = doubles.reshape(ncoords, 12)

            # origins = doubles[:, 3]
            # zaxes = origins + doubles[:, 9:12]
            # xzplanes = origins + doubles[:, 6:9]

            #print('ints =', ints.tolist())
            if ncoords == ndoubles // 12:
                values = doubles.reshape(ncoords, 12)
                #print('doubles =', doubles.tolist())
            else:
                values = np.frombuffer(blocks[1], dtype='float32').reshape(ncoords // 4, 12)
                #print('floats =', floats.tolist())

            op2.op2_results.cstm.data = np.concatenate([ints, values], axis=1)
            if not is_geometry:
                return

            for intsi, valuesi in zip(ints, values):
                cid, cid_type, unused_int_index, unused_double_index = intsi
                assert len(valuesi) == 12, valuesi

                origin = valuesi[:3]
                i = valuesi[3:6]
                unused_j = valuesi[6:9]
                k = valuesi[9:12]
                zaxis = origin + k
                xzplane = origin + i
                assert len(origin) == 3, origin
                assert len(zaxis) == 3, zaxis
                assert len(xzplane) == 3, xzplane
                if cid_type == 1:
                    coord = op2.add_cord2r(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                elif cid_type == 2:
                    coord = op2.add_cord2c(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                elif cid_type == 3:
                    coord = op2.add_cord2s(cid, rid=0,
                                           origin=origin, zaxis=zaxis, xzplane=xzplane,
                                           comment='')
                else:  # pragma: no cover
                    raise NotImplementedError(f'cid_type={cid_type}')
                str(coord)
        else:  # pragma: no cover
            raise NotImplementedError(f'nCSTM blocks={nblocks} (not 1 or 2)')

        #print(self.op2.coords)

        #while i < len(ints):
            #break
            #cid = ints[i]
            #coord_type_int = ints[i + 1]

            #if coord_type_int in coord_type_map:
                #unused_coord_type = coord_type_map[coord_type_int]
            #else:  # pragma: no cover
                #msg = 'cid=%s coord_type_int=%s is not supported\n' % (cid, coord_type_int)
                #if hasattr(self, 'coords'):
                    #print(op2.coords)
                #raise RuntimeError(msg)

            #if coord_type_int in [1, 2, 3]: # rectangular, cylindrical, spherical
                #translations_cosines = floats[i+3:i+14]
                #print(len(translations_cosines))
                #i += 14
            #else:
                #raise NotImplementedError('coord_type_int = %i' % coord_type_int)
        #if 0:
            ## vectorized
            #assert nints % 14 == 0, 'nints=%s' % (nints)
            #ncstm = get_table_size_from_ncolumns('CSTM', nints, 14)
            #ints = ints.reshape(ncstm, 14)[:, :2]
            #floats = floats.reshape(ncstm, 14)[:, 2:]
            ##assert ncstm == 1, 'ncoords = %s' % ncstm
            ##print(self.coords)
            #for i, unused_coord in enumerate(ints):
                #cid = ints[i, 0]
                #coord_type_int = ints[i, 1]
            #if coord_type_int in coord_type_map:
                #unused_coord_type = coord_type_map[coord_type_int]
            #else:  # pragma: no cover
                #msg = 'cid=%s coord_type_int=%s is not supported\n' % (cid, coord_type_int)
                #if hasattr(self, 'coords'):
                    #print(op2.coords)
                #raise RuntimeError(msg)

            #print(self.coords)
            #print('cid = ', cid)
            #print('coord_type = ', coord_type)
            #print('myints =', ints)
            #print('floats =', floats)

    def read_trmbd(self):
        """
        Reads the TRMBD table, which defines the transform from material (deformed) to basic.
        """
        op2 = self.op2
        is_geometry = op2.is_geometry
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        #print(self.show_data(data, types='ifsqd'))
        # 101, 466286, 15,      1, 1, 180, 0
        # 101, 466286, ncoords, 1, 1, 180, 0
        header_int = np.frombuffer(data, dtype=op2.idtype8)
        #header_float = np.frombuffer(data, dtype=op2.fdtype8)
        ncoords = header_int[2]
        #print('data =', header_int)
        #print('ncoords =', ncoords)
        factor = self.factor
        #if self.size == 4:
            #idtype = 'int32'
            #fdtype = 'float32'
        #else:
            #idtype = 'int64'
            #fdtype = 'float64'
        assert len(data) == 28 * factor, len(data)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record() # CSTM
        # print(self.show_data(data, types='s'))
        # assert len(data) == 8 * factor, len(data)

        self.read_3_markers([-3, 1, 0])

        itable = -4

        blocks = []
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            data = self._read_record()
            blocks.append(data)
            self.read_markers([itable, 1, 0])
            itable -= 1
        markers = self.get_nmarkers(1, rewind=False)

        if self.read_mode == 1:
            return
        nblocks = len(blocks)

        assert(nblocks % 2 == 0)  # Should be even, first block is IDENT, second block is DATA

        # Get time steps per subcase
        time_steps = {}  # {subcase_id: {t1, t2, ...}
        for i in range(0, nblocks, 2):
            identifiers_int = np.frombuffer(blocks[i], dtype=op2.idtype8)
            identifiers_float = np.frombuffer(blocks[i], dtype=op2.fdtype8)

            acode = identifiers_int[0]
            tcode = identifiers_int[1]
            #eltype = identifiers_int[2]
            subcase = identifiers_int[3]
            time_step = identifiers_float[4]

            if subcase not in time_steps:
                time_steps[subcase] = {time_step, }
            else:
                time_steps[subcase].add(time_step)

        # Create time step to tid mapper per subcase
        for subcase, tsteps in time_steps.items():
            time_steps[subcase] = np.array(list(tsteps))  # convert set to numpy array
            tstep_to_index_mapper = {time_steps[subcase][x]: x for x in range(0, time_steps[subcase].shape[0])}

        trmbd = op2.trmbd
        # Read data
        #print('nblocks =', nblocks)
        for i in range(0, nblocks, 2):
            identifiers_int = np.frombuffer(blocks[i], dtype=op2.idtype8)
            identifiers_float = np.frombuffer(blocks[i], dtype=op2.fdtype8)

            acode = identifiers_int[0]
            tcode = identifiers_int[1]
            element_type = identifiers_int[2]
            subcase = identifiers_int[3]
            time_step = identifiers_float[4]

            op2.data_code = {}
            op2.subtable_name = ''
            approach_code = acode
            tCode = tcode
            int3 = element_type
            isubcase = subcase
            op2.element_type = element_type
            op2._set_approach_code(approach_code, tCode, int3, isubcase)
            #------------------------------------------------------------------------
            int_data = np.frombuffer(blocks[i+1], dtype=op2.idtype8)
            float_data = np.frombuffer(blocks[i+1], dtype=op2.fdtype8)

            #grids = int_data[:, np.array([1, 5, 9, 13, 17, 21])]
            if element_type == 363:
                element = 'CROD'
                nnodes = 2
            elif element_type == 347:
                nnodes = 2  # ???
                element = 'CBAR'
            elif element_type == 348:
                element = 'CBEAM'
                nnodes = 2  # ???

            elif element_type == 301:
                element = 'CPENTA'
                nnodes = 6
                #grids = int_data[:, np.array([1, 5, 9, 13, 17, 21])]
                #eulers_x = float_data[:, np.array([2, 6, 10, 14, 18, 22])]
                #eulers_y = float_data[:, np.array([3, 7, 11, 15, 19, 23])]
                #eulers_z = float_data[:, np.array([4, 8, 12, 16, 20, 24])]
            elif element_type == 302:
                element = 'CTETRA'
                nnodes = 4
                #grids = int_data[:, np.array([1, 5, 9, 13])]
                #eulers_x = float_data[:, np.array([2, 6, 10, 14])]
                #eulers_y = float_data[:, np.array([3, 7, 11, 15])]
                #eulers_z = float_data[:, np.array([4, 8, 12, 16])]
            elif element_type == 303:
                element = 'CPYRAM'
                nnodes = 5

            elif element_type == 343:
                element = 'CTRIA6'
                nnodes = 3
            elif element_type == 344:
                element = 'CQUAD8'
                nnodes = 4
            elif element_type == 346:
                element = 'CQUADR'
                nnodes = 4

            elif element_type == 300: # or eltype == 67:
                element = 'CHEXA8'
                nnodes = 8

            elif element_type == 345:
                nnodes = 3
                element = 'CTRIAR'
            elif element_type == 346:
                element = 'CQUADR'
                nnodes = 4

            elif element_type == 313:
                element = 'CQUADX4'
                nnodes = 4
            elif element_type == 315:
                element = 'CQUADX8'
                nnodes = 4

            elif element_type == 316:
                nnodes = 3
                element = 'CPLSTN3'
            elif element_type == 317:
                nnodes = 4
                element = 'CPLSTN4'
            elif element_type == 318:
                nnodes = 3
                element = 'CPLSTN6'
            elif element_type == 319:
                nnodes = 4
                element = 'CPLSTN8'

            #elif element_type == 323:
                #nnodes = 6
                #element = 'CPLSTS6'
            elif element_type == 321:
                nnodes = 4
                element = 'CPLSTS4'
            elif element_type == 323:
                nnodes = 4
                element = 'CPLSTS8'

            else:
                #print(int_data[:20])
                print(op2.code_information())
                raise NotImplementedError(f"Element type {element_type} is not implemented.\n")
            n_elements, int_data, float_data = reshape_trmbd(element, nnodes, int_data, float_data)

            eids = (int_data[:, 0] - op2.device_code) // 10
            #nelements = len(eids)
            #print('ndatai =', int_data.size)
            #print('nelements =', nelements)
            index = np.arange(1, 1+nnodes*4, 4)
            grids = int_data[:, index]
            eulers_x = float_data[:, index+1]
            eulers_y = float_data[:, index+2]
            eulers_z = float_data[:, index+3]

            if subcase not in trmbd.eulersx:
                trmbd.eulersx[subcase] = {}
                trmbd.eulersy[subcase] = {}
                trmbd.eulersz[subcase] = {}
            if element not in trmbd.eulersx[subcase]:
                #isubcase = time_steps[subcase]
                trmbd.eulersx[subcase][element] = np.empty([time_steps[subcase].shape[0], n_elements, eulers_x.shape[1]])
                trmbd.eulersy[subcase][element] = np.empty([time_steps[subcase].shape[0], n_elements, eulers_y.shape[1]])
                trmbd.eulersz[subcase][element] = np.empty([time_steps[subcase].shape[0], n_elements, eulers_z.shape[1]])
            if element not in trmbd.nodes:
                nodes = np.empty([n_elements, grids.shape[1] + 1], dtype='int32')
                nodes[:, 0] = eids
                nodes[:, 1:] = grids
                trmbd.nodes[element] = nodes

            itime = tstep_to_index_mapper[time_step]
            trmbd.eulersx[subcase][element][itime, :, :] = eulers_x
            trmbd.eulersy[subcase][element][itime, :, :] = eulers_y
            trmbd.eulersz[subcase][element][itime, :, :] = eulers_z

        return

    def read_trmbu(self):
        """
        Reads the TRMBU table, which defines the transform from material (undeformed) to basic.
        """
        op2 = self.op2
        #is_geometry = op2.is_geometry
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        # print(self.show_data(data, types='ifsqd'))
        # 101, 466286, 15,      1, 1, 180, 0
        # 101, 466286, ncoords, 1, 1, 180, 0
        factor = self.factor
        #if self.size == 4:
            #idtype = 'int32'
            #fdtype = 'float32'
        #else:
            #idtype = 'int64'
            #fdtype = 'float64'
        assert len(data) == 28 * factor, len(data)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record() # CSTM
        # print(self.show_data(data, types='s'))
        # assert len(data) == 8 * factor, len(data)

        self.read_3_markers([-3, 1, 0])
        itable = -4

        blocks = []
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            data = self._read_record()
            blocks.append(data)
            self.read_markers([itable, 1, 0])
            itable -= 1
        markers = self.get_nmarkers(1, rewind=False)

        if self.read_mode == 1:
            return
        nblocks = len(blocks)

        assert(nblocks % 2 == 0)  # Should be even, first block is IDENT, second block is DATA

        # Get time steps per subcase
        time_steps = {}  # {subcase_id: {t1, t2, ...}
        for i in range(0, nblocks, 2):
            identifiers_int = np.frombuffer(blocks[i], dtype=op2.idtype8)
            identifiers_float = np.frombuffer(blocks[i], dtype=op2.fdtype8)

            acode = identifiers_int[0]
            tcode = identifiers_int[1]
            #eltype = identifiers_int[2]
            subcase = identifiers_int[3]
            time_step = identifiers_float[4]

            if subcase not in time_steps:
                time_steps[subcase] = {time_step, }
            else:
                time_steps[subcase].add(time_step)

        # Create time step to tid mapper per subcase
        for subcase, tsteps in time_steps.items():
            time_steps[subcase] = np.array(list(tsteps))  # convert set to numpy array
            tstep_to_index_mapper = {time_steps[subcase][x]: x for x in range(0, time_steps[subcase].shape[0])}

        # Read data
        trmbu = op2.trmbu
        for i in range(0, nblocks, 2):
            identifiers_int = np.frombuffer(blocks[i], dtype=op2.idtype8)
            identifiers_float = np.frombuffer(blocks[i], dtype=op2.fdtype8)

            acode = identifiers_int[0]
            tcode = identifiers_int[1]
            element_type = identifiers_int[2]
            subcase = identifiers_int[3]
            time_step = identifiers_float[4]

            op2.data_code = {}
            op2.subtable_name = ''
            approach_code = acode
            tCode = tcode
            int3 = element_type
            isubcase = subcase
            op2.element_type = element_type
            op2._set_approach_code(approach_code, tCode, int3, isubcase)
            #------------------------------------------------------------------------

            int_data = np.frombuffer(blocks[i+1], dtype=op2.idtype8)
            float_data = np.frombuffer(blocks[i+1], dtype=op2.fdtype8)
            numwide = 4
            if element_type == 301: # or eltype == 67:
                element = 'CPENTA'
            elif element_type == 302: # or eltype == 67:
                element = 'CTETRA'
            elif element_type == 303: # or eltype == 67:
                element = 'CPYRAM'

            elif element_type == 300: # or eltype == 67:
                element = 'CHEXA8'

                numwide = 4
                n_elements = int_data.shape[0] // numwide  # elid + 3 euler angles
                #int_data = int_data.reshape(n_elements, numwide)
                #float_data = float_data.reshape(n_elements, numwide)

                #eids = (int_data[:, 0] - 2) // 10
                # grids = int_data[:, np.array([1, 5, 9, 13, 17, 21, 25, 29])]

                #eulers = float_data[:, np.array([1, 2, 3])]
                # eulers_y = float_data[:, np.array([3, 7, 11, 15, 19, 23, 27, 31])]
                # eulers_z = float_data[:, np.array([4, 8, 12, 16, 20, 24, 28, 32])]

            elif element_type == 363:
                numwide = 4
                element = 'CROD'
            elif element_type == 347:
                numwide = 4
                element = 'CBAR'
            elif element_type == 348:
                numwide = 4
                element = 'CBEAM'

            elif element_type == 345:
                element = 'CTRIAR'
                numwide = 4
            elif element_type == 346:
                element = 'CQUADR'
                numwide = 4

            elif element_type == 343:
                element = 'CTRIA6'
                numwide = 4
            elif element_type == 344:
                element = 'CQUAD8'
                numwide = 4

            elif element_type == 316:
                numwide = 4
                element = 'CPLSTN3'
            elif element_type == 317:
                numwide = 4
                element = 'CPLSTN4'
            elif element_type == 318:
                element = 'CPLSTN6'
            elif element_type == 319:
                element = 'CPLSTN8'

            #elif element_type == 316:
                #numwide = 4
                #element = 'CPLSTS3'
            elif element_type == 321:
                numwide = 4
                element = 'CPLSTS4'
            #elif element_type == 318:
                #nnodes = 6
                #element = 'CPLSTS6'
            elif element_type == 323:
                numwide = 4
                element = 'CPLSTS8'

            elif element_type == 313:
                element = 'CQUADX4'
                numwide = 4
            elif element_type == 315:
                element = 'CQUADX8'
                numwide = 4
            else:
                print(int_data[:10])
                print(op2.code_information())
                raise NotImplementedError(f"Element type {element_type} is not implemented.\n")

            assert numwide == 4, numwide
            n_elements = int_data.shape[0] // numwide  # element_id + 3 euler angles
            int_data = int_data.reshape(n_elements, numwide)
            float_data = float_data.reshape(n_elements, numwide)

            eids = (int_data[:, 0] - op2.device_code) // 10
            eulers = float_data[:, 1:]

            if subcase not in trmbu.eulers:
                trmbu.eulers[subcase] = {}
            if element not in trmbu.eulers[subcase]:
                trmbu.eulers[subcase][element] = np.empty([time_steps[subcase].shape[0], n_elements, 4])

            itime = tstep_to_index_mapper[time_step]
            trmbu.eulers[subcase][element][itime, :, 0] = eids
            trmbu.eulers[subcase][element][itime, :, 1:] = eulers

        return

    def read_obc1(self):
        op2 = self.op2
        unused_table_name = self._read_table_name(rewind=False)

        #(102, 0, 0, 2, 1, 0, 0)
        self.read_markers([-1])
        data = self._read_record()

        self.read_3_markers([-2, 1, 0])

        table_names = ['OBC1', 'OBC1X', 'OBG1']
        subtable_name = self._read_subtable_name(table_names)

        self.read_3_markers([-3, 1, 0])

        op2.subtable_name = subtable_name
        op2.data_code = {
            'subtable_name': subtable_name,
        }
        #print(len(data))
        #print(op2.data_code)
        #str(op2.code_information())

        itable = -4

        if 0:
            import struct
            markers = self.get_nmarkers(1, rewind=True)
            while markers[0] != 0:
                self.read_3_markers([itable, 1, 0])

                try:
                    self.read_3_markers([itable - 1, 1, 0])
                except struct.error:
                    data = self._read_record()
                    print(itable, len(data))
                itable -= 1
        else:
            while 1:
                #print(itable)
                data = self._read_record() # 584
                ndata = len(data)
                #print('ndata A =', ndata)
                op2._op2_readers.reader_oqg._read_obc1_3(data, ndata)
                #print('---------------')

                self.read_3_markers([itable, 1, 0])  # -4
                itable -= 1

                marker = self.get_nmarkers(1, rewind=True)[0]
                #print('marker A1 =', itable, marker)
                if marker == itable:
                    self.read_3_markers([itable, 1, 0]) # -5
                    break
                #self.show(100, types='isdq')
                #self.show(100, types='isf')

                data = self._read_record()
                #print('ndata B =', len(data))

                #self.show(100, types='isdq')
                #print('---------------')

                marker = self.get_nmarkers(1, rewind=True)[0]
                #print('marker B =', marker, itable)
                if marker != itable:
                    aaa
                self.read_3_markers([itable, 1, 0]) # -5
                itable -= 1

                marker = self.get_nmarkers(1, rewind=True)[0]
                #print('marker C =', marker)
                if marker == 0:
                    break
                #print('---------------')

                #self.show(100, types='isdq')
            marker = self.get_nmarkers(1, rewind=False)[0]
            assert marker == 0, marker

        #(4, -5, 4, 4, 1, 4, 4, 0, 4, 4, 0, 4, 4, 0, 4)
        #data = self._read_record()

        #self.show_data(data)
        #self.show(2000)
        #aaa

    def _read_subtable_name(self, table_names: list[str]):
        data, ndata = self._read_record_ndata()
        if ndata == 8:
            table_name2, = self.op2.struct_8s.unpack(data)
            utable_name2 = table_name2.decode('utf-8').strip()
        elif self.size == 8 and ndata == 56:
            table_name2, day, month, year, onea, oneb = Struct(b'16s 5q').unpack(data)
            utable_name2 = table_name2.decode('utf-8').strip()
            #self.show_data(data[16:], types='qsd')
        else:
            self.show_data(data)
            raise NotImplementedError(data)
        assert utable_name2 in table_names, utable_name2
        return utable_name2

    def read_qualinfo(self):
        r"""
        Reads the QUALINFO table

        -100001 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;PVALID=0;DESINC=0;DISCRETE=FALSE;MASSID=0;ARBMID=0;PARTNAME=' ';TRIMID=0;MODULE=0)
        -100000 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;PVALID=0;DESINC=0;ARBMID=0;PARTNAME=' ';DISCRETE=FALSE;TRIMID=0)
         -99999 (AUXMID=0;AFPMID=0;HIGHQUAL=0;PVALID=0;DESINC=0;PRESEQP=TRUE;ARBMID=0;PARTNAME=' ';TRIMID=0;FLXBDYID=0;DFPHASE=' ')
         -99998 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;DESINC=0;DISCRETE=FALSE;ARBMID=0;MASSID=0;PARTNAME=' ';TRIMID=0)
           1431 (HIGHQUAL=0;AUXMID=0;AFPMID=0;DESINC=0;ARBMID=0;PARTNAME=' ';TRIMID=0;FLXBDYID=0)
           1459 (PEID=0;DESITER=0;PVALID=0;NL99=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;AUXMID=0;DESINC=0;DISCRETE=FALSE;MASSID=0;PARTNAME=' ';MODULE=0)
           1461 (PEID=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;AUXMID=0;DESINC=0;PARTNAME=' ';MODULE=0)
           1541 (SEID=0;PEID=0;MTEMP=0;DESITER=0;PVALID=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;P2G=' ';K2GG=' ';M2GG=' ';DELTA=FALSE;AUXMID=0;BNDSHP=FALSE;ADJOINT=FALSE;DESINC=0;DISCRETE=FALSE;CASEF06=' ';ISOLAPP=1;SUBCID=0;OSUBID=1;STEPID=0;RGYRO=0;PARTNAME=' ';SSTEPID=0;MODULE=0)

        Word Name Type Description
        1 NAME(2) CHAR4 Datablock Name

        Word Name Type Description
        1 DBKEY   I     database KEY associated with qualifiers
        2 QLEN(C) I     length in words of qualifiers string
        3 QUALSTR CHAR4 Qualifier information string
        Word 3 repeats QLEN times

        Word Name Type Description
        1 FUNIT   I Fortran unit op2 file was written to
        2 NUMKEYS I Number of keys
        3 BIT(5)  I ,{

        """
        # we read the table on the first pass, so if we ever see a
        # 64-bit table, the error message makes a bit more sense
        self.op2.set_as_msc()
        read_record_ndata = self.get_skip_read_record_ndata()

        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])

        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')

        # (301, 1, 8, 0, 0, 0, 0)
        # (???, ?, n, ?, ?, ?, ?)
        data, ndata = read_record_ndata()
        assert ndata == 28, self.show_data(data) # 7*4

        self.read_3_markers([-2, 1, 0])

        # QUALINFO
        data, ndata = read_record_ndata()
        assert ndata == 8, self.show_data(data)

        itable = -3
        qual_strings = {}
        while 1:
            self.read_3_markers([itable, 1, 0])
            stop_marker = self.get_marker1(rewind=True)
            if stop_marker == 0:
                break

            data, ndata = read_record_ndata()
            if self.read_mode == 1:
                db_key, qlen = unpack(self._endian + b'2i', data[:8])
                fmt = self._endian + b'%is' % (ndata - 8)
                qual_str = unpack(fmt, data[8:])[0].decode('latin1')
                qual = QualInfo.from_str(db_key, qual_str.rstrip(), self.log)
                qual_strings[db_key] = qual_str
            itable -= 1
        stop_marker = self.get_marker1(rewind=False)

    def read_extdb(self):
        r"""
        fails if a streaming block:
         - nx_spike\extse04c_0.op2

        DB NAME   op2 Name  Description
        =======   ========  ===========
        XSOP2DIR  XSOP2DIR  Table of contents of the op2 file.
        GEOM1X    GEOM1X    GRID point geometry.
        GEOM2X    GEOM2X    SPOINTs.
        GEOM4X    GEOM4X    ASET/ASET1 and QSET/QSET1 entries.
        EXTDB     MATK      Boundary stiffness.
        EXTDB     MATM      Boundary mass.
        EXTDB     MATP      Boundary loads.
        EXTDB     MATV      Boundary fluid-structure partitioning vector.
        EXTDB     TUG1      Displacement OTM table.
        EXTDB     MUG1      Displacement OTM matrix.
        EXTDB     TES1      Stress OTM table.
        EXTDB     MES1      Stress OTM matrix.
        EXTDB     TEF1      Force OTM table.
        EXTDB     MEF1      Force OTM matrix.
        EXTDB     MUG1B     Displacement OTM matrix in basic coordinates system.

        https://help.hexagonmi.com/bundle/MSC_Nastran_2022.2/page/Nastran_Combined_Book/Superelements_User_s_Guide/superOTM/TOC.OUTPUT2.Files.xhtml
        """
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\extse04c_cnv1_0.op2
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        self.read_3_markers([-2, 1, 0])
        marker = -3
        if self.read_mode == 1 or op2.make_geom is False:
            if self.read_mode == 1 and op2.make_geom is False:
                self.log.warning('reading the EXTRN tables requires the read_op2_geom')
            data = self._skip_record()
            while 1:
                #print('====================')
                #print(f'***reading {marker}')
                try:
                    self.read_markers([marker, 1])
                except FortranMarkerError:
                    #op2.show_ndata(100)
                    raise
                nfields = self.get_marker1(rewind=True)
                if nfields == 0:
                    break
                data, ndata = self._read_record_ndata()
            return

        # drop XSOP2DIR and PVT0
        iextdb = op2.table_count[b'EXTDB'] + 0
        xsop2dir_name = self.xsop2dir_names[iextdb]
        if self.read_mode == 2:
            data, ndata = self._read_record_ndata()
            name = ''
            name1 = ''
            name2 = ''
            if self.size == 4:
                if ndata == 8:
                    #self.show_data(data, types='ifs', endian=None)
                    name, = Struct(self._endian + b'8s').unpack(data)
                    name = name.decode('latin1').rstrip()
                    self.log.info(f'A: name={name!r} -> {xsop2dir_name!r}')
                    #print(name, 8)
                elif ndata == 16:
                    name, int1, int2 = Struct(self._endian + b'8s 2i').unpack(data)
                    name = name.decode('latin1').rstrip()
                    #name = name.decode(self._encoding)
                    #print(name, int1, int2, 16)
                    self.log.info(f'B: name={name!r} -> {xsop2dir_name!r} int1={int1} int2={int2}')
                elif ndata == 28:
                    #self.show_data(data)
                    name1, int1, name2, int2 = Struct(self._endian + b'8s i 12s i').unpack(data)
                    name1 = name1.decode('latin1').rstrip()
                    name2 = name2.decode('latin1').rstrip()
                    self.log.info(f'C DMI: name1={name1!r} -> {xsop2dir_name!r} int1={int1} name2={name2!r} int2={int2}')
                    # m.add_dmi(name, form, tin, tout, nrows, ncols, GCj, GCi, Real, Complex=None, comment='')
                    #print(name1, int1, name2, int2, 28)
                else:
                    self.show_data(data, types='ifs')
                    raise NotImplementedError(ndata)

            else:
                if ndata == 16:
                    name, = Struct(self._endian + b'16s').unpack(data)
                    name = reshape_bytes_block(name).decode('latin1').strip()
                    self.log.info(f'A: name={name!r} -> {xsop2dir_name}')
                elif ndata == 32:
                    name, int1, int2 = Struct(self._endian + b'16s 2q').unpack(data)
                    name = reshape_bytes_block(name).decode('latin1').strip()
                    self.log.info(f'B: name={name!r} -> {xsop2dir_name} int1={int1} int2={int2}')
                elif ndata == 56:
                    # (PHIP, 7, REAL, 1)
                    #self.show_data(data, types='ifsdq')
                    name1, int1, name2, int2 = Struct(self._endian + b'16s q 24s q').unpack(data)
                    name1 = reshape_bytes_block(name1).decode('latin1').strip()
                    name2 = reshape_bytes_block(name2).decode('latin1').strip()
                    self.log.info(f'C: name1={name1!r} -> {xsop2dir_name} int1={int1} name2={name2!r} int2={int2}')
                else:
                    self.show_data(data, types='ifsdq')
                    raise NotImplementedError(ndata)

        #self.show(200)
        while 1:
            #print('====================')
            #print(f'***reading {marker}')
            try:
                self.read_markers([marker, 1])
            except FortranMarkerError:
                #op2.show_ndata(100)
                raise

            nfields1 = self.get_marker1(rewind=True)
            if nfields1 == 0:
                nfields1 = self.get_marker1(rewind=False)
            elif nfields1 == 1:
                #data, ndata = self._read_record_ndata()
                nfields1 = self.read_markers([1])

                nfields_test = self.get_marker1(rewind=True)
                while nfields_test > 0:
                    nfields = self.get_marker1(rewind=False)
                    block = self._read_block()
                    nblock = len(block)
                    ndouble = (nblock - 4) // 8
                    fmt = mapfmt(self._endian + b'i%dd' % (ndouble), self.size)
                    out = Struct(fmt).unpack(block)
                    #print(out, nblock)
                    nfields_test = self.get_marker1(rewind=True)
                #print('-------')
                #print(f'end of marker={marker}')
                marker -= 1
                #marker = self.get_marker1(rewind=True)
                continue
            else:
                raise RuntimeError('EXTDB error')

            #op2.show_ndata(100)
            nfields = self.get_marker1(rewind=True)
            #print('nfields =', nfields)
            if nfields == 0:
                #if self.read_mode == 2:
                    #self.log.warning('breaking...')
                #self.show(200)
                break
            # ----------------------------------------------------------------------

            data, ndata = self._read_record_ndata()
            if self.read_mode == 2:
                if name not in ['GEOM1', 'GEOM2', 'GEOM2X', 'IGEOM2X', 'GEOM4', 'EXTDB']:
                    if ndata != 12:
                        self.log.warning(f'--B; ndata={ndata}--')

                if self.size == 4:
                    if ndata == 12:
                        name, int1, int2 = Struct(self._endian + b'4s 2i').unpack(data)
                        # b'\xff\xff\x00\x00' 65535 25535 12 ???
                        #print(name, int1, int2, 12)
                        #self.show_data(data)

                    elif name == 'GEOM1':
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom1.geom1_map)
                    elif name in ['GEOM2', 'GEOM2X', 'IGEOM2X']:
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom2.geom2_map)

                    elif name == 'GEOM4':
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom4.geom4_map)

                    elif name == 'EXTDB':
                        _read_extdb_extdb(self, xsop2dir_name, data, self._endian, 'int32')
                    elif name1 in ['TES', 'PHIP']:
                        _read_extdb_phip(self, marker, data, self._endian, 'i', op2.idtype8)
                    elif name1 in ['TQP', 'TEF']:
                        _read_extdb_phip(self, marker, data, self._endian, 'i', op2.idtype8)
                    else:
                        raise NotImplementedError(f'EXTDB; name={name!r} name1={name1!r} name2={name2!r}')
                else:
                    if ndata == 24:
                        ints = Struct(self._endian + b'3q').unpack(data)
                        assert ints == (65535, 65535, 65535), ints
                    elif name == 'GEOM1':
                        # _read_extdb_geom1(self, data, self._endian)
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom1.geom1_map)
                    elif name in ['GEOM2', 'IGEOM2X']:
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom2.geom2_map)
                    elif name == 'GEOM4':
                        _read_extdb_geomx(self, data, self._endian, op2.reader_geom4.geom4_map)

                    elif name == 'EXTDB':
                        _read_extdb_extdb(self, xsop2dir_name, data, self._endian, 'int64')
                    elif name1 == 'PHIP':
                        _read_extdb_phip(self, marker, data, self._endian, 'q', op2.idtype8)
                    elif name1 == 'TES':
                        _read_extdb_phip(self, marker, data, self._endian, 'q', op2.idtype8)

                    else:
                        self.log.warning(f'EXTDB; name={name!r} name1={name1!r} ndata={ndata}')
                        raise RuntimeError(f'EXTDB; name={name!r} name1={name1!r} ndata={ndata}')
                        #self.show_data(data, types='sqd')
            marker -= 1
            #print('--------------------')
        unused_marker_end = self.get_marker1(rewind=False)
        #aa
        #if self.read_mode == 2:
            #self.log.warning('returning...')
        return

    def read_descyc(self):
        """reads the DESCYC table"""
        op2 = self.op2

        # TODO: I think this is used to handle optimization
        op2._count += 1

        #op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        fmt = mapfmt(self._endian + b'7i', self.size)
        unused_ints = Struct(fmt).unpack(data)
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        if self.size == 4:
            name, = Struct(self._endian + b'8s').unpack(data)
        else:
            name, = Struct(self._endian + b'16s').unpack(data)
            name = reshape_bytes_block(name)
        assert name == b'DESCYC  ', name

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        if self.size == 4:
            design_cycle, design_cycle_type_bytes = Struct(self._endian + b'i8s').unpack(data)
        else:
            design_cycle, design_cycle_type_bytes = Struct(self._endian + b'q16s').unpack(data)
            design_cycle_type_bytes = reshape_bytes_block(design_cycle_type_bytes)

        #Design cycle type; 'D' for discretized design cycle
        #Blank for continuous design cycle
        if design_cycle_type_bytes == b'        ':
            design_cycle_type_str = 'continuous'
        elif design_cycle_type_bytes == b'D       ':
            design_cycle_type_str = 'discrete'
        else:
            raise NotImplementedError(design_cycle_type_bytes)

        self.read_markers([-4, 1, 0, 0])
        #print('descyc')
        #print('  ints =', ints)
        #print('  name =', name)
        #print(f'  design_cycle={design_cycle} type={design_cycle_type_str!r} count={op2._count}')

    def read_dbcopt(self):
        """reads the DBCOPT table, which is a design variable history table"""
        #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cc577.op2
        # ints    = (101, 11, 10, 3, 4, 0, 0)
        # objective_function = [1175749.5, 711181.875, 369194.03125, 112453.1406, 229.4625, 50.9286, 50.8863, 50.8316, 50.7017, 49.7571, 49.3475]
        # approx = [0.0, 651307.3125, 388093.0625, 150222.453125, 5894.2626, 50.9092, 50.8834, 50.8252, 49.4544, 49.6455, 49.2533]
        # max_value_of_constraint = [nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan]
        # desvar_ids = [1, 2, 3]
        # cycle_1_values = [1.0, 1.0, 1.0]
        # cycle_n_values = [1.4, 0.6, 1.4,
        #                   1.96, 0.36, 1.96,
        #                   2.744, 0.216, 2.744,
        #                   3.8416, 0.1296, 3.8416, ...]

        # 1 NFEA I Number of finite element analyses
        # 2 NAOP I Number of optimization cycles w.r.t. approximate model
        # 3 NDV I Number of design variables
        # 4 NCC I Convergence criterion
        # 5 UNDEF(2 ) None
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        fmt = mapfmt(self._endian + b'7i', self.size)
        num, nopt, napprox, nvars, one, zeroa, zerob = Struct(fmt).unpack(data)
        assert num == 101, num
        assert zeroa == 0, zeroa
        assert zerob == 0, zerob

        # (101, 11, 10, 3, 4, 0, 0)
        #self.show_data(data)
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()

        if self.size == 4:
            name, = Struct(self._endian + b'8s').unpack(data)
        else:
            name, = Struct(self._endian + b'16s').unpack(data)
            name = reshape_bytes_block(name)
        assert name == b'DBCOPT  ', name

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        #ndata = len(data) // 4
        if self.size == 4:
            fdtype = 'float32'
            idtype = 'int32'
        else:
            fdtype = 'float64'
            idtype = 'int64'
        objective_function = np.frombuffer(data, dtype=fdtype) # .tolist()
        #print(f'  objective_function = {objective_function}; n={len(objective_function)}')
        assert len(objective_function) == nopt, f'len(objective_function)={len(objective_function)} nopt={nopt}'

        self.read_3_markers([-4, 1, 0])
        data = self._read_record()
        approx = np.frombuffer(data, dtype=fdtype).copy()# .tolist()
        napprox_actual = len(approx)
        if approx[0] == 0.0:
            approx[0] = np.nan
            napprox_actual -= 1
        #print(approx.tolist())
        #assert napprox_actual == napprox, f'napprox_actual={napprox_actual} napprox={napprox}'
        #print(f'  approx = {approx}; n={len(approx)}')

        self.read_3_markers([-5, 1, 0])
        data = self._read_record()
        max_value_of_constraint = np.frombuffer(data, dtype=fdtype).tolist()
        #print(f'  max_value_of_constraint = {max_va/lue_of_constraint}; n={len(max_value_of_constraint)}')


        self.read_3_markers([-6, 1, 0])
        data = self._read_record()
        desvar_ids = np.frombuffer(data, dtype=idtype).tolist()
        assert len(desvar_ids) == nvars, f'len(desvars)={len(desvars)} nvars={nvars}'

        self.read_3_markers([-7, 1, 0])
        data = self._read_record()
        cycle_1_values = np.frombuffer(data, dtype=fdtype).tolist()

        self.read_3_markers([-8, 1, 0])
        marker0 = self.get_marker1(rewind=True)
        if marker0 == 0:
            self.read_markers([0])
            return
        data = self._read_record()
        cycle_n_values = np.frombuffer(data, dtype=fdtype)
        cycle_n_values2 = cycle_n_values.reshape(len(cycle_n_values) // nvars, nvars)
        cycle_n_values3 = np.vstack([cycle_1_values, cycle_n_values2])

        approx_obj_constraint = np.vstack([approx, objective_function, max_value_of_constraint]).T
        #print(f'  desvar_ids = {desvar_ids}')
        #print(f'  approx_obj_constraint; {approx_obj_constraint.shape}:\n{approx_obj_constraint}')
        #print(f'  cycle_n_values {cycle_n_values3.shape}:\n{cycle_n_values3}')

        self.read_markers([-9, 1, 0, 0])
        #data = self._read_record()
        #self.show_data(data)
        #self.show_ndata(100)

    def read_dscmcol(self):
        """reads the DSCMCOL table, which defines the columns? for the DSCM2 table"""
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        #fmt = mapfmt(self._endian + b'7i', self.size)
        #num, ndesvars, one_zero, zeroa, zerob, zeroc, zerod = Struct(fmt).unpack(data)
        #print(num, ndesvars, one_zero, zeroa, zerob, zeroc, zerod)
        # (101, 3, 1, 0, 0, 0, 0)
        #self.show_data(data)
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        if self.size == 4:
            name, = Struct(self._endian + b'8s').unpack(data)
        else:
            name, = Struct(self._endian + b'16s').unpack(data)
            name = reshape_bytes_block(name)
        assert name == b'DSCMCOL ', name

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()

        responses = {}
        if self.read_mode == 2:
            ints = np.frombuffer(data, dtype=op2.idtype8)
            floats = np.frombuffer(data, dtype=op2.fdtype8)
            nresponses_dresp1 = len(ints) // 9
            dscmcol_dresp1(responses, nresponses_dresp1, ints, floats)

        #self.show_data(data[4*idata:])
        self.read_3_markers([-4, 1, 0])
        nfields = self.get_marker1(rewind=True)
        if nfields == 0:
            self._save_dscmcol_response(responses)
            self.read_markers([0])
            return

        data = self._read_record()
        if self.read_mode == 2:
            # read the DRESP2 columns
            ints = np.frombuffer(data, dtype=op2.idtype8)
            floats = np.frombuffer(data, dtype=op2.fdtype8)
            nresponses_dresp2 = len(ints) // 6
            dscmcol_dresp2(responses, nresponses_dresp2, ints, floats)

        self._save_dscmcol_response(responses)
        self.read_markers([-5, 1, 0, 0])

    def _save_dscmcol_response(self, responses):
        """saves the DSCMCOL dictionary"""
        if self.read_mode == 2:
            assert len(responses) > 0
        if responses:
            if self.op2.op2_results.responses.dscmcol is not None:
                self.log.warning('overwriting DSCMCOL')
            respi = DSCMCOL(responses)
            str(respi)
            self.op2.op2_results.responses.dscmcol = respi

    def read_fol(self):
        """
        Reads the FOL table
        Frequency response frequency output list

        tested by TestOP2.test_monpnt3

        +------+---------+-------+-----------------+
        | Word |  Name   | Type  |   Description   |
        +======+=========+=======+=================+
        |  1   | NAME(2) | CHAR4 | Data block name |
        +------+---------+-------+-----------------+
        |  3   |  FREQ   |  RS   |   Frequency     |
        +------+---------+-------+-----------------+
        | Word 3 repeats until End of Record       |
        +------------------------------------------+

        +------+----------+------+-----------------------------+
        | Word |  Name    | Type |   Description               |
        +======+==========+======+=============================+
        |  1   |  WORD1   |  I   | Number of frequencies       |
        +------+----------+------+-----------------------------+
        |  2   |  WORD2   |  I   | Frequency set record number |
        +------+----------+------+-----------------------------+
        |  3   |  WORD3   |  I   | Number of loads             |
        +------+----------+------+-----------------------------+
        |  4   | UNDEF(3) | None | Not used                    |
        +------+----------+------+-----------------------------+

        """
        op2 = self.op2
        #op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        ndata = len(data)
        subtable_name_raw, = op2.struct_8s.unpack(data[:8])
        subtable_name = subtable_name_raw.strip()
        assert subtable_name == b'FOL', 'subtable_name=%r' % subtable_name

        nfloats = (ndata - 8) // 4
        assert nfloats * 4 == (ndata - 8)
        fmt = self._endian + b'%if' % nfloats
        freqs = np.array(list(unpack(fmt, data[8:])), dtype='float32')

        if self.read_mode == 2:
            if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                msg = (
                    'Cannot overwrite op2._frequencies...\n'
                    'op2._frequencies = %s\n'
                    'new_freqs = %s\n' % (op2._frequencies, freqs))
                raise RuntimeError(msg)
            op2._frequencies = freqs
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, freqs]\n'  % (subtable_name_raw))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
                self.binary_debug.write('  freqs = %s' % freqs)
        self._read_subtables()

    def read_frl(self):
        """
        reads the FRL (Frequency Response List) table

        tested by TestOP2.test_op2_good_sine_01

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        fmt1 = mapfmt(self._endian + b'7i', self.size)
        idata = unpack(fmt1, data)
        assert idata[0] == 101, f'idata[0]={idata[0]}; idata={idata}'
        assert idata[1] in [1, 2, 3, 4], f'idata[1]={idata[1]}; idata={idata}'
        assert idata[2] == 0, f'idata[2]={idata[2]}; idata={idata}'
        assert idata[3] == 0, f'idata[3]={idata[3]}; idata={idata}'
        assert idata[4] == 0, f'idata[4]={idata[4]}; idata={idata}'
        assert idata[5] == 0, f'idata[5]={idata[5]}; idata={idata}'
        assert idata[6] == 0, f'idata[6]={idata[6]}; idata={idata}'
        #print(self.show_data(data))


        self.read_3_markers([-2, 1, 0])
        data = self._read_record()

        if len(data) == 12:
            subtable_name_raw, = op2.struct_8s.unpack(data[:8])
            subtable_name = subtable_name_raw.strip()
            assert subtable_name in [b'FRL', b'FRL0'], 'subtable_name=%r' % subtable_name
        elif len(data) == 16:
            #(FRL, 200, 201)
            subtable_name_raw, = op2.struct_8s.unpack(data[:8])
            subtable_name = subtable_name_raw.strip()
            assert subtable_name in [b'FRL', b'FRL0'], 'subtable_name=%r' % subtable_name
        elif len(data) == 20:
            #(FRL, 70, 71, 72)
            subtable_name_raw, = op2.struct_8s.unpack(data[:8])
            subtable_name = subtable_name_raw.strip()
            assert subtable_name in [b'FRL', b'FRL0'], 'subtable_name=%r' % subtable_name
        elif len(data) == 24:
            #(FRL, 71, 72, 73, 74)
            subtable_name_raw, = op2.struct_8s.unpack(data[:8])
            subtable_name = subtable_name_raw.strip()
            assert subtable_name in [b'FRL', b'FRL0'], 'subtable_name=%r' % subtable_name
        else:
            self.show_data(data, types='ifsd')
            raise RuntimeError('bad length...')

        self.read_3_markers([-3, 1, 0])
        isubtable = -3
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            if self.read_mode == 1:
                self._skip_record()
            else:
                data = self._read_record()
                #self.show_data(data)
                freqs = np.frombuffer(data, dtype=op2.fdtype).copy()
                #print('read_mode=%s itable=%s freqs=%s' % (self.read_mode, isubtable, freqs.tolist()))
                if isubtable == -3:
                    if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                        msg = (
                            'Cannot overwrite op2._frequencies...\n'
                            'op2._frequencies = %s\n'
                            'new_freqs = %s\n' % (op2._frequencies, freqs))
                        raise RuntimeError(msg)
                    op2._frequencies = freqs
                else:
                    #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\rtr_mfreq41kf.op2
                    if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                        msg = (
                            'Cannot overwrite op2._frequencies...\n'
                            'op2._frequencies = %s\n'
                            'new_freqs = %s\n' % (op2._frequencies, freqs))
                        self.log.warning(msg)

            isubtable -= 1
            self.read_markers([isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        del isubtable
        self.read_markers([0])

    def read_rst(self):
        """
        reads the RST table (restart file?)

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\k402rerun3d2ss.op2
        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        header_data = self._read_record()  # (101, 1, 0, 0, 0, 0, 0) - longs
        #self.show_data(header_data, types='qds', endian=None, force=False)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        #print(data[:16])

        table_name, = Struct('<16s').unpack(data[:16])
        #date
        #self.show_data(data[16:], types='qds', endian=None, force=False)

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data, types='qd', endian=None, force=False)


        self.read_3_markers([-4, 1, 0])
        data = self._read_record()
        #self.show_data(data, types='ifqds', endian=None, force=False)

        if 0:  # pragma: no cover
            #(1, 10, 4617315517961601024, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536, 2314885530818453536,
            # -1, -1)
            ndata = len(data) - 5 * 8
            #print(len(data))
            aint, bint, cfloat, spaces, m1a, m1b = unpack('<2qd ' + str(ndata) +'s' + ' 2q', data)
            #print((aint, bint, cfloat, spaces, m1a, m1b))
            spaces = spaces.decode('latin1').strip()
            print((aint, bint, cfloat, spaces, m1a, m1b))


        ints = np.frombuffer(data, op2.idtype8)
        #print('ints', ints)
        i1 = np.where(ints == -1)[0]
        print('i1', i1)
        i0 = i1[:-1] + 1
        i0 = np.hstack([[0], i0])

        print(f'i0={i0}')
        print(f'i1={i1}')

        # header = (101, 1, 0, 0, 0, 0, 0)
        # num = 402
        # data = [
        #    (1, 16, 0.8000000000000002, '', -1); nspaces=1024
        #    (2, 24, 1.2, '', -1); nspaces=1024
        # ]
        for i0i, i1i in zip(i0, i1):
            if i0i == i1i:
                continue
            datai = data[i0i*8:i1i*8+8]
            ndatai = len(datai) - 4 * 8
            fmt = '<2qd ' + str(ndatai) +'s' + ' q'
            #print(fmt)
            aint, bint, cfloat, spaces, m1a = unpack(fmt, datai)
            nspaces = len(spaces)
            spaces = spaces.decode('latin1').strip()
            print((aint, bint, cfloat, spaces, m1a), nspaces)
            #self.show_data(datai, types='qds', endian=None, force=False)
        #aaa
        #ints2 = [val if val != 2314885530818453536 else ' '
                 #for val in ints]
        #ints2[2] = cfloat
        #print(ints2)
        #ints = np.frombuffer(data, op2.idtype8).tolist()
        #floats = np.frombuffer(data, op2.fdtype8).tolist()

        #self.show_data(data, types='qds', endian=None, force=False)
        self.read_3_markers([-5, 1, 0, 0])

        #self.show(32, types='ifqds', endian=None, force=False)
        #asdf

    def read_gpl(self):
        """
        reads the GPL table (grid point list?)

        tested by TestOP2.test_beam_modes

        """
        if self.read_mode == 1:
            read_record = self._skip_record
        else:
            read_record = self._read_record

        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)

        self.read_markers([-1])
        header_data = self._read_record()  # (102, 117, 0, 0, 0, 0, 0)
        ints = np.frombuffer(header_data, op2.idtype)

        #seid = ints[0] # ???
        nnodes = ints[1]

        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        #self.show_data(unused_data)
        #print('--------------------')

        self.read_3_markers([-2, 1, 0])
        self.read_table_name(['GPL', 'GPLOUT'])
        #else ndata == 12:  # TestOP2Matrix.test_gpspc
        #print('--------------------')

        self.read_3_markers([-3, 1, 0])
        unused_data = read_record() # nids 1-117

        self.read_3_markers([-4, 1, 0])
        data = read_record()
        if self.read_mode == 2 and self.size == 4:
            # nids 1-117 (column 1) with nid*1000 (column 2)
            #
            # External grid or scalar identification number = node_id
            # Sequence number = 1000 * external identification number
            unused_nid_seq = np.frombuffer(data, op2.idtype).reshape(nnodes, 2)
        self.read_markers([-5, 1, 0, 0])

    def read_gpls(self):
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)

        self.read_markers([-1])
        data = self._read_record()  # (101, 139, 0, 0, 0, 0, 0)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        ndata = len(data)
        method = 0
        if self.size == 4:
            if ndata == 8:
                # 'GPL     ' in GPLS table
                # osmpf1_results.op2
                gpl_gpls, = Struct(self._endian + b'8s').unpack(data)
            else:
                gpl_gpls, method = Struct(self._endian + b'8si').unpack(data)
        else:
            if ndata == 16:
                gpl_gpls, = Struct(self._endian + b'16s').unpack(data)
            else:
                gpl_gpls, method = Struct(self._endian + b'16sq').unpack(data)
        assert gpl_gpls.strip() in [b'GPL', b'GPLS'], gpl_gpls.strip()
        assert method in [0, 1, 2, 3, 4, 5, 6, 7, 10, 12, 13, 15, 20, 30, 40, 99,
                          101, 201], f'GPLS method={method}'

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        ints = np.frombuffer(data, op2.idtype)
        #print(ints)

        self.read_3_markers([-4, 1, 0])
        data = self._read_record()
        ints = np.frombuffer(data, op2.idtype)
        nints = len(ints)
        ints = ints.reshape(nints//2, 2)
        #print(ints)

        self.read_markers([-5, 1, 0, 0])
        #data = self._read_record()
        #self.show_data(data, types='ifs', endian=None)

    def read_table_name(self, table_names: list[bytes]) -> str:
        if self.size == 4:
            return self.read_table_name4(table_names)
        return self.read_table_name8(table_names)

    def read_table_name4(self, table_names: list[bytes]) -> str:
        assert isinstance(table_names, list), table_names
        data, ndata = self._read_record_ndata4() # GPL
        if ndata == 8:
            table_name_bytes, = self.op2.struct_8s.unpack(data)
        elif ndata == 12:
            table_name_bytes, zero = self.op2.struct_8s_i.unpack(data)
            assert zero == 0, self.show_data(data)
        else:
            self.show_data(data)
            raise SubTableReadError('cannot read table_name=%r' % table_names)

        table_name_str = table_name_bytes.decode('utf-8').strip()
        assert table_name_str in table_names, f'actual={table_name_str} allowed={table_names}'
        return table_name_str

    def read_table_name8(self, table_names: list[bytes]) -> str:
        assert isinstance(table_names, list), table_names
        data, ndata = self._read_record_ndata8() # GPL
        if ndata == 16:
            table_name_bytes, = self.op2.struct_16s.unpack(data)
        elif ndata == 24:
            table_name_bytes, zero = self.op2.struct_16s_q.unpack(data)
            assert zero == 0, self.show_data(data)
        else:
            print(ndata)
            self.show_data(data, types='ifsq')
            #self.show_data(data[16:], types='ifsq')
            raise SubTableReadError(f'cannot read table_name={table_names}')

        table_name_bytes = reshape_bytes_block(table_name_bytes)
        table_name_str = table_name_bytes.decode('utf-8').strip()
        assert table_name_str in table_names, f'actual={table_name_str} allowed={table_names}'
        return table_name_str

    def read_gpdt(self):
        """
        reads the GPDT table

        tested by ???

        """
        #if self.read_mode == 1:
            #read_record = self._skip_record
        #else:
        read_record = self._read_record
        skip_record = self._skip_record

        op2 = self.op2
        table_name = self._read_table_name(rewind=False)
        op2.table_name = table_name
        #self.log.debug('table_name = %r' % table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_gpdt - %s\n' % table_name)

        self.read_markers([-1])
        header_data = self._read_record()  # (103, 117, 0, 0, 0, 0, 0)
        header_ints = np.frombuffer(header_data, op2.idtype8)
        header_floats = np.frombuffer(header_data, op2.fdtype8)

        # what a mess...
        # there are 3 ways to define the header...
        # for each of the 3, there are size 7/10 and single/double precision...
        #[ 102 3       0    0    0    0    0]  # 3 nodes....7 words   C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\extse11s.op2
        #[ 102 28      0    0    0    0    0]  # 16 nodes...28 bytes  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\boltld06.op2
        #[ 102 5523    7    0    1    0    0]  # 5523 nodes...7 words C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\acssnbena1.op2
        #print('gpdt ints =', header_ints)
        #print('gpdt floats =', header_floats)
        #seid = ints[0] # ??? is this a table number>
        #nnodes = ints[1]

        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        #print('--------------------')

        self.read_3_markers([-2, 1, 0])
        self.read_table_name(['GPDT', 'GPDTS', 'SAGPDT'])

        #print('--------------------')

        self.read_3_markers([-3, 1, 0])


        ## TODO: no idea how this works...
        read_mode = self.read_mode
        result_name = 'gpdt'
        #if op2._results.is_not_saved(result_name):
            #read_mode = 2 #  fake the read mode to skip the table
        #else:
            #op2._results._found_result(result_name)

        if read_mode == 1 and op2._results.is_saved(result_name):
            op2._results._found_result(result_name)
            i = 0
            data = read_record() # nid,cp,x,y,z,cd,ps
            ndata = len(data)



            # assume nodes
            #_get_gpdt_nnodes(self.size, ndata, header_ints, self.log)
            nnodes, numwide = _get_gpdt_nnodes2(ndata, header_ints, self.size)
            if (self.size, numwide) == (8, 14):
                # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\z402cdamp1_04.op2
                nid_cp_cd_ps, xyz = _read_gpdt_8_14(op2, data, nnodes)
            elif (self.size, numwide) == (4, 7):
                nid_cp_cd_ps, xyz = _read_gpdt_4_7(op2, data, nnodes)
            elif (self.size, numwide) == (4, 10):
                nid_cp_cd_ps, xyz = _read_gpdt_4_10(op2, data, nnodes)
            else:
                self.show_data(data, types='if')
                #ndata = 84:
                #[102   3   0   0   0   0   0]
                #ints    = (1,   0, 150.0,   0,   0,   0,   0,
                           #2, 0, 1125515264, 1103626240, 0, 0, 0,
                           #3, 0, 1125515264, -1043857408, 0, 0, 0)
                #floats  = (1, 0.0, 150.0, 0.0, 0.0, 0.0, 0.0,
                           #2, 0.0, 150.0, 25.0, 0.0, 0.0, 0.0,
                           #3, 0.0, 150.0, -25.0, 0.0, 0.0, 0.0)
                raise NotImplementedError((self.size, numwide))

            if 0:  # pragma: no cover
                if nvalues % 7 == 0:
                    # mixed ints, floats
                    #  0   1   2   3   4   5   6
                    # id, cp, x1, x2, x3, cd, ps
                    nrows = get_table_size_from_ncolumns('GPDT', nvalues, 7)
                    #print('self.size =', self.size)
                    ntotal = 28
                    structi = Struct(self._endian + b'2i 4f i')
                    for j in range(10):
                        edata = data[i:i+ntotal]
                        self.show_data(edata, types='ifqd')
                        out = structi.unpack(edata)
                        i += ntotal
                        #print(out)
                    asdf
                    ints = np.frombuffer(data, op2.idtype8).reshape(nnodes, 7).copy()
                    floats = np.frombuffer(data, op2.fdtype8).reshape(nnodes, 7).copy()
                    iints = [0, 1, 5, 6] # [1, 2, 6, 7] - 1
                    nid_cp_cd_ps = ints[:, iints]
                    xyz = floats[:, 2:5]
                elif nvalues % 10 == 0:
                    # mixed ints, doubles
                    nrows = get_table_size_from_ncolumns('GPDT', nvalues, 10)
                    iints = [0, 1, 8, 9]
                    #ifloats = [2, 3, 4, 5, 6, 7]
                    idoubles = [1, 2, 3]
                    unused_izero = [1, 8, 9]

                    #print('ints:')
                    #print(ints)
                    # nid cp, x,    y,    z,    cd, ps
                    # [0, 1,  2, 3, 4, 5, 6, 7, 8,   9]
                    # [ ,  ,  1, 1, 2, 2, 3, 3,  ,    ]
                    if self.read_mode == 1:
                        ints = np.frombuffer(data, op2.idtype).reshape(nnodes, 10).copy()
                        #floats = np.frombuffer(data, op2.fdtype).reshape(nrows, 10).copy()
                        doubles = np.frombuffer(data, 'float64').reshape(nnodes, 5).copy()

                        nid_cp_cd_ps = ints[:, iints]
                        xyz = doubles[:, idoubles]
                else:
                    raise NotImplementedError(nvalues)

            self.op2.op2_results.gpdt = GPDT(nid_cp_cd_ps, xyz)
        else:
            unused_data = skip_record() # nid,cp,x,y,z,cd,ps


        # 1. Scalar points are identified by CP=-1 and words X1 through
        #    PS are zero.
        # 3. or fluid grid points, CD=-1.
        #print(nid_cp_cd_ps)
        #print(xyz)


        isubtable = -4
        markers = self.get_nmarkers(1, rewind=True)

        if markers[0] != isubtable:
            self.read_markers([markers[0], 1, 0, 0])
            self.show(200)
            self.log.error('unexpected GPDT marker marker=%s; expected=%s' % (
                markers[0], isubtable))
            #markers = self.get_nmarkers(1, rewind=False)
            return

        self.read_3_markers([isubtable, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            #markers = self.get_nmarkers(1, rewind=True)
            #self.log.debug('GPDT record; markers=%s' % str(markers))
            if self.read_mode == 1:
                self._skip_record()
            else:
                #self.log.debug('unexpected GPDT record; markers=%s' % str(markers))
                data = self._read_record()
                #print('read_mode=%s freqs=%s' % (self.read_mode, freqs.tolist()))

            markers = self.get_nmarkers(1, rewind=True)
            self.read_3_markers([isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
            isubtable -= 1
        del isubtable
        self.read_markers([0])

    def read_bgpdt(self):
        """
        reads the BGPDT, BGPDTS, BGPDTOLD tables

        tested by TestOP2Matrix.test_gpspc

        """
        #if self.read_mode == 1:
            #read_record = self._skip_record
        #else:
        read_record = self._read_record

        op2 = self.op2
        table_name = self._read_table_name(rewind=False)
        op2.table_name = table_name
        #self.log.debug('table_name = %r' % table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_bgpdt - %s\n' % table_name)

        self.read_markers([-1])
        header_data = self._read_record()  # (105, 51, 0, 0, 0, 0, 0)
        header_ints = np.frombuffer(header_data, op2.idtype8)
        header_floats = np.frombuffer(header_data, op2.fdtype8)
        #print('bgpdt ints =', header_ints)
        #print('bgpdt floats =', header_floats)

        #seid = ints[0] # ??? is this a table number?
        nnodes = header_ints[1]  # validated
        #print('nnodes =', nnodes)
        #print(ints)
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        #print('--------------------')

        self.read_3_markers([-2, 1, 0])
        self.read_table_name(['BGPDT', 'BGPDTS', 'BGPDTOLD', 'BGPDTOUT'])

        #print('--------------------')

        self.read_3_markers([-3, 1, 0])

        #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\s402_sphere_03.op2
        #GRID 1 0 0.0    0.0 0.0      0
        #GRID 2 0 0.0    0.0 0.0      0
        #GRID 3 0 0.0871 0.0 -0.99619 0
        #D = (0.0, 5e-324,
        #     5e-324, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5e-323,
        #     1e-323, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.4e-323,
        #     1.5e-323, 3e-322, 0.0, 0.0, 0.0871, 0.0, -0.99619)
        #L = (0, 1,
        #     1, 61, 0, 0, 0, 0, 0, 0, 7,
        #     2, 61, 0, 0, 0, 0, 0, 0, 13,
        #     3, 61, 0, 0, 0.0871, 0, -0.99619)

        #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\s402_flxslddriver_05.op2
        #doubles (float64) = (0.0, 5e-324, 5e-324, 3e-322, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5e-323, 1e-323, 3e-322, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 6.4e-323, 1.5e-323, 3e-322, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 9.4e-323, 2e-323, 3e-322, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 1.24e-322, 2.5e-323, 3e-322, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 1.53e-322, 3e-323, 3e-322, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 1.83e-322, 3.5e-323, 3e-322, 0.0, 0.0, 60.0, 0.0, 0.0, 0.0, 2.1e-322, 4e-323, 3e-322, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 2.4e-322, 4.4e-323, 3e-322, 0.0, 0.0, 80.0, 0.0, 0.0, 0.0, 2.7e-322, 5e-323, 3e-322, 0.0, 0.0, 90.0, 0.0, 0.0, 0.0, 3e-322, 5.4e-323, 3e-322, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 3.3e-322, 6e-323, 3e-322, 0.0, 0.0, 5.0, 0.0, 20.0, 0.0, 3.6e-322, 6.4e-323, 3e-322, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 3.9e-322, 7e-323, 3e-322, 0.0, 0.0, 10.0, 0.0, 0.0)
        #long long (int64) = (
            #0, 1,
            #1, 61, 0, 0, 0, 0, 0, 0, 7,
            #2, 61, 0, 0, 4621819117588971520, 0, 0, 0, 13,
            #3, 61, 0, 0, 4626322717216342016, 0, 0, 0, 19,
            #4, 61, 0, 0, 4629137466983448576, 0, 0, 0, 25,
            #5, 61, 0, 0, 4630826316843712512, 0, 0, 0, 31,
            #6, 61, 0, 0, 4632233691727265792, 0, 0, 0, 37,
            #7, 61, 0, 0, 4633641066610819072, 0, 0, 0, 43,
            #8, 61, 0, 0, 4634626229029306368, 0, 0, 0, 49,
            #9, 61, 0, 0, 4635329916471083008, 0, 0, 0, 55,
            #10, 61, 0, 0, 4636033603912859648, 0, 0, 0, 61,
            #11, 61, 0, 0, 4636737291354636288, 0, 0, 0, 67,
            #12, 61, 0, 0, 4617315517961601024, 0, 4626322717216342016, 0, 73,
            #13, 61, 0, 0, 4617315517961601024, 0, 0, 0, 79,
            #14, 61, 0, 0, 4621819117588971520, 0, 0)

        if self.read_mode == 1:
            self._skip_record()

        elif self.read_mode == 2:
            i = 0
            data = read_record() # cd,x,y,z
            # xword = 4 * self.factor
            nvalues4 = len(data) // 4
            assert len(data) % 4 == 0, len(data) % 4
            # assert len(data) % xword == 0

            numwide = nvalues4 // nnodes
            result_name = 'bgpdt'
            if op2._results.is_saved(result_name):
                op2._results._found_result(result_name)

                if self.size == 4:
                    #C:\NASA\m4\formats\git\examples\x33_blog\blog_materials\pressure_vessel_fem1_sim1-my_limit_load.bdf
                    #GRID*                  1               0-1.792343211E+003.4539249252E+00+
                    #*       3.6302008282E-01               0
                    #GRID*                  2               0-1.792313297E+003.3970433059E+00+
                    #*       7.2220293673E-01               0
                    #GRID*                  3               0-1.792337984E+003.3029109367E+00+
                    #*       1.0733895770E+00               0
                    #GRID*                  4               0-1.792341147E+003.1725968908E+00+
                    #*       1.4128346960E+00               0
                    #GRID*                  5               0-1.792333634E+003.0076863630E+00+
                    #*       1.7364594482E+00               0

                    assert nvalues4 % nnodes == 0, nvalues4 % nnodes
                    #print(self.size, numwide, nvalues4, len(data))

                    if numwide == 4:
                        assert numwide == 4, numwide
                        ntotal = 16 # 4*4
                        assert len(data) % ntotal == 0, len(data) % ntotal
                        structi = Struct(self._endian + b'i3f')  # 16
                        for j in range(nnodes):
                            out = structi.unpack(data[i:i+ntotal])
                            cd, x, y, z = out
                            outs = f'nid={j+1} cd={cd} xyz=({x:g},{y:g},{z:g})'
                            #print(outs)
                            i += ntotal
                        # [cd, x, y, z]
                        cd = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 4)[:, :-3]
                        xyz = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 4)[:, -3:]

                    elif numwide == 12:
                        ntotal = 48 # (9+3)*4
                        structi = Struct(self._endian + b'6i3d')
                        # len(data) // 4 / ngrids = 4
                        # 30768 // 4 / 1923
                        nrows = get_table_size_from_ncolumns('BGPDT', nvalues4, 12)
                        for j in range(nnodes):
                            edata = data[i:i+ntotal]
                            #self.show_data(edata, 'ifqd')
                            out = structi.unpack(edata)

                            #self.show_data(data[i:i+ntotal], types='qd')
                            #print(out)
                            cd, sil, nid, sixtyone, ps, zero_c, x, y, z = out
                            outs = f'cd={cd} sil={sil} nid={nid} sixtyone={sixtyone} ps={ps} zero_c={zero_c} xyz=({x:g},{y:g},{z:g})'
                            assert nid > 0, outs
                            #print(outs)
                            #assert zero_a in [0, 225, 362, 499], out
                            #assert zero_c == 0, outs
                            #11: C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\gluedg01f.op2
                            assert sixtyone in [11, 12, 61], f'sixtyone={sixtyone} outs={outs}'
                            i += ntotal

                        # [cd, x, y, z]
                        #self.show_data(data, types='if', endian=None, force=False)
                        #[0, num, nid, 61, 0, 0, x, _, y, _, z, _]
                        ints = np.frombuffer(data, op2.idtype).reshape(nrows, 12).copy()
                        floats = np.frombuffer(data, 'float64').reshape(nrows, 6).copy()
                        cd = ints[:, 0]
                        xyz = floats[:, 1:]
                        #print(ints[:-6])
                        #print(xyz)
                    else:  # pragma: no cover
                        raise NotImplementedError((self.size, numwide))
                    op2.op2_results.bgpdt = BGPDT(cd, xyz)
                else:
                    #bad = []
                    #nvalues = len(data) // 4
                    #for i in [2, 3, 6]: # 2-16 checked
                        #if nvalues % i != 0:
                            #bad.append(i)
                    #if bad:
                        #print(nvalues, bad)
                    #self.show_data(data, types='ifqd')
                    #print(nvalues4)
                    # 112 / 28.
                    #nrows9 = nvalues4 // 18
                    #nrows4 = nvalues4 // 8
                    #nvalues_9 = nvalues4 % 18
                    #nvalues_4 = nvalues4 % 8
                    #print('nvalues_9', nvalues_9, nvalues_9)
                    #print('nvalues_4', nvalues_4, nvalues_4)
                    if numwide == 18: # nvalues_9 == 0 and nvalues_4 != 0:
                        #assert nvalues4 % 18 == 0, nvalues4 % 18
                        #print(nrows)
                        #i = np.arange(0, nrows)
                        #i1 = i * 4
                        #i2 = (i + 1) * 4
                        #zero_a, sil, nid, sixtyone, ps, zero_c, x, y, z
                        ints = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 9)[:, :-3]
                        floats = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 9)[:, -3:]
                        nid = ints[:, 2]
                        cd = ints[:, [4, 5]]
                        xyz = floats

                        #GRID           1       0    27.5    20.0     0.0       0
                        #GRID           2       0 25.30329 25.30329     0.0       0
                        #GRID           3       0    20.0    27.5     0.0       0
                        #GRID           4       0 14.69671 25.30329     0.0       0

                        #GRID           1       0    27.5    20.0     0.0       0
                        #GRID           2       025.3032925.30329     0.0       0
                        #GRID           3       0    20.0    27.5     0.0       0
                        #GRID           4       014.6967125.30329     0.0       0
                        #GRID           5       0    12.5    20.0     0.0       0
                        #GRID           6       014.6967114.69671     0.0       0
                        #GRID           7       0    20.0    12.5     0.0       0
                        #GRID           8       025.3032914.69671     0.0       0
                        #GRID           9       0 22.3122 7.30685     0.0       0
                        #GRID          10       017.8372532.93714     0.0       0
                        #GRID          11       07.22962217.75039     0.0       0
                        #GRID          12       032.86626 17.7885     0.0       0
                        #GRID          13       0 30.56189.511215     0.0       0
                        #GRID          14       024.9442132.70425     0.0       0
                        #GRID          15       032.4224632.42249     0.0       0
                        #GRID          16       032.6858324.92913     0.0       0
                        #GRID          17       09.68575330.94779     0.0       0
                        #GRID          18       07.90302224.50037     0.0       0
                        #GRID          19       09.0897749.672205     0.0       0
                        #GRID          20       015.524677.918942     0.0       0
                        #GRID          21       0     0.0    32.0     0.0       0
                        #GRID          22       0     0.0    24.0     0.0       0
                        #GRID          23       0     0.0    16.0     0.0       0
                        #GRID          24       0     0.0     8.0     0.0       0
                        #GRID          25       0     8.0     0.0     0.0       0
                        #GRID          26       0    16.0     0.0     0.0       0
                        #GRID          27       0    24.0     0.0     0.0       0
                        #GRID          28       0    32.0     0.0     0.0       0
                        #GRID          29       0    32.0    40.0     0.0       0
                        #GRID          30       0    24.0    40.0     0.0       0
                        #GRID          31       0    16.0    40.0     0.0       0
                        #GRID          32       0     8.0    40.0     0.0       0
                        #GRID          33       0    40.0     8.0     0.0       0
                        #GRID          34       0    40.0    16.0     0.0       0
                        #GRID          35       0    40.0    24.0     0.0       0
                        #GRID          36       0    40.0    32.0     0.0       0
                        #GRID          37       0     0.0     0.0     0.0       0
                        #GRID          38       0     0.0    40.0     0.0       0
                        #GRID          39       0    40.0     0.0     0.0       0
                        #GRID          40       0    40.0    40.0     0.0       0
                        #GRID          41       0    12.5    20.0    20.1       0
                        #(0, 1, 1, 61, 0, 0, 27.5, 20.0, 0.0)
                        #(0, 7, 2, 61, 0, 0, 25.30329, 25.30329, 0.0)
                        #(0, 13, 3, 61, 0, 0, 20.0, 27.5, 0.0)
                        #(0, 19, 4, 61, 0, 0, 14.69671, 25.30329, 0.0)
                        #(0, 25, 5, 61, 0, 0, 12.5, 20.0, 0.0)
                        #(0, 31, 6, 61, 0, 0, 14.69671, 14.69671, 0.0)
                        #(0, 37, 7, 61, 0, 0, 20.0, 12.5, 0.0)
                        #(0, 43, 8, 61, 0, 0, 25.30329, 14.69671, 0.0)
                        #(0, 49, 9, 61, 0, 0, 22.3122, 7.30685, 0.0)
                        #(0, 55, 10, 61, 0, 0, 17.83725, 32.93714, 0.0)
                        #(0, 61, 11, 61, 0, 0, 7.229622, 17.75039, 0.0)
                        #(0, 67, 12, 61, 0, 0, 32.86626, 17.7885, 0.0)
                        #(0, 73, 13, 61, 0, 0, 30.5618, 9.511215, 0.0)
                        #(0, 79, 14, 61, 0, 0, 24.94421, 32.70425, 0.0)
                        #(0, 85, 15, 61, 0, 0, 32.42246, 32.42249, 0.0)
                        #(0, 91, 16, 61, 0, 0, 32.68583, 24.92913, 0.0)
                        #(0, 97, 17, 61, 0, 0, 9.685753, 30.94779, 0.0)
                        #(0, 103, 18, 61, 0, 0, 7.903022, 24.50037, 0.0)
                        #(0, 109, 19, 61, 0, 0, 9.089774, 9.672205, 0.0)
                        #(0, 115, 20, 61, 0, 0, 15.52467, 7.918942, 0.0)
                        #(0, 121, 21, 61, 0, 0, 0.0, 32.0, 0.0)
                        #(0, 127, 22, 61, 0, 0, 0.0, 24.0, 0.0)
                        #(0, 133, 23, 61, 0, 0, 0.0, 16.0, 0.0)
                        #(0, 139, 24, 61, 0, 0, 0.0, 8.0, 0.0)
                        #(0, 145, 25, 61, 0, 0, 8.0, 0.0, 0.0)
                        #(0, 151, 26, 61, 0, 0, 16.0, 0.0, 0.0)
                        #(0, 157, 27, 61, 0, 0, 24.0, 0.0, 0.0)
                        #(0, 163, 28, 61, 0, 0, 32.0, 0.0, 0.0)
                        #(0, 169, 29, 61, 0, 0, 32.0, 40.0, 0.0)
                        #(0, 175, 30, 61, 0, 0, 24.0, 40.0, 0.0)
                        #(0, 181, 31, 61, 0, 0, 16.0, 40.0, 0.0)
                        #(0, 187, 32, 61, 0, 0, 8.0, 40.0, 0.0)
                        #(0, 193, 33, 61, 0, 0, 40.0, 8.0, 0.0)
                        #(0, 199, 34, 61, 0, 0, 40.0, 16.0, 0.0)
                        #(0, 205, 35, 61, 0, 0, 40.0, 24.0, 0.0)
                        #(0, 211, 36, 61, 0, 0, 40.0, 32.0, 0.0)
                        #(0, 217, 37, 61, 0, 0, 0.0, 0.0, 0.0)
                        #(0, 223, 38, 61, 0, 0, 0.0, 40.0, 0.0)
                        #(0, 229, 39, 61, 0, 0, 40.0, 0.0, 0.0)
                        #(0, 235, 40, 61, 0, 0, 40.0, 40.0, 0.0)
                        #(0, 241, 41, 61, 0, 0, 12.5, 20.0, 20.1)
                        #doubles (float64) = (0, 7, 2, 61, 0.0, 0.0, 25.30329, 25.30329, 0.0)
                        #long long (int64) = (0, 7, 2, 61, 0, 0, 25.30329, 25.30329, 0)
                        #self.show_data(data[:80], types='qd')
                        ntotal = 72 # 9*8
                        structi = Struct(self._endian + b'6q3d')
                        for j in range(nnodes):
                            out = structi.unpack(data[i:i+ntotal])
                            #self.show_data(data[i:i+ntotal], types='qd')
                            #print(out)
                            cd, sil, nid, sixtyone, ps, zero_c, x, y, z = out
                            outs = f'cd={cd} sil={sil} nid={nid} sixtyone={sixtyone} ps={ps} zero_c={zero_c} xyz=({x:g},{y:g},{z:g})'
                            assert nid > 0, outs

                            assert zero_c == 0, outs
                            assert sixtyone == 61, outs
                            i += ntotal


                    # elif nvalues_9 != 0 and nvalues_4 == 0:
                    elif numwide == 8:
                        #assert numwide == 8, numwide
                        ntotal = 32 # 24+4 = 28
                        assert len(data) % ntotal == 0, len(data) % ntotal
                        structi = Struct(self._endian + b'q3d')  # 28
                        for j in range(nnodes):
                            out = structi.unpack(data[i:i+ntotal])
                            #self.show_data(data[i:i+ntotal], types='ifqd')
                            #print(out)
                            cd, x, y, z = out
                            i += ntotal
                        # [cd, x, y, z]
                        cd = np.frombuffer(data, op2.idtype8).copy().reshape(nnodes, 4)[:, :-3]
                        xyz = np.frombuffer(data, op2.fdtype8).copy().reshape(nnodes, 4)[:, -3:]
                    else:
                        raise RuntimeError((self.size, numwide))
                    #cd = ints[::7]
                    # [cd, x, _, y, _, z, _]

                    #print(ints)
                    #print(floats)
                    #print(nrows*7, len(floats))
                    #ints = ints[2:].reshape(nrows, 7)
                    #floats = floats[2:].reshape(nrows, 7)
                    #for inti, floati in zip(ints, floats):
                        #print(inti[:-3], floats[-3:])
                    #print(xyz, xyz.shape)
                    op2.op2_results.bgpdt = BGPDT(cd, xyz)
                #print('cd = %s' % cd.tolist())
                #print('xyz:\n%s' % xyz)

        self.read_3_markers([-4, 1, 0])
        marker = self.get_nmarkers(1, rewind=True)[0]
        if marker == 0:
            self.read_markers([0])
            return

        ## TODO: why is this needed??? (it is, but dmap is not clear)
        data = self._read_record()
        #self.show_data(data, types='i')
        isubtable = -5
        while 1:
            self.read_3_markers([isubtable, 1, 0])
            marker = self.get_nmarkers(1, rewind=True)[0]
            if marker == 0:
                break
            data = self._read_record()
            isubtable -= 1

        self.read_markers([0])

    def read_hisadd(self):
        """optimization history (SOL200) table"""
        op2 = self.op2

        result_name = 'responses.convergence_data'
        is_saved_result = op2._results.is_saved(result_name)
        #read_mode = self.read_mode
        if is_saved_result:
            op2._results._found_result(result_name)
        else:
            self._skip_table('HISADD', warn=False)
            return
        is_not_saved_result = not is_saved_result

        #slot = op2.get_result(result_name)

        responses = op2.op2_results.responses
        op2.table_name = self._read_table_name(rewind=False)

        if self.read_mode == 1 or is_not_saved_result:
            self.read_markers([-1])
            self._skip_record()
            self.read_3_markers([-2, 1, 0])
            self._skip_record()
            self.read_3_markers([-3, 1, 0])

            if is_saved_result and responses.convergence_data is None:
                data = self._read_record()
                ndvs = len(data) // 4 - 7
                responses.convergence_data = Convergence(ndvs)
            #elif is_not_saved_result:
                #self._skip_record()
            else:
                self._skip_record()
                responses.convergence_data.n += 1

            self.read_markers([-4, 1, 0, 0])
            return

        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        #self.log.info('----marker1----')
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()  # ()102, 303, 0, 0, 0, 0, 0) date???
        #print('hisadd data1')
        #self.show_data(data)

        #self.log.info('----marker2----')
        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()  # ('HISADD', )
        #print('hisadd data2')
        #self.show_data(data)

        #self.log.info('----marker3----')
        self.read_3_markers([-3, 1, 0])
        data = self._read_record()

        fmt = mapfmt(self._endian + b'3i3fi', self.size)
        (design_iter, iconvergence, conv_result, obj_intial, obj_final,
         constraint_max, row_constraint_max) = unpack(fmt, data[:28 * self.factor])

        if iconvergence == 1:
            iconvergence = 'soft'
        elif iconvergence == 2:
            iconvergence = 'hard'
        elif iconvergence == 6:
            self.log.warning('HISADD iconverge=6')
            iconvergence = '???'
        else:  # pragma: no cover
            msg = 'iconvergence=%s\n' % iconvergence
            self.show_data(data, types='ifs', endian=None)
            raise NotImplementedError(msg)

        if conv_result == 0:
            conv_result = 'no'
        elif conv_result == 1:
            conv_result = 'soft'
        elif conv_result == 2:
            conv_result = 'hard'
        elif conv_result in [3, 4]:
            #self.log.warning('HISADD conv_result=%s' % conv_result)
            # not sure why this happens, but the field is wrong
            # it seems to apply to one step before this one
            conv_result = 'best_design'
        else:
            self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                           'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                               design_iter, iconvergence, conv_result, obj_intial,
                               obj_final, constraint_max, row_constraint_max))
            raise NotImplementedError('conv_result=%s' % conv_result)
        #self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                       #'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                           #design_iter, iconvergence, conv_result, obj_intial,
                           #obj_final, constraint_max, row_constraint_max))

        ndvs = len(data) // 4 - 7
        desvar_values = unpack('%sf' % ndvs, data[28:])

        responses.convergence_data.append(
            design_iter, iconvergence, conv_result, obj_intial,
            obj_final, constraint_max, row_constraint_max, desvar_values)
        self.read_markers([-4, 1, 0, 0])

    def read_ibulk(self):
        """
        tested by TestOP2.test_ibulk

        read_mode = 1 (array sizing)
        read_mode = 1 (reading)

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()
        #print(self.show_data(data))

        unused_markers = self.get_nmarkers(1, rewind=True)

        if self._dump_deck:
            write_deck = True
            bulk_filename = 'bulk.test_op2.bdf'
            save_lines = False
            self._read_deck_section(bulk_filename, save_lines,
                                    write_deck, mode='w',
                                    read_mode=1)
        else:
            save_lines = False
            write_deck = False
            bulk_filename = 'bulk.test_op2.bdf'
            self._read_deck_section(bulk_filename, save_lines,
                                    write_deck, mode='a',
                                    read_mode=2)


    def _read_deck_section(self, deck_filename: str, save_lines: bool,
                           write_deck: bool, mode='w', read_mode: int=1) -> list[str]:
        """helper for ``read_ibulk`` and ``read_icase``"""
        marker = -2
        if save_lines and write_deck:
            raise RuntimeError(f'save_lines={save_lines} write_deck={write_deck}; '
                               'one or more must be False')

        lines = []
        size = self.size
        if write_deck and self.read_mode == read_mode:
            with open(deck_filename, mode) as bdf_file:  # pragma: no cover
                while 1:
                    self.read_3_markers([marker, 1, 0])
                    nfields = self.get_marker1(rewind=True)
                    if nfields > 0:
                        # we're not reading the record because the IBULK/ICASE
                        # table is literally just an unsorted echo of the
                        # BULK/CASE data table in the BDF
                        data = self._read_record()
                        line = reshape_bytes_block_size(data.replace(b'\xff', b' '), size=size)
                        bdf_file.write(line + '\n')
                    elif nfields == 0:
                        #op2.show_ndata(100, types='ifs')
                        break
                    else:
                        raise RuntimeError('nfields=%s' % nfields)
                    marker -= 1

        elif save_lines and self.read_mode == read_mode:
            while 1:
                self.read_3_markers([marker, 1, 0])
                nfields = self.get_marker1(rewind=True)
                if nfields > 0:
                    # we're not reading the record because the IBULK
                    # table is literally just an unsorted echo of the
                    # BULK data table in the BDF
                    data = self._read_record()
                    line = reshape_bytes_block_size(data.replace(b'\xff', b' '), size=size)
                    lines.append(line + '\n')
                elif nfields == 0:
                    #op2.show_ndata(100, types='ifs')
                    break
                else:
                    raise RuntimeError('nfields=%s' % nfields)
                marker -= 1
        else:
            while 1:
                self.read_3_markers([marker, 1, 0])
                nfields = self.get_marker1(rewind=True)
                if nfields > 0:
                    # we're not reading the record because the IBULK
                    # table is literally just an unsorted echo of the
                    # BULK data table in the BDF
                    unused_data = self._skip_record()
                elif nfields == 0:
                    #op2.show_ndata(100, types='ifs')
                    break
                else:
                    raise RuntimeError('nfields=%s' % nfields)
                marker -= 1
        #print("marker = ", marker)
        unused_marker_end = self.get_marker1(rewind=False)
        return lines

    def read_icase(self):
        """
        tested by ???

        read_mode = 1 (array sizing)
        read_mode = 1 (reading)

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

        unused_markers = self.get_nmarkers(1, rewind=True)

        save_lines = True
        write_deck = False
        bulk_filename = 'bulk.test_op2.bdf'
        lines = self._read_deck_section(bulk_filename, save_lines,
                                        write_deck=write_deck, mode='w',
                                        read_mode=1)
        if self._dump_deck:
            with open(bulk_filename, 'a') as bdf_file:
                bdf_file.writelines(lines)
        self._case_control_lines = lines

    def read_casecc(self):
        """reads the CASECC table"""
        op2 = self.op2
        size = self.size
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)

        #print('reading -1')
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')

        # (101, 1, 0, 1237, 0, 0, 0)
        data = self._read_record()
        #self.show_data(data, types='iq', endian=None, force=False)

        #print('reading -2, 1, 0')
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        #self.show_data(data, types='dqs', endian=None, force=False)

        #self.read_3_markers([-3, 1, 0])
        #data = self._read_record()

        itable = -3
        #print(f'reading {itable}, 1, 0')
        self.read_3_markers([itable, 1, 0])
        marker = self.get_marker1(rewind=True, macro_rewind=False)

        from pyNastran.op2.tables.geom.subcase import set_casecc
        while marker != 0:
            itable -= 1
            data = self._read_record()
            #print(itable, len(data))
            try:
                subcase = set_casecc(self, data, op2.idtype8, op2.fdtype8, size=size,
                                     nastran_format=self.op2._nastran_format)
                self.op2.case_control_deck.subcases[subcase.id] = subcase
                #print(subcase)
            except Exception:
                pass #raise
            self.read_3_markers([itable, 1, 0])
            marker = self.get_marker1(rewind=True, macro_rewind=False)

        marker = self.get_marker1(rewind=False, macro_rewind=False)

    def read_xcasecc(self):
        r"""
        Poorly reads the XCASECC table.
        It's somehow aero related...

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\z402cbush1d_02a.op2
        word1 = b'XCASECC '
        word2 = (
            b'                                                                                                                                                                                                                                                                 '
            b'SUBCASE - NONLINEAR IMPLICIT                                                                           '
            b'SUBCASE 1               '
        word3 = b'AEROSG2D'
        word4 = b'YES '
        """
        op2 = self.op2
        size = self.size
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')

        if self.read_mode == 1:
            #(103, 1, 0, 1200, 0, 0, 0)
            data = self._skip_record()

            self.read_3_markers([-2, 1, 0])
            data = self._skip_record()

            self.read_3_markers([-3, 1, 0])
            data = self._skip_record()
        else:
            from pyNastran.op2.tables.geom.subcase import set_casecc
            #(103, 1, 0, 1200, 0, 0, 0)
            data = self._read_record()

            self.read_3_markers([-2, 1, 0])
            data = self._read_record()
            #word1 = reshape_bytes_block(data)
            #print('word1 =', word1)

            self.read_3_markers([-3, 1, 0])
            data = self._read_record()
            subcase = set_casecc(self, data, op2.idtype8, op2.fdtype8, size=size,
                                 nastran_format=self.op2._nastran_format)
            self.op2.case_control_deck.subcases[subcase.id] = subcase
            #print(subcase)
            if size == 8:
                word2 = reshape_bytes_block(data[38*size:134*size])
                #print('word2 =', word2)

                word3 = reshape_bytes_block(data[283*size:285*size])
                #print('word3 =', word3)

                word4 = reshape_bytes_block(data[518*size:519*size])
                #print('word4 =', word4)

        self.read_3_markers([-4, 1, 0, 0])

    def read_cddata(self):
        """Cambell diagram summary"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._read_record()
        #(101, 15, 5, 7, 0, 2, 0)
        #(101, 118, 4, 7, 0, 2, 0)

        # CDDATA, 2
        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        cddata, method = Struct(self._endian + b'8si').unpack(data)

        marker = -3
        cddata_list = []
        if method == 1:
            while 1:
                #print(f'read marker={marker}...')
                self.read_3_markers([marker, 1, 0])
                nfields = self.get_marker1(rewind=True)
                if nfields == 0:
                    break
                #print(nfields)
                data, ndata = self._read_record_ndata4()
                if self.read_mode == 1:
                    marker -= 1
                    continue
                #self.show_data(data, types='if', endian=None)
                marker -= 1
            #self.show(200)

        elif method == 2:
            while 1:
                self.read_3_markers([marker, 1, 0])
                nfields = self.get_marker1(rewind=True)
                if nfields == 0:
                    break
                data = self._read_record()
                if self.read_mode == 1:
                    marker -= 1
                    continue

                ints = np.frombuffer(data, op2.idtype)
                nints = len(ints)
                floats = np.frombuffer(data, op2.fdtype) # .reshape(nints//7, 7)

                # 1 NVAL I Number of values
                # 2 NCRV I Number of sections/solution (i.e., number of curves)
                # 3 KEYW I Keyword=10000+(SOLN*10)+DATTYP,
                # where:
                #   SOLN=solution number
                #   DATTYP=1 for list of rotor speeds in user-defined units
                #   DATTYP=2 for list of eigenfrequencies in the analysis system
                #   DATTYP=3 for list of Lehr damping values
                #   DATTYP=4 for list of real part of eigenvalues
                #   DATTYP=5 for list of imaginary part of eigenvalues
                #   DATTYP=6 for list of whirl direction codes (2.0=backwards, 3.0=forward, 4.0=linear)
                #   DATTYP=7 for list of converted frequencies in analysis system
                #   DATTYP=8 for list of whirl directions codes for converted solution (2.0=backwards, 3.0=forward, 4.0=linear)
                # 4 VALS(NVAL) RS list of values
                # Words 1–4 repeat for NCRV curves. For DATTYP≠1, NVAL and NCRV=0

                dict_map = {
                    1: 'RPM',
                    2: 'eigenfreq',
                    3: 'Lehr',
                    4: 'real eig',
                    5: 'imag eig',
                    6: 'whirl_dir',
                    7: 'converted_freq',
                    8: 'whirl_code',
                }
                i = 0
                data_out = {}
                while i < nints:
                    nvalues = ints[i]
                    #ncurves = ints[i+1]
                    keyword = ints[i+2]  # 10000+(SOLN*10)+DATTYP
                    base = keyword - 10000 # (SOLN*10)+DATTYP
                    #solution = base // 10
                    datatype = base % 10
                    assert datatype in [1, 2, 3, 4, 5, 6, 7, 8], datatype
                    values = floats[i+3:i+3+nvalues]
                    if datatype in [6, 8]:
                        values = values.astype('int32')
                    #print(f'{nvalues}, {ncurves}, {solution}, {datatype}, {dict_map[datatype]:14s} {values}')
                    data_out[datatype] = values
                    i += 3 + nvalues
                marker -= 1
                cddata_list.append(data_out)
                #import matplotlib.pyplot as plt

                #dict_map = {
                    #1: 'RPM',
                    #2: 'eigenfreq',
                    #3: 'Lehr',
                    #4: 'real eig',
                    #5: 'imag eig',
                    #6: 'whirl_dir',
                    #7: 'converted_freq',
                    #8: 'whirl_code',
                #}

                #plt.figure(2)
                #plt.plot(data_out[1], data_out[2]) # RPM vs. eigenfreq
                #plt.grid(True)

                #plt.figure(3)
                #plt.plot(data_out[1], data_out[3]) # RPM vs. Lehr
                #plt.grid(True)

                #plt.figure(4)
                #plt.plot(data_out[1], data_out[4]) # RPM vs. real_eig
                #plt.grid(True)

                #plt.figure(5)
                #plt.plot(data_out[1], data_out[5]) # RPM vs. imag_eig
                #plt.grid(True)

                #plt.figure(7)
                #plt.plot(data_out[1], data_out[7]) # RPM vs. converted_freq
                #plt.grid(True)

            #print('------------------------------------')
        else:
            raise NotImplementedError(f'CDDATA method={method}')
        if self.read_mode == 2:
            #plt.grid(True)
            #plt.show()
            self.op2.op2_results.cddata = cddata_list
        nfields = self.get_marker1(rewind=False)

    def read_stdisp(self):
        """reads the STDISP table"""
        #C:\NASA\m4\formats\git\examples\backup\aeroelasticity\loadf.op2
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        #(101, 1, 27, 0, 3, 1, 0)
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        assert Struct('8s').unpack(data)[0] == b'STDISP  '

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()
        a, b1, b2, b3, b4, b5, b6, b7, b8, b9 = Struct('64s 5i 12s 3i').unpack(data)
        #a 345 0 0 0 0 b'DISPSTWING  ' 1 0 123
        out = (b1, b2, b3, b4, b5, b6, b7, b8, b9)
        assert b1 == 345, out
        assert b2 == 0, out
        assert b3 == 0, out
        assert b4 == 0, out
        assert b5 == 0, out
        assert b6 == b'DISPSTWING  ', out
        assert b7 == 1, out
        assert b8 == 0, out
        assert b9 == 123, out
        #print(a, b1, b2, b3, b4, b5, b6, b7, b8, b9)
        #'STWING  COMPLETE STRUCTURAL WING                                '
        #(345, 0, 0, 0, 0, 1347635524, 1230459987, 538986318, 1, 0, 123)
        #(4.834479701920619e-43, 0.0, 0.0, 0.0, 0.0, 14179176448.0, 881989.1875, 1.35761196e-19, 1.4012e-45, 0.0, 1.72359)
        #self.show_data(data[96:], types='ifs', endian=None)

        self.read_3_markers([-4, 1, 0])
        self.read_markers([0])

    def read_omm2(self):
        """reads the OMM2 table"""
        op2 = self.op2
        #op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        else:
            raise NotImplementedError(self.show_data(data))
        self._read_subtables()

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???

        """
        self._skip_pcompts()
        return
        #if self.read_mode == 1:
            #return
        #op2.log.debug("table_name = %r" % op2.table_name)
        #table_name = self._read_table_name(rewind=False)

        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_3_markers([-2, 1, 0])
        #data = self._read_record()
        #table_name, = op2.struct_8s.unpack(data)
        ##print "table_name = %r" % table_name

        #self.read_3_markers([-3, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [-4]:
            #data = self._read_record()

        #self.read_3_markers([-4, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [0]:
            #data = self._read_record()
        #else:
            #self.read_markers([0])
            #return

        #self.read_3_markers([-5, 1, 0])
        #data = self._read_record()

        #self.read_3_markers([-6, 1, 0])
        #self.read_markers([0])

    def _skip_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???

        """
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        table_name = self._read_table_name(rewind=False)

        #read_record = self._skip_record if self.read_mode == 2 else self._read_record
        self.read_markers([-1])

        # (104, 0,   1,   3, 0, 0,   0) - tpl\qrcomp.op2 PCOMPTS (length 3?)
        # (104, 8,   0,   0, 0, 0,   0) - tpl\lmtas1.op2 PCOMPTS (return after -4)
        # (104, 0, 103, 412, 0, 0, 103) - output.op2     PCOMPT  (return at end)
        data_header = self._read_record()
        #self.show_data(data_header, types='ifsd', endian=None)
        a, bi, n4a, n5, e, f, n4b = Struct(b'<7i').unpack(data_header)
        self.log.debug('a=%s b=%s n4a=%s n5=%s e=%s f=%s n4b=%s' % (
            a, bi, n4a, n5, e, f, n4b))
        #assert n4a == n4b, 'n4a=%s n4b=%s' % (n4a, n4b)

        if table_name == b'PCOMPT':
            n4words = n4a * 4
            n5words = n5 * 32
        elif table_name == b'PCOMPTS':
            n4words = n4a * 3
            n5words = n5
        else:  # pragma: no cover
            raise NotImplementedError(table_name)

        self.read_3_markers([-2, 1, 0])
        # 'IPCOMPT '
        unused_data = self._read_record()
        #table_name, = op2.struct_8s.unpack(data)

        isubtable = -3
        self.read_3_markers([isubtable, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)

        if markers != [-4]:
            unused_data = self._read_record()

        self.read_3_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]: # n4a=0
            #self.show_data(data_header, types='ifsd')
            # (1, 4, 0, '    ')
            #(421, 4, 128, '    ')
            #(814, 4, 256, '    ')
            data = self._read_record()

            assert len(data) == n4words*4, 'n4words=%s len(data)=%s n4words*4=%s'  % (n4words, len(data), n4words*4)

            if table_name == b'PCOMPTS':
                # (11, 3, 0)
                pass
                #self.show_data(data, types='ifs', endian=None)
            elif 0: # pramga: no cover
                i = 0
                j = 0
                s1 = Struct(b'< 3i 4s')
                while i < len(data):
                    datai = data[i:i+16]
                    out = s1.unpack(datai)
                    print(out)
                    i += 16
                    j += 1
                print("j (-4) = %s" % j)
                assert i == len(data), '-4'
        else:
            self.read_markers([0])
            assert n4a == 0, n4a
            return

        self.read_3_markers([-5, 1, 0])
        data = self._read_record()
        assert len(data) == n5words * 4, 'n5words=%s len(data)=%s n5words*4=%s'  % (n5words, len(data), n5words*4)

        if table_name == b'PCOMPTS':
            # (101, 102, 103)
            pass
            #self.show_data(data, types='ifs', endian=None)
        elif 0: # pramga: no cover
            i = 0
            j = 0
            nbytes = 128
            s1 = Struct(b'< 6i fi 96s')
            while i < len(data):
                datai = data[i:i+nbytes]
                a, b, c, d, e, f, g, h, blank = s1.unpack(datai)
                print('%r %r %r %r %r %r %r %r %r ' % (
                    a, b, c, d, e, f, g, h, blank.rstrip()))
                i += nbytes
                j += 1
            print("j (-5) = %s" % j)
            assert i == len(data), '-5'

        self.read_3_markers([-6, 1, 0])
        self.read_markers([0])

    def read_meff(self):
        """reads the MEFF table"""
        asdf
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_3_markers([-2, 1, 0])
        unused_data = self._read_record()

        for n in [-3, -4, -5, -6, -7, -8]:
            self.read_3_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            unused_data = op2.f.read(nbytes)
            op2.n += nbytes
        n = -9
        self.read_markers([n, 1, 0, 0])

    def read_intmod(self):
        """
        reads the INTMOD table

        tested by TestNastranGUI.test_femap_rougv1_01

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()
        #print('intmod data1')
        #self.show_data(data)

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_3_markers([-2, 1, 0])
        unused_data = self._read_record()
        #print('intmod data2')
        #self.show_data(data)

        for n in [-3, -4, -5, -6, -7, -8,]:
            self.read_3_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            unused_data = op2.f.read(nbytes)
            #print('intmod data%i' % n)
            #self.show_data(data)
            op2.n += nbytes

        n = -9
        self.read_markers([n, 1, 0, 0])
        #self.show(50)
        #raise NotImplementedError(op2.table_name)

    def read_r1tabrg(self):
        """Reads the R1TABRG design response optimization table"""
        op2 = self.op2

        # TODO: I think this is used to handle optimization
        op2._count += 1

        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])

        read_record_ndata = self.get_skip_read_record_ndata()

        # (101, 221355, 0, 0, 0, 0, 0)
        # (???, nvalues,?, ?, ?, ?, ?)
        data = self._read_record()
        values = unpack(mapfmt(self._endian + b'7i', self.size), data)
        unused_nvalues = values[1]
        #print(values)
        #self.show_data(data, types='i', endian=None)

        #'R1TAB   '
        self.read_3_markers([-2, 1, 0])
        data, ndata = read_record_ndata()
        assert ndata == 8 * self.factor, ndata

        itable = -3
        while 1:
            self.read_3_markers([itable, 1, 0])
            stop_marker = self.get_marker1(rewind=True)
            if stop_marker == 0:
                break

            data, ndata = read_record_ndata()
            self._read_r1tabrg(data, ndata)
            itable -= 1
        stop_marker = self.get_marker1(rewind=False)

    def get_skip_read_record_ndata(self):
        """selects the read_record or skip_record depending on read_mode"""
        if self.read_mode == 1:
            read_record_ndata = self._read_record_ndata
        else:
            read_record_ndata = self._skip_record_ndata
        return read_record_ndata

    def _read_r1tabrg(self, data: bytes, ndata: int):
        """
        Design Responses:
          - Weight
          - Flutter Speed
          - Stress
          - Strain
          - Displacement

        """
        op2 = self.op2
        result_name = 'responses'
        if op2._results.is_not_saved(result_name):
            return ndata

        #op2._results._found_result(result_name)
        responses = op2.op2_results.responses

        # create the result object
        if self.read_mode == 1:
            assert data is not None, data
            assert len(data) > 12, len(data)
            if self.size == 4:
                response_type, = op2.struct_i.unpack(data[8:12])
            else:
                response_type, = op2.struct_q.unpack(data[16:24])

            #assert response_type in [1, 6, 10, 84], response_type
            if response_type == 1:
                if responses.weight_response is None:
                    responses.weight_response = WeightResponse()
                else:
                    responses.weight_response.n += 1
            elif response_type == 4:
                #TYPE =4 EIGN or FREQ
                #8 MODE I Mode number
                #9 APRX I Approximation code
                pass
            elif response_type == 5:
                #TYPE =5 DISP
                #8 COMP I Displacement component
                #9 UNDEF None
                #10 GRID I Grid identification number
                if responses.displacement_response is None:
                    responses.displacement_response = DisplacementResponse()
                else:
                    responses.displacement_response.n += 1
            elif response_type == 6:
                if responses.stress_response is None:
                    responses.stress_response = StressResponse()
                else:
                    responses.stress_response.n += 1
            elif response_type == 7:
                if responses.strain_response is None:
                    responses.strain_response = StrainResponse()
                else:
                    responses.strain_response.n += 1

            elif response_type == 8:
                if responses.force_response is None:
                    responses.force_response = ForceResponse()
                else:
                    responses.force_response.n += 1
            elif response_type == 15:
                # CEIG
                #8 MODE I Mode number
                #9 ICODE I 1: Real component or 2: Imaginary component
                pass
            elif response_type == 23:
                if responses.flutter_response is None:
                    responses.fractional_mass_response = FractionalMassResponse()
                else:
                    responses.fractional_mass_response.n += 1
            elif response_type == 84:
                if responses.flutter_response is None:
                    responses.flutter_response = FlutterResponse()
                else:
                    responses.flutter_response.n += 1
            else:
                self.log.warning('skipping response_type=%d' % response_type)
            return ndata
            #else: # response not added...
                #pass

        # fill the result object
        read_r1tabrg = True
        if read_r1tabrg:
            #self.show_data(data, types='ifs', endian=None)
            fmt = self._endian + b'3i 8s 4i i 5i' if self.size == 4 else b'3q 16s 4q q 5q'
            out = unpack(fmt, data)
            # per the R1TAB DMAP page:
            #   all indicies are downshift by 1
            #   indices above out[3] are off by +2 because of the 2 field response_label
            internal_id = out[0]
            dresp_id = out[1]
            response_type = out[2]
            response_label = out[3].strip()
            # -1 for 2 field wide response_label
            region = out[4]
            subcase = out[5]
            type_flag = out[12]  # no meaning per MSC DMAP 2005
            seid = out[13]

            response_label = response_label.decode(self._encoding)

            if response_type == 1:
                responses.weight_response.add_from_op2(out, self.log)
            elif response_type == 5:  # DISP
                # out = (1, 101, 5, 'DISP1   ', 101, 1, 3, 0, 1, 0, 0, 0, 0, 0)
                comp = out[6]
                nid = out[8]
                #msg = f'DISP - label={response_label!r} region={region} subcase={subcase} nid={nid}'
                #print(msg)
                #print(out[6:])
                # (3,   0,  1,    0,   0,   0,   0,   0)
                # (???, NA, comp, ???, ???, ???, ???, ???)
                responses.displacement_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    nid, comp)
            elif response_type == 6:  # STRESS
                #                              -----   STRESS RESPONSES   -----
                #  -------------------------------------------------------------------------------------------
                #   INTERNAL  DRESP1  RESPONSE  ELEMENT   VIEW   COMPONENT  LOWER   INPUT    OUTPUT    UPPER
                #      ID       ID     LABEL       ID    ELM ID     NO.     BOUND   VALUE     VALUE    BOUND
                #  -------------------------------------------------------------------------------------------
                #         21      209  S09L      144747             17       N/A   4.85E+04  5.00E+04  5.00E+04
                # (21, 209, 6, 'S09L    ', 30, 1011, 17, 0, 144747, 0, 0, 0, 0, 0)
                stress_code = out[6]
                pid = out[8]
                #msg = ('STRESS - response_type=%r label=%r region=%s subcase=%s '
                       #'stress_code=%s pid=%s' % (
                           #response_type, response_label, region, subcase,
                           #stress_code, pid))
                responses.stress_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    stress_code, pid)

            elif response_type == 7:  # STRAIN
                strain_code = out[6]
                pid = out[8]
                responses.strain_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    strain_code, pid)
            elif response_type == 8:  # FORCE
                #print('internal_id=%s dresp_id=%s response_type=%s response_label=%s'
                      #' region=%s subcase=%s type_flag=%s seid=%s' % (
                    #internal_id, dresp_id, response_type, response_label,
                    #region, subcase, type_flag, seid
                #))
                force_code = out[6]
                pid = out[8]
                #msg = 'FORCE - label=%r region=%s subcase=%s force_code=%s pid=%s' % (
                    #response_label, region, subcase, force_code, pid)
                #print(msg)
                #print(out)
                responses.force_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    force_code, pid)

            elif response_type == 10:  # CSTRESS
                stress_code = out[6]
                #ply = out[7]
                #pid = out[8]  # is this element id?
                #msg = 'CSTRESS - label=%r region=%s subcase=%s stress_code=%s ply=%s pid=%s' % (
                    #response_label, region, subcase, stress_code, ply, pid)
                #print(msg)
            #elif response_type == 10:  # CSTRAIN
                #pass
            elif response_type == 23:  # fractional mass response
                #msg = f'??? - label={response_label!r} region={region} subcase={subcase}'
                #print(msg)
                #print(out[6:])
                responses.fractional_mass_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid)

            elif response_type == 24:  # FRSTRE
                #8 ICODE I Stress item code
                #9 UNDEF None
                #10 ELID I Element identification number
                #11 FREQ RS Frequency
                #12 IFLAG I Integrated response flag. See Remark 20 of DRESP1.
                #Value is -1 to -6, for SUM, AVG, SSQ,
                pass
            elif response_type == 28:  # RMSACCL
                #8 COMP I RMS Acceleration component
                #9 RANDPS I RANDPS entry identification number
                #10 GRID I Grid identification number
                #11 DMFREQ RS Dummy frequency for internal use
                pass
            elif response_type == 84:
                # FLUTTER  (iii, label, mode, (Ma, V, rho), flutter_id, fff)
                out = unpack(self._endian + b'iii 8s iii fff i fff', data)
                mode = out[6]
                mach = out[7]
                velocity = out[8]
                density = out[9]
                flutter_id = out[10]
                #msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                       #'mach=%s velocity=%s density=%s flutter_id=%s' % (
                           #op2._count, response_label, region, subcase, mode,
                           #mach, velocity, density, flutter_id))
                responses.flutter_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    mode, mach, velocity, density, flutter_id)
                #print(msg)
                #self.log.debug(msg)
            #else:
                #self.log.debug('R1TABRG response response_type=%s not supported' % response_type)
                #raise NotImplementedError(response_type)
            assert len(out) == 14, len(out)
        #self.response1_table[self._count] = out
        return ndata

    def read_sdf(self):
        """reads the SDF table"""
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_3_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 16:
            subtable_name, dummy_a, dummy_b = unpack(self._endian + b'8sii', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n'  % (
                    subtable_name, dummy_a, dummy_b))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
                assert dummy_a == 170, dummy_a
                assert dummy_b == 170, dummy_b
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_3_markers([-3, 1, 1])

        unused_markers0 = self.get_nmarkers(1, rewind=False)
        record = self._read_block()
        self.show_data(record, types='ifs', endian=None, force=False)

        #data = self._read_record()
        self.read_markers([-4, 1, 0])
        return
        self.read_markers([-4, 1, 1])
        record, ndata = self._read_block_ndata()
        #ndata, = op2.struct_i.unpack(record)
        nbytes = ndata * self.size

        print(f'*sdf ndata={ndata} nbytes={nbytes}')
        self.show_data(record, types='ifs', endian=None, force=False)

        self.show_ndata(100, types='ifs', force=False, endian=None)

        #markers1 = self.get_nmarkers4(1, rewind=True)
        #print('markers1', markers1)
        #self.show_ndata(100, types='ifs', force=False, endian=None)
        sys.stdout.flush()
        self.read_markers([0])
        #self._read_subtables()

    def read_tol(self):
        """
        This is probably broken for MSC Nastran

        TOL
        ---
        -2 - nitimes?
        -3 - list of times?

        """
        unused_table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        unused_data = self._read_record()
        #self.show_data(data)

        self.read_3_markers([-2, 1, 0])
        #op2.show_ndata(440, types='if')
        unused_data = self._read_record()
        #print('----')
        self.read_3_markers([-3, 1, 0])
        #op2.show_ndata(440, types='if')
        #print('----')
        self.read_markers([0])
        #data = self._read_record()


        #op2.show_ndata(440, types='ifs')

        #self.show_data(data)

    def _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout):
        """
        +------+---------------------------+
        | Type | Meaning                   |
        +------+---------------------------+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+

        """
        if tout == 1:
            nfloats = nvalues
            nterms = nvalues
            fmt = self._endian + b'i %if' % nfloats
        elif tout == 2:
            nfloats = nvalues // 2
            nterms = nvalues // 2
            fmt = self._endian + b'i %id' % nfloats
        elif tout == 3:
            nfloats = nvalues
            nterms = nvalues // 2
            fmt = self._endian + b'i %if' % nfloats
        elif tout == 4:
            nfloats = nvalues // 2
            nterms = nvalues // 4
            fmt = self._endian + b'i %id' % nfloats
        else:
            raise RuntimeError('tout = %s' % tout)
        return fmt, nfloats, nterms

    def _read_matrix_mat(self):
        """
        Reads a matrix in "standard" form.  The forms are::
            standard:
                Return a matrix that looks similar to a matrix found
                in the OP4.  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``OUTPUT2 KGG//0/12``
            matpool:
                Return a matrix that looks similar to a DMIG matrix
                (e.g., it contains the node id and DOF).  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``TODO: add the magic keyword...``
                ``OUTPUT2 KGG//0/12``

        Matrix Trailer:
        +------+---------------------------------------------------+
        | Word | Contents                                          |
        +======+===================================================+
        |  1   | Number of columns in matrix                       |
        |  2   | Number of rows in matrix                          |
        |  3   | Form of the matrix                                |
        |  4   | Type of matrix                                    |
        |  5   | Largest number of nonzero words among all columns |
        |  6   | Density of the matrix multiplied by 10000         |
        |  7   | Size in blocks                                    |
        |  8   | Maximum string length over all strings            |
        |  9   | Number of strings                                 |
        |  10  | Average bandwidth                                 |
        |  11  | Maximum bandwidth                                 |
        |  12  | Number of null columns                            |
        +------+---------------------------------------------------+

        +------+--------------------------------+
        | Form | Meaning                        |
        +======+================================+
        |  1   | Square                         |
        |  2   | Rectangular                    |
        |  3   | Diagonal                       |
        |  4   | Lower triangular factor        |
        |  5   | Upper triangular factor        |
        |  6   | Symmetric                      |
        |  8   | Identity                       |
        |  9   | Pseudo identity                |
        |  10  | Cholesky factor                |
        |  11  | Trapezoidal factor             |
        |  13  | Sparse lower triangular factor |
        |  15  | Sparse upper triangular factor |
        +------+--------------------------------+

        +------+---------------------------+
        | Type | Meaning                   |
        +======+===========================+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+

        """
        op2 = self.op2
        log = self.log
        allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
        #self.log.debug('----------------------------------------------------------------')
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()

        # old-bad
        #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(self._endian + b'7i', data)

        fmt1 = mapfmt(self._endian + b'7i', self.size)
        #           good   good   good  good  ???    ???
        matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(fmt1, data)
        #print('g =', g)

        utable_name = table_name.decode('utf-8')
        m = Matrix(utable_name, form=form)
        op2.matrices[utable_name] = m

        # matrix_num is a counter (101, 102, 103, ...)
        # 101 will be the first matrix 'A' (matrix_num=101),
        # then we'll read a new matrix 'B' (matrix_num=102),
        # etc.
        #
        # the matrix is Mrows x Ncols
        #
        # it has nvalues in it
        #
        # tout is the precision of the matrix
        # 0 - set precision by cell
        # 1 - real, single precision (float32)
        # 2 - real, double precision (float64)
        # 3 - complex, single precision (complex64)
        # 4 - complex, double precision (complex128)

        # form (bad name)
        # 1 - column matrix
        # 2 - factor matrix
        # 3 - factor matrix
        if tout == 1:
            dtype = 'float32'
        elif tout == 2:
            dtype = 'float64'
        elif tout == 3:
            dtype = 'complex64'
        elif tout == 4:
            dtype = 'complex128'
        else:
            dtype = '???'
            msg = ('unexpected tout for %s: matrix_num=%s form=%s '
                   'mrows=%s ncols=%s tout=%s nvalues=%s g=%s'  % (
                       table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            log.warning(msg)
            raise RuntimeError(msg)

        #self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
        #               'ncols=%s tout=%s nvalues=%s g=%s' % (
        #                   table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
        if form == 1:
            if ncols != mrows:
                log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                    table_name, form, mrows, ncols))
        elif form not in allowed_forms:
            log.error('name=%r matrix_num=%s form=%s mrows=%s '
                      'ncols=%s tout=%s nvalues=%s g=%s' % (
                          table_name, matrix_num, form, mrows, ncols,
                          tout, nvalues, g))
            raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))
        if self.size == 4:
            log.debug('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s '
                      'nvalues=%s g=%s' % (
                          table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
        else:
            #if tout == 1:
                #tout = 2
            log.info('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s '
                     'nvalues=%s g=%s' % (
                         table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        if self.size == 4:
            if len(data) == 16:
                unused_name, ai, bi = unpack(self._endian + b'8s 2i', data)
                assert ai == 170, ai
                assert bi == 170, bi
            else:
                self.log.warning('unexpected matrix length=%s' % len(data))
                #self.log.warning(self.show_data(data, types='if'))
        elif self.size == 8:
            if len(data) == 32:
                unused_name, ai, bi = unpack(self._endian + b'16s 2q', data)
                # name isn't mapped
                assert ai == 170, ai
                assert bi == 170, bi
            else:
                self.log.warning('unexpected matrix length=%s' % len(data))
                #self.log.warning(self.show_data(data, types='ifsqd', endian=None))
        else:
            raise RuntimeError(self.size)

        itable = -3
        unused_j = None

        niter = 0
        niter_max = 100000000

        GCi = []
        GCj = []
        reals = []
        jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)

                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    fmt, unused_nfloats, nterms = self._get_matrix_row_fmt_nterms_nfloats(
                        nvalues, tout)
                    GCjj = [jj] * nterms
                    GCj += GCjj

                    #-----------
                    data = self._read_block()
                    if self.size == 8:
                        #self.log.warning('skipping matrix')
                        #self.show_data(data, types='ifqd')
                        fmt = mapfmt(fmt, self.size)
                        #self.log.warning(fmt)
                        #print('***itable=%s nvalues=%s fmt=%r' % (itable, nvalues, fmt))
                        #continue
                    out = unpack(fmt, data)
                    #print(out)
                    ii = out[0]
                    values = out[1:]

                    #list(range(2, 10))
                    #[2, 3, 4, 5, 6, 7, 8, 9]
                    GCii = list(range(ii, ii + nterms))
                    GCi += GCii
                    reals += values
                    nvalues = self.get_marker1(rewind=True)
                    if self.debug_file:
                        self.binary_debug.write('  GCi = %s\n' % GCii)
                        self.binary_debug.write('  GCj = %s\n' % GCjj)
                        self.binary_debug.write('  reals/imags = %s\n' % str(values))
                assert len(GCi) == len(GCj), 'nGCi=%s nGCj=%s' % (len(GCi), len(GCj))
                if self.size == 4:
                    if tout in [1, 2]:
                        assert len(GCi) == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi), len(reals))
                    else:
                        assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi)*2, len(reals))
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues

                matrix = _cast_matrix_mat(GCi, GCj, mrows, ncols, reals, tout, dtype, log)
                if table_name in DENSE_MATRICES:
                    matrix = matrix.toarray()
                m.data = matrix
                if matrix is not None:
                    op2.matrices[table_name.decode('utf-8')] = m
                #nvalues = self.get_marker1(rewind=True)
                return
            itable -= 1
            niter += 1
        raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

    def _cast_matrix_mat(self, GCi, GCj, mrows, ncols, reals, tout, dtype):
        """helper method for _read_matrix_mat"""
        op2 = self.op2
        #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
        #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)

        # we subtract 1 to the indicides to account for Fortran
        GCi = np.array(GCi, dtype='int32') - 1
        GCj = np.array(GCj, dtype='int32') - 1
        try:
            if dtype == '???':
                matrix = None
                self.log.warning('what is the dtype?')
            elif tout in [1, 2]:
                # real
                real_array = np.array(reals, dtype=dtype)
                matrix = scipy.sparse.coo_matrix(
                    (real_array, (GCi, GCj)),
                    shape=(mrows, ncols), dtype=dtype)
                #self.log.info('created %s (real)' % self.table_name)
            elif tout in [3, 4]:
                # complex
                real_array = np.array(reals, dtype=dtype)
                nvalues_matrix = real_array.shape[0] // 2
                real_complex = real_array.reshape((nvalues_matrix, 2))
                real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
                if self.binary_debug:
                    #self.binary_debug.write('reals = %s' % real_complex[:, 0])
                    #self.binary_debug.write('imags = %s' % real_complex[:, 1])
                    self.binary_debug.write('real_imag = %s' % real_imag)
                matrix = scipy.sparse.coo_matrix(
                    (real_imag, (GCi, GCj)),
                    shape=(mrows, ncols), dtype=dtype)
                #msg = 'created %s (complex)' % self.table_name
                #self.log.debug(msg)
                #raise RuntimeError(msg)
            else:
                raise RuntimeError('this should never happen')
        except ValueError:
            self.log.warning('shape=(%s, %s)' % (mrows, ncols))
            self.log.warning('cant make a coo/sparse matrix...trying dense')

            if dtype == '???':
                matrix = None
                self.log.warning('what is the dtype?')
            else:
                real_array = np.array(reals, dtype=dtype)
                self.log.debug('shape=%s mrows=%s ncols=%s' % (
                    str(real_array.shape), mrows, ncols))
                if len(reals) == mrows * ncols:
                    real_array = real_array.reshape(mrows, ncols)
                    self.log.info('created %s' % op2.table_name)
                else:
                    self.log.warning('cant reshape because invalid sizes : created %s' %
                                     op2.table_name)

                matrix = real_array
        return matrix

    def _skip_matrix_mat(self):
        """
        Reads a matrix in "standard" form.

        Notes
        -----
        see read_matrix_mat

        """
        unused_table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        unused_data = self._skip_record()

        self.read_3_markers([-2, 1, 0])
        unused_data = self._skip_record()

        itable = -3
        niter = 0
        niter_max = 100000000

        #jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            #print('nvalues4a =', nvalues)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)
                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    unused_data = self._skip_block()
                    nvalues = self.get_marker1(rewind=True)
                #jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                return
            itable -= 1
            niter += 1
        raise RuntimeError('this should never happen; n=%s' % niter_max)

    def read_matrix(self, table_name):
        """
        General method for reading matrices and MATPOOL matrices

        Note
        ----
        Matrices are read on read_mode = 1

        .. todo:: Doesn't support checking matrices vs. MATPOOLs

        """
        read_mode_to_read_matrix = 1
        op2 = self.op2
        i = op2.f.tell()
        # if we skip on read_mode=1, we don't get debugging
        # if we just use read_mode=2, some tests fail
        #
        if self.read_mode != read_mode_to_read_matrix and not self.debug_file:
            try:
                self._skip_matrix_mat()  # doesn't work for matpools
            except MemoryError:
                raise
            except(RuntimeError, AssertionError, ValueError):
                self._goto(i)
                self._skip_table(table_name)
            return

        try:
            self._read_matrix_mat()
        except MemoryError:
            raise
        except(RuntimeError, AssertionError, ValueError):
            # read matpool matrix
            self._goto(i)
            try:
                self._read_matrix_matpool()
            except(RuntimeError, AssertionError, ValueError):
                raise
                #self._goto(i)
                #self._skip_table(op2.table_name)

    def _read_mef1(self):
        op2 = self.op2
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        #print(utable_name)
        self.read_markers([-1])

        #(101, 15, 12, 2, 2, 4, 1333)
        data = self._read_record()
        ints = unpack(self._endian + b'7i', data)
        print(ints)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        name, int1, int2 = unpack(self._endian + b'8s ii', data)
        print(name, int1, int2)
        #self.show_data(data, types='ifs', endian=None, force=False)

        self.read_3_markers([-3, 1])
        read_dofs(op2, size=4)
        #data = self._read_record()
        #(4, 2, 4, 12, 1, -901, -1074790401, 12,
         #4, 2, 4, 12, 7, -901, -1074790401, 12)
        #self.show_data(data, types='ifsqd', endian=None, force=False)
        #self.show_data(data[4:], types='ifsqd', endian=None, force=False)

        assert op2.f.tell() == op2.n

        ni = -4
        marker = 1
        while marker != 0:
            #print(f'ni = {ni}')
            self.read_3_markers([ni, 1])
            read_dofs(op2, size=4)
            ni -= 1

            marker = self.get_marker1(rewind=True, macro_rewind=False)
            #if marker == 0:
                #break
        #print('marker =', marker)
        marker = self.get_marker1(rewind=False, macro_rewind=False)
        #self.read_3_markers([ni, 1])
        #self.show_ndata(200, types='ifs')

        #self.read_3_markers([-9, 1])
        #read_dofs(op2, size=4)

    # (104, 32768, 0, 0, 0, 0, 0)
        #aaa

    def _read_matrix_matpool(self):
        """
        Reads a MATPOOL matrix

        MATPOOL matrices are always sparse

        +------+-----------------+
        | Form | Meaning         |
        +======+=================+
        |  1   | Square          |
        |  2   | Rectangular     |
        |  6   | Symmetric       |
        |  9   | Pseudo identity |
        +------+-----------------+


        Record 2  - BNDFL(9614,96,0)
        Record 3  - DMIAX(214,2,221)
        Record 4  - DMIG(114,1,120)
        Record 5  - DMIJ(514,5,578)
        Record 6  - DMIJI(614,6,579)
        Record 7  - DMIK(714,7,580)
        Record 8  - ELIST(314,3,279)
        Record 9  - MFLUID(414,4,284)
        Record 10 - RADCAV(2509,25,418)
        Record 11 - RADSET(8602,86,421)
        Record 12 - RADLST(2014,20,243)
        Record 13 - RADMTX(3014,30,244)
        """
        #print('-------------------------------------')

        op2 = self.op2
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        #print(utable_name)
        self.read_markers([-1])

        # (104, 32768, 0, 0, 0, 0, 0)
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None, force=False)

        self.read_3_markers([-2, 1, 0])

        # MATPOOL
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None, force=False)
        if self.size == 8:
            data = reshape_bytes_block(data)
        #self.show_data(data)
        ndata = len(data)
        if ndata == 8:
            table_name2, = op2.struct_8s.unpack(data)
            utable_name = table_name2.decode('utf-8').strip()
            #assert utable_name == utable_name2, f'utable_name={utable_name} utable_name2={utable_name2}'

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        if self.size == 4:
            struct_3i = op2.struct_3i
        else:
            struct_3i = op2.struct_3q
        itable = -4
        while 1:
            #self.show_data(data[:12*self.factor], types='ifsqd')
            code = struct_3i.unpack(data[:12*self.factor])
            #self.show(36)

            #if utable_name == 'DELTAK':
                #self.show_data(data)

            #nvalues = len(data) // 4
            assert len(data) % 4 == 0, len(data) / 4

            if self.read_mode == 1:
                self.read_matpool_result(code, op2, data, utable_name)

            self.read_3_markers([itable, 1, 0])
            #self.log.debug(f'  read [{itable},1,0]')
            expected_marker = itable - 1
            data, ndatas = self.read_long_block(expected_marker)
            if ndatas == 0:
                itable -= 1
                #self.log.debug(f'  read [{itable},1,0]')
                self.read_3_markers([itable, 1, 0])
                break
            elif self.read_mode == 1:
                #self.show_data(data, types='ifs', endian=None, force=False)
                code = struct_3i.unpack(data[:12*self.factor])
                self.read_matpool_result(code, op2, data, utable_name)
                #self.log.info('showing data...')
                #data, ndatas = self.read_long_block(expected_marker)
            itable -= 1

        #self.show(100, types='ifsqd')
        self.read_markers([0])
        #raise RuntimeError('failed on _read_matpool_matrix')

    def read_long_block(self, expected_marker: int):
        op2 = self.op2
        if self.size == 4:
            struct_3i = op2.struct_3i
        else:
            struct_3i = op2.struct_3q

        marker1 = None
        datas = []
        ndatas = 0
        while marker1 != expected_marker:
            data, ndata = self._read_record_ndata()
            if ndata == 12 * self.factor:
                out = struct_3i.unpack(data)
                #print('****', out)
                #self.show(100, types='ifsq')
                marker1 = self.get_marker1(rewind=True, macro_rewind=False)
                #self.log.info(f'marker1 = {marker1}')
                continue
            #self.log.info('adding data')
            datas.append(data)
            ndatas += ndata
            #self.show(500, types='ifsq')
            marker1 = self.get_marker1(rewind=True, macro_rewind=False)
            #self.log.info(f'marker1 = {marker1}')

        if len(datas) != 1:
            #print('ndatas =', len(datas))
            data = b''.join(datas)
        #self.show(200, types='ifsq')
        return data, ndatas


    def read_matpool_result(self, code, op2, data, utable_name):
        if code == (114, 1, 120):
            self.log.debug(f'  code = {code}')
            try:
                if self.size == 4:
                    read_matpool_dmig_4(op2, data, utable_name, debug=False)
                else:
                    read_matpool_dmig_8(op2, data, utable_name, debug=False)
            #try:
                #read_matpool_dmig(op2, data, utable_name, debug=False)
            except Exception as excep:
                self.log.error(str(excep))
                self.log.warning('  skipping MATPOOL-DMIG')
                #raise
        elif code == (314, 3, 279):
            # geom
            self._read_matpool_elist(op2, data, utable_name, debug=False)
        elif code == (414, 4, 284):
            # geom
            self._read_matpool_mfluid(op2, data, utable_name, debug=False)
        elif code == (2509, 25, 418):
            # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
            self.log.warning('  skipping MATPOOL-RADCAV')
        elif code == (3014, 30, 244):
            # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
            self.log.warning('  skipping MATPOOL-RADMTX')
        elif code == (8602, 86, 421):
            # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15305.op2
            self.log.warning('  skipping MATPOOL-RADSET')
        elif code == (2014, 20, 243):
            # C:\NASA\m4\formats\git\examples\pyNastran_examples\demo_sort2_post_m2\hd15306.op2
            self.log.warning('  skipping MATPOOL-RADLST')

        elif code == (9614, 96, 0):
            # some axisymmetric matrix
            self._read_matpool_bndfl(op2, data, utable_name, debug=False)
        else:
            print(f'  code = {code}')
            self.show_data(data, types='ifs', endian=None, force=False)
            raise NotImplementedError(code)

    def _read_matpool_mfluid(self, op2: OP2, data, utable_name: str, debug: bool=False):
        r"""
        Word Name Type Description
        1 SID       I
        2 CID       I
        3 ZFR      RS
        4 RHO      RS
        5 ELIST1    I
        6 ELIST2    I
        7 PLANE1    I
        8 PLANE2    I
        9 RMAX     RS
        10 FMEXACT RS

        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\cms01.op2
        """
        assert len(data) == 12 + 40, len(data)
        structi = Struct(self._endian + b'2i 2f 4i 2f')
        (unused_sid, unused_cid, unused_zfr, unused_rho, unused_elist1, unused_elist2,
         unused_plane1, unused_plane2, unused_rmax, unused_fmexact) = structi.unpack(data[12:])
        self.log.warning('skipping MATPOOL-MFLUID table')

    def _read_matpool_elist(self, op2: OP2, data: bytes, utable_name: str, debug: bool=False):
        """
        Word Name Type Description
        1 LID I
        2 E1 I
        Word 2 repeats until End of Record
        """
        n = 12 * self.factor
        datai = data[n:]
        ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
        assert ints[-1] == 0, ints
        elist_id = ints[0]
        element_ids = ints[1:-1].tolist()
        if not hasattr(self, 'elist'):
            self.elist = {}
        self.elist[elist_id] = element_ids

    def _read_matpool_bndfl(self, op2: OP2, data: bytes, utable_name: str, debug: bool=False):
        r"""
        Word Name Type Description
        1 CSF      I
        2 G        RS
        3 RHO      RS
        4 B        RS
        5 NOSYM    I
        6 M        I
        7 S1       I
        8 S2       I
        9 NHARM(C) I
        10 NI      I
        Word 10 repeats NHARM times
        11 IDFL I
        12 R    RS
        13 Z    RS
        14 L    RS
        15 C    RS
        16 S    RS
        17 RHOI RS
        18 G    I
        19 PHI  RS
        Words 18 through 19 repeat until (-1,-1) occurs
        Words 11 through 19 repeat until End of Record

        ints = [7, 26, 27, 43, 44, 60, 61]
        C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1072x.op2
        """
        endian = self._endian
        struct1 = Struct(endian + b'i 3f 5i')
        #          c   g    rho   b   nosym m s1  s2 nharm n1, n2
        #ints    = (2, 32.2, 0.03, 0.0, 0,   4, 1, -1, 2,    4, 8,
        #           # idfl r    z     l    c    s    rho
        #           2,     8.0, 10.0, 2.5, 1.0, 0,   0.03,
        #           # g  phi
        #           3, 0,
        #           4, 30.0,
        #           5, 60.0,
        #           6, 90.0,
        #           -1, -1,
        #
        #           # idfl  r    z    l      c     s    rho
        #           8,      8.0, 5.0, 5.0, 0.93, -0.35, 0.03,
        #           # g phi
        #           9, 0,
        #           10, 30.0,
        #           11, 60.0,
        #           12, 90.0,
        #           -1, -1,
        #           ...)

        size = self.size
        factor = self.factor
        log = self.log
        ndata = len(data)
        log.warning('skipping MATPOOL-BNDFL table')
        n = 12 * factor
        datai = data[n:]
        ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
        floats = np.frombuffer(datai, dtype=op2.fdtype8).copy()
        iminus1 = np.where(ints == -1)[0]
        #print(iminus1)

        iharm = 9
        nharm = ints[iharm]
        #idfl = ints[iharm + nharm]
        #print(f'nharm = {nharm}')

        #b_idfl = (iharm + nharm) * size + 12
        #self.show_data(datai, types='if')
        #self.show_data(datai[:b_idfl], types='ifs')

        ntotal1 = 9 * size
        #print('i =', i)

        edata = data[n:n+ntotal1]
        #print('len(edata)', len(edata))
        c, g, rho, b, nosym, m, s1, s2, nharm = struct1.unpack(edata)
        log.debug(f'  c={c} g={g:g} rho={rho:g} nosym={nosym} m={m} s1={s1} s2={s2} nharm={nharm}')
        n += ntotal1
        i = 9

        ntotal2 = nharm * size
        #ni = Struct(b'%di' % nharm).unpack(data[n:n+ntotal2])
        ni = ints[i:i+nharm]
        #print(f'n={n} i={i} ni={ni}')
        n += ntotal2
        i += nharm

        while n < ndata:
            #print('------------------------------------')
            #print(ints[i-2:i+10])
            #print(floats[i-2:i+10])
            # 11 IDFL I
            # 12 R    RS
            # 13 Z    RS
            # 14 L    RS
            # 15 C    RS
            # 16 S    RS
            # 17 RHOI RS
            idfl = ints[i]
            r, z, l, c, s, rhoi = floats[i+1:i+7]
            log.debug(f'    idfl={idfl} r={r} z={z} L={l:g} c={c:g} s={s:g} rhoi={rhoi:g}')
            i += 7
            n += 7 * size
            assert idfl < 1000

            ints2 = ints[i+1::2]
            ints1 = ints[i::2][:len(ints2)]
            iminus1i = np.where((ints1 == -1) & (ints2 == -1))[0][0] - 1

            #intsi = ints[i:i+iminus1i+4:2]  # good
            intsi = ints1[:iminus1i+1]
            #floatsi = floats[i+1::2][:iminus1i+1]  #good
            floatsi = floats[i+1:i+2*(iminus1i+1):2]
            log.debug(f'    {intsi} {floatsi}')
            assert intsi.min() >= 1, intsi
            assert intsi.max() <= 1000, intsi
            #n += (iminus1i + 4) * size# ???
            #i += iminus1i + 1
            i += len(intsi) * 2 + 2 # ???
            n2 = (i + 3) * size
            n = n2
            #print(f'n={n} i={i} -> n2={n2}')
        return


    def _read_units(self):
        r"""models/msc/units_mass_spring_damper"""
        op2 = self.op2
        assert self.factor == 1, '64-bit is not supported'

        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])

        #read_record_ndata = self.get_skip_read_record_ndata()

        #(101, 32767, 32767, 32767, 32767, 32767, 32767)
        data = self._read_record()
        #self.show_data(data, types='ifs', endian=None, force=False)

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        assert data == b'UNITS   ', data

        self.read_3_markers([-3, 1, 0])
        data = self._read_record()

        # per MSC DMAP 2020
        #
        # Word Name Type Description
        # 1 MASS(2)   CHAR4 Units assumed for mass
        # 3 FORCE(2)  CHAR4 Units assumed for force
        # 5 LENGTH(2) CHAR4 Units assumed for length
        # 7 TIME(2)   CHAR4 Units for assumed time
        # 9 STRESS(2) CHAR4 Units for assumed stress

        ndata = len(data)
        if op2.is_geometry:
            if ndata == 32:
                #'MGG     N       MM      S       '
                out = Struct(self._endian + b'8s 8s 8s 8s').unpack(data)
                mass_bytes, force_bytes, length_bytes, time_bytes = out
                mass = mass_bytes.decode(self._encoding)
                force = force_bytes.decode(self._encoding)
                length = length_bytes.decode(self._encoding)
                time = time_bytes.decode(self._encoding)
                print(mass, force, length, time)
                fields = {
                    'mass' : mass,
                    'force' : force,
                    'length' : length,
                    'time' : time, }
                op2.add_dti('UNITS', fields)
            else:
                raise RuntimeError(f'ndata={len(data)} (expected 40); data={data!r}')
        self.read_3_markers([-4, 1, 0])
        self.read_markers([0])
    #---------------------------------------------------------------------------

    def _get_marker_n(self, nmarkers: int) -> list[int]:
        """
        Gets N markers

        A marker is a flag that is used.  It's a series of 3 ints (4, n, 4)
        where n changes from marker to marker.

        Parameters
        ----------
        nmarkers : int
            the number of markers to read

        Returns
        -------
        markers : list[int, int, int]
            a list of nmarker integers

        """
        op2 = self.op2
        markers = []
        struc = Struct('3i')
        for unused_i in range(nmarkers):
            block = op2.f.read(12)
            marker = struc.unpack(block)
            markers.append(marker)
        return markers

    def get_nmarkers(self, n: int, rewind=True, macro_rewind=False):
        if self.size == 4:
            return self.get_nmarkers4(n, rewind=rewind, macro_rewind=macro_rewind)
        return self.get_nmarkers8(n, rewind=rewind, macro_rewind=macro_rewind)

    def get_nmarkers4(self, n: int, rewind=True, macro_rewind=False):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        Parameters
        ----------
        n : int
            number of markers to get
        rewind : bool
            should the file be returned to the starting point

        Returns
        -------
        markers : list[int]
            list of [1, 2, 3, ...] markers

        """
        op2 = self.op2
        ni = op2.n
        markers = []
        for i in range(n):
            data = self.read_block4()
            marker, = op2.struct_i.unpack(data)
            markers.append(marker)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
            #for i in range(n):
                #self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (
                    #i, macro_rewind or rewind))
        else:
            #if not macro_rewind:
            if self.is_debug_file:
                for i in range(n):
                    self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (
                        i, macro_rewind or rewind))
        return markers

    def get_nmarkers8(self, n: int, rewind=True, macro_rewind=False):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        Parameters
        ----------
        n : int
            number of markers to get
        rewind : bool
            should the file be returned to the starting point

        Returns
        -------
        markers : list[int]
            list of [1, 2, 3, ...] markers

        """
        op2 = self.op2
        ni = op2.n
        markers = []
        for i in range(n):
            data = self.read_block8()
            marker, = op2.struct_q.unpack(data)
            markers.append(marker)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
            #for i in range(n):
                #self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (
                    #i, macro_rewind or rewind))
        else:
            #if not macro_rewind:
            if self.is_debug_file:
                for i in range(n):
                    self.binary_debug.write('get_nmarkers- [8, %i, 8]; macro_rewind=%s\n' % (
                        i, macro_rewind or rewind))
        return markers

    def read_markers(self, markers: list[int], macro_rewind: bool=True) -> None:
        if self.size == 4:
            return self.read_markers4(markers, macro_rewind=macro_rewind)
        return self.read_markers8(markers, macro_rewind=macro_rewind)

    def read_markers4(self, markers: list[int], macro_rewind: bool=True) -> None:
        """
        Gets specified markers, where a marker has the form of [4, value, 4].
        The "marker" corresponds to the value, so 3 markers takes up 9 integers.
        These are used to indicate position in the file as well as the number
        of bytes to read.

        Because we're checking the markers vs. what we expect, we just throw
        the data away.

        Parameters
        ----------
        markers : list[int]
            markers to get; markers = [-10, 1]

        Raises
        ------
        FortranMarkerError
            if the expected table number is not found

        """
        op2 = self.op2
        for i, marker in enumerate(markers):
            #self.log.debug('markers[%i] = %s' % (i, marker))
            data = self.read_block4()
            imarker, = op2.struct_i.unpack(data)
            if marker != imarker:
                #self.show_data(data)
                msg = 'marker=%r imarker=%r; markers=%s; i=%s; table_name=%r; iloc=%s/%s' % (
                    marker, imarker, markers, i, op2.table_name,
                    op2.f.tell(), os.path.getsize(op2.op2_filename))
                raise FortranMarkerError(msg)
            if self.is_debug_file:
                self.binary_debug.write(f'  read_markers -> [4, {marker:d}, 4]\n')

    def read_markers8(self, markers: list[int], macro_rewind: int=True) -> None:
        """
        Gets specified markers, where a marker has the form of [4, value, 4].
        The "marker" corresponds to the value, so 3 markers takes up 9 integers.
        These are used to indicate position in the file as well as the number
        of bytes to read.

        Because we're checking the markers vs. what we expect, we just throw
        the data away.

        Parameters
        ----------
        markers : list[int]
            markers to get; markers = [-10, 1]

        Raises
        ------
        FortranMarkerError
            if the expected table number is not found

        """
        op2 = self.op2
        for i, marker in enumerate(markers):
            #self.log.debug('markers[%i] = %s' % (i, marker))
            data = self.read_block8()
            imarker, = op2.struct_q.unpack(data)
            if marker != imarker:
                #self.show_data(data, types='iq')
                msg = 'marker=%r imarker=%r; markers=%s; i=%s; table_name=%r; iloc=%s/%s' % (
                    marker, imarker, markers, i, op2.table_name,
                    op2.f.tell(), os.path.getsize(op2.op2_filename))
                raise FortranMarkerError(msg)
            if self.is_debug_file:
                self.binary_debug.write(f'  read_markers -> [8, {marker:d}, 8]\n')

    def _skip_table(self, table_name: str, warn: bool=True) -> None:
        """bypasses the next table as quickly as possible"""
        #if table_name in ['DIT', 'DITS']:  # tables
            #self._read_dit()
        if table_name in ['PCOMPTS', 'PCOMPTS']:
            self._read_pcompts()
        else:
            self._skip_table_helper(warn=warn)

    def _print_month(self, month: int, day: int, year: int, zero: int, one: int) -> None:
        """
        Creates the self.date attribute from the 2-digit year.

        Parameters
        ----------
        month : int
            the month (integer <= 12)
        day :  int
            the day (integer <= 31)
        year : int
            the day (integer <= 99)
        zero : int
            a dummy integer (???)
        one : int
            a dummy integer (???)

        """
        if year > 2000:
            # would you believe they changed the date format?
            month, day = day, month
        else:
            year += 2000
        month, day, year = self._set_op2_date(month, day, year)

        #self.log.debug("%s/%s/%4i zero=%s one=%s" % (month, day, year, zero, one))
        #if self.is_debug_file:
        if self.is_debug_file:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, '
                                    'zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        #assert zero == 0, zero  # is this the RTABLE indicator???
        assert one in [0, 1], one  # 0, 50

    def _set_op2_date(self, month: int, day: int, year: int) -> tuple[int, int, int]:
        """sets the date the job was run"""
        date = (month, day, year)
        self.op2.date = date
        return date

    #----------------------------------------------------------------------------------------
    def _read_record(self, debug=True, macro_rewind=False) -> bytes:
        """
        Reads a record.

        A record is defined N blocks.  Blocks are split every 2^12 bytes,
        which is an oddity of the OP2, which is a "Fortran formatted" file.
        You can think of a block as a partial result.  A record is a full
        result.

        If a block is small enough, it will fit into 2^12 bytes and a record
        is a block.

        """
        if self.size == 4:
            return self._read_record_ndata4(debug=debug, macro_rewind=macro_rewind)[0]
        return self._read_record_ndata8(debug=debug, macro_rewind=macro_rewind)[0]

    def _read_record_ndata(self, debug=True, macro_rewind=False) -> tuple[bytes, int]:
        """reads a record and the length of the record"""
        if self.size == 4:
            return self._read_record_ndata4(debug=debug, macro_rewind=macro_rewind)
        return self._read_record_ndata8(debug=debug, macro_rewind=macro_rewind)

    def _read_record_ndata4(self, debug=True, macro_rewind=False) -> tuple[bytes, int]:
        """reads a record and the length of the record for size=4"""
        op2 = self.op2
        marker0 = self.get_marker1_4(rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (
                marker0, macro_rewind))
        na = op2.n
        record, nrecord = self._read_block_ndata4()

        if self.is_debug_file and debug:
            msg = 'read_record - record = [%i, recordi, %i]; macro_rewind=%s\n' % (
                nrecord, nrecord, macro_rewind)
            self.binary_debug.write(msg)

        if marker0*4 != len(record):
            op2.f.seek(na)
            op2.n = na
            if nrecord == 4:
                self.log.error(f'EmptyRecordError: marker0={marker0} nrecord={nrecord}')
                raise EmptyRecordError('nrecord=4')
            self.log.error(f'marker0={marker0} nrecord={nrecord}')
            raise FortranMarkerError('marker0=%s*4 len(record)=%s; table_name=%r' % (
                marker0*4, len(record), op2.table_name))

        #markers1 = self.get_marker1_4(rewind=True)
        markers1 = self.get_nmarkers4(1, rewind=True)
        if markers1[0] > 0:
            #nloop = 0
            records = [record]
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                recordi, nrecordi = self._read_block_ndata()
                nrecord += nrecordi
                records.append(recordi)
                #record += recordi
                markers1 = self.get_nmarkers(1, rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                #nloop += 1

            # if nloop == 0:
                # record = records[0]
            # elif nloop == 1:
                # record = records[0] + records[1]
            # else:
            record = b''.join(records)
        return record, nrecord

    def _read_record_ndata8(self, debug=True, macro_rewind=False) -> tuple[bytes, int]:
        """reads a record and the length of the record for size=8"""
        op2 = self.op2
        markers0 = self.get_nmarkers8(1, rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [8, %i, 8]; macro_rewind=%s\n' % (
                markers0[0], macro_rewind))
        record, nrecord = self._read_block_ndata8()

        if self.is_debug_file and debug:
            msg = 'read_record - record = [%i, recordi, %i]; macro_rewind=%s\n' % (
                nrecord, nrecord, macro_rewind)
            self.binary_debug.write(msg)
        if markers0[0]*8 != len(record):
            self.show_data(record, types='ifsqd')
            raise FortranMarkerError('markers0=%s*8 len(record)=%s; table_name=%r' % (
                markers0[0]*8, len(record), op2.table_name))

        markers1 = self.get_nmarkers8(1, rewind=True)
        if markers1[0] > 0:
            #nloop = 0
            records = [record]
            while markers1[0] > 0:
                markers1 = self.get_nmarkers8(1, rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [8, %i, 8]\n' % markers1[0])
                recordi, nrecordi = self._read_block_ndata8()
                nrecord += nrecordi
                records.append(recordi)
                #record += recordi
                markers1 = self.get_nmarkers8(1, rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [8, %i, 8]\n' % markers1[0])
                #nloop += 1

            # if nloop == 0:
                # record = records[0]
            # elif nloop == 1:
                # record = records[0] + records[1]
            # else:
            record = b''.join(records)
        #self.show(200, types='ifq')
            #if out == (65535, 65535, 65535):
                #continue
                #self.show_data(data, types='ifsqd')
                #self.log.info('break....')
                #break

        return record, nrecord

    def _read_block_ndata4(self) -> tuple[bytes, int]:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary
        ndata : int
            len(data)

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = op2.f.read(ndata)
        data = op2.f.read(4)
        op2.n += 8 + ndata
        return data_out, ndata

    def _read_block_ndata(self) -> tuple[bytes, int]:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary
        ndata : int
            len(data)

        """
        if self.size == 4:
            return self._read_block_ndata4()
        return self._read_block_ndata8()

    def _read_block_ndata8(self) -> tuple[bytes, int]:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary
        ndata : int
            len(data)

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = op2.f.read(ndata)
        data = op2.f.read(4)
        op2.n += 8 + ndata
        return data_out, ndata
    #------------------------------------------------------------------
    def unpack_table_name(self, data: bytes) -> bytes:
        if self.size == 4:
            return self.unpack_table_name4(data)
        return self.unpack_table_name8(data)

    def unpack_table_name4(self, data: bytes) -> bytes:
        """table names can apparently be 8 or 32 characters"""
        if len(data) == 8:
            # 'GEOM4   '
            structi = self.op2.struct_8s
            table_name, = structi.unpack(data)
        elif len(data) == 32:
            # 'GEOM1   20140   0   _y\xfe\xffGEOM1   '
            # ['GEOM1   ', '20140   ', '0   ', -10000, ' GEOM1   ']
            structi = self.op2.struct_8s
            table_name, = structi.unpack(data[:8])
        else:
            raise NotImplementedError(f'data={data!r}; n={len(data)}')
        return table_name.strip()

    def unpack_table_name8(self, data: bytes) -> bytes:
        """table names can apparently be 8 or 32 characters"""
        if len(data) == 8:
            # 'GEOM4   '
            structi = self.op2.struct_8s
            table_name, = structi.unpack(data)
        elif len(data) == 16:
            # 'GEOM4   '
            #b'PVT0            '
            structi = self.op2.struct_16s
            table_name, = structi.unpack(data)
        #elif len(data) == 32:
            # 'GEOM1   20140   0   _y\xfe\xffGEOM1   '
            # ['GEOM1   ', '20140   ', '0   ', -10000, ' GEOM1   ']
            #structi = self.op2.struct_8s
            #table_name, = structi.unpack(data[:8])
        else:
            raise NotImplementedError(f'data={data!r}; n={len(data)}')
        return table_name.strip()

    def _read_table_name(self, last_table_name: Optional[bytes]=None,
                         rewind: bool=False, stop_on_failure: bool=True) -> bytes:
        """
        Reads the next OP2 table name (e.g. OUG1, OES1X1)

        Parameters
        ----------
        last_table_name : bytes; default=Noen
            the last table name

        Returns
        -------
        table_name : bytes
            the table name

        """
        op2 = self.op2
        table_name = None
        data = None
        if self.is_debug_file:
            self.binary_debug.write('_read_table_name - rewind=%s\n' % rewind)
        ni = op2.n
        if stop_on_failure:
            data = self._read_record(debug=False, macro_rewind=rewind)
            if self.size == 8:
                data = reshape_bytes_block(data)
            table_name = self.unpack_table_name(data)

            if self.is_debug_file and not rewind:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write('table_header = [8, %r, 8]\n\n' % table_name)
            table_name = table_name.strip()
        else:
            try:
                data = self._read_record(macro_rewind=rewind)
                if self.size == 8:
                    data = reshape_bytes_block(data)
                table_name = self.unpack_table_name(data)
            except (NameError, MemoryError):
                raise
            except Exception: # struct_error:
                # we're done reading
                op2.n = ni
                op2.f.seek(op2.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0], macro_rewind=rewind)
                except Exception: #struct_error:
                    # if we hit this block, we have a FATAL error
                    is_special_nastran = op2._nastran_format.lower().startswith(('imat', 'autodesk'))
                    if not is_special_nastran and op2.post != -4:
                        op2.f.seek(op2.n)
                        self.show(1000)
                        if last_table_name:
                            self.log.error(f'finished table_name = {last_table_name}')
                        raise FatalError('There was a Nastran FATAL Error.  Check the F06.\n'
                                         f'last table={op2.table_name!r}; post={op2.post} '
                                         f'version={self.op2._nastran_format!r}')
                table_name = None

                # we're done reading, so we're going to ignore the rewind
                rewind = False

        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
        if table_name is not None:
            assert len(table_name) <= 8, f'table_name={table_name}; n={len(table_name)}'
        #if self.size == 8:
            #print(table_name)
        return table_name

    def _read_block(self):
        if self.size == 4:
            return self.read_block4()
        return self.read_block8()

    def read_string_block(self) -> bytes:
        block = self._read_block()
        if self.size == 4:
            return block
        return reshape_bytes_block(block)

    def read_block4(self) -> bytes:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary

        Notes
        -----
        see read_3_blocks

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = op2.f.read(ndata)
        data = op2.f.read(4)
        op2.n += 8 + ndata
        return data_out

    def read_block8(self) -> bytes:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary

        Notes
        -----
        see read_3_blocks

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)
        data_out = op2.f.read(ndata)
        data = op2.f.read(4)
        op2.n += 8 + ndata
        return data_out

    def read_3_blocks4(self) -> bytes:
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]
            [nbytes, data, nbytes]
            [nbytes, data, nbytes]

        This is intended to be used for reading marker triples

        Returns
        -------
        data : bytes
            the data in binary

        """
        op2 = self.op2
        data_out = b''
        for unused_i in range(3):
            data = op2.f.read(4)
            ndata, = op2.struct_i.unpack(data)

            data_out += op2.f.read(ndata)
            data = op2.f.read(4)
            op2.n += 8 + ndata
        return data_out

    def read_3_markers(self, markers, macro_rewind=True) -> None:
        """Micro-optimizes ``read_markers`` for 3 markers."""
        if self.size == 4:
            self.read_3_markers4(markers, macro_rewind=macro_rewind)
        else:
            self.read_markers8(markers, macro_rewind=macro_rewind)

    def read_3_markers4(self, markers, macro_rewind=True) -> None:
        """
        Micro-optimizes ``read_markers`` for 3 markers.

        Parameters
        ----------
        markers : list[int, int, int]
            markers to get; markers = [-10, 1, 0]

        Raises
        ------
        FortranMarkerError
            if the expected table number is not found

        """
        #data1 = self.read_block4()
        #data2 = self.read_block4()
        #data3 = self.read_block4()
        #data = data1 + data2 + data3
        op2 = self.op2
        data = self.read_3_blocks4()
        markers_actual = op2.struct_3i.unpack(data)
        for imarker, marker, marker_actual in zip(count(), markers, markers_actual):
            if marker != marker_actual:
                raise FortranMarkerError(f'imarker={imarker}; markers={markers}; '
                                         f'marker_actual={markers_actual} table_name={op2.table_name!r}')
            if self.is_debug_file:
                self.binary_debug.write('  read_markers -> [4, %i, 4]\n' % marker)

    def get_marker1(self, rewind: bool=True, macro_rewind: bool=False) -> int:
        if self.size == 4:
            return self.get_marker1_4(rewind=rewind, macro_rewind=macro_rewind)
        return self.get_marker1_8(rewind=rewind, macro_rewind=macro_rewind)

    def get_marker1_4(self, rewind: bool=True, macro_rewind: bool=False) -> int:
        """
        Gets 1 marker
        See get_n_markers(...)

        Parameters
        ----------
        rewind : bool
            should the file be returned to the starting point
        macro_rewind : bool
            ???

        Returns
        -------
        markers : int
            a single marker

        """
        op2 = self.op2
        ni = op2.n
        data = self.read_block4()
        marker, = op2.struct_i.unpack(data)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
        else:
            if self.is_debug_file:
                msg = 'get_marker - [4, %i, 4]; macro_rewind=%s\n' % (
                    marker, macro_rewind or rewind)
                self.binary_debug.write(msg)
        return marker

    def get_marker1_8(self, rewind=True, macro_rewind=False) -> int:
        """
        Gets 1 marker
        See get_n_markers(...)

        Parameters
        ----------
        rewind : bool
            should the file be returned to the starting point
        macro_rewind : bool
            ???

        Returns
        -------
        markers : int
            a single marker

        """
        op2 = self.op2
        ni = op2.n
        data = self.read_block8()
        marker, = op2.struct_q.unpack(data)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
        else:
            if self.is_debug_file:
                msg = 'get_marker - [8, %i, 8]; macro_rewind=%s\n' % (
                    marker, macro_rewind or rewind)
                self.binary_debug.write(msg)
        return marker

    #------------------------------------------------------------------
    @property
    def is_debug_file(self) -> bool:
        """interface to the op2 object"""
        return self.op2.is_debug_file

    @property
    def binary_debug(self):
        """interface to the op2 object"""
        return self.op2.binary_debug

    @property
    def log(self):
        """interface to the op2 object"""
        return self.op2.log

    @property
    def _endian(self) -> bytes:
        """interface to the op2 object"""
        return self.op2._endian

    @property
    def _uendian(self) -> str:
        """interface to the op2 object"""
        return self.op2._uendian

    @property
    def _encoding(self) -> str:
        """interface to the op2 object"""
        return self.op2._encoding

    @property
    def read_mode(self) -> int:
        """interface to the op2 object"""
        return self.op2.read_mode

    @property
    def debug_file(self):
        """interface to the op2 object"""
        return self.op2.debug_file

    #------------------------------------------------------------------
    # skip methods

#class OP2Skip:
    #def __init__(self, op2):
        #self.op2 = op2

    #@property
    #def is_debug_file(self):
        #return self.op2.is_debug_file

    #@property
    #def binary_debug(self):
        #return self.op2.binary_debug

    #@property
    #def log(self):
        #return self.op2.log

    def _skip_table_helper(self, warn=True):
        """
        Skips the majority of geometry/result tables as they follow a very standard format.
        Other tables don't follow this format.

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write(f'skipping table...{op2.table_name!r}\n')
        if warn:
            self.log.warning(f'    skipping table_helper = {op2.table_name}')

        self.read_markers([-1])
        unused_data = self._skip_record()
        self.read_3_markers([-2, 1, 0])
        unused_data = self._skip_record()
        self._skip_subtables()

    def _skip_subtables(self):
        """skips a set of subtables"""
        op2 = self.op2
        op2.isubtable = -3
        self.read_3_markers([-3, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            unused_data = self._skip_record()
            if self.is_debug_file:
                desc = self.desc_map.get(op2.table_name, '???')
                #assert desc != '???', self.table_name
                msgi = "skipping table_name = %r ({desc})".rstrip('(?)')
                self.log.debug(msgi)
            #if len(data) == 584:
                #self._parse_results_table3(data)
            #else:
                #data = self._parse_results_table4(data)

            op2.isubtable -= 1
            self.read_3_markers([op2.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])

    def _skip_record(self):
        """
        the skip version of ``_read_record``

        Returns
        -------
        record : None
            a record of None indicates a skipped block

        """
        unused_markers0 = self.get_nmarkers(1, rewind=False)
        record = self._skip_block()

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        while markers1[0] > 0:
            markers1 = self.get_nmarkers(1, rewind=False)
            record = self._skip_block()
            markers1 = self.get_nmarkers(1, rewind=True)
        return record

    def _skip_record_ndata(self, debug=True, macro_rewind=False):
        if self.size == 4:
            return self._skip_record_ndata4(debug=debug, macro_rewind=macro_rewind)
        return self._skip_record_ndata8(debug=debug, macro_rewind=macro_rewind)

    def _skip_record_ndata4(self, debug=True, macro_rewind=False):
        """the skip version of ``_read_record_ndata``"""
        op2 = self.op2
        marker0 = self.get_marker1_4(rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (
                marker0, macro_rewind))

        na = op2.n
        record, nrecord = self._skip_block_ndata()

        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - record = [%i, recordi, %i]; '
                                    'macro_rewind=%s\n' % (nrecord, nrecord, macro_rewind))

        if marker0*4 != nrecord:
            self.log.debug(str(marker0))
            op2.f.seek(na)
            op2.n = na
            if nrecord == 4:
                #self.read_3_markers([1, 0], macro_rewind=False)
                #self.show(500, types='ifs', endian=None, force=False)
                self.log.error(f'EmptyRecordError: marker0={marker0} nrecord={nrecord}')
                raise EmptyRecordError('nrecord=4')
            self.log.error(f'marker0={marker0} nrecord={nrecord}')
            #self.show(500, types='ifs', endian=None, force=False)
            #self.show(record, types='ifs', endian=None, force=False)
            self.log.error('returning to before data block is skipped...')
            self.show(500, types='ifs', endian=None, force=False)
            msg = f'marker0={marker0*4}*4 len(record)={nrecord}; table_name={op2.table_name!r}'
            raise FortranMarkerError(msg)

        #self.log.debug(f'marker0={marker0} nrecord={nrecord}')
        marker1 = self.get_marker1_4(rewind=True)

        if marker1 > 0:
            while marker1 > 0:
                marker1 = self.get_marker1_4(rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write(f'read_record - marker1 = [4, {marker1}, 4]\n')
                unused_recordi, nrecordi = self._skip_block_ndata()
                nrecord += nrecordi

                marker1 = self.get_marker1_4(rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write(f'read_record - marker1 = [4, {marker1}, 4]\n')
        return record, nrecord

    def _skip_record_ndata8(self, debug=True, macro_rewind=False):
        """the skip version of ``_read_record_ndata``"""
        op2 = self.op2
        marker0 = self.get_marker1_8(rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [8, %i, 8]; macro_rewind=%s\n' % (
                marker0, macro_rewind))
        record, nrecord = self._skip_block_ndata()

        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - record = [%i, recordi, %i]; '
                                    'macro_rewind=%s\n' % (nrecord, nrecord, macro_rewind))
        if marker0*8 != nrecord:
            msg = 'marker0=%s*8 len(record)=%s; table_name=%r' % (
                marker0*8, nrecord, op2.table_name)
            raise FortranMarkerError(msg)

        marker1 = self.get_marker1_8(rewind=True)

        if marker1 > 0:
            while marker1 > 0:
                marker1 = self.get_marker1_8(rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write(f'read_record - marker1 = [8, {marker1}, 8]\n')
                unused_recordi, nrecordi = self._skip_block_ndata()
                nrecord += nrecordi

                marker1 = self.get_marker1_8(rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write(f'read_record - marker1 = [8, {marker1}, 8]\n')
        return record, nrecord

    def _get_record_length(self):
        """
        The record length helps us figure out data block size, which is used
        to quickly size the arrays.  We just need a bit of meta data and can
        jump around quickly.

        Returns
        -------
        record_length : int
            the length of the data block

        """
        op2 = self.op2
        if self.is_debug_file:
            self.binary_debug.write('_get_record_length\n')
        len_record = 0
        n0 = op2.n
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('  markers0=%s\n' % markers0)

        n = op2.n
        unused_record = self._skip_block()
        len_record += op2.n - n - 8  # -8 is for the block
        if self.is_debug_file:
            self.binary_debug.write('  len_record=%s\n' % len_record)

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        while markers1[0] > 0:
            markers1 = self.get_nmarkers(1, rewind=False)
            if self.is_debug_file:
                self.binary_debug.write('  markers1=%s\n' % markers1)
            n = op2.n
            unused_record = self._skip_block()
            len_record += op2.n - n - 8  # -8 is for the block
            markers1 = self.get_nmarkers(1, rewind=True)
        self._goto(n0)
        return len_record

    def _skip_block(self):
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data :  since data can never be None, a None value
                indicates something bad happened.

        """
        return self._skip_block_ndata()[0]

    def _skip_block_ndata(self) -> tuple[None, int]:
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : None
            since data can never be None, a None value
            indicates something bad happened.
        ndata : int
            the length of the data block that was skipped

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)
        op2.n += 8 + ndata
        self._goto(op2.n)
        return None, ndata

    #---------------------------------------------------------------------------
    def _goto(self, n: int) -> None:
        """
        Jumps to position n in the file

        Parameters
        ----------
        n : int
            the position to goto

        """
        self.op2.n = n
        self.op2.f.seek(n)

    def is_valid_subcase(self) -> bool:
        """
        Lets the code check whether or not to read a subcase

        Returns
        -------
        is_valid : bool
            should this subcase defined by self.isubcase be read?

        """
        op2 = self.op2
        if not op2.is_all_subcases:
            if hasattr(op2, 'isubcase') and op2.isubcase in op2.valid_subcases:
                return True
            return False
        return True

    def read_results_table(self):
        """Reads a results table"""
        if self.size == 4:
            self.read_results_table4()
        else:
            self.read_results_table8()

    def read_results_table4(self) -> None:
        """Reads a results table"""
        op2 = self.op2
        if self.is_debug_file:
            self.binary_debug.write(f'read_results_table - {op2.table_name}\n')
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
            #self.binary_debug.write('marker = [4, -1, 4]\n')
        data = self._read_record_ndata4()[0]
        # (101, 0, 1, 0, 0, 0, 3)
        #self.show_data(data)  # TODO: what is this???

        self.read_3_markers([-2, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-2, 1, 0]---\n')
        data, ndata = self._read_record_ndata4()

        subtable_name = self.get_subtable_name4(op2, data, ndata)
        op2.subtable_name = subtable_name
        self._read_subtables()

    def read_results_table8(self) -> None:
        """Reads a results table"""
        op2 = self.op2
        if self.is_debug_file:
            self.binary_debug.write(f'read_results_table - {op2.table_name}\n')
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers8([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
            #self.binary_debug.write('marker = [4, -1, 4]\n')
        data = self._read_record_ndata8()[0]

        self.read_markers8([-2, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-2, 1, 0]---\n')
        data, ndata = self._read_record_ndata8()

        subtable_name = self.get_subtable_name8(op2, data, ndata)
        subtable_name = reshape_bytes_block(subtable_name)
        op2.subtable_name = subtable_name
        self._read_subtables()

    def get_subtable_name8(self, op2, data: bytes, ndata: int) -> bytes:
        if ndata == 16: # 8*2
            subtable_name = op2.struct_16s.unpack(data)
            subtable_name = reshape_bytes_block(subtable_name)
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 32: # 16*2
            #(name1, name2, 170, 170)
            subtable_name, = op2.struct_16s.unpack(data[:16])
            assert len(subtable_name) == 16, len(subtable_name)
            subtable_name = reshape_bytes_block(subtable_name)
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 56: # 28*2
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'16s5q', data)
            subtable_name = reshape_bytes_block(subtable_name)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        else:
            self.show_data(data, types='ifsqd', endian=None)
            raise NotImplementedError(len(data))
        return subtable_name

    def get_subtable_name4(self, op2, data: bytes, ndata: int) -> bytes:
        if ndata == 8:
            subtable_name = op2.struct_8s.unpack(data)
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 12:
            subtable_name, unused_ten = unpack(self._endian + b'8si', data)  # type: tuple[bytes, int]
            subtable_name = subtable_name.strip().decode(self._encoding)
            #assert ten == 10, self.show_data(data, types='ifs', endian=None)
            assert subtable_name in ['GPL', 'GPLS'], subtable_name
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        elif ndata == 612: # ???
            strings, ints, floats = self.show_data(data)
            msg = f'len(data) = {ndata:d}\n'
            #msg += 'strings  = %r\n' % strings
            #msg += 'ints     = %r\n' % str(ints)
            #msg += 'floats   = %r' % str(floats)
            print(msg)
            subtable_name, = op2.struct_8s.unpack(data[:8])
            print('subtable_name = %r' % subtable_name.strip())
        elif ndata == 24 and op2.table_name == b'ICASE':
            subtable_name = 'CASE CONTROL SECTION'
            #strings = 'CASE CONTROL SECTION\xff\xff\xff\xff'
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        if hasattr(self, 'subtable_name'):
            raise RuntimeError('the file has not been cleaned up; subtable_name_old=%s new=%s' % (
                op2.subtable_name, subtable_name))
        return subtable_name

    def generic_stop_table(self, data: bytes, ndata: int):  # pragma: no cover
        """print table data when things get weird"""
        strings, ints, floats = self.show_data(data)
        msg = 'Unhandled table length error\n'
        msg += f'table_name = {self.op2.table_name}\n'
        msg += f'len(data) = {ndata:d}\n'
        msg += 'strings  = %r\n' % strings
        msg += 'ints     = %r\n' % str(ints)
        msg += 'floats   = %r' % str(floats)
        raise NotImplementedError(msg)

    def read_geom_table(self):
        """Reads a geometry table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)

        if self.is_debug_file:
            self.binary_debug.write(f'read_geom_table - {op2.table_name}\n')
        self.read_markers([-1])
        data = self._read_record() # length=28=7*4

        self.read_3_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if self.size == 4:
            if ndata == 8:
                subtable_name, = self.op2.struct_8s.unpack(data)
            elif ndata == 28:
                # date = august 2, 2019
                #(b'OGPFB1  ', 8, 2, 19, 0, 1)
                #(b'OGPFB1  ', 8, 2, 19, 1, 1)
                fmt = self._endian + b'8s 3i 2i'
                subtable_name, month, day, year, zero, one = Struct(fmt).unpack(data)
                if zero != 0 or one != 1:  # pragma: no cover
                    self.log.warning(f'  subtable_name={subtable_name} possible error; zero={zero} one={one}')
                    #self.generic_stop_table(data, ndata)
                self._set_op2_date(month, day, year)
            else:
                self.generic_stop_table(data, ndata)
        else:
            if ndata == 16:
                subtable_name, = self.op2.struct_16s.unpack(data)
            else:
                self.generic_stop_table(data, ndata)
        op2.subtable_name = subtable_name.rstrip()
        self._read_subtables()

    def _read_subtables(self) -> None:
        """reads a series of subtables"""
        # this parameters is used for numpy streaming
        op2 = self.op2
        op2._table4_count = 0
        op2.is_table_1 = True
        op2._data_factor = 1

        #nstart = op2.n
        op2.isubtable = -3
        self.read_3_markers([-3, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write(f'***isubtable = {op2.isubtable:d}\n')
            self.binary_debug.write('---markers = [-3, 1, 0]---\n')

        # get the parsing functions (table3_parser, table4_parser)
        # or find out we're going to be skipping the tables
        #
        # table3 - the table with the meta data (e.g. subcase_id, time, is_stress/strain)
        # table4 - the actual results data
        #
        # we indicate table3/4 by isubtable, which starts from -3 (so table3) and counts
        # down (yes down) to 4 to indicate table4.  If we count down again, we end up
        # back at table 3 (with isubtable=-5), which will occur in the case of multiple
        # times/element types/results in a single macro table (e.g. OUG, OES).
        #table_mapper = op2._get_table_mapper()
        table_mapper = op2.table_mapper
        table_name = op2.table_name
        desc = '???'
        if table_name in table_mapper:
            #if self.read_mode == 2:
                #self.log.debug("table_name = %r" % table_name)
            try:
                table3_parser, table4_parser = table_mapper[table_name]
            except:
                table3_parser, table4_parser, desc = table_mapper[table_name]
            passer = False
        else:
            if table_name in op2.op2_reader.desc_map:
                desc = op2.op2_reader.desc_map[table_name]

            if self.read_mode == 2:
                self.log.info(f'skipping table_name = {table_name!r} ({desc})')
                #assert desc != '???', table_name
                    #raise NotImplementedError(table_name)
            table3_parser = None
            table4_parser = None
            passer = True

        # we need to check the marker, so we read it and rewind, so we don't
        # screw up our positioning in the file
        markers = self.get_nmarkers(1, rewind=True)
        #if markers[0] == 0:
            #self.log.debug('    returning early')
            #return

        if self.is_debug_file:
            self.binary_debug.write(f'---marker0 = {markers}---\n')

        # while the subtables aren't done
        while markers[0] != 0:
            #print(markers)
            op2.is_start_of_subtable = True
            if self.is_debug_file:
                self.binary_debug.write(f'***isubtable = {op2.isubtable:d}\n')

            try:
                self._read_subtable_3_4(table3_parser, table4_parser, passer)
            except EmptyRecordError:
                self.log.error('catching EmptyRecordError')
                self.read_markers([1, 0], macro_rewind=False)
                #n = op2.n
                #try:
                marker146 = self.get_marker1(rewind=True)
                #except AssertionError:
                    #self.log.debug('resetting n!')
                    #op2.f.seek(n)
                    #op2.n = n
                    #raise
                op2.isubtable -= 1
                if marker146 == 146:
                    continue
                break
            except Exception:  # pragma: no cover
                print(f'failed reading {table_name} isubtable={op2.isubtable:d}')
                raise
            #force_table4 = self._read_subtable_3_4(table3_parser, table4_parser, passer)
            op2.isubtable -= 1

            iloc = op2.f.tell()
            if op2._nastran_format == 'optistruct' and 0:
                #  (4, -4, 4, 4, 1, 4, 4, 0, 4)
                nbytes = self.size * 9
                datai = op2.f.read(nbytes)  # 48=9*4 bytes
                structi = Struct(self._endian + mapfmt(b'9i', self.size))
                outi = structi.unpack(datai)
                op2.n += nbytes
                if not outi[1] == op2.isubtable:
                    self.log.warning(f'outi={outi} isubtable={op2.isubtable:d}')
            else:
                #self.show_ndata(36, types='i', force=False, endian=None)
                try:
                    self.read_3_markers([op2.isubtable, 1, 0])
                    #self.log.debug('markers=%s' % [op2.isubtable, 1, 0])
                    df = op2.f.tell() - iloc
                    if op2.size == 8:
                        assert df == 48, df
                    else:
                        assert df == 36, df
                except FortranMarkerError:
                    self.log.error(f'isubtable={op2.isubtable:d}')
                    op2.f.seek(iloc)
                    op2.n = iloc

                    self.show(4*3*3, types='i')
                    self.show(500)

                    marker0 = self.get_nmarkers(1, rewind=True)[0]
                    self.log.debug(f'marker0 = {marker0}')
                    if marker0 < op2.isubtable:
                        raise RuntimeError('marker0 < isubtable; marker0=%s isubtable=%s' % (
                            marker0, op2.isubtable))
                        #self.read_3_markers([marker0, 1, 0])
                        ##self.log.debug('markers=%s' % [marker0, 1, 0])
                        #self.show(200)
                        #break
                    raise
            markers = self.get_nmarkers(1, rewind=True)

        if self.is_debug_file:
            self.binary_debug.write(f'breaking on marker={markers}\n')

        # we've finished reading all subtables, but have one last marker to read
        marker = self.get_marker1(rewind=False, macro_rewind=False)
        assert marker == 0, marker
        op2._finish()

    def _read_subtable_3_4(self,
                           table3_parser: Optional[Callable],
                           table4_parser: Optional[Callable],
                           passer: Optional[Callable]) -> Optional[bool]:
        """
        Reads a series of subtable 3/4

        Parameters
        ----------
        table3_parser : bool / function
            None : just to break the code if we're on the array sizing step
            function : the table 3 reading function
        table4_parser :bool / function
            None : just to break the code if we're on the array sizing step
            function : the table 4 reading function
        passer : bool
            flag to see if we're skipping tables

        Returns
        -------
        flag : bool
            True : ???
            False : failed???
            None : passed???

        """
        op2 = self.op2
        IS_TESTING = op2.IS_TESTING
        if self.binary_debug:
            self.binary_debug.write('-' * 60 + '\n')
        # this is the length of the current record inside table3/table4
        record_len = self._get_record_length()
        if self.is_debug_file:
            self.binary_debug.write(f'record_length = {record_len:d}\n')

        oes_nl = [b'OESNLXD', b'OESNL1X', b'OESNLXR'] # 'OESCP'?
        factor = self.factor
        #print('record_len =', record_len)
        table_name = op2.table_name
        if record_len == 584 * factor:  # table3 has a length of 584
            if table_name in oes_nl and hasattr(op2, 'num_wide') and op2.num_wide == 146:
                data_code_old = deepcopy(op2.data_code)

            if self.load_as_h5:
                assert self.h5_file is not None, self.h5_file
            op2.data_code = {
                '_encoding' : self._encoding,
                'load_as_h5' : self.load_as_h5,
                'h5_file' : self.h5_file,
                'size' : self.size,
            }
            op2.obj = None
            data, ndata = self._read_record_ndata()
            if not passer:
                try:
                    table3_parser(data, ndata)
                except SortCodeError:
                    if self.is_debug_file:
                        self.binary_debug.write('except SortCodeError!\n')
                    if table_name in oes_nl:
                        update_op2_datacode(op2, data_code_old)

                        n = table4_parser(data, ndata)
                        #print(data_code_old)
                        if not isinstance(n, integer_types):
                            msg = 'n is not an integer; table_name=%s n=%s table4_parser=%s' % (
                                self.op2.table_name, n, table4_parser)
                            raise TypeError(msg)
                        if IS_TESTING:
                            self._run_checks(table4_parser)

                        if self.read_mode == 1:
                            #op2_reader._goto(n)
                            #n = op2_reader._skip_record()
                            #if hasattr(self.op2, 'table_name'):
                                #print('***_init_vector_counter', self.op2.table_name)
                            #print('record_len', record_len)
                            self.op2._init_vector_counter(record_len)
                        else:
                            self.op2._reset_vector_counter()

                        #print('except...')
                        return False
                    raise RuntimeError(op2.code_information())
                #if hasattr(op2, 'isubcase'):
                    #print("code = ", op2._get_code())
        else:
            if table_name in GEOM_TABLES:
                if passer:
                    data = self._skip_record()
                else:
                    data, ndata = self._read_record_ndata()
                    unused_n = table4_parser(data, ndata)

            elif passer or not self.is_valid_subcase():
                data = self._skip_record()
            else:
                if hasattr(op2, 'num_wide'):
                    # num_wide is the result size and is usually found in
                    # table3, but some B-list tables don't have it
                    unused_n = op2._read_subtable_results(table4_parser, record_len)
                else:
                    data, ndata = self._read_record_ndata()
                    unused_n = table4_parser(data, ndata)
                    if IS_TESTING:
                        self._run_checks(table4_parser)
                #del n
        return None

    def _run_checks(self, table4_parser):
        """helper method"""
        if table4_parser != self.op2._table_passer:
            str(self.op2.code_information())
        if hasattr(self, 'obj'):
            str(self.obj.get_stats())

    def show(self, n: int, types: str='ifs', endian: Optional[str]=None,
             force: bool=False):  # pragma: no cover
        """
        shows the next N bytes

        Parameters
        ----------
        n : int
            the number of bytes to show
        types : str; default='ifs'
            the data types to show
        endian : str; default=None -> active endian
            the data endian
        force : bool; default=False
            overwrite the n=2000 limiter

        """
        op2 = self.op2
        assert op2.n == op2.f.tell()
        data = op2.f.read(n)
        strings, ints, floats = self.show_data(data, types=types, endian=endian, force=force)
        op2.f.seek(op2.n)
        return strings, ints, floats

    def show_data(self, data: bytes, types: str='ifs', endian: Optional[str]=None,
                  force: bool=False):  # pragma: no cover
        """
        Shows a data block as various types

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}
        force : bool; default=False
            overwrite the n=2000 limiter

        .. warning:: 's' is apparently not Python 3 friendly

        """
        #ifsdqlILQ
        return self._write_data(sys.stdout, data, types=types, endian=endian, force=force)

    def _write_data(self, f, data: bytes, types: str='ifs',
                    endian: Optional[str]=None, force: bool=False):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}
        force : bool; default=False
            overwrite the n=2000 limiter

        """
        n = len(data)
        if not force and n > 2000:
            self.log.warning(f'limiting n={n} to 2000')
            n = 2000
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._uendian
            assert endian is not None, endian

        f.write(f'\nndata = {n}:\n')
        for typei in types:
            assert typei in 'sifdq lILQ', f'type={typei!r} is invalid; use sifdq lILQ'

        data4 = data[:nints * 4]
        #data8 = data[:ndoubles * 8]
        if 's' in types:
            strings = unpack('%s%is' % (endian, n), data[:n])
            f.write(f"  strings = {strings}\n")
        if 'i' in types:
            ints = unpack('%s%ii' % (endian, nints), data4)
            f.write(f'  ints    = {ints}\n')
        if 'f' in types:
            floats = unpack('%s%if' % (endian, nints), data4)
            f.write(f'  floats  = {floats}\n')
        if 'd' in types:
            doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
            f.write(f'  doubles (float64) = {doubles}\n')

        if 'l' in types:
            longs = unpack('%s%il' % (endian, nints), data4)
            f.write(f'  long  = {longs}\n')
        if 'I' in types:
            ints2 = unpack('%s%iI' % (endian, nints), data4)
            f.write(f'  unsigned int = %s\n' % str(ints2))
        if 'L' in types:
            longs2 = unpack('%s%iL' % (endian, nints), data4)
            f.write(f'  unsigned long = {longs2}\n')
        if 'q' in types:
            longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write(f'  long long (int64) = {longs}\n')
        if 'Q' in types:
            longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write(f'  unsigned long long (int64) = {longs}\n')
        f.write('\n')
        return strings, ints, floats

    def show_ndata(self, n: int, types: str='ifs', force: bool=False, endian=None):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types, force=force, endian=endian)

    def _write_ndata(self, f, n: int, types: str='ifs', force: bool=False, endian=None):  # pragma: no cover
        """Useful function for seeing what's going on locally when debugging."""
        op2 = self.op2
        nold = op2.n
        data = op2.f.read(n)
        op2.n = nold
        op2.f.seek(op2.n)
        return self._write_data(f, data, types=types, force=force, endian=endian)

def eqexin_to_nid_dof_doftype(eqexin1, eqexin2) -> tuple[Any, Any, Any]:
    """assemble dof table"""
    dof = eqexin2[1::2] // 10
    doftype = eqexin2[1::2] - 10 * dof
    nid = eqexin2[::2]

    # eqexin is in external sort, so sort it
    i = eqexin1[1::2].argsort()
    dof = dof[i]
    doftype = doftype[i]
    nid = nid[i]
    return nid, dof, doftype

def get_table_size_from_ncolumns(table_name: str, nvalues: int, ncolumns: int) -> int:
    nrows = nvalues // ncolumns
    if nvalues % ncolumns != 0:
        msg = 'table=%s: nrows=nvalues/ncolumns=%s/%s=%s; nrows=%s must be an int' % (
            table_name, nvalues, ncolumns, nrows, nvalues / ncolumns)
        raise RuntimeError(msg)
    return nrows

def dscmcol_dresp1(responses: dict[int, dict[str, Any]],
                   nresponses_dresp1: int,
                   ints, floats) -> None:
    """helper for DSCMCOL"""
    idata = 0
    if nresponses_dresp1 == 0:
        return

    for iresp in range(nresponses_dresp1):
        # Record – Type 1 Responses
        #   1 IRID  I Internal response identification number
        #   2 RID   I External response identification number
        #   3 RTYPE I Response Type
        # RTYPE=8 Static force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 VIEWID I View element identification number
        #   9 SEID I Superelement identification number
        # RTYPE=10 Composite stress
        #   4 EID I Element identification number
        #   5 COMP I Stress component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 PLY I Ply number
        #   9 SEID I Superelement identification number
        # RTYPE=13 Static SPC force
        #   4 GRID I Grid identification number
        #   5 COMP I SPC force component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=14 Element static strain energy
        #   4 EID I Element identification number
        #   5 COMP I Strain energy component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=17 Compliance
        #   4 UNDEF(2) None
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF(3) None
        # RTYPE=21 Frequency response velocity
        #   4 GRID I Grid identification number
        #   5 COMP I Velocity component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=22 Frequency response acceleration
        #   4 GRID I Grid identification number
        #   5 COMP I Acceleration component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=23 Frequency response SPC Force
        #   4 GRID I Grid identification number
        #   5 COMP I SPC Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=24 Frequency response stress
        #   4 EID I Element identification number
        #   5 COMP I Stress component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=25 Frequency response force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 FREQ RS Frequency
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=26 RMS displacement
        #   4 GRID I Grid identification number
        #   5 COMP I RMS displacement component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=27 RMS velocity
        #   4 GRID I Grid identification number
        #   5 COMP I RMS velocity component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=28 RMS acceleration
        #   4 GRID I Grid identification number
        #   5 COMP I RMS acceleration component number
        #   6 SUBCASE I Subcase identification number
        #   7 UNDEF None
        #   8 RANDPS I RANDPS ID
        #   9 SEID I Superelement identification number
        # RTYPE=30 PSD velocity
        # 4 GRID I Grid identification number
        # 5 COMP I PSD velocity component number
        # 6 SUBCASE I Subcase identification number
        # 7 FREQ RS Frequency
        # 8 RANDPS I RANDPS ID
        # 9 SEID I Superelement identification number
        # RTYPE=60 Transient response displacement
        # 4 GRID I Grid identification number
        # 5 COMP I Displacement component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=61 Transient response velocity
        # 4 GRID I Grid identification number
        # 5 COMP I Velocity component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=62 Transient response acceleration
        # 4 GRID I Grid identification number
        # 5 COMP I Acceleration component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=63 Transient response SPC Force
        # 4 GRID I Grid identification number
        # 5 COMP I SPC force component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=64 Transient response stress
        # 4 EID I Element identification number
        # 5 COMP I Stress component number
        # 6 SUBCASE I Subcase identification number
        # 7 TIME RS Time
        # 8 UNDEF None
        # 9 SEID I Superelement identification number
        # RTYPE=65 Transient response force
        #   4 EID I Element identification number
        #   5 COMP I Force component number
        #   6 SUBCASE I Subcase identification number
        #   7 TIME RS Time
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=81 Aeroelastic divergence
        #   4 SUBCASE I Subcase identification number
        #   5 UNDEF None
        #   6 ROOT I Root
        #   7 MACH RS Mach number
        #   8 UNDEF None
        #   9 SEID I Superelement identification number
        # RTYPE=82 Aeroelastic trim
        #   4 SUBCASE I Subcase identification number
        #   5 UNDEF None
        #   6 XID I XID
        #   7 UNDEF(2) None
        #   9 SEID I Superelement identification number
        # RTYPE=83 Aeroelastic stability derivative
        #   4 SUBCASE I Subcase identification number
        #   5 RU I R/U
        #   6 COMP I Component number
        #   7 UNDEF None
        #   8 XID I XID
        #   9 SEID I Superelement identification number
        # RTYPE=84 Aeroelastic flutter damping
        #   4 SUBCASE I Subcase identification number
        #   5 MODE I Mode number
        #   6 DENSITY RS Density
        #   7 MACH RS Mach number
        #   8 VEL RS Velocity
        #   9 SEID I Superelement identification number
        # End RTYPE

        internal_response_id = ints[idata]
        external_response_id = ints[idata+1]
        response_type = ints[idata+2]
        #print(f'internal_response_id={internal_response_id} '
              #f'external_response_id={external_response_id} response_type={response_type}')

        if response_type == 1:
            # RTYPE=1 Weight
            #   4 UNDEF(5) None
            #   9 SEID I Superelement identification number
            seid = ints[idata+8]
            response = {'name': 'weight', 'seid': seid}
            #print(f'  seid={seid} (weight)')
        elif response_type == 2:
            # RTYPE=2 Volume
            #   4 UNDEF(5) None
            #   9 SEID I Superelement identification number
            seid = ints[idata+8]
            response = {'name': 'volume', 'seid': seid}
            #print(f'  seid={seid} (volume)')
        elif response_type == 3:
            # RTYPE=3 Buckling
            #   4 MODE    I Mode number
            #   5 UNDEF None
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            mode = ints[idata+3]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'buckling', 'mode_num': mode, 'subcase': subcase, 'seid': seid}
            #print(f'  mode={mode} subcase={subcase} seid={seid} (buckling)')

        elif response_type == 4:
            # RTYPE=4 Normal modes
            #   4 MODE    I Mode number
            #   5 UNDEF None
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            mode_num = ints[idata+3]
            # blank
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'normal modes', 'mode_num': mode_num, 'subcase': subcase, 'seid': seid}
            #print(f'  mode_num={mode_num} subcase={subcase} seid={seid} (normal modes)')
        elif response_type == 5:
            # RTYPE=5 Static displacement
            #   4 GRID    I Grid identification number
            #   5 COMP    I Displacement component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF(2) None
            #   9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'static displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} seid={seid} (static displacement)')
        elif response_type == 6:
            # RTYPE=6 Static stress
            # 4 EID     I Element identification number
            # 5 COMP    I Stress component number
            # 6 SUBCASE I Subcase identification number
            # 7 UNDEF(2) None
            # 9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            response = {'name': 'static stress', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} seid={seid} (static stress)')
        elif response_type == 7:
            # RTYPE=7 Static strain
            # 4 EID     I Element identification number
            # 5 COMP    I Strain component number
            # 6 SUBCASE I Subcase identification number
            # 7 UNDEF None
            # 8 VIEWID  I View element identification number
            # 9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            view_id = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'static strain', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'view id': view_id, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} view_id={view_id} seid={seid} (static strain)')
        elif response_type == 9:
            # RTYPE=9 Composite failure
            #   4 EID     I Element identification number
            #   5 COMP    I Failure component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF None
            #   8 PLY     I Ply number
            #   9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'composite failure', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} ply={ply} seid={seid} (composite failure)')
        elif response_type == 11:
            # RTYPE=11 Composite strain
            #   4 EID     I Element identification number
            #   5 COMP    I Strain component number
            #   6 SUBCASE I Subcase identification number
            #   7 UNDEF None
            #   8 PLY     I Ply number
            #   9 SEID    I Superelement identification number
            eid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'composite strain', 'eid': eid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  eid={eid} comp={comp} subcase={subcase} ply={ply} seid={seid} (composite strain)')
        elif response_type == 15:
            # CEIG - complex eigenvalues
            mode_num = ints[idata+3]
            subcase = ints[idata+5]
            seid = ints[idata+8]
            #print(ints[idata+4:idata+9])
            #print(floats[idata+4:idata+9])
            #print(f'internal_response_id={internal_response_id} '
                  #f'external_response_id={external_response_id} response_type={response_type}')
            #print(f'  mode_num={mode_num} subcase={subcase} seid={seid} (CEIG)')
            response = {'name': 'ceig', 'mode_num': mode_num, 'subcase': subcase, 'seid': seid}

        elif response_type == 19:
            # RTYPE=19 Equivalent radiated power
            # 4 PANEL      I Element identification number
            # 5 FLAG       I A subcase ID based code. +: magnitude; –:density
            # 6 SUBCASE    I Subcase identification number
            # 7 FREQUENCY RS Frequency
            # 8 UNDEF None
            # 9 SEID       I Superelement identification number
            panel = ints[idata+3]
            flag = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            seid = ints[idata+8]
            response = {'name': 'equivalent radiated power', 'panel': panel, 'flag': flag,
                        'subcase': subcase, 'freq': freq, 'seid': seid}
            #print(f'  panel={panel} flag={flag} subcase={subcase} freq={freq} seid={seid} '
                  #'(equivalent radiated power)')
        elif response_type == 20:
            # RTYPE=20 Frequency response displacement
            # 4 GRID    I Grid identification number
            # 5 COMP    I Displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 UNDEF None
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            ply = ints[idata+7]
            unused_freq8 = floats[idata+7] # also freq despite DMAP saying this is unused...
            seid = ints[idata+8]
            response = {'name': 'frequency response displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'ply={ply} seid={seid} (frequency response displacement)')
        elif response_type == 21:
            # RTYPE=21 frequency response stress???
            # 4 GRID    I Grid identification number
            # 5 COMP    I Displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 UNDEF None
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            #sixI = ints[idata+6]
            ply = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'frequency response stress?', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'ply': ply, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'ply={ply} seid={seid} (frequency response stress?)')
        elif response_type == 29:
            # RTYPE=29 PSD displacement
            # 4 GRID    I Grid identification number
            # 5 COMP    I PSD displacement component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ    RS Frequency
            # 8 RANDPS  I RANDPS ID
            # 9 SEID I   Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            randps = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'psd displacement', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'freq': freq, 'randps': randps, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'randps={randps} seid={seid} (PSD displacement)')
        elif response_type == 31:
            # RTYPE=31 PSD acceleration
            # 4 GRID    I Grid identification number
            # 5 COMP    I PSD acceleration component number
            # 6 SUBCASE I Subcase identification number
            # 7 FREQ   RS Frequency
            # 8 RANDPS  I RANDPS ID
            # 9 SEID    I Superelement identification number
            grid = ints[idata+3]
            comp = ints[idata+4]
            subcase = ints[idata+5]
            freq = floats[idata+6]
            randps = ints[idata+7]
            seid = ints[idata+8]
            response = {'name': 'psd acceleration', 'grid': grid, 'component': comp,
                        'subcase': subcase, 'freq': freq, 'randps': randps, 'seid': seid}
            #print(f'  grid={grid} comp={comp} subcase={subcase} freq={freq} '
                  #f'randps={randps} seid={seid} (PSD acceleration)')
        else:  # pragma: no cover
            print(f'internal_response_id={internal_response_id} '
                  f'external_response_id={external_response_id} response_type={response_type}')
            raise NotImplementedError(response_type)
        response['internal_response_id'] = internal_response_id
        response['external_response_id'] = external_response_id
        response['response_type'] = response_type
        response['iresponse'] = iresp
        response['response_number'] = 1
        responses[iresp] = response
        idata += 9
    return

def dscmcol_dresp2(responses: dict[int, dict[str, Any]],
                   nresponses_dresp2: int,
                   ints, floats) -> None:
    """helper for DSCMCOL"""
    if nresponses_dresp2 == 0:
        return
    nresponses = len(responses)

    idata = 0
    for iresp in range(nresponses_dresp2):
        # Word Name Type Description
        # 1 IRID      I Internal response identification number
        # 2 RID       I External response identification number
        # 3 SUBCASE   I Subcase identification number
        # 4 DFLAG     I Dynamic response flag (See Note)
        # 5 FREQTIME RS Frequency or time step
        # 6 SEID      I Superelement identification number
        internal_response_id = ints[idata]
        external_response_id = ints[idata+1]
        subcase = ints[idata+2]
        dflag = ints[idata+3]
        freqtime = floats[idata+4]
        seid = ints[idata+5]
        iresp2 = iresp + nresponses
        response = {
            'iresponse': iresp2,
            'response_number': 2,
            'internal_response_id': internal_response_id,
            'external_response_id': external_response_id,
            'subcase': subcase, 'dflag': dflag,
            'freq': freqtime, 'seid': seid}
        responses[iresp2] = response
        #print(f'internal_response_id={internal_response_id} '
              #f'external_response_id={external_response_id} '
              #f'subcase={subcase} dflag={dflag} freq/time={freqtime} seid={seid}')
        idata += 6
    return

def update_op2_datacode(op2, data_code_old):
    op2.data_code = data_code_old
    for key, value in data_code_old.items():
        if key == 'size':
            continue
        setattr(op2, key, value)

def _cast_matrix_mat(GCi: np.ndarray, GCj: np.ndarray,
                     mrows: int, ncols: int,
                     reals: np.ndarray,
                     tout: int,
                     dtype: str, log: SimpleLogger):
    """helper method for _read_matrix_mat"""
    #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
    #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)

    # we subtract 1 to the indicides to account for Fortran
    GCi = np.array(GCi, dtype='int32') - 1
    GCj = np.array(GCj, dtype='int32') - 1
    try:
        if dtype == '???':
            matrix = None
            log.warning('what is the dtype?')
        elif tout in {1, 2}:
            # real
            real_array = np.array(reals, dtype=dtype)
            matrix = scipy.sparse.coo_matrix(
                (real_array, (GCi, GCj)),
                shape=(mrows, ncols), dtype=dtype)
            #log.info(f'created {self.table_name} (real)')
        elif tout in {3, 4}:
            # complex
            real_array = np.array(reals, dtype=dtype)
            nvalues_matrix = real_array.shape[0] // 2
            real_complex = real_array.reshape((nvalues_matrix, 2))
            real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
            #if self.binary_debug:
                #self.binary_debug.write('reals = %s' % real_complex[:, 0])
                #self.binary_debug.write('imags = %s' % real_complex[:, 1])
                #self.binary_debug.write('real_imag = %s' % real_imag)
            matrix = scipy.sparse.coo_matrix(
                (real_imag, (GCi, GCj)),
                shape=(mrows, ncols), dtype=dtype)
            #msg = 'created %s (complex)' % self.table_name
            #log.debug(msg)
            #raise RuntimeError(msg)
        else:
            raise RuntimeError('this should never happen')
    except ValueError:
        log.warning('shape=(%s, %s)' % (mrows, ncols))
        log.warning('cant make a coo/sparse matrix...trying dense')

        if dtype == '???':
            matrix = None
            log.warning('what is the dtype?')
        else:
            real_array = np.array(reals, dtype=dtype)
            log.debug('shape=%s mrows=%s ncols=%s' % (
                str(real_array.shape), mrows, ncols))
            if len(reals) == mrows * ncols:
                real_array = real_array.reshape(mrows, ncols)
                log.info(f'created {op2.table_name}')
            else:
                log.warning(f'cant reshape because invalid sizes : created {op2.table_name}')

            matrix = real_array
    return matrix



def read_dofs(op2: OP2, size: int=4) -> None:
    op2.log.debug('read_dofs')
    dofs = []
    while 1:
        #self.show_ndata(32, types='if')
        tell = op2.f.tell()
        data = op2.f.read(12)
        out = unpack(b'<3i', data)
        if out[1] != 2:
            #print(out)
            op2.f.seek(tell)
            break
        op2.n += 16
        #if
        #tell = op2.f.tell()
        data = op2.f.read(20)
        a_int, b_int, c_int, d_int, e_float = unpack(op2._endian + b'3i f i', data)
        dofs.append((a_int, b_int, c_int, d_int, e_float))
        #out = op2.struct_3i.unpack(data)
        op2.n += 16
        #op2.show_data(data, types='if')
        #op2.show_ndata(64, types='if')
        #bbb
    for dof in dofs:
        print('dof', dof)
    print()
    return dofs

def reshape_trmbd(element_name: str, nnodes: int, int_data, float_data):
    ndata_per_element = 1 + nnodes + 3 * nnodes  # 1+4*(nnnodes) = 1+4*2 = 9
    n_elements = int_data.shape[0] // ndata_per_element  # elid + 2 grid + 3*2 euler angles
    assert int_data.shape[0] % ndata_per_element == 0, f'{element_name}: nnodes={nnodes} int_data.shape={int_data.shape}'

    int_data = int_data.reshape(n_elements, ndata_per_element)
    float_data = float_data.reshape(n_elements, ndata_per_element)
    return n_elements, int_data, float_data
