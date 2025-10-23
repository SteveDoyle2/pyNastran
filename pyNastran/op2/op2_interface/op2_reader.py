"""
Defines various tables that don't fit in other sections:
  - OP2Reader
    - read_aemonpt(self)
    - read_monitor(self)

    - read_cstm(self)
    - read_dit(self)
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
  - functions
    - _get_matrix_row_fmt_nterms_nfloats(nvalues, tout, endian)

"""
from __future__ import annotations
import os
import sys
from copy import deepcopy
from itertools import count
from functools import partial
from struct import unpack, Struct  # , error as struct_error
from typing import Optional, Callable, TYPE_CHECKING

import numpy as np

from cpylog import SimpleLogger
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.f06.errors import FatalError
from pyNastran.f06.flutter_response import FlutterResponse
from pyNastran.f06.f06_tables.trim import (
    MonitorLoads, TrimResults, ControllerState,
    AeroPressure, AeroForce,
    TrimVariables, TrimDerivatives,
    HingeMomentDerivatives,
    ControlSurfacePositionHingeMoment,)

from pyNastran.op2.errors import FortranMarkerError, SortCodeError, EmptyRecordError
from pyNastran.op2.result_objects.eqexin import EQEXIN
#from pyNastran.op2.result_objects.matrix import Matrix
from pyNastran.op2.result_objects.matrix_dict import MatrixDict
from pyNastran.op2.result_objects.qualinfo import QualInfo
from pyNastran.op2.result_objects.op2_results import CSTM
from pyNastran.op2.result_objects.campbell import CampbellData
from pyNastran.op2.op2_interface.msc_tables import MSC_GEOM_TABLES
from pyNastran.op2.op2_interface.nx_tables import NX_GEOM_TABLES

from pyNastran.op2.op2_interface.utils import (
    mapfmt, reshape_bytes_block,
    reshape_bytes_block_size)

from .version import parse_nastran_version
from .gpdt import (read_gpdt, read_bgpdt,
                   get_table_size_from_ncolumns)
from pyNastran.op2.op2_interface.read_extdb import (
    read_extdb, read_tug1, read_mef1)
from pyNastran.op2.op2_interface.read_external_superelement import read_cmodext
from pyNastran.op2.op2_interface.read_matrix_matpool import read_matrix_matpool
from pyNastran.op2.result_objects.monpnt import MONPNT1, MONPNT3

#from pyNastran.op2.op2_interface.read_matrix import (
    #read_matrix_mat)
from pyNastran.op2.op2_interface.read_optimization import (
    read_dbcopt, read_descyc, read_destab, read_dscmcol,
    read_hisadd, read_r1tabrg)
from pyNastran.op2.op2_interface.read_trmbu_trmbd import read_trmbu, read_trmbd

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

#class MinorTables:
    #def __init__(self, op2_reader):
        #self.op2_reader = op2_reader


class SubTableReadError(Exception):
    pass


GEOM_TABLES = MSC_GEOM_TABLES + NX_GEOM_TABLES


# https://pyyeti.readthedocs.io/en/latest/modules/nastran/generated/pyyeti.nastran.bulk.wtextseout.html
EXTSEOUT = [
    'KAA', 'MAA', 'BAA', 'K4XX', 'PA', 'GPXX', 'GDXX', 'RVAX',
    'VA', 'MUG1', 'MUG1O', 'MES1', 'MES1O', 'MEE1', 'MEE1O', 'MGPF',
    'MGPFO', 'MEF1O', 'MQG1', 'MQG1O', 'MQMG1', 'MQMG1O',
    'MEF1',
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

        self.op2: OP2 = op2

        fread_gpdt = partial(read_gpdt, self)
        fread_bgpdt = partial(read_bgpdt, self)
        fread_extdb = partial(read_extdb, self)
        fread_tug1 = partial(read_tug1, self)
        fread_mef1 = partial(read_mef1, self)
        fread_matrix_matpool = partial(read_matrix_matpool, self)
        self.mapped_tables = {
            b'RST': (self.read_rst, 'restart file?'),
            b'GPL': (self.read_gpl, 'grid point list'),
            b'GPLS': (self.read_gpls, 'grid point list (superelement)'),

            # GPDT  - Grid point definition table
            b'GPDT': (fread_gpdt, 'grid point locations'),
            b'GPDTS': (fread_gpdt, 'grid point locations (superelement)'),

            # BGPDT - Basic grid point definition table.
            b'BGPDT': (fread_bgpdt, 'grid points in cid=0 frame'),
            b'BGPDTS': (fread_bgpdt, 'grid points in cid=0 (superelement)'),
            b'BGPDTOLD': (fread_bgpdt, 'grid points in cid=0 frame'),
            b'BGPDTVU': (fread_bgpdt, 'VU grid points in cid=0 frame'),

            # optimization
            b'DESCYC': (partial(read_descyc, self), 'design iteration'),
            b'DBCOPT': (partial(read_dbcopt, self), 'design variable history table'),
            b'DSCMCOL': (partial(read_dscmcol, self), 'creates op2_results.responses.dscmcol'),
            b'DESTAB':  (partial(read_destab, self), 'creates op2_results.responses.desvars'),

            #b'MEFF': self.read_meff,

            # flutter
            b'OVG': (partial(read_ovg, self), 'aeroelastic velocity; creates op2_results.vg_vf_responses'),

            # trim
            b'OAEROTV': (partial(read_oaerotv, self), 'trim variables'),
            b'OAEROSCD': (partial(read_oaeroscd, self), 'stability and control derivatives'),
            b'OAERCSHM': (partial(read_oaercshm, self), 'control surface position & hinge moment'),
            b'OAEROHMD': (partial(read_oaerohmd, self), 'hinge moment derivatives'),
            b'OAEROP': (partial(read_oaerop, self), 'aero pressures'),
            b'OAEROF': (partial(read_oaerof, self), 'aero forces'),

            b'INTMOD': (self.read_intmod, '???'),
            b'HISADD': (partial(read_hisadd, self), 'optimization history; op2_results.responses.convergence_data'),
            #b'MEF1': (self.read_extdb, 'external superlelement matrix'),
            b'EXTDB': (fread_extdb, 'external superlelements'),
            b'OMM2': (self.read_omm2, 'max/min table'),
            b'STDISP': (self.read_stdisp, 'aero-structural displacement?'),
            b'TOL': (self.read_tol, 'time output list?'),
            b'PCOMPT': (self._read_pcompts, 'NX: LAM option input from the PCOMP bulk entry'),
            b'PCOMPTS': (self._read_pcompts, 'NX: LAM option input from the PCOMP bulk entry (superelement)'),
            b'MONITOR': (self.read_monitor, 'MONITOR point output'),
            b'AEMONPT': (self.read_aemonpt, 'aero matrix'),
            b'FOL': (self.read_fol, 'frequency output list'),
            b'FRL': (self.read_frl, 'frequency response list'),
            b'SDF': (self.read_sdf, 'aero-structural displacement?'),
            b'IBULK': (self.read_ibulk, 'explicit bulk data'),
            b'ICASE': (self.read_icase, 'explicit case control'),
            b'CASECC': (self.read_casecc, 'case control'),
            b'XCASECC': (self.read_xcasecc, 'case control'),

            b'CDDATA': (self.read_cddata, 'Cambell diagram'),
            b'CMODEXT': (partial(read_cmodext, self), 'Component mode synthesis - external'),

            #MSC
            #msc / units_mass_spring_damper
            b'UNITS': (self._read_units, 'units'),
            #b'CPHSF': self._read_cphsf,

            # element matrices
            #b'KELM': self._read_element_matrix,
            #b'MELM': self._read_element_matrix,
            #b'BELM': self._read_element_matrix,
            #b'KELMP': self._read_element_matrix,
            #b'MELMP': self._read_element_matrix,

            # element dictionaries
            b'KDICT': (self._read_dict, 'matrix'),
            b'MDICT': (self._read_dict, 'matrix'),
            b'BDICT': (self._read_dict, 'matrix'),
            b'KDICTP': (self._read_dict, 'matrix'),
            b'MDICTP': (self._read_dict, 'matrix'),
            b'KDICTDS': (self._read_dict, 'matrix'),
            b'KDICTX': (self._read_dict, 'matrix'),
            b'XDICT': (self._read_dict, 'matrix'),
            b'XDICTB': (self._read_dict, 'matrix'),
            b'XDICTDS': (self._read_dict, 'matrix'),
            b'XDICTX': (self._read_dict, 'matrix'),

            # coordinate system transformation matrices
            b'CSTM': (self.read_cstm, 'coordinate transforms'),
            b'CSTMS': (self.read_cstm, 'coordinate transforms (superelement)'),
            b'TRMBD': (partial(read_trmbd, self), 'euler angles for transforming from material to (deformed) basic csys'),
            b'TRMBU': (partial(read_trmbu, self), 'euler angles for transforming from material to (undeformed) basic csys'),

            b'R1TABRG': (partial(read_r1tabrg, self), 'DRESP1 optimization table'),
            # Qualifier info table???
            b'QUALINFO': (self.read_qualinfo, 'Qualifier info table'),

            # Equivalence between external and internal grid/scalar numbers
            b'EQEXIN': (self.read_eqexin, 'internal/external ids'),
            b'EQEXINS': (self.read_eqexin, 'internal/external ids (superelement)'),

            b'XSOP2DIR': (self.read_xsop2dir, 'list of external superelement matrices?'),
            b'TUG1': (fread_tug1, 'table Displacement g-set sort 1'),
            b'TEF1': (fread_tug1, 'table element forces sort 1'),
            b'TES1': (fread_tug1, 'table stress sort 1'),
            b'TEE1': (fread_tug1, 'table strain energy sort 1'),
            b'TQMG1': (fread_tug1, 'table mpc forces sort 1'),

            b'MEF1': (fread_mef1, 'external superelement'),
            b'MES1': (fread_mef1, 'external superelement'),
            b'MEE1': (fread_mef1, 'external superelement'),
            b'MEE1O': (fread_mef1, 'external superelement'),
            b'MES1O': (fread_mef1, 'external superelement'),
            b'MEF1O': (fread_mef1, 'external superelement'),
            b'MUG1B': (fread_mef1, 'external superelement'),
            b'MUG1OB': (fread_mef1, 'external superelement'),
            b'MKQG1': (fread_mef1, 'external superelement'),

            b'MMQMG1': (fread_mef1, 'external superelement'),
            b'MBQMG1': (fread_mef1, 'external superelement'),
            b'MKQMG1': (fread_mef1, 'external superelement'),
            b'MK4QMG1': (fread_mef1, 'external superelement'),

            b'MMQG1': (fread_mef1, 'external superelement'),
            b'MBQG1': (fread_mef1, 'external superelement'),
            b'MK4QG1': (fread_mef1, 'external superelement'),

            b'MATPOOL': (fread_matrix_matpool, 'matrices'),

            #b'OBC1': (self.read_obc1, 'Contact pressures and tractions at grid points'),
            #b'OBG1': (self.read_obc1, 'Glue normal and tangential tractions at grid point in cid=0 frame'),
            b'PTMIC': (self._read_ptmic, 'property of VATV microphone points'),

            # OVG: Table of aeroelastic x-y plot data for V-g or V-f curves
            b'MKLIST': (self._read_mklist, 'M/K aero pairs'),
        }
        desc_map = {
            b'PERF': 'aero matrix',
            b'META': 'string matrix',
            b'TUG1': 'external superelement',
            b'TQG1': 'external superelement',
            b'MKQG1': 'external superelement',
            b'MATRV': 'external superelement',
            b'MUG1B': 'external superelement',
            #b'MEF1': 'external superelement',
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
        version_str = ''
        op2: OP2 = self.op2
        read_mode = op2.read_mode
        encoding = self._encoding
        markers = self.get_nmarkers(1, rewind=True)
        #except Exception:
            #self._goto(0)
            #try:
                #self.f.read(4)
            #except Exception:
                #raise FatalError("The OP2 is empty.")
            #raise
        if self.is_debug_file:
            if read_mode == 1:
                self.binary_debug.write(f'read_mode = {read_mode:d} (vectorized; 1st pass)\n')
            elif read_mode == 2:
                self.binary_debug.write(f'read_mode = {read_mode:d} (vectorized; 2nd pass)\n')

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
            data = self.read_string_block(is_interlaced_block=True)  # 'NASTRAN FORT TAPE ID CODE - '
            if data == b'NASTRAN FORT TAPE ID CODE - ':
                macro_version = 'nastran'
            elif data == b'HAJIF FORT TAPE ID CODE -   ':
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
            #version_str = version.decode(encoding)
            #print('version = %r' % version_str)

            if macro_version == 'nastran':
                mode, version_str = parse_nastran_version(
                    data, version, encoding, op2.log)
                if mode == 'msc' and op2.size == 8:
                    op2.is_interlaced = False

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
                self.binary_debug.write(data.decode(encoding) + '\n')
            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -2
            if self.is_debug_file:
                self.binary_debug.write('marker = 2 -> PARAM,POST,-2?\n')
            if op2.post is None:
                op2.post = -2
        else:  # pragma: no cover
            raise NotImplementedError(markers)

        #print(f'mode = {mode!r} fmt={op2._nastran_format!r}')
        if mode == 'autodesk' or op2._nastran_format == 'autodesk':
            op2.post = -4
            mode = 'autodesk'
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
        if read_mode == 1:
            self.log.debug(f'mode={mode!r} version={version_str!r}')
        self.op2.set_mode(mode)
        self.op2.set_table_type()

    def _read_mklist(self):
        """
        reads the MKLIST table and puts it in:
           - op2.op2_results.mklist
        """
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
        if op2.read_mode == 1:
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
        op2: OP2 = self.op2
        self.xsop2dir_names = []
        if op2.read_mode == 1 or op2.make_geom is False:
            table_name = self._read_table_name(rewind=True)
            self._skip_table(table_name, warn=False)
            return

        #op2: OP2 = self.op2
        unused_table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        # (101, 14, 0, 0, 0, 0, 0)
        # (101, 17, 0, 0, 0, 0, 0) # modes_elements_dmig.op2
        data = self._read_record()
        values = np.frombuffer(data, dtype=op2.idtype8)
        self.log.debug(f'xsop2dir: header_values = {values}')
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
        op2: OP2 = self.op2
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        fmt1 = mapfmt(self._endian + b'7i', self.size)

        idata = unpack(fmt1, data)
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
        op2: OP2 = self.op2
        #self.log.debug("table_name = %r" % op2.table_name)
        unused_table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if op2.read_mode == 2:
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
                print('  word=%r (e, f, g)=(%s, %s, %s)' % (word, e, f, g))  # (1, 2, 0)
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
        op2: OP2 = self.op2
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
        if op2.read_mode == 2:
            op2.monitor_data = []

        while markers[0] != 0:

            if op2.read_mode == 1:
                self._skip_record()
            else:
                assert op2.read_mode == 2, op2.read_mode
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
                    'name': name,
                    'cp': cp,
                    'cd': cd,
                    'xyz': [x, y, z],
                    'comps': comps,
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
        op2: OP2 = self.op2
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

    def _read_dict(self) -> None:
        """testing the KDICT"""
        op2: OP2 = self.op2
        if op2.read_mode == 1:
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
        op2: OP2 = self.op2
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
        data = self._read_record()  # CSTM
        #print(self.show_data(data, types='s'))
        assert len(data) == 8 * factor, len(data)

        self.read_3_markers([-3, 1, 0])

        coord_type_map = {
            1: 'CORD2R',
            2: '???',
            3: 'CORD2S',
            5: 'GMSURF',
            7: '???',
            8: '???',
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

        if not op2.read_mode == 1:  # or b'GEOM1' in op2.table_names:
        #if not is_geometry: # or op2.read_mode == 1 or b'GEOM1' in op2.table_names:
            return
        nblocks = len(blocks)
        if nblocks == 1:
            # vectorized
            ints = np.frombuffer(blocks[0], dtype=op2.idtype8)
            floats = np.frombuffer(blocks[0], dtype=op2.fdtype8)
            #doubles = np.frombuffer(blocks[1], dtype='float64')
            nints = len(ints)

            assert nints % 14 == 0, f'nints={nints:d}'
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

            cstm = CSTM()
            cstm.data = np.concatenate([ints, values], axis=1)

            op2.op2_results.cstm = cstm
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

        -100001 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;PVALID=0;DESINC=0;DISCRETE=FALSE;MASSID=0;ARBMID=0;PARTNAME=' ';
                 TRIMID=0;MODULE=0)
        -100000 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;PVALID=0;DESINC=0;ARBMID=0;PARTNAME=' ';DISCRETE=FALSE;TRIMID=0)
         -99999 (AUXMID=0;AFPMID=0;HIGHQUAL=0;PVALID=0;DESINC=0;PRESEQP=TRUE;ARBMID=0;PARTNAME=' ';TRIMID=0;FLXBDYID=0;
                 DFPHASE=' ')
         -99998 (AUXMID=0;AFPMID=0;DESITER=0;HIGHQUAL=0;DESINC=0;DISCRETE=FALSE;ARBMID=0;MASSID=0;PARTNAME=' ';TRIMID=0)
           1431 (HIGHQUAL=0;AUXMID=0;AFPMID=0;DESINC=0;ARBMID=0;PARTNAME=' ';TRIMID=0;FLXBDYID=0)
           1459 (PEID=0;DESITER=0;PVALID=0;NL99=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;AUXMID=0;DESINC=0;DISCRETE=FALSE;
                 MASSID=0;PARTNAME=' ';MODULE=0)
           1461 (PEID=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;AUXMID=0;DESINC=0;PARTNAME=' ';MODULE=0)
           1541 (SEID=0;PEID=0;MTEMP=0;DESITER=0;PVALID=0;APRCH=' ';QCPLD=' ';HIGHQUAL=0;P2G=' ';K2GG=' ';M2GG=' ';
                 DELTA=FALSE;AUXMID=0;BNDSHP=FALSE;ADJOINT=FALSE;DESINC=0;DISCRETE=FALSE;CASEF06=' ';ISOLAPP=1;
                 SUBCID=0;OSUBID=1;STEPID=0;RGYRO=0;PARTNAME=' ';SSTEPID=0;MODULE=0)

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

        op2: OP2 = self.op2
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
        assert ndata == 28, self.show_data(data)  # 7*4

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
            if op2.read_mode == 1:
                db_key, qlen = unpack(self._endian + b'2i', data[:8])
                fmt = self._endian + b'%is' % (ndata - 8)
                qual_str = unpack(fmt, data[8:])[0].decode('latin1')
                qual = QualInfo.from_str(db_key, qual_str.rstrip(), self.log)
                qual_strings[db_key] = qual_str
            itable -= 1
        stop_marker = self.get_marker1(rewind=False)

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
        op2: OP2 = self.op2
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

        if op2.read_mode == 2:
            if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                msg = (
                    'Cannot overwrite op2._frequencies...\n'
                    'op2._frequencies = %s\n'
                    'new_freqs = %s\n' % (op2._frequencies, freqs))
                raise RuntimeError(msg)
            op2._frequencies = freqs
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name_raw!r}, freqs]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
                self.binary_debug.write('  freqs = %s' % freqs)
        self._read_subtables()

    def read_frl(self):
        """
        reads the FRL (Frequency Response List) table

        tested by TestOP2.test_op2_good_sine_01

        """
        op2: OP2 = self.op2
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
            if op2.read_mode == 1:
                self._skip_record()
            else:
                data = self._read_record()
                #self.show_data(data)
                freqs = np.frombuffer(data, dtype=op2.fdtype).copy()
                #print('read_mode=%s itable=%s freqs=%s' % (op2.read_mode, isubtable, freqs.tolist()))
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
        r"""
        reads the RST table (restart file?)

        """
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\k402rerun3d2ss.op2
        op2: OP2 = self.op2
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
            aint, bint, cfloat, spaces, m1a, m1b = unpack('<2qd ' + str(ndata) + 's' + ' 2q', data)
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
            fmt = '<2qd ' + str(ndatai) + 's' + ' q'
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
        op2: OP2 = self.op2
        read_record = get_read_skip_record(self, imode_skip=1)

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
        unused_data = read_record()  # nids 1-117

        self.read_3_markers([-4, 1, 0])
        data = read_record()
        if op2.read_mode == 2 and self.size == 4:
            # nids 1-117 (column 1) with nid*1000 (column 2)
            #
            # External grid or scalar identification number = node_id
            # Sequence number = 1000 * external identification number
            unused_nid_seq = np.frombuffer(data, op2.idtype).reshape(nnodes, 2)
        self.read_markers([-5, 1, 0, 0])

    def read_gpls(self):
        op2: OP2 = self.op2
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
        data, ndata = self._read_record_ndata4()  # GPL
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
        data, ndata = self._read_record_ndata8()  # GPL
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

        is_interlaced_block = self.op2.is_interlaced  # is_nx
        table_name_bytes = reshape_bytes_block(table_name_bytes,
                                               is_interlaced_block=is_interlaced_block)
        table_name_str = table_name_bytes.decode('utf-8').strip()
        assert table_name_str in table_names, f'actual={table_name_str} allowed={table_names}'
        return table_name_str

    def read_ibulk(self):
        """
        tested by TestOP2.test_ibulk

        read_mode = 1 (array sizing)
        read_mode = 1 (reading)

        """
        op2: OP2 = self.op2
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

        op2: OP2 = self.op2
        size = self.size
        lines = []
        if write_deck and op2.read_mode == read_mode:
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

        elif save_lines and op2.read_mode == read_mode:
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
                pass  # raise
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
        op2: OP2 = self.op2
        size = self.size
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')

        if op2.read_mode == 1:
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
        op2: OP2 = self.op2
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
                if op2.read_mode == 1:
                    marker -= 1
                    continue
                #self.show_data(data, types='if', endian=None)
                marker -= 1
            #self.show(200)

        elif method == 2:
            dict_map = {
                1: 'RPM',
                2: 'eigenfreq',
                3: 'Lehr',
                4: 'real_eig',
                5: 'imag_eig',
                6: 'whirl_dir',
                7: 'converted_freq',
                8: 'whirl_code',
            }
            while 1:
                self.read_3_markers([marker, 1, 0])
                nfields = self.get_marker1(rewind=True)
                if nfields == 0:
                    break

                if op2.read_mode == 1:
                    data = self._skip_record()
                    marker -= 1
                    continue

                data = self._read_record()
                ints = np.frombuffer(data, op2.idtype)
                nints = len(ints)
                floats = np.frombuffer(data, op2.fdtype)  # .reshape(nints//7, 7)

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

                i = 0
                data_out = {}
                solution0 = -1
                while i < nints:
                    nvalues = ints[i]
                    ncurves = ints[i+1]
                    keyword = ints[i+2]  # 10000+(SOLN*10)+DATTYP
                    base = keyword - 10000  # (SOLN*10)+DATTYP
                    solution = base // 10
                    datatype = base % 10
                    assert datatype in [1, 2, 3, 4, 5, 6, 7, 8], datatype
                    values = floats[i+3:i+3+nvalues]
                    if datatype in [6, 8]:
                        values = values.astype('int32')
                    datatype_str = dict_map[datatype]
                    #print(f'{nvalues}, {ncurves}, {solution}, {datatype}, {datatype_str:14s} {values}')
                    data_out[datatype_str] = values
                    if solution0 == -1:
                        solution0 = solution
                    assert solution == solution0, (solution, solution0)
                    i += 3 + nvalues
                marker -= 1
                cddata_list.append(data_out)

            #print('------------------------------------')
        else:  # pragma: no cover
            raise NotImplementedError(f'CDDATA method={method}')
        if op2.read_mode == 2:
            #plt.grid(True)
            #plt.show()
            cambell = CampbellData(solution, cddata_list)
            self.op2.op2_results.cddata[solution] = cambell
        nfields = self.get_marker1(rewind=False)

    def read_stdisp(self):
        """reads the STDISP table"""
        #C:\NASA\m4\formats\git\examples\backup\aeroelasticity\loadf.op2
        op2: OP2 = self.op2
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
        #(345, 0,   0,   0,   0,   1347635524,   1230459987,   538986318,      1, 0, 123)
        #(345, 0.0, 0.0, 0.0, 0.0, 14179176448.0, 881989.1875, 1.35761196e-19, 1, 0.0, 1.72359)
        #self.show_data(data[96:], types='ifs', endian=None)

        self.read_3_markers([-4, 1, 0])
        self.read_markers([0])

    def read_omm2(self):
        """reads the OMM2 table"""
        op2: OP2 = self.op2
        #op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_3_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n' % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        else:  # pragma: no cover
            raise NotImplementedError(self.show_data(data))
        self._read_subtables()

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???

        """
        self._skip_pcompts()
        return
        #if op2.read_mode == 1:
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
        op2: OP2 = self.op2
        log = op2.log
        log.debug("table_name = %r" % op2.table_name)
        table_name = self._read_table_name(rewind=False)

        #read_record = self._skip_record if op2.read_mode == 2 else self._read_record
        self.read_markers([-1])

        # (104, 0,   1,   3, 0, 0,   0) - tpl\qrcomp.op2 PCOMPTS (length 3?)
        # (104, 8,   0,   0, 0, 0,   0) - tpl\lmtas1.op2 PCOMPTS (return after -4)
        # (104, 0, 103, 412, 0, 0, 103) - output.op2     PCOMPT  (return at end)
        data_header = self._read_record()
        #self.show_data(data_header, types='ifsd', endian=None)
        a, bi, n4a, n5, e, f, n4b = Struct(b'<7i').unpack(data_header)
        log.debug('a=%s b=%s n4a=%s n5=%s e=%s f=%s n4b=%s' % (
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
        if markers != [0]:  # n4a=0
            #self.show_data(data_header, types='ifsd')
            # (1, 4, 0, '    ')
            #(421, 4, 128, '    ')
            #(814, 4, 256, '    ')
            #---------------------------------
            # (201, 4, 0,   '    ')
            data = self._read_record()
            # self.show_data(data)

            #assert len(data) == n4words*4, 'n4words=%s len(data)=%s n4words*4=%s' % (n4words, len(data), n4words*4)

            if table_name == b'PCOMPTS':
                # (11, 3, 0)
                pass
                #self.show_data(data, types='ifs', endian=None)
            elif 0:  # pramga: no cover
                i = 0
                j = 0
                s1 = Struct(b'< 3i 4s')
                while i < len(data):
                    datai = data[i:i+16]
                    out = s1.unpack(datai)
                    log.debug(str(out))
                    i += 16
                    j += 1
                log.debug("j (-4) = %s" % j)
                # assert i == len(data), '-4'
        else:
            self.read_markers([0])
            assert n4a == 0, n4a
            return

        self.read_3_markers([-5, 1, 0])
        data = self._read_record()
        # assert len(data) == n5words * 4, 'n5words=%s len(data)=%s n5words*4=%s' % (n5words, len(data), n5words*4)

        if table_name == b'PCOMPTS':
            # (101, 102, 103)
            pass
            #self.show_data(data, types='ifs', endian=None)
        elif 0:  # pramga: no cover
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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

    def get_skip_read_record_ndata(self):
        """selects the read_record or skip_record depending on read_mode"""
        if self.op2.read_mode == 1:
            read_record_ndata = self._read_record_ndata
        else:
            read_record_ndata = self._skip_record_ndata
        return read_record_ndata

    def read_sdf(self):
        """reads the SDF table"""
        op2: OP2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_3_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 16:
            subtable_name, dummy_a, dummy_b = unpack(self._endian + b'8sii', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n' % (
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

    def read_long_block(self, expected_marker: int) -> tuple[bytes, int]:
        op2: OP2 = self.op2
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

    def _read_units(self):
        r"""models/msc/units_mass_spring_damper"""
        op2: OP2 = self.op2
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
                encoding = op2._encoding
                out = Struct(self._endian + b'8s 8s 8s 8s').unpack(data)
                mass_bytes, force_bytes, length_bytes, time_bytes = out
                mass = mass_bytes.decode(encoding)
                force = force_bytes.decode(encoding)
                length = length_bytes.decode(encoding)
                time = time_bytes.decode(encoding)
                #print(mass, force, length, time)
                fields = {
                    'mass': mass,
                    'force': force,
                    'length': length,
                    'time': time, }
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
        op2: OP2 = self.op2
        markers = []
        struc = Struct('3i')
        for unused_i in range(nmarkers):
            block = op2.f.read(12)
            marker = struc.unpack(block)
            markers.append(marker)
        return markers

    def get_nmarkers(self, n: int, rewind: bool=True, macro_rewind: bool=False):
        if self.size == 4:
            return self.get_nmarkers4(n, rewind=rewind, macro_rewind=macro_rewind)
        return self.get_nmarkers8(n, rewind=rewind, macro_rewind=macro_rewind)

    def get_nmarkers4(self, n: int, rewind: bool=True, macro_rewind: bool=False):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        Parameters
        ----------
        n : int
            number of markers to get
        rewind : bool; default=True
            should the file be returned to the starting point
        macro_rewind : bool; default=False
            print flag for debugging

        Returns
        -------
        markers : list[int]
            list of [1, 2, 3, ...] markers

        """
        op2: OP2 = self.op2
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

    def get_nmarkers8(self, n: int, rewind: bool=True, macro_rewind: bool=False):
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        #self.log.debug(f'{month}/{day}/{year:d} zero={zero} one={one}')
        month, day, year = self._set_op2_date(month, day, year)

        #if self.is_debug_file:
        if self.is_debug_file:
            self.binary_debug.write(f'  [subtable_name, month={month:d}, day={day:d}, year={year:d}, '
                                    f'zero={zero:d}, one={one:d}]\n\n')
        #assert zero == 0, zero  # is this the RTABLE indicator???
        assert one in [0, 1], one  # 0, 50

    def _set_op2_date(self, month: int, day: int, year: int) -> tuple[int, int, int]:
        """sets the date the job was run"""
        date = (month, day, year)
        self.op2.date = date
        return date

    #----------------------------------------------------------------------------------------
    def _read_record(self, debug: bool=True,
                     macro_rewind: bool=False) -> bytes:
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

    def _read_record_ndata(self, debug: bool=True, macro_rewind: bool=False) -> tuple[bytes, int]:
        """reads a record and the length of the record"""
        if self.size == 4:
            return self._read_record_ndata4(debug=debug, macro_rewind=macro_rewind)
        return self._read_record_ndata8(debug=debug, macro_rewind=macro_rewind)

    def _read_record_ndata4(self, debug: bool=True, macro_rewind: bool=False) -> tuple[bytes, int]:
        """reads a record and the length of the record for size=4"""
        op2: OP2 = self.op2
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
                if op2.allow_empty_records:
                    self.log.error(f'EmptyRecordError: marker0={marker0} nrecord={nrecord}')
                    raise EmptyRecordError('nrecord=4')
                else:
                    self.log.error(f'EmptyRecordError: marker0={marker0} nrecord={nrecord}')
                    raise FortranMarkerError('EmptyRecord: nrecord=4')
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

    def _read_record_ndata8(self, debug: bool=True, macro_rewind: bool=False) -> tuple[bytes, int]:
        """reads a record and the length of the record for size=8"""
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        last_table_name : bytes; default=None
            the last table name

        Returns
        -------
        table_name : bytes
            the table name

        """
        op2: OP2 = self.op2
        is_interlaced_block = op2.is_interlaced

        table_name = None
        data = None
        if self.is_debug_file:
            self.binary_debug.write('_read_table_name - rewind=%s\n' % rewind)
        ni = op2.n
        if stop_on_failure:
            data = self._read_record(debug=False, macro_rewind=rewind)
            if self.size == 8:
                data = reshape_bytes_block(data, is_interlaced_block=is_interlaced_block)
            table_name = self.unpack_table_name(data)

            if self.is_debug_file and not rewind:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write(f'table_header = [8, {table_name!r}, 8]\n\n')
            table_name = table_name.strip()
        else:
            try:
                #data = self.read_string_block(is_interlaced_block=False)
                data = self._read_record(macro_rewind=rewind)
                if self.size == 8:
                    data = reshape_bytes_block(data, is_interlaced_block=is_interlaced_block)
                table_name = self.unpack_table_name(data)
            except (NameError, MemoryError):
                raise
            except Exception:  # struct_error:
                # we're done reading
                op2.n = ni
                op2.f.seek(op2.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0], macro_rewind=rewind)
                except Exception:  # struct_error:
                    # if we hit this block, we have a FATAL error
                    is_special_nastran = op2._nastran_format.lower().startswith(('imat', 'autodesk'))
                    if not is_special_nastran and op2.post != -4:
                        op2.f.seek(op2.n)
                        self.show(1000)
                        if last_table_name:
                            self.log.error(f'finished table_name = {last_table_name}')

                        msg = ('There was a Nastran FATAL Error.  Check the F06.\n'
                               f'last table={op2.table_name!r}; post={op2.post} '
                               f'version={self.op2._nastran_format!r}')
                        self.log.error(msg)
                        if op2.stop_on_unclosed_file:
                            raise FatalError(msg)
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

    def _read_block(self) -> bytes:
        if self.size == 4:
            return self.read_block4()
        return self.read_block8()

    def read_string_block(self, is_interlaced_block: bool) -> bytes:
        block = self._read_block()
        if self.size == 4:
            return block
        return reshape_bytes_block(block, is_interlaced_block)

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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
        data_out = b''
        for unused_i in range(3):
            data = op2.f.read(4)
            ndata, = op2.struct_i.unpack(data)

            data_out += op2.f.read(ndata)
            data = op2.f.read(4)
            op2.n += 8 + ndata
        return data_out

    def read_3_markers(self, markers, macro_rewind: bool=True) -> None:
        """Micro-optimizes ``read_markers`` for 3 markers."""
        if self.size == 4:
            self.read_3_markers4(markers, macro_rewind=macro_rewind)
        else:
            self.read_markers8(markers, macro_rewind=macro_rewind)

    def read_3_markers4(self, markers, macro_rewind: bool=True) -> None:
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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
    def log(self) -> SimpleLogger:
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

    def _skip_table_helper(self, warn: bool=True) -> None:
        """
        Skips the majority of geometry/result tables as they follow a
        very standard format.  Other tables don't follow this format.

        """
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
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

    def _skip_record(self) -> Optional[bytes]:
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

    def _skip_record_ndata(self, debug: bool=True, macro_rewind: bool=False) -> tuple[Optional[bytes], int]:
        if self.size == 4:
            return self._skip_record_ndata4(debug=debug, macro_rewind=macro_rewind)
        return self._skip_record_ndata8(debug=debug, macro_rewind=macro_rewind)

    def _skip_record_ndata4(self, debug: bool=True, macro_rewind: bool=False) -> tuple[Optional[bytes], int]:
        """the skip version of ``_read_record_ndata``"""
        op2: OP2 = self.op2
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
                if op2.allow_empty_records:
                    raise EmptyRecordError('nrecord=4')
                else:
                    raise RuntimeError('EmptyRecord: nrecord=4')
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

    def _skip_record_ndata8(self, debug: bool=True, macro_rewind: bool=False) -> tuple[Optional[bytes], int]:
        """the skip version of ``_read_record_ndata``"""
        op2: OP2 = self.op2
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

    def _get_record_length(self) -> int:
        """
        The record length helps us figure out data block size, which is used
        to quickly size the arrays.  We just need a bit of metadata and can
        jump around quickly.

        Returns
        -------
        record_length : int
            the length of the data block

        """
        op2: OP2 = self.op2
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

    def _skip_block(self) -> Optional[bytes]:
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
        op2: OP2 = self.op2
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
        op2: OP2 = self.op2
        if not op2.is_all_subcases:
            if hasattr(op2, 'isubcase') and op2.isubcase in op2.valid_subcases:
                return True
            return False
        return True

    def read_results_table(self) -> None:
        """Reads a results table"""
        if self.size == 4:
            self.read_results_table4()
        else:
            self.read_results_table8()

    def read_results_table4(self) -> None:
        """Reads a results table"""
        op2: OP2 = self.op2
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
        #try:
        op2._results.log = op2.log
        self._read_subtables()
        #except Exception:
            #self.log.error('\n' + str(self.op2.code_information()))
            #raise

    def read_results_table8(self) -> None:
        """Reads a results table"""
        op2: OP2 = self.op2
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

        is_interlaced_block = self.op2.is_interlaced  # is_nx
        subtable_name = self.get_subtable_name8(op2, data, ndata)
        subtable_name = reshape_bytes_block(
            subtable_name, is_interlaced_block=is_interlaced_block)
        op2.subtable_name = subtable_name
        try:
            self._read_subtables()
        except Exception:
            self.log.error('\n' + str(self.op2.code_information()))
            raise

    def get_subtable_name8(self, op2, data: bytes, ndata: int) -> bytes:
        is_interlaced_block = self.op2.is_nx
        if ndata == 16:  # 8*2
            subtable_name, = op2.struct_16s.unpack(data)
            subtable_name = reshape_bytes_block(
                subtable_name, is_interlaced_block)
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 32:  # 16*2
            #(name1, name2, 170, 170)
            subtable_name, = op2.struct_16s.unpack(data[:16])
            assert len(subtable_name) == 16, len(subtable_name)
            subtable_name = reshape_bytes_block(
                subtable_name, is_interlaced_block=is_interlaced_block)
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 56:  # 28*2
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'16s5q', data)
            subtable_name = reshape_bytes_block(
                subtable_name, is_interlaced_block=is_interlaced_block)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n' % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        else:  # pragma: no cover
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
            subtable_name_bytes, unused_ten = unpack(self._endian + b'8si', data)
            subtable_name = subtable_name_bytes.strip().decode(self._encoding)
            #assert ten == 10, self.show_data(data, types='ifs', endian=None)
            assert subtable_name in ['GPL', 'GPLS'], subtable_name
            if self.is_debug_file:
                self.binary_debug.write(f'  recordi = [{subtable_name!r}]\n')
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
        elif ndata == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n' % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write(f'  subtable_name={subtable_name!r}\n')
            self._print_month(month, day, year, zero, one)
        elif ndata == 612:  # ???
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
        else:  # pragma: no cover
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
        op2: OP2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)

        if self.is_debug_file:
            self.binary_debug.write(f'read_geom_table - {op2.table_name}\n')
        self.read_markers([-1])
        data = self._read_record()  # length=28=7*4

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
        op2: OP2 = self.op2

        # _table4_count is used for numpy streaming
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
        # table3 - the table with the metadata (e.g. subcase_id, time, is_stress/strain)
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
            #if op2.read_mode == 2:
                #self.log.debug("table_name = %r" % table_name)
            try:
                table3_parser, table4_parser, desc = table_mapper[table_name]
            except:
                table3_parser, table4_parser = table_mapper[table_name]
            passer = False
        else:
            if table_name in op2.op2_reader.desc_map:
                desc = op2.op2_reader.desc_map[table_name]

            if op2.read_mode == 2:
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
                raise
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
            #self.log.debug(f'op2.isubtable = {op2.isubtable}')

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
                #self.show(4 * 3 * 3, types='i')
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
            #self.show_ndata(100, types='ifs')
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
        op2: OP2 = self.op2
        table_name = op2.table_name
        IS_TESTING = op2.IS_TESTING
        if self.binary_debug:
            self.binary_debug.write('-' * 60 + '\n')

        # if table_name == b'GEOM2':
        #     self.show(100)
        # this is the length of the current record inside table3/table4
        record_len = self._get_record_length()
        if self.is_debug_file:
            self.binary_debug.write(f'record_length = {record_len:d}\n')

        oes_nl = [b'OESNLXD', b'OESNL1X', b'OESNLXR']  # 'OESCP'?
        factor = self.factor
        #print('record_len =', record_len)
        if record_len == 584 * factor:  # table3 has a length of 584
            if table_name in oes_nl and hasattr(op2, 'num_wide') and op2.num_wide == 146:
                data_code_old = deepcopy(op2.data_code)

            if self.load_as_h5:
                assert self.h5_file is not None, self.h5_file
            op2.data_code = {
                '_encoding': self._encoding,
                'load_as_h5': self.load_as_h5,
                'h5_file': self.h5_file,
                'size': self.size,
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

                        if op2.read_mode == 1:
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
        op2: OP2 = self.op2
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

    def show_ndata(self, n: int, types: str='ifs', force: bool=False, endian=None):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types, force=force, endian=endian)

    def _write_ndata(self, f, n: int, types: str='ifs', force: bool=False, endian=None):  # pragma: no cover
        """Useful function for seeing what's going on locally when debugging."""
        op2: OP2 = self.op2
        nold = op2.n
        data = op2.f.read(n)
        op2.n = nold
        op2.f.seek(op2.n)
        return self._write_data(f, data, types=types, force=force, endian=endian)

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


def eqexin_to_nid_dof_doftype(eqexin1: np.ndarray,
                              eqexin2: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def update_op2_datacode(op2: OP2, data_code_old):
    op2.data_code = data_code_old
    for key, value in data_code_old.items():
        if key == 'size':
            continue
        setattr(op2, key, value)


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
    #for dof in dofs:
        #print('dof', dof)
    #print()
    return dofs


def get_read_skip_record(op2_reader: OP2Reader,
                         imode_skip: int) -> Callable:
    if op2_reader.read_mode == imode_skip:
        read_record = op2_reader._skip_record
    else:
        read_record = op2_reader._read_record
    return read_record


def read_ovg(op2_reader: OP2Reader) -> None:
    """
    V-g or V-f curves
    NX specific table
    """
    op2: OP2 = op2_reader.op2
    #log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size

    result_name = 'flutter_response'
    if 0:
        is_saved_result = op2._results.is_saved(result_name)
        #read_mode = self.read_mode
        if is_saved_result:
            op2._results._found_result(result_name)
        else:
            op2_reader._skip_table('OVG', warn=False)
            return
        is_not_saved_result = not is_saved_result

    #responses = op2.op2_results.responses
    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OVGX    '
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OVGX    ', data

    #print('id,mach,rho,velocity???')
    #print('id?,kfreq?,damping,velocity')
    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    imode_old = 1
    modes = [imode_old]
    datafs = []
    dataf = []
    structi = Struct(endian + b'5i 3f 3i 156s 128s 128s 128s')
    while 1:
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        data = op2_reader._read_record(debug=False)  # table 3
        #self.show_data(data)
        out = structi.unpack(data)

        #1  ACODE(C)    I Device code + 10*Approach Code
        #2  TCODE(C)    I 2002
        #3  METHOD      I Method flag; 1=K, 2=KE, 3=PK, 4=PKNL
        #4  SUBCASE     I Subcase identification number
        #5  POINTID     I Device code + 10*Point identification number
        #6  MACHNUMBER RS Mach number   – METHOD = K, KE, PK, or PKNL
        #7  DENSRATIO  RS Density ratio – METHOD = K, KE, PK, or PKNL
        #8  KFREQ      RS Reduced frequency – METHOD = K
        #9  FCODE       I Format Code = '1'
        #10 NUMWDE      I Number of words per entry in DATA record, set to 4
        #11 MODENUM     I Mode number – METHOD = KE, PK, or PKNL
        # 12 UNDEF(39)    None
        #(60, 2002, 4, 1, 10, 1039199643, 1067257355, 0, 1, 4, 1)
        (acode, tcode, method_int, subcase_id, point_device, mach, rho, kfreq, fcode, numwide, imode,
         b, title, subtitle, subcase) = out

        device_code = acode % 10
        imode10 = point_device
        assert device_code == 0, (acode, device_code)
        #point_id = point_device // 10

        b = b.strip()
        title = title.strip()
        subtitle = subtitle.strip()
        subcase = subcase.strip()

        #print(f'a = {a!r}')
        #print(f'b = {b!r}')
        #print(f'title = {title!r}')
        #print(f'subtitle = {subtitle!r}')
        #print(f'subcase = {subcase!r}')
        assert len(b) == 0, b

        assert acode == 60, acode
        assert tcode == 2002, tcode
        if method_int == 1:
            method = 'K'
        elif method_int == 2:
            method = 'KE'
        elif method_int == 3:
            method = 'PK'
        elif method_int == 4:
            method = 'PKNL'
        else:  # pragma: no cover
            raise NotImplementedError(method_int)

        assert fcode == 1, fcode
        assert numwide == 4, numwide
        try:
            assert imode10 == imode * 10, (imode10, imode)
            assert kfreq == 0.0, kfreq
        except AssertionError:  # pragma: no cover
            print(out)
            raise
        #print(f'method={method}; imode={imode} point_id={point_id}; mach={mach:.3f} rho={rho:.5f}')

        op2_reader.read_3_markers([itable-1, 1, 0])
        data = op2_reader._read_record(debug=False)  # table 4
        fdata = np.frombuffer(data, dtype=op2.fdtype8)
        idata = np.frombuffer(data, dtype=op2.idtype8)
        nvalues = len(fdata) // 4

        if imode != imode_old:
            datafs.append(dataf)
            modes.append(imode)
            imode_old = imode
            dataf = []

        if nvalues == 4:
            velocity, is_complexf, damping, freq = fdata
            is_complex = idata[1]
            assert is_complex in {0, 1}, f'fdata={fdata} idata={idata} is_complex={is_complex}'
            dataf.append([rho, mach, velocity, damping, freq])
        else:
            assert len(fdata) % 4 == 0, len(fdata)
            fdata = fdata.reshape(nvalues, 4)
            idata = idata.reshape(nvalues, 4)
            velocity = fdata[:, 0]
            is_complexf = idata[:, 1]  # flag for if complex modes were output; 0=NO, 1=YES
            damping = fdata[:, 2]      # g
            freq = fdata[:, 3]
            assert is_complexf.min() in {0, 1}, is_complexf.min()
            assert is_complexf.max() in {0, 1}, is_complexf.max()

            #rho = np.ones(numwide, dtype=velocity.dtype) * rho
            #mach = np.ones(numwide, dtype=velocity.dtype) * mach
            for velocityi, dampingi, freqi in zip(velocity, damping, freq):
                dataf.append([rho, mach, velocityi, dampingi, freqi])
        itable -= 2
    datafs.append(dataf)
    op2_reader.read_markers([0])

    # (20, 11, 4) = (nmodes, npoints, [mach, rho, velocity, damping, freq])
    fdata2 = np.array(datafs)

    cref = 1.0
    is_xysym = False
    is_xzsym = False
    if hasattr(op2, 'aero') and op2.aero is not None:
        aero = op2.aero
        cref = aero.cref
        is_xysym = aero.is_symmetric_xy
        is_xzsym = aero.is_symmetric_xz

    if not hasattr(op2, 'in_units'):
        op2.in_units = {
            'density': 'slinch/in^3', 'velocity': 'in/s',
            'altitude': 'ft', 'eas': 'in/s', 'dynamic_pressure': 'psi',
        }
    resp = FlutterResponse.from_nx(
        method, fdata2, op2.op2_filename,
        subcase_id=subcase_id, subtitle=subtitle, cref=cref,
        is_xysym=is_xysym, is_xzsym=is_xzsym,
        in_units=op2.in_units)

    op2.op2_results.vg_vf_response[subcase_id] = resp
    return


def _skip_table(op2_reader: OP2Reader, itable: int) -> None:
    while 1:
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        data = op2_reader._skip_record()  # table 3
        op2_reader.read_3_markers([itable - 1, 1, 0])

        next_marker = op2_reader.get_marker1(rewind=True)
        if next_marker == itable-2:
            op2_reader.read_3_markers([itable-2, 1, 0])
            next_marker2 = op2_reader.get_marker1(rewind=False)
            #log.debug(f'{op2.table_name}; exit on marker={next_marker2}')
            return

        data = op2_reader._read_record()  # table 4
        itable -= 2
    op2_reader.read_markers([0])
    return


def read_oaerotv(op2_reader: OP2Reader) -> None:
    """trim variables"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size
    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROTV '
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAEROTV ', data

    #print('id,mach,rho,velocity???')
    #print('id?,kfreq?,damping,velocity')
    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    #6: f
    #7: f
    #                              Ma q aero ? ? ? c.b.Sref ?            subcase title subtitle
    structi = Struct(endian + b'5i f  f 8s   i i i fff      35i 128s    128s  128s')
    structi2 = Struct(endian + b'8s ii f')
    trim_type_map = {
        0: ('Intercept', ''),
        1: ('Rigid body', 'load/rate'),
        2: ('Control surface', 'RADIANS'),
        3: ('General control', 'RADIANS'),
    }
    trim_status_map = {
        1: 'Free',
        2: 'Free',
        3: 'Fixed',
        4: 'Scheduled',
        5: 'Linked',
        6: 'Scheduled',
    }
    trim = op2.op2_results.trim
    result_name = 'trim.variables'
    is_saved = op2._results.is_saved(result_name)

    while 1:
        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 3
            #ni = -128*3
            #print(op2.show_data(data[:12*4]))
            #print(op2.show_data(data[15*4:ni]))
            out = structi.unpack(data)

            #1  ACODE(C)    I Device code + 10*Approach Code
            #2  TCODE(C)    I 2002
            # 3 DATCOD I Data code = 0
            # 4 SUBCASE I Subcase identification number
            # 5 UNDEF None
            # 6 MACHNUM RS Mach number
            # 7 Q RS Dynamic pressure
            # 8 CONFIG(2) CHAR4 Aerodynamic configuration name
            # 10 NUMWDE I Number of words per entry in DATA, set to 8
            # 11 SYMXY I Aerodynamic configuration XY symmetry
            #     -1 = SYMMETRIC
            #      0 = ASYMMETRIC
            #      1 = ANTISYMMETRIC
            # 12 SYMXZ I Aerodynamic configuration XZ symmetry
            #      -1 = ANTISYMMETRIC
            #       0 = ASYMMETRIC
            #       1 = SYMMETRIC
            # 13 CHORD RS Reference chord length
            # 14 SPAN  RS Reference span length
            # 15 AREA  RS Reference area
            # 16 UNDEF(35) None
            (acode, tcode, method_int, subcase_id,
             point_device, mach, q, aerosg2d, numwide, symxy, symxz,
             cref, bref, sref, *outi,
             title, subtitle, subcase) = out

            op2.subtable_name = ''
            op2.parse_approach_code(data)

            log.debug(f'mach={mach:g} q={q:g} aerosg2d={aerosg2d!r} symxy={symxy} symxz={symxz}\n'
                      f'  cbs_ref=[{cref:g},{bref:g},{sref:g}]')
            if max(outi) != 0 or min(outi) != 0:
                log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

            # device_code = acode % 10
            #imode10 = point_device
            #assert device_code == 0, (acode, device_code)
            #point_id = point_device // 10

            op2.isubcase = subcase_id
            data_code = op2._read_title_helper(data)
            title = data_code['title']
            subtitle = data_code['subtitle']
            label = data_code['label']

            #print(f'title = {title!r}')
            #print(f'subtitle = {subtitle!r}')
            #print(f'subcase = {subcase!r}')
            assert acode == 12, acode
            assert tcode == 103, tcode
            assert numwide == 5, numwide

            op2.tCode = tcode   # trim
            op2.sort_code = 0   # SORT1, real, not-random
            subcase_key = op2._get_code()
        else:
            data = op2_reader._skip_record()  # table 3

        op2_reader.read_3_markers([itable-1, 1, 0])
        next_marker = op2_reader.get_marker1(rewind=True)
        if next_marker == itable-2:
            op2_reader.read_3_markers([itable-2, 1, 0])
            next_marker2 = op2_reader.get_marker1(rewind=False)
            #log.debug(f'{op2.table_name}; exit on marker={next_marker2}')
            return

        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 4
            idata = 0
            trim_values_list = []
            name_trimtype_trimstatus_units_list = []
            while idata*4 < len(data):
                datai = data[idata*4:(idata+5)*4]
                #print(op2.show_data(datai))
                name, trim_type_int, trim_status_int, value = structi2.unpack(datai)
                trim_type, units = trim_type_map[trim_type_int]
                trim_status = trim_status_map[trim_status_int]
                name = name.rstrip().decode('latin1')
                if name == 'INTERCPT':
                    name = 'INTERCEPT'
                    units = ''   # trim_type_int is wrong...
                    # assert units is None, (name, trim_type, units, trim_status)

                elif units == 'load/rate':
                    if name in ['ROLL', 'PITCH', 'YAW']:
                        units = 'NONDIMEN. RATE'
                    elif name in ['ANGLEA', 'SIDES']:
                        units = 'RADIANS'
                    else:
                        units = 'LOAD FACTOR'

                name_trimtype_trimstatus_units_list.append((name, trim_type, trim_status, units))
                trim_values_list.append(value)
                #print(f'name={name!r} ai={ai} bi={bi} value={value}')
                # (b'INTERCPT', 1, 3, 1.0)
                # (b'ANGLEA  ', 1, 1, 0.104)
                # (b'PITCH   ', 1, 3, 0.0)
                # (b'URDD3   ', 1, 3, 2.5)
                # (b'URDD5   ', 1, 3, 0.0)
                # (b'PITCH   ', 1, 3, 0.0)
                # (b'TFLAP   ', 2, 2, -0.45418)
                idata += 5

            op2._results._found_result(result_name)
            name_type_status_units = np.array(name_trimtype_trimstatus_units_list, dtype='U16')
            trim_values_array = np.array(trim_values_list)

            names = name_type_status_units[:, 0]
            names_list = names.tolist()
            nnames = len(names)
            ids = np.full(nnames, -1, dtype='int32')
            if hasattr(op2, 'aestats'):
                for aestat_id, aestat in op2.aestats.items():
                    index = names_list.index(aestat.label)
                    ids[index] = aestat_id
                for aesurf_id, aesurf in op2.aesurf.items():
                    index = names_list.index(aesurf.label)
                    ids[index] = aesurf_id

            trim_vars = TrimVariables(
                mach, q, cref, bref, sref,
                name_type_status_units, trim_values_array,
                ids=ids,
                subcase=subcase_id, title=title,
                subtitle=subtitle, label=label)
            trim_vars.print_f06()
            assert subcase_key not in trim.variables, subcase_key
            log.debug(f'trim.variables: {subcase_key}')
            trim.variables[subcase_key] = trim_vars
        else:
            data = op2_reader._skip_record()  # table 4
        itable -= 2
    op2_reader.read_markers([0])

    #if hasattr(op2, 'aeros'):
        #op2.add_trim(trim_id, mach, q, cref=cref, bref=bref, sref=sref)
        # is_xysym = aero.is_symmetric_xy
        # is_xzsym = aero.is_symmetric_xz
    return


def read_oaerof(op2_reader: OP2Reader) -> None:
    """aero box forces"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size
    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROF  '
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAEROF  ', data

    #print('id,mach,rho,velocity???')
    #print('id?,kfreq?,damping,velocity')
    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    #6: f
    #7: f
    #                              Ma q aero ? ? ? c.b.Sref ?            subcase title subtitle
    structi = Struct(endian + b'5i f  f 8s   i i i fff      35i 128s    128s  128s')
    structi2 = Struct(endian + b'i 4s 6f')

    result_name = 'trim.aero_force'
    is_saved = op2._results.is_saved(result_name)
    while 1:

        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 3
            #ni = -128*3
            #print(op2.show_data(data[:12*4]))
            #print(op2.show_data(data[15*4:ni]))
            out = structi.unpack(data)

            #1  ACODE(C)    I Device code + 10*Approach Code
            #2  TCODE(C)    I 2002
            #3  METHOD      I Method flag; 1=K, 2=KE, 3=PK, 4=PKNL
            #4  SUBCASE     I Subcase identification number
            #5  POINTID     I Device code + 10*Point identification number
            #6  MACH       RS
            #8  Q          RS Dynamic pressure
            #9  config     8d aerodynamic configuration name
            #10 NUMWDE      I Number of words per entry in DATA record, set to 4
            # 11 SYMXY I Aerodynamic configuration XY symmetry
            #  -1 = SYMMETRIC
            #   0 = ASYMMETRIC
            #   1 = ANTISYMMETRIC
            # 12 SYMXZ I Aerodynamic configuration XZ symmetry
            #  -1 = ANTISYMMETRIC
            #   0 = ASYMMETRIC
            #   1 = SYMMETRIC
            # 13 CHORD RS Reference chord length
            # 14 SPAN RS Reference span length
            # 15 AREA RS Reference area
            # 16 UNDEF(39)     None
            #(60, 2002, 4, 1, 10, 1039199643, 1067257355, 0, 1, 4, 1)
            (acode, tcode, method_int, subcase_id,
             point_device, mach, q, aerosg2d, numwide, zero, coord,
             cref, bref, sref, *outi,
             title_bytes, subtitle_bytes, subcase_bytes) = out

            op2._results._found_result(result_name)
            log.debug(f'mach={mach:g} q={q:g} aerosg2d={aerosg2d!r} coord={coord}\n'
                      f'  cbs_ref=[{cref:g},{bref:g},{sref:g}]')
            assert zero == 0, zero
            if max(outi) != 0 or min(outi) != 0:
                log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

            device_code = acode % 10
            #imode10 = point_device
            #assert device_code == 0, (acode, device_code)
            #point_id = point_device // 10

            op2.isubcase = subcase_id
            data_code = op2._read_title_helper(data)
            title = data_code['title'].rstrip()
            subtitle = data_code['subtitle'].rstrip()
            label = data_code['label'].rstrip()

            assert acode == 12, acode
            assert tcode == 102, tcode
            assert numwide == 8, numwide

            op2.tCode = tcode  # trim
            op2.sort_code = 0  # SORT1, real, not-random
            subcase_key = op2._get_code()

            op2_reader.read_3_markers([itable-1, 1, 0])
            data = op2_reader._read_record(debug=False)  # table 4
            idata = 0
            grid_list = []
            label_list = []
            force_list = []
            while idata*4 < len(data):
                datai = data[idata*4:(idata+numwide)*4]
                #print(op2.show_data(datai))
                nid, labeli, *forcei = structi2.unpack(datai)
                #print(nid, label, *force)
                grid_list.append(nid)
                label_list.append(labeli.rstrip().decode('latin1'))
                force_list.append(forcei)
                idata += numwide

            nodes = np.array(grid_list, dtype='int32')
            force = np.array(force_list, dtype='float64')
            force_labels = np.array(label_list)
            assert isinstance(title, str), title
            assert isinstance(subtitle, str), subtitle
            assert isinstance(label, str), label
            aforce = AeroForce(
                subcase_id,
                mach, q, cref, bref, sref,
                nodes, force, force_labels,
                title=title, subtitle=subtitle, label=label)

            assert subcase_key not in op2.op2_results.trim.aero_force
            op2.op2_results.trim.aero_force[subcase_key] = aforce
        else:
            data = op2_reader._skip_record()  # table 3
            op2_reader.read_3_markers([itable-1, 1, 0])
            data = op2_reader._skip_record()  # table 4

        itable -= 2

    op2_reader.read_markers([0])
    # trim_forces = {
    #     'nid': np.array(grid_ids, dtype='int32'),
    #     'label': np.array(labels),
    #     'force': np.array(forces, dtype='float64'),
    # }


def read_oaerop(op2_reader: OP2Reader) -> None:
    """aero box pressures"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size
    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROP  '
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAEROP  ', data

    #print('id,mach,rho,velocity???')
    #print('id?,kfreq?,damping,velocity')
    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    #6: f
    #7: f
    #                              Ma q aero ? ? ? c.b.Sref ?            subcase title subtitle
    structi = Struct(endian + b'5i f  f 8s   i i i fff      35i 128s    128s  128s')
    structi2 = Struct(op2._endian + b'i 4s ff')

    result_name = 'trim.aero_pressure'
    is_saved = op2._results.is_saved(result_name)
    while 1:

        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 3
            # ni = -128*3
            # print(op2.show_data(data[:12*4]))
            # print(op2.show_data(data[15*4:ni]))

            out = structi.unpack(data)

            # 1  ACODE(C)    I Device code + 10*Approach Code
            # 2  TCODE(C)    I 2002
            # 3  METHOD      I Method flag; 1=K, 2=KE, 3=PK, 4=PKNL
            # 4  SUBCASE     I Subcase identification number
            # (60, 2002, 4, 1, 10, 1039199643, 1067257355, 0, 1, 4, 1)
            (acode, tcode, method_int, subcase_id,
             point_device, mach, q, aerosg2d, numwide, zero, coord,
             cref, bref, sref, *outi,
             title_bytes, subtitle_bytes, subcase_bytes) = out
            op2.subtable_name = ''
            op2.parse_approach_code(data)

            log.debug(f'mach={mach:g} q={q:g} aerosg2d={aerosg2d!r} coord={coord}\n'
                      f'  cbs_ref=[{cref:g},{bref:g},{sref:g}]')
            assert zero == 0, zero
            if max(outi) != 0 or min(outi) != 0:
                log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

            device_code = acode % 10
            #imode10 = point_device
            #assert device_code == 0, (acode, device_code)
            #point_id = point_device // 10

            op2.isubcase = subcase_id
            data_code = op2._read_title_helper(data)
            title = data_code['title'].rstrip()
            subtitle = data_code['subtitle'].rstrip()
            label = data_code['label'].rstrip()

            assert acode == 12, acode
            assert tcode == 101, tcode
            assert numwide == 4, numwide

            op2.tCode = tcode  # trim
            op2.sort_code = 0  # SORT1, real, not-random
            subcase_key = op2._get_code()
            log.debug(f'aero_pressure.subcase_key = {subcase_key}')

            op2_reader.read_3_markers([itable-1, 1, 0])
            data = op2_reader._read_record(debug=False)  # table 4

            op2._results._found_result(result_name)

            grid_list = []
            label_list = []
            cp_list = []
            pressure_list = []
            idata = 0
            while idata * 4 < len(data):
                datai = data[idata * 4:(idata + numwide) * 4]
                # print(op2.show_data(datai))
                nid, labeli, aero_pressure_coeff, aero_pressure = structi2.unpack(datai)
                grid_list.append(nid)
                label_list.append(labeli.rstrip().decode('latin1'))
                cp_list.append(aero_pressure_coeff)
                pressure_list.append(aero_pressure)
                idata += numwide

            nodes = np.array(grid_list, dtype='int32')
            cp = np.array(cp_list, dtype='float64')
            pressure = np.array(pressure_list, dtype='float64')
            labels = np.array(label_list)
            assert isinstance(title, str), title
            assert isinstance(subtitle, str), subtitle
            assert isinstance(label, str), label
            # assert isinstance(subcase_id, integer_types), subcase_id
            apress = AeroPressure(
                subcase_id,
                mach, q, cref, bref, sref,
                nodes, cp, pressure,  # labels,
                title=title, subtitle=subtitle, label=label)

            assert subcase_key not in op2.op2_results.trim.aero_force
            log.debug(f'trim.aero_pressure: {subcase_key}')
            op2.op2_results.trim.aero_pressure[subcase_key] = apress
        else:
            data = op2_reader._skip_record()  # table 3
            op2_reader.read_3_markers([itable - 1, 1, 0])
            data = op2_reader._skip_record()  # table 4

        itable -= 2

    op2_reader.read_markers([0])

    #if hasattr(op2, 'aeros'):
        #op2.add_trim(trim_id, mach, q, cref=cref, bref=bref, sref=sref)
        # is_xysym = aero.is_symmetric_xy
        # is_xzsym = aero.is_symmetric_xz
    return


def read_oaeroscd(op2_reader: OP2Reader) -> None:
    """stability & control derivatives"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size

    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROSCD'
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAEROSCD', data

    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    # TODO: QRG incorrectly casts (chord, span, sref) to floats instead of integers
    #       I think this is related to defaults for the AESURF (blank vs. specified)
    #                              Ma q cofnig numwide symxy symxz chord,span,sref zero subcase title subtitle
    #structf = Struct(endian + b'5i f  f 8s     i       i     i     3f              35i  128s    128s  128s')
    structi = Struct(endian + b'5i f  f 8s     i       i     i     3i              35i  128s    128s  128s')
    structi2 = Struct(endian + b'8s 36f')
    result_name = 'trim.derivatives'
    is_saved = op2._results.is_saved(result_name)
    op2 = op2_reader.op2

    # title = ''
    # subtitle = ''
    # label = ''
    while 1:
        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        # $       acsid   rcsid   chord   bref    sref    symxz   symxy
        # AEROS          0       0 1.6137 7.74596       1.       0       0

        # $             id    mach       q
        # TRIM           1     0.5  50000.   URDD3    -50.    Slat      0.
        #             Flap      0.AlrnLeft      0.AlrnRite      0.ElvnLeft      0.
        #         ElvnRite      0.

        #           (acode=12 tcode=106  0  subcase 0  mach q        config      numwide sym1 sym2
        # ints    = (12,      106,       0, 1,      0, 0.5, 50000.0, 'AEROSG2D', 38,     -1, 1,     0, 1, 0)
        # floats  = (12,      106,       0, 1,      0, 0.5, 50000.0, 'AEROSG2D', 38,     -1, 1,     0.0, 1.401298464324817e-45, 0.0)
        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 3
            ni = -128 * 3
            # op2.show_data(data[:15*4], types='ifs')
            # op2.show_data(data[15*4:ni], types='if')
            out = structi.unpack(data)

            # 1 ACODE(C) I Device code + 10* Approach Code = 12
            # 2 TCODE(C) I Table code = 106
            # 3 DATCOD   I Data code = 0
            # 4 SUBCASE  I Subcase identification number
            # 5 UNDEF None
            # 6 MACHNUM      RS Mach number
            # 7 Q            RS Dynamic pressure
            # 8 CONFIG(2) CHAR4 Aerodynamic configuration name
            # 10 NUMWDE       I Number of words per entry in DATA, set to 8   ## QRG is wrong...it's 38
            # 11 SYMXY        I Aerodynamic configuration XY symmetry
            #   -1 = SYMMETRIC
            #    0 = ASYMMETRIC
            #    1 = ANTISYMMETRIC
            # 12 SYMXZ  I Aerodynamic configuration XZ symmetry
            #   -1 = ANTISYMMETRIC
            #    0 = ASYMMETRIC
            #    1 = SYMMETRIC
            # 13 CHORD RS Reference chord length  ## QRG is wrong...this is an int
            # 14 SPAN  RS Reference span length   ## QRG is wrong...this is an int
            # 15 AREA  RS Reference area          ## QRG is wrong...this is an int
            # 16 UNDEF(35) None
            #
            (acode, tcode, method_int, subcase_id,
             point_device, mach, q, aerosg2d, numwide, symxy, symxz,
             chord, span, sref, *outi,
             title_bytes, subtitle_bytes, subcase) = out

            op2.isubcase = subcase_id
            data_code = op2._read_title_helper(data)
            title = data_code['title']
            subtitle = data_code['subtitle']
            label = data_code['label']

            allowed_cbs = [
                (0, 0, 0),
                (0, 1, 0),
            ]
            if (chord, span, sref) not in allowed_cbs:
                log.error(f'Expected {op2.table_name} (chord,span,sref) flags can be {allowed_cbs}; got ({chord},{span},{sref})')

            log.debug(f'mach={mach:g} q={q:g} aerosg2d={aerosg2d!r} symxy={symxy}; symxz={symxz}')
            assert numwide == 38, numwide
            if max(outi) != 0 or min(outi) != 0:
                log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

            device_code = acode % 10
            #imode10 = point_device
            #assert device_code == 0, (acode, device_code)
            #point_id = point_device // 10

            # title = title.strip()
            # subtitle = subtitle.strip()
            # subcase = subcase.strip()
            # label = ''

            #print(f'title = {title!r}')
            #print(f'subtitle = {subtitle!r}')
            #print(f'subcase = {subcase!r}')

            assert acode == 12, acode
            assert tcode == 106, tcode
            assert numwide == 38, numwide

            op2.tCode = tcode  # trim
            op2.sort_code = 0  # SORT1, real, not-random
            subcase_key = op2._get_code()

            op2_reader.read_3_markers([itable - 1, 1, 0])
            data = op2_reader._read_record(debug=False)  # table 4

            op2._results._found_result(result_name)
            names = []
            all_values = []
            idata = 0
            while idata*4 < len(data):
                datai = data[idata*4:(idata+numwide)*4]
                #print(op2.show_data(datai))
                name_bytes, *data_list = structi2.unpack(datai)
                name = name_bytes.rstrip(b' ').decode(op2._encoding)
                data_values = np.array(data_list, dtype='float32')
                values = data_values.reshape(6, 6)

                names.append(name)
                all_values.append(values)
                # print(values)
                # name = [REFCOEFF, ANGLEA, ELVNLEFT, ELVNRITE, URDD3, ALRNLEFT, SLAT, FLAP, ALRNRITE]
                idata += numwide

            # create the output---
            nnames = len(names)
            derivatives_array = np.array(all_values).reshape(nnames, 6, 6)
            # print(derivatives_array[0, :, :])
            assert derivatives_array.shape == (nnames, 6, 6), (nnames, derivatives_array.shape)

            names_array = np.array(names)
            trim_derivatives = TrimDerivatives(
                mach, q, chord, span, sref,
                names_array, derivatives_array,
                subcase=subcase_id, title=title, subtitle=subtitle, label=label)
            trim = op2.op2_results.trim
            assert subcase_key not in trim.derivatives, subcase_key
            trim.derivatives[subcase_key] = trim_derivatives
        else:
            data = op2_reader._skip_record()  # table 3
            op2_reader.read_3_markers([itable-1, 1, 0])
            data = op2_reader._skip_record()  # table 4
        itable -= 2

    op2_reader.read_markers([0])

    #if hasattr(op2, 'aeros'):
        #op2.add_trim(trim_id, mach, q, cref=cref, bref=bref, sref=sref)
        # is_xysym = aero.is_symmetric_xy
        # is_xzsym = aero.is_symmetric_xz
    return


def read_oaercshm(op2_reader: OP2Reader) -> None:
    """control surface position and hinge moment"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size

    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROSCD'
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAERCSHM', data

    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    #                              Ma q aero ? ? ? cbs_ref zero subcase title subtitle
    structi = Struct(endian + b'5i f  f 8s   i i i 3f      35i  128s    128s  128s')

    structi2 = Struct(endian + b'8s 6f')
    trim = op2.op2_results.trim

    result_name = 'trim.control_surface_position_hinge_moment'
    is_saved = op2._results.is_saved(result_name)
    while 1:
        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        if is_saved:
            data = op2_reader._read_record(debug=False)  # table 3
            #ni = -128*3
            #print(op2.show_data(data[:12*4]))
            #print(op2.show_data(data[15*4:ni]))
            out = structi.unpack(data)

            #1  ACODE(C)    I Device code + 10*Approach Code
            #2  TCODE(C)    I 2002
            #3  METHOD      I Method flag; 1=K, 2=KE, 3=PK, 4=PKNL
            #4  SUBCASE     I Subcase identification number
            #5  POINTID     I Device code + 10*Point identification number
            #6  MACH       RS
            #8  KFREQ      RS Reduced frequency – METHOD = K
            #9  FCODE       I Format Code = '1'
            #10 NUMWDE      I Number of words per entry in DATA record, set to 4
            #11 MODENUM     I Mode number – METHOD = KE, PK, or PKNL
            # 12 UNDEF(39)    None
            (acode, tcode, method_int, subcase_id,
             point_device, mach, q, aerosg2d, numwide, zero, coord,
             cref, bref, sref, *outi,
             title, subtitle, subcase) = out

            op2.subtable_name = ''
            op2.parse_approach_code(data)

            log.debug(f'mach={mach:g} q={q:g} aerosg2d={aerosg2d!r} coord={coord}')
            log.debug(f'  cbs_ref=[{cref:g},{bref:g},{sref:g}]')
            #assert zero == 0, zero
            if max(outi) != 0 or min(outi) != 0:
                log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

            device_code = acode % 10
            #imode10 = point_device
            #assert device_code == 0, (acode, device_code)
            #point_id = point_device // 10

            op2.isubcase = subcase_id
            data_code = op2._read_title_helper(data)
            title = data_code['title']
            subtitle = data_code['subtitle']
            label = data_code['label']

            assert acode == 12, acode
            assert tcode == 104, tcode
            assert numwide == 8, numwide

            op2.tCode = tcode  # trim
            op2.sort_code = 0  # SORT1, real, not-random
            subcase_key = op2._get_code()
        else:
            data = op2_reader._skip_record()  # table 3

        op2_reader.read_3_markers([itable-1, 1, 0])
        next_marker = op2_reader.get_marker1(rewind=True)
        if next_marker == itable-2:  # -5
            op2_reader.read_3_markers([itable-2, 1, 0], macro_rewind=False)
            next_marker2 = op2_reader.get_marker1(rewind=False)
            #log.debug(f'{op2.table_name}; exit on marker={next_marker2}')
            return

        if is_saved:
            #op2_reader.show(80, types='ifs')
            data = op2_reader._read_record(debug=False)  # table 4
            idata = 0
            name_map = {'INTERCPT': 'INTERCEPT'}

            names_list = []
            trim_values_list = []
            data_list = []
            while idata*4 < len(data):
                # TFLAP [-4.5418206e-01 -1.5707964e+00  1.5707964e+00  1.6729131e+06
                #        -1.0000000e+10  1.0000000e+10]
                datai = data[idata*4:(idata+numwide)*4]
                name, trim_value, *data_listi = structi2.unpack(datai)
                name = name.rstrip().decode(op2.encoding)
                if name == 'INTERCPT':
                    name = 'INTERCEPT'
                # name = name_map.get(name, name)
                names_list.append(name)
                trim_values_list.append(trim_value)
                data_list.append(data_listi)
                # log.debug(f'{name}={trim_value:g} values={values.round(4)}')
                idata += numwide

            op2._results._found_result(result_name)

            names = np.array(names_list, dtype='U10')
            trim_values = np.array(trim_values_list)
            data_array = np.array(data_list)
            trim_control_surface_position_hinge_moment = ControlSurfacePositionHingeMoment(
                mach, q, cref, bref, sref,
                names, trim_values, data_array,
                subcase=subcase_id, title=title,
                subtitle=subtitle, label=label)
            assert subcase_key not in trim.control_surface_position_hinge_moment, subcase_key
            trim.control_surface_position_hinge_moment[subcase_key] = trim_control_surface_position_hinge_moment
        else:
            data = op2_reader._skip_record()  # table 4
        itable -= 2
    op2_reader.read_markers([0])

    #if hasattr(op2, 'aeros'):
        #op2.add_trim(trim_id, mach, q, cref=cref, bref=bref, sref=sref)
        # is_xysym = aero.is_symmetric_xy
        # is_xzsym = aero.is_symmetric_xz
    return


def read_oaerohmd(op2_reader: OP2Reader) -> None:
    """hinge moment derivatives"""
    op2: OP2 = op2_reader.op2
    log = op2.log
    endian = op2_reader._endian
    size = op2_reader.size
    #factor: int = op2_reader.factor
    assert size == 4, size

    op2.table_name = op2_reader._read_table_name(rewind=False)

    #if op2.read_mode == 1 or is_not_saved_result:
    op2_reader.read_markers([-1])

    #(101, 1, 0, 0, 0, 0, 0)
    data = op2_reader._read_record(debug=False)

    #n = 8
    #'OAEROHMD'
    #op2_reader._skip_record()
    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record(debug=False)
    assert data == b'OAEROHMD', data

    itable = -3
    if op2_reader.read_mode == 1:
        _skip_table(op2_reader, itable)
        return

    #                              Ma q aero ? ? ? name one one zero subcase title subtitle
    structi = Struct(endian + b'5i f  f 8s   i i i 8s   i   i   34i  128s    128s  128s')
    structi2 = Struct(endian + b'8s 5f')
    trim = op2.op2_results.trim

    result_name = 'trim.hinge_moment_derivatives'
    is_saved = op2._results.is_saved(result_name)

    while 1:
        #       trimid    coord
        # AEROS ACSID RCSID       REFC      REFB      REFS SYMXZ SYMXY
        # AEROS   1       1       131.0   2556.4  734000.01       0
        op2_reader.read_3_markers([itable, 1, 0])
        itablei = op2_reader.get_marker1()
        if itablei == 0:
            break

        data = op2_reader._read_record(debug=False)  # table 3
        #ni = -128*3
        #op2.show_data(data[:12*4])
        #op2.show_data(data[15*4:ni])
        out = structi.unpack(data)

        #1  ACODE(C)      I Device code + 10*Approach Code
        #2  TCODE(C)      I 105
        # 3 DATCOD        I Data code = 0
        # 4 SUBCASE       I Subcase identification number
        # 5 UNDEF None
        # 6 MACHNUM      RS Mach number
        # 7 Q            RS Dynamic pressure
        # 8 CONFIG(2) CHAR4 Aerodynamic configuration name
        # 10 NUMWDE I Number of words per entry in DATA, set to 8
        # 11 SYMXY  I Aerodynamic configuration XY symmetry
        # -1 = SYMMETRIC
        #  0 = ASYMMETRIC
        #  1 = ANTISYMMETRIC
        # 12 SYMXZ  I Aerodynamic configuration XZ symmetry
        # -1 = ANTISYMMETRIC
        #  0 = ASYMMETRIC
        #  1 = SYMMETRIC
        # 13 CHORD          RS Reference chord length
        # 14 SPAN           RS Reference span length
        # 15 AREA           RS Reference area
        # 16 CNTLSURF(2) CHAR4 Control surface
        # 18 REFCORDL       RS Reference chord length
        # 19 REFAREA        RS Reference area
        # 20 UNDEF(31) None

        nzero = 34
        # op2.show_data(data[:-nzero*4-128*3])

        # op2.show_data(data[4*10:-nzero*4-128*3])
        (acode, tcode, method_int, subcase_id,
         point_device, mach, q, aerosg2d_bytes, numwide, symxy, symxz,
         cs_name_bytes, one_a, one_b, *outi,
         title, subtitle, subcase) = out
        # print(f'cs_name = {cs_name_bytes}')

        op2.subtable_name = ''
        op2.parse_approach_code(data)

        #op2.show_data(data[14*4:15*4])
        cs_name = cs_name_bytes.decode('latin1').rstrip()
        aerosg2d = aerosg2d_bytes.decode('latin1').rstrip()
        log.debug(f'mach={mach:g} q={q:.3f} cs_name={cs_name!r} aerosg2d={aerosg2d!r} symxy={symxy}; symxz={symxz}')
        #log.debug(f'  name=[{name}]')

        allowed = [
            (1, 1),
            (1, 0),
        ]
        assert (one_a, one_b) in allowed, f'Allowed={allowed}; actual=({one_a:d} ,{one_b:d})'

        if max(outi) != 0 or min(outi) != 0:
            log.error(f'Expected all 0s in {op2.table_name}; outi={outi}')

        # device_code = acode % 10
        #imode10 = point_device
        #assert device_code == 0, (acode, device_code)
        #point_id = point_device // 10

        op2.isubcase = subcase_id
        data_code = op2._read_title_helper(data)
        title = data_code['title']
        subtitle = data_code['subtitle']
        label = data_code['label']

        assert acode == 12, acode
        assert tcode == 105, tcode
        assert numwide == 7, numwide

        op2.tCode = tcode  # trim
        op2.sort_code = 0  # SORT1, real, not-random
        subcase_key = op2._get_code()

        op2_reader.read_3_markers([itable-1, 1, 0])
        data = op2_reader._read_record(debug=False)  # table 4
        idata = 0

        # 1 LABEL(2) CHAR4 Trim variable name
        # 3 RIGID       RS Rigid hinge moment derivative
        # 4 ELSTRES     RS Elastic    restrained hinge moment derivative
        # 5 ELSTURSTN   RS Elastic  unrestrained hinge moment derivative
        # 6 INRLRES     RS Inertial   restrained hinge moment derivative
        # 7 INRLURSTN   RS Inertial unrestrained hinge moment derivative
        # trim_hinge_moment_derivatives = {}
        if is_saved:
            op2._results._found_result(result_name)

            names_list = []
            values_list = []
            while idata*4 < len(data):
                # TFLAP [-4.5418206e-01 -1.5707964e+00  1.5707964e+00  1.6729131e+06
                #        -1.0000000e+10  1.0000000e+10]
                datai = data[idata*4:(idata+numwide)*4]
                name, *data_listi = structi2.unpack(datai)
                name = name.rstrip().decode(op2.encoding)
                if name == 'INTERCPT':
                    name = 'INTERCEPT'
                names_list.append(name)
                values_list.append(data_listi)
                # log.debug(f'HMD {name}: values={values.round(6)}')
                idata += numwide

            nnames = len(names_list)
            names = np.array(names_list, dtype='U8')
            values = np.array(values_list)
            assert values.shape == (nnames, 5), (nnames, values.shape)
            chord = 0.
            span = 0.
            sref = 0.
            derivs = HingeMomentDerivatives(
                mach, q, chord, span, sref,
                cs_name, names, values, subcase=subcase_id,
                title=title, subtitle=subtitle, label=label)

            # print(cs_name, names)
            # print(f'subcase_key = {subcase_key}')
            assert subcase_key[-1] == '', subcase_key
            subcase_list = list(subcase_key)
            subcase_list[-1] = cs_name
            subcase_key2 = tuple(subcase_list)

            assert subcase_key2 not in trim.hinge_moment_derivatives, subcase_key2
            #print(subcase_key2)
            # trim.hinge_moment_derivatives[subcase_key] = trim_hinge_moment_derivatives
            trim.hinge_moment_derivatives[subcase_key2] = derivs

        itable -= 2
    op2_reader.read_markers([0])
    #if hasattr(op2, 'aeros'):
        #op2.add_trim(trim_id, mach, q, cref=cref, bref=bref, sref=sref)
        # is_xysym = aero.is_symmetric_xy
        # is_xzsym = aero.is_symmetric_xz
    return
