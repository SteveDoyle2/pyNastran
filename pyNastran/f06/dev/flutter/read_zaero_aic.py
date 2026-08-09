"""Reader and writer for ZAERO binary AIC (Aerodynamic Influence Coefficient) files.

Binary format (Fortran unformatted sequential access):
  - Each record: [4-byte marker (nbytes)] [data] [4-byte marker (nbytes)]
  - All numeric data is little-endian double (float64) or int64.

File layout:
  1. Header record: [unused, mach, nk_int64, freq_0, freq_1, ..., freq_{nk-1}]
  2. For each symmetry condition (symmetric, then antisymmetric if sym_flag=2):
     For each reduced frequency k_i:
       a. Dim record (5 int64): [sym_flag, nrows, ncols, method, nmach]
       b. nrows data records: [row_idx, ncols, re_0, im_0, re_1, im_1, ...]
       c. Map header (5 int64): [3, nboxes, ndof, method, nmach]
       d. nboxes map records: DOF mapping (40 or 64 bytes depending on k)
       e. Zone header (5 int64): [3, nboxes, ndof, method, nmach]
       f. nboxes zone records (16 bytes): [1, 0]
"""

from __future__ import annotations

import os
import struct
from dataclasses import dataclass, field
from typing import BinaryIO

import numpy as np

from pyNastran.utils import PathLike, print_bad_path


@dataclass
class ZaeroAicResult:
    mach: float
    reduced_freqs: np.ndarray
    aic_matrices: dict[str, np.ndarray]
    nrows: int
    ncols: int
    sym_flag: int
    method: int
    nmach: int
    map_records: dict[str, list[bytes]] = field(repr=False)
    zone_records: dict[str, list[bytes]] = field(repr=False)
    map_headers: dict[str, np.ndarray] = field(repr=False)
    zone_headers: dict[str, np.ndarray] = field(repr=False)
    footer_records: list[bytes] = field(default_factory=list, repr=False)


def _read_record(data: bytes, pos: int) -> tuple[bytes, int]:
    nbytes = struct.unpack_from("<i", data, pos)[0]
    rec = data[pos + 4 : pos + 4 + nbytes]
    return rec, pos + 4 + nbytes + 4


def _write_record(f: BinaryIO, rec: bytes) -> None:
    marker = struct.pack("<i", len(rec))
    f.write(marker)
    f.write(rec)
    f.write(marker)


def read_zaero_aic(aic_filename: PathLike) -> ZaeroAicResult:
    """Read a ZAERO binary AIC file.

    Parameters
    ----------
    aic_filename : PathLike
        Path to the ZAERO AIC binary file.

    Returns
    -------
    ZaeroAicResult
        Parsed AIC data including Mach number, reduced frequencies,
        and complex AIC matrices keyed by symmetry label.
    """
    assert os.path.exists(aic_filename), print_bad_path(aic_filename)
    with open(aic_filename, "rb") as f:
        data = f.read()

    pos = 0

    rec_header, pos = _read_record(data, pos)
    header = np.frombuffer(rec_header, dtype="<f8")
    mach = header[1]
    nk = np.frombuffer(rec_header[16:24], dtype="<i8")[0]
    reduced_freqs = header[3 : 3 + nk].copy()

    aic_matrices: dict[str, np.ndarray] = {}
    map_records: dict[str, list[bytes]] = {}
    zone_records: dict[str, list[bytes]] = {}
    map_headers: dict[str, np.ndarray] = {}
    zone_headers: dict[str, np.ndarray] = {}

    sym_labels = ["symmetric", "antisymmetric"]
    sym_idx = 0
    nrows = ncols = sym_flag = method = nmach = 0
    nrows_expected = -1

    footer_records: list[bytes] = []
    while pos < len(data):
        peek_rec, peek_pos = _read_record(data, pos)
        peek_dims = np.frombuffer(peek_rec, dtype="<i8")

        if len(peek_dims) == 5 and nrows_expected > 0 and int(peek_dims[1]) != nrows_expected:
            pos = peek_pos
            footer_records.append(peek_rec)
            while pos < len(data):
                rec_extra, pos = _read_record(data, pos)
                footer_records.append(rec_extra)
            break

        label = sym_labels[sym_idx] if sym_idx < len(sym_labels) else f"sym_{sym_idx}"
        aic_3d = None

        for ik in range(nk):
            rec_dim, pos = _read_record(data, pos)
            dims = np.frombuffer(rec_dim, dtype="<i8")
            sym_flag, nrows, ncols, method, nmach = dims[0], dims[1], dims[2], dims[3], dims[4]

            if ik == 0:
                nrows_expected = nrows
                aic_3d = np.zeros((nk, nrows, ncols), dtype=complex)

            for irow in range(nrows):
                rec_data, pos = _read_record(data, pos)
                vals = np.frombuffer(rec_data, dtype="<f8")
                real_parts = vals[2::2]
                imag_parts = vals[3::2]
                aic_3d[ik, irow, :] = real_parts[:ncols] + 1j * imag_parts[:ncols]

            rec_map_hdr, pos = _read_record(data, pos)
            map_hdr = np.frombuffer(rec_map_hdr, dtype="<i8").copy()

            nboxes = int(map_hdr[1])
            map_recs: list[bytes] = []
            for _ in range(nboxes):
                rec_map, pos = _read_record(data, pos)
                map_recs.append(rec_map)

            rec_zone_hdr, pos = _read_record(data, pos)
            zone_hdr = np.frombuffer(rec_zone_hdr, dtype="<i8").copy()

            zone_recs: list[bytes] = []
            for _ in range(nboxes):
                rec_zone, pos = _read_record(data, pos)
                zone_recs.append(rec_zone)

            freq_key = f"{label}_k{ik}"
            map_records[freq_key] = map_recs
            zone_records[freq_key] = zone_recs
            map_headers[freq_key] = map_hdr
            zone_headers[freq_key] = zone_hdr

        aic_matrices[label] = aic_3d
        sym_idx += 1

    return ZaeroAicResult(
        mach=mach,
        reduced_freqs=reduced_freqs,
        aic_matrices=aic_matrices,
        nrows=nrows,
        ncols=ncols,
        sym_flag=sym_flag,
        method=method,
        nmach=nmach,
        map_records=map_records,
        zone_records=zone_records,
        map_headers=map_headers,
        zone_headers=zone_headers,
        footer_records=footer_records,
    )


def write_zaero_aic(aic_filename: PathLike, result: ZaeroAicResult) -> None:
    """Write a ZAERO binary AIC file.

    Parameters
    ----------
    aic_filename : PathLike
        Output path for the binary AIC file.
    result : ZaeroAicResult
        AIC data to write.
    """
    nk = len(result.reduced_freqs)

    header_vals = np.zeros(3 + nk, dtype="<f8")
    header_vals[0] = 0.0
    header_vals[1] = result.mach
    header_vals[2] = 0.0
    header_vals[3 : 3 + nk] = result.reduced_freqs
    header_bytes = bytearray(header_vals.tobytes())
    struct.pack_into("<q", header_bytes, 16, nk)

    with open(aic_filename, "wb") as f:
        _write_record(f, bytes(header_bytes))

        for sym_label, aic_3d in result.aic_matrices.items():
            for ik in range(nk):
                dims = np.array(
                    [result.sym_flag, result.nrows, result.ncols, result.method, result.nmach],
                    dtype="<i8",
                )
                _write_record(f, dims.tobytes())

                for irow in range(result.nrows):
                    row_data = np.zeros(2 + 2 * result.ncols, dtype="<f8")
                    row_header = np.array([irow + 1, result.ncols], dtype="<i8")
                    row_data[:2] = np.frombuffer(row_header.tobytes(), dtype="<f8")
                    row_data[2::2] = aic_3d[ik, irow, :].real
                    row_data[3::2] = aic_3d[ik, irow, :].imag
                    _write_record(f, row_data.tobytes())

                freq_key = f"{sym_label}_k{ik}"
                _write_record(f, result.map_headers[freq_key].tobytes())
                for rec in result.map_records[freq_key]:
                    _write_record(f, rec)

                _write_record(f, result.zone_headers[freq_key].tobytes())
                for rec in result.zone_records[freq_key]:
                    _write_record(f, rec)

        for rec in result.footer_records:
            _write_record(f, rec)
