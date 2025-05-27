"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import Any, BinaryIO, TYPE_CHECKING

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues
from .edt_writer import remove_unsupported_cards

if TYPE_CHECKING:  # pragma: no cover
    #from cpylog import SimpleLogger
    #from pyNastran.bdf.cards.aero.aero import CAERO1, CAERO2, PAERO1, PAERO2, AESURF, AESURFS
    #from pyNastran.bdf.cards.aero.static_loads import AEROS, AESTAT # , CSSCHD, DIVERG, TRIM, TRIM2
    from pyNastran.bdf.cards.aero.dynamic_loads import GUST
    from pyNastran.op2.op2_geom import OP2Geom, BDF

def write_dit(op2_file: BinaryIO, op2_ascii,
              model: BDF | OP2Geom,
              endian: bytes=b'<', nastran_format: str='nx') -> None:
    """writes the DIT/DITS table"""
    if not hasattr(model, 'loads'):  # OP2
        return
    card_types = [
        'GUST',
        'TABDMP1',
        'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
        'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
        'TABRNDG',
        'TABLES1', 'TABLEST',
        'TABLEH1',
        #'TABLEHT',
    ]

    cards_to_skip = []
    out = defaultdict(list)
    # geometry
    for table_id, table in sorted(model.tables.items()):
        out[table.type].append(table_id)
    for table_id, table in sorted(model.tables_d.items()):
        out[table.type].append(table_id)
    for table_id, table in sorted(model.tables_m.items()):
        out[table.type].append(table_id)
    for table_id, table in sorted(model.tables_sdamping.items()):
        out[table.type].append(table_id)
    for table_id, table in sorted(model.random_tables.items()):
        out[table.type].append(table_id)
    for gust_id, gust in sorted(model.gusts.items()):
        out[gust.type].append(gust_id)

    remove_unsupported_cards(out, card_types, model.log)
    # other
    if len(out) == 0:
        return

    write_geom_header(b'DIT', op2_file, op2_ascii, endian=endian)
    itable = -3

    for name, ids in sorted(out.items()):
        model.log.debug('DIT %s %s' % (name, ids))
        ncards = len(ids)
        assert ncards > 0, ncards
        if name in cards_to_skip:
            model.log.warning('skipping EDT-%s' % name)
            continue

        #if nmaterials == 0:
            #continue
        #model.log.info(f'EDT: {name}; n={len(ids)}')
        try:
            func = DIT_MAP[name]
        except KeyError:  # pragma: no cover
            #continue
            raise NotImplementedError(name)

        nbytes = func(model, name, ids, ncards, op2_file, op2_ascii,
                      endian, nastran_format)
        op2_file.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2_file.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2_file, op2_ascii, itable)


def write_tabdmp1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file: BinaryIO,
                  op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    TABDMP1(15, 21, 162)

    1 ID    I  Table identification number
    9 F     RS Natural frequency
    10 G    RS Damping
    Words 9 through 10 repeat until (-1,-1) occurs
    """
    key = (15, 21, 162)
    nvalues = 0
    datas = []

    fmt = endian
    for table_id in table_ids:
        table = model.tables_sdamping[table_id]
        nx = len(table.x)
        fmt += b'8i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)
        data = [table_id, 0, 0, 0, 0, 0, 0, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  TABDMP1 data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tabrndg(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file: BinaryIO, op2_ascii,
                  endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    TABRNDG(56, 26, 303)
    Power spectral density for gust loads in aeroelastic analysis.

    1 ID        I   Table identification number
    2 TYPE      I   Power spectral density type
    3 LU        RS  Scale of turbulence divided by velocity
    4 WG        RS  Root-mean-square gust velocity
    5 UNDEF(4) none Not used
    Words 1 through 8 repeat until (-1,-1) occurs

    """
    key = (56, 26, 303)
    nvalues = 0
    datas = []

    fmt = endian + b'2i 2f 4i 2i' * ncards
    for table_id in table_ids:
        table = model.random_tables[table_id]

        data = (table_id, table.Type,
                table.LU, table.WG,
                0, 0, 0, 0, -1, -1)
        assert None not in data, data
        #assert (unused_dunno_a, unused_dunno_b, unused_dunno_c, unused_dunno_d) == (0, 0, 0, 0), out
        #if tid > 100000000:
            #tid = -(tid - 100000000)
        #op2.add_tabrndg(tid, table_type, lu, wg, comment='')
        assert table_id > 0, table
        datas.extend(data)
        op2_ascii.write(f'  TABDMP1 data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tabled1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table1(
        model, name, table_ids, ncards,
        model.tables_d, 'TABLED1',
        (1105, 11, 133),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tablem1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table1(
        model, name, table_ids, ncards,
        model.tables_m, 'TABLEM1',
        (105, 1, 93),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tabled2(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table2(
        model, name, table_ids, ncards,
        model.tables_d, 'TABLED2',
        (1205, 12, 134),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tablem2(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table2(
        model, name, table_ids, ncards,
        model.tables_m, 'TABLEM2',
        (205, 2, 94),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tabled3(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table3(
        model, name, table_ids, ncards,
        model.tables_d, 'TABLED3',
        (1305, 13, 140),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tablem3(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table3(
        model, name, table_ids, ncards,
        model.tables_m, 'TABLEM3',
        (305, 3, 95),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tabled4(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table4(
        model, name, table_ids, ncards,
        model.tables_d, 'TABLED4',
        (1405, 14, 141),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_tablem4(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    nbytes = _write_table4(
        model, name, table_ids, ncards,
        model.tables_m, 'TABLEM4',
        (405, 4, 96),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def _write_table1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  table_dict: dict[int, Any], table_name: str,
                  key: tuple[int, int, int],
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    TABLED1(1105,11,133)

    Word Name Type Description
    1 ID     I Table identification number
    2 CODEX  I Type of interpolation for the x-axis
    3 CODEY  I Type of interpolation for the y-axis
    4 FLAG   I Extrapolation on/off flag
    5 LOCUT RS Low cutoff value
    6 HICUT RS High cutoff value
    7 UNDEF(2) None
    9  X    RS X tabular value
    10 Y    RS Y tabular value
    Words 9 through 10 repeat until (-1,-1) occurs
    """
    nvalues = 0
    datas = []

    fmt = endian
    locut = 0.
    hicut = 0.
    for table_id in table_ids:
        table = table_dict[table_id]
        nx = len(table.x)
        fmt += b'4i 2f 2i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        data = [table_id]
        for axis_type in [table.xaxis, table.yaxis]:
            if axis_type == 'LINEAR':
                axis = 0
            elif axis_type == 'LOG':
                axis = 1
            else: # pragma: no cover
                raise ValueError('axis=%r' % axis_type)
            data.append(axis)

        data += [table.extrap, locut, hicut, 0, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  {table_name} data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def _write_table2(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  table_dict: dict[int, Any], table_name: str,
                  key: tuple[int, int, int],
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    1 ID       I  Table identification number
    2 X1      RS X-axis shift
    3 EXTRAP   I  Extrapolation on/off flag
    4 UNDEF    I  None
    9  X RS    X  value
    10 Y RS    Y  value
    Words 9 through 10 repeat until (-1,-1) occurs
    """
    nvalues = 0
    datas = []

    fmt = endian
    #locut = 0.
    #hicut = 0.
    for table_id in table_ids:
        table = table_dict[table_id]
        nx = len(table.x)
        fmt += b'if 2i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        data = [table_id, table.x1, table.extrap, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  {table_name} data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def _write_table3(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  table_dict: dict[int, Any], table_name: str,
                  key: tuple[int, int, int],
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    Record – TABLEM3(305,3,95)

    Word Name Type Description
    1 ID      I Table identification number
    2 X1     RS X-axis shift
    3 X2     RS X-axis normalization
    4 EXTRAP  I Extrapolation on/off flag
    5 UNDEF(4)  None
    9  X     RS X value
    10 Y     RS Y value
    Words 9 through 10 repeat until (-1,-1) occurs
    """
    nvalues = 0
    datas = []

    fmt = endian
    #locut = 0.
    #hicut = 0.
    for table_id in table_ids:
        table = table_dict[table_id]
        nx = len(table.x)
        fmt += b'i2fi 4i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        data = [table_id, table.x1, table.x2, table.extrap, 0, 0, 0, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  {table_name} data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def _write_table4(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  table_dict: dict[int, Any], table_name: str,
                  key: tuple[int, int, int],
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    TABLED4(1405,14,141)

    Word Name Type Description
    1 ID        I Table identification number
    2 X1       RS X-axis shift
    3 X2       RS X-axis normalization
    4 X3       RS X value when x is less than X3
    5 X4       RS X value when x is greater than X4
    6 UNDEF(3)    None
    9 A        RS
    10 MINUS1   I End of record flag
    """
    nvalues = 0
    datas = []

    fmt = endian
    #locut = 0.
    #hicut = 0.
    for table_id in table_ids:
        table = table_dict[table_id]
        na = len(table.a)
        fmt += b'i 4f 3i ' + b'f' * na + b' i'

        data = [table_id, table.x1, table.x2, table.x3, table.x4, 0, 0, 0] + table.a.tolist() + [-1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  {table_name} data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tables1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    TABLES1(3105, 31, 97)
    Word Name Type Description
    1 ID   I   Table identification number
    2 UNDEF(7) None
    9  X  RS X value
    10 Y  RS Y value
    Words 9 through 10 repeat until (-1,-1) occurs

    TYPE is unique to MSC and is not documented...
    """
    key = (3105, 31, 97)
    nvalues = 0
    datas = []

    fmt = endian
    #locut = 0.
    #hicut = 0.
    for table_id in table_ids:
        table = model.tables[table_id]
        assert table.Type == 1, table.get_stats()
        nx = len(table.x)
        fmt += b'8i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        data = [table_id, 0, 0, 0, 0, 0, 0, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  TABLES1 data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tablest(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """
    Record – TABLEST(1905,19,178)
    Word Name Type Description
    1 ID        I Table identification number
    2 EXTRAP    I Extrapolation on/off flag
    3 UNDEF(6)    None
    9 TI       RS Temperature
    10 TIDI     I TABLES1 Bulk Data entry identification number
    Words 9 through 10 repeat until (-1,-1) occurs

    EXTRAP is unique to NX
    """
    key = (1905, 19, 178)
    nvalues = 0
    datas = []

    fmt = endian
    for table_id in table_ids:
        table = model.tables[table_id]
        nx = len(table.x)
        fmt += b'8i ' + b'fi' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        extrap = table.extrap
        data = [table_id, extrap, 0, 0, 0, 0, 0, 0] + xys + [-1, -1]
        #print(data)
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  TABLEST data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tableh1(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:
    """MSC-specific card"""
    key = (14605, 146, 617)
    nvalues = 0
    datas = []

    fmt = endian
    #locut = 0.
    #hicut = 0.
    for table_id in table_ids:
        table = model.tables[table_id]
        nx = len(table.x)
        fmt += b'8i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)

        data = [table_id, 0, 0, 0, 0, 0, 0, 0] + xys + [-1, -1]
        assert None not in data, data
        #print(data)
        datas.extend(data)
        op2_ascii.write(f'  TABLEH1 data={data}\n')
        nvalues += len(data)

    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    data_bytes = Struct(fmt).pack(*datas)
    op2_file.write(data_bytes)
    return nbytes

def write_tableht(model: BDF | OP2Geom, name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes,
                  nastran_format: str='nx') -> int:  # pragma: no cover
    """(14705, 147, 618)"""
    nbytes = _write_tablet_sh(
        model, name, table_ids, ncards,
        model.tables, 'TABLEHT',
        (14705, 147, 618),
        op2_file, op2_ascii, endian,
        nastran_format=nastran_format)
    return nbytes

def write_gust(model: BDF | OP2Geom, name: str,
               gust_ids: list[GUST], ncards: int,
               op2_file, op2_ascii, endian: bytes,
               nastran_format: str='nx') -> int:
    """
    GUST(1005,10,174) - the marker for Record 1
    """
    key = (1005, 10, 174)
    nfields = 5
    structi = Struct(endian + b'ii3f')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for gust_id in gust_ids:
        gust: GUST = model.gusts[gust_id]
        # (sid, dload, wg, x0, V) = out
        data = [gust.sid, gust.dload, gust.wg, gust.x0, gust.V]
        #print(gust.get_stats())
        assert None not in data, data
        op2_ascii.write(f'  GUST data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes


DIT_MAP = {
    'TABDMP1': write_tabdmp1,
    'TABLED1': write_tabled1,
    'TABLED2': write_tabled2,
    'TABLED3': write_tabled3,
    'TABLED4': write_tabled4,

    'TABLEM1': write_tablem1,
    'TABLEM2': write_tablem2,
    'TABLEM3': write_tablem3,
    'TABLEM4': write_tablem4,

    'TABRNDG': write_tabrndg,
    'TABLES1': write_tables1,
    'TABLEST': write_tablest,

    'TABLEH1': write_tableh1,
    #'TABLEHT': write_tableht,
    'GUST': write_gust,
}
