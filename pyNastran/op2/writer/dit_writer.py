"""writes the MPT/MPTS table"""
from __future__ import annotations
from copy import deepcopy
from collections import defaultdict
from struct import pack, Struct
from typing import Union, TYPE_CHECKING

import numpy as np

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues
from .edt_writer import remove_unsupported_cards

from pyNastran.utils.numpy_utils import integer_types
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    #from pyNastran.bdf.cards.aero.aero import CAERO1, CAERO2, PAERO1, PAERO2, AESURF, AESURFS
    #from pyNastran.bdf.cards.aero.static_loads import AEROS, AESTAT # , CSSCHD, DIVERG, TRIM, TRIM2
    from pyNastran.bdf.cards.aero.dynamic_loads import GUST
    from pyNastran.op2.op2_geom import OP2Geom, BDF

def write_dit(op2_file, op2_ascii, model: Union[BDF, OP2Geom],
              endian: bytes=b'<', nastran_format: str='nx') -> None:
    """writes the DIT/DITS table"""
    if not hasattr(model, 'loads'):  # OP2
        return
    card_types = [
        'GUST',
        # 'TABDMP1', # op2 is bugged?
    ]

    cards_to_skip = [
        'TABDMP1', # op2 is bugged?
    ]
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
    #-------------------------------------

#def remove_unsupported_cards(card_dict: dict[str, Any],
                             #card_types: list[str],
                             #log: SimpleLogger):

    #for card_type in list(card_dict):
        #if card_type not in card_types:
            #del card_dict[card_type]
            #log.warning(f"removing {card_type} in OP2 writer because it's unsupported")

def write_tabdmp1(model: Union[BDF, OP2Geom], name: str,
                  table_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
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
        fmt += b'i ' + b'2f' * nx + b' 2i'
        xys = []
        for x, y in zip(table.x, table.y):
            xys.append(x)
            xys.append(y)
        data = [table_id] + xys + [-1, -1]
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

def write_gust(model: Union[BDF, OP2Geom], name: str,
               gust_ids: list[GUST], ncards: int,
               op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    GUST(1005,10,174) - the marker for Record 1
    """
    key = (1005, 10, 174)
    nfields = 5
    structi = Struct(endian + b'ii3f')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for gust_id in gust_ids:
        gust = model.gusts[gust_id]
        # (sid, dload, wg, x0, V) = out
        data = [gust.sid, gust.dload, gust.wg, gust.x0, gust.V]
        #flutter = model.loads[flutter_id]  # type: FLUTTER
        #print(flutter.get_stats())
        assert None not in data, data
        op2_ascii.write(f'  GUST data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes


DIT_MAP = {
    'TABDMP1': write_tabdmp1,
    'GUST': write_gust,
}
