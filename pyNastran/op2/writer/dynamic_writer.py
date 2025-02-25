"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import Any, BinaryIO, TYPE_CHECKING

from pyNastran.op2.tables.geom.edom import (
    DSCREEN_RTYPE_TO_INT, DRESP_FLAG_TO_RESP_MSC, DRESP_FLAG_TO_RESP_NX)
from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header_nvalues # write_header,

if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.bdf.cards.optimization import DVPREL2, DVMREL2, DTABLE
    #from pyNastran.bdf.cards.aero.static_loads import AEROS # , AESTAT, CSSCHD, DIVERG, TRIM, TRIM2
    #from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, FLUTTER # , FLFACT, MKAERO2
    from pyNastran.op2.op2_geom import OP2Geom, BDF
from pyNastran.utils.numpy_utils import integer_types, float_types

#DRESP_MSC_TO_FLAG = {value: key for key, value in DRESP_FLAG_TO_RESP_MSC.items()}
#DRESP_NX_TO_FLAG = {value: key for key, value in DRESP_FLAG_TO_RESP_NX.items()}

def write_dynamic(op2_file: BinaryIO,
                  op2_ascii, model: BDF | OP2Geom,
                  endian: bytes=b'<',
                  nastran_format: str='nx') -> None:
    """writes the DYNAMIC table"""
    # if not hasattr(model, 'loads'):  # OP2
    #     return
    #card_types = [
        #'DESVAR', 'DCONSTR',
        #'DVMREL2', 'DVPREL2',
        #'DTABLE', 'DSCREEN',
    #]
    card_types = list(DYNAMIC_MAP.keys())
    #print(card_types)

    cards_to_skip = [
        #'DVCREL2',  # disabled b/c no reader
        # 'DOPTPRM',
        # 'DCONADD',
        #
        # 'DLINK',
    ]
    out = defaultdict(list)

    # for rtype, dscreen in sorted(model.dscreen.items()):
    #     out[dscreen.type].append(dscreen)
    # for desvar_id, desvar in sorted(model.desvars.items()):
    #     out[desvar.type].append(desvar_id)
    # for oid, dresp in sorted(model.dresps.items()):
    #     out[dresp.type].append(oid)
    # for oid, dconadd in sorted(model.dconadds.items()):
    #     out[dconadd.type].append(oid)
    # for oid, dconstrs in sorted(model.dconstrs.items()):
    #     for dconstr in dconstrs:
    #         out[dconstr.type].append(dconstr)
    #
    # for dvprel_id, dvprel in sorted(model.dvprels.items()):
    #     out[dvprel.type].append(dvprel_id)
    # for dvmrel_id, dvmrel in sorted(model.dvmrels.items()):
    #     out[dvmrel.type].append(dvmrel_id)
    # for dvcrel_id, dvcrel in sorted(model.dvcrels.items()):
    #     out[dvcrel.type].append(dvcrel_id)
    for delay_id, delay in sorted(model.delays.items()):
        out[delay.type].append(delay_id)
    for dphase_id, dphase in sorted(model.dphases.items()):
        out[dphase.type].append(dphase_id)
    # for dvgrid_id, dvgrid in sorted(model.dvgrids.items()):
    #     out['DVGRID'].append(dvgrid_id)
    # if model.doptprm:
    #     out[model.doptprm.type].append(model.doptprm)
    # if model.dtable:
    #     out[model.dtable.type].append(model.dtable)

    for card_type in list(out):
        if card_type not in card_types:
            del out[card_type]
            model.log.warning(f'removing {card_type} in OP2 writer')
    # other
    if len(out) == 0:
        return

    write_geom_header(b'DYNAMIC', op2_file, op2_ascii, endian=endian)
    itable = -3

    for name, ids in sorted(out.items()):
        model.log.debug('DYNAMIC %s %s' % (name, ids))
        #print('EDOM %s %s' % (name, ids))
        ncards = len(ids)
        assert ncards > 0, ncards
        if name in cards_to_skip:
            model.log.warning('skipping DYNAMIC-%s' % name)
            continue

        #if nmaterials == 0:
            #continue
        try:
            func = DYNAMIC_MAP[name]
        except KeyError:  # pragma: no cover
            raise NotImplementedError(name)

        nbytes = func(model, name, ids, ncards, op2_file, op2_ascii, endian,
                      nastran_format=nastran_format)
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


def write_delay(model: BDF | OP2Geom, name: str,
                delay_ids: list[int], ncards: int,
                op2_file: BinaryIO, op2_ascii, endian: bytes,
                nastran_format: str='nx') -> int:
    """
    DELAY(37,18,183) - Record 3

    1 SID I  Set identification number
    2 P   I  Grid, scalar, or extra point identification number
    3 C   I  Component number
    4 T   RS Time delay

    """
    key = (37, 18, 183)
    structi = Struct(endian + b'3if')

    ncards = 0
    for delay_id in delay_ids:
        delay = model.delays[delay_id]
        ncards += len(delay.nodes)
    nvalues = 4 * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    for delay_id in delay_ids:
        delay = model.delays[delay_id]
        #sid, nodes, components, delays = out
        sid = delay.sid
        for nid, comp, delay in zip(delay.nodes, delay.components, delay.delays):
            data = [sid, nid, comp, delay]
            assert None not in data, data
            op2_ascii.write(f'  DELAY data={data}\n')
            op2_file.write(structi.pack(*data))
    return nbytes

def write_dphase(model: BDF | OP2Geom, name: str,
                dphase_ids: list[int], ncards: int,
                op2_file: BinaryIO, op2_ascii, endian: bytes,
                nastran_format: str='nx') -> int:
    """
    DPHASE(77,19,184) - Record 5

    1 SID I Load set identification number
    2 P   I Grid, scalar, or extra point identification number
    3 C   I Component number
    4 TH  RS Phase lead

    """
    key = (77, 19, 184)
    structi = Struct(endian + b'3if')

    ncards = 0
    for dphase_id in dphase_ids:
        dphase = model.dphases[dphase_id]
        ncards += len(dphase.nodes)
    nvalues = 4 * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    for dphase_id in dphase_ids:
        dphase = model.dphases[dphase_id]
        sid = dphase.sid
        for nid, comp, phase_lead in zip(dphase.nodes, dphase.components, dphase.phase_leads):
            data = [sid, nid, comp, phase_lead]
            assert None not in data, data
            op2_ascii.write(f'  DPHASE data={data}\n')
            op2_file.write(structi.pack(*data))
    return nbytes

DYNAMIC_MAP = {
    'DELAY': write_delay,
    'DPHASE': write_dphase,
}
