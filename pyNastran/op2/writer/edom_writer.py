"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import List, Union, Any, TYPE_CHECKING

from pyNastran.op2.tables.geom.edom import DSCREEN_RTYPE_TO_INT
from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.cards.optimization import DVPREL2, DVMREL2, DTABLE
    #from pyNastran.bdf.cards.aero.static_loads import AEROS # , AESTAT, CSSCHD, DIVERG, TRIM, TRIM2
    #from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, FLUTTER # , FLFACT, MKAERO2
    from pyNastran.op2.op2_geom import OP2Geom, BDF
integer_types = int

def write_edom(op2_file, op2_ascii, model: Union[BDF, OP2Geom], endian: bytes=b'<') -> None:
    """writes the EDOM table"""
    if not hasattr(model, 'loads'):  # OP2
        return
    card_types = [
        'DESVAR', 'DCONSTR',
        'DVMREL2', 'DVPREL2',
        'DTABLE', 'DSCREEN',
    ]

    cards_to_skip = [
        'DVCREL1', 'DVMREL1', 'DVPREL1',
        'DVCREL2',
        'DOPTPRM',
        'DCONADD',

        'DLINK',
        'DVGRID',
    ]
    out = defaultdict(list)

    for rtype, dscreen in sorted(model.dscreen.items()):
        out[dscreen.type].append(dscreen)
    for desvar_id, desvar in sorted(model.desvars.items()):
        out[desvar.type].append(desvar_id)
    for oid, dresp in sorted(model.dresps.items()):
        out[dresp.type].append(oid)
    for oid, dconadd in sorted(model.dconadds.items()):
        out[dconadd.type].append(oid)
    for oid, dconstrs in sorted(model.dconstrs.items()):
        for dconstr in dconstrs:
            out[dconstr.type].append(dconstr)

    for dvprel_id, dvprel in sorted(model.dvprels.items()):
        out[dvprel.type].append(dvprel_id)
    for dvmrel_id, dvmrel in sorted(model.dvmrels.items()):
        out[dvmrel.type].append(dvmrel_id)
    for dvcrel_id, dvcrel in sorted(model.dvcrels.items()):
        out[dvcrel.type].append(dvcrel_id)
    for dlink_id, dlink in sorted(model.dlinks.items()):
        out[dlink.type].append(dlink_id)
    for dvgrid_id, dvgrid in sorted(model.dvgrids.items()):
        out['DVGRID'].append(dvgrid_id)
    if model.doptprm:
        out[model.doptprm.type].append(model.doptprm)
    if model.dtable:
        out[model.dtable.type].append(model.dtable)

    for card_type in list(out):
        if card_type not in card_types:
            del out[card_type]
            model.log.warning(f'removing {card_type} in OP2 writer')
    # other
    if len(out) == 0:
        return

    write_geom_header(b'EDOM', op2_file, op2_ascii, endian=endian)
    itable = -3

    for name, ids in sorted(out.items()):
        model.log.debug('EDOM %s %s' % (name, ids))
        #print('EDOM %s %s' % (name, ids))
        ncards = len(ids)
        assert ncards > 0, ncards
        if name in cards_to_skip:
            model.log.warning('skipping EDOM-%s' % name)
            continue

        #if nmaterials == 0:
            #continue
        try:
            func = EDOM_MAP[name]
        except KeyError:  # pragma: no cover
            raise NotImplementedError(name)

        nbytes = func(model, name, ids, ncards, op2_file, op2_ascii, endian)
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

def _write_dscreen(model: Union[BDF, OP2Geom], name: str,
                   dscreens: List[Any], ncards: int,
                   op2_file, op2_ascii, endian: bytes) -> int:
    """
    DSCREEN(4206, 42, 363)
    Design constraint screening data.

    Word Name Type Description
    1 RTYPE I Response type for which the screening criteria apply
    2 TRS  RS Truncation threshold
    3 NSTR  I Maximum number of constraints to be retained per region per load case

    data = (5, -0.70, 10)
    """
    key = (4206, 42, 363)
    structi = Struct(endian + b'ifi')

    nvalues = 3 * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    for dscreen in dscreens:
        rtype_int = DSCREEN_RTYPE_TO_INT[dscreen.rtype]
        data = [rtype_int, dscreen.trs, dscreen.nstr]
        assert None not in data, data
        op2_ascii.write(f'  DSCREEN data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def _write_desvar(model: Union[BDF, OP2Geom], name: str,
                  desvar_ids: List[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
    """
    (3106, 31, 352)
    NX 2019.2

    Word Name  Type  Description
    1 ID       I     Unique design variable identification number
    2 LABEL(2) CHAR4 User-supplied name for printing purposes
    4 XINIT    RS    Initial value
    5 XLB      RS    Lower bound
    6 XUB      RS    Upper bound
    7 DELXV    RS    Fractional change allowed for the design variable
                     during approximate optimization
    8 DDVAL    I     ID of a DDVAL entry that provides a set of allowable
                     discrete values

    """
    key = (3106, 31, 352)
    structi = Struct(endian + b'i8s ffff i')

    nvalues = 8 * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    for desvar_id in desvar_ids:
        desvar = model.desvars[desvar_id]

        label = desvar.label
        label_bytes = ('%-8s' % label).encode('ascii')

        xinit = desvar.xinit
        xlb = desvar.xlb
        xub = desvar.xub
        delx = desvar.delx
        ddval = desvar.ddval
        if delx is None:
            delx = 0
        if ddval is None:
            ddval = 0

        data = [desvar_id, label_bytes, xinit, xlb, xub, delx, ddval]
        assert None not in data, data
        #print(data)
        op2_ascii.write(f'  DESVAR data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def _write_dconstr(model: Union[BDF, OP2Geom], name: str,
                   dconstrs: List[int], ncards: int,
                   op2_file, op2_ascii, endian: bytes) -> int:
    """
    Record – DCONSTR(4106,41,362)
    Design constraints.
    Word Name Type Description
    1 DCID    I Design constraint set identification number
    2 RID     I DRESPi entry identification number
    3 LALLOW RS Lower bound on the response quantity. Undefined if
                LTID is nonzero
    4 UALLOW RS Upper bound on the response quantity. Undefined if
                UTID is nonzero
    5 LOWFQ  RS Low end of frequency range in Hz
    6 HIGHFQ RS High end of frequency range in Hz
    7 LTID    I Identification number of TABLEDi entry giving lower
                bound on the response quantity as a function of
                frequency or 0 if not specified
    8 UTID    I Identification number of TABLEDi entry giving upper
                bound on the response quantity as a function of
                frequency or 0 if not specified

    data  = (50, 2, 0.0016, 0.0018, 0.0, 1.0e+20, 0, 0)

    """
    key = (4106, 41, 362)
    structi = Struct(endian + b'ii 4f ii')

    nvalues = 8 * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    for dconstr in dconstrs:
        if isinstance(dconstr.uid, integer_types):
            utid = dconstr.uid
            uallow = 0.0
        else:
            utid = 0
            uallow = dconstr.uid

        if isinstance(dconstr.lid, integer_types):
            ltid = dconstr.lid
            lallow = 0.0
        else:
            ltid = 0
            lallow = dconstr.lid

        data = [dconstr.oid, dconstr.dresp_id, lallow, uallow, dconstr.lowfq, dconstr.highfq, ltid, utid]
        assert None not in data, data
        #print(data)
        op2_ascii.write(f'  DCONSTR data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def _write_dvprel2(model: Union[BDF, OP2Geom], name: str,
                   dvprel_ids: List[int], ncards: int,
                   op2_file, op2_ascii, endian: bytes):
    """
    Record – DVPREL2(3406,34,355)
    Design variable to property relation based on a user-supplied equation.

    Word Name Type Description
    1 ID          I Unique identification number
    2 TYPE(2) CHAR4 Name of a property entry
    4 PID         I Property entry identification number
    5 FID      RS/I FID number input. Otherwise, either 0 if property
                    name is input, or frequency (RS) if entry is for
                    frequency dependent property. (See Words 9 and 10)
    6 PMIN       RS Minimum value allowed for this property
    7 PMAX       RS Maximum value allowed for this property
    8 EQID        I DEQATN entry identification number
    9 PNAME1  CHAR4 First word of property name, if any, or blanks if
                    FID number is nonzero in Word 5
    10 PNAME2 CHAR4 Second word of property name, if any. Otherwise,
                    either blanks if FID number is nonzero in Word 5,
                    or frequency (RS) if entry is for frequency
                    dependent property. (See Word 5)
    11 FLAG I DESVAR/DTABLE
    FLAG = 1000 DESVAR
      12 DVIDi I A DESVAR entry identification number
      Word 12 repeats until -1000
    FLAG = 2000 DTABLE
      12 LABLi(2) CHAR4 Label for a constant on the DTABLE entry
      Words 12 and 13 repeat until -2000
    End flag when -1 occurs

    data = (2, PROD, 101, 4, -1.0e+35, 1.0e+20, 2, '', 1000, 2, -1000,
                                                       2000, L1, -2000)
    """
    key = (3406, 34, 355)
    structi = Struct(endian + b'ii 4f ii')

    ndata = 9 * ncards
    nvalues = 11 * ncards
    for dvprel_id in dvprel_ids:
        dvprel = model.dvprels[dvprel_id]
        if len(dvprel.dvids):
            nvalues += len(dvprel.dvids) + 2
            ndata += len(dvprel.dvids) + 2
        if len(dvprel.labels):
            nvalues += 2 * len(dvprel.labels) + 2
            ndata += len(dvprel.labels) + 2

    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    data_all = []
    for dvprel_id in dvprel_ids:
        dvprel = model.dvprels[dvprel_id]

        # TODO: doesn't handle fid = float
        fid = 0
        if isinstance(dvprel.pname_fid, str):
            fmt = b'i 8s 2i 2fi 8s'
            fid = 0
            pname = dvprel.pname_fid
        else:
            fmt = b'i 8s 2i 2fi 8s'
            fid = dvprel.pname_fid
            pname = ''

        p_min = dvprel.p_min
        if p_min is None:
            # NX
            # Minimum value allowed for this property. If FID references a stress
            # recovery location field, then the default value for PMIN is -1.0+35.
            # PMIN must be explicitly set to a negative number for properties that
            # may be less than zero (for example, field ZO on the PCOMP entry).
            # See Remark 10. (Real; Default = 1.E-20)
            #
            # MSC; default=1e-15
            p_min = 1e-20

        prop_type_bytes = ('%-8s' % dvprel.prop_type).encode('ascii')
        pname_bytes = ('%-8s' % pname).encode('ascii')
        data = [dvprel_id, prop_type_bytes, dvprel.pid, fid,
                p_min, dvprel.p_max, dvprel.dequation, pname_bytes]
        fmt += _write_dvxrel2_flag(dvprel, data)

        assert None not in data, data
        structi = Struct(endian + fmt)
        op2_ascii.write(f'  DVPREL2 data={data}\n')
        op2_file.write(structi.pack(*data))
        data_all += data
    #assert len(data_all) == nvalues, f'ndata={len(data_all)} nvalues={nvalues}'
    assert len(data_all) == ndata, f'ndata={len(data_all)} nvalues={ndata}'
    return nbytes

def _write_dvmrel2(model: Union[BDF, OP2Geom], name: str,
                   dvmrel_ids: List[int], ncards: int,
                   op2_file, op2_ascii, endian: bytes):
    """
    Record – DVMREL2(6400,64,432)
    Design variable to material relation based on a user-supplied equation.
    Word Name Type Description
    1 ID            I Unique identification number
    2 TYPE(2)   CHAR4 Name of a material property entry
    4 MID           I Material identification number
    5 FID           I Entry is 0
    6 MPMIN        RS Minimum value allowed for this property
    7 MPMAX        RS Maximum value allowed for this property
    8 EQID          I DEQATN entry identification number
    9 MPNAME(2) CHAR4 Name of material property
    11 FLAG         I DESVAR/DTABLE
    FLAG = 1000 DESVAR
      12 DVIDi I A DESVAR entry identification number
      Word 12 repeats until -1000
    FLAG = 2000 DTABLE
      12 LABLi(2) CHAR4 Label for a constant on the DTABLE entry
      Words 12 and 13 repeat until -2000
    End flag when -1 occurs

    """
    key = (6400, 64, 432)
    structi = Struct(endian + b'ii 4f ii')

    ndata = 9 * ncards
    nvalues = 11 * ncards
    for dvmrel_id in dvmrel_ids:
        dvmrel = model.dvmrels[dvmrel_id]
        if len(dvmrel.dvids):
            nvalues += len(dvmrel.dvids) + 2
            ndata += len(dvmrel.dvids) + 2
        if len(dvmrel.labels):
            nvalues += 2 * len(dvmrel.labels) + 2
            ndata += len(dvmrel.labels) + 2

    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    data_all = []
    for dvmrel_id in dvmrel_ids:
        dvmrel = model.dvmrels[dvmrel_id]
        #print(dvmrel.get_stats())
        # TODO: doesn't handle fid = float
        #fid = 0
        #if isinstance(dvprel.pname_fid, str):
        fmt = b'i 8s 2i 2fi 8s'
        fid = 0
        #else:
            #fmt = 'i 8s 2i 2fi ii'
            #fid = dvprel.pname_fid
            #pname = ''

        mp_min = dvmrel.mp_min
        if mp_min is None:
            # TODO: not supported
            #
            # Minimum value allowed for this property. If MPNAME references a
            # material property that can only be positive, then the default value for
            # MPMIN is 1.0E-15. Otherwise, it is -1.0E35. (Real)
            mp_min = 1e-20

        mat_type_bytes = ('%-8s' % dvmrel.mat_type).encode('ascii')
        mp_name_bytes = ('%-8s' % dvmrel.mp_name).encode('ascii')
        data = [dvmrel_id, mat_type_bytes, dvmrel.mid, fid,
                mp_min, dvmrel.mp_max, dvmrel.dequation, mp_name_bytes]
        fmt += _write_dvxrel2_flag(dvmrel, data)

        assert None not in data, data
        structi = Struct(endian + fmt)
        op2_ascii.write(f'  DVMREL2 data={data}\n')
        op2_file.write(structi.pack(*data))
        data_all += data
    #assert len(data_all) == nvalues, f'ndata={len(data_all)} nvalues={nvalues}'
    assert len(data_all) == ndata, f'ndata={len(data_all)} nvalues={ndata}'
    return nbytes

def _write_dvxrel2_flag(dvxrel2: Union[DVPREL2, DVMREL2], data: List[Any]) -> bytes:
    """writes the DVxREL2 flag table"""
    fmt = b''
    ndesvars = len(dvxrel2.dvids)
    if ndesvars:
        fmt += b' %ii' % (2 + ndesvars)
        data += [1000] + dvxrel2.dvids + [-1000]

    nlabels = len(dvxrel2.labels)
    if nlabels:
        fmt += b' i%isi' % (nlabels*8)
        labels_bytes = [('%-8s' % label).encode('ascii') for label in dvxrel2.labels]
        data += [2000] + labels_bytes + [-2000]
    fmt += b' i'
    data.append(-1)
    return fmt

def _write_dtable(model: Union[BDF, OP2Geom], name: str,
                  dtables: List[DTABLE], ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
    """
    Record – DTABLE(3706,37,358)
    Table constants.
    Word Name Type Description
    1 LABLi(2) CHAR4 Label for the constant
    3 VALUi       RS Value of the constant
    Words 1 thru 3 repeat until -1 occurs

    """
    assert ncards == 1, ncards
    #if self.size == 4:
    #else:
        #aaa

    dtable = model.dtable
    nkeys = len(dtable.default_values)
    fmt = endian + b'8s f'*nkeys + b'i'
    structi = Struct(fmt)

    key = (3706, 37, 358)
    nvalues = 3 * nkeys + 1
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    data = []
    for key, value in dtable.default_values.items():
        key_bytes = ('%-8s' % key).encode('ascii')
        data += [key_bytes, value]
    data.append(-1)

    assert None not in data, data
    op2_ascii.write(f'  DTABLE data={data}\n')
    out = structi.pack(*data)
    op2_file.write(out)

    return nbytes


EDOM_MAP = {
    #'DVPREL1': _write_dvprel1, 'DVMREL1': _write_dvmrel1, 'DVCREL1': _write_dvcrel1,
    'DVPREL2': _write_dvprel2, 'DVMREL2': _write_dvmrel2, # 'DVCREL2': _write_dvcrel2,

    'DESVAR': _write_desvar,
    #'DOPTPRM': _write_doptprm,
    'DTABLE': _write_dtable,
    'DCONSTR': _write_dconstr,
    #'DCONADD': _write_dconadd,

    #'DLINK': _write_dlink,
    #'DESVAR': _write_desvar,
    #'DVGRID': _write_dvgrid,
    'DSCREEN': _write_dscreen,

    #(3206, 32, 353) : ['DLINK', self._read_fake],
    #(3306, 33, 354) : ['DVPREL1', self._read_dvprel1],
    #(3806, 38, 359) : ['DRESP1', self._read_fake],
    #(3906, 39, 360) : ['DRESP2', self._read_fake],
    #(4206, 42, 363) : ['DSCREEN', self._read_fake],
    #(4306, 43, 364) : ['DOPTPRM', self._read_doptprm],
    #(4406, 44, 372) : ['DVGRID', self._read_dvgrid],
}
