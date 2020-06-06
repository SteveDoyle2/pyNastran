from struct import pack, Struct
from collections import defaultdict
from typing import List, Dict, Tuple, Union, Any

from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.op2.errors import SixtyFourBitError
from .geom1_writer import write_geom_header, close_geom_table

def write_geom4(op2, op2_ascii, obj, endian: bytes=b'<', nastran_format: str='nx') -> None:
    if not hasattr(obj, 'rigid_elements'):
        return
    loads_by_type = defaultdict(list)  # type: Dict[str, Any]
    for unused_id, rigid_element in obj.rigid_elements.items():
        loads_by_type[rigid_element.type].append(rigid_element)
    for aset in obj.asets:
        assert not isinstance(aset, int), obj.asets
        # ASET
        loads_by_type[aset.type].append(aset)
    for bset in obj.bsets:
        assert not isinstance(bset, int), obj.bsets
        loads_by_type[bset.type].append(bset)
    for cset in obj.csets:
        assert not isinstance(cset, int), obj.csets
        loads_by_type[cset.type].append(cset)
    for qset in obj.qsets:
        assert not isinstance(qset, int), obj.qsets
        loads_by_type[qset.type].append(qset)
    for unused_name, usets in obj.usets.items():
        for uset in usets:
            loads_by_type[uset.type].append(uset)
    for omit in obj.omits:
        assert not isinstance(omit, int), obj.omits
        loads_by_type[omit.type].append(omit)

    for suport in obj.suport:
        loads_by_type[suport.type].append(suport)
    for unused_idi, suport in obj.suport1:
        loads_by_type[suport.type].append(suport)

    for unused_id, spcs in obj.spcs.items():
        for spc in spcs:
            loads_by_type[spc.type].append(spc)
    for unused_id, mpcs in obj.mpcs.items():
        for mpc in mpcs:
            loads_by_type[mpc.type].append(mpc)

    for unused_id, spcadds in obj.spcadds.items():
        for spcadd in spcadds:
            loads_by_type[spcadd.type].append(spcadd)
    for unused_id, mpcadds in obj.mpcadds.items():
        for mpcadd in mpcadds:
            loads_by_type[mpcadd.type].append(mpcadd)

    #for unused_load_id, load in obj.tempds.items():
        #loads_by_type[load.type].append(load)

    # return if no supported cards are found
    skip_cards = {
        #'SUPORT', 'SUPORT1', # suport

        # spcs
        'GMSPC',

        # rigid elements
        'RROD', 'RSSCON',

        # sets
        'CSET1',
        'USET', 'USET1',
    }
    # not defined in DMAP
    not_defined_cards = {'RBAR1'}
    supported_cards = {
        'SUPORT', 'SUPORT1', # suport

        # sets
        'ASET', 'BSET', 'CSET', 'QSET', 'OMIT',
        'ASET1', 'BSET1', 'QSET1', 'OMIT1',
        # rigid
        'RBAR', 'RBE1', 'RBE2', 'RBE3',
        # constraints
        'MPC', 'SPC', 'SPC1', 'SPCADD', 'MPCADD',
    }

    is_constraints = False
    for card_type in sorted(loads_by_type.keys()):
        if card_type in skip_cards:
            obj.log.warning('skipping GEOM4-%s' % card_type)
            continue
        if card_type in not_defined_cards:
            continue
        if card_type in supported_cards:
            is_constraints = True
            continue
            #break
        obj.log.warning('skipping GEOM4-%s' % card_type)
    #else:
        #return

    if not is_constraints:
        return
    write_geom_header(b'GEOM4', op2, op2_ascii)
    itable = -3
    for card_type, cards in sorted(loads_by_type.items()):
        #if card_type in ['SPCD']: # not a GEOM3 load
            #continue
        if card_type in skip_cards or card_type in not_defined_cards:
            continue

        try:
            nbytes = write_card(op2, op2_ascii, card_type, cards, endian,
                                nastran_format=nastran_format)
        except:  # pragma: no cover
            obj.log.error('failed GEOM4-%s' % card_type)
            raise
        op2.write(pack('i', nbytes))

        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)

    #-------------------------------------

def write_card(op2, op2_ascii, card_type: str, cards, endian: bytes,
               nastran_format: str='nx') -> int:
    ncards = len(cards)
    if card_type in ['ASET1', 'BSET1', 'CSET1', 'QSET1', 'OMIT1']:
        nbytes = _write_xset1(card_type, cards, ncards, op2, op2_ascii,
                              endian)
    elif card_type in ['ASET', 'BSET', 'CSET', 'OMIT', 'QSET']:
        nbytes = _write_xset(card_type, cards, ncards, op2, op2_ascii,
                             endian)
    elif card_type == 'SUPORT':
        key = (5601, 56, 14)
        data = []
        fmt = endian
        for suport in cards:
            datai = []
            nnodes = len(suport.Cs)
            for nid, ci in zip(suport.node_ids, suport.Cs):
                assert isinstance(nid, int), suport.get_stats()
                assert isinstance(ci, str), suport.get_stats()
                datai.extend([nid, int(ci)])
            fmt += b'%ii' % (nnodes * 2)
            data.extend(datai)
            op2_ascii.write('  SUPORT data=%s\n' % str(datai))
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        op2.write(pack(fmt, *data))
        del data, fmt

    elif card_type == 'SUPORT1':
        key = (10100, 101, 472)
        data = []
        fmt = endian
        for suport1 in cards:
            suport1i = [suport1.conid]
            nnodes = len(suport1.Cs)
            for nid, ci in zip(suport1.node_ids, suport1.Cs):
                assert isinstance(nid, int), suport1.get_stats()
                assert isinstance(ci, int), suport1.get_stats()
                suport1i.extend([nid, ci])
            suport1i.append(-1)
            op2_ascii.write('  SUPORT1 data=%s\n' % str(suport1i))
            fmt += b'%ii' % (2 * nnodes + 2)
            data.extend(suport1i)
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        op2.write(pack(fmt, *data))
        del data, fmt

    elif card_type == 'MPC':
        key = (4901, 49, 17)
        data = []
        fmt = endian
        for mpc in cards:
            datai = [mpc.conid, ]
            fmt += b'i' + b'ifi' * len(mpc.coefficients) + b'3i'
            for nid, coeff, comp in zip(mpc.node_ids, mpc.coefficients, mpc.components):
                datai += [nid, coeff, int(comp)]
            datai += [-1, -1, -1]
            op2_ascii.write('  MPC data=%s\n' % str(datai))
            data += datai

        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        op2.write(pack(fmt, *data))
        del data, fmt

    elif card_type == 'RBE1':
        nbytes = _write_rbe1(card_type, cards, ncards, op2, op2_ascii,
                             endian)
    elif card_type == 'RBE2':
        nbytes = _write_rbe2(card_type, cards, ncards, op2, op2_ascii,
                             endian)
    elif card_type == 'RBE3':
        nbytes = _write_rbe3(card_type, cards, ncards, op2, op2_ascii,
                             endian)

    elif card_type == 'RBAR':
        nbytes = _write_rbar(card_type, cards, ncards, op2, op2_ascii,
                             endian, nastran_format=nastran_format)

    elif card_type == 'SPC1':
        nbytes = _write_spc1(card_type, cards, ncards, op2, op2_ascii,
                             endian)

    elif card_type in ['SPCADD', 'MPCADD']:
        if card_type == 'SPCADD':
            key = (5491, 59, 13)
        elif card_type == 'MPCADD':
            key = (4891, 60, 83)
        else:  # pragma: no cover
            raise NotImplementedError(card_type)
        #SPCADD(5491,59,13)
        #MPCADD(4891,60,83)
        #[2  1 10 -1]
        #[3  1 -1]
        data = []
        for spcadd in cards:
            datai = [spcadd.conid] + spcadd.ids + [-1]
            op2_ascii.write('  %s data=%s\n' % (card_type, str(datai)))
            data += datai
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        spack = Struct(endian + b'%ii' % nfields)
        op2.write(spack.pack(*data))
    elif card_type == 'SPC':
        nbytes = _write_spc(card_type, cards, ncards, op2, op2_ascii, endian,
                            nastran_format=nastran_format)
    #elif card_type == 'TEMPD':
        #key = (5641, 65, 98)
        #nfields = 6
        #spack = Struct(endian + b'if')
        #nbytes = write_header(card_type, nfields, ncards, key, op2, op2_ascii)
        #for load in cards:
            #print(load.get_stats())
            ##sid, T = data
            #data = [load.sid, load.temperature]
            #op2_ascii.write('  TEMPD data=%s\n' % str(data))
            #op2.write(spack.pack(*data))
    else:  # pragma: no cover
        card0 = cards[0]
        raise NotImplementedError(card0)
    return nbytes

def write_header_nvalues(name: str, nvalues: int, key: Tuple[int, int, int], op2, op2_ascii):
    """a more precise version of write_header for when card lengths can vary"""
    nvalues += 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes

def write_header(name: str, nfields: int, ncards: int, key: Tuple[int, int, int],
                 op2, op2_ascii) -> int:
    """writes the op2 card header given the number of cards and the fields per card"""
    nvalues = nfields * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2, op2_ascii)
    return nbytes

def _write_spc(card_type: str, cards, ncards: int, op2, op2_ascii,
               endian: bytes, nastran_format: str='nx') -> int:
    """writes an SPC"""
    key = (5501, 55, 16)
    #nastran_format = 'msc'
    max_spc_id = max([spc.conid for spc in cards])
    max_nid = max([max(spc.node_ids) for spc in cards])
    if max_spc_id > 99999999:
        raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max spc_id={max_spc_id}')
    if max_nid > 99999999:
        raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max SPC nid={max_nid}')

    data = []  # type: List[Union[int, float]]
    if nastran_format == 'msc':
        # MSC
        # SPC(5501,55,16) - Record 44
        #
        # 1 SID   I    Set identification number
        # 2 ID    I    Grid or scalar point identification number
        # 3 C     I    Component numbers
        # 4 UNDEF none Not used
        # 5 D     RX   Enforced displacement
        nfields = 5
        nbytes = write_header(card_type, nfields, ncards, key, op2, op2_ascii)
        for spc in cards:
            node_ids = spc.node_ids
            for nid, comp, enforcedi in zip(node_ids, spc.components, spc.enforced):
                datai = [spc.conid, nid, int(comp), 0, enforcedi]
            op2_ascii.write('  SPC data=%s\n' % str(datai))
            data += datai
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        fmt = endian + b'4if' * (nfields // 5)
        op2.write(pack(fmt, *data))

    elif nastran_format == 'nx':
        # NX
        # SPC(5501,55,16) - Record 44
        #
        # 1 SID I  Set identification number
        # 2 ID  I  Grid or scalar point identification number
        # 3 C   I  Component numbers
        # 4 D   RS Enforced displacement
        for spc in cards:
            node_ids = spc.node_ids
            #assert len(node_ids) == 1, spc.get_stats()
            datai = []
            for nid, comp, enforcedi in zip(node_ids, spc.components, spc.enforced):
                datai = [spc.conid, nid, int(comp), enforcedi]
                op2_ascii.write('  SPC data=%s\n' % str(datai))

            data += datai
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        fmt = endian + b'3if' * (nfields // 4)
        op2.write(pack(fmt, *data))
    else:  # pragma: no cover
        raise RuntimeError(f'nastran_format={nastran_format} not msc, nx')
    return nbytes

def _write_rbe1(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                endian: bytes) -> int:
    """
    RBE1(6801,68,294) - Record 23

    MSC/NX
    Word Name Type Description
    1 EID I Element identification number
    2 GN  I Grid point identification number for independent degrees-of-freedom
    3 CN  I Component numbers of independent degrees-of-freedom
    Words 2 through 3 repeat until (-2,-2) occurs

    4 GM  I Grid point identification number for dependent degrees-of-freedom
    5 CM  I Component numbers of dependent degreesof-freedom
    Words 4 through 5 repeat until (-1,-1) occurs

    6 ALPHA RS Thermal expansion coefficient
    7 UNDEF none Not used

    """
    # TODO: neither reader or writer considers alpha; no current examples
    key = (6801, 68, 294)
    fields = []  # type: List[int]
    fmt = endian
    for rbe1 in cards:
        fieldsi = [rbe1.eid]
        ngn = len(rbe1.Gni)
        ngm = len(rbe1.Gmi)
        #fmt += b'i %if 2i %if 2i fi' % (ngn * 2, ngm * 2)
        fmt += b'i %if 2i %if 2i' % (ngn * 2, ngm * 2)
        for gn, cn in zip(rbe1.Gni, rbe1.Cni):
            fieldsi += [gn, int(cn)]
        fieldsi += [-2, -2]
        for gm, cm in zip(rbe1.Gmi, rbe1.Cmi):
            fieldsi += [gm, int(cm)]
        #fieldsi += [-1, -1, rbe1.alpha, 0]
        fieldsi += [-1, -1]
        fields += fieldsi
    nfields = len(fields)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *fields))
    del fields, fmt
    return nbytes

def _write_rbe2(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                endian: bytes) -> int:
    """
    RBE2(6901,69,295) - Record 24

    Word Name Type Description
    1 EID I Element identification number
    2  GN I Grid point identification number for independent degrees-of-freedom
    3  CM I Component numbers of dependent degrees of-freedom
    4  GM I Grid point identification number for dependent degrees-of-freedom
    Word 4 repeats until End of Record
    5 ALPHA RS Thermal expansion coefficient

    """
    is_alpha = False
    for rbe2 in cards:
        if rbe2.alpha != 0.:
            is_alpha = True
            break

    key = (6901, 69, 295)
    fields = []  # type: List[Union[int, float]]
    fmt = endian
    if is_alpha:
        for rbe2 in cards:
            ngm = len(rbe2.Gmi)
            fmt += b'3i %ii' % ngm
            fields += [rbe2.eid, rbe2.gn, int(rbe2.cm)] + rbe2.Gmi
    else:
        for rbe2 in cards:
            ngm = len(rbe2.Gmi)
            fmt += b'3i %ii f' % ngm
            fields += [rbe2.eid, rbe2.gn, int(rbe2.cm)] + rbe2.Gmi + [rbe2.alpha]

    nfields = len(fields)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *fields))
    return nbytes

def _write_rbe3(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                endian: bytes) -> int:
    """
    1 EID   I Element identification number
    2 REFG  I Reference grid point identification number
    3 REFC  I Component numbers at the reference grid point

    4 WT1  RS Weighting factor for components of motion at G
    5 C     I Component numbers
    6 G     I Grid point identification number
    Word 6 repeats until -1 occurs
    Words 4 through 6 repeat until -2 occurs

    7 GM    I Grid point identification number for dependent degrees-of-freedom
    8 CM    I Component numbers of dependent degrees-of-freedom

    Words 7 through 8 repeat until -3 occurs

    """
    # TODO: alpha is not supported in the writer
    key = (7101, 71, 187)
    fields = []  # type: List[Union[int, float]]
    fmt = endian
    for rbe3 in cards:
        fieldsi = [rbe3.eid, rbe3.refgrid, int(rbe3.refc)]
        fmti = b'3i'
        for weight, cg, group in rbe3.wt_cg_groups:
            fieldsi += [weight, int(cg)] + group + [-1]
            ngroup = len(group)
            fmti += b'fi %ii i' % ngroup
        fmti += b'i'
        fieldsi.append(-2)
        ngmi = len(rbe3.Gmi)
        if ngmi:
            ncmi = len(rbe3.Cmi)
            assert ncmi == ngmi, rbe3.get_stats()
            fmti += b'2i' * ncmi
            for cmi, gmi in zip(rbe3.Cmi, rbe3.Gmi):
                fieldsi += [cmi, gmi]
            #op2.log.warning('UM is not implemented')
        fmti += b'i'
        fieldsi += [-3]

        fmt += fmti
        fields += fieldsi
    nfields = len(fields)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *fields))
    del fields, fmt
    return nbytes

def _write_rbar(card_type: str, cards, ncards: int, op2, op2_ascii,
                endian: bytes, nastran_format: str='nx') -> int:
    """writes an RBAR"""
    # MSC
    key = (6601, 66, 292)
    fields = []  # type: List[Union[int, float]]
    if nastran_format == 'msc':
        fmt = endian + b'7if' * ncards
        for rbar in cards:
            cna = int(rbar.cna)
            cma = int(rbar.cma)
            cnb = int(rbar.cnb)# if rbar.cnb != '' else 0
            cmb = int(rbar.cmb)
            fields += [
                rbar.eid, rbar.ga, rbar.gb,
                cna, cnb,
                cma, cmb,
                rbar.alpha]
    elif nastran_format == 'nx':
        fmt = endian + b'7i' * ncards
        for rbar in cards:
            cna = int(rbar.cna)
            cma = int(rbar.cma) if rbar.cma != '' else 0
            cnb = int(rbar.cnb) if rbar.cnb != '' else 0
            cmb = int(rbar.cmb)
            fields += [
                rbar.eid, rbar.ga, rbar.gb,
                cna, cnb, cma, cmb]
    else:  # pragma: no cover
        raise NotImplementedError(nastran_format)
    nfields = len(fields)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *fields))
    del fields, fmt
    return nbytes

def _write_spc1(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                endian: bytes) -> int:
    key = (5481, 58, 12)
    #sid, components = out[:2]
    #thru_flag = out[2]
    #if thru_flag == 0:  # repeat 4 to end
        #nids = out[3:].tolist()
        #thru_check = False
    #elif thru_flag == 1:
        #n1 = out[3]
        #n2 = out[4]
        #nids = list(range(n1, n2+1))
        #thru_check = True
    #else:
        #raise NotImplementedError('SPC1; thru_flag=%s' % thru_flag)

    max_spc_id = max([spc.conid for spc in cards])
    try:
        nids = [max(spc.node_ids) for spc in cards]
    except ValueError:
        for spc in cards:
            try:
                max(spc.node_ids)
            except ValueError:
                print(spc)
                raise
    max_nid = max(nids)

    if max_spc_id > 99999999:
        raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max spc_id={max_spc_id}')
    if max_nid > 99999999:
        raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max SPC nid={max_nid}')

    nfields = 0
    fields = []  # type: List[int]
    data = defaultdict(list)  # type: Dict[Tuple[int, int], List[int]]
    for spc in cards:
        data_key = (spc.conid, int(spc.components))
        ids = spc.node_ids
        data[data_key] += ids

    for data_key, nodes in data.items():
        conid, components = data_key
        singles, doubles = collapse_thru_packs(nodes)
        nsingles = len(singles)
        ndoubles = len(doubles)
        if nsingles:
            nfields += 4 + nsingles # conid, components, thru_flag, singles, -1
            fields += [conid, components, 0, ] + singles
            fields.append(-1)
        if ndoubles:
            for double in doubles:
                fields += [conid, components, 1, double[0], double[2], -1] # + doubles
            nfields += 6 * ndoubles
        assert len(fields) == nfields

    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(endian + b'%ii' % nfields, *fields))
    return nbytes

def _write_xset(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                endian: bytes) -> int:
    """
    Word Name Type Description
    1 ID I Grid or scalar point identification number
    2  C I Component numbers

    """
    if card_type == 'ASET':
        key = (5561, 76, 215)
    elif card_type == 'BSET':
        key = (110, 1, 311)
    elif card_type == 'CSET':
        key = (310, 3, 313)
    elif card_type == 'QSET':
        key = (510, 5, 315)
    elif card_type == 'OMIT':
        key = (5001, 50, 15)
    else:  # pragma: no cover
        raise NotImplementedError(card_type)

    data = []  # type: List[int]
    fmt = endian
    for set_obj in cards:
        nnodes = len(set_obj.components)
        fmt += b'%ii' % (nnodes * 2)
        for nid, comp in zip(set_obj.node_ids, set_obj.components):
            data += [nid, int(comp)]
    nfields = len(data)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *data))
    del data, fmt
    return nbytes

def _write_xset1(card_type: str, cards, unused_ncards: int, op2, op2_ascii,
                 endian: bytes) -> int:
    """
    Word Name Type Description
    1 C        I Component numbers
    2 THRUFLAG I Thru range flag
    THRUFLAG=0 No
      3 ID I Grid or scalar point identification number
      Word 3 repeats until End of Record
    THRUFLAG=1 Yes
      3 ID1 I First grid or scalar point identification number
      4 ID2 I Second grid or scalar point identification number
    End THRUFLAG

    """
    if card_type == 'ASET1':
        key = (5571, 77, 216)
    elif card_type == 'BSET1':
        key = (410, 4, 314)
    elif card_type == 'CSET1':
        key = (210, 2, 312)
    elif card_type == 'QSET1':
        key = (610, 6, 316)
    elif card_type == 'OMIT1':
        key = (4951, 63, 92)
    else:  # pragma: no cover
        raise NotImplementedError(card_type)

    data = []  # type: List[int]
    fmt = endian
    for set_obj in cards:
        nodes = set_obj.node_ids
        singles, doubles = collapse_thru_packs(nodes)
        nsingles = len(singles)
        ndoubles = len(doubles)
        if nsingles == 0 and ndoubles == 1:
            min_node, unused_thru, max_node = doubles[0]
            data += [int(set_obj.components), 1, min_node, max_node]
            fmt += b'4i'
        else:
            data += [int(set_obj.components), 0, ] + nodes
            nnodes = len(nodes)
            fmt += b'%ii' % (nnodes + 2)

    nfields = len(data)
    nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)

    op2.write(pack(fmt, *data))
    del data, fmt
    return nbytes
