from __future__ import absolute_import
from struct import pack, Struct
from collections import defaultdict

from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from .geom1_writer import write_geom_header, close_geom_table

def write_geom4(op2, op2_ascii, obj, endian=b'<', nastran_format='nx'):
    if not hasattr(obj, 'rigid_elements'):
        return
    loads_by_type = defaultdict(list)
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
    for name, usets in obj.usets.items():
        for uset in usets:
            loads_by_type[uset.type].append(uset)
    for omit in obj.omits:
        assert not isinstance(omit, int), obj.omits
        loads_by_type[omit.type].append(omit)

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
    skip_cards = [

        # rigid elements
        'RBE1', 'RBE2', 'RBE3', 'RBAR', 'RROD', 'RSSCON',

        # sets
        'ASET', 'ASET1', 'BSET', 'BSET1', 'CSET', 'CSET1',
        'QSET', 'QSET1', 'OMIT1', 'USET', 'USET1',
    ]
    supported_cards = [
        'MPC', 'SPC', 'SPC1', 'SPCADD', 'MPCADD',
    ]
    is_constraints = False
    for card_type in sorted(loads_by_type.keys()):
        if card_type in skip_cards:
            obj.log.warning('skipping GEOM4-%s' % card_type)
            continue
        if card_type in supported_cards:
            is_constraints = True
            continue
            #break
        obj.log.warning('skipping GEOM4-%s' % card_type)
    #else:
        #return

    write_geom_header(b'GEOM4', op2, op2_ascii)
    itable = -3
    for card_type, cards in sorted(loads_by_type.items()):
        #if card_type in ['SPCD']: # not a GEOM3 load
            #continue
        if card_type in skip_cards:
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

def write_card(op2, op2_ascii, card_type, cards, endian, nastran_format='nx'):
    ncards = len(cards)
    if card_type == 'MPC':
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
        #while i < nfields:
            #sid, grid, comp = ints[i:i+3]
            #coeff = floats[i+3]
            #mpc_data = [sid, grid, comp, coeff]
            #nodes = [grid]
            #components = [comp]
            #coefficients = [coeff]

            #intsi = ints[i+4:i+7]
            #assert len(intsi) == 3, intsi
            #while intsi != (-1, -1, -1):
                #gridi, compi, coeffi = ints[i+4], ints[i+5], floats[i+6]
                #mpc_data.extend([gridi, compi, coeffi])
                #nodes.append(gridi)
                #components.append(compi)
                #coefficients.append(coeffi)
                #i += 3
                #intsi = ints[i+4:i+7]
            #mpc_data.extend([-1, -1, -1])
            #i += 7 # 3 + 4 from (-1,-1,-1) and (sid,grid,comp,coeff)
            #if self.is_debug_file:
                #self.binary_debug.write('  MPC=%s\n' % str(mpc_data))
            #mpci = MPC.add_op2_data((sid, nodes, components, coefficients))
            #self._add_constraint_mpc_object(mpci)
        #raise NotImplementedError('MPC')

        #for load in cards:
            #data = [load.sid, load.node_id, load.Cid(), load.mag] + list(load.xyz)
            #op2_ascii.write('  MPC data=%s\n' % str(data))
            #op2.write(spack.pack(*data))

    elif card_type == 'SPC1':
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

        nfields = 0
        fields = []
        data = defaultdict(list)
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
        del fields

    elif card_type in ['SPCADD', 'MPCADD']:
        if card_type == 'SPCADD':
            key = (5491, 59, 13)
        else:
            # MPCADD
            key = (4891, 60, 83)
        #SPCADD(5491,59,13)
        #MPCADD(4891,60,83)
        #[2  1 10 -1]
        #[3  1 -1]
        data = []
        for spcadd in cards:
            #print(spcadd.get_stats())
            datai = [spcadd.conid] + spcadd.ids + [-1]
            op2_ascii.write('  %s data=%s\n' % (card_type, str(datai)))
            data += datai
        nfields = len(data)
        nbytes = write_header_nvalues(card_type, nfields, key, op2, op2_ascii)
        spack = Struct(endian + b'%ii' % nfields)
        op2.write(spack.pack(*data))
    elif card_type == 'SPC':
        key = (5501, 55, 16)
        #nastran_format = 'msc'
        if nastran_format == 'msc':
            nfields = 5
            nbytes = write_header(card_type, nfields, ncards, key, op2, op2_ascii)
            # MSC
            # SPC(5501,55,16) - Record 44
            #
            # 1 SID   I    Set identification number
            # 2 ID    I    Grid or scalar point identification number
            # 3 C     I    Component numbers
            # 4 UNDEF none Not used
            # 5 D     RX   Enforced displacement
            data = []
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
            data = []
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

def write_header_nvalues(name, nvalues, key, op2, op2_ascii):
    """a more precise version of write_header for when card lengths can vary"""
    nvalues += 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes

def write_header(name, nfields, ncards, key, op2, op2_ascii):
    """writes the op2 card header given the number of cards and the fields per card"""
    nvalues = nfields * ncards
    nbytes = write_header_nvalues(name, nvalues, key, op2, op2_ascii)
    return nbytes
