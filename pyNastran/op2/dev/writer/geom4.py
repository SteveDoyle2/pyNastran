from __future__ import absolute_import
from struct import pack, Struct
from collections import defaultdict

from .geom1 import write_geom_header, close_geom_table

def write_geom4(op2, op2_ascii, obj, endian=b'<'):
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
        'MPC', 'SPC1', 'SPCADD', 'SPC',

        # rigid elements
        'RBE1', 'RBE2', 'RBE3', 'RBAR', 'RROD', 'RSSCON',

        # sets
        'ASET', 'ASET1', 'BSET', 'BSET1', 'CSET', 'CSET1',
        'QSET', 'QSET1', 'OMIT1', 'USET', 'USET1',
    ]
    supported_cards = [
    ]
    is_constraints = False
    for card_type in sorted(loads_by_type.keys()):
        if card_type in skip_cards:
            obj.log.warning('skipping GEOM4-%s' % card_type)
            continue
        if card_type in supported_cards:
            is_constraints = True
            break
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
            nbytes = write_card(op2, op2_ascii, card_type, cards, endian)
        except:
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

def write_card(op2, op2_ascii, card_type, cards, endian):
    ncards = len(cards)
    if card_type == 'MPC':
        key = (4901, 49, 17)
        nfields = 7
        #print(cards)
        #spack = Struct(endian + b'3i 4f')
        nbytes = write_header(card_type, nfields, ncards, key, op2, op2_ascii)
        raise NotImplementedError('MPC')

        #for load in cards:
            #data = [load.sid, load.node_id, load.Cid(), load.mag] + list(load.xyz)
            #op2_ascii.write('  MPC data=%s\n' % str(data))
            #op2.write(spack.pack(*data))
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

def write_header(name, nfields, ncards, key, op2, op2_ascii):
    nvalues = nfields * ncards + 3 # +3 comes from the keys
    nbytes = nvalues * 4
    op2.write(pack('3i', *[4, nvalues, 4]))
    op2.write(pack('i', nbytes)) #values, nbtyes))

    op2.write(pack('3i', *key))
    op2_ascii.write('%s %s\n' % (name, str(key)))
    return nbytes
