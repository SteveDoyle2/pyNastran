"""
defines:
 - model, pin_flag_map = split_cbars_by_pin_flag(
       bdf_filename, pin_flags_filename=None, bdf_filename_out=None)

"""
from pyNastran.bdf.bdf import read_bdf


def split_cbars_by_pin_flag(bdf_filename,
                            pin_flags_filename=None, bdf_filename_out=None,
                            debug=False):
    """
    Splits bar elements if they have a pin flag.  That way you can
    each side of the element (A/B) a unique color based on the pin flag.
    This doesn't split non-pin flagged bars.

    Parameters
    ----------
    bdf_filename : str; BDF
        str : use the ``read_bdf`` method
        BDF : assume it's a model
    pin_flags_filename : str; default=None
        the pin flag file to write
        this is optional as you may need to define your own map
    bdf_filename_out : str; default=None
        write the updated deck
    debug : bool/None; default=True
        used to set the logger if no logger is passed in
            True:  logs debug/info/error messages
            False: logs info/error messages
            None:  logs error messages

    Returns
    -------
    model : BDF()
        the BDF object
    pin_flag_map : dict[eid] : pin_flag
        eid : int
            the element id
        pin_flag : int
            the bar pin flag

    """
    if isinstance(bdf_filename, str):
        model = read_bdf(bdf_filename, xref=False, debug=debug)
    else:
        model = bdf_filename
    card_types = ['CBAR', 'CBEAM']
    out = model.get_card_ids_by_card_types(card_types=card_types)
    bar_beam_eids = out['CBAR'] + out['CBEAM']

    pin_flag_map = {}

    nid_new = max(model.nodes) + 1
    eid_new = max(model.elements) + 1
    for eid in sorted(list(model.elements.keys())):
        elem = model.elements[eid]
        etype = elem.type
        if etype in ['CBAR', 'CBEAM']:
            pa = elem.pa
            pb = elem.pb
            if pa == 0 and pa == pb:
                pin_flag_map[eid] = 0
                continue

            ga, gb = elem.node_ids
            xyz_a = model.nodes[ga].xyz
            xyz_b = model.nodes[gb].xyz
            xyz_mid = (xyz_a + xyz_b) / 2.
            model.add_grid(nid_new, xyz_mid, cp=0, cd=0, ps='', seid=0,
                           comment='')

            comment = elem.comment
            ga = elem.ga
            gb = elem.gb
            pid = elem.pid
            x = elem.x
            g0 = elem.g0
            offt = elem.offt
            wa = elem.wa
            wb = elem.wb
            wc = (wa + wb) / 2.

            del model.elements[eid]
            if etype == 'CBAR':
                model.add_cbar(eid, pid, [ga, nid_new], x, g0, offt=offt, pa=0, pb=0,
                               wa=wa, wb=wc, comment=comment)
                model.add_cbar(eid_new, pid, [nid_new, gb], x, g0, offt=offt, pa=0, pb=0,
                               wa=wc, wb=wb, comment='')
            else:  # CBEAM
                model.add_cbeam(eid, pid, [ga, nid_new], x, g0, offt=offt, pa=0, pb=0,
                                wa=wa, wb=wc, sa=0, sb=0, bit=None, comment=comment)
                model.add_cbeam(eid_new, pid, [nid_new, gb], x, g0, offt=offt, pa=0, pb=0,
                                wa=wc, wb=wb, sa=0, sb=0, bit=None, comment='')
            pin_flag_map[eid] = pa
            pin_flag_map[eid_new] = pb
            nid_new += 1
            eid_new += 1
        else:
            pin_flag_map[eid] = 0
            #eid_new += 1
    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    assert len(model.elements) == len(pin_flag_map)
    #nid_release_map = make_release_map(model, bar_beam_eids)

    write_pin_flag_map(pin_flag_map, pin_flags_filename)
    return model, pin_flag_map

def write_pin_flag_map(pin_flag_map, pin_flags_filename):
    """writes the pin flag map"""
    if pin_flags_filename is not None:
        pin_flag_remap = {
            0 : 0,
            6 : 1,
            45 : 2,
            46 : 3,
            56 : 4,
            456 : 5,
            16 : 6,
            1456 : 7,
        }
        with open(pin_flags_filename, 'w') as pin_flags_file:
            pin_flags_file.write('# Eid(%i), flag(%i)\n')
            for eid, pin_flag in sorted(pin_flag_map.items()):
                pin_flag2 = pin_flag_remap[pin_flag]
                pin_flags_file.write('%s,%s\n' % (eid, pin_flag2))
