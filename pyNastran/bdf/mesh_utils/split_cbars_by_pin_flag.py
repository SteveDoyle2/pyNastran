"""
defines:
 - model, pin_flag_map = split_cbars_by_pin_flag(
       bdf_filename, pin_flags_filename=None, bdf_filename_out=None)
"""
from __future__ import print_function
from collections import defaultdict
from six import iteritems, string_types, integer_types
import numpy as np
from pyNastran.bdf.bdf import read_bdf

def make_release_map(model, bar_beam_eids, debug=False):  # pragma: no cover
    """may be used later"""
    bar_types = {
        # PBAR
        'bar' : [],

        # PBEAML/PBARL
        "ROD": [],
        "TUBE": [],
        "TUBE2" : [],
        "I": [],
        "CHAN": [],
        "T": [],
        "BOX": [],
        "BAR": [],
        "CROSS": [],
        "H": [],
        "T1": [],
        "I1": [],
        "CHAN1": [],
        "Z": [],
        "CHAN2": [],
        "T2": [],
        "BOX1": [],
        "HEXA": [],
        "HAT": [],
        "HAT1": [],
        "DBOX": [],  # was 12

        # PBEAM
        'beam' : [],

        # PBEAML specfic
        "L" : [],
    }  # for GROUP="MSCBML0"
    allowed_types = [
        'BAR', 'BOX', 'BOX1', 'CHAN', 'CHAN1', 'CHAN2', 'CROSS', 'DBOX',
        'H', 'HAT', 'HAT1', 'HEXA', 'I', 'I1', 'L', 'ROD',
        'T', 'T1', 'T2', 'TUBE', 'TUBE2', 'Z', 'bar', 'beam',
    ]

    # bar_types['bar'] = [ [...], [...], [...] ]
    #bar_types = defaultdict(lambda : defaultdict(list))

    found_bar_types = set([])
    #neids = len(self.element_ids)
    for bar_type, data in iteritems(bar_types):
        eids = []
        lines_bar_y = []
        lines_bar_z = []
        bar_types[bar_type] = (eids, lines_bar_y, lines_bar_z)


        nid_release_map = defaultdict(list)

        #debug = True
        bar_nids = set([])
        #print('bar_beam_eids = %s' % bar_beam_eids)
        for eid in bar_beam_eids:
            if eid not in self.eid_map:
                model.log.error('eid=%s is not a valid bar/beam element...' % eid)
                if debug:
                    print('eid=%s is not a valid bar/beam element...' % eid)
                continue
            ieid = self.eid_map[eid]
            elem = model.elements[eid]
            pid = elem.pid
            assert not isinstance(pid, integer_types), elem
            if pid.type in ['PBAR', 'PBEAM']:
                bar_type = 'bar'
            elif pid.type in ['PBEAM']:
                bar_type = 'beam'
            elif pid.type in ['PBARL', 'PBEAML']:
                bar_type = pid.Type
            else:
                if debug:
                    print('NotImplementedError(pid)')
                raise NotImplementedError(pid)
            #print('bar_type =', bar_type)

            if debug:
                print('%s' % elem)
                print('  bar_type =', bar_type)
            found_bar_types.add(bar_type)

            (nid1, nid2) = elem.node_ids
            bar_nids.update([nid1, nid2])
            node1 = model.nodes[nid1]
            node2 = model.nodes[nid2]
            n1 = node1.get_position()
            n2 = node2.get_position()
            #centroid = (n1 + n2) / 2.
            i = n2 - n1
            Li = np.linalg.norm(i)
            ihat = i / Li

            if elem.pa != 0:
                #assert elem.pa in [], elem.pa
                nid_release_map[nid1].append((eid, elem.pa))
            if elem.pb != 0:
                nid_release_map[nid2].append((eid, elem.pb))
        return nid_release_map


def split_cbars_by_pin_flag(bdf_filename,
                            pin_flags_filename=None, bdf_filename_out=None):
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
    if isinstance(bdf_filename, string_types):
        model = read_bdf(bdf_filename, xref=False)
    else:
        model = bdf_filename
    card_types = ['CBAR', 'CBEAM']
    out = model.get_card_ids_by_card_types(card_types=card_types)
    bar_beam_eids = out['CBAR'] + out['CBEAM']

    pin_flag_map = {}
    #with open('modes_old.bdf', 'w') as bdf_file:
    nid_new = max(model.nodes) + 1
    eid_new = max(model.elements) + 1
    for eid in sorted(list(model.elements.keys())):
        elem = model.elements[eid]
        if elem.type == 'CBAR':
            pa = elem.pa
            pb = elem.pb
            if pa == 0 and pa == pb:
                pin_flag_map[eid] = 0
                continue

            ga, gb = elem.node_ids
            xyz_a = model.nodes[ga].xyz
            xyz_b = model.nodes[gb].xyz
            xyz_mid = (xyz_a + xyz_b) / 2.
            model.add_grid(nid_new, cp=0, xyz=xyz_mid, cd=0, ps='', seid=0,
                           comment='')

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
            model.add_cbar(eid, pid, ga, nid_new, x, g0, offt=offt, pa=0, pb=0,
                           wa=wa, wb=wc, comment='')
            model.add_cbar(eid_new, pid, nid_new, gb, x, g0, offt=offt, pa=0, pb=0,
                           wa=wc, wb=wb, comment='')
            pin_flag_map[eid] = pa
            pin_flag_map[eid_new] = pb
            nid_new += 1
            eid_new += 1
        elif elem.type in ['CELAS1', 'CQUAD4']:
            pin_flag_map[eid] = 0
            #eid_new += 1
        else:
            raise RuntimeError(elem.type)
    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    assert len(model.elements) == len(pin_flag_map)

    if pin_flags_filename is not None:
        pin_flag_remap = {
            0 : 0,
            6 : 1,
            45 : 2,
            46 : 3,
            56 : 4,
            456 : 5,
            16 : 6,
        }
        with open(pin_flags_filename, 'w') as pin_flags_file:
            pin_flags_file.write('# Eid(%i), flag(%i)\n')
            for eid, pin_flag in sorted(iteritems(pin_flag_map)):
                pin_flag2 = pin_flag_remap[pin_flag]
                pin_flags_file.write('%s,%s\n' % (eid, pin_flag2))

    #nid_release_map = make_release_map(model, bar_beam_eids)
    return model, pin_flag_map

def main():  # pragma: no cover
    #bdf_filename = 'modes_old.bdf'
    #bdf_filename_out = 'modes_split_old.bdf'

    bdf_filename = 'merged_model_newgeom_fatigue.bdf'
    bdf_filename_out = 'split_newgeom_fatigue.bdf'

    pin_flags_filename = 'pin_flags.csv'
    split_cbars_by_pin_flag(bdf_filename, pin_flags_filename, bdf_filename_out=bdf_filename_out)

if __name__ == '__main__':  # pragma: no cover
    main()
