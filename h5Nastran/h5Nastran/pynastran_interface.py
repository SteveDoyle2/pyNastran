"""
Interfaces to pynastran bdf file.  Returns a dict of all card objects where the key is the card name.
"""

from __future__ import print_function, absolute_import

from six import iteritems


def get_bdf_cards(bdf):
    """

    :type bdf: pyNastran.bdf.bdf.BDF
    :return:
    """

    cards = {}

    def _add_cards_from_dict(obj):
        for key, value in iteritems(obj):
            if isinstance(value, (list, tuple)):
                values = value

                for value in values:
                    card_type = value.type

                    try:
                        cards[card_type]
                    except KeyError:
                        cards[card_type] = {}

                    try:
                        cards[card_type][key].append(value)
                    except KeyError:
                        cards[card_type][key] = [value]
            else:
                card_type = value.type

                try:
                    cards[card_type]
                except KeyError:
                    cards[card_type] = {}

                cards[card_type][key] = value

    def _add_cards_from_list(obj):
        for value in obj:
            if isinstance(value, (list, tuple)):
                values = value

                for value in values:
                    card_type = value.type

                    try:
                        data = cards[card_type]
                    except KeyError:
                        data = cards[card_type] = []

                    for _ in value:
                        data.append(value)

            else:
                card_type = value.type

                try:
                    cards[card_type].append(value)
                except KeyError:
                    cards[card_type] = [value]

    attrs = [
        'nodes', 'coords', 'elements', 'properties', 'rigid_elements', 'plotels', 'masses', 'properties_mass',
        'materials', 'thermal_materials', 'MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'creep_materials', 'hyperelastic_materials',
        'load_combinations', 'loads', 'tics', 'dloads', 'dload_entries',
        'nlpcis', 'nlparms', 'rotors', 'tsteps', 'tstepnls', 'transfer_functions', 'delays',
        'aeros', 'caeros', 'paeros', 'splines',
        'aecomps', 'aefacts', 'aelinks', 'aeparams', 'aesurf', 'aesurfs', 'aestats', 'trims', 'divergs', 'csschds', 'mkaeros', 'monitor_points',
        'aero', 'flfacts', 'flutters', 'gusts',
        'bcs', 'phbdys', 'convection_properties', 'tempds',
        'bcrparas', 'bctadds', 'bctparas', 'bctsets', 'bsurf', 'bsurfs',
        'suport1', 'suport', 'se_suport',
        'spcadds', 'spcs', 'mpcadds', 'mpcs',
        'dareas', 'dphases', 'pbusht', 'pdampt', 'pelast', 'frequencies',
        'dmis', 'dmigs', 'dmijs', 'dmijis', 'dmiks',
        'sets', 'usets', 'asets', 'bsets', 'csets', 'qsets', 'se_sets', 'se_usets', 'se_bsets', 'se_csets', 'se_qsets',
        'tables', 'tables_d', 'tables_m', 'random_tables', 'tables_sdamping',
        'methods', 'cMethods',
        'dconadds', 'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps', 'dtable', 'doptprm', 'dequations', 'dvprels', 'dvmrels', 'dvcrels', 'dscreen', 'dvgrids',
    ]

    for attr in attrs:
        _attr = getattr(bdf, attr)
        if attr is None:
            continue
        if isinstance(_attr, dict):
            _add_cards_from_dict(_attr)
        elif isinstance(_attr, list):
            _add_cards_from_list(_attr)

    return cards


def get_bdf_add_methods():
    add_methods = {}

    from pyNastran.bdf.bdf_interface.add_card import AddCards

    methods = dir(AddCards)

    for method in methods:
        if method.startswith('add_') and callable(getattr(AddCards, method)):
            card_id = method[4:].upper()
            add_methods[card_id] = method

    return add_methods



if __name__ == '__main__':
    from pyNastran.bdf.bdf import BDF

    # bdf = BDF()
    # bdf.read_bdf(r'../new_bdf.bdf')
    #
    # cards = get_bdf_cards(bdf)
    #
    # pbeam = cards['PBEAM'][200007]
    #
    # print(pbeam.repr_fields())
    #
    # print(cards['MOMENT'])

    add_methods = get_bdf_add_methods()
    print(add_methods)




