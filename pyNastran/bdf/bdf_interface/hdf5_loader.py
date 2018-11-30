"""Defines various helper functions for loading a HDF5 BDF file"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from collections import defaultdict
from six import string_types, StringIO, text_type, binary_type
import numpy as np
from pyNastran.utils.dict_to_h5py import _add_list_tuple, _cast, integer_types, float_types, string_types
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP
from pyNastran.utils import object_attributes

# dict[key] : [value1, value2, ...]
dict_int_list_obj_attrs = [
    'spcs', 'spcadds',
    'mpcs', 'mpcadds',
    'loads', 'load_combinations',
    'dloads', 'dload_entries',
    #'usets',
    'nsms', 'nsmadds',
    #'dconstrs',
    'frequencies',
    'bcs', 'transfer_functions',
    'dvgrids',
]

# dict[key] : value
dict_int_obj_attrs = [
    # are handled explictly
    #'elements',
    #'nodes',
    #'coords',
    #'materials',
    #'properties',
    #'masses',
    #'tables',
    #'methods',
    #'creep_materials', 'csschds',
    #'flutters',
    #'gusts',
    #'trims',
    'MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',

    # TODO: don't work
    #'reject_count',
    #'dresps',

    'rigid_elements',
    'aecomps', 'aefacts', 'aelinks', 'aelists', 'aeparams',
    'aestats', 'aesurf', 'aesurfs', 'ao_element_flags', 'bconp', 'bcrparas', 'bctadds',
    'bctparas', 'bctsets', 'blseg', 'bsurf', 'bsurfs', 'cMethods',
    'convection_properties',
    'csuper', 'csupext', 'dareas',
    'dconadds', 'ddvals', 'delays', 'dequations', 'desvars', 'divergs', 'dlinks',
    'dmigs', 'dmijis', 'dmijs', 'dmiks', 'dmis', 'dphases',
    'dscreen', 'dti', 'dvcrels', 'dvmrels', 'dvprels',
    'epoints', 'flfacts',
    'gridb',
    'nlparms', 'nlpcis',
    'normals',
    'nxstrats', 'paeros',
    'pbusht', 'pdampt', 'pelast', 'phbdys', 'plotels', 'points',
    'properties_mass',
    'radcavs', 'radmtx', 'random_tables',
    'ringaxs', 'ringfl',
    'rotors',
    'se_sets', 'se_usets', 'sebndry', 'sebulk', 'seconct', 'seelt',
    'seexcld', 'selabel', 'seload', 'seloc', 'sempln', 'senqset', 'setree', 'sets',
    'spcoffs',
    'spoints',
    'suport1',
    'tables_d', 'tables_m', 'tables_sdamping', 'tempds',
    'tics',
    'tstepnls', 'tsteps',
    'view3ds', 'views',
]

scalar_obj_keys = [
    # required
    'aero', 'aeros', 'axic', 'axif', 'baror', 'beamor',
    'doptprm', 'dtable', 'grdset', 'radset', 'seqgp',
    #'case_control_deck',
    #'zona',
]

scalar_keys = [
    # basic types
    'bdf_filename',
    '_auto_reject', '_encoding', '_iparse_errors', '_is_axis_symmetric', '_is_cards_dict',
    '_is_dynamic_syntax', '_is_long_ids', '_ixref_errors', '_nastran_format', '_nparse_errors',
    '_nxref_errors', '_sol', '_stop_on_duplicate_error', '_stop_on_parsing_error',
    '_stop_on_xref_error',
    #'_unique_bulk_data_cards',
    #'cards_to_read',
    '_xref', 'active_filename',
    'dumplines', 'echo', 'force_echo_off', 'include_dir',
    #'is_bdf_vectorized',
    #'is_long_ids',
    'is_msc', 'is_nx', 'punch',
    'read_includes', 'save_file_structure', 'sol', 'sol_iline', 'nastran_format',
    'is_superelements', 'is_zona', 'sol_method', 'debug',

    # not sure
    #'nid_cp_cd', 'xyz_cid0',

    # removed
    #'ncaeros', 'ncoords', 'nnodes', 'nelements', 'nproperties', 'nmaterials',
    #'point_ids', 'wtmass', 'log',
]

list_keys = [
    # required
    'active_filenames', 'executive_control_lines', 'case_control_lines', #'initial_superelement_models',
    #'reject_cards',

    # maybe...
    #'_duplicate_coords', '_duplicate_elements', '_duplicate_masses', '_duplicate_materials',
    #'_duplicate_nodes', '_duplicate_properties', '_duplicate_thermal_materials', '_stored_parse_errors',
    #'_stored_xref_errors',
    'system_command_lines', #'units', 'xyz_limits',

    # removed
    #'coord_ids', 'element_ids', 'material_ids', 'node_ids', 'property_ids', 'caero_ids',
    #'special_cards',
]

list_obj_keys = [
    # TODO: required
    'asets', 'bsets', 'csets', 'omits', 'qsets',
    'mkaeros',
    'monitor_points',
    'suport',
    'se_bsets', 'se_csets', 'se_qsets', 'se_suport',
]

dict_attrs = [
    # required
    'params',

    # removed
    #'_solmap_to_value',
    #'card_count',
    #'_card_parser',
    #'_card_parser_prepare',
    #'_slot_to_type_map',
    #'_type_to_id_map',
    #'_type_to_slot_map',
]

def export_to_hdf5_file(hdf5_file, model, exporter=None):
    attrs = object_attributes(model, mode='both', keys_to_skip=None)
    encoding = model.get_encoding()

    if 'GRID' in model.card_count:
        model.log.debug('exporting nodes to hdf5')
        node_group = hdf5_file.create_group('nodes')
        grid_group = node_group.create_group('GRID')
        nids = model._type_to_id_map['GRID']
        if len(nids) == 0:
            assert len(model.nodes) == 0, len(model.nodes)
        CARD_MAP['GRID'].export_to_hdf5(grid_group, model, nids)

    _hdf5_export_group(hdf5_file, model, 'coords', encoding, debug=False)
    _hdf5_export_elements(hdf5_file, model, encoding)

    # explicit groups
    #
    # these are broken down by card type
    # they came from dict_int_obj_attrs
    groups_to_export = [
        'properties', 'masses',

        # materials
        'materials', 'thermal_materials', 'creep_materials', 'hyperelastic_materials',
        #'MATS1',
        #'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',

        # aero
        'caeros', 'splines', 'flutters', 'trims', 'csschds', 'gusts',

        # other
        'methods', 'tables',
    ]
    for group_name in groups_to_export:
        _hdf5_export_group(hdf5_file, model, group_name, encoding)


    dict_int_attrs = [
        # required
        '_dmig_temp',
        'include_filenames',
        'rsolmap_to_str',
        'superelement_models',
        'values_to_skip',

        # removed
        #'nid_map',
        #'subcases',
    ]

    _export_dict_int_obj_attrs(model, hdf5_file, encoding)
    _export_dict_int_list_obj_attrs(model, hdf5_file, encoding)

    #for key in scalar_obj_keys:
        #value = getattr(model, key)
        #hdf5_file.create_dataset(key, value)

    scalars_keys_to_analyze = []
    for key in scalar_keys:
        if hasattr(model, key):
            scalars_keys_to_analyze.append(key)

    if scalars_keys_to_analyze:
        scalar_group = hdf5_file.create_group('minor_attributes')
        encoding = model.get_encoding()
        scalar_group.create_dataset('encoding', data=encoding)
        for key in scalars_keys_to_analyze:
            value = getattr(model, key)

            if value is None:
                continue
                #print('None: %s %s' % (key, value))

            #elif isinstance(value, bool):
                #print('bool: %s %s' % (key, value))
            elif isinstance(value, (integer_types, float_types, string_types, np.ndarray)):
                try:
                    scalar_group.create_dataset(key, data=value)
                except TypeError:  # pragma: no cover
                    print('key=%r value=%s type=%s' % (key, str(value), type(value)))
                    raise
            #elif isinstance(value, set):
                #_add_list_tuple(hdf5_file, key, value, 'set', model.log)
            elif isinstance(value, StringIO):
                pass
            else:
                #print(key, value)
                scalar_group.create_dataset(key, data=value)

    if model.params:
        skip_attrs = ['comment', '_field_map']
        group = hdf5_file.create_group('params')
        for key, param in model.params.items():
            _h5_export_class(group, model, key, param, skip_attrs, debug=False)

    for key in ['case_control_lines', 'executive_control_lines', 'system_command_lines', 'active_filenames']:
        list_str = getattr(model, key)
        if not len(list_str):
            continue
        list_bytes = [line.encode(encoding) for line in list_str]
        scalar_group.create_dataset(key, data=list_bytes)

    if len(model.reject_lines):
        #print('reject_lines', model.reject_lines)
        for i, list_str in enumerate(model.reject_lines):
            #print(list_str, type(list_str))
            list_bytes = [line.encode(encoding) for line in list_str]
            scalar_group.create_dataset(str(i), data=list_bytes)

    skip_attrs = ['comment', '_field_map']
    for key in scalar_obj_keys:
        value = getattr(model, key)
        if value is None:
            #print('None: %s %s' % (key, value))
            pass
        else:
            _h5_export_class(hdf5_file, model, key, value, skip_attrs, debug=False)

    _export_list_keys(model, hdf5_file, list_keys)
    _export_list_obj_keys(model, hdf5_file, list_obj_keys, encoding)

    cards_to_read = [key.encode(encoding) for key in list(model.cards_to_read)]
    cards_to_read = list(cards_to_read)
    cards_to_read.sort()
    hdf5_file.create_dataset('cards_to_read', data=cards_to_read)
    #dict_keys2 = []
    #list_keys2 = []
    #other_keys2 = []
    #for key in attrs:
        #value = getattr(model, key)
        #if isinstance(value, dict):
            #dict_keys2.append(key)
        #elif isinstance(value, list):
            #list_keys2.append(key)
        #else:
            #other_keys2.append(key)

    #print('dict_keys2 = %s' % (set(dict_keys) - set(dict_keys2)))
    #print('list_keys2 = %s' % (set(list_keys) - set(list_keys2)))
    #print('other_keys2 = %s' % (set(other_keys) - set(other_keys2)))
    #asd

def _export_dict_int_obj_attrs(model, hdf5_file, encoding):
    cards = set(list(CARD_MAP.keys()))
    for attr in dict_int_obj_attrs:
        dict_obj = getattr(model, attr)
        if not len(dict_obj):
            continue

        #model.log.info(attr)
        try:
            group = hdf5_file.create_group(attr) # 'gusts'
        except ValueError:  # pragma: no cover
            model.log.error('cant create %r' % attr)
            raise
        _hdf5_export_object_dict(group, model, attr, dict_obj, dict_obj.keys(), encoding)

def _export_dict_int_list_obj_attrs(model, hdf5_file, encoding):
    for attr in dict_int_list_obj_attrs:
        dict_obj = getattr(model, attr) # spcs
        if not len(dict_obj):
            continue

        #model.log.info(attr)
        try:
            group = hdf5_file.create_group(attr) # 'spcs'
        except ValueError:  # pragma: no cover
            model.log.error('cant create %r' % attr)
            raise

        keys = list(dict_obj.keys())
        keys.sort()
        group.create_dataset('keys', data=keys)
        for spc_id, spcs in sorted(dict_obj.items()):
            id_group = group.create_group(str(spc_id))

            card_types = defaultdict(list)
            for spc in spcs:
                card_types[spc.type].append(spc)
            for card_type, spcs in card_types.items():
                card_group = id_group.create_group(card_type)

                class_obj = spcs[0]
                if hasattr(class_obj, 'export_to_hdf5'):
                    class_obj.export_to_hdf5(card_group, model, spcs)
                else:
                    indices = list(range(len(spcs)))
                    _hdf5_export_object_dict(card_group, model,
                                             '%s/id=%s/%s' % (attr, spc_id, card_type),
                                             spcs, indices, encoding)

def _export_list_keys(model, hdf5_file, list_keys):
    for attr in list_keys:
        #print('list_key: %s' % attr)
        list_obj = getattr(model, attr) # active_filenames
        if not len(list_obj):
            continue
        #model.log.info(attr)

        #try:
            #group = hdf5_file.create_group(attr) # 'active_filenames'
        #except ValueError:
            #model.log.error('cant create %r' % attr)
            #raise

        if isinstance(list_obj, list):
            Type = 'list'
        elif isinstance(list_obj, tuple):
            Type = 'tuple'
        else:
            raise NotImplementedError(type(list_obj))
        #indices = list(range(len(list_keys)))
        #group.create_dataset('keys', data=keys)

        if isinstance(list_obj[0], (int, float, string_types)):
            try:
                _add_list_tuple(hdf5_file, attr, list_obj, Type, model.log)
            except TypeError:  # pragma: no cover
                print(list_obj)
                raise
        #elif isinstance(list_obj[0], list):
            #group = hdf5_file.create_group(attr)
            #group.attrs['type'] = Type
            #for keyi, valuei in enumerate(list_obj):
                ##sub_group = hdf5_file.create_group(str(keyi))
                ##group
                #_add_list_tuple(group, str(keyi), valuei, Type, model.log)
        else:
            raise NotImplementedError(type(list_obj[0]))
        #_hdf5_export_object_dict(group, model, attr, list_obj, indices, encoding)


def _export_list_obj_keys(model, hdf5_file, list_obj_keys, encoding):
    for attr in list_obj_keys:
        #print('list_key: %s' % attr)
        list_obj = getattr(model, attr) # active_filenames
        if not len(list_obj):
            continue
        model.log.info(attr)

        try:
            group = hdf5_file.create_group(attr) # 'active_filenames'
        except ValueError:  # pragma: no cover
            model.log.error('cant create %r' % attr)
            raise

        #if isinstance(list_obj, list):
            #Type = 'list'
        #elif isinstance(list_obj, tuple):
            #Type = 'tuple'
        #else:
            #raise NotImplementedError(type(list_obj))

        indices = list(range(len(list_obj)))
        #group.create_dataset('keys', data=indices)
        #try:
            #_add_list_tuple(hdf5_file, attr, list_obj, Type, model.log)
        #except TypeError:
            #print(list_obj)
            #raise
        _hdf5_export_object_dict(group, model, attr, list_obj, indices, encoding)


def _h5_export_class(sub_group, model, key, value, skip_attrs, debug=True):
    #sub_groupi = sub_group.create_group('values')
    class_group = sub_group.create_group(str(key))
    try:
        class_group.attrs['type'] = value.type
    except:  # pragma: no cover
        model.log.error('ERROR: key=%s value=%s' % (key, value))
        raise

    #if hasattr(value, 'get_h5attrs'):
        #getattrs
    if hasattr(value, 'object_attributes'):
        keys_to_skip = []
        if hasattr(value, '_properties'):
            keys_to_skip = value._properties
        h5attrs = value.object_attributes(mode='both', keys_to_skip=keys_to_skip)
        if hasattr(value, '_properties'):
            h5attrs.remove('_properties')
        #sub_group = hdf5_file.create_group(key)
    else:
        raise NotImplementedError(value)

    #if hasattr(value, '_properties'):
        #print(value.type, value._properties)
        #if debug:
            #print(h5attrs)
        #for prop in value._properties:
            #try:
                #h5attrs.remove(prop)
            #except:
                #print('cant remove %s' % prop)
                #print(value)
                #raise
        #h5attrs.remove('_properties')

    if debug:
        model.log.debug(value)
    for h5attr in h5attrs:
        if '_ref' in h5attr or h5attr in skip_attrs:
            continue
        class_value = getattr(value, h5attr)
        if class_value is None:
            continue

        is_none = False
        if debug:
            model.log.info('%s %s %s' % (key, h5attr, class_value))
        if isinstance(class_value, list):
            if len(class_value) == 0: # empty list
                continue
            for class_valuei in class_value:
                if class_valuei is None:  # PAERO2 : lth
                    is_none = True
                    break

        if is_none:
            model.log.warning('skipping %s attribute: %s %s %s' % (value.type, key, h5attr, class_value))
        else:
            try:
                class_group.create_dataset(h5attr, data=class_value)
            except ValueError:  # pragma: no cover
                print(h5attr, class_group)
                raise
            except TypeError:
                # contains unicode
                class_group.attrs['type'] = 'list'
                param_group = class_group.create_group(h5attr)
                for i, valuei in enumerate(class_value):
                    if isinstance(valuei, text_type):
                        param_group.create_dataset(str(i), data=valuei.encode('ascii'))
                    else:
                        param_group.create_dataset(str(i), data=valuei)
                #if isinstance(class_value, list):
                    #print('type(value[0] =', class_value, type(class_value[0]))
                    #raise

    #assert isinstance(key, int), 'key=%s value=%s' % (key, value)

    if isinstance(value, list):
        raise NotImplementedError('list: %s' % value)
        #for valuei in value:
            #if valuei.type not in cards:
                #msg = 'key=%s type=%s value=%s=' % (key, valuei.type, value)
                #print(msg)
        #continue

    #if attr in ['elements']:
        #continue
    #if value.type not in cards:
        #msg = 'key=%s type=%s value=%s=' % (key, value.type, value)
        #print(msg)

def _hdf5_export_elements(hdf5_file, model, encoding):
    """
    exports the elements to an hdf5_file

    TODO: not done
    """
    etypes_actual = []
    etypes = model._slot_to_type_map['elements']
    for card_name in model.card_count:
        #CTRIA3, CQUAD4
        #CONROD
        #CBUSH
        #CBEAM
        #CPENTA, CHEXA
        if card_name in etypes:
            #model.log.debug(card_name)
            etypes_actual.append(card_name)
            continue

    if etypes_actual:
        elements_group = hdf5_file.create_group('elements')
        def save_solids(etype, slot_name):
            element_group = elements_group.create_group(etype)
            eids = model._type_to_id_map[etype]
            CARD_MAP[slot_name].export_to_hdf5(element_group, model, eids)

        solids = [
            ('CTETRA', 'CTETRA4'),
            ('CPENTA', 'CPENTA6'),
            ('CPYRAM', 'CPYRAM5'),
            ('CHEXA', 'CHEXA20'),
        ]
        for card_name, slot_name in solids:
            if card_name in model.card_count:
                save_solids(card_name, slot_name)
                etypes_actual.remove(card_name)

        for card_type in etypes_actual:
            element_group = elements_group.create_group(card_type)
            eids = model._type_to_id_map[card_type]
            class_obj = CARD_MAP[card_type]

            if hasattr(class_obj, 'export_to_hdf5'):
                class_obj.export_to_hdf5(element_group, model, eids)
            else:
                _hdf5_export_object_dict(element_group, model, card_type,
                                         model.elements, eids, encoding)

def _hdf5_export_group(hdf5_file, model, group_name, encoding, debug=False):
    """
    exports the properties to an hdf5_file
    """
    data_dict = getattr(model, group_name) # self.properties
    if debug:
        model.log.debug('%s %s' % (group_name, data_dict))

    types_actual = []
    types = model._slot_to_type_map[group_name]
    if debug:
        model.log.debug('card_count = %s' % model.card_count)
        model.log.debug('types = %s' % types)
    for card_name in model.card_count:
        #PSHELL
        if card_name in types:
            types_actual.append(card_name)
            continue

    if types_actual:
        model.log.debug('exporting %s to hdf5' % group_name)
        if debug:
            print('types_actual =', types_actual)
        group = hdf5_file.create_group(group_name)
        for card_type in types_actual:
            sub_group = group.create_group(card_type)
            ids = model._type_to_id_map[card_type]
            if debug:
                print(ids)
            assert len(ids) > 0, '%s : %s' % (card_type, ids)
            class_obj = CARD_MAP[card_type]
            #print(class_obj)
            if hasattr(class_obj, 'export_to_hdf5'):
                class_obj.export_to_hdf5(sub_group, model, ids)
            else:
                _hdf5_export_object_dict(sub_group, model, card_type, data_dict, ids, encoding)

def _hdf5_export_object_dict(group, model, name, obj_dict, keys, encoding):
    i = 0
    sub_group = group.create_group('values')
    assert isinstance(name, string_types), 'name=%s; type=%s' % (name, type(name))

    skip_attrs = ['comment', '_field_map']

    keys_write = list(keys)
    if isinstance(keys_write[0], text_type):
        keys_write = list([key.encode(encoding) if isinstance(key, text_type) else key for key in list(keys_write)])

    for key in keys:
        value = obj_dict[key]
        #if isinstance(value, text_type):
            #value = value.encode(encoding)

        try:
            _h5_export_class(sub_group, model, key, value, skip_attrs, debug=False)
        except:  # pragma: no cover
            # for debugging
            sub_group2 = group.create_group('values2')
            _h5_export_class(sub_group2, model, key, value, skip_attrs, debug=True)
        i += 1

    #group.attrs['type'] = class_name
    #print('%s keys = %s' % (name, keys))
    try:
        group.create_dataset('keys', data=keys_write)
    except TypeError:  # pragma: no cover
        print('name =', name)
        print('encoding =', encoding)
        print('keys =', keys)
        raise

# exporter
#-------------------------------------------------------------------------------------
# importer
def load_hdf5_file(h5_file, model):
    encoding = _cast(h5_file['minor_attributes']['encoding'])
    keys = h5_file.keys()

    mapper = {
        'elements' : hdf5_load_elements,
        'properties' : hdf5_load_properties,
        'coords' : hdf5_load_coords,
        'tables' : hdf5_load_tables,
        'methods' : hdf5_load_methods,
        'masses' : hdf5_load_masses,
        'materials' : hdf5_load_materials,

        'spcs' : hdf5_load_spcs,
        'spcadds' : hdf5_load_spcadds,
        'mpcs' : hdf5_load_mpcs,
        'mpcadds' : hdf5_load_mpcadds,

        'loads' : hdf5_load_loads,
        'load_combinations' : hdf5_load_load_combinations,
        'dloads' : hdf5_load_dloads,
        'dload_entries' : hdf5_load_dload_entries,
        'bcs' : hdf5_load_bcs,
        'transfer_functions' : hdf5_load_transfer_functions,
        'dvgrids': hdf5_load_dvgrids,

        'nsms' : hdf5_load_nsms,
        'nsmadds' : hdf5_load_nsmadds,
        'frequencies' : hdf5_load_frequencies,
    }
    generic_mapper = {
        'flutters' : hdf5_load_generic,
        'trims' : hdf5_load_generic,
        'csschds' : hdf5_load_generic,
        'gusts' : hdf5_load_generic,
        'caeros' : hdf5_load_generic,
        'splines' : hdf5_load_generic,
        'thermal_materials' : hdf5_load_generic,
        'creep_materials' : hdf5_load_generic,
        'hyperelastic_materials' : hdf5_load_generic,
        #'MATS1' : hdf5_load_generic,
        #'MATT1' : hdf5_load_generic,
        #'MATT2' : hdf5_load_generic,
        #'MATT3' : hdf5_load_generic,
        #'MATT4' : hdf5_load_generic,
        #'MATT5' : hdf5_load_generic,
        #'MATT8' : hdf5_load_generic,
        #'MATT9' : hdf5_load_generic,
    }
    for key in keys:
        group = h5_file[key]
        if key == 'nodes':
            nodes = {}
            grids = group['GRID']
            nids = _cast(grids['nid'])
            xyz = _cast(grids['xyz'])
            cp = _cast(grids['cp'])
            cd = _cast(grids['cd'])
            ps = _cast(grids['ps'])
            seid = _cast(grids['seid'])
            for nid, xyzi, cpi, cdi, psi, seidi in zip(nids, xyz, cp, cd, ps, seid):
                model.add_grid(nid, xyzi, cp=cpi, cd=cdi, ps=psi, seid=seidi, comment='')

        elif key in mapper:
            func = mapper[key]
            func(model, group, encoding)
        elif key in generic_mapper:
            func = generic_mapper[key]
            func(model, group, key, encoding)
        elif key in dict_int_obj_attrs:
            load_cards_from_keys_values(key, group)
        elif key in ['info', 'matrices'] or key.startswith('Subcase'): # op2
            continue
        elif key in ['cards_to_read']: # handled separaaately
            continue
        elif key == 'params':
            keys = list(group.keys())
            _load_cards_from_keys_values('params', group, keys)
        elif key == 'minor_attributes':
            keys_attrs = group.keys()
            for keyi in keys_attrs:
                sub_group = group[keyi]
                value = _cast(sub_group)
                try:
                    setattr(model, keyi, value)
                except AttributeError:  # pragma: no cover
                    model.log.warning('cant set minor_attributes/%s as %s' % (keyi, value))
                    raise

        elif key in ['case_control_lines', 'executive_control_lines', 'system_command_lines']:
            lst = _load_indexed_list_str(model, key, group, encoding)
            setattr(model, key, lst)

        elif key in ['active_filenames']:
            if 'value' not in group:
                lst = _load_indexed_list_str(model, key, group, encoding)
                continue

            lst = _cast(group['value']).tolist()
            #else:
            #except KeyError:  # pragma: no cover
                #print('group', group)
                #print('group.keys()', list(group.keys()))
                #raise

            if isinstance(lst[0], text_type):
                pass
            else:
                lst = [line.encode(encoding) for line in lst]
            setattr(model, key, lst)

        elif key == 'reject_lines':
            #print(group)
            lst = []
            ddd
            for key in group.keys():
                value = _cast(group[key])
                print('value', value)
                #assert isinstance(value, text_type), type(value)
                #lst.append(value)

        elif key in list_obj_keys:
            #model.log.info('key = %s' % key)
            #model.log.info('group = %s' % group)
            #model.log.info('group.keys() = %s' % list(group.keys()))
            keys = _cast(group['keys'])
            values = group['values']
            lst = [None] * len(keys)
            for key in values.keys():
                ikey = int(key)
                class_obj_hdf5 = values[key]
                card_type = _cast(class_obj_hdf5['type'])
                class_instance = _load_from_class(class_obj_hdf5, card_type)
                lst[ikey] = class_instance
            #model.log.info('keys = %s' % keys)
            #model.log.info('values = %s' % values)
            #model.log.info('values.keys() = %s' % values.keys())

        elif key in scalar_obj_keys:
            keys = list(group.keys())
            keys.remove('type')
            card_type = _cast(group['type'])
            class_instance = _load_from_class(group, card_type)
            setattr(model, key, class_instance)
        #elif key in scalar_keys:
            #value = _cast(group)
            #try:
                #setattr(model, key, value)
            #except AttributeError:
                #model.log.warning('cant set %r as %s' % (key, value))
                #raise

        #elif key in list_keys:
            #value = _cast(group)
            #try:
                #setattr(model, key, value)
            #except AttributeError:
                #model.log.warning('cant set %r as %s' % (key, value))
                #raise
        else:
            model.log.warning('skipping hdf5 load for %s' % key)

    cards_to_read = _cast(h5_file['cards_to_read'])
    cards_to_read = [key.decode(encoding) for key in cards_to_read]
    model.cards_to_read = set(list(cards_to_read))

def _load_indexed_list_str(model, key, group, encoding):
    lst = []
    for key in group.keys():
        value = _cast(group[key])
        lst.append(value)

    #try:
        #value0 = value[0]
    #except IndexError:  # pragma: no cover
        #print('key =', key)
        #print('value = %r' % value)
        #print('group =', group)
        #print('group.keys() =', list(group.keys()))
        #raise

    if isinstance(value, text_type):
        pass
    else:
        lst = [line.encode(encoding) for line in lst]
        assert isinstance(lst[0], text_type), type(lst[0])
    return lst

def hdf5_load_coords(model, coords_group, encoding):
    """loads the coords from an HDF5 file"""
    for card_type in coords_group.keys():
        coords = coords_group[card_type]
        if card_type in ['CORD2R', 'CORD2C', 'CORD2S']:
            if card_type == 'CORD2R':
                func = model.add_cord2r
            elif card_type == 'CORD2C':
                func = model.add_cord2c
            elif card_type == 'CORD2S':
                func = model.add_cord2s

            cids = _cast(coords['cid'])
            rids = _cast(coords['rid'])
            e1s = _cast(coords['e1'])
            e2s = _cast(coords['e2'])
            e3s = _cast(coords['e3'])
            for cid, rid, origin, zaxis, xzplane in zip(
                    cids, rids, e1s, e2s, e3s):
                func(cid, origin, zaxis, xzplane, rid=rid, comment='')
        elif card_type in ['CORD1R', 'CORD1C', 'CORD1S']:
            if card_type == 'CORD1R':
                func = model.add_cord1r
            elif card_type == 'CORD1C':
                func = model.add_cord1c
            elif card_type == 'CORD1S':
                func = model.add_cord1s

            cids = _cast(coords['cid'])
            nodes = _cast(coords['nodes'])
            for cid, (n1, n2, n3) in zip(cid, nodes):
                func(cid, n1, n2, n3, comment='')
        else:
            load_cards_from_keys_values('coords/%s' % card_type, coords)

def hdf5_load_tables(unused_model, group, encoding):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        load_cards_from_keys_values('tables/%s' % card_type, sub_group)

def hdf5_load_methods(unused_model, group, encoding):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'EIGRL':
            #pass
        load_cards_from_keys_values('methods/%s' % card_type, sub_group)

def hdf5_load_masses(model, group, encoding):
    for card_type in group.keys():
        masses = group[card_type]
        if card_type == 'CONM2':
            eid = _cast(masses['eid'])
            nid = _cast(masses['nid'])
            cid = _cast(masses['cid'])
            X = _cast(masses['X'])
            I = _cast(masses['I'])
            mass = _cast(masses['mass'])
            for eidi, nidi, cidi, Xi, Ii, massi in zip(eid, nid, cid, X, I, mass):
                model.add_conm2(eidi, nidi, massi, cid=cidi, X=Xi, I=Ii, comment='')
        elif card_type == 'CMASS2':
            eids = _cast(masses['eid'])
            mass = _cast(masses['mass'])
            nodes = _cast(masses['nodes']).tolist()
            components = _cast(masses['components'])
            for eid, massi, nids, (c1, c2) in zip(eids, mass, nodes, components):
                model.add_cmass2(eid, massi, nids, c1, c2, comment='')

        else:
            #model.add_cmass1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cmass3(eid, pid, nids, comment='')
            #model.add_cmass4(eid, mass, nids, comment='')
            #model.add_conm1(eid, nid, mass_matrix, cid=0, comment='')
            load_cards_from_keys_values('masses/%s' % card_type, masses)


def hdf5_load_materials(model, group, encoding):
    for card_type in group.keys():
        sub_group = group[card_type]
        if card_type == 'MAT1':
            mid = _cast(sub_group['mid'])
            E = _cast(sub_group['E'])
            G = _cast(sub_group['G'])
            nu = _cast(sub_group['nu'])
            rho = _cast(sub_group['rho'])
            a = _cast(sub_group['A'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            St = _cast(sub_group['St'])
            Sc = _cast(sub_group['Sc'])
            Ss = _cast(sub_group['Ss'])
            mcsid = _cast(sub_group['mcsid'])
            for midi, Ei, Gi, nui, rhoi, ai, trefi, gei, Sti, Sci, Ssi, mcsidi in zip(
                    mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid):
                model.add_mat1(midi, Ei, Gi, nui, rho=rhoi, a=ai, tref=trefi,
                               ge=gei, St=Sti, Sc=Sci, Ss=Ssi, mcsid=mcsidi, comment='')

        elif card_type == 'MAT2':
            mid = _cast(sub_group['mid'])
            G = _cast(sub_group['G'])
            rho = _cast(sub_group['rho'])
            a = _cast(sub_group['A'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            St = _cast(sub_group['St'])
            Sc = _cast(sub_group['Sc'])
            Ss = _cast(sub_group['Ss'])
            mcsid = _cast(sub_group['mcsid'])

            for (midi, (G11, G22, G33, G12, G13, G23), rhoi, (a1i, a2i, a3i),
                 trefi, gei, Sti, Sci, Ssi, mcsidi) in zip(
                     mid, G, rho, a, tref, ge, St, Sc, Ss, mcsid):
                if mcsidi == -1:
                    mcsidi = None
                model.add_mat2(midi, G11, G12, G13, G22, G23, G33, rho=rhoi,
                               a1=a1i, a2=a2i, a3=a3i, tref=trefi, ge=gei,
                               St=Sti, Sc=Sci, Ss=Ssi, mcsid=mcsidi, comment='')

        elif card_type == 'MAT3':
            mid = _cast(sub_group['mid'])
            ex = _cast(sub_group['Ex'])
            eth = _cast(sub_group['Eth'])
            ez = _cast(sub_group['Ez'])

            nuxth = _cast(sub_group['Nuxth'])
            nuzx = _cast(sub_group['Nuzx'])
            nuthz = _cast(sub_group['Nuthz'])
            gxz = _cast(sub_group['Gzx'])

            ax = _cast(sub_group['Ax'])
            ath = _cast(sub_group['Ath'])
            az = _cast(sub_group['Az'])

            rho = _cast(sub_group['rho'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            for (midi, exi, ethi, ezi, nuxthi, nuzxi, nuthzi,
                 rhoi, gzxi, axi, athi, azi, trefi, gei) in zip(
                     mid, ex, eth, ez, nuxth, nuzx, nuthz, rho, gxz, ax, ath, az, tref, ge):
                model.add_mat3(midi, exi, ethi, ezi, nuxthi, nuthzi, nuzxi, rho=rhoi,
                               gzx=gzxi, ax=axi, ath=athi, az=azi, tref=trefi, ge=gei, comment='')

        elif card_type == 'MAT8':
            mid = _cast(sub_group['mid'])
            e11 = _cast(sub_group['E11'])
            e22 = _cast(sub_group['E22'])
            nu12 = _cast(sub_group['Nu12'])
            g12 = _cast(sub_group['G12'])
            g1z = _cast(sub_group['G1z'])
            g2z = _cast(sub_group['G2z'])

            a1 = _cast(sub_group['A1'])
            a2 = _cast(sub_group['A2'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            rho = _cast(sub_group['rho'])

            xt = _cast(sub_group['Xt'])
            xc = _cast(sub_group['Xc'])
            yt = _cast(sub_group['Yt'])
            yc = _cast(sub_group['Yc'])
            s = _cast(sub_group['S'])

            f12 = _cast(sub_group['F12'])
            strn = _cast(sub_group['strn'])
            for (midi, e11i, e22i, nu12i, g12i, g1zi, g2zi, rhoi, a1i, a2i, trefi,
                 xti, xci, yti, yci, si, gei, f12i, strni) in zip(
                     mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, tref, xt, xc, yt, yc, s, ge, f12, strn):
                model.add_mat8(midi, e11i, e22i, nu12i, g12=g12i, g1z=g1zi, g2z=g2zi, rho=rhoi,
                               a1=a1i, a2=a2i, tref=trefi, Xt=xti, Xc=xci, Yt=yti, Yc=yci,
                               S=si, ge=gei, F12=f12i, strn=strni, comment='')
        elif card_type == 'MAT9':
            ## TODO: add G
            mid = _cast(sub_group['mid'])
            a = _cast(sub_group['A'])
            tref = _cast(sub_group['tref'])
            ge = _cast(sub_group['ge'])
            rho = _cast(sub_group['rho'])
            for midi, ai, trefi, gei, rhoi in zip(mid, a, tref, ge, rho):
                model.add_mat9(
                    midi,
                    G11=0., G12=0., G13=0., G14=0., G15=0., G16=0.,
                    G22=0., G23=0., G24=0., G25=0., G26=0.,
                    G33=0., G34=0., G35=0., G36=0.,
                    G44=0., G45=0., G46=0.,
                    G55=0., G56=0.,
                    G66=0.,
                    rho=rhoi, A=ai, tref=trefi, ge=gei, comment='')

        elif card_type in ['MAT8', 'MAT9']:
            model.log.warning('skipping materials/%s because its vectorized '
                              'and needs a loader' % card_type)
        else:
            #model.add_mat4(mid, k, cp=0.0, rho=1.0, H=None, mu=None, hgen=1.0,
                           #ref_enthalpy=None, tch=None, tdelta=None, qlat=None, comment='')
            #model.add_mat5(mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0., kzz=0.,
                           #cp=0., rho=1., hgen=1., comment='')
            #model.add_mat10(mid, bulk, rho, c, ge=0.0, gamma=None,
                            #table_bulk=None, table_rho=None, table_ge=None,
                            #table_gamma=None, comment='')
            #model.add_mat11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
                            #rho=0.0, a1=0.0, a2=0.0, a3=0.0, tref=0.0, ge=0.0, comment='')
            load_cards_from_keys_values('materials/%s' % card_type, sub_group)

def hdf5_load_spcs(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for spc_id in keys:
        cards_group = group[spc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'SPC1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('spcs/%s/%s' % (spc_id, card_type), sub_group)

def hdf5_load_spcadds(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for spc_id in keys:
        cards_group = group[spc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'SPC1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('spcadds/%s/%s' % (spc_id, card_type), sub_group)

def hdf5_load_mpcs(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #mpc_ids = _cast(group['keys'])
    for mpc_id in keys:
        cards_group = group[mpc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MPC':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('mpcs/%s/%s' % (mpc_id, card_type), sub_group)

def hdf5_load_mpcadds(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for mpc_id in keys:
        cards_group = group[mpc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MPCADD':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('mpcadds/%s/%s' % (mpc_id, card_type), sub_group)

def hdf5_load_loads(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    spc_ids = _cast(group['keys'])
    for load_id in keys:
        cards_group = group[load_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            if card_type in ['FORCE', 'MOMENT']:
                if card_type == 'FORCE':
                    func = model.add_force
                else:
                    func = model.add_moment
                sid = _cast(sub_group['sid'])
                node = _cast(sub_group['node'])
                cid = _cast(sub_group['cid'])
                mag = _cast(sub_group['mag'])
                xyz = _cast(sub_group['xyz'])
                for (sidi, nodei, magi, xyzi, cidi) in zip(sid, node, mag, xyz, cid):
                    func(sidi, nodei, magi, xyzi, cid=cidi, comment='')
            else:
                #model.add_force1(sid, node, mag, g1, g2, comment='')
                load_cards_from_keys_values('loads/%s/%s' % (load_id, card_type), sub_group)

def hdf5_load_load_combinations(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for load_id in keys:
        cards_group = group[load_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'LOAD':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('load_combinations/%s/%s' % (load_id, card_type), sub_group)

def hdf5_load_nsms(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for nsm_id in keys:
        cards_group = group[nsm_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'NSM':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('nsms/%s/%s' % (nsm_id, card_type), sub_group)

def hdf5_load_nsmadds(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for nsm_id in keys:
        cards_group = group[nsm_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'NSMADD':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('nsmadds/%s/%s' % (nsm_id, card_type), sub_group)

def hdf5_load_frequencies(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for freq_id in keys:
        cards_group = group[freq_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'FREQ':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('frequencies/%s/%s' % (freq_id, card_type), sub_group)

def hdf5_load_dloads(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #dload_ids = _cast(group['keys'])
    for dload_id in keys:
        cards_group = group[dload_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'DLOAD':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('dloads/%s/%s' % (dload_id, card_type), sub_group)

def hdf5_load_dload_entries(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #dload_ids = _cast(group['keys'])
    for dload_id in keys:
        cards_group = group[dload_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'TLOAD1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('dload_entries/%s/%s' % (dload_id, card_type), sub_group)

def hdf5_load_bcs(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #dload_ids = _cast(group['keys'])
    for bc_id in keys:
        cards_group = group[bc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('bcs/%s/%s' % (bc_id, card_type), sub_group)

def hdf5_load_transfer_functions(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #dload_ids = _cast(group['keys'])
    for bc_id in keys:
        cards_group = group[bc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('transfer_functions/%s/%s' % (bc_id, card_type), sub_group)

def hdf5_load_dvgrids(model, group, encoding):
    keys = list(group.keys())
    keys.remove('keys')
    #dload_ids = _cast(group['keys'])
    for opt_id in keys:
        cards_group = group[opt_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            load_cards_from_keys_values('dvgrids/%s/%s' % (opt_id, card_type), sub_group)

def hdf5_load_generic(unused_model, group, name, encoding):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        load_cards_from_keys_values('%s/%s' % (name, card_type), sub_group)



def hdf5_load_properties(model, properties_group, encoding):
    """loads the properties from an HDF5 file"""
    for card_type in properties_group.keys():
        properties = properties_group[card_type]
        if card_type == 'PSHELL':
            pid = _cast(properties['pid'])
            mids = _cast(properties['mids'])
            z = _cast(properties['z'])
            t = _cast(properties['t'])
            twelveIt3 = _cast(properties['twelveIt3'])
            tst = _cast(properties['tst'])
            nsm = _cast(properties['nsm'])
            for pidi, (mid1, mid2, mid3, mid4), (z1, z2), ti, twelveIt3i, tsti, nsmi in zip(
                    pid, mids, z, t, twelveIt3, tst, nsm):
                if np.isnan(ti):
                    ti = None
                    raise RuntimeError('Differential shell thickness is not supported')
                if np.isnan(z1):
                    z1 = None
                if np.isnan(z2):
                    z2 = None
                model.add_pshell(pidi, mid1=mid1, t=ti, mid2=mid2, twelveIt3=twelveIt3i,
                                 mid3=mid3, tst=tsti, nsm=nsmi, z1=z1, z2=z2, mid4=mid4,
                                 comment='')
        elif card_type in ['PSOLID', 'PIHEX']:
            func = model.add_psolid if card_type == 'PSOLID' else model.add_pihex
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            cordm = _cast(properties['cordm'])
            integ = _cast(properties['integ'])
            isop = _cast(properties['isop'])
            stress = _cast(properties['stress'])
            fctn = _cast(properties['fctn'])
            for pidi, midi, cordmi, integi, stressi, isopi, fctni in zip(
                    pid, mid, cordm, integ, stress, isop, fctn):
                func(pidi, midi, cordm=cordmi, integ=integi, stress=stressi,
                     isop=isopi, fctn=fctni, comment='')

        elif card_type == 'PROD':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            A = _cast(properties['A'])
            j = _cast(properties['J'])
            c = _cast(properties['c'])
            nsm = _cast(properties['nsm'])
            for pidi, midi, Ai, ji, ci, nsmi in zip(
                    pid, mid, A, j, c, nsm):
                model.add_prod(pidi, midi, Ai, j=ji, c=ci, nsm=nsm, comment='')

        elif card_type == 'PTUBE':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            OD = _cast(properties['OD'])
            t = _cast(properties['t'])
            nsm = _cast(properties['nsm'])
            for pidi, midi, (OD1, OD2), ti, nsmi in zip(
                    pid, mid, OD, t, nsm):
                model.add_ptube(pidi, midi, OD1, t=ti, nsm=nsmi, OD2=OD2, comment='')

        elif card_type == 'PBAR':
            pid = _cast(properties['pid'])
            mid = _cast(properties['mid'])
            A = _cast(properties['A'])
            J = _cast(properties['J'])
            I = _cast(properties['I'])

            c = _cast(properties['c'])
            d = _cast(properties['d'])
            e = _cast(properties['e'])
            f = _cast(properties['f'])
            k = _cast(properties['k'])

            nsm = _cast(properties['nsm'])
            for (pidi, midi, Ai, Ji, (i1, i2, i12),
                 (c1, c2), (d1, d2), (e1, e2), (f1, f2), (k1, k2), nsmi) in zip(
                     pid, mid, A, J, I,
                     c, d, e, f, k, nsm):
                if k1 == np.nan:
                    k1 = None
                if k2 == np.nan:
                    k2 = None
                model.add_pbar(pidi, midi, A=Ai, i1=i1, i2=i2, i12=i12, j=Ji, nsm=nsmi,
                               c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                               f1=f1, f2=f2, k1=k1, k2=k2, comment='')

        else:
            load_cards_from_keys_values('properties/%s' % card_type, properties)
            #model.add_pshear(pid, mid, t, nsm=0., f1=0., f2=0., comment='')
            #model.add_pvisc(pid, ce, cr, comment='')
            #model.add_pelas(pid, k, ge=0., s=0., comment='')
            #model.add_pdamp(pid, b, comment='')
            #model.add_pcomp(pid, mids, thicknesses, thetas=None, souts=None, nsm=0., sb=0.,
                            #ft=None, tref=0., ge=0., lam=None, z0=None, comment='')
            #model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                             #nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None,
                             #comment='')


def load_cards_from_keys_values(name, properties):
    try:
        keys = _cast(properties['keys'])
    except KeyError:  # pragma: no cover
        print('name = %s' % name)
        print(properties)
        raise
    #except TypeError:  # pragma: no cover
        #print('name = %s' % name)
        #print(properties)
        #print(properties['keys'])
        #raise
    values = properties['values']
    _load_cards_from_keys_values(name, values, keys)

def _load_cards_from_keys_values(name, values, keys):
    for key, keyi in zip(keys, values.keys()):
        value = values[keyi]
        keys_to_read = list(value.keys())
        card_type = _cast(value['type'])
        class_obj = CARD_MAP[card_type]
        if hasattr(class_obj, '_init_from_empty'):
            class_instance = class_obj._init_from_empty()
        else:
            try:
                class_instance = class_obj()
            except TypeError:  # pragma: no cover
                print('error loading %r' % card_type)
                print(class_obj)
                raise

        _properties = []
        if hasattr(class_obj, '_properties'):
            _properties = class_obj._properties
        for key_to_cast in keys_to_read:
            if key_to_cast in _properties:
                continue

            try:
                valuei = _cast(value[key_to_cast])
            except AttributeError:
                valuei = None

            try:
                setattr(class_instance, key_to_cast, valuei)
            except AttributeError:  # pragma: no cover
                print('error loading %r' % card_type)
                print(_properties)
                print(key, key_to_cast, valuei)
                raise

def _load_from_class(value, card_type):
    keys_to_read = list(value.keys())
    class_obj = CARD_MAP[card_type]
    if hasattr(class_obj, '_init_from_empty'):
        class_instance = class_obj._init_from_empty()
    else:
        try:
            class_instance = class_obj()
        except TypeError:  # pragma: no cover
            print('error loading %r' % card_type)
            print(class_obj)
            raise

    _properties = []
    if hasattr(class_obj, '_properties'):
        _properties = class_obj._properties

    for key_to_cast in keys_to_read:
        if key_to_cast in _properties:
            continue
        try:
            valuei = _cast(value[key_to_cast])
        except AttributeError:
            valuei = None

        try:
            setattr(class_instance, key_to_cast, valuei)
        except AttributeError:  # pragma: no cover
            print('error loading %r' % card_type)
            print(_properties)
            print(key_to_cast, valuei)
            raise
    return class_instance


def hdf5_load_elements(model, elements_group, encoding):
    """loads the elements from an HDF5 file"""
    for card_type in elements_group.keys():
        elements = elements_group[card_type]
        if card_type == 'CTETRA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctetra(eid, pid, nids, comment='')
        elif card_type == 'CPENTA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_chexa(eid, pid, nids, comment='')
        elif card_type == 'CPYRAM':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_cpyram(eid, pid, nids, comment='')
        elif card_type == 'CHEXA':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_chexa(eid, pid, nids, comment='')

        elif card_type == 'CROD':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_crod(eid, pid, nids, comment='')
        elif card_type == 'CTUBE':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctube(eid, pid, nids, comment='')
        elif card_type == 'CONROD':
            eids = _cast(elements['eid'])
            mids = _cast(elements['mid'])
            nodes = _cast(elements['nodes']).tolist()
            A = _cast(elements['A'])
            J = _cast(elements['J'])
            c = _cast(elements['c'])
            nsm = _cast(elements['nsm'])
            for eid, mid, nids, ai, ji, ci, nsmi in zip(eids, mids, nodes, A, J, c, nsm):
                model.add_conrod(eid, mid, nids, A=ai, j=ji, c=ci, nsm=nsmi, comment='')

        elif card_type == 'CBAR':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            offt = _cast(elements['offt'])
            wa = _cast(elements['wa'])
            wb = _cast(elements['wb'])
            pa = _cast(elements['pa'])
            pb = _cast(elements['pb'])
            for eid, pid, nids, xi, g0i, offti, pai, pbi, wai, wbi in zip(
                    eids, pids, nodes, x, g0, offt, pa, pb, wa, wb):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                model.add_cbar(eid, pid, nids, xi, g0i, offt=offti.decode(encoding),
                               pa=pai, pb=pbi, wa=wai, wb=wbi, comment='')

        elif card_type == 'CBEAM':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            bit = _cast(elements['bit'])
            offt = _cast(elements['offt'])
            sa = _cast(elements['sa'])
            sb = _cast(elements['sb'])
            wa = _cast(elements['wa'])
            wb = _cast(elements['wb'])
            pa = _cast(elements['pa'])
            pb = _cast(elements['pb'])
            for eid, pid, nids, xi, g0i, offti, biti, pai, pbi, wai, wbi, sai, sbi in zip(
                    eids, pids, nodes, x, g0, offt, bit, pa, pb, wa, wb, sa, sb):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                if biti == np.nan:
                    offti = offti.decode(encoding)
                else:
                    offti = None
                model.add_cbeam(eid, pid, nids, xi, g0i, offt=offti, bit=biti,
                                pa=pai, pb=pbi, wa=wai, wb=wbi, sa=sai, sb=sbi, comment='')

        elif card_type in ['CELAS1', 'CDAMP1']:
            func = model.add_celas1 if card_type == 'CELAS1' else model.add_cdamp1
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, pid, nids, (c1, c2) in zip(eids, pids, nodes, components):
                func(eid, pid, nids, c1=c1, c2=c2, comment='')
        elif card_type == 'CELAS2':
            eids = _cast(elements['eid'])
            k = _cast(elements['K'])
            ge = _cast(elements['ge'])
            s = _cast(elements['s'])
            nodes = _cast(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, ki, nids, (c1, c2), gei, si in zip(eids, k, nodes, components, ge, s):
                model.add_celas2(eid, ki, nids, c1=c1, c2=c2, ge=gei, s=si, comment='')
        elif card_type == 'CDAMP2':
            eids = _cast(elements['eid'])
            b = _cast(elements['B'])
            nodes = _cast(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, bi, nids, (c1, c2) in zip(eids, b, nodes, components):
                model.add_cdamp2(eid, bi, nids, c1=c1, c2=c2, comment='')

        elif card_type in ['CELAS3', 'CDAMP3', 'CDAMP5', 'CVISC']:
            if card_type == 'CELAS3':
                func = model.add_celas3
            elif card_type == 'CDAMP3':
                func = model.add_cdamp3
            elif card_type == 'CDAMP5':
                func = model.add_cdamp5
            elif card_type == 'CVISC':
                func = model.add_cvisc
            else:
                raise NotImplementedError(card_type)
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_celas3(eid, pid, nids, comment='')
        elif card_type == 'CELAS4':
            eids = _cast(elements['eid'])
            k = _cast(elements['K'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, ki, nids in zip(eids, k, nodes):
                model.add_celas4(eid, ki, nids, comment='')
        elif card_type == 'CDAMP4':
            eids = _cast(elements['eid'])
            b = _cast(elements['B'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, bi, nids in zip(eids, b, nodes):
                model.add_cdamp4(eid, bi, nids, comment='')


        elif card_type == 'CBUSH':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            cid = _cast(elements['cid'])
            ocid = _cast(elements['ocid'])
            s = _cast(elements['s'])
            si = _cast(elements['si'])
            for eid, pid, nids, xi, g0i, cidi, s2, ocidi, si2 in zip(
                    eids, pids, nodes, x, g0, cid, s, ocid, si):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                if cidi == -1:
                    cidi = None

                if si2[0] == np.nan:
                    si2 = [None, None, None]
                model.add_cbush(eid, pid, nids, x, g0, cid=cidi, s=s2, ocid=ocidi, si=si2,
                                comment='')

        elif card_type == 'CGAP':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            g0 = _cast(elements['g0'])
            x = _cast(elements['x'])
            cid = _cast(elements['cid'])
            for eid, pid, nids, xi, g0i, cidi in zip(
                    eids, pids, nodes, x, g0, cid):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                if cidi == -1:
                    cidi = None
                model.add_cgap(eid, pid, nids, xi, g0i, cid=cidi, comment='')

        elif card_type == 'CBUSH1D':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            cid = _cast(elements['cid'])
            for eid, pid, nids, cidi in zip(eids, pids, nodes, cid):
                if cidi == -1:
                    cidi = None
                model.add_cbush1d(eid, pid, nids, cid=cidi, comment='')

        elif card_type in ['CTRIA3', 'CTRIAR']:
            func = model.add_ctria3 if card_type == 'CTRIA3' else model.add_ctriar
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_ctria3(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, comment='')
        elif card_type in ['CQUAD4', 'CQUADR']:
            func = model.add_cquad4 if card_type == 'CQUAD4' else model.add_cquadr
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                func(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                     tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')

        elif card_type == 'CTRIA6':
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_ctria6(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, comment='')
        elif card_type == 'CQUAD8':
            # TODO: doesn't support tflag
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            zoffsets = _cast(elements['zoffset'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_cquad8(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')

        elif card_type == 'CQUAD':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta in zip(eids, pids, nodes, mcids, thetas):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_cquad(eid, pid, nids, theta_mcid=theta_mcid, comment='')

        elif card_type == 'CSHEAR':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_cshear(eid, pid, nids, comment='')

        elif card_type == 'CTRIAX':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta in zip(eids, pids, nodes, mcids, thetas):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_ctriax(eid, pid, nids, theta_mcid=theta_mcid, comment='')
        elif card_type == 'CTRAX3':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                model.add_ctrax3(eid, pid, nids, theta=theta, comment='')
        elif card_type == 'CTRAX6':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                model.add_ctrax6(eid, pid, nids, theta=theta, comment='')
        elif card_type == 'CTRIAX6':
            eids = _cast(elements['eid'])
            mids = _cast(elements['mid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, mid, nids, theta in zip(eids, mids, nodes, thetas):
                model.add_ctriax6(eid, mid, nids, theta=theta, comment='')

        elif card_type == 'CQUADX':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            mcids = _cast(elements['mcid'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta, mcid in zip(eids, pids, nodes, thetas, mcids):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                nids = [None if nid == 0 else nid
                        for nid in nids]
                model.add_cquadx(eid, pid, nids, theta_mcid=theta_mcid, comment='')

        elif card_type == 'CQUADX4':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                model.add_cquadx4(eid, pid, nids, theta=theta, comment='')
        elif card_type == 'CQUADX8':
            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                model.add_cquadx8(eid, pid, nids, theta=theta, comment='')

        elif card_type in ['CPLSTN3', 'CPLSTN4']:
            func = model.add_cplstn3 if card_type == 'CPLSTN3' else model.add_cplstn4

            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')
        elif card_type in ['CPLSTN6', 'CPLSTN8']:
            func = model.add_cplstn6 if card_type == 'CPLSTN6' else model.add_cplstn8

            eids = _cast(elements['eid'])
            pids = _cast(elements['pid'])
            thetas = _cast(elements['theta'])
            nodes = _cast(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')
        else:
            load_cards_from_keys_values('elements/%s' % card_type, elements)
            #model.add_cdamp4(eid, b, nids, comment='')
            #model.add_cbush2d(eid, pid, nids, cid=0, plane='XY', sptid=None, comment='')
            #model.add_cfast(eid, pid, Type, ida, idb, gs=None, ga=None, gb=None,
                            #xs=None, ys=None, zs=None, comment='')
            #model.add_cmass1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cmass2(eid, mass, nids, c1, c2, comment='')
            #model.add_cmass3(eid, pid, nids, comment='')
            #model.add_cmass4(eid, mass, nids, comment='')
            #model.log.debug(card_type)
