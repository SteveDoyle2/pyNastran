"""Defines various helper functions for exporting a HDF5 BDF file"""
from __future__ import annotations
from collections import defaultdict
from typing import List, Any, TYPE_CHECKING
from io import StringIO
import numpy as np

from pyNastran.utils.dict_to_h5py import (
    _add_list_tuple, integer_types, float_types, string_types)
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP
from pyNastran.utils import object_attributes

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

# dict[key] : [value1, value2, ...]
dict_int_list_obj_attrs = [
    'spcs', 'spcadds',
    'mpcs', 'mpcadds',
    'loads', 'load_combinations',
    'dloads', 'dload_entries',
    #'usets', # has string keys
    'nsms', 'nsmadds',
    'frequencies',
    'bcs', 'transfer_functions',
    'dvgrids',

    # parametric
    'pval',
]

# dict[key] : value
dict_int_obj_attrs = [
    # are handled explictly----
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
    #'plotels',
    'MATS1', 'MATS3', 'MATS8', 'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',

    # TODO: don't work
    #'reject_count',
    'dresps',

    # optimization
    'dconadds', 'ddvals', 'dequations', 'dscreen',
    'dvcrels', 'dvmrels', 'dvprels',
    # aero
    'aecomps', 'aefacts', 'aelists', 'aeparams',
    'aestats', 'aesurf', 'aesurfs',
    'divergs', 'dlinks',
    'flfacts', 'paeros',

    # contact
    'bconp', 'bcrparas', 'bctadds',
    'bctparas', 'bctsets', 'blseg', 'bsurf', 'bsurfs',
    'bfric',

    # other
    'ao_element_flags', 'cMethods',
    'convection_properties',
    'dareas',
    'dmig', 'dmiji', 'dmij', 'dmik', 'dmi', 'dmiax',
    'dti',
    'dphases', 'delays',
    'epoints', 'gridb',
    'nlparms', 'nlpcis',
    'normals',
    'nxstrats',
    'pbusht', 'pdampt', 'pelast', 'phbdys', 'points',
    'properties_mass',
    'radcavs', 'radmtx', 'random_tables',
    'rotors',
    'sets',
    'spcoffs',
    'spoints',
    'suport1',
    'tables_d', 'tables_m', 'tables_sdamping', 'tempds',
    'tics',
    'tstepnls', 'tsteps',
    'view3ds', 'views',

    # superelements
    'csuper', 'csupext',
    'se_sets', 'se_usets', 'sebndry', 'sebulk', 'seconct', 'seelt',
    'seexcld', 'selabel', 'seload', 'seloc', 'sempln', 'senqset', 'setree',
    'release',

    # axisymmetric
    'ringaxs', 'ringfl',

    # parametric
    'pset', 'gmcurv', 'feedge', 'feface', 'gmsurf',

    # cyclic
    'cyjoin',
]

scalar_obj_keys = [
    # required----
    'aero', 'aeros', 'axic', 'axif', 'cyax', 'baror', 'beamor',
    'acmodl', 'modtrak',
    'doptprm',
    'dtable', 'grdset', 'radset', 'seqgp',
    'case_control_deck',
    #'zona',
]

scalar_keys = [
    # handled separately----
    #'cards_to_read',

    # basic types----
    'bdf_filename',
    '_auto_reject', '_encoding', '_iparse_errors', '_is_axis_symmetric', '_is_cards_dict',
    '_is_dynamic_syntax', '_is_long_ids', '_ixref_errors', '_nastran_format', '_nparse_errors',
    '_nxref_errors', '_sol', '_stop_on_duplicate_error', '_stop_on_parsing_error',
    '_stop_on_xref_error',
    '_xref', 'active_filename',
    'dumplines', 'echo', 'force_echo_off', 'include_dir',
    'is_msc', 'is_nx', 'punch',
    'read_includes', 'save_file_structure', 'sol', 'sol_iline', 'nastran_format',
    'is_superelements', 'is_zona', 'sol_method', 'debug',
    #'_unique_bulk_data_cards',
    #'is_bdf_vectorized',
    #'is_long_ids',


    # not sure----
    #'nid_cp_cd', 'xyz_cid0',

    # removed-----
    #'ncaeros', 'ncoords', 'nnodes', 'nelements', 'nproperties', 'nmaterials',
    #'point_ids', 'wtmass', 'log',
]

LIST_KEYS = [
    # handled in minor_attributes----
    #'active_filenames', 'executive_control_lines', 'case_control_lines',
    #'system_command_lines',
    #'reject_cards',

    # required------------------
    #'initial_superelement_models',

    # maybe...
    #'_duplicate_coords', '_duplicate_elements', '_duplicate_masses', '_duplicate_materials',
    '_duplicate_nodes', '_duplicate_properties',
    '_duplicate_thermal_materials', '_stored_parse_errors',
    #'_stored_xref_errors',
    #'units', 'xyz_limits',

    # removed
    #'coord_ids', 'element_ids', 'material_ids', 'node_ids', 'property_ids', 'caero_ids',
    #'special_cards',
]

LIST_OBJ_KEYS = [  ## TODO: not done
    # TODO: required
    'asets', 'bsets', 'csets', 'omits', 'qsets',
    'mkaeros',
    'monitor_points',
    'suport',
    'se_bsets', 'se_csets', 'se_qsets', 'se_suport',
]


def export_bdf_to_hdf5_file(hdf5_file, model: BDF, exporter=None):
    """
    Converts the BDF objects into hdf5 object

    Parameters
    ----------
    hdf5_file : H5File()
        an h5py object
    exporter : HDF5Exporter; default=None
        unused

    """
    unused_attrs = object_attributes(model, mode='both', keys_to_skip=None,
                                     filter_properties=True)
    encoding = model.get_encoding()

    if 'GRID' in model.card_count:
        model.log.debug('exporting nodes')
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
        'properties', 'masses', 'rigid_elements', 'plotels',

        # materials
        'materials', 'thermal_materials', 'creep_materials', 'hyperelastic_materials',
        #'MATS1',
        #'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',

        # aero
        'caeros', 'splines', 'flutters', 'trims', 'csschds', 'gusts',

        # other
        'methods', 'tables', 'desvars', 'topvar',
    ]
    for group_name in groups_to_export:
        _hdf5_export_group(hdf5_file, model, group_name, encoding)


    unused_dict_int_attrs = [  # TODO: not used...
        # required
        '_dmig_temp',
        'include_filenames',
        'superelement_models',
        'values_to_skip',

        # removed
        #'rsolmap_to_str',
        #'nid_map',
        #'subcases',
    ]

    _export_dict_int_obj_attrs(model, hdf5_file, encoding)
    _export_dict_int_list_obj_attrs(model, hdf5_file, encoding)
    _export_dconstrs(hdf5_file, model, encoding)

    #for key in scalar_obj_keys:
        #value = getattr(model, key)
        #hdf5_file.create_dataset(key, value)
    if model.params:
        model.log.debug('exporting params')
        skip_attrs = ['comment', '_field_map']
        group = hdf5_file.create_group('params')
        for key, param in model.params.items():
            _h5_export_class(group, model, key, param, skip_attrs, encoding, debug=False)

    if model.aelinks:
        model.log.debug('exporting aelinks')
        skip_attrs = ['comment', '_field_map']
        group = hdf5_file.create_group('aelinks')
        for aelink_id, aelinks in model.aelinks.items():
            groupi = group.create_group(str(aelink_id))
            for j, aelinki in enumerate(aelinks):
                key = str(j)
                _h5_export_class(groupi, model, key, aelinki, skip_attrs, encoding, debug=False)

    if model.usets:
        model.log.debug('exporting usets')
        skip_attrs = ['comment', '_field_map']
        group = hdf5_file.create_group('usets')
        for name, usets in model.usets.items():
            groupi = group.create_group(name)
            #print(usets)
            for i, uset in enumerate(usets):
                #print(uset.get_stats())
                key = str(i)
                _h5_export_class(groupi, model, key, uset, skip_attrs, encoding, debug=False)

    _export_scalar_group(hdf5_file, model, encoding)

    skip_attrs = ['comment', '_field_map']
    for key in scalar_obj_keys:
        value = getattr(model, key)
        if value is None:
            #print('None: %s %s' % (key, value))
            pass
        else:
            model.log.debug('exporting %s' % key)
            _h5_export_class(hdf5_file, model, key, value, skip_attrs, encoding, debug=False)

    _export_list_keys(model, hdf5_file, LIST_KEYS)
    _export_list_obj_keys(model, hdf5_file, LIST_OBJ_KEYS, encoding)

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


def _export_dconstrs(hdf5_file, model: BDF, encoding):
    """exports the dconstrs, which includes DCONSTRs and DCONADDs"""
    if model.dconstrs:
        dconstrs_group = hdf5_file.create_group('dconstrs')
        unused_keys = list(model.dconstrs.keys())

        dconstrsi = []
        dconadds = []
        for unused_key, dconstrs in model.dconstrs.items():
            #group = dconstrs_group.create_group(str(key))
            for dconstr in dconstrs:
                Type = dconstr.type
                if Type == 'DCONSTR':
                    dconstrsi.append(dconstr)
                elif Type == 'DCONADD':
                    dconadds.append(dconstr)

        ndconstrs = len(dconstrsi)
        if ndconstrs:
            dconstr_group = dconstrs_group.create_group('DCONSTR')
            unused_keys = np.arange(ndconstrs, dtype='int32')
            dconstr0 = dconstrsi[0]
            dconstr0.export_to_hdf5(dconstr_group, dconstrsi, encoding)
            ndconstrs = len(dconstrsi)

        ndconadds = len(dconadds)
        if ndconadds:
            dconadd_group = dconstrs_group.create_group('DCONADD')
            unused_keys = np.arange(ndconadds, dtype='int32')
            dconadds0 = dconadds[0]
            dconadds0.export_to_hdf5(dconadd_group, dconadds, encoding)

def _export_scalar_group(hdf5_file, model: BDF, encoding):
    scalar_group = _export_minor_attributes(hdf5_file, model, encoding)

    for key in ['case_control_lines', 'executive_control_lines',
                'system_command_lines', 'active_filenames']:
        # these are exported to the minor_attributes group
        list_str = getattr(model, key)
        if not len(list_str):
            continue
        list_bytes = [line.encode(encoding) for line in list_str]
        #print(scalar_group)
        scalar_group.create_dataset(key, data=list_bytes)

    if len(model.reject_lines):
        reject_group = scalar_group.create_group('reject_lines')
        for i, list_str in enumerate(model.reject_lines):
            list_bytes = [line.encode(encoding) for line in list_str]
            reject_group.create_dataset(str(i), data=list_bytes)

    if len(model.reject_cards):
        reject_group = scalar_group.create_group('reject_cards')
        for i, reject_card in enumerate(model.reject_cards):
            fields = reject_card.fields()
            list_bytes = [field.encode(encoding) if field is not None else b''
                          for field in fields]
            reject_group.create_dataset(str(i), data=list_bytes)

def _export_minor_attributes(hdf5_file, model: BDF, encoding):
    """
    Minor atributes include:
     - encoding
     - include_dir
     - is_enddata
     - is_msc
     - is_nx
     - punch
     - read_includes
     - active_filenames
     - bdf_filename
     - executive_control_lines
     - case_control_lines
     - sol
     - sol_iline
    """
    scalar_group = None
    scalars_keys_to_analyze = []
    for key in scalar_keys:
        if hasattr(model, key):
            scalars_keys_to_analyze.append(key)

    if scalars_keys_to_analyze:
        scalar_group = hdf5_file.create_group('minor_attributes')
        scalar_group.create_dataset('encoding', data=encoding)
        for key in sorted(scalars_keys_to_analyze):
            value = getattr(model, key)
            #model.log.debug('  minor %s' % key)

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
        #del scalar_group

    scalar_group.create_dataset('is_enddata', data='ENDDATA' in model.card_count)
    return scalar_group

def _export_dict_int_obj_attrs(model: BDF, hdf5_file, encoding):
    unused_cards = set(list(CARD_MAP.keys()))
    for attr in dict_int_obj_attrs:
        dict_obj = getattr(model, attr)
        #print(attr, dict_obj)
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

        model.log.debug('exporting %s' % attr)
        try:
            group = hdf5_file.create_group(attr) # 'spcs'
        except ValueError:  # pragma: no cover
            model.log.error('cant create %r' % attr)
            raise

        keys = list(dict_obj.keys())
        keys.sort()
        #model.log.debug('keys = %s' % keys)
        if attr in ['dmig', 'dmij', 'dmi', 'dmik', 'dmiji', 'dmiax']:
            #print('keys =', keys)
            key0 = keys[0]
            value = dict_obj[key0]
            group.attrs['type'] = value.type
            #print('setting type', group, value.type)
            model.log.debug('type = %s' % value.type)
            model.log.debug('export 364')
            value.export_to_hdf5(group, model, encoding)
            return

        group.create_dataset('keys', data=keys)
        for spc_id, spcs_obj in sorted(dict_obj.items()):
            id_group = group.create_group(str(spc_id))

            card_types = defaultdict(list)
            for spc in spcs_obj:
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

def _export_list_keys(model: BDF, hdf5_file, list_keys):
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


def _export_list_obj_keys(model: BDF, hdf5_file, list_obj_keys, encoding):
    for attr in list_obj_keys:
        list_obj = getattr(model, attr) # active_filenames
        if not len(list_obj):
            #model.log.debug('skipping list_key: %s' % attr)
            continue
        #model.log.debug('exporting %s' % attr)

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


def _h5_export_class(sub_group: Any, model: BDF, key: str, value: Any,
                     skip_attrs: List[str], encoding: str, debug: bool=True) -> None:
    #model.log.debug('exporting %s to hdf5' % key)
    #sub_groupi = sub_group.create_group('values')
    class_group = sub_group.create_group(str(key))
    try:
        class_group.attrs['type'] = value.type
    except:  # pragma: no cover
        print('key = %r' % key)
        print('value', value)
        model.log.error('ERROR: key=%s value=%s' % (key, value))
        raise

    #if hasattr(value, 'get_h5attrs'):
        #getattrs

    if hasattr(value, 'export_to_hdf5'):
        #print('value =', value, type(value))
        #print('class_group', class_group)
        #model.log.debug('  export 477 - %s' % class_group)
        value.export_to_hdf5(class_group, model, encoding)
        return
    elif hasattr(value, 'object_attributes'):
        keys_to_skip = []
        if hasattr(value, '_properties'):
            keys_to_skip = value._properties
        h5attrs = value.object_attributes(mode='both', keys_to_skip=keys_to_skip,
                                          filter_properties=True)
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

        #model.log.info('%s %s %s' % (key, h5attr, class_value))
        if debug:
            model.log.info('%s %s %s' % (key, h5attr, class_value))

        if isinstance(class_value, dict):
            class_group.attrs['type'] = 'dict'
            param_group = class_group.create_group(h5attr)
            keysi = []
            valuesi = []
            for i, (keyi, valuei) in enumerate(class_value.items()):
                keysi.append(keyi)
                valuesi.append(valuei)
                #if isinstance(valuei, str):
                    #param_group.create_dataset(str(i), data=valuei.encode('ascii'))
                #elif valuei is None:
                    #param_group.create_dataset(str(i), data=np.nan)
                #else:
                    #param_group.create_dataset(str(i), data=valuei)
            model.log.debug('  exporting dict as keys/values for %s (%s)' % (h5attr, value.type))
            _export_list(param_group, '%s/%s/%s' % (value.type, key, h5attr),
                         'keys', keysi, encoding)
            _export_list(param_group, '%s/%s/%s' % (value.type, key, h5attr),
                         'values', valuesi, encoding)
            continue
        elif isinstance(class_value, (list, np.ndarray)):
            if len(class_value) == 0: # empty list
                continue

            is_nones = []
            for class_valuei in class_value:
                is_none = False
                if class_valuei is None:  # PAERO2 : lth
                    is_none = True
                    #break
                is_nones.append(is_none)
        elif isinstance(class_value, (integer_types, float_types, string_types, bool)):
            is_nones = [False]
        #elif isinstance(class_value, dict) and len(class_value) == 0:
            #pass
        else:
            raise NotImplementedError('%s %s; class_value=%s type=%s' % (
                getattr(value, 'type'), key, class_value, type(class_value)))

        is_none = any(is_nones)

        if all(is_nones):
            model.log.warning('skipping %s attribute: %s %s %s' % (
                value.type, key, h5attr, class_value))
        elif all([not is_nonei for is_nonei in is_nones]): # no Nones
            # no Nones
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
                    if isinstance(valuei, str):
                        param_group.create_dataset(str(i), data=valuei.encode('ascii'))
                    else:
                        param_group.create_dataset(str(i), data=valuei)
                #if isinstance(class_value, list):
                    #print('type(value[0] =', class_value, type(class_value[0]))
                    #raise
        else:
            # mixed Nones and values
            class_group.attrs['type'] = 'list'
            param_group = class_group.create_group(h5attr)
            for i, valuei in enumerate(class_value):
                if isinstance(valuei, str):
                    param_group.create_dataset(str(i), data=valuei.encode('ascii'))
                elif valuei is None:
                    param_group.create_dataset(str(i), data=np.nan)
                else:
                    param_group.create_dataset(str(i), data=valuei)

            #raise RuntimeError('mixed %s attribute: %s %s %s' % (
                #value.type, key, h5attr, class_value))

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

#def _export_lists(h5_group, attr, name, values, encoding):
    #print(name, attr, values)

def _export_list(h5_group, attr, name, values, encoding):
    """
    exports a list of:
     - constant type to a dataset
     - variable type to a numbered list

    """
    values2 = [value.encode(encoding) if isinstance(value, str) else value
               for value in values]
    types = {type(value) for value in values}
    if len(types) == 1:
        #print('types =', types)
        #if isinstance(values[0], list):
            #return _export_lists(h5_group, attr, name, values, encoding)

        if not isinstance(values[0], (integer_types, float_types, string_types)):
            raise TypeError('not a base type; %s; %s' % (attr, values2))
        try:
            h5_group.create_dataset(name, data=values2)
        except TypeError:  # pragma: no cover
            print(attr, name, values2)
            raise
    else:
        sub_group = h5_group.create_group(name)
        sub_group.attrs['type'] = 'list'
        for i, value in enumerate(values2):
            if value is None:
                sub_group.create_dataset(str(i), data=np.nan)
            else:
                try:
                    sub_group.create_dataset(str(i), data=value)
                except TypeError:  # pragma: no cover
                    print(attr, name, values2, i)
                    raise

        #print('%s2' % name, values2)
        #h5_group.create_dataset(name, data=values2)
        #raise

def _hdf5_export_elements(hdf5_file, model: BDF, encoding):
    """
    exports the elements to an hdf5_file

    TODO: not done
    """
    etypes_actual = []
    if len(model.elements) == 0:
        return
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

        for card_type in sorted(etypes_actual):
            element_group = elements_group.create_group(card_type)
            eids = model._type_to_id_map[card_type]
            #print(card_type, eids)
            if len(eids) == 0:
                continue
            #for eid in eids:
                #elem = model.elements[eid]
                #print(elem)
            class_obj = CARD_MAP[card_type]

            if hasattr(class_obj, 'export_to_hdf5'):
                class_obj.export_to_hdf5(element_group, model, eids)
            else:
                _hdf5_export_object_dict(element_group, model, card_type,
                                         model.elements, eids, encoding)

def _hdf5_export_group(hdf5_file, model: BDF, group_name, encoding, debug=False):
    """
    exports the properties to an hdf5_file
    """
    data_dict = getattr(model, group_name) # model.properties

    #if group_name in ['splines']:
        #debug = True
    if debug:
        model.log.debug('%s %s' % (group_name, data_dict))

    types_actual = []
    types = model._slot_to_type_map[group_name]
    if debug:
        model.log.debug('card_count = %s' % model.card_count)
        model.log.debug('types = %s' % types)

    for card_name in types:
        #PSHELL
        if card_name in model.card_count:
            types_actual.append(card_name)
            continue

    if types_actual:
        model.log.debug('exporting %s' % group_name)
        if debug:  # pragma: no cover
            print('types_actual =', types_actual)

        group = hdf5_file.create_group(group_name)
        for card_type in types_actual:
            sub_group = group.create_group(card_type)
            ids = model._type_to_id_map[card_type]
            if debug:  # pragma: no cover
                print(ids)
            assert len(ids) > 0, '%s : %s' % (card_type, ids)
            class_obj = CARD_MAP[card_type]
            #print(class_obj)
            if hasattr(class_obj, 'export_to_hdf5'):
                class_obj.export_to_hdf5(sub_group, model, ids)
            else:
                _hdf5_export_object_dict(sub_group, model, card_type, data_dict, ids, encoding)
    #else:
        #model.log.debug('skipping %s to hdf5' % group_name)

def _hdf5_export_object_dict(group, model: BDF, name, obj_dict, keys, encoding):
    #i = 0
    skip_attrs = ['comment', '_field_map']

    keys_write = list(keys)
    if name in ['dmig', 'dmij', 'dmik', 'dmiji', 'dmiax', 'dmi', 'dresps']:
        keys = list(keys)
        #print(group)
        key0 = keys_write[0]
        value = obj_dict[key0]
        group.attrs['type'] = value.type
        #print('group...', group)
        model.log.debug('exporting %s' % name)
        value.export_to_hdf5(group, model, encoding)
        return

    if isinstance(keys_write[0], str):
        keys_write = list([key.encode(encoding) if isinstance(key, str) else key
                           for key in list(keys_write)])

    sub_group = group.create_group('values')
    assert isinstance(name, string_types), 'name=%s; type=%s' % (name, type(name))

    for key in keys:
        value = obj_dict[key]
        #if isinstance(value, str):
            #value = value.encode(encoding)

        #try:
        _h5_export_class(sub_group, model, key, value, skip_attrs, encoding, debug=False)
        #except:  # pragma: no cover
            #raise
            # for debugging
            #sub_group2 = group.create_group('values2')
            #_h5_export_class(sub_group2, model, key, value, skip_attrs, encoding, debug=True)
        #i += 1

    #group.attrs['type'] = class_name
    #print('%s keys = %s' % (name, keys))
    try:
        group.create_dataset('keys', data=keys_write)
    except TypeError:  # pragma: no cover
        print('name =', name)
        print('encoding =', encoding)
        print('keys =', keys)
        raise
