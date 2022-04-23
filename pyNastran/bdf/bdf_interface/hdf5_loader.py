"""Defines various helper functions for loading a HDF5 BDF file"""
from __future__ import annotations
from itertools import count
from typing import List, Any, TYPE_CHECKING
import numpy as np
import h5py

from pyNastran.bdf.bdf import DMIAX, MDLPRM
from pyNastran.utils.dict_to_h5py import _cast, _cast_array, cast_string, cast_strings
from pyNastran.bdf.bdf_interface.encoding import decode_lines
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP
from pyNastran.bdf.bdf_interface.hdf5_exporter import (
    dict_int_obj_attrs, scalar_obj_keys, LIST_OBJ_KEYS)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

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

def load_bdf_from_hdf5_file(h5_file, model):
    """
    Loads an h5 file object into an OP2 object

    Parameters
    ----------
    h5_file : H5File()
        an h5py file object
    model : BDF()
        the BDF file to put the data into

    """
    encoding = cast_string(h5_file['minor_attributes']['encoding'], 'latin1')
    model._encoding = encoding
    assert isinstance(encoding, str), f'encoding={encoding!r}; type={type(encoding)}'
    model.get_encoding()
    keys = h5_file.keys()

    mapper = {
        'elements' : hdf5_load_elements,
        'plotels' : hdf5_load_plotels,
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

        'pval' : hdf5_load_pval,

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
        'aelinks' : hdf5_load_aelinks,
        'desvars' : hdf5_load_desvars,

        'dmig' : hdf5_load_dmigs,
        'dmiax' : hdf5_load_dmigs,
        'dmij' : hdf5_load_dmigs,
        'dmik' : hdf5_load_dmigs,
        'dmiji' : hdf5_load_dmigs,
        'dmi' : hdf5_load_dmigs,
        'dti' : hdf5_load_dti,

        'dconstrs' : hdf5_load_dconstrs,
        'dresps' : hdf5_load_dresps,
        'usets' : hdf5_load_usets,
    }
    generic_mapper = {
        'rigid_elements' : hdf5_load_generic,
        'thermal_materials' : hdf5_load_generic,
        'creep_materials' : hdf5_load_generic,
        'hyperelastic_materials' : hdf5_load_generic,

        'flutters' : hdf5_load_generic,
        'trims' : hdf5_load_generic,
        'csschds' : hdf5_load_generic,
        'gusts' : hdf5_load_generic,
        'caeros' : hdf5_load_generic,
        'splines' : hdf5_load_generic,
        'mdlparm' : hdf5_load_generic,
        #'MATS1' : hdf5_load_generic,
        #'MATT1' : hdf5_load_generic,
        #'MATT2' : hdf5_load_generic,
        #'MATT3' : hdf5_load_generic,
        #'MATT4' : hdf5_load_generic,
        #'MATT5' : hdf5_load_generic,
        #'MATT8' : hdf5_load_generic,
        #'MATT9' : hdf5_load_generic,
    }
    #print('keys =', list(keys))
    for key in keys:
        #model.log.debug('loading %s' % key)
        group = h5_file[key]
        if key == 'nodes':
            grids = group['GRID']
            nids = _cast_array(grids['nid'])
            xyz = _cast_array(grids['xyz'])
            cp = _cast_array(grids['cp'])
            cd = _cast_array(grids['cd'])
            ps = _cast(grids['ps'])
            seid = _cast_array(grids['seid'])
            for nid, xyzi, cpi, cdi, psi, seidi in zip(nids, xyz, cp, cd, ps, seid):
                model.add_grid(nid, xyzi, cp=cpi, cd=cdi, ps=psi, seid=seidi, comment='')
            model.card_count['GRID'] = len(nids)

        elif key in mapper:
            func = mapper[key]
            func(model, group, encoding)
        elif key in generic_mapper:
            func = generic_mapper[key]
            func(model, group, key, encoding)
        elif key in dict_int_obj_attrs:
            #model.log.debug('  dict_int_obj')
            dkeys, values = load_cards_from_keys_values(
                key, group, encoding, model.log)
            _put_keys_values_into_dict(model, key, dkeys, values)
            card_type = values[0].type
            model.card_count[card_type] = len(dkeys)

        elif key in ['info', 'matrices'] or key.startswith('Subcase'): # op2
            continue
        elif key in ['cards_to_read']: # handled separately
            continue
        elif key == 'params':
            keys = list(group.keys())
            values = _load_cards_from_keys_values('params', group, keys, encoding, model.log)
            _put_keys_values_into_dict(model, 'params', keys, values, cast_int_keys=False)
            model.card_count['PARAM'] = len(keys)
        elif key == 'bcparas':
            keys = list(group.keys())
            values = _load_cards_from_keys_values('bcparas', group, keys, encoding, model.log)
            _put_keys_values_into_dict(model, 'bcparas', keys, values, cast_int_keys=False)
            model.card_count['BCPARA'] = len(keys)

        elif key == 'minor_attributes':
            _load_minor_attributes(key, group, model, encoding)
        #elif key in ['case_control_lines', 'executive_control_lines', 'system_command_lines']:
            #lst = _load_indexed_list_str(keyi, sub_group, encoding)

        elif key == 'active_filenames':
            if 'value' not in group:
                lst = _load_indexed_list_str(key, group, encoding)
                continue

            lst = _cast_array(group['value']).tolist()
            #else:
            #except KeyError:  # pragma: no cover
                #print('group', group)
                #print('group.keys()', list(group.keys()))
                #raise

            if isinstance(lst[0], str):
                pass
            else:
                lst = [line.encode(encoding) for line in lst]
            setattr(model, key, lst)

        elif key in LIST_OBJ_KEYS:
            #model.log.debug('  list_obj')
            #model.log.info('  key = %s' % key)
            #model.log.info('  group = %s' % group)
            #model.log.info('  group.keys() = %s' % list(group.keys()))
            keys = _cast(group['keys'])
            values = group['values']
            lst = [None] * len(keys)
            for keyi in values.keys():
                ikey = int(keyi)
                class_obj_hdf5 = values[keyi]
                card_type = cast_string(class_obj_hdf5['type'], encoding)
                class_instance = _load_from_class(class_obj_hdf5, card_type, encoding)
                lst[ikey] = class_instance
            _put_keys_values_into_list(model, key, keys, lst)
            #model.log.info('keys = %s' % keys)
            #model.log.info('values = %s' % values)
            #model.log.info('values.keys() = %s' % values.keys())

        elif key in 'case_control_deck':
            lines = []
            model.case_control_deck = CaseControlDeck(lines, log=model.log)
            model.case_control_deck.load_hdf5_file(group, encoding)
            str(model.case_control_deck)

        elif key in scalar_obj_keys: # these only have 1 value
            #model.log.debug('  scalar_obj')
            keys = list(group.keys())
            keys.remove('type')
            card_type = cast_string(group['type'], encoding)
            class_instance = _load_from_class(group, card_type, encoding)
            write_card(class_instance)
            setattr(model, key, class_instance)
            model.card_count[card_type] = 1
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
        elif key in 'mdlprm':
            lines = []
            model.mdlprm = MDLPRM({'HDF5': 1})
            model.mdlprm.load_hdf5_file(group, encoding)
            model.card_count['MDLPRM'] = 1
        else:
            model.log.warning('skipping hdf5 load for %s' % key)
            raise RuntimeError('skipping hdf5 load for %s' % key)

    cards_to_read = _cast(h5_file['cards_to_read'])
    cards_to_read = [key.decode(encoding) for key in cards_to_read]
    model.cards_to_read = set(list(cards_to_read))

def _load_minor_attributes(unused_key: str, group, model: BDF,
                           encoding: str) -> None:
    keys_attrs = group.keys()
    list_attrs = {'case_control_lines', 'executive_control_lines',
                  'system_command_lines', 'active_filenames'}
    str_attrs = {'nastran_format', 'include_dir'}
    #skip_attrs = []
    for keyi in keys_attrs:
        sub_group = group[keyi]
        #model.log.debug('  %s' % keyi)

        if keyi in list_attrs:
            lst = _cast(sub_group)
            if isinstance(lst[0], str):
                pass
            else:
                lst = decode_lines(lst, encoding)
                assert isinstance(lst[0], str), type(lst[0])
            setattr(model, keyi, lst)
            continue
        elif keyi == 'reject_lines':
            reject_keys = list(sub_group.keys())
            lst = [None] * len(reject_keys)

            for reject_key in reject_keys:
                reject_key_int = int(reject_key)
                h5_value = sub_group[reject_key]
                value = _cast(h5_value)
                lst[reject_key_int] = value
                comment = value[0].decode(encoding)
                card_lines = value[1:]
                card_lines = decode_lines(card_lines, encoding)
                try:
                    line0 = card_lines[0]
                except IndexError:
                    # C:\Program Files\Siemens\NX 12.0\NXNASTRAN\nxn12\nast\del\gentim1.dat
                    print(value)
                    print(card_lines)
                    raise
                card_name_field0 = line0.split(',', 1)[0].split('\t', 1)[0]
                card_name = card_name_field0[:8].rstrip().upper().rstrip('*')
                assert isinstance(comment, str), type(comment)

                ## TODO: swap out
                #model.add_card(card_lines, card_name, comment=comment,
                               #ifile=None, is_list=True, has_none=True)
                model.reject_card_lines(card_name, card_lines, comment=comment)
            continue
        elif keyi == 'reject_cards':
            reject_keys = list(sub_group.keys())
            for ireject in sub_group.keys():
                reject_card = _cast(sub_group[ireject])
                if not isinstance(reject_card, list):
                    reject_card = reject_card.tolist()
                fields = decode_lines(reject_card, encoding)
                #fields = [field if field != 'nan' else None for field in fields]
                card_name = fields[0]
                model.add_card(fields, card_name, comment='', ifile=None,
                               is_list=True, has_none=True)
            continue
        elif keyi in str_attrs:
            value = cast_string(sub_group, encoding)
            #print(f'adding key={keyi!r} value={value!r}')
            assert isinstance(value, str), value

            try:
                setattr(model, keyi, value)
            except RuntimeError:  # pragma: no cover
                model.log.error('cant set minor_attributes/%s as %s' % (keyi, value))
            except AttributeError:  # pragma: no cover
                model.log.warning('cant set minor_attributes/%s as %s' % (keyi, value))
                raise
            continue
        elif keyi == 'is_enddata':
            model.card_count['ENDDATA'] = 1
            continue

        value = _cast(sub_group)
        try:
            setattr(model, keyi, value)
        except AttributeError:  # pragma: no cover
            model.log.warning('cant set minor_attributes/%s as %s' % (keyi, value))
            raise

    return


def _load_indexed_list(key, group, unused_encoding):
    lst = []
    for key in group.keys():
        value = _cast(group[key])
        lst.append(value)
    #print('_load_indexed_list: %s' % lst)
    return lst

def _load_indexed_list_str(key, group, encoding):
    lst = _load_indexed_list(key, group, encoding)

    #try:
        #value0 = value[0]
    #except IndexError:  # pragma: no cover
        #print('key =', key)
        #print('value = %r' % value)
        #print('group =', group)
        #print('group.keys() =', list(group.keys()))
        #raise

    if isinstance(lst, str):
        pass
    else:
        lst = decode_lines(lst, encoding)
        assert isinstance(lst[0], str), type(lst[0])
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

            cids = _cast_array(coords['cid'])
            rids = _cast_array(coords['rid'])
            e1s = _cast_array(coords['e1'])
            e2s = _cast_array(coords['e2'])
            e3s = _cast_array(coords['e3'])
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

            cids = _cast_array(coords['cid'])
            nodes = _cast_array(coords['nodes'])
            for cid, (n1, n2, n3) in zip(cids, nodes):
                func(cid, n1, n2, n3, comment='')
        else:
            cids, values = load_cards_from_keys_values(
                'coords/%s' % card_type,
                coords, encoding, model.log)
            _put_keys_values_into_dict(model, 'coords', cids, values)
        model.card_count[card_type] = len(cids)

def hdf5_load_tables(model: BDF, group, encoding: str) -> None:
    """loads the tables"""
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        keys, values = load_cards_from_keys_values(
            'tables/%s' % card_type,
            sub_group, encoding, model.log)
        _put_keys_values_into_dict(model, 'tables', keys, values)
        model.card_count[card_type] = len(keys)

def hdf5_load_methods(model: BDF, group, encoding: str) -> None:
    """loads the methods"""
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'EIGRL':
            #pass
        keys, values = load_cards_from_keys_values(
            'methods/%s' % card_type,
            sub_group, encoding, model.log)
        _put_keys_values_into_dict(model, 'methods', keys, values)
        model.card_count[card_type] = len(keys)

def hdf5_load_masses(model: BDF, group, encoding: str) -> None:
    """loads the masses"""
    for card_type in group.keys():
        masses = group[card_type]
        if card_type == 'CONM2':
            eid = _cast_array(masses['eid'])
            nid = _cast_array(masses['nid'])
            cid = _cast_array(masses['cid'])
            X = _cast_array(masses['X'])
            I = _cast_array(masses['I'])
            mass = _cast_array(masses['mass'])
            for eidi, nidi, cidi, Xi, Ii, massi in zip(eid, nid, cid, X, I, mass):
                model.add_conm2(eidi, nidi, massi, cid=cidi, X=Xi, I=Ii, comment='')
        elif card_type == 'CMASS2':
            eid = _cast_array(masses['eid'])
            mass = _cast_array(masses['mass'])
            nodes = _cast_array(masses['nodes']).tolist()
            components = _cast(masses['components'])
            for eidi, massi, nids, (c1, c2) in zip(eid, mass, nodes, components):
                model.add_cmass2(eidi, massi, nids, c1, c2, comment='')

        else:
            #model.add_cmass1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cmass3(eid, pid, nids, comment='')
            #model.add_cmass4(eid, mass, nids, comment='')
            #model.add_conm1(eid, nid, mass_matrix, cid=0, comment='')
            eid, values = load_cards_from_keys_values(
                'masses/%s' % card_type,
                masses, encoding, model.log)
            _put_keys_values_into_dict(model, 'masses', eid, values)
        model.card_count[card_type] = len(eid)


def hdf5_load_materials(model: BDF, group, encoding: str) -> None:
    """loads the materials"""
    for card_type in group.keys():
        sub_group = group[card_type]
        if card_type == 'MAT1':
            mid = _cast_array(sub_group['mid'])
            E = _cast_array(sub_group['E'])
            G = _cast_array(sub_group['G'])
            nu = _cast_array(sub_group['nu'])
            rho = _cast_array(sub_group['rho'])
            a = _cast_array(sub_group['A'])
            tref = _cast_array(sub_group['tref'])
            ge = _cast_array(sub_group['ge'])
            St = _cast_array(sub_group['St'])
            Sc = _cast_array(sub_group['Sc'])
            Ss = _cast_array(sub_group['Ss'])
            mcsid = _cast_array(sub_group['mcsid'])
            for midi, Ei, Gi, nui, rhoi, ai, trefi, gei, Sti, Sci, Ssi, mcsidi in zip(
                    mid, E, G, nu, rho, a, tref, ge, St, Sc, Ss, mcsid):
                model.add_mat1(midi, Ei, Gi, nui, rho=rhoi, a=ai, tref=trefi,
                               ge=gei, St=Sti, Sc=Sci, Ss=Ssi, mcsid=mcsidi, comment='')

        elif card_type == 'MAT2':
            mid = _cast(sub_group['mid'])
            G = _cast_array(sub_group['G'])
            rho = _cast_array(sub_group['rho'])
            a = _cast_array(sub_group['A'])
            tref = _cast_array(sub_group['tref'])
            ge = _cast_array(sub_group['ge'])
            St = _cast_array(sub_group['St'])
            Sc = _cast_array(sub_group['Sc'])
            Ss = _cast_array(sub_group['Ss'])
            mcsid = _cast_array(sub_group['mcsid'])

            for (midi, (G11, G22, G33, G12, G13, G23), rhoi, (a1i, a2i, a3i),
                 trefi, gei, Sti, Sci, Ssi, mcsidi) in zip(
                     mid, G, rho, a, tref, ge, St, Sc, Ss, mcsid):
                if mcsidi == -1:
                    mcsidi = None
                model.add_mat2(midi, G11, G12, G13, G22, G23, G33, rho=rhoi,
                               a1=a1i, a2=a2i, a3=a3i, tref=trefi, ge=gei,
                               St=Sti, Sc=Sci, Ss=Ssi, mcsid=mcsidi, comment='')

        elif card_type == 'MAT3':
            mid = _cast_array(sub_group['mid'])
            ex = _cast_array(sub_group['Ex'])
            eth = _cast_array(sub_group['Eth'])
            ez = _cast_array(sub_group['Ez'])

            nuxth = _cast_array(sub_group['Nuxth'])
            nuzx = _cast_array(sub_group['Nuzx'])
            nuthz = _cast_array(sub_group['Nuthz'])
            gxz = _cast_array(sub_group['Gzx'])

            ax = _cast_array(sub_group['Ax'])
            ath = _cast_array(sub_group['Ath'])
            az = _cast_array(sub_group['Az'])

            rho = _cast_array(sub_group['rho'])
            tref = _cast_array(sub_group['tref'])
            ge = _cast_array(sub_group['ge'])
            for (midi, exi, ethi, ezi, nuxthi, nuzxi, nuthzi,
                 rhoi, gzxi, axi, athi, azi, trefi, gei) in zip(
                     mid, ex, eth, ez, nuxth, nuzx, nuthz, rho, gxz, ax, ath, az, tref, ge):
                model.add_mat3(midi, exi, ethi, ezi, nuxthi, nuthzi, nuzxi, rho=rhoi,
                               gzx=gzxi, ax=axi, ath=athi, az=azi, tref=trefi, ge=gei, comment='')

        elif card_type == 'MAT8':
            mid = _cast_array(sub_group['mid'])
            e11 = _cast_array(sub_group['E11'])
            e22 = _cast_array(sub_group['E22'])
            nu12 = _cast_array(sub_group['Nu12'])
            g12 = _cast_array(sub_group['G12'])
            g1z = _cast_array(sub_group['G1z'])
            g2z = _cast_array(sub_group['G2z'])

            a1 = _cast_array(sub_group['A1'])
            a2 = _cast_array(sub_group['A2'])
            tref = _cast_array(sub_group['tref'])
            ge = _cast_array(sub_group['ge'])
            rho = _cast_array(sub_group['rho'])

            xt = _cast_array(sub_group['Xt'])
            xc = _cast_array(sub_group['Xc'])
            yt = _cast_array(sub_group['Yt'])
            yc = _cast_array(sub_group['Yc'])
            s = _cast_array(sub_group['S'])

            f12 = _cast_array(sub_group['F12'])
            strn = _cast_array(sub_group['strn'])
            for (midi, e11i, e22i, nu12i, g12i, g1zi, g2zi, rhoi, a1i, a2i, trefi,
                 xti, xci, yti, yci, si, gei, f12i, strni) in zip(
                     mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, tref,
                     xt, xc, yt, yc, s, ge, f12, strn):
                model.add_mat8(midi, e11i, e22i, nu12i, g12=g12i, g1z=g1zi, g2z=g2zi, rho=rhoi,
                               a1=a1i, a2=a2i, tref=trefi, Xt=xti, Xc=xci, Yt=yti, Yc=yci,
                               S=si, ge=gei, F12=f12i, strn=strni, comment='')
        elif card_type == 'MAT9':
            ## TODO: add G
            mid = _cast_array(sub_group['mid'])
            a = _cast_array(sub_group['A'])
            tref = _cast_array(sub_group['tref'])
            ge = _cast_array(sub_group['ge'])
            rho = _cast_array(sub_group['rho'])
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
            mid, values = load_cards_from_keys_values(
                'materials/%s' % card_type,
                sub_group, encoding, model.log)
            _put_keys_values_into_dict(model, 'materials', mid, values)
        model.card_count[card_type] = len(mid)

def hdf5_load_spcs(model: BDF, group, encoding: str) -> None:
    """loads the spcs"""
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for spc_id in keys:
        ispc_id = int(spc_id)
        cards_group = group[spc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'SPC1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'spcs/%s/%s' % (spc_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'spcs', ispc_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_spcadds(model: BDF, group, encoding: str) -> None:
    """loads the spcadds"""
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for spc_id in keys:
        ispc_id = int(spc_id)
        cards_group = group[spc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'SPC1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'spcadds/%s/%s' % (spc_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'spcadds', ispc_id, lkeys, values)

def hdf5_load_mpcs(model: BDF, group, encoding: str) -> None:
    """loads the mpcs"""
    keys = list(group.keys())
    keys.remove('keys')
    #mpc_ids = _cast(group['keys'])
    for mpc_id in keys:
        impc_id = int(mpc_id)
        cards_group = group[mpc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MPC':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'mpcs/%s/%s' % (mpc_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'mpcs', impc_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_mpcadds(model: BDF, group, encoding: str) -> None:
    """loads the mpcadds"""
    keys = list(group.keys())
    keys.remove('keys')
    #spc_ids = _cast(group['keys'])
    for mpc_id in keys:
        unused_impc_id = int(mpc_id)
        cards_group = group[mpc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MPCADD':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'mpcadds/%s/%s' % (mpc_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'mpcadds', mpc_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_pval(model: BDF, group, encoding: str) -> None:
    """loads the pval"""
    keys = list(group.keys())
    keys.remove('keys')
    for adapt_id in keys:
        adapt_idi = int(adapt_id)
        cards_group = group[adapt_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'TEMP':  # this has a weird dictionary structure
                #sid = sub_group.keys()
                #for index in sid:
                    #cardi = sub_group[index]
                    #nodes = _cast(cardi['node']).tolist()
                    #temp = _cast(cardi['temperature']).tolist()
                    #temperatures = {nid : tempi for (nid, tempi) in zip(nodes, temp)}
                    #model.add_temp(iload_id, temperatures, comment='')
            #else:
            sid, values = load_cards_from_keys_values(
                'pval/%s/%s' % (adapt_idi, card_type),
                sub_group, encoding, model.log)
            #for value in values:
                #print(value)
            _put_keys_values_into_dict_list(model, 'pval', adapt_idi, sid, values)
            model.card_count[card_type] = len(sid)

def hdf5_load_loads(model: BDF, group, encoding: str) -> None:
    """loads the loads"""
    keys = list(group.keys())
    keys.remove('keys')
    for load_id in keys:
        iload_id = int(load_id)
        cards_group = group[load_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            if card_type in ['FORCE', 'MOMENT']:
                if card_type == 'FORCE':
                    func = model.add_force
                else:
                    func = model.add_moment
                sid = _cast_array(sub_group['sid'])
                node = _cast_array(sub_group['node'])
                cid = _cast_array(sub_group['cid'])
                mag = _cast_array(sub_group['mag'])
                xyz = _cast_array(sub_group['xyz'])
                for (sidi, nodei, magi, xyzi, cidi) in zip(sid, node, mag, xyz, cid):
                    func(sidi, nodei, magi, xyzi, cid=cidi, comment='')
            elif card_type == 'TEMP':  # this has a weird dictionary structure
                sid = sub_group.keys()
                for index in sid:
                    cardi = sub_group[index]
                    nodes = _cast_array(cardi['node']).tolist()
                    temp = _cast_array(cardi['temperature']).tolist()
                    temperatures = {nid : tempi for (nid, tempi) in zip(nodes, temp)}
                    model.add_temp(iload_id, temperatures, comment='')
            else:
                #model.add_force1(sid, node, mag, g1, g2, comment='')
                sid, values = load_cards_from_keys_values(
                    'loads/%s/%s' % (load_id, card_type),
                    sub_group, encoding, model.log)
                #for value in values:
                    #print(value)
                _put_keys_values_into_dict_list(model, 'loads', iload_id, sid, values)
            model.card_count[card_type] = len(sid)

def hdf5_load_load_combinations(model: BDF, group, encoding: str) -> None:
    """loads the load_combinations"""
    keys = list(group.keys())
    keys.remove('keys')
    for load_id in keys:
        iload_id = int(load_id)
        cards_group = group[load_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'LOAD':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'load_combinations/%s/%s' % (load_id, card_type),
                sub_group, encoding, model.log)
            #for value in values:
                #print(value)
            _put_keys_values_into_dict_list(model, 'load_combinations', iload_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_nsms(model: BDF, group, encoding: str) -> None:
    """loads the nsms"""
    keys = list(group.keys())
    keys.remove('keys')
    for nsm_id in keys:
        insm_id = int(nsm_id)
        cards_group = group[nsm_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'NSM':
                #mid = _cast(sub_group['mid'])
            #else:
            keys, values = load_cards_from_keys_values(
                'nsms/%s/%s' % (nsm_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'nsms', insm_id, keys, values)
            model.card_count[card_type] = len(keys)

def hdf5_load_nsmadds(model: BDF, group, encoding: str) -> None:
    """loads the nsmadds"""
    keys = list(group.keys())
    keys.remove('keys')
    for nsm_id in keys:
        insm_id = int(nsm_id)
        cards_group = group[nsm_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'NSMADD':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'nsmadds/%s/%s' % (nsm_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'nsmadds', insm_id, lkeys, values)
            model.card_count[card_type] = len(keys)

def hdf5_load_frequencies(model: BDF, group, encoding: str) -> None:
    """loads the frequencies"""
    keys = list(group.keys())
    keys.remove('keys')
    for freq_id in keys:
        ifreq_id = int(freq_id)
        cards_group = group[freq_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'FREQ':
                #mid = _cast(sub_group['mid'])
            #else:
            fkeys, values = load_cards_from_keys_values(
                'frequencies/%s/%s' % (freq_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'frequencies', ifreq_id, fkeys, values)
            model.card_count[card_type] = len(fkeys)

def hdf5_load_aelinks(model: BDF, group, encoding: str) -> None:
    """loads the aelinks"""
    keys = group.keys()
    naelinks = 0
    add_methods = model._add_methods
    for aelink_id in keys:
        unused_iaelink_id = int(aelink_id)
        jlinks_group = group[aelink_id]
        keys = jlinks_group.keys()
        aelink = [None] * len(keys)
        for jlink in keys:
            j_int = int(jlink)
            aelinki_group = jlinks_group[jlink]
            value = aelinki_group
            aelinki = _load_class(jlink, value, 'AELINK', encoding)
            aelink[j_int] = aelinki
            naelinks += 1
        for aelinki in aelink:
            add_methods._add_aelink_object(aelinki)
    model.card_count['AELINK'] = naelinks

def hdf5_load_dloads(model: BDF, group, encoding: str) -> None:
    """loads the dloads"""
    keys = list(group.keys())
    keys.remove('keys')
    for dload_id in keys:
        idload_id = int(dload_id)
        cards_group = group[dload_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'DLOAD':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'dloads/%s/%s' % (dload_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'dloads', idload_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_dload_entries(model: BDF, group, encoding: str) -> None:
    """loads the dload_entries"""
    keys = list(group.keys())
    keys.remove('keys')
    for dload_id in keys:
        idload_id = int(dload_id)
        cards_group = group[dload_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'TLOAD1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'dload_entries/%s/%s' % (dload_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'dload_entries', idload_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_bcs(model: BDF, group, encoding: str) -> None:
    """loads the bcs"""
    keys = list(group.keys())
    keys.remove('keys')
    for bc_id in keys:
        ibc_id = int(bc_id)
        cards_group = group[bc_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'bcs/%s/%s' % (bc_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'bcs', ibc_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_transfer_functions(model: BDF, group, encoding: str) -> None:
    """loads the transfer_functions"""
    keys = list(group.keys())
    keys.remove('keys')
    for tf_id in keys:
        itf_id = int(tf_id)
        cards_group = group[tf_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'transfer_functions/%s/%s' % (tf_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'transfer_functions', itf_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_dvgrids(model: BDF, group, encoding: str) -> None:
    """loads the dvgrids"""
    keys = list(group.keys())
    keys.remove('keys')
    for opt_id in keys:
        iopt_id = int(opt_id)
        cards_group = group[opt_id]
        for card_type in cards_group.keys():
            sub_group = cards_group[card_type]
            #if card_type == 'MAT1':
                #mid = _cast(sub_group['mid'])
            #else:
            lkeys, values = load_cards_from_keys_values(
                'dvgrids/%s/%s' % (opt_id, card_type),
                sub_group, encoding, model.log)
            _put_keys_values_into_dict_list(model, 'dvgrids', iopt_id, lkeys, values)
            model.card_count[card_type] = len(lkeys)

def hdf5_load_desvars(model: BDF, group, encoding: str) -> None:
    """loads the desvars"""
    for card_type in group.keys():
        sub_group = group[card_type]
        if card_type == 'DESVAR':
            desvar = _cast_array(sub_group['desvar'])
            label = _cast(sub_group['label'])
            xinit = _cast_array(sub_group['xinit'])
            xlb = _cast_array(sub_group['xlb'])
            xub = _cast_array(sub_group['xub'])
            delx = _cast_array(sub_group['delx'])
            ddval = _cast_array(sub_group['ddval'])
            for desvari, labeli, xiniti, xlbi, xubi, delxi, ddvali in zip(
                    desvar, label, xinit, xlb, xub, delx, ddval):
                labeli = labeli.decode(encoding)
                assert isinstance(labeli, str), labeli
                model.add_desvar(desvari, labeli, xiniti, xlb=xlbi, xub=xubi,
                                 delx=delxi, ddval=ddvali, comment='')
        else:  # pragma: no cover
            raise RuntimeError('card_type=%s in hdf5_load_desvars' % card_type)
        model.card_count[card_type] = len(desvar)

def hdf5_load_dmigs(model: BDF, group, unused_encoding: str) -> None:
    """loads the dmigs"""
    keys = group.keys()
    if len(keys) == 0:
        #model.log.warning('skipping loading %s' % group)
        raise RuntimeError('error loading %s' % group)
        #return

    for name in keys:
        sub_group = group[name]
        #print('group', group)
        #print('sub_group', sub_group)

        class_type = group.attrs['type']
        if class_type == 'DMIG' and name == 'UACCEL':
            _load_dmig_uaccel(model, sub_group)
        elif class_type == 'DMI':
            _load_dmi(model, name, sub_group)
        elif class_type == 'DMIAX':
            _load_dmiax(model, name, sub_group)
        else:
            _load_dmig(model, name, sub_group, class_type)
    model.card_count[class_type] = len(keys)


def _load_dmig_uaccel(model: BDF, sub_group):
    """loads the DMIG,UACCEL"""
    keysi = list(sub_group.keys())
    tin = _cast(sub_group['tin'])
    keysi.remove('tin')
    ncol = None
    if 'ncol' in keysi:
        keysi.remove('ncol')
        ncol = _cast(sub_group['ncol'])

    load_sequences = {}
    for idi in keysi:
        lseq = int(idi)
        sub_groupi = sub_group[idi]
        dofs = _cast(sub_groupi['dofs'])
        nids = _cast(sub_groupi['nids'])
        values = _cast(sub_groupi['values'])
        load_sequences[lseq] = list([
            [nid, dof, value] for (nid, dof, value)
            in zip(nids, dofs, values)])
    dmig_uaccel = model.add_dmig_uaccel(tin, ncol, load_sequences, comment='')
    str(dmig_uaccel)

def _load_dmi(model: BDF, name, sub_group):
    """loads the DMI"""
    ncols = _cast(sub_group['ncols'])
    nrows = _cast(sub_group['nrows'])
    #polar = _cast(sub_group['polar'])
    matrix_form = _cast(sub_group['matrix_form'])
    tin = _cast(sub_group['tin'])
    tout = _cast(sub_group['tout'])
    GCi = _cast(sub_group['GCi'])
    GCj = _cast(sub_group['GCj'])
    Real = _cast(sub_group['Real'])
    Complex = None
    if 'Complex' in sub_group:
        Complex = _cast(sub_group['Complex'])

    #ifo = matrix_form
    form = matrix_form
    model.add_dmi(name, form, tin, tout, nrows, ncols, GCj, GCi,
                  Real, Complex=Complex, comment='')

def _load_dmig(model, name, sub_group, class_type):
    """loads the DMIG, DMIJ, DMIJI, DMIK"""
    class_obj = CARD_MAP[class_type]
    ncols = None
    if 'ncols' in sub_group:
        ncols = _cast(sub_group['ncols'])
    polar = _cast(sub_group['polar'])
    matrix_form = _cast(sub_group['matrix_form'])
    tin = _cast(sub_group['tin'])
    tout = _cast(sub_group['tout'])
    #dmig_group.create_dataset('tin_dtype', data=dmig.tin_dtype)
    #dmig_group.create_dataset('tout_dtype', data=dmig.tout_dtype)

    #dmig_group.create_dataset('matrix_type', data=dmig.matrix_type)
    #dmig_group.create_dataset('is_complex', data=dmig.is_complex)
    #dmig_group.create_dataset('is_real', data=dmig.is_real)
    #dmig_group.create_dataset('is_polar', data=dmig.is_polar)

    GCi = _cast(sub_group['GCi'])
    GCj = _cast(sub_group['GCj'])
    Real = _cast(sub_group['Real'])
    Complex = None
    if 'Complex' in sub_group:
        Complex = _cast(sub_group['Complex'])

    ifo = matrix_form
    dmig = class_obj(name, ifo, tin, tout, polar, ncols,
                     GCj, GCi, Real, Complex=Complex, comment='', finalize=True)
    assert class_type in ['DMIG', 'DMIK', 'DMIJ', 'DMIJI'], class_type
    slot_name = class_type.lower() + 's'
    slot = getattr(model, slot_name)
    slot[name] = dmig
    str(dmig)
    #model.dmigs[name] = dmig

def _load_dmiax(model, name, sub_group):
    """loads the DMIAX"""
    class_obj = CARD_MAP['DMIAX']
    ncols = None
    if 'ncols' in sub_group:
        ncols = _cast(sub_group['ncols'])
    matrix_form = _cast(sub_group['matrix_form'])
    tin = _cast(sub_group['tin'])
    tout = _cast(sub_group['tout'])

    gcni = _cast(sub_group['GCNi_j'])
    gcnj = _cast(sub_group['GCNj'])
    i_none_flags = _cast(sub_group['i_none_flags'])
    j_none_flags = _cast(sub_group['j_none_flags'])
    dmiax_GCNi = []
    dmiax_GCNj = []

    k = 0
    for GCNj, is_none_flag_j in zip(gcnj, j_none_flags):
        gj, cj, nj = GCNj
        if is_none_flag_j:
            nj = None
        dmiax_GCNj.append((gj, cj, nj))
    del GCNj, is_none_flag_j, nj

    j_old = -1
    gcni_group = []
    for GCNi_j, is_none_flag_j in zip(gcni, i_none_flags):
        gi, ci, ni, j = GCNi_j
        is_none_flag_i = i_none_flags[k]
        if is_none_flag_i:
            ni = None
        if j != j_old:
            j_old = j
            gcni_group = []
            dmiax_GCNi.append(gcni_group)
        gcni_group.append((gi, ci, ni))
    #print('GCNj =', dmiax_GCNj)
    #print('GCNi =', dmiax_GCNi)

    Real = _cast(sub_group['Real'])
    Complex = None
    if 'Complex' in sub_group:
        Complex = _cast(sub_group['Complex'])

    ifo = matrix_form
    dmiax = DMIAX(name, matrix_form, tin, tout, ncols,
                  dmiax_GCNj, dmiax_GCNi, Real, Complex=Complex)
    model.dmiax[name] = dmiax
    str(dmiax)
    #print(dmiax)

def hdf5_load_dconstrs(model, group, encoding):
    """loads the dconstrs"""
    keys = group.keys()
    if len(keys) == 0:
        #model.log.warning('skipping loading %s' % group)
        raise RuntimeError('error loading %s' % group)
        #return
    add_methods = model._add_methods
    for card_type in keys:
        sub_group = group[card_type]
        #print('group', group)
        #print('sub_group', sub_group)

        if card_type == 'DCONSTR':
            #keys_group = list(sub_group.keys())
            oid = _cast(sub_group['oid'])
            dresp_id = _cast(sub_group['dresp_id'])
            lid = _cast(sub_group['lid'])
            uid = _cast(sub_group['uid'])
            lowfq = _cast(sub_group['lowfq'])
            highfq = _cast(sub_group['highfq'])

            for oidi, dresp_idi, lidi, uidi, lowfqi, highfqi in zip(
                    oid, dresp_id, lid, uid, lowfq, highfq):
                model.add_dconstr(oidi, dresp_idi, lid=lidi, uid=uidi,
                                  lowfq=lowfqi, highfq=highfqi, comment='')

        elif card_type == 'DCONADD':
            keys = sub_group.keys()
            #print('keys_group', keys_group)
            #debug = False
            unused_name = 'dconstrs/%s' % card_type

            for key in keys:
                value = sub_group[key]
                dconadd = _load_class(key, value, card_type, encoding)
                add_methods._add_dconstr_object(dconadd)
                #model.add_dconadd(oid, dconstrs, comment='')
        else:
            raise RuntimeError('error loading %s' % card_type)

        model.card_count[card_type] = len(keys)

def hdf5_load_dti(model, group, encoding):
    """loads the dti"""
    group_keys = group.keys()
    if len(group_keys) == 0:
        #model.log.warning('skipping loading %s' % group)
        raise RuntimeError('error loading %s' % group)

    names = cast_strings(group['keys'], encoding)
    values = group['values']
    for name in names:
        sub_group = values[name]
        records = sub_group.keys()

        fields = {}
        #print('records', records)
        for irecord in records:
            sub_groupi = sub_group[irecord]
            #print(sub_group, sub_groupi)
            if 'keys' in sub_groupi:
                lst = _load_indexed_list(irecord, sub_groupi, encoding)
                lst2 = [val.decode(encoding) if isinstance(val, bytes) else val for val in lst]
            else:
                if isinstance(sub_groupi, h5py._hl.dataset.Dataset):
                    #print(sub_group, sub_groupi)
                    lst = _cast(sub_groupi).tolist()
                    lst2 = [val.decode(encoding) if isinstance(val, bytes) else val for val in lst]
                else:
                    #print(sub_group, sub_groupi, len(sub_groupi.keys()))
                    keys = sub_groupi.keys()
                    lst = []
                    for key in keys:
                        sub_groupii = sub_groupi[key]
                        if len(sub_groupii.shape) == 0:

                            # h5py between 2.8.0 and 3.1.0 (probably 3.0)
                            # changed str to bytes
                            scalar_value = np.array(sub_groupii).tolist()

                            # str/bytes/float
                            if isinstance(scalar_value, str):
                                pass
                            elif isinstance(scalar_value, bytes):
                                scalar_value = scalar_value.decode(encoding)
                            elif np.isnan(scalar_value):
                                scalar_value = None
                            lst.append(scalar_value)
                        else:
                            lsti = _cast(sub_groupii)
                            #assert isinstance(lsti, int, float, str), lsti
                            lst.append(lsti)
                    #lst = _cast(sub_groupi)
                    #print(lst)
                    lst2 = lst
            if name == 'UNITS':
                fields[irecord] = lst2[0]
            else:
                fields[irecord] = lst2
        assert len(fields) > 0, fields
        model.add_dti(name, fields)
    model.card_count['DTI'] = len(names)

def hdf5_load_usets(model, group, encoding):
    """loads the usets"""
    keys = group.keys()
    if len(keys) == 0:
        #model.log.warning('skipping loading %s' % group)
        raise RuntimeError('error loading %s' % group)
    add_methods = model._add_methods
    for name in keys:
        sub_group = group[name]
        keys = sub_group.keys()
        unused_lst = [None] * len(keys)

        for key in keys:
            sub_groupi = sub_group[key]
            unused_keys2 = sub_groupi.keys()

            value = sub_groupi
            card_type = _cast(sub_groupi['type'])
            class_obj = _load_class(key, value, card_type, encoding)
            add_methods._add_uset_object(class_obj)
            if card_type not in model.card_count:
                model.card_count[card_type] = 1
            else:
                model.card_count[card_type] += 1

def hdf5_load_dresps(model, group, encoding):
    """loads the dresps"""
    keys = list(group.keys())
    if len(keys) == 0:
        #model.log.warning('skipping loading %s' % group)
        raise RuntimeError('error loading %s' % group)

    for class_type in group.keys():
        sub_group = group[class_type]

        if class_type == 'DRESP1':
            unused_keys_group = list(sub_group.keys())
            #print('keys_group', keys_group)

            #'atta', u'attb', u'dresp_id', u'label', u'region', u'response_type'
            dresp_id = _cast_array(sub_group['dresp_id'])
            atta = _cast(sub_group['atta'])
            #print('atta =', atta)
            attb = _cast(sub_group['attb'])
            label = _cast_array(sub_group['label'])
            region = _cast_array(sub_group['region'])
            response_type = _cast_array(sub_group['response_type'])
            property_type = _cast_array(sub_group['property_type'])
            atti = []
            for (i, dresp_idi, labeli, response_typei, property_typei, regioni,
                 attai, attbi) in zip(count(), dresp_id, label, response_type, property_type,
                                      region, atta, attb):
                drespi_group = sub_group[str(i)]

                labeli = labeli.decode(encoding)
                response_typei = response_typei.decode(encoding)

                if property_typei == b'':
                    property_typei = None
                elif property_typei.isdigit():
                    property_typei = int(property_typei)
                else:
                    property_typei = property_typei.decode(encoding)

                if regioni == -1:
                    regioni = None
                #else:
                    #regioni = regioni.decode(encoding)

                # int, float, str, blank
                if attai == b'':
                    attai = None
                elif b'.' in attai:
                    attai = float(attai)
                elif attai.isdigit():
                    attai = int(attai)
                else:
                    attai = attai.decode(encoding)

                # int, float, str, blank
                if attbi == b'':
                    attbi = None
                elif b'.' in attbi:
                    attbi = float(attbi)
                elif attbi.isdigit():
                    attbi = int(attbi)
                else:
                    attbi = attbi.decode(encoding)

                atti = []
                if 'atti' in drespi_group:
                    atti = _cast_array(drespi_group['atti']).tolist()

                model.add_dresp1(dresp_idi, labeli, response_typei, property_typei, regioni,
                                 attai, attbi, atti, validate=False, comment='')

        elif class_type == 'DRESP2':
            dresp_id = _cast_array(sub_group['dresp_id'])
            label = _cast(sub_group['label'])
            dequation = _cast(sub_group['dequation'])
            dequation_str = _cast(sub_group['func'])
            #dequation_str = _cast(sub_group['dequation_str'])
            region = _cast_array(sub_group['region'])
            method = _cast(sub_group['method'])
            c123 = _cast(sub_group['c123'])

            for (i, dresp_idi, labeli, dequationi, dequation_stri, regioni, methodi, c123i) in zip(
                    count(), dresp_id, label, dequation, dequation_str, region, method, c123):
                c1, c2, c3 = c123i
                if regioni == -1:
                    regioni = None
                #paramsi = {(0, u'DESVAR'): [1, 2, 3]}
                paramsi = {}
                dresp_groupi = sub_group[str(i)]
                param_keys = _cast(dresp_groupi['param_keys'])
                #print('param_keys', param_keys)

                for j, param_key_j in enumerate(param_keys):
                    param_values = _cast(dresp_groupi[str(j)]['values'])
                    param_key = param_key_j.decode(encoding)
                    #print('  param_values', (i, j), param_values)
                    param_values2 = [val.decode(encoding) if isinstance(val, bytes) else val
                                     for val in param_values]
                    paramsi[(j, param_key)] = param_values2
                model.log.debug('DRESP2 params = %s' % paramsi)

                if dequationi == -1:
                    dequationi = dequation_stri.decode(encoding)

                labeli = labeli.decode(encoding)
                methodi = methodi.decode(encoding)
                model.add_dresp2(dresp_idi, labeli, dequationi, regioni, paramsi,
                                 method=methodi, c1=c1, c2=c2, c3=c3,
                                 validate=False, comment='')
        else:
            raise RuntimeError('error loading %s' % class_type)

        model.card_count[class_type] = len(dresp_id)

def hdf5_load_generic(model, group, name, encoding):
    for card_type in group.keys():
        sub_group = group[card_type]
        #if card_type == 'TABLES1':
            #pass
        lkeys, values = load_cards_from_keys_values(
            '%s/%s' % (name, card_type),
            sub_group, encoding, model.log)
        _put_keys_values_into_dict(model, name, lkeys, values)
        model.card_count[card_type] = len(lkeys)


def hdf5_load_properties(model, properties_group, encoding):
    """loads the properties from an HDF5 file"""
    for card_type in properties_group.keys():
        properties = properties_group[card_type]
        if card_type == 'PSHELL':
            pid = _cast_array(properties['pid'])
            mids = _cast_array(properties['mids'])
            z = _cast_array(properties['z'])
            t = _cast_array(properties['t'])
            twelveIt3 = _cast_array(properties['twelveIt3'])
            tst = _cast_array(properties['tst'])
            nsm = _cast_array(properties['nsm'])
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
            pid = _cast_array(properties['pid'])
            mid = _cast_array(properties['mid'])
            cordm = _cast_array(properties['cordm'])
            integ = _cast_array(properties['integ'])
            isop = _cast_array(properties['isop'])
            stress = _cast_array(properties['stress'])
            fctn = _cast_array(properties['fctn'])
            for pidi, midi, cordmi, integi, stressi, isopi, fctni in zip(
                    pid, mid, cordm, integ, stress, isop, fctn):
                integi = integi.decode(encoding)
                fctni = fctni.decode(encoding)
                isopi = isopi.decode(encoding)
                stressi = stressi.decode(encoding)
                if integi == '':
                    integi = None
                if fctni == '':
                    fctni = None
                if isopi == '':
                    isopi = None
                if stressi == '':
                    stressi = None
                func(pidi, midi, cordm=cordmi, integ=integi, stress=stressi,
                     isop=isopi, fctn=fctni, comment='')

        elif card_type == 'PROD':
            pid = _cast_array(properties['pid'])
            mid = _cast_array(properties['mid'])
            A = _cast_array(properties['A'])
            j = _cast_array(properties['J'])
            c = _cast_array(properties['c'])
            nsm = _cast_array(properties['nsm'])
            for pidi, midi, Ai, ji, ci, nsmi in zip(
                    pid, mid, A, j, c, nsm):
                model.add_prod(pidi, midi, Ai, j=ji, c=ci, nsm=nsmi, comment='')

        elif card_type == 'PTUBE':
            pid = _cast_array(properties['pid'])
            mid = _cast_array(properties['mid'])
            OD = _cast_array(properties['OD'])
            t = _cast_array(properties['t'])
            nsm = _cast_array(properties['nsm'])
            for pidi, midi, (OD1, OD2), ti, nsmi in zip(
                    pid, mid, OD, t, nsm):
                model.add_ptube(pidi, midi, OD1, t=ti, nsm=nsmi, OD2=OD2, comment='')

        elif card_type == 'PBAR':
            pid = _cast_array(properties['pid'])
            mid = _cast_array(properties['mid'])
            A = _cast_array(properties['A'])
            J = _cast_array(properties['J'])
            I = _cast_array(properties['I'])

            c = _cast_array(properties['c'])
            d = _cast_array(properties['d'])
            e = _cast_array(properties['e'])
            f = _cast_array(properties['f'])
            k = _cast_array(properties['k'])

            nsm = _cast_array(properties['nsm'])
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

            #if card_type == 'PCOMP':
                #debug = True
            pid, values = load_cards_from_keys_values(
                'properties/%s' % card_type,
                properties, encoding, model.log)
            _put_keys_values_into_dict(model, 'properties', pid, values)

            #model.add_pshear(pid, mid, t, nsm=0., f1=0., f2=0., comment='')
            #model.add_pvisc(pid, ce, cr, comment='')
            #model.add_pelas(pid, k, ge=0., s=0., comment='')
            #model.add_pdamp(pid, b, comment='')
            #model.add_pcomp(pid, mids, thicknesses, thetas=None, souts=None, nsm=0., sb=0.,
                            #ft=None, tref=0., ge=0., lam=None, z0=None, comment='')
            #model.add_pcompg(pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                             #nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None,
                             #comment='')
        model.card_count[card_type] = len(pid)

    for prop in model.properties.values():
        write_card(prop)

def _put_keys_values_into_dict(model, name: str, keys, values, cast_int_keys: bool=True) -> None:
    """add something like an element to a dictionary"""
    for value in values:
        write_card(value)

    slot = getattr(model, name)
    card_count = model.card_count

    # 'dmig', 'dmik', 'dmij', 'dmiji', 'dmi', 'dmiax'
    if cast_int_keys and name not in ['dscreen', 'dti', 'aecomps', 'seconct', 'sebndry']:
        #print('keys =', keys, cast_int_keys, name)
        try:
            keys = [int(key) for key in keys]
        except ValueError:  # pragma: no cover
            # If this hits, you need probably have a non-integer key
            # (e.g., a tuple of 2 ints) and need to skip the above
            # caster and figure out the right way to cast it.
            #
            # This could be a string (in which case you just pass the
            # initial check above and then use the normal adder below)
            # similar to 'dscreen'.
            #
            # Another possibility is you have a (int_a, int_b) tuple key.
            # Follow the pattern for 'seconct'.
            print('name =', name)
            print('keys = ', keys)
            print('values = ', values)
            raise

    tuple_integer_casts = ('seconct', 'sebndry')
    if name in tuple_integer_casts:
        for key, value in zip(keys, values):
            key = tuple(key)
            slot[key] = value
            #print('  *%s %s' % (value.type, key))
            card_type = value.type
            if card_type not in card_count:
                card_count[card_type] = 0
            card_count[card_type] += 1
            model._type_to_id_map[card_type].append(key)
    else:
        for key, value in zip(keys, values):
            slot[key] = value
            #print('  *%s %s' % (value.type, key))
            card_type = value.type
            if card_type not in card_count:
                card_count[card_type] = 0
            card_count[card_type] += 1
            model._type_to_id_map[card_type].append(key)

def _put_keys_values_into_list(model, name, keys, values):
    """add something like an MKAERO1 to a list"""
    for value in values:
        try:
            write_card(value)
        except RuntimeError:
            print(value)
            raise

    slot = getattr(model, name)
    card_count = model.card_count
    for key, value in zip(keys, values):
        slot.append(value)
        #print('  *%s %s' % (value.type, key))
        Type = value.type
        if Type not in card_count:
            card_count[Type] = 0
        card_count[Type] += 1
        model._type_to_id_map[Type].append(key)

def _put_keys_values_into_dict_list(model: Any, name: str, idi: int,
                                    keys: np.ndarray,
                                    values: List[Any]):
    """add something like an SPC into a dictionary that has a list"""
    for value in values:
        #print(value)
        write_card(value)

    slot = getattr(model, name)
    idi = int(idi)
    #assert isinstance(idi, int), 'idi=%s type=%s' % (idi, type(idi))
    if idi in slot:
        slot_list = slot[idi]
    else:
        slot_list = []
        slot[idi] = slot_list

    card_count = model.card_count
    for key, value in zip(keys, values):
        slot_list.append(value)
        #print('  *%s %s' % (value.type, key))
        Type = value.type
        if Type not in card_count:
            card_count[Type] = 0
        card_count[Type] += 1
        model._type_to_id_map[Type].append(key)

def load_cards_from_keys_values(name, properties, encoding, log):
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
    value_objs = _load_cards_from_keys_values(name, values, keys, encoding, log)
    return keys, value_objs

def _load_cards_from_keys_values(unused_name, values, keys, encoding, unused_log):
    value_objs = []
    for key, keyi in zip(keys, values.keys()):
        #print('%s - %s' % (name, key))
        value = values[keyi]
        card_type = _cast(value['type'])
        class_instance = _load_class(key, value, card_type, encoding)
        value_objs.append(class_instance)
    return value_objs

def _load_class(key: str, value, card_type: str, encoding: str):
    if isinstance(card_type, bytes):
        card_type = card_type.decode(encoding)
    keys_to_read = list(value.keys())
    class_obj = CARD_MAP[card_type]  # see add_card.py ~line 200
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

    #print('  keys_to_read = ', keys_to_read)
    for key_to_cast in keys_to_read:
        if key_to_cast in _properties:
            continue

        try:
            valuei = _get_casted_value(value, key_to_cast, encoding)
        except AssertionError:
            print('error loading %r' % card_type)
            print(_properties)
            print(key, key_to_cast)
            raise

        if isinstance(valuei, np.ndarray):
            valuei = valuei.tolist()
            if isinstance(valuei, list) and isinstance(valuei[0], bytes):
                valuei = [valueii.decode(encoding) for valueii in valuei]
        elif isinstance(valuei, list) and isinstance(valuei[0], bytes):
            valuei = [valueii.decode(encoding) for valueii in valuei]

        elif isinstance(valuei, bytes):
            raise TypeError(f'class={card_type} key={key_to_cast} value={valuei} must be a string (not bytes)')

        try:
            setattr(class_instance, key_to_cast, valuei)
            #print('  set %s to %s' % (key_to_cast, valuei))
        except AttributeError:  # pragma: no cover
            print('error loading %r' % card_type)
            print(_properties)
            print(key, key_to_cast, valuei)
            raise
    #if debug:
        #print(class_instance.get_stats())
        #print(class_instance)
    if hasattr(class_instance, '_finalize_hdf5'):
        class_instance._finalize_hdf5(encoding)
    #else:
        #print('no %s' % class_instance.type)
    str(class_instance)
    return class_instance

def _cast_encoding(value_h5, encoding: str):
    valuei = _cast(value_h5)
    if isinstance(valuei, (int, float, np.ndarray)):
        pass
    elif isinstance(valuei, bytes):
        valuei = valuei.decode(encoding)
    elif isinstance(valuei, list):
        valuei = [val.decode(encoding) if isinstance(val, bytes) else val
                  for val in valuei]
    #else:
        #print(type(valuei))
    return valuei

def _get_casted_value(value, key_to_cast: str, encoding: str) -> Any:
    value_h5 = value[key_to_cast]
    if isinstance(value_h5, h5py._hl.dataset.Dataset):
        valuei = _cast_encoding(value_h5, encoding)
    else:
        h5_keys = list(value_h5.keys())
        if len(h5_keys) == 0:
            valuei = _cast_encoding(value_h5, encoding)
        else:
            #print('h5_keys =', h5_keys)
            lst = []
            for h5_key in h5_keys:
                slot_h5 = value_h5[h5_key]

                if isinstance(slot_h5, h5py._hl.dataset.Dataset):
                    valueii = _cast(slot_h5)
                elif isinstance(slot_h5, h5py._hl.group.Group):
                    valueii = _load_indexed_list(h5_key, slot_h5, encoding)
                else:  # pragma: no cover
                    print(key_to_cast, h5_key)
                    print(slot_h5, type(slot_h5))
                    raise NotImplementedError()

                #print(f'key={key_to_cast}; valueii={valueii}; type={type(valueii)}')
                valueii = to_list_int_float_str(valueii, encoding)
                lst.append(valueii)
            valuei = lst

        #valuei = None
    #else:
    #try:
        #valuei = _cast(value_h5)
    #except AttributeError:
        #print(value, key_to_cast, value.keys())
        #print(value_h5, value_h5.keys())
        #raise
        #valuei = None
    #print(f'key={key_to_cast}; valuei={valuei}; type={type(valuei)}')
    #assert not isinstance(valuei, bytes), f'key={key_to_cast}; valuei={valuei}; type={type(valuei)}'
    return valuei

def to_list_int_float_str(valueii: Any, encoding: str) -> Any:
    if isinstance(valueii, (int, float, str)):
        pass
    elif isinstance(valueii, bytes):
        valueii = valueii.decode(encoding)
    elif isinstance(valueii, np.ndarray):
        valueii = valueii.tolist()
        if isinstance(valueii[0], bytes):
            valueii = [val.decode(encoding) if isinstance(val, bytes) else val
                       for val in valueii]
    elif isinstance(valueii, list):
        if len(valueii) > 0 and isinstance(valueii[0], bytes):
            valueii = [val.decode(encoding) if isinstance(val, bytes) else val
                       for val in valueii]
    else:
        print(valueii)
        raise NotImplementedError(type(valueii))
    return valueii

def _load_from_class(value, card_type: str, encoding: str):
    """generic loader that only requires an ``_init_from_empty`` method"""
    keys_to_read = list(value.keys())
    assert isinstance(card_type, str), card_type
    class_obj = CARD_MAP[card_type]  # see add_card.py ~line 200
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

        valuei = _get_casted_value(value, key_to_cast, encoding)
        #print('%s set to %r' % (key_to_cast, valuei))
        #h5_value = value[key_to_cast]
        #try:
            #valuei = _cast(h5_value)
        #except AttributeError:
            #print('key =', key)
            #raise
            #valuei = None

        try:
            setattr(class_instance, key_to_cast, valuei)
        except AttributeError:  # pragma: no cover
            print('error loading %r' % card_type)
            print(_properties)
            print(key_to_cast, valuei)
            raise

    if hasattr(class_instance, '_finalize_hdf5'):
        class_instance._finalize_hdf5(encoding)
    return class_instance


def hdf5_load_elements(model, elements_group, encoding):
    """loads the elements from an HDF5 file"""
    for card_type in elements_group.keys():
        elements = elements_group[card_type]
        if card_type == 'CTETRA':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctetra(eid, pid, nids, comment='')
        elif card_type == 'CPENTA':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_cpenta(eid, pid, nids, comment='')
        elif card_type == 'CPYRAM':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_cpyram(eid, pid, nids, comment='')
        elif card_type == 'CHEXA':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_chexa(eid, pid, nids, comment='')

        elif card_type == 'CROD':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_crod(eid, pid, nids, comment='')
        elif card_type == 'CTUBE':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_ctube(eid, pid, nids, comment='')
        elif card_type == 'CONROD':
            eids = _cast_array(elements['eid'])
            mids = _cast_array(elements['mid'])
            nodes = _cast_array(elements['nodes']).tolist()
            A = _cast_array(elements['A'])
            J = _cast_array(elements['J'])
            c = _cast_array(elements['c'])
            nsm = _cast_array(elements['nsm'])
            for eid, mid, nids, ai, ji, ci, nsmi in zip(eids, mids, nodes, A, J, c, nsm):
                model.add_conrod(eid, mid, nids, A=ai, j=ji, c=ci, nsm=nsmi, comment='')

        elif card_type == 'CBAR':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            g0 = _cast_array(elements['g0'])
            x = _cast_array(elements['x'])
            offt = _cast_array(elements['offt'])
            wa = _cast_array(elements['wa'])
            wb = _cast_array(elements['wb'])
            pa = _cast_array(elements['pa'])
            pb = _cast_array(elements['pb'])
            for eid, pid, nids, xi, g0i, offti, pai, pbi, wai, wbi in zip(
                    eids, pids, nodes, x, g0, offt, pa, pb, wa, wb):
                if g0i == -1:
                    g0i = None
                if xi[0] == np.nan:
                    xi = [None, None, None]
                model.add_cbar(eid, pid, nids, xi, g0i, offt=offti.decode(encoding),
                               pa=pai, pb=pbi, wa=wai, wb=wbi, comment='')

        elif card_type == 'CBEAM':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            g0 = _cast_array(elements['g0'])
            x = _cast_array(elements['x'])
            bit = _cast_array(elements['bit'])
            offt = _cast_array(elements['offt'])
            sa = _cast_array(elements['sa'])
            sb = _cast_array(elements['sb'])
            wa = _cast_array(elements['wa'])
            wb = _cast_array(elements['wb'])
            pa = _cast_array(elements['pa'])
            pb = _cast_array(elements['pb'])
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
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, pid, nids, (c1, c2) in zip(eids, pids, nodes, components):
                func(eid, pid, nids, c1=c1, c2=c2, comment='')
        elif card_type == 'CELAS2':
            eids = _cast_array(elements['eid'])
            k = _cast_array(elements['K'])
            ge = _cast_array(elements['ge'])
            s = _cast_array(elements['s'])
            nodes = _cast_array(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, ki, nids, (c1, c2), gei, si in zip(eids, k, nodes, components, ge, s):
                model.add_celas2(eid, ki, nids, c1=c1, c2=c2, ge=gei, s=si, comment='')
        elif card_type == 'CDAMP2':
            eids = _cast_array(elements['eid'])
            b = _cast_array(elements['B'])
            nodes = _cast_array(elements['nodes']).tolist()
            components = _cast(elements['components'])
            for eid, bi, nids, (c1, c2) in zip(eids, b, nodes, components):
                nids = list([nid if nid != 0 else None for nid in nids])
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
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_celas3(eid, pid, nids, comment='')
        elif card_type == 'CELAS4':
            eids = _cast_array(elements['eid'])
            k = _cast_array(elements['K'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, ki, nids in zip(eids, k, nodes):
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_celas4(eid, ki, nids, comment='')
        elif card_type == 'CDAMP4':
            eids = _cast_array(elements['eid'])
            b = _cast_array(elements['B'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, bi, nids in zip(eids, b, nodes):
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_cdamp4(eid, bi, nids, comment='')


        elif card_type == 'CBUSH':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            g0 = _cast_array(elements['g0'])
            x = _cast_array(elements['x']).tolist()
            cid = _cast_array(elements['cid'])
            ocid = _cast_array(elements['ocid'])
            s = _cast_array(elements['s'])
            si = _cast_array(elements['si']).tolist()
            for eid, pid, nids, xi, g0i, cidi, s2, ocidi, si2 in zip(
                    eids, pids, nodes, x, g0, cid, s, ocid, si):
                nids = list([nid if nid != 0 else None for nid in nids])
                if g0i == -1:
                    g0i = None
                #if xi[0] == np.nan:
                    #xi = [None, None, None]
                if cidi == -1:
                    cidi = None

                if si2[0] == np.nan:
                    si2 = [None, None, None]
                elem = model.add_cbush(eid, pid, nids, xi, g0i, cid=cidi, s=s2, ocid=ocidi, si=si2,
                                       comment='')
                write_card(elem)

        elif card_type == 'CGAP':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            g0 = _cast_array(elements['g0'])
            x = _cast_array(elements['x']).tolist()
            cid = _cast_array(elements['cid'])
            for eid, pid, nids, xi, g0i, cidi in zip(
                    eids, pids, nodes, x, g0, cid):
                nids = list([nid if nid != 0 else None for nid in nids])
                if g0i == -1:
                    g0i = None
                #if xi[0] == np.nan:
                    #xi = [None, None, None]
                if cidi == -1:
                    cidi = None
                elem = model.add_cgap(eid, pid, nids, xi, g0i, cid=cidi, comment='')
                #write_card(elem)

        elif card_type == 'CBUSH1D':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            cid = _cast_array(elements['cid'])
            for eid, pid, nids, cidi in zip(eids, pids, nodes, cid):
                nids = list([nid if nid != 0 else None for nid in nids])
                if cidi == -1:
                    cidi = None
                model.add_cbush1d(eid, pid, nids, cid=cidi, comment='')
        elif card_type == 'CBUSH2D':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            cid = _cast_array(elements['cid'])
            sptid = _cast_array(elements['sptid'])
            plane = _cast_array(elements['plane']).tolist()
            for eid, pid, nids, cidi, planei, sptidi in zip(eids, pids, nodes, cid, plane, sptid):
                planei = planei.decode(encoding)
                model.add_cbush2d(eid, pid, nids, cid=cidi, plane=planei, sptid=sptidi, comment='')

        elif card_type in ['CTRIA3', 'CTRIAR']:
            func = model.add_ctria3 if card_type == 'CTRIA3' else model.add_ctriar
            # TODO: doesn't support tflag
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            zoffsets = _cast_array(elements['zoffset'])
            nodes = _cast_array(elements['nodes']).tolist()
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
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            zoffsets = _cast_array(elements['zoffset'])
            nodes = _cast_array(elements['nodes']).tolist()
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
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            zoffsets = _cast_array(elements['zoffset'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_ctria6(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, comment='')
        elif card_type == 'CQUAD8':
            # TODO: doesn't support tflag
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            zoffsets = _cast_array(elements['zoffset'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta, zoffset in zip(
                    eids, pids, nodes, mcids, thetas, zoffsets):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_cquad8(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                                 tflag=0, T1=None, T2=None, T3=None, T4=None, comment='')

        elif card_type == 'CQUAD':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta in zip(eids, pids, nodes, mcids, thetas):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_cquad(eid, pid, nids, theta_mcid=theta_mcid, comment='')

        elif card_type == 'CSHEAR':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids in zip(eids, pids, nodes):
                model.add_cshear(eid, pid, nids, comment='')

        elif card_type == 'CTRIAX':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, mcid, theta in zip(eids, pids, nodes, mcids, thetas):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                model.add_ctriax(eid, pid, nids, theta_mcid=theta_mcid, comment='')
        elif card_type in ['CTRAX3', 'CTRAX6']:
            if card_type == 'CTRAX3':
                func = model.add_ctrax3
            else:
                func = model.add_ctrax6
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')
        elif card_type == 'CTRIAX6':
            eids = _cast_array(elements['eid'])
            mids = _cast_array(elements['mid'])
            thetas = _cast_array(elements['theta'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, mid, nids, theta in zip(eids, mids, nodes, thetas):
                nids = list([nid if nid != 0 else None for nid in nids])
                model.add_ctriax6(eid, mid, nids, theta=theta, comment='')

        elif card_type == 'CQUADX':
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            mcids = _cast_array(elements['mcid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, theta, mcid in zip(eids, pids, nodes, thetas, mcids):
                if mcid == -1:
                    theta_mcid = theta
                else:
                    theta_mcid = mcid
                nids = [None if nid == 0 else nid
                        for nid in nids]
                model.add_cquadx(eid, pid, nids, theta_mcid=theta_mcid, comment='')

        elif card_type in ['CQUADX4', 'CQUADX8']:
            if card_type == 'CQUADX4':
                func = model.add_cquadx4
            else:
                func = model.add_cquadx8
            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')

        elif card_type in ['CPLSTN3', 'CPLSTN4',
                           'CPLSTS3', 'CPLSTS4']:
            func_map = {
                'CPLSTN3' : model.add_cplstn3,
                'CPLSTN4' : model.add_cplstn4,
                'CPLSTS3' : model.add_cplsts3,
                'CPLSTS4' : model.add_cplsts4,
            }
            func = func_map[card_type]

            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')
        elif card_type in ['CPLSTN6', 'CPLSTN8',
                           'CPLSTS6', 'CPLSTS8']:
            func_map = {
                'CPLSTN6' : model.add_cplstn6,
                'CPLSTN8' : model.add_cplstn8,
                'CPLSTS6' : model.add_cplsts6,
                'CPLSTS8' : model.add_cplsts8,
            }
            func = func_map[card_type]

            eids = _cast_array(elements['eid'])
            pids = _cast_array(elements['pid'])
            thetas = _cast_array(elements['theta'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, pid, nids, theta in zip(eids, pids, nodes, thetas):
                func(eid, pid, nids, theta=theta, comment='')
        else:
            eids, values = load_cards_from_keys_values(
                'elements/%s' % card_type,
                elements, encoding, model.log)
            _put_keys_values_into_dict(model, 'elements', eids, values)
            #model.add_cdamp4(eid, b, nids, comment='')
            #model.add_cbush2d(eid, pid, nids, cid=0, plane='XY', sptid=None, comment='')
            #model.add_cfast(eid, pid, Type, ida, idb, gs=None, ga=None, gb=None,
                            #xs=None, ys=None, zs=None, comment='')
            #model.add_cmass1(eid, pid, nids, c1=0, c2=0, comment='')
            #model.add_cmass2(eid, mass, nids, c1, c2, comment='')
            #model.add_cmass3(eid, pid, nids, comment='')
            #model.add_cmass4(eid, mass, nids, comment='')
            #model.log.debug(card_type)
        model.card_count[card_type] = len(eids)

def hdf5_load_plotels(model, elements_group, unused_encoding):
    """loads the plotels from an HDF5 file"""
    for card_type in elements_group.keys():
        elements = elements_group[card_type]
        if card_type == 'PLOTEL':
            eids = _cast_array(elements['eid'])
            nodes = _cast_array(elements['nodes']).tolist()
            for eid, nids in zip(eids, nodes):
                model.add_plotel(eid, nids, comment='')
        else:  # pragma: no cover
            raise RuntimeError('card_type=%s in hdf5_load_plotels' % card_type)
        model.card_count[card_type] = len(eids)

def write_card(elem):  # pragma: no cover
    """verifies that the card was built correctly near where the card was made"""
    try:
        elem.write_card(size=8, is_double=False)
    except RuntimeError:
        elem.write_card(size=16, is_double=False)
    except Exception:  # pragma: no cover
        print(elem.get_stats())
        raise
