from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import cast_encoding_strip, cast_encoding_strip0, get_group_name
if TYPE_CHECKING: # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def _load_h5_parameter_pvt(model: BDF, group: Group, encoding: str):
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            group_name = get_group_name(group)
            raise RuntimeError(group_name)
            #msg += f'{indent}Group: {group_name}\n'
            #indent += '    '
            #print_group(group, indent=indent)
            #continue
        elif class_id == 'TABLE':
            pass
        else: # pragma: no cover
            raise RuntimeError(class_id)
        datatype = h5_element.name
        data = h5_element.read()
        name_bytes = data['NAME']
        name = cast_encoding_strip(name_bytes, encoding)
        if datatype == 'INT':
            # NAME, VALUE
            value = data['VALUE']
            for key, valuei in zip(name, value):
                model.add_param(key, [valuei])
        elif datatype == 'CHAR':
            # NAME, VALUE
            value_bytes = data['VALUE']
            value = cast_encoding_strip(value_bytes, encoding)
            for key, valuei in zip(name, value):
                model.add_param(key, [valuei])
        elif datatype == 'DOUBLE':
            # NAME, VALUE
            value = data['VALUE']
            for key, valuei in zip(name, value):
                model.add_param(key, [valuei])
        else: # pragma: no cover
            raise NotImplementedError(datatype)
        x = 1

def _load_h5_parameter_group(model: BDF, group: Group, encoding: str):
    group_name = get_group_name(group)
    skip_names = {'CASECC', 'TSTEP'} # , 'PVT',
    if group_name in skip_names:
        model.log.warning(f'skipping {group_name} in _load_h5_parameter_group')
    elif group_name == 'PVT':
        _load_h5_parameter_pvt(model, group, encoding)
    elif group_name == 'AECOMP':
        #('NAME', 'LISTTYPE', 'LISTIDS_POS', 'LISTIDS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        list_type = cast_encoding_strip(identity['LISTTYPE'], encoding)
        nlists = identity['LISTIDS_LEN']
        lists = group['LISTIDS'].read()
        model.aecomp._save(name, list_type, nlists, lists)
    elif group_name == 'AECOMPL':
        #('NAME', 'LABELS_POS', 'LABELS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        nlabels = identity['LABELS_LEN']
        labels_bytes = group['LABELS'].read()
        labels = cast_encoding_strip0(labels_bytes, encoding)
        model.aecompl._save(name, nlabels, labels)
    elif group_name == 'AEFACT':
        #('SID', 'DIS_POS', 'DIS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        aefact_id = identity['SID']
        nfractions = identity['DIS_LEN']
        fractions = group['DIS'].read()
        model.aefact._save(aefact_id, nfractions, fractions)
    else:
        raise NotImplementedError(group_name)

def load_h5_parameter(model: BDF, group: Group, encoding: str='latin1'):
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            _load_h5_parameter_group(model, h5_element, encoding)
            continue
        elif class_id == 'TABLE':
            pass
        else: # pragma: no cover
            raise RuntimeError(class_id)
        name = h5_element.name
        data = h5_element.read()
        #if name in {'PVT'}:
            #model.log.warning(f'skipping {element_name} in load_h5_parameter')
        if name == 'TSTEPNL':
            model.log.warning(f'skipping {name} in load_h5_parameter')
        elif name == 'AEROS':
            #('ACSID', 'RCSID', 'REFC', 'REFB', 'REFS', 'SYMXZ', 'SYMXY', 'DOMAIN_ID')
            acsid = data['ACSID'][0]
            rcsid = data['RCSID'][0]
            bref = data['REFB'][0]
            cref = data['REFC'][0]
            sref = data['REFS'][0]
            sym_xz = data['SYMXZ'][0]
            sym_xy = data['SYMXY'][0]
            aeros = model.add_aeros(cref, bref, sref, acsid=acsid,
                                    rcsid=rcsid, sym_xz=sym_xz, sym_xy=sym_xy, comment='')
            str(aeros)
        elif name == 'AESTAT':
            #('ID', 'LABEL', 'DOMAIN_ID')
            aestat_ids = data['ID']
            labels = data['LABEL']
            for aestat_id, label in zip(aestat_ids, labels):
                label_str = label.strip().decode('latin1')
                model.add_aestat(aestat_id, label_str, comment='')
        elif name == 'MDLPRM':
            #('NAME', 'VALUE', 'DOMAIN_ID')
            names = data['NAME']
            values = data['VALUE']
            mdlprm_dict = {}
            for name, value in zip(names, values):
                name_str = name.strip().decode('latin1')
                mdlprm_dict[name_str] = value
            assert len(mdlprm_dict) == 1, mdlprm_dict
            assert 'HDF5' in mdlprm_dict, mdlprm_dict
            model.add_mdlprm(mdlprm_dict, comment='')
        elif name == 'MONPNT1':
            name = data['NAME']
            label = data['LABEL']
            axes = data['AXES']
            comp = data['COMP']
            cp = data['CP']
            xyz = np.stack(
                [data['X'], data['Y'], data['Z']],
                axis=1)
            assert xyz.shape[1] == 3, xyz.shape
            cd = data['CD']
            model.monpnt1._save(
                name, label, comp, axes, cp, xyz, cd)
        else: # pragma: no cover
            raise NotImplementedError(name)
