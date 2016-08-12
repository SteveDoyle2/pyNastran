from __future__ import print_function
import os
import numpy as np
from six import string_types, iteritems
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.bdf import BDF


def bdf_merge(bdf_filenames, bdf_filename_out=None, renumber=True, encoding=None,
              size=8, is_double=False, cards_to_skip=None):
    """
    Merges multiple BDF into one file

    Parameters
    ----------
    bdf_filenames : List[str]
        list of bdf filenames
    bdf_filename_out : str / None
        the output bdf filename (default=None; None -> no writing)
    renumber : bool
        should the bdf be renumbered (default=True)
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write

    Supports
    --------
      nodes:      GRID
      coords:     CORDx
      elements:   CQUAD4, CTRIA3, CTETRA, CPENTA, CHEXA, CELASx, CBAR, CBEAM
                  CONM1, CONM2, CMASS
      properties: PSHELL, PCOMP, PSOLID, PMASS
      materials:  MAT1, MAT8

    .. todo:: doesn't support SPOINTs/EPOINTs
    .. warning:: still very preliminary
    """
    if not isinstance(bdf_filenames, (list, tuple)):
        raise TypeError('bdf_filenames is not a list/tuple...%s' % str(bdf_filenames))

    if not len(bdf_filenames) > 1:
        raise RuntimeError("You can't merge one BDF...bdf_filenames=%s" % str(bdf_filenames))
    for bdf_filename in bdf_filenames:
        if not isinstance(bdf_filename, string_types):
            raise TypeError('bdf_filenames is not a string...%s' % bdf_filename)
        #bdf_filenames = [bdf_filenames]

    #starting_id_dict_default = {
        #'cid' : max(model.coords.keys()),
        #'nid' : max(model.nodes.keys()),
        #'eid' : max([
            #max(model.elements.keys()),
            #max(model.masses.keys()),
        #]),
        #'pid' : max([
            #max(model.properties.keys()),
            #max(model.properties_mass.keys()),
        #]),
        #'mid' : max(model.material_ids),
    #}
    model = BDF(debug=False)
    model.disable_cards(cards_to_skip)
    bdf_filename0 = bdf_filenames[0]
    model.read_bdf(bdf_filename0, encoding=encoding)
    model.log.info('primary=%s' % bdf_filename0)

    data_members = [
        'coords', 'nodes', 'elements', 'masses', 'properties', 'properties_mass',
        'materials',
    ]
    for bdf_filename in bdf_filenames[1:]:
        model.log.info('model.masses = %s' % model.masses)
        starting_id_dict = {
            'cid' : max(model.coords.keys()) + 1,
            'nid' : max(model.nodes.keys()) + 1,
            'eid' : max([
                max(model.elements.keys()),
                0 if len(model.masses) == 0 else max(model.masses.keys()),
            ]) + 1,
            'pid' : max([
                max(model.properties.keys()),
                0 if len(model.properties_mass) == 0 else max(model.properties_mass.keys()),
            ]) + 1,
            'mid' : max(model.material_ids) + 1,
        }
        #for param, val in sorted(iteritems(starting_id_dict)):
            #print('  %-3s %s' % (param, val))

        model.log.info('secondary=%s' % bdf_filename)
        model2 = BDF(debug=False)
        model2.disable_cards(cards_to_skip)
        bdf_dump = 'bdf_merge_temp.bdf'
        #model2.read_bdf(bdf_filename, xref=False)

        bdf_renumber(bdf_filename, bdf_dump, starting_id_dict=starting_id_dict,
                     size=size, is_double=is_double, cards_to_skip=cards_to_skip)
        model2 = BDF(debug=False)
        model2.disable_cards(cards_to_skip)
        model2.read_bdf(bdf_dump)
        os.remove(bdf_dump)

        model.log.info('model2.node_ids = %s' % np.array(model2.node_ids))
        for data_member in data_members:
            data1 = getattr(model, data_member)
            data2 = getattr(model2, data_member)
            if isinstance(data1, dict):
                model.log.info('  working on %s' % (data_member))
                for key, value in iteritems(data2):
                    if data_member in 'coords' and key == 0:
                        continue
                    if isinstance(value, list):
                        raise NotImplementedError(type(value))
                    else:
                        assert key not in data1, key
                        data1[key] = value
                        #print('   %s' % key)
            else:
                raise NotImplementedError(type(data1))
    #if bdf_filenames_out:
        #model.write_bdf(bdf_filenames_out, size=size)

    if renumber:
        model.log.info('final renumber...')
        starting_id_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(model, bdf_filename_out, starting_id_dict=starting_id_dict,
                     size=size, is_double=is_double, cards_to_skip=cards_to_skip)
    elif bdf_filename_out:
        model.write_bdf(out_filename=bdf_filename_out, encoding=None,
                        size=size, is_double=is_double,
                        interspersed=True,
                        enddata=None)
    return model
