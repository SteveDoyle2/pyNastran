from __future__ import annotations
import sys
from typing import TextIO, Optional, Any, cast, TYPE_CHECKING
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.write_mesh_utils import (
    find_aero_location,
    write_dict, write_list,
    write_bdfs_dict, write_bdfs_list, write_bdfs_dict_list, write_xpoints_file,
    get_properties_by_element_type, _ifile,
)
from pyNastran.bdf.cards.nodes import write_xpoints
from pyNastran.bdf.bdf_interface.utils import sorteddict
try:
    from natsort import natsorted
except ModuleNotFoundError:
    natsorted = sorted

if TYPE_CHECKING:
    from pathlib import PathLike
    from pyNastran.bdf.bdf import BDF, DESVAR
    from io import StringIO
    TextFile = StringIO | TextIO


class Writer:
    def __init__(self, model: BDF):
        self.model = model

    #def __deepcopy__(self, memo):
    #    # Return the original instance unmodified/uncopied
    #    # TODO: update self.model...
    #    return self

    def write_h5(self, h5file, nastran_group):
        pass

    def get_encoding(self, encoding: Optional[str]=None) -> str:
        if encoding is not None:
            pass
        else:
            model = self.model
            encoding = model._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
            elif isinstance(encoding, bytes):
                # needed for hdf5 loader for some reason...
                encoding = encoding.decode('latin1')
                model._encoding = encoding
        encoding = cast(str, encoding)  # just for typing
        assert isinstance(encoding, str), encoding
        return encoding

    def get_long_ids(self, size: int) -> tuple[bool, int]:
        """calculate is_long_ids"""
        # required for MasterModelTaxi
        model = self.model
        is_long_ids = (
            any((len(dicti) and max(dicti) > 100000000 for dicti in (
                model.nodes, model.coords, model.elements, model.properties,
                model.materials, model.thermal_materials,
                #self.loads,
                model.load_combinations, model.masses,
                model.nsms, model.nsmadds,
            ))))
        if not is_long_ids:
            for loads in model.loads.values():
                for load in loads:
                    if load.type == 'TEMP':
                        if max(load.temperatures) > 100000000:
                            is_long_ids = True
                            break
                        #print(load.get_stats())
        if is_long_ids:
            size = 16
        return is_long_ids, size

    def write_mesh_long_ids_size(self, size: int, is_long_ids: bool) -> tuple[int, bool]:
        """helper method"""
        if is_long_ids and size == 16 or is_long_ids is False:
            return size, is_long_ids

        if size == 16 and is_long_ids is None or self.model.is_long_ids:
            size = 16
            is_long_ids = True
        else:
            is_long_ids = False
        return size, is_long_ids

    def write_bulk_data(self, bdf_file,
                        size: int=8, is_double: bool=False,
                        interspersed: bool=False,
                        enddata: Optional[bool]=None, close: bool=True,
                        coords_size: Optional[int]=None,
                        nodes_size: Optional[int]=None,
                        elements_size: Optional[int]=None,
                        loads_size: Optional[int]=None,
                        table_size: Optional[int]=None,
                        flfact_size: Optional[int]=None,
                        sort_cards: bool=True,
                        is_long_ids: bool=False,
                        is_csv: bool=False) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        bdf_file : varies
            file       - a file object
            StringIO() - a StringIO object
        size : int; {8, 16}
            the field size
        is_double : bool; default=False
            False : small field
            True : large field
        interspersed : bool; default=True
            Writes a bdf with properties & elements
            interspersed like how Patran writes the bdf.  This takes
            slightly longer than if interspersed=False, but makes it
            much easier to compare to a Patran-formatted bdf and is
            more clear.
        enddata : bool; default=None
            bool - enable/disable writing ENDDATA
            None - depends on input BDF
        sort_cards : bool; default=True
            sort the nodes, elements, ... to make finding things easier
        close : bool; default=True
            should the output file be closed

        .. note:: is_long_ids is only needed if you have ids longer
                  than 8 characters. It's an internal parameter, but if
                  you're calling the new sub-function, you might need
                  it.  Chances are you won't.
        """
        model = self.model
        size, coords_size, nodes_size, elements_size, loads_size, table_size, flfact_size = _fix_sizes(
            size, coords_size, nodes_size, elements_size, loads_size, table_size, flfact_size)

        args = {
            'is_double': is_double,
            'is_long_ids': is_long_ids,
            # 'is_csv': is_csv,
            'sort_cards': sort_cards,
        }
        self.write_params(bdf_file, size, **args)
        self.write_model_groups(bdf_file, sort_cards=sort_cards)
        self.write_nodes(bdf_file, nodes_size, is_csv=is_csv, **args)

        if interspersed:
            self.write_elements_interspersed(bdf_file, elements_size, **args)
        else:
            self.write_elements(bdf_file, elements_size, **args)
            self.write_properties(bdf_file, size, **args)
            #self.write_properties_by_element_type(bdf_file, size, is_double, is_long_ids, sort_cards=sort_cards)

        for cards in (model.bolt, model.boltseq, model.boltfor, model.boltfrc, model.boltld):
            for key, card in sorteddict(cards, sort_cards):
                bdf_file.write(card.write_card(size, is_double))

        self.write_materials(bdf_file, size, **args)
        self.write_masses(bdf_file, size, **args)

        # split out for write_bdf_symmetric
        self.write_rigid_elements(bdf_file, size, is_csv=is_csv, **args)
        self.write_aero(bdf_file, size, **args)

        self.write_common(bdf_file, loads_size,
                          coords_size=coords_size,
                          table_size=table_size, flfact_size=flfact_size,
                          is_csv=is_csv, **args)
        if (enddata is None and 'ENDDATA' in model.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def write_header(self, bdf_file: TextIO, encoding: str,
                     write_header: bool=True) -> None:
        """Writes the executive and case control decks."""
        model = self.model
        model._set_punch()

        if model.nastran_format and write_header:
            bdf_file.write(f'$pyNastran: version={model.nastran_format}\n')
            bdf_file.write(f'$pyNastran: punch={model.punch}\n')
            bdf_file.write(f'$pyNastran: encoding={encoding}\n')
            #bdf_file.write(f'$pyNastran: nnodes={len(model.nodes):d}\n')
            #bdf_file.write(f'$pyNastran: nelements={len(model.elements):d}\n')

        if not model.punch:
            self.write_executive_control_deck(bdf_file)
            self.write_case_control_deck(bdf_file)

    def write_executive_control_deck(self, bdf_file: TextIO) -> None:
        """Writes the executive control deck."""
        msg = ''
        model = self.model
        for line in model.system_command_lines:
            msg += line + '\n'

        if model.executive_control_lines:
            msg += '$EXECUTIVE CONTROL DECK\n'

            if model.sol_iline is not None:
                if model.sol == 600 and model.sol_method:
                    new_sol = f'SOL 600,{model.sol_method}'
                else:
                    new_sol = f'SOL {model.sol}'
                model.executive_control_lines[model.sol_iline] = new_sol

            for line in model.executive_control_lines:
                msg += line + '\n'
            bdf_file.write(msg)
            if 'CEND' not in msg:
                bdf_file.write('CEND\n')

    def write_case_control_deck(self, bdf_file: TextIO) -> None:
        """Writes the Case Control Deck."""
        model = self.model
        if model.case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            if model.superelement_models:
                msg += model.case_control_deck.write(write_begin_bulk=False)
            else:
                msg += str(model.case_control_deck)
                assert 'BEGIN BULK' in msg, msg
            bdf_file.write(''.join(msg))
        else:
            # if you run:
            #   model.BDF()
            #   model.sol = 101
            #   ... # add stuff
            #   model.write_bdf(...)
            #
            #  without this line, you'll get a CEND, but not BEGIN BULK
            bdf_file.write('BEGIN BULK\n')

    def write_params(self, bdf_file: TextFile,
                     size: int=8, is_double: bool=False,
                     sort_cards: bool=True,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the PARAM cards"""
        size, is_long_ids = self.write_mesh_long_ids_size(
            size, is_long_ids)

        model = self.model
        if model.params or model.dti or model.mdlprm:
            bdf_file.write('$PARAMS\n')
            is_csv = False
            write_dict(bdf_file, model.dti, size, is_double, is_csv, is_long_ids, sort_cards)
            write_dict(bdf_file, model.params, size, is_double, is_csv, is_long_ids, sort_cards)
            # for unused_name, dti in sorteddict(model.dti, sort_cards):
            #     bdf_file.write(dti.write_card(size=size, is_double=is_double))
            # for (unused_key, param) in sorteddict(model.params, sort_cards):
            #     bdf_file.write(param.write_card(size, is_double))
            if model.mdlprm:
                bdf_file.write(model.mdlprm.write_card(size, is_double))

    def write_model_groups(self, bdf_file: TextFile,
                           sort_cards: bool=True) -> None:
        if self.model.model_groups:
            #bdf_file.write('$ MODELGROUPS\n')
            for group in self.model.model_groups.values():
                #bdf_file.write(f'$ {group}\n')
                print(group)
            #x = 1

    def write_nodes(self, bdf_file: TextFile,
                    size: int=8, is_double: bool=False,
                    sort_cards: bool=True,
                    is_long_ids: Optional[bool]=None,
                    is_csv: bool=False) -> None:
        """Writes the NODE-type cards"""
        model = self.model
        if model.spoints:
            bdf_file.write('$SPOINTS\n')
            bdf_file.write(write_xpoints('SPOINT', model.spoints))
        if model.epoints:
            bdf_file.write('$EPOINTS\n')
            bdf_file.write(write_xpoints('EPOINT', model.epoints))
        if model.points:
            bdf_file.write('$POINTS\n')
            for unused_point_id, point in sorteddict(model.points, sort_cards):
                bdf_file.write(point.write_card(size, is_double))

        if model.cyax:
            bdf_file.write(model.cyax.write_card(size, is_double))

        self.write_grids(
            bdf_file, size=size, is_double=is_double,
            sort_cards=sort_cards, is_csv=is_csv)
        if model.seqgp:
            bdf_file.write(model.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def write_grids(self, bdf_file: TextFile,
                    size: int=8, is_double: bool=False,
                    sort_cards: bool=True,
                    is_long_ids: Optional[bool]=None,
                    write_as_cid0: bool=False,
                    is_csv: bool=False) -> None:
        """Writes the GRID-type cards"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.nodes or model.grdset:
            bdf_file.write('$NODES\n')
            if model.grdset:
                bdf_file.write(model.grdset.write_card(size))

            if write_as_cid0:
                for (unused_nid, node) in sorteddict(model.nodes, sort_cards):
                    if node.cp != 0:
                        xyz = node.get_position()
                        node.uncross_reference()
                        node.cp = 0
                        node.xyz = xyz
                        bdf_file.write(node.write_card(size, is_double))
                    else:
                        bdf_file.write(node.write_card(size, is_double))
            else:
                write_dict(bdf_file, model.nodes, size, is_double,
                           is_csv, is_long_ids, sort_cards)

    def write_coords(self, bdf_file: TextFile,
                     size: int=8, is_double: bool=False,
                     sort_cards: bool=True,
                     is_long_ids: Optional[bool]=None,
                     breakout_coords: bool=False) -> None:
        """Writes the coordinate cards in a sorted order"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if len(model.coords) > 1:
            bdf_file.write('$COORDS\n')

        if breakout_coords:
            for (coord_id, coord) in sorteddict(model.coords, sort_cards):
                if coord_id == 0:
                    continue
                if coord.origin is not None:
                    bdf_file.write(f'$ cid={coord_id}\n')
                    bdf_file.write(f'$   origin={coord.origin}\n')
                    bdf_file.write(f'$   i={coord.i}\n')
                    bdf_file.write(f'$   j={coord.j}\n')
                    bdf_file.write(f'$   k={coord.k}\n')
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card_16(is_double))
        else:
            for (coord_id, coord) in sorteddict(model.coords, sort_cards):
                if coord_id == 0:
                    continue
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card_16(is_double))

    def write_matcids(self, bdf_file: TextFile,
                      size: int=8, is_double: bool=False,
                      sort_cards: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the MATCID cards in a sorted order"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)

        if len(model.matcid):
            bdf_file.write('$MATCID\n')
        for (cid, matcids) in sorteddict(model.matcid, sort_cards):
            for matcid in matcids:
                try:
                    bdf_file.write(matcid.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(matcid.write_card_16(is_double))

    def write_elements(self, bdf_file: TextIO,
                       size: int=8, is_double: bool=False,
                       sort_cards: bool=True,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.elements:
            bdf_file.write('$ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorteddict(model.elements, sort_cards):
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorteddict(model.elements, sort_cards):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
        if model.ao_element_flags:
            for (eid, element) in sorteddict(model.ao_element_flags, sort_cards):
                bdf_file.write(element.write_card(size, is_double))
        if model.normals:
            for (unused_nid, snorm) in sorteddict(model.normals, sort_cards):
                bdf_file.write(snorm.write_card(size, is_double))
        self.write_nsm(bdf_file, size, is_double, sort_cards)

    def write_elements_interspersed(self, bdf_file: TextIO,
                                    size: int=8, is_double: bool=False,
                                    sort_cards: bool=True,
                                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements and properties in and interspersed order"""
        model = self.model
        # from pyNastran.bdf.bdf import BDF
        # assert isinstance(model, BDF), type(model)
        # model.log
        missing_properties = []
        if model.properties:
            bdf_file.write('$ELEMENTS_WITH_PROPERTIES\n')

        eids_written: list[int] = []
        pids = sorted(model.properties.keys())
        pid_eids = model.get_element_ids_dict_with_pids(
            pids, stop_if_no_eids=False)

        #failed_element_types = set()
        for (pid, eids) in sorted(pid_eids.items()):
            prop = model.properties[pid]
            if eids:
                bdf_file.write(prop.write_card(size, is_double))
                eids.sort()
                for eid in eids:
                    element = model.elements[eid]
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type!r} eid={eid}')
                        raise
                eids_written += eids
            else:
                missing_properties.append(prop.write_card(size, is_double))

        eids_missing = set(model.elements.keys()).difference(set(eids_written))
        if eids_missing:
            bdf_file.write('$ELEMENTS_WITH_NO_PROPERTIES '
                           '(PID=0 and unanalyzed properties)\n')
            for eid in sorted(eids_missing):
                element = model.elements[eid]
                try:
                    bdf_file.write(element.write_card(size, is_double))
                except Exception:
                    print(f'failed printing element...type={element.type} eid={eid}')
                    raise

        if missing_properties or model.pdampt or model.pbusht or model.pelast:
            bdf_file.write('$UNASSOCIATED_PROPERTIES\n')
            for pid, card in sorteddict(model.pbusht, sort_cards):
                bdf_file.write(card.write_card(size, is_double))
            for pid, card in sorteddict(model.pdampt, sort_cards):
                bdf_file.write(card.write_card(size, is_double))
            for pid, card in sorteddict(model.pelast, sort_cards):
                bdf_file.write(card.write_card(size, is_double))
            for card in missing_properties:
                # this is a string...
                #print("missing_property = ", card
                bdf_file.write(card)

        if model.ao_element_flags:
            for (eid, element) in sorted(model.ao_element_flags.items()):
                bdf_file.write(element.write_card(size, is_double))
        if model.normals:
            for (unused_nid, snorm) in sorted(model.normals.items()):
                bdf_file.write(snorm.write_card(size, is_double))
        self.write_nsm(bdf_file, size, is_double)

    def write_nsm(self, bdf_file: TextIO,
                  size: int=8, is_double: bool=False,
                  sort_cards: bool=True,
                  is_long_ids: Optional[bool]=None) -> None:
        """Writes the nsm in a sorted order"""
        model = self.model
        if model.nsms or model.nsmadds:
            bdf_file.write('$NSM\n')
            for (unused_id, nsmadds) in sorteddict(model.nsmadds, sort_cards):
                for nsmadd in nsmadds:
                    bdf_file.write(str(nsmadd))
            for (key, nsms) in sorteddict(model.nsms, sort_cards):
                for nsm in nsms:
                    try:
                        bdf_file.write(nsm.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing nsm...type={nsm.type} key={key!r}')
                        raise

    def write_masses(self, bdf_file: TextFile,
                     size: int=8, is_double: bool=False,
                     sort_cards: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the mass cards sorted by ID"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.properties_mass:
            bdf_file.write('$PROPERTIES_MASS\n')
            for (pid, mass) in sorteddict(model.properties_mass, sort_cards):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except Exception:
                    print(f'failed printing mass property...type={mass.type} pid={pid}')
                    raise

        if model.masses:
            bdf_file.write('$MASSES\n')
            for (eid, mass) in sorteddict(model.masses, sort_cards):
                try:
                    bdf_file.write(mass.write_card(size, is_double))
                except Exception:
                    print(f'failed printing masses...type={mass.type} eid={eid}')
                    raise

    def write_properties(self, bdf_file: TextFile,
                         size: int=8, is_double: bool=False,
                         sort_cards: bool=True,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the properties in a sorted order"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        is_big_properties = hasattr(model, 'big_properties') and model.big_properties
        is_properties = (model.properties or model.pelast or
                         model.pdampt or model.pbusht or is_big_properties)
        is_csv = False
        if not is_properties:
            return
        bdf_file.write('$PROPERTIES\n')
        prop_groups = (model.properties, model.pelast, model.pdampt, model.pbusht)
        if is_long_ids:
            for prop_group in prop_groups:
                write_dict(bdf_file, prop_group, size, is_double, is_csv, is_long_ids, sort_cards)
                # for unused_pid, prop in sorteddict(prop_group, sort_cards):
                #     bdf_file.write(prop.write_card_16(is_double))
            #except Exception:
                #print('failed printing property type=%s' % prop.type)
                #raise
        else:
            for prop_group in prop_groups:
                write_dict(bdf_file, prop_group, size, is_double, is_csv, is_long_ids, sort_cards)
                # for unused_pid, prop in sorteddict(prop_group, sort_cards):
                #     bdf_file.write(prop.write_card(size, is_double))

        if is_big_properties:
            for unused_pid, prop in sorteddict(model.big_properties, sort_cards):
                bdf_file.write(prop.write_card_16(is_double))

    def write_properties_by_element_type(self, bdf_file: TextFile, size: int=8,
                                         is_double: bool=False,
                                         is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the properties in a sorted order by property type grouping

        TODO: Missing some property types.
        """
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        is_properties = model.properties or model.pelast or model.pdampt or model.pbusht
        if not is_properties:
            return

        out = get_properties_by_element_type(model)
        prop_class_to_property_types, prop_type_to_property_class, props_by_class = out

        bdf_file.write('$PROPERTIES\n')
        for prop_class, prop_types in prop_class_to_property_types.items():
            # print(prop_class, prop_types)
            # for prop_type in prop_types:
            #     if prop_type not in properties_by_class:
            #         continue
            #     print('  ', prop_type)
            props = props_by_class[prop_class]
            if not props:
                continue
            bdf_file.write('$' + '-' * 80 + '\n')
            bdf_file.write('$ %s\n' % prop_class)

            for prop in props:
                bdf_file.write(prop.write_card(size, is_double))
        bdf_file.write('$' + '-' * 80 + '\n')

    def write_materials(self, bdf_file: TextFile,
                        size: int=8, is_double: bool=False,
                        sort_cards: bool=True,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the materials in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        is_big_materials = hasattr(model, 'big_materials') and model.big_materials
        is_materials = (model.materials or model.hyperelastic_materials or model.creep_materials or
                        model.MATS1 or model.MATS3 or model.MATS8 or model.MATT1 or
                        model.MATT2 or model.MATT3 or model.MATT4 or model.MATT5 or
                        model.MATT8 or model.MATT9 or model.nxstrats or is_big_materials)
        if not is_materials:
            return
        bdf_file.write('$MATERIALS\n')
        for (unused_mid, material) in sorteddict(model.materials, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.hyperelastic_materials, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.creep_materials, sort_cards):
            bdf_file.write(material.write_card(size, is_double))

        for (unused_mid, material) in sorteddict(model.MATS1, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATS3, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATS8, sort_cards):
            bdf_file.write(material.write_card(size, is_double))

        for (unused_mid, material) in sorteddict(model.MATT1, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT2, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT3, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT4, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT5, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT8, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATT9, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_mid, material) in sorteddict(model.MATDMG, sort_cards):
            bdf_file.write(material.write_card(size, is_double))
        for (unused_sid, nxstrat) in sorteddict(model.nxstrats, sort_cards):
            bdf_file.write(nxstrat.write_card(size, is_double))

        if is_big_materials:
            for unused_mid, mat in sorteddict(model.big_materials, sort_cards):
                bdf_file.write(mat.write_card_16(is_double))

    def write_rigid_elements(self, bdf_file: TextFile,
                             size: int=8, is_double: bool=False,
                             sort_cards: bool=True,
                             is_long_ids: Optional[bool]=None,
                             is_csv: bool=False) -> None:
        """Writes the rigid elements in a sorted order"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.rigid_elements:
            bdf_file.write('$RIGID ELEMENTS\n')
            if is_long_ids:
                for (eid, element) in sorteddict(model.rigid_elements, sort_cards):
                    try:
                        bdf_file.write(element.write_card_16(is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
            else:
                for (eid, element) in sorteddict(model.rigid_elements, sort_cards):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing element...type={element.type} eid={eid}')
                        raise
        if model.plotels:
            bdf_file.write('$PLOT ELEMENTS\n')
            write_dict(bdf_file, model.plotels, size, is_double, is_csv, is_long_ids, sort_cards)

    def write_constraints(self, bdf_file: TextFile,
                          size: int=8, is_double: bool=False,
                          sort_cards: bool=True,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the constraint cards sorted by ID"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.suport or model.suport1:
            bdf_file.write('$CONSTRAINTS\n')
            for suport in model.suport:
                bdf_file.write(suport.write_card(size, is_double))
            for unused_suport_id, suport in sorteddict(model.suport1, sort_cards):
                bdf_file.write(suport.write_card(size, is_double))

        if model.spcs or model.spcadds or model.spcoffs:
            # bdf_file.write('$SPCs\n')
            # str_spc = str(self.spcObject) # old
            # if str_spc:
            #     bdf_file.write(str_spc)
            # else:
            bdf_file.write('$SPCs\n')
            for (unused_id, spcadds) in sorteddict(model.spcadds, sort_cards):
                for spcadd in spcadds:
                    bdf_file.write(str(spcadd))
            for (unused_id, spcs) in sorteddict(model.spcs, sort_cards):
                for spc in spcs:
                    bdf_file.write(str(spc))
            for (unused_id, spcoffs) in sorteddict(model.spcoffs, sort_cards):
                for spc in spcoffs:
                    bdf_file.write(str(spc))

        if model.mpcs or model.mpcadds:
            bdf_file.write('$MPCs\n')
            for (unused_id, mpcadds) in sorteddict(model.mpcadds, sort_cards):
                for mpcadd in mpcadds:
                    bdf_file.write(str(mpcadd))
            for (unused_id, mpcs) in sorteddict(model.mpcs, sort_cards):
                for mpc in mpcs:
                    bdf_file.write(mpc.write_card(size, is_double))

    def write_loads(self, bdf_file: TextFile,
                    size: int=8, is_double: bool=False,
                    is_csv: bool=False, sort_cards: bool=True,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.load_combinations or model.loads or model.tempds or model.cyjoin:
            bdf_file.write('$LOADS\n')
            for (key, load_combinations) in sorteddict(model.load_combinations, sort_cards):
                for load_combination in load_combinations:
                    try:
                        bdf_file.write(load_combination.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load_combination.type} key={key!r}')
                        raise

            if is_long_ids:
                for (key, loadcase) in sorteddict(model.loads, sort_cards):
                    for load in loadcase:
                        try:
                            bdf_file.write(load.write_card_16(is_double))
                        except Exception:
                            print(f'failed printing load...type={load.type} key={key!r}')
                            raise
            else:
                for (key, loadcase) in sorteddict(model.loads, sort_cards):
                    for load in loadcase:
                        # if load.type == 'PLOAD2':
                        #     try:
                        #         load.raw_fields()
                        #     except Exception:
                        #         bdf_file.write(load.write_card_separate(self, size, is_double))
                        #         continue
                        try:
                            bdf_file.write(load.write_card(size, is_double))
                        except Exception:
                            print(f'failed printing load...type={load.type} key={key!r}')
                            raise

            for unused_key, tempd in sorteddict(model.tempds, sort_cards):
                bdf_file.write(tempd.write_card(size, is_double))
            for unused_key, cyjoin in sorteddict(model.cyjoin, sort_cards):
                bdf_file.write(cyjoin.write_card(size, is_double))
        self.write_dloads(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids,
                           is_csv=is_csv, sort_cards=sort_cards)

    def write_dloads(self, bdf_file: TextFile,
                     size: int=8, is_double: bool=False,
                     is_csv: bool=False,
                     sort_cards: bool=True,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the dload cards sorted by ID"""
        model = self.model
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        if model.dloads or model.dload_entries:
            bdf_file.write('$DLOADS\n')
            for (key, loadcase) in sorteddict(model.dloads, sort_cards):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load.type} key={key!r}')
                        raise

            for (key, loadcase) in sorteddict(model.dload_entries, sort_cards):
                for load in loadcase:
                    try:
                        bdf_file.write(load.write_card(size, is_double))
                    except Exception:
                        print(f'failed printing load...type={load.type} key={key!r}')
                        raise

    def write_dynamic(self, bdf_file: TextFile,
                      size: int=8, is_double: bool=False,
                      sort_cards: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the dynamic cards sorted by ID"""
        model = self.model
        is_dynamic = (model.dareas or model.dphases or model.nlparms or model.frequencies or
                      model.methods or model.cMethods or model.tsteps or model.tstepnls or
                      model.transfer_functions or model.delays or model.rotors or model.tics or
                      model.nlpcis or model.acmodl is not None or
                      model.acplnw or model.amlreg or model.micpnt)

        if not is_dynamic:
            return
        bdf_file.write('$DYNAMIC\n')
        for (unused_id, method) in sorteddict(model.methods, sort_cards):
            bdf_file.write(method.write_card(size, is_double))
        for (unused_id, cmethod) in sorteddict(model.cMethods, sort_cards):
            bdf_file.write(cmethod.write_card(size, is_double))
        for (unused_id, darea) in sorteddict(model.dareas, sort_cards):
            bdf_file.write(darea.write_card(size, is_double))
        for (unused_id, dphase) in sorteddict(model.dphases, sort_cards):
            bdf_file.write(dphase.write_card(size, is_double))
        for (unused_id, nlparm) in sorteddict(model.nlparms, sort_cards):
            bdf_file.write(nlparm.write_card(size, is_double))
        for (unused_id, nlpci) in sorteddict(model.nlpcis, sort_cards):
            bdf_file.write(nlpci.write_card(size, is_double))
        for (unused_id, tstep) in sorteddict(model.tsteps, sort_cards):
            bdf_file.write(tstep.write_card(size, is_double))
        for (unused_id, tstepnl) in sorteddict(model.tstepnls, sort_cards):
            bdf_file.write(tstepnl.write_card(size, is_double))
        for (unused_id, freqs) in sorteddict(model.frequencies, sort_cards):
            for freq in freqs:
                bdf_file.write(freq.write_card(size, is_double))
        for (unused_id, delay) in sorteddict(model.delays, sort_cards):
            bdf_file.write(delay.write_card(size, is_double))
        for (unused_id, rotor) in sorteddict(model.rotors, sort_cards):
            bdf_file.write(rotor.write_card(size, is_double))
        for (unused_id, tic) in sorteddict(model.tics, sort_cards):
            bdf_file.write(tic.write_card(size, is_double))

        for (unused_id, tfs) in sorteddict(model.transfer_functions, sort_cards):
            for transfer_function in tfs:
                bdf_file.write(transfer_function.write_card(size, is_double))
        if model.acmodl:
            bdf_file.write(model.acmodl.write_card(size, is_double))
        for key, acplnw in model.acplnw.items():
            bdf_file.write(acplnw.write_card(size, is_double))
        for key, amlreg in model.amlreg.items():
            bdf_file.write(amlreg.write_card(size, is_double))
        for key, micpnt in model.micpnt.items():
            bdf_file.write(micpnt.write_card(size, is_double))

    def write_aero(self, bdf_file: TextFile,
                   size: int=8, is_double: bool=False,
                   sort_cards: bool=True,
                   is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero cards"""
        model = self.model
        if model.caeros or model.paeros or model.monitor_points or model.splines:
            bdf_file.write('$AERO\n')
            for (unused_id, caero) in sorteddict(model.caeros, sort_cards):
                bdf_file.write(caero.write_card(size, is_double))
            for (unused_id, paero) in sorteddict(model.paeros, sort_cards):
                bdf_file.write(paero.write_card(size, is_double))
            for (unused_id, spline) in sorteddict(model.splines, sort_cards):
                bdf_file.write(spline.write_card(size, is_double))

        if model.monitor_points or model.group:
            for monitor_point in model.monitor_points:
                bdf_file.write(monitor_point.write_card(size, is_double))
            for (unused_id, group) in sorteddict(model.group, sort_cards):
                bdf_file.write(group.write_card(size, is_double))
        model.zaero.write_bdf(bdf_file, size=8, is_double=False)

    def write_aero_control(self, bdf_file: TextFile,
                           size: int=8, is_double: bool=False,
                           sort_cards: bool=True,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero control surface cards"""
        model = self.model
        is_aero = (model.aecomps or model.aefacts or model.aeparams or model.aelinks or
                   model.aelists or model.aestats or model.aesurf or model.aesurfs)
        if not is_aero:
            return
        bdf_file.write('$AERO CONTROL SURFACES\n')
        for (unused_id, aelinks) in sorteddict(model.aelinks, sort_cards):
            for aelink in aelinks:
                bdf_file.write(aelink.write_card(size, is_double))

        for (unused_id, aecomp) in sorteddict(model.aecomps, sort_cards):
            bdf_file.write(aecomp.write_card(size, is_double))
        for (unused_id, aeparam) in sorteddict(model.aeparams, sort_cards):
            bdf_file.write(aeparam.write_card(size, is_double))
        for (unused_id, aestat) in sorteddict(model.aestats, sort_cards):
            bdf_file.write(aestat.write_card(size, is_double))

        for (unused_id, aelist) in sorteddict(model.aelists, sort_cards):
            bdf_file.write(aelist.write_card(size, is_double))
        for (unused_id, aesurf) in sorteddict(model.aesurf, sort_cards):
            bdf_file.write(aesurf.write_card(size, is_double))
        for (unused_id, aesurfs) in sorteddict(model.aesurfs, sort_cards):
            bdf_file.write(aesurfs.write_card(size, is_double))
        for (unused_id, aefact) in sorteddict(model.aefacts, sort_cards):
            bdf_file.write(aefact.write_card(size, is_double))

    def write_static_aero(self, bdf_file: TextFile,
                          size: int=8, is_double: bool=False,
                          sort_cards: bool=True,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the static aero cards"""
        model = self.model
        is_aero = (
            model.aeros or model.trims or model.divergs or model.uxvec or
            model.aeforce or model.aepress or model.aedw)
        if not is_aero:
            return
        bdf_file.write('$STATIC AERO\n')
        # static aero
        if model.aeros:
            bdf_file.write(model.aeros.write_card(size, is_double))
        for (unused_id, trim) in sorteddict(model.trims, sort_cards):
            bdf_file.write(trim.write_card(size, is_double))
        for (unused_id, diverg) in sorteddict(model.divergs, sort_cards):
            bdf_file.write(diverg.write_card(size, is_double))
        for (unused_id, uxvec) in sorteddict(model.uxvec, sort_cards):
            bdf_file.write(uxvec.write_card(size, is_double))
        for aedw in model.aedw:
            bdf_file.write(aedw.write_card(size, is_double))
        for aeforce in model.aeforce:
            bdf_file.write(aeforce.write_card(size, is_double))
        for aepress in model.aepress:
            bdf_file.write(aepress.write_card(size, is_double))

    def write_flutter(self, bdf_file: TextFile, size: int=8,
                      flfact_size: int=8,
                      is_double: bool=False,
                      sort_cards: bool=True,
                      write_aero_in_flutter: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the flutter cards"""
        model = self.model
        if (write_aero_in_flutter and model.aero) or model.flfacts or model.flutters or model.mkaeros:
            bdf_file.write('$FLUTTER\n')
            if write_aero_in_flutter and model.aero is not None:
                bdf_file.write(model.aero.write_card(size, is_double))
            for (unused_id, flutter) in sorteddict(model.flutters, sort_cards):
                bdf_file.write(flutter.write_card(size, is_double))
            for (unused_id, flfact) in sorteddict(model.flfacts, sort_cards):
                bdf_file.write(flfact.write_card(flfact_size, is_double))
            for mkaero in model.mkaeros:
                bdf_file.write(mkaero.write_card(size, is_double))

    def write_gust(self, bdf_file: TextFile, size: int=8, is_double: bool=False,
                    write_aero_in_gust: bool=True,
                    sort_cards: bool=True,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the gust cards"""
        model = self.model
        if (write_aero_in_gust and model.aero) or model.gusts:
            bdf_file.write('$GUST\n')
            if write_aero_in_gust:
                if model.aero is not None:
                    bdf_file.write(model.aero.write_card(size, is_double))
            for (unused_id, gust) in sorteddict(model.gusts, sort_cards):
                bdf_file.write(gust.write_card(size, is_double))

    def write_sets(self, bdf_file: TextFile,
                   size: int=8, is_double: bool=False,
                   is_csv: bool=False,
                   sort_cards: bool=True,
                   is_long_ids: Optional[bool]=None) -> None:
        """Writes the SETx cards sorted by ID"""
        model = self.model
        is_sets = (model.sets or model.asets or model.omits or model.bsets
                   or model.csets or model.qsets or model.usets)
        if not is_sets:
            return
        bdf_file.write('$SETS\n')
        write_dict(bdf_file, model.sets, size, is_double, is_csv, is_long_ids, sort_cards)
        # for (unused_id, set_obj) in sorteddict(model.sets, sort_cards):  # dict
        #     bdf_file.write(set_obj.write_card(size, is_double))
        write_list(bdf_file, model.asets, size, is_double, is_csv, is_long_ids)
        write_list(bdf_file, model.omits, size, is_double, is_csv, is_long_ids)
        write_list(bdf_file, model.bsets, size, is_double, is_csv, is_long_ids)
        write_list(bdf_file, model.csets, size, is_double, is_csv, is_long_ids)
        write_list(bdf_file, model.qsets, size, is_double, is_csv, is_long_ids)

        for unused_name, usets in sorted(model.usets.items()):  # dict
            for set_obj in usets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))

    def write_dmigs(self, bdf_file: TextFile,
                    size: int=8, is_double: bool=False,
                    sort_cards: bool=True,
                    is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        model = self.model
        for (unused_name, dmig) in natsorted(model.dmig.items()):
            bdf_file.write(dmig.write_card(size, is_double))
        for (unused_name, dmi) in natsorted(model.dmi.items()):
            bdf_file.write(dmi.write_card(size, is_double))
        for (unused_name, dmij) in natsorted(model.dmij.items()):
            bdf_file.write(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in natsorted(model.dmiji.items()):
            bdf_file.write(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in natsorted(model.dmik.items()):
            bdf_file.write(dmik.write_card(size, is_double))
        for (unused_name, dmiax) in natsorted(model.dmiax.items()):
            bdf_file.write(dmiax.write_card(size, is_double))

    def write_contact(self, bdf_file: TextFile,
                      size: int=8, is_double: bool=False,
                      sort_cards: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the contact cards sorted by ID"""
        model = self.model
        is_contact = (model.bcrparas or model.bctadds or model.bctparas
                      or model.bctsets or model.bsurf or model.bsurfs
                      or model.bconp or model.blseg or model.bfric
                      or model.bgadds or model.bgsets or model.bctparms
                      or model.bcbodys or model.bcparas)
        if not is_contact:
            return
        bdf_file.write('$CONTACT\n')
        for (unused_id, bcbody) in sorteddict(model.bcbodys, sort_cards):
            bdf_file.write(bcbody.write_card(size, is_double))
        for (unused_id, bcpara) in sorteddict(model.bcparas, sort_cards):
            bdf_file.write(bcpara.write_card(size, is_double))

        for (unused_id, bcrpara) in sorteddict(model.bcrparas, sort_cards):
            bdf_file.write(bcrpara.write_card(size, is_double))
        for (unused_id, bctparam) in sorteddict(model.bctparms, sort_cards):
            bdf_file.write(bctparam.write_card(size, is_double))
        for (unused_id, bctadds) in sorteddict(model.bctadds, sort_cards):
            bdf_file.write(bctadds.write_card(size, is_double))
        for (unused_id, bctpara) in sorteddict(model.bctparas, sort_cards):
            bdf_file.write(bctpara.write_card(size, is_double))

        for (unused_id, bctset) in sorteddict(model.bctsets, sort_cards):
            bdf_file.write(bctset.write_card(size, is_double))
        for (unused_id, bsurfi) in sorteddict(model.bsurf, sort_cards):
            bdf_file.write(bsurfi.write_card(size, is_double))
        for (unused_id, bsurfsi) in sorteddict(model.bsurfs, sort_cards):
            bdf_file.write(bsurfsi.write_card(size, is_double))
        for (unused_id, bconp) in sorteddict(model.bconp, sort_cards):
            bdf_file.write(bconp.write_card(size, is_double))
        for (unused_id, blseg) in sorteddict(model.blseg, sort_cards):
            bdf_file.write(blseg.write_card(size, is_double))
        for (unused_id, bfric) in sorteddict(model.bfric, sort_cards):
            bdf_file.write(bfric.write_card(size, is_double))
        for (unused_id, bgadd) in sorteddict(model.bgadds, sort_cards):
            bdf_file.write(bgadd.write_card(size, is_double))
        for (unused_id, bgset) in sorteddict(model.bgsets, sort_cards):
            bdf_file.write(bgset.write_card(size, is_double))

    def write_superelements(self, bdf_file: TextFile,
                            size: int=8, is_double: bool=False,
                            sort_cards: bool=True,
                            is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the Superelement cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        model = self.model
        is_sets = (model.se_sets or model.se_bsets or model.se_csets
                   or model.se_qsets or model.se_usets)
        if is_sets:
            bdf_file.write('$SUPERELEMENTS\n')
            for set_obj in model.se_bsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in model.se_csets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for set_obj in model.se_qsets:  # list
                bdf_file.write(set_obj.write_card(size, is_double))
            for (unused_set_id, set_obj) in sorted(model.se_sets.items()):  # dict
                bdf_file.write(set_obj.write_card(size, is_double))
            for unused_name, usets in sorted(model.se_usets.items()):  # dict
                for set_obj in usets:  # list
                    bdf_file.write(set_obj.write_card(size, is_double))
            for suport in model.se_suport:  # list
                bdf_file.write(suport.write_card(size, is_double))

        for unused_seid, csuper in sorted(model.csuper.items()):
            bdf_file.write(csuper.write_card(size, is_double))
        for unused_seid, csupext in sorted(model.csupext.items()):
            bdf_file.write(csupext.write_card(size, is_double))

        for unused_seid, sebulk in sorted(model.sebulk.items()):
            bdf_file.write(sebulk.write_card(size, is_double))
        for unused_seid, seconct in sorted(model.seconct.items()):
            bdf_file.write(seconct.write_card(size, is_double))

        for unused_seid, sebndry in sorted(model.sebndry.items()):
            bdf_file.write(sebndry.write_card(size, is_double))
        for unused_seid, seelt in sorted(model.seelt.items()):
            bdf_file.write(seelt.write_card(size, is_double))
        for unused_seid, seexcld in sorted(model.seexcld.items()):
            bdf_file.write(seexcld.write_card(size, is_double))

        for unused_seid, selabel in sorted(model.selabel.items()):
            bdf_file.write(selabel.write_card(size, is_double))
        for unused_seid, seloc in sorted(model.seloc.items()):
            bdf_file.write(seloc.write_card(size, is_double))
        for unused_seid, seload in sorted(model.seload.items()):
            bdf_file.write(seload.write_card(size, is_double))
        for unused_seid, sempln in sorted(model.sempln.items()):
            bdf_file.write(sempln.write_card(size, is_double))
        for unused_setid, senqset in sorted(model.senqset.items()):
            bdf_file.write(senqset.write_card(size, is_double))
        for unused_seid, setree in sorted(model.setree.items()):
            bdf_file.write(setree.write_card(size, is_double))
        for unused_seid, release in sorted(model.release.items()):
            bdf_file.write(release.write_card(size, is_double))

    def write_rejects(self, bdf_file: TextFile, size: int=8,
                      is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the rejected (processed) cards and the rejected unprocessed
        cardlines

        """
        model = self.model
        if size == 8:
            print_func = print_card_8
        else:
            print_func = print_card_16

        if model.reject_cards:
            bdf_file.write('$REJECT_CARDS\n')
            for reject_card in model.reject_cards:
                try:
                    bdf_file.write(print_func(reject_card))
                except RuntimeError:
                    if len(reject_card) > 0:
                        line0 = reject_card[0].upper()
                        if line0.startswith('ADAPT'):
                            for line in reject_card:
                                assert isinstance(line, str), line
                                bdf_file.write(line+'\n')
                            continue
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if model.reject_lines:
            bdf_file.write('$REJECT_LINES\n')
            for reject_lines in model.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_file.write('%s\n' % reject2)
                elif isinstance(reject_lines, str):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_file.write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def write_parametric(self, bdf_file: TextFile, size: int=8, is_double: bool=False,
                         sort_cards: bool=True,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        model = self.model
        is_parametric = model.pset or model.pval or model.gmcurv or model.feedge or model.feface
        if not is_parametric:
            return
        for (unused_id, pset) in sorted(model.pset.items()):
            bdf_file.write(pset.write_card(size, is_double))
        for (unused_adapt_id, pvals) in sorted(model.pval.items()):
            for pval in pvals:
                bdf_file.write(pval.write_card(size, is_double))
        for (unused_id, gmcurv) in sorted(model.gmcurv.items()):
            bdf_file.write(gmcurv.write_card(size, is_double))
        for (unused_id, gmsurf) in sorted(model.gmsurf.items()):
            bdf_file.write(gmsurf.write_card(size, is_double))
        for (unused_id, feedge) in sorted(model.feedge.items()):
            bdf_file.write(feedge.write_card(size, is_double))
        for (unused_id, feface) in sorted(model.feface.items()):
            bdf_file.write(feface.write_card(size, is_double))

    def write_tables(self, bdf_file: TextFile,
                     size: int=8, is_double: bool=False,
                     sort_cards: bool=True,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the TABLEx cards sorted by ID"""
        model = self.model
        if model.tables or model.tables_d or model.tables_m or model.tables_sdamping:
            bdf_file.write('$TABLES\n')
            is_csv = False
            write_dict(bdf_file, model.tables, size, is_double, is_csv, is_long_ids, sort_cards)
            write_dict(bdf_file, model.tables_d, size, is_double, is_csv, is_long_ids, sort_cards)
            write_dict(bdf_file, model.tables_m, size, is_double, is_csv, is_long_ids, sort_cards)
            write_dict(bdf_file, model.tables_sdamping, size, is_double, is_csv, is_long_ids, sort_cards)
            # for (unused_id, table) in sorteddict(model.tables, sort_cards):
            #     bdf_file.write(table.write_card(size, is_double))
            # for (unused_id, table) in sorteddict(model.tables_d, sort_cards):
            #     bdf_file.write(table.write_card(size, is_double))
            # for (unused_id, table) in sorteddict(model.tables_m, sort_cards):
            #     bdf_file.write(table.write_card(size, is_double))
            # for (unused_id, table) in sorteddict(model.tables_sdamping, sort_cards):
            #     bdf_file.write(table.write_card(size, is_double))

        if model.random_tables:
            bdf_file.write('$RANDOM TABLES\n')
            for (unused_id, table) in sorteddict(model.random_tables, sort_cards):
                bdf_file.write(table.write_card(size, is_double))

    def write_optimization(self, bdf_file: TextFile,
                           size: int=8, is_double: bool=False,
                           sort_cards: bool=True,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        model = self.model
        is_optimization = (
            model.dconadds or model.dconstrs or model.desvars or model.ddvals or
            model.dresps or
            model.dvprels or model.dvmrels or model.dvcrels or model.doptprm or
            model.dlinks or model.dequations or model.dtable is not None or
            model.dvgrids or model.dscreen or model.topvar or model.modtrak or
            # nx optimization
            model.dvtrels or model.dmncon
        )
        if not is_optimization:
            return
        bdf_file.write('$OPTIMIZATION\n')
        for (unused_id, dconadd) in sorteddict(model.dconadds, sort_cards):
            bdf_file.write(dconadd.write_card(size, is_double))
        for (unused_id, dconstrs) in sorteddict(model.dconstrs, sort_cards):
            for dconstr in dconstrs:
                bdf_file.write(dconstr.write_card(size, is_double))
        for (unused_id, desvar) in sorteddict(model.desvars, sort_cards):
            bdf_file.write(desvar.write_card(size, is_double))
        for (unused_id, topvar) in sorteddict(model.topvar, sort_cards):
            bdf_file.write(topvar.write_card(size, is_double))
        for (unused_id, ddval) in sorteddict(model.ddvals, sort_cards):
            bdf_file.write(ddval.write_card(size, is_double))
        for (unused_id, dlink) in sorteddict(model.dlinks, sort_cards):
            bdf_file.write(dlink.write_card(size, is_double))
        for (unused_id, dresp) in sorteddict(model.dresps, sort_cards):
            bdf_file.write(dresp.write_card(size, is_double))

        for (unused_id, dvcrel) in sorteddict(model.dvcrels, sort_cards):
            bdf_file.write(dvcrel.write_card(size, is_double))
        for (unused_id, dvmrel) in sorteddict(model.dvmrels, sort_cards):
            bdf_file.write(dvmrel.write_card(size, is_double))
        for (unused_id, dvprel) in sorteddict(model.dvprels, sort_cards):
            bdf_file.write(dvprel.write_card(size, is_double))
        for (unused_id, dvgrids) in sorteddict(model.dvgrids, sort_cards):
            for dvgrid in dvgrids:
                bdf_file.write(dvgrid.write_card(size, is_double))
        for (unused_id, dscreen) in sorteddict(model.dscreen, sort_cards):
            bdf_file.write(str(dscreen))

        for (unused_id, equation) in sorteddict(model.dequations, sort_cards):
            bdf_file.write(str(equation))

        if model.dtable is not None:
            bdf_file.write(model.dtable.write_card(size, is_double))
        if model.doptprm is not None:
            bdf_file.write(model.doptprm.write_card(size, is_double))
        if model.modtrak is not None:
            bdf_file.write(model.modtrak.write_card(size, is_double))

        # nx optimization
        for (unused_id, dvtrel) in sorteddict(model.dvtrels, sort_cards):
            bdf_file.write(dvtrel.write_card(size, is_double))
        for (unused_id, dmncon) in sorteddict(model.dmncon, sort_cards):
            bdf_file.write(dmncon.write_card(size, is_double))

    def write_thermal(self, bdf_file: TextFile,
                      size: int=8, is_double: bool=False,
                      sort_cards: bool=True,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal cards"""
        # PHBDY
        model = self.model
        is_thermal = (model.phbdys or model.convection_properties or model.bcs or
                      model.views or model.view3ds or model.radset or model.radcavs)
        if not is_thermal:
            return
        bdf_file.write('$THERMAL\n')
        is_csv = False
        write_dict(bdf_file, model.phbdys, size, is_double, is_csv, is_long_ids, sort_cards)

        #for unused_key, prop in sorted(model.thermal_properties.items()):
        #    bdf_file.write(str(prop))
        write_dict(bdf_file, model.convection_properties, size, is_double, is_csv, is_long_ids, sort_cards)

        # BCs
        for (unused_key, bcs) in sorteddict(model.bcs, sort_cards):
            for boundary_condition in bcs:  # list
                bdf_file.write(boundary_condition.write_card(size, is_double))

        write_dict(bdf_file, model.views, size, is_double, is_csv, is_long_ids, sort_cards)
        write_dict(bdf_file, model.view3ds, size, is_double, is_csv, is_long_ids, sort_cards)
        if model.radset:
            bdf_file.write(model.radset.write_card(size, is_double))
        write_dict(bdf_file, model.radcavs, size, is_double, is_csv, is_long_ids, sort_cards)

    def write_thermal_materials(self, bdf_file: TextFile,
                                size: int=8, is_double: bool=False,
                                sort_cards: bool=True,
                                is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal materials in a sorted order"""
        model = self.model
        if model.thermal_materials:
            bdf_file.write('$THERMAL MATERIALS\n')
            is_csv = False
            write_dict(bdf_file, model.thermal_materials, size, is_double, is_csv, is_long_ids, sort_cards)

    def write_common(self, bdf_file: TextFile,
                     size: int=8,
                     coords_size: int=8,
                     table_size: int=8,
                     flfact_size: int=8,
                     is_double: bool=False,
                     is_csv: bool=False,
                     sort_cards: bool=True,
                     is_long_ids: Optional[bool]=None) -> None:
        """
        Write the common outputs so none get missed...

        Parameters
        ----------
        bdf_file : file
            the file object
        size : int; default=8
            the field width
        is_double : bool; default=False
            is this double precision
         coords_size: int; default=8
            takes priority over size
         table_size: int; default=8
            takes priority over size
         flfact_size: int; default=8
            takes priority over size

        """
        model = self.model
        args = {
            'is_double': is_double,
            'is_long_ids': is_long_ids,
            # 'is_csv': is_csv,
            'sort_cards': sort_cards,
        }
        self.write_dmigs(bdf_file, size, **args)
        self.write_loads(bdf_file, size, **args)
        self.write_dynamic(bdf_file, size, **args)
        self.write_aero_control(bdf_file, size, **args)
        self.write_static_aero(bdf_file, size, **args)

        write_aero_in_flutter, write_aero_in_gust = find_aero_location(model)
        self.write_flutter(bdf_file, size=size, flfact_size=flfact_size,
                            write_aero_in_flutter=write_aero_in_flutter,
                            **args)
        self.write_gust(bdf_file, size, is_double, write_aero_in_gust, is_long_ids=is_long_ids)

        self.write_thermal(bdf_file, size, **args)
        self.write_thermal_materials(bdf_file, size, **args)
        self.write_constraints(bdf_file, size, **args)
        self.write_optimization(bdf_file, size, **args)
        self.write_tables(bdf_file, table_size, **args)
        self.write_sets(bdf_file, size, is_csv=is_csv, **args)
        self.write_superelements(bdf_file, size, **args)
        self.write_contact(bdf_file, size, **args)
        self.write_parametric(bdf_file, size, **args)
        self.write_rejects(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self.write_coords(bdf_file, coords_size, **args)
        self.write_matcids(bdf_file, size, **args)

    # def write_nodes_associated(self, bdf_file, size=8, is_double=False):
    #     """
    #     Writes the NODE-type in associated and unassociated groups.
    #
    #     .. warning:: Sometimes crashes, probably on invalid BDFs.
    #     """
    #     associated_nodes = set()
    #     for (eid, element) in model.elements.items():
    #         associated_nodes = associated_nodes.union(set(element.node_ids))
    #
    #     all_nodes = set(model.nodes.keys())
    #     unassociated_nodes = list(all_nodes.difference(associated_nodes))
    #     #missing_nodes = all_nodes.difference(
    #
    #     # TODO: this really shouldn't be a list...???
    #     associated_nodes = list(associated_nodes)
    #
    #     if associated_nodes:
    #         bdf_file.write('$ASSOCIATED NODES\n')
    #         if model.grdset:
    #             bdf_file.write(model.grdset.write_card(size, is_double))
    #         # TODO: this really shouldn't be a dictionary...???
    #         for key, node in sorted(associated_nodes.items()):
    #             bdf_file.write(node.write_card(size, is_double))
    #
    #     if unassociated_nodes:
    #         bdf_file.write('$UNASSOCIATED NODES\n')
    #         if model.grdset and not associated_nodes:
    #             v(model.grdset.write_card(size, is_double))
    #         for key, node in sorted(unassociated_nodes.items()):
    #             if key in model.nodes:
    #                 bdf_file.write(node.write_card(size, is_double))
    #             else:
    #                 bdf_file.write('$ Missing NodeID=%s' % key)

    def write_elements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the elements in a sorted order
        """
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)

        model = self.model
        if model.elements:
            write_bdfs_dict(bdf_files, model.elements, size, is_double, is_long_ids)

        if model.ao_element_flags:
            write_bdfs_dict(bdf_files, model.ao_element_flags, size, is_double, is_long_ids)
        if model.normals:
            write_bdfs_dict(bdf_files, model.normals, size, is_double, is_long_ids)
        self.write_nsm_file(bdf_files, size, is_double)

    def write_nsm_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the nsm in a sorted order"""
        model = self.model
        if model.nsms or model.nsmadds:
            write_bdfs_dict_list(bdf_files, model.nsmadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.nsms, size, is_double, is_long_ids)

    def write_aero_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero cards"""
        model = self.model
        if model.caeros or model.paeros or model.monitor_points or model.splines:
            write_bdfs_dict(bdf_files, model.caeros, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.paeros, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.splines, size, is_double, is_long_ids)
            for monitor_point in model.monitor_points:
                bdf_files[monitor_point.ifile].write(monitor_point.write_card(size, is_double))
        model.zaero.write_bdf(bdf_files[0], size=8, is_double=False)

    def write_aero_control_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero control surface cards"""
        model = self.model
        is_aero = (
            model.aecomps or model.aefacts or model.aeparams or model.aelinks or
            model.aelists or model.aestats or model.aesurf or model.aesurfs)
        if is_aero:
            return
        write_bdfs_dict_list(bdf_files, model.aecomps, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.aecomps, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aeparams, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aestats, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aelists, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aesurf, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aesurfs, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.aefacts, size, is_double, is_long_ids)

    def write_static_aero_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                               is_long_ids: Optional[bool]=None) -> None:
        """Writes the static aero cards"""
        model = self.model
        if model.aeros or model.trims or model.divergs:
            # static aero
            if model.aeros:
                bdf_files[model.aeros.ifile].write(model.aeros.write_card(size, is_double))

            write_bdfs_dict(bdf_files, model.trims, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.divergs, size, is_double, is_long_ids)

    def write_flutter_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           write_aero_in_flutter: bool=True,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the flutter cards"""
        model = self.model
        if (write_aero_in_flutter and model.aero) or model.flfacts or model.flutters or model.mkaeros:
            if write_aero_in_flutter and model.aero is not None:
                file_obj = bdf_files[_ifile(model.aero)]
                file_obj.write(model.aero.write_card(size, is_double))
            write_bdfs_dict(bdf_files, model.flutters, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.flfacts, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, model.mkaeros, size, is_double, is_long_ids)

    def write_gust_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                        write_aero_in_gust: bool=True, is_long_ids: Optional[bool]=None) -> None:
        """Writes the gust cards"""
        model = self.model
        if (write_aero_in_gust and model.aero) or model.gusts:
            if write_aero_in_gust:
                for (unused_id, aero) in sorted(model.aero.items()):
                    bdf_files[_ifile(aero)].write(aero.write_card(size, is_double))
            write_bdfs_dict(bdf_files, model.gusts, size, is_double, is_long_ids)

    def write_common_file(self, bdf_files: TextFile,
                          size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """
        Write the common outputs so none get missed...

        Parameters
        ----------
        bdf_file : file
            the file object
        size : int (default=8)
            the field width
        is_double : bool (default=False)
            is this double precision

        """
        self.write_dmigs_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_loads_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_dynamic_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_aero_control_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_static_aero_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

        write_aero_in_flutter, write_aero_in_gust = find_aero_location(self.model)
        self.write_flutter_file(bdf_files, size, is_double, write_aero_in_flutter,
                                 is_long_ids=is_long_ids)
        self.write_gust_file(bdf_files, size, is_double, write_aero_in_gust,
                             is_long_ids=is_long_ids)

        self.write_thermal_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_thermal_materials_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_constraints_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_optimization_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_tables_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_sets_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_superelements_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_contact_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_rejects_file(bdf_files, size, is_double, is_long_ids=is_long_ids)
        self.write_coords_file(bdf_files, size, is_double, is_long_ids=is_long_ids)

    def write_constraints_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                               is_long_ids: Optional[bool]=None) -> None:
        """Writes the constraint cards sorted by ID"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.suport or model.suport1:
            for suport in model.suport:
                bdf_files[suport.ifile].write(suport.write_card(size, is_double))
            for unused_suport_id, suport in sorted(model.suport1.items()):
                bdf_files[suport.ifile].write(suport.write_card(size, is_double))

        if model.spcs or model.spcadds or model.spcoffs:
            # bdf_file.write('$SPCs\n')
            # str_spc = str(model.spcObject) # old
            # if str_spc:
            #     bdf_file.write(str_spc)
            # else:
            write_bdfs_dict_list(bdf_files, model.spcadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.spcs, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.spcoffs, size, is_double, is_long_ids)

        if model.mpcs or model.mpcadds:
            write_bdfs_dict_list(bdf_files, model.mpcadds, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.mpcs, size, is_double, is_long_ids)

    def write_contact_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the contact cards sorted by ID"""
        model = self.model
        is_contact = (model.bcrparas or model.bctadds or model.bctparas or model.bctparms
                      or model.bctsets or model.bsurf or model.bsurfs
                      or model.bconp or model.blseg or model.bfric
                      or model.bgsets or model.bgadds)
        if is_contact:
            return
        write_bdfs_dict(bdf_files, model.bcrparas, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bctadds, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bctparas, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bctparms, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bctsets, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bsurf, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bsurfs, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.bconp, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.blseg, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bfric, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bgadds, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.bgsets, size, is_double, is_long_ids)

    def write_coords_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the coordinate cards in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)

        model = self.model
        for (unused_id, coord) in sorted(model.coords.items()):
            if unused_id != 0:
                bdf_file = bdf_files[coord.ifile]
                try:
                    bdf_file.write(coord.write_card(size, is_double))
                except RuntimeError:
                    bdf_file.write(coord.write_card(16, is_double))

    def write_dmigs_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        model = self.model
        write_bdfs_dict(bdf_files, model.dmig, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dmi, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dmij, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dmiji, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dmik, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dmiax, size, is_double, is_long_ids)

    def write_dynamic_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the dynamic cards sorted by ID"""
        model = self.model
        is_dynamic = (model.dareas or model.dphases or model.nlparms or model.frequencies or
                      model.methods or model.cMethods or model.tsteps or model.tstepnls or
                      model.transfer_functions or model.delays or model.rotors or model.tics or
                      model.nlpcis)
        if is_dynamic:
            return
        write_bdfs_dict(bdf_files, model.methods, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.cMethods, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dareas, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dphases, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.nlparms, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.nlpcis, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.tsteps, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.tstepnls, size, is_double, is_long_ids)

        write_bdfs_dict_list(bdf_files, model.frequencies, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.delays, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.rotors, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.tics, size, is_double, is_long_ids)

        write_bdfs_dict_list(bdf_files, model.transfer_functions, size, is_double, is_long_ids)

    def write_loads_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.load_combinations or model.loads or model.tempds:
            write_bdfs_dict_list(bdf_files, model.load_combinations, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.loads, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.tempds, size, is_double, is_long_ids)
        self.write_dloads_file(bdf_files, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def write_dloads_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the dload cards sorted by ID"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.dloads or model.dload_entries:
            write_bdfs_dict_list(bdf_files, model.dloads, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.dload_entries, size, is_double, is_long_ids)

    def write_masses_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the mass cards sorted by ID"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.properties_mass:
            write_bdfs_dict(bdf_files, model.properties_mass, size, is_double, is_long_ids)
        if model.masses:
            write_bdfs_dict(bdf_files, model.masses, size, is_double, is_long_ids)

    def write_materials_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                             is_long_ids: Optional[bool]=None) -> None:
        """Writes the materials in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        is_materials = (model.materials or model.hyperelastic_materials or model.creep_materials or
                        model.MATS1 or model.MATS3 or model.MATS8 or model.MATT1 or
                        model.MATT2 or model.MATT3 or model.MATT4 or model.MATT5 or
                        model.MATT8 or model.MATT9 or model.nxstrats)
        if not is_materials:
            return
        write_bdfs_dict(bdf_files, model.materials, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.hyperelastic_materials, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.creep_materials, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATS1, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATS3, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATS8, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT1, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT2, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT3, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT4, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT5, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT8, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.MATT9, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.nxstrats, size, is_double, is_long_ids)

    def write_nodes_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the NODE-type cards"""
        model = self.model
        if model.spoints:
            write_xpoints_file(bdf_files, 'SPOINT', model.spoints)
        if model.epoints:
            write_xpoints_file(bdf_files, 'EPOINT', model.epoints)
        if model.points:
            write_bdfs_dict(bdf_files, model.points, size, is_double, is_long_ids)

        # if model._is_axis_symmetric:
        #     if model.axic:
        #         bdf_files[model.axic.ifile].write(model.axic.write_card(size, is_double))
        #     if model.axif:
        #         bdf_files[model.axif.ifile].write(model.axif.write_card(size, is_double))
        #     write_bdfs_dict(bdf_files, model.ringaxs, size, is_double, is_long_ids)
        #     write_bdfs_dict(bdf_files, model.ringfl, size, is_double, is_long_ids)
        #     write_bdfs_dict(bdf_files, model.gridb, size, is_double, is_long_ids)

        self.write_grids_file(bdf_files, size=size, is_double=is_double)
        if model.seqgp:
            bdf_files[model.seqgp.ifile].write(model.seqgp.write_card(size, is_double))

        #if 0:  # not finished
            #self.write_nodes_associated(bdf_file, size, is_double)

    def write_grids_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the GRID-type cards"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.nodes:
            if model.grdset:
                bdf_files[model.grdset.ifile].write(model.grdset.write_card(size))
            write_bdfs_dict(bdf_files, model.nodes, size, is_double, is_long_ids)

    def write_optimization_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        model = self.model
        is_optimization = (model.dconadds or model.dconstrs or model.desvars or model.ddvals or
                           model.dresps or
                           model.dvprels or model.dvmrels or model.dvcrels or model.doptprm or
                           model.dlinks or model.dequations or model.dtable is not None or
                           model.dvgrids or model.dscreen or model.topvar)
        if is_optimization:
            return
        write_bdfs_dict(bdf_files, model.dconadds, size, is_double, is_long_ids)
        write_bdfs_dict_list(bdf_files, model.dconadds, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.desvars, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.topvar, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.ddvals, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dlinks, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dresps, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dvcrels, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dvmrels, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.dvprels, size, is_double, is_long_ids)

        write_bdfs_dict_list(bdf_files, model.dvgrids, size, is_double, is_long_ids)

        for (unused_id, dscreen) in sorted(model.dscreen.items()):
            bdf_files[dscreen.ifile].write(str(dscreen))

        for (unused_id, equation) in sorted(model.dequations.items()):
            bdf_files[equation.ifile].write(str(equation))

        if model.dtable is not None:
            bdf_files[model.dtable.ifile].write(model.dtable.write_card(size, is_double))
        if model.doptprm is not None:
            bdf_files[model.doptprm.ifile].write(model.doptprm.write_card(size, is_double))
        if model.modtrak is not None:
            bdf_files[model.modtrak.ifile].write(model.modtrak.write_card(size, is_double))

    def write_params_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the PARAM cards
        """
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.params or model.dti:
            write_bdfs_dict(bdf_files, model.params, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.dti, size, is_double, is_long_ids)

    def write_properties_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                              is_long_ids: Optional[bool]=None) -> None:
        """Writes the properties in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        is_properties = model.properties or model.pelast or model.pdampt or model.pbusht
        if is_properties:
            write_bdfs_dict(bdf_files, model.properties, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.pelast, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.pdampt, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.pbusht, size, is_double, is_long_ids)

    def write_rejects_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the rejected (processed) cards and the rejected unprocessed
        cardlines
        """
        model = self.model
        if size == 8:
            print_func = print_card_8
        else:
            print_func = print_card_16

        if model.reject_cards:
            for reject_card in model.reject_cards:
                try:
                    bdf_files[0].write(print_func(reject_card))
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if model.reject_lines:
            #print(model.reject_lines)
            for reject_lines in model.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_files[0].write('%s\n' % reject2)
                elif isinstance(reject_lines, str):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_files[0].write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def write_rigid_elements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                  is_long_ids: Optional[bool]=None) -> None:
        """Writes the rigid elements in a sorted order"""
        size, is_long_ids = self.write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.rigid_elements:
            write_bdfs_dict(bdf_files, model.rigid_elements, size, is_double, is_long_ids)

        if model.plotels:
            write_bdfs_dict(bdf_files, model.plotels, size, is_double, is_long_ids)

    def write_sets_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the SETx cards sorted by ID"""
        model = self.model
        is_sets = (model.sets or model.asets or model.omits or model.bsets or model.csets or
                   model.qsets or model.usets)
        if is_sets:
            return
        write_bdfs_dict(bdf_files, model.sets, size, is_double, is_long_ids)
        write_bdfs_list(bdf_files, model.asets, size, is_double, is_long_ids)
        write_bdfs_list(bdf_files, model.omits, size, is_double, is_long_ids)
        write_bdfs_list(bdf_files, model.bsets, size, is_double, is_long_ids)
        write_bdfs_list(bdf_files, model.csets, size, is_double, is_long_ids)
        write_bdfs_list(bdf_files, model.qsets, size, is_double, is_long_ids)

        write_bdfs_dict_list(bdf_files, model.usets, size, is_double, is_long_ids)

    def write_superelements_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                 is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the Superelement cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        model = self.model
        is_sets = (model.se_sets or model.se_bsets or model.se_csets or
                   model.se_qsets or model.se_usets)
        if is_sets:
            write_bdfs_list(bdf_files, model.se_bsets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, model.se_csets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, model.se_qsets, size, is_double, is_long_ids)

            write_bdfs_dict(bdf_files, model.se_sets, size, is_double, is_long_ids)
            write_bdfs_dict_list(bdf_files, model.se_usets, size, is_double, is_long_ids)
            write_bdfs_list(bdf_files, model.se_suport, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.csuper, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.csupext, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.sebndry, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.sebulk, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.seconct, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.seelt, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.seexcld, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.seloc, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.seload, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.sempln, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.senqset, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.setree, size, is_double, is_long_ids)


    def write_tables_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the TABLEx cards sorted by ID"""
        model = self.model
        if model.tables or model.tables_d or model.tables_m or model.tables_sdamping:
            write_bdfs_dict(bdf_files, model.tables, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.tables_d, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.tables_m, size, is_double, is_long_ids)
            write_bdfs_dict(bdf_files, model.tables_sdamping, size, is_double, is_long_ids)

        if model.random_tables:
            write_bdfs_dict(bdf_files, self.random_tables, size, is_double, is_long_ids)

    def write_thermal_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal cards"""
        # PHBDY
        model = self.model
        is_thermal = (model.phbdys or model.convection_properties or model.bcs or
                      model.views or model.view3ds or model.radset or model.radcavs)
        if not is_thermal:
            return
        write_bdfs_dict(bdf_files, model.phbdys, size, is_double, is_long_ids)

        #for unused_key, prop in sorted(model.thermal_properties.items()):
        #    bdf_file.write(str(prop))
        write_bdfs_dict(bdf_files, model.convection_properties, size, is_double, is_long_ids)

        # BCs
        write_bdfs_dict_list(bdf_files, model.bcs, size, is_double, is_long_ids)

        write_bdfs_dict(bdf_files, model.views, size, is_double, is_long_ids)
        write_bdfs_dict(bdf_files, model.view3ds, size, is_double, is_long_ids)
        if model.radset:
            bdf_files[model.radset.ifile].write(model.radset.write_card(size, is_double))
        write_bdfs_dict(bdf_files, model.radcavs, size, is_double, is_long_ids)

    def write_thermal_materials_file(self, bdf_files: Any, size: int=8, is_double: bool=False,
                                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal materials in a sorted order"""
        model = self.model
        if model.thermal_materials:
            write_bdfs_dict(bdf_files, model.thermal_materials, size, is_double, is_long_ids)


def _fix_sizes(size: int,
               coords_size: Optional[int],
               nodes_size: Optional[int],
               elements_size: Optional[int],
               loads_size: Optional[int],
               table_size: Optional[int],
               flfact_size: Optional[int],
               ) -> tuple[int, int, int, int, int, int, int]:
    if coords_size is None:
        coords_size = size
    if nodes_size is None:
        nodes_size = size
    if elements_size is None:
        elements_size = size
    if loads_size is None:
        loads_size = size
    if table_size is None:
        table_size = size
    if flfact_size is None:
        flfact_size = size
    return size, coords_size, nodes_size, elements_size, loads_size, table_size, flfact_size


def get_optimization_include(model: BDF) -> tuple[list[int], list[int]]:
    """gets the properties and materials refereced by DVPRELx/DVMRELx"""
    property_types = {'PELAS', 'PDAMP', 'PGAP', 'PBUSH', 'PBUSH1D', 'PVISC', 'PWELD',
                      'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBMSECT',
                      'PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG', }
    material_types = {'MAT1', 'MAT8', 'MAT9'}

    pids_to_remove = set([])
    for dvprel_id, dvprel in model.dvprels.items():
        if dvprel.prop_type in property_types:
            pid = dvprel.pid
            pids_to_remove.add(pid)

    mids_to_remove = set([])
    for dvmrel_id, dvmrel in model.dvmrels.items():
        if dvmrel.mat_type in material_types:
            mid = dvmrel.mid
            mids_to_remove.add(mid)
    return list(pids_to_remove), list(mids_to_remove)


def delete_optimization_data(model: BDF) -> tuple[list[Any], list[Any], list[DESVAR]]:
    """removes optimization referenced data (that will be in the PCH file)"""
    pids_to_remove, mids_to_remove = get_optimization_include(model)
    properties_to_write, materials_to_write, desvars = _delete_optimization_data(
        model, pids_to_remove, mids_to_remove)
    return properties_to_write, materials_to_write, desvars


def _delete_optimization_data(model: BDF,
                              pids_to_remove: list[int],
                              mids_to_remove: list[int],
                              ) -> tuple[list[Any], list[Any], list[DESVAR]]:
    """heper method"""
    properties_to_write = []
    materials_to_write = []

    #--------------------------------------------------------------
    prop_types = set([])
    for pid in sorted(pids_to_remove):
        try:
            prop = model.properties[pid]
        except KeyError:
            continue
        prop_types.add(prop.type)

    pids_to_delete = []
    for pid, prop in model.properties.items():
        if prop.type in prop_types:
            properties_to_write.append(prop)
            pids_to_delete.append(pid)
    for pid in pids_to_delete:
        del model.properties[pid]

    # ----------------------------------------------------------------
    #mids_to_remove = list(mids_to_remove)
    mat_types = set([])
    for mid in sorted(mids_to_remove):
        try:
            mat = model.materials[mid]
        except KeyError:
            continue
        mat_types.add(mat.type)

    mids_to_delete = []
    for mid, mat in model.materials.items():
        if mat.type in mat_types:
            materials_to_write.append(mat)
            del model.materials[mid]
            mids_to_delete.append(mid)
    for mid in mids_to_delete:
        del model.mids_to_delete[mid]

    # ----------------------------------------------------------------
    desvars = []
    for desvar_id, desvar in sorted(model.desvars.items()):
        desvars.append(desvar)
    model.desvars = {}
    return properties_to_write, materials_to_write, desvars


def write_optimization_include(model: BDF, pch_include_filename: PathLike,
                               size: int=8) -> None:
    """writes an optimization deck"""
    properties_to_write, materials_to_write, desvars = delete_optimization_data(
        model)
    encoding = model.writer.get_encoding()
    with open(pch_include_filename, 'w', encoding=encoding) as pch_file:
        for prop in properties_to_write:
            pch_file.write(prop.write_card(size=size))
        for mat in materials_to_write:
            pch_file.write(mat.write_card(size=size))

        for desvar in desvars:
            pch_file.write(desvar.write_card(size=size))
    return
