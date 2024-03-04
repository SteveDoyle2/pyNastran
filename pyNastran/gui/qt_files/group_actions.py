from __future__ import annotations
import traceback
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.gui.menus.groups_modify.groups_modify import Group
from pyNastran.bdf.bdf_interface.model_group import ModelGroup
from pyNastran.bdf.utils import parse_patran_syntax
if TYPE_CHECKING:
    from pyNastran.gui.main_window import MainWindow

class GroupActions:
    def __init__(self, gui: MainWindow):
        self.gui: MainWindow = gui

    def create_groups_by_visible_result(self, nlimit: int=50) -> int:
        """
        Creates group by the active result

        This should really only be called for integer results < 50-ish.
        """
        gui: MainWindow = self.gui
        log = gui.log
        try:
            #self.scalar_bar.title
            case_key = gui.case_keys[gui.icase] # int for object
            obj, (i, res_name) = gui.result_cases[case_key]
            default_title = obj.get_default_legend_title(i, res_name)
            location = obj.get_location(i, res_name)
            if obj.data_format != '%i':
                log.error(f'not creating result={res_name!r}; must be an integer result')
                return 0
            if location != 'centroid':
                log.error(f'not creating result={res_name!r}; must be a centroidal result')
                return 0

            word = default_title
            prefix = default_title
            ngroups = self._create_groups_by_name(word, prefix, nlimit=nlimit)
            gui.log_command('self.create_groups_by_visible_result()'
                                 f' # created {ngroups:d} groups for result_name={res_name!r}')
        except Exception as error:
            gui.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            gui.log_error(str(error))
            if gui.IS_GUI_TESTING:
                raise
        return ngroups

    def create_groups_by_property_id(self, nlimit: int=500) -> None:
        """
        Creates a group for each Property ID.

        As this is somewhat Nastran specific,
        create_groups_by_visible_result exists as well.
        """
        self._create_groups_by_name('PropertyID', 'property', nlimit=nlimit)
        self.gui.log_command('self.create_groups_by_property_id()')

    def _create_groups_by_name(self, name: str, prefix: str,
                               nlimit: int=500) -> int:
        """
        Helper method for `create_groups_by_visible_result` and
        `create_groups_by_property_id`
        """
        #eids = self.find_result_by_name('ElementID', restype='fringe')
        #elements_pound = eids.max()
        gui: MainWindow = self.gui
        log = gui.log
        try:
            group = gui.groups['main']
            model = gui.models['main']
        except Exception as error:  # pragma: no cover
            log.error('Cannot create groups as there is no main geometry')
            log.error(str(error))
            if gui.IS_GUI_TESTING:
                raise
            return 0

        try:
            eids = group.element_ids
            nids = group.node_ids
            elements_pound = group.elements_pound
            nodes_pound = group.nodes_pound
        except Exception as error:  # pragma: no cover
            log.error('Cannot create groups as there are no elements in the model')
            log.error(str(error))
            if gui.IS_GUI_TESTING:
                raise
            return 0

        result = gui.find_result_by_name(name, restype='fringe')
        ures = np.unique(result)
        ngroups = len(ures)
        if ngroups > nlimit:
            gui.log.error(f'not creating result; {ngroups:d} new groups would be created; '
                          f'increase nlimit={nlimit:d} if you really want to')
            return 0

        for uresi in ures:
            ids = np.where(uresi == result)[0]
            eidsi = eids[ids]
            nidsi = model.get_node_ids_with_elements(eidsi, msg='', return_set=False)

            name = f'{prefix} {uresi}'
            element_str = ''
            node_str = ''
            group = Group(
                name,
                element_str, elements_pound,
                node_str, nodes_pound,
                editable=True)
            group.element_ids = eidsi
            group.node_ids = nidsi
            gui.log_info(f'creating group={name!r}')
            gui.groups[name] = group
        return ngroups

    def create_groups_by_model_group(self, model_groups: dict[int, ModelGroup],
                                     ) -> int:
        """
        Uses the model_groups object
        """
        #eids = self.find_result_by_name('ElementID', restype='fringe')
        #elements_pound = eids.max()
        gui = self.gui
        log = gui.log
        group = gui.groups['main']
        try:
            eids = group.element_ids
            nids = group.node_ids
            elements_pound = group.elements_pound
            nodes_pound = group.nodes_pound
        except Exception as error:
            log.error('Cannot create groups as there are no elements in the model')
            log.error(str(error))
            if gui.IS_GUI_TESTING:
                raise
            return 0

        #result = gui.find_result_by_name(name, restype='fringe')
        ngroups = 0
        for name, group in model_groups.items():
            if group.elements is None:
                continue
            elements_str = group.get_patran_syntax(group.elements)
            elements = parse_patran_syntax(elements_str)

            nodes_str = ''
            nodes = np.array([], dtype='int32')
            if group.nodes is not None:
                nodes_str = group.get_patran_syntax(group.nodes)
                nodes = parse_patran_syntax(nodes_str)

            ids = np.searchsorted(eids, elements)
            actual_eids = eids[ids]
            assert np.all(actual_eids == elements)

            _node_str = ''
            _element_str = ''
            group = Group(
                name,
                _element_str, elements_pound,
                _node_str, nodes_pound,
                editable=True)
            group.element_ids = actual_eids
            group.node_ids = nodes
            gui.log_info(f'creating group={name!r}')
            gui.groups[name] = group
            ngroups += 1
        return ngroups
