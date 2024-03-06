import traceback
import numpy as np

from pyNastran.gui.menus.groups_modify.groups_modify import Group
from pyNastran.bdf.bdf_interface.model_group import ModelGroup
from pyNastran.bdf.utils import parse_patran_syntax


class GroupActions:
    def __init__(self, gui):
        self.gui = gui

    def create_groups_by_visible_result(self, nlimit: int=50) -> int:
        """
        Creates group by the active result

        This should really only be called for integer results < 50-ish.
        """
        log = self.gui.log
        try:
            #self.scalar_bar.title
            case_key = self.gui.case_keys[self.gui.icase] # int for object
            obj, (i, res_name) = self.gui.result_cases[case_key]
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
            self.gui.log_command('self.create_groups_by_visible_result()'
                                 f' # created {ngroups:d} groups for result_name={res_name!r}')
        except Exception as error:
            self.gui.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.gui.log_error(str(error))
            if self.gui.IS_GUI_TESTING:
                raise
        return ngroups

    def create_groups_by_property_id(self, nlimit: int=500):
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
        log = self.gui.log
        try:
            eids = self.gui.groups['main'].element_ids
            elements_pound = self.gui.groups['main'].elements_pound
        except Exception as error:
            log.error('Cannot create groups as there are no elements in the model')
            log.error(str(error))
            if self.gui.IS_GUI_TESTING:
                raise
            return 0

        result = self.gui.find_result_by_name(name, restype='fringe')
        ures = np.unique(result)
        ngroups = len(ures)
        if ngroups > nlimit:
            self.gui.log.error(f'not creating result; {ngroups:d} new groups would be created; '
                               f'increase nlimit={nlimit:d} if you really want to')
            return 0

        for uresi in ures:
            ids = np.where(uresi == result)[0]

            name = '%s %s' % (prefix, uresi)
            element_str = ''
            group = Group(
                name, element_str, elements_pound,
                editable=True)
            group.element_ids = eids[ids]
            self.gui.log_info('creating group=%r' % name)
            self.gui.groups[name] = group
        return ngroups

    def create_groups_by_model_group(self, model_groups: dict[int, ModelGroup],
                                     ) -> int:
        """
        Uses the model_groups object
        """
        #eids = self.find_result_by_name('ElementID', restype='fringe')
        #elements_pound = eids.max()
        log = self.gui.log
        try:
            eids = self.gui.groups['main'].element_ids
            elements_pound = self.gui.groups['main'].elements_pound
        except Exception as error:
            log.error('Cannot create groups as there are no elements in the model')
            log.error(str(error))
            if self.gui.IS_GUI_TESTING:
                raise
            return 0

        #result = self.gui.find_result_by_name(name, restype='fringe')
        ngroups = 0
        for name, group in model_groups.items():
            if group.elements is None:
                continue
            elements_str = group.get_patran_syntax(group.elements)
            elements = parse_patran_syntax(elements_str)

            ids = np.searchsorted(eids, elements)
            actual_eids = eids[ids]
            assert np.all(actual_eids == elements)

            element_str = ''
            group = Group(
                name, element_str, elements_pound,
                editable=True)
            group.element_ids = actual_eids
            self.gui.log_info(f'creating group={name!r}')
            self.gui.groups[name] = group
            ngroups += 1
        return ngroups
