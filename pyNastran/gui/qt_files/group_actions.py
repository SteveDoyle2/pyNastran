import traceback

import numpy as np

from pyNastran.gui.menus.groups_modify.groups_modify import Group


class GroupActions:
    def __init__(self, gui):
        self.gui = gui

    def create_groups_by_visible_result(self, nlimit=50):
        """
        Creates group by the active result

        This should really only be called for integer results < 50-ish.
        """
        try:
            #self.scalar_bar.title
            case_key = self.gui.case_keys[self.gui.icase] # int for object
            obj, (i, res_name) = self.gui.result_cases[case_key]
            default_title = obj.get_default_title(i, res_name)
            location = obj.get_location(i, res_name)
            if obj.data_format != '%i':
                self.gui.log.error('not creating result=%r; must be an integer result' % res_name)
                return 0
            if location != 'centroid':
                self.gui.log.error('not creating result=%r; must be a centroidal result' % res_name)
                return 0

            word = default_title
            prefix = default_title
            ngroups = self._create_groups_by_name(word, prefix, nlimit=nlimit)
            self.gui.log_command('create_groups_by_visible_result()'
                                 ' # created %i groups for result_name=%r' % (ngroups, res_name))
        except Exception as error:
            self.gui.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.gui.log_error(str(error))
            if self.gui.IS_GUI_TESTING:
                raise
        return ngroups

    def create_groups_by_property_id(self, nlimit=500):
        """
        Creates a group for each Property ID.

        As this is somewhat Nastran specific, create_groups_by_visible_result exists as well.
        """
        self._create_groups_by_name('PropertyID', 'property', nlimit=nlimit)
        self.gui.log_command('create_groups_by_property_id()')

    def _create_groups_by_name(self, name, prefix, nlimit=500):
        """
        Helper method for `create_groups_by_visible_result` and `create_groups_by_property_id`
        """
        #eids = self.find_result_by_name('ElementID')
        #elements_pound = eids.max()
        try:
            eids = self.gui.groups['main'].element_ids
            elements_pound = self.gui.groups['main'].elements_pound
        except Exception as error:
            self.gui.log.error('Cannot create groups as there are no elements in the model')
            self.gui.log.error(str(error))
            if self.gui.IS_GUI_TESTING:
                raise
            return 0

        result = self.gui.find_result_by_name(name)
        ures = np.unique(result)
        ngroups = len(ures)
        if ngroups > nlimit:
            self.gui.log.error('not creating result; %i new groups would be created; '
                               'increase nlimit=%i if you really want to' % (ngroups, nlimit))
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
