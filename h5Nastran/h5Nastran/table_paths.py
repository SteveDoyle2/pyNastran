from __future__ import print_function, absolute_import


class TablePaths(object):
    bdf_lines_path = r'/H5NASTRAN/NASTRAN/INPUT'
    bdf_lines_table = 'BDF_LINES'

    unsupported_cards_path = r'/H5NASTRAN/NASTRAN/INPUT'
    unsupported_cards_table = 'UNSUPPORTED_CARDS'

    about_path = r'/H5NASTRAN/INFO'
    about_table = 'ABOUT'

    defaults_path = r'/H5NASTRAN/INFO'
    defaults_table = 'DEFAULTS'

    unsupported_result_tables_path = r'/H5NASTRAN/NASTRAN/RESULT'
    unsupported_result_tables_table = 'UNSUPPORTED_RESULT_TABLES'

    shell_element_info_path = r'/H5NASTRAN/NASTRAN/INPUT'
    shell_element_info_table = 'SHELL_ELEMENT_INFO'

    private_index_path = r'/H5NASTRAN/INDEX'

    def __getattr__(self, attr):
        dict = self.__class__.__dict__
        return dict[attr + '_path'] + '/' + dict[attr + '_table']
