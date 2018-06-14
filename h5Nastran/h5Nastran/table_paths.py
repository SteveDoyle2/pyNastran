from __future__ import print_function, absolute_import


class TablePaths(object):
    def __init__(self):
        # TODO: this should be saved to h5 file, and read back in when reading
        self.bdf_lines_path = r'/H5NASTRAN/NASTRAN/INPUT'
        self.bdf_lines_table = 'BDF_LINES'

        self.bdf_file_path = r'/H5NASTRAN/NASTRAN/INPUT'
        self.bdf_file_table = 'BDF_FILE'

        self.unsupported_cards_path = r'/H5NASTRAN/NASTRAN/INPUT'
        self.unsupported_cards_table = 'UNSUPPORTED_CARDS'

        self.about_path = r'/H5NASTRAN/INFO'
        self.about_table = 'ABOUT'

        self.versioning_path = r'/H5NASTRAN/INFO'
        self.versioning_table = 'VERSIONING'

        self.defaults_path = r'/H5NASTRAN/INFO'
        self.defaults_table = 'DEFAULTS'

        self.unsupported_result_tables_path = r'/H5NASTRAN/NASTRAN/RESULT'
        self.unsupported_result_tables_table = 'UNSUPPORTED_RESULT_TABLES'

        self.shell_element_info_path = r'/H5NASTRAN/NASTRAN/INPUT'
        self.shell_element_info_table = 'SHELL_ELEMENT_INFO'

        self.subcase_path = r'/H5NASTRAN/NASTRAN/RESULT'
        self.subcase_table = 'SUBCASES'

        self.private_index_path = r'/H5NASTRAN/INDEX'

    def __getattr__(self, attr):
        dict = self.__dict__
        return dict[attr + '_path'] + '/' + dict[attr + '_table']
