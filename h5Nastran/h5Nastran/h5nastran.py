from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

from pyNastran.bdf.bdf import BDF

import tables
import numpy as np

from .input import Input
from .result import Result

from .pynastran_interface import get_bdf_cards
from .punch import PunchReader
from .f06 import F06Reader
from .exceptions import pyNastranReadBdfError, pyNastranWriteBdfError
from .table_paths import TablePaths
from .data_helper import DataHelper


class H5Nastran(object):
    h5n_version_str = '0.1.0a0'
    h5n_version = (0, 1, 0)
    nastran_type = None
    nastran_version = (0, 0, 0)

    default_driver = None

    def __init__(self, h5filename, mode='r', nastran_type=None, nastran_version=None, in_memory=False):
        self.file_mode = mode

        if mode == 'r':
            assert nastran_type is None and nastran_version is None
        else:
            if nastran_type is not None:
                self.nastran_type = nastran_type

            if nastran_version is not None:
                assert isinstance(nastran_version, tuple) and len(nastran_version) == 3
                self.nastran_version = nastran_version

        if in_memory:
            driver = 'H5FD_CORE'
        else:
            driver = self.default_driver

        filters = tables.Filters(complib='zlib', complevel=5)
        self.h5f = tables.open_file(h5filename, mode=mode, filters=filters, driver=driver)

        self._card_tables = {}
        self._result_tables = {}

        self.input = Input(self)
        self.result = Result(self)

        self.bdf = None  # pyNastran bdf file

        self._tables = set()
        self._unsupported_tables = set()
        self._unsupported_bdf_cards = set()

        self._bdf = None
        self._punch = None
        self._f06 = None
        self._op2 = None

        self._element_results_tables = {}

        self.table_paths = TablePaths()

        self.defaults = DataHelper()

        if mode == 'w':
            self._write_info()
        else:
            self._update()

    def close(self):
        self.h5f.close()

    def visualize(self):
        from .gui.visualization import to_vtk

        if self.bdf is None:
            self.load_bdf()

        vtk_data = to_vtk(self.bdf)
        vtk_data.visualize()

    def load_bdf(self, filename=None):
        if self._bdf is not None:
            raise Exception('BDF already loaded!')

        if filename is None:
            self._load_bdf()
            return self.bdf

        self._bdf = filename

        self.bdf = BDF(debug=False)
        try:
            self.bdf.read_bdf(filename)  # allow xref, could catch bdf errors
        except Exception:
            raise pyNastranReadBdfError("h5Nastran: pyNastran is unable to load the bdf '%s' for some reason." % filename)

        bdf = self.bdf

        assert bdf is not None

        cards = get_bdf_cards(bdf)

        # TODO: I'm not quite sure how to interface with card_tables
        #       but using this should be quite a bit faster
        # cards = self.bdf._read_bdf_cards(filename)
        # TODO: figure out how to use bdf._read_bdf_cards

        tables = set()
        unsupported = []

        card_names = sorted(cards.keys())

        for card_name in card_names:
            table = self._card_tables.get(card_name, None)

            print(card_name)

            if table is None:
                print(card_name, 'not supported')
                unsupported.append(card_name)
                continue

            try:
                table.write_data(cards[card_name])
            except NotImplementedError:
                print(card_name, 'not supported')
                unsupported.append(card_name)

            tables.add(table)

        for table in tables:
            table.finalize()

        self._unsupported_cards(unsupported)

        self._save_bdf()
        
        self.input.element.write_shell_element_info(self.bdf, cards)

        self.input.update()

        return self.bdf

    def load_f06(self, f06file):
        if self._bdf is None:
            raise Exception('BDF must be loaded first!')

        if self._punch is not None:
            raise Exception('Punch file has already been loaded.  Cannot load f06 file after punch.')

        self._f06 = f06file

        reader = F06Reader(f06file)
        reader.register_callback(self._load_result_table)
        reader.read()

        for table in self._tables:
            table.finalize()

        self._tables.clear()

    def load_punch(self, filename):
        if self._bdf is None:
            raise Exception('BDF must be loaded first!')

        if self._f06 is not None:
            raise Exception('F06 has already been loaded.  Cannot load punch file after f06.')

        self._punch = filename

        reader = PunchReader(filename)
        reader.register_callback(self._load_punch_table)
        reader.read()

        self.h5f.flush()

        for table in self._tables:
            table.finalize()

        self._tables.clear()
        self._write_unsupported_tables()

    def load_op2(self, filename):
        if self._bdf is None:
            raise Exception('BDF must be loaded first!')

        if self._punch is not None:
            raise Exception('Punch file has already been loaded.  Cannot load op2 file after punch.')

        if self._f06 is not None:
            raise Exception('F06 file has already been loaded.  Cannot load op2 file after f06.')

        self._op2 = filename

        # TODO: need OP2Reader
        reader = OP2Reader(filename)
        reader.register_callback(self._load_op2_table)
        reader.read()

        self.h5f.flush()

        for table in self._tables:
            table.finalize()

        self._tables.clear()
        self._write_unsupported_tables()

    def path(self):
        return ['', 'NASTRAN']

    def register_card_table(self, card_table):
        assert card_table.card_id not in self._card_tables
        self._card_tables[card_table.card_id] = card_table

    def register_result_table(self, result_table):
        result_type = result_table.result_type

        if isinstance(result_type, str):
            result_type = [result_type]

        for _result_type in result_type:
            assert _result_type not in self._result_tables

            if 'ELEMENT' in _result_type:
                tmp = _result_type.split()
                try:
                    int(tmp[-3])
                    del tmp[-2]
                except ValueError:
                    pass
                _result_type = ' '.join(tmp)

            self._result_tables[_result_type] = result_table

    def supported_from_bdf_cards(self):
        cards = []

        from .input.input_table import InputTable

        keys = sorted(self._card_tables.keys())

        for key in keys:
            table = self._card_tables[key]
            if table.from_bdf_implemented():
                cards.append(key)

        return cards

    def supported_to_bdf_cards(self):
        cards = []

        from .input.input_table import InputTable

        keys = sorted(self._card_tables.keys())

        for key in keys:
            table = self._card_tables[key]
            if table.to_bdf_implemented():
                cards.append(key)

        return cards

    def _load_bdf(self):
        from zlib import decompress

        bdf_lines = decompress(self.h5f.get_node(self.table_paths.bdf_lines).read()).decode()

        from six import StringIO

        class DummyIO(StringIO):
            # pyNastran expects StringIO to have a readlines method
            def readlines(self):
                return self.getvalue().split('\n')

        data = DummyIO()

        data.write(bdf_lines)

        bdf = BDF(debug=False)

        for card in self.supported_to_bdf_cards():
            bdf.cards_to_read.remove(card)

        bdf.read_bdf(data, xref=False)
        data.close()

        self.input.to_bdf(bdf)

        bdf.cross_reference()

        self.bdf = bdf

    def _load_result_table(self, table_data):
        print(table_data.header)

        results_type = table_data.header.results_type
        results_type_ = table_data.header.results_type_no_options

        print(results_type)
        print(results_type_)

        table = self._result_tables.get(results_type_, None)

        if table is None:
            table = self._find_element_result_table(results_type_)

        if table is None:
            return self._unsupported_table(table_data)

        table.results_type = results_type

        table.write_data(table_data)

        self._tables.add(table)

    def _load_op2_table(self, table_data):
        result_type = table_data.result_type

        table = self._result_tables.get(result_type, None)

        if table is None:
            return self._unsupported_table(table_data)

        table.write_op2_data(table_data)

        self._tables.add(table)

    def _load_punch_table(self, table_data):
        print(table_data.header)

        results_type = table_data.header.results_type_basic

        table = self._result_tables.get(results_type, None)

        if table is None:
            return self._unsupported_table(table_data)

        table.write_punch_data(table_data)

        self._tables.add(table)

    # TODO: remove element type from results tables, since it's hard to know what they should be
    # item code is good enough
    def _find_element_result_table(self, results_type):
        if not results_type.startswith('ELEMENT'):
            return None

        tmp = results_type.split(' ')
        del tmp[2]
        results_type = ' '.join(tmp)

        if len(self._element_results_tables) == 0:
            keys = self._result_tables.keys()
            for key in keys:
                if not key.startswith('ELEMENT'):
                    continue
                tmp = key.split(' ')
                del tmp[2]
                tmp = ' '.join(tmp)
                self._element_results_tables[tmp] = key

        _results_type = self._element_results_tables.get(results_type, None)

        return self._result_tables.get(_results_type, None)

    def _save_bdf(self):
        from six import StringIO

        out = StringIO()

        try:
            self.bdf.write_bdf(out, close=False)
        except Exception:
            raise pyNastranWriteBdfError("h5Nastran: pyNastran is unable to write bdf '%s' for some reason." % self._bdf)

        from zlib import compress

        self.h5f.create_array(self.table_paths.bdf_lines_path, self.table_paths.bdf_lines_table,
                              obj=compress(out.getvalue().encode()), title='BDF LINES', createparents=True)

        self.h5f.create_table(self.table_paths.defaults_path, self.table_paths.defaults_table, self.defaults.Format,
                              'DATA DEFAULTS', expectedrows=1, createparents=True)

        table = self.h5f.get_node(self.table_paths.defaults)
        table.append(self.defaults.save())

        self.h5f.flush()

    def _unsupported_cards(self, cards):
        cards = np.array(cards, dtype='S8')
        self.h5f.create_array(self.table_paths.unsupported_cards_path, self.table_paths.unsupported_cards_table,
                              obj=cards, title='UNSUPPORTED BDF CARDS', createparents=True)
        
        self._unsupported_bdf_cards.clear()
        self._unsupported_bdf_cards.update(set(cards))

    def _unsupported_table(self, table_data):
        print('Unsupported table %s' % table_data.header.results_type)
        self._unsupported_tables.add(table_data.header.results_type)

    def _write_info(self):
        import pyNastran

        info = 'h5Nastran version %s\nPowered by pyNastran version %s' % (self.h5n_version_str, pyNastran.__version__)

        self.h5f.create_array(self.table_paths.about_path, self.table_paths.about_table, obj=info.encode(),
                              title='h5Nastran Info', createparents=True)

        versioning_dtype = np.dtype([
            ('H5NASTRAN_VERSION_STR', 'S8'),
            ('H5NASTRAN_VERSION', '<i8', (3,)),
            ('NASTRAN_TYPE', 'S8'),
            ('NASTRAN_VERSION', '<i8', (3,))
        ])

        format = tables.descr_from_dtype(versioning_dtype)[0]

        self.h5f.create_table(self.table_paths.versioning_path, self.table_paths.versioning_table, format,
                              'VERSIONING', expectedrows=1, createparents=True)

        table = self.h5f.get_node(self.table_paths.versioning)

        data = np.zeros(1, dtype=versioning_dtype)

        data['H5NASTRAN_VERSION_STR'][0] = self.h5n_version_str
        data['H5NASTRAN_VERSION'][0] = self.h5n_version

        nastran_type = self.nastran_type

        if nastran_type is None:
            nastran_type = ''

        data['NASTRAN_TYPE'][0] = nastran_type
        data['NASTRAN_VERSION'][0] = self.nastran_version

        table.append(data)

        self.h5f.flush()

    def _write_unsupported_tables(self):
        headers = list(sorted(self._unsupported_tables))
        data = np.array(headers, dtype='S256')

        self.h5f.create_array(self.table_paths.unsupported_result_tables_path,
                              self.table_paths.unsupported_result_tables_table, obj=data,
                              title='UNSUPPORTED RESULT TABLES', createparents=True)

    def _update(self):
        self.input.update()
        defaults = self.h5f.get_node(self.table_paths.defaults).read()
        self.defaults.load(defaults)
