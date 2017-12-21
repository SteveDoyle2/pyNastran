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


class H5Nastran(object):

    version = '0.1.0'

    def __init__(self, h5filename, mode='r'):
        filters = tables.Filters(complib='zlib', complevel=5)
        self.h5f = tables.open_file(h5filename, mode=mode, filters=filters)

        self._card_tables = {}
        self._result_tables = {}

        self.input = Input(self)
        self.result = Result(self)

        self.bdf = None  # pyNastran bdf file

        self._tables = set()
        self._unsupported_tables = set()

        self._bdf = None
        self._punch = None
        self._f06 = None

        self._bdf_domain = 1

        if mode == 'w':
            self._write_info()

    def close(self):
        self.h5f.close()

    def load_bdf(self, filename=None):
        if self._bdf is not None:
            raise Exception('BDF already loaded!')

        if filename is None:
            self._load_bdf()
            return self.bdf

        self._bdf = filename

        self.bdf = BDF(debug=False)
        self.bdf.read_bdf(filename)

        bdf = self.bdf

        assert bdf is not None

        cards = get_bdf_cards(bdf)

        # TODO: I'm not quite sure how to interface with card_tables
        #       but using this should be quite a bit faster
        #cards = self.bdf._read_bdf_cards(filename)

        tables = set()
        unsupported = []

        card_names = sorted(cards.keys())

        for card_name in card_names:
            table = self._card_tables.get(card_name, None)

            if table is None:
                print(card_name)
                unsupported.append(card_name)
                continue

            try:
                table.write_data(cards[card_name], self._bdf_domain)
            except NotImplementedError:
                print(card_name)
                unsupported.append(card_name)

            tables.add(table)

        for table in tables:
            table.finalize()

        self._unsupported_cards(unsupported)

        self._bdf_domain += 1

        self._save_bdf()

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
        reader.register_callback(self._load_result_table)
        reader.read()

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
            self._result_tables[_result_type] = result_table

    def _load_bdf(self):
        from zlib import decompress

        bdf_lines = decompress(self.h5f.get_node('/PRIVATE/NASTRAN/INPUT/BDF_LINES').read()).decode()

        from six import StringIO

        class DummyIO(StringIO):
            # pyNastran expects StringIO to have a readlines method
            def readlines(self):
                return self.getvalue().split('\n')

        data = DummyIO()

        data.write(bdf_lines)

        bdf = BDF(debug=False)
        bdf.read_bdf(data)

        data.close()

        self.bdf = bdf

    def _load_result_table(self, table_data):
        print(table_data.header)

        results_type = table_data.header.results_type

        table = self._result_tables.get(results_type, None)

        if table is None:
            return self._unsupported_table(table_data)

        table.results_type = results_type

        table.write_data(table_data)

        self._tables.add(table)

    def _save_bdf(self):
        from six import StringIO

        out = StringIO()

        self.bdf.write_bdf(out, close=False)

        from zlib import compress

        self.h5f.create_array('/PRIVATE/NASTRAN/INPUT', 'BDF_LINES', obj=compress(out.getvalue().encode()), title='BDF LINES',
                              createparents=True)

    def _unsupported_cards(self, cards):
        cards = np.array(cards, dtype='S8')
        self.h5f.create_array('/PRIVATE/NASTRAN/INPUT', 'UNSUPPORTED_CARDS', obj=cards, title='UNSUPPORTED BDF CARDS',
                         createparents=True)

    def _unsupported_table(self, table_data):
        print('Unsupported table %s' % table_data.header.results_type)
        self._unsupported_tables.add(table_data.header.results_type)

    def _write_info(self):
        import pyNastran

        info = 'h5Nastran version %s\nPowered by pyNastran version %s' % (self.version, pyNastran.__version__)

        self.h5f.create_array('/PRIVATE/h5Nastran', 'h5Nastran', obj=info.encode(), title='h5Nastran Info',
                              createparents=True)

    def _write_unsupported_tables(self):
        headers = list(sorted(self._unsupported_tables))
        data = np.array(headers, dtype='S256')

        self.h5f.create_array('/PRIVATE/NASTRAN/RESULT', 'UNSUPPORTED_RESULT_TABLES', obj=data, title='UNSUPPORTED RESULT TABLES',
                         createparents=True)
