from __future__ import print_function, absolute_import

from ._base import H5NastranBase

from h5Nastran.exceptions import pyNastranReadBdfError, pyNastranWriteBdfError
from h5Nastran.pynastran_interface import get_bdf_cards

from pyNastran.bdf.bdf import BDF

import numpy as np


class H5NastranBDF(H5NastranBase):
    def __init__(self, *args, **kwargs):
        super(H5NastranBDF, self).__init__(*args, **kwargs)
        self.bdf = None  # pyNastran bdf file
        self._unsupported_bdf_cards = set()
        self._bdf = None
        
    # TODO: rename to register_input_table
    def register_card_table(self, card_table):
        assert card_table.card_id not in self._card_tables
        self._card_tables[card_table.card_id] = card_table
        
    def supported_from_bdf_cards(self):
        cards = []

        keys = sorted(self._card_tables.keys())

        for key in keys:
            table = self._card_tables[key]
            if table.from_bdf_implemented():
                cards.append(key)

        return cards

    def supported_to_bdf_cards(self):
        cards = []

        keys = sorted(self._card_tables.keys())

        for key in keys:
            table = self._card_tables[key]
            if table.to_bdf_implemented():
                cards.append(key)

        return cards

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
            raise pyNastranReadBdfError(
                "h5Nastran: pyNastran is unable to load the bdf '%s' for some reason." % filename)

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

            # print(card_name)

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

        self.nastran.input.element.write_shell_element_info(self.bdf, cards)

        self.nastran.input.update()

        return self.bdf
    
    def _load_bdf(self):
        from zlib import decompress

        def decompress(d):
            return d

        self._bdf = self.h5f.get_node(self.table_paths.bdf_file).read().decode()

        bdf_lines = decompress(self.h5f.get_node(self.table_paths.bdf_lines).read()).decode()

        from six import StringIO

        class DummyIO(StringIO):
            # pyNastran expects StringIO to have a readlines method
            def readlines(self):
                return self.getvalue().split('\n')

        data = DummyIO()

        data.write(bdf_lines)

        bdf = BDF(debug=False)

        _cards = self.supported_to_bdf_cards()

        for card in _cards:
            bdf.cards_to_read.remove(card)

        bdf.read_bdf(data, xref=False)
        data.close()

        for card in _cards:
            bdf.cards_to_read.add(card)
            
        del bdf.reject_lines[:]

        self.nastran.input.to_bdf(bdf)

        bdf.cross_reference()

        self.bdf = bdf

        self.bdf.bdf_filename = self._bdf
        
    def _save_bdf(self):
        from six import StringIO

        out = StringIO()

        try:
            self.bdf.write_bdf(out, close=False)
        except Exception:
            raise pyNastranWriteBdfError("h5Nastran: pyNastran is unable to write bdf '%s' for some reason." % self._bdf)

        from zlib import compress

        def compress(d):
            return d

        self.h5f.create_array(self.table_paths.bdf_file_path, self.table_paths.bdf_file_table,
                              obj=self._bdf.encode(), title='BDF FILE', createparents=True)

        self.h5f.create_array(self.table_paths.bdf_lines_path, self.table_paths.bdf_lines_table,
                              obj=compress(out.getvalue().encode()), title='BDF LINES', createparents=True)

        self.h5f.flush()

    def _unsupported_cards(self, cards):
        cards = np.array(cards, dtype='S8')
        self.h5f.create_array(self.table_paths.unsupported_cards_path, self.table_paths.unsupported_cards_table,
                              obj=cards, title='UNSUPPORTED BDF CARDS', createparents=True)

        self._unsupported_bdf_cards.clear()
        self._unsupported_bdf_cards.update(set(cards))
