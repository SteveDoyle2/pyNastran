from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

from ._file_reader import FileReader
from ._table_data import PunchTableData


def _default_callback(table_data):
    print(table_data.header)


class OP2Reader(object):
    def __init__(self, filename):
        self.file = None  # fortran format reader
        self._done_reading = False
        self._callback = _default_callback

    def register_callback(self, callback):
        assert callable(callback)
        self._callback = callback

    def close(self):
        self.file.close()

    def read(self):
        self._done_reading = True  # remove this when done writing code

        while not self._done_reading:
            data = self._read_data_block()  # reads a table, matrix, whatever
            # data is a class of OP2Table, OP2Matrix, whatever... classes that don't exist yet
            # which are subclasses of ResultTableData

            # ResultTableData.data should be a numpy array, can be in memory if possible
            # but in general should be a numpy.memmap array
            # table/matrix formats are determined by ident tables
            # the data from the ident tables is used in a lookup table that will return the appropriate dtype
            # the lookup data can be defined in a comma delimited file

            # if data is result data, then it should be sorted by subcase if not already SORT1
            # this would be tricky though since each result table would be an element, for example,
            # and all element result tables haven't been read yet... they're passed over to the h5n db;
            # might be able to do this in pytables when done reading op2... just need to be able to
            # write the correct domains for each element's subcase

            # TODO: in future if geometry tables are to be read, then it doesn't make sense
            #       to subclass from ResultTableData; maybe should be moved and renamed to
            #       TableData or something
            self._callback(data)

    def _read_header(self):
        pass

    def _read_data_block(self):
        # if self.file.tell() == 0:
        #     return self._read_header()
        # otherwise read some data block
        pass
