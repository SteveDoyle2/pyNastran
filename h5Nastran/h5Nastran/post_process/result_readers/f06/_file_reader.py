from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import os


class FileReader(object):
    def __init__(self, filename):
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)

        self.f = open(self.filename, 'rb')

        some_data = self.f.read(1024)

        if b'\r\n' in some_data:
            self.separator = b'\r\n'
        else:
            self.separator = b'\n'

        self.f.seek(0)

        self.chunksize = int(self.filesize / 100)

        maxchunk = 1024 * 1000

        if self.chunksize > maxchunk:
            self.chunksize = maxchunk

        self._counter = 0

        self._data = []
        self._old_data = None
        self._new_data = None

        self._data_read = 0
        self._old_data_read = 0

        self._line_number = 0

    def __del__(self):
        self.close()

    def close(self):
        try:
            self.f.close()
        except AttributeError:
            pass

        self.f = None

    def next_line(self):
        try:
            tmp = self._data[self._counter]
            self._counter += 1
            self._line_number += 1
            return tmp
        except IndexError:
            if self._new_data is not None:
                self._data = self._new_data
                self._new_data = None
                self._counter = 0
                self._line_number = 0

                return self._data[self._counter]

            self._old_data = self._data

            self._old_data_read = self._data_read

            # print('reading data...')

            _data = self.f.read(self.chunksize)

            _data_len = _data.rfind(self.separator) + 1

            rewind = len(_data) - _data_len

            self.f.seek(self.f.tell() - rewind)

            _data = _data[:_data_len]

            self._data_read += len(_data)

            if len(_data) == 0:
                return None

            self._data = _data.split(self.separator)
            self._data.pop()
            self._counter = 0

            try:
                tmp = self._data[self._counter]
                self._counter += 1
                self._line_number += 1
                return tmp
            except IndexError:
                return None

    def previous_line(self):
        self._counter -= 1
        self._line_number -= 1

        try:
            return self._data[self._counter]
        except IndexError:
            if self._old_data is not None:
                self._new_data = self._data
                self._data = self._old_data
                self._old_data = None
                self._counter = len(self._data) - 1
                self._line_number = self._counter

                return self._data[self._counter]
            else:
                return None

    def line_number(self):
        return self._line_number
