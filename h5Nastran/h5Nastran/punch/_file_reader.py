from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import os


class FileReader(object):
    def __init__(self, filename):
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)

        if self.filesize % 82 == 0:
            self.separator = b'\r\n'
            self.linesize = 82
        elif self.filesize % 81 == 0:
            self.separator = b'\n'
            self.linesize = 81
        else:
            raise Exception('%s is not a valid punch file!' % self.filename)

        self.f = open(self.filename, 'rb')

        self.chunksize = 10000 * self.linesize

        self._counter = 0

        self._data = []
        self._old_data = None
        self._new_data = None

        self._data_read = 0
        self._old_data_read = 0

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
            return tmp[:72]
        except IndexError:
            if self._new_data is not None:
                self._data = self._new_data
                self._new_data = None
                self._counter = 0

                return self._data[self._counter][:72]

            self._old_data = self._data

            self._old_data_read = self._data_read

            _data = self.f.read(self.chunksize)

            self._data_read += len(_data)

            if len(_data) == 0:
                return None

            self._data = _data.split(self.separator)
            self._data.pop()
            self._counter = 0

            try:
                tmp = self._data[self._counter]
                self._counter += 1
                return tmp[:72]
            except IndexError:
                return None

    def previous_line(self):
        self._counter -= 1

        try:
            return self._data[self._counter][:72]
        except IndexError:
            if self._old_data is not None:
                self._new_data = self._data
                self._data = self._old_data
                self._old_data = None
                self._counter = len(self._data) - 1

                return self._data[self._counter][:72]
            else:
                return None

    def line_number(self):
        return self._data_read / self.linesize + self._counter
