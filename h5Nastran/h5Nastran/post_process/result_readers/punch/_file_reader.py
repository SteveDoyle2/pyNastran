from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import os


class FileReader(object):
    def __init__(self, filename):
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)

        self.f = open(self.filename, 'rb')

        tmp = self.f.read(100)

        if b'\r\n' in tmp:
            self.separator = b'\r\n'
            self.linesize = 82
        elif b'\n':
            self.separator = b'\n'
            self.linesize = 81
        else:
            raise Exception('%s is not a valid punch file!' % self.filename)

        assert self.filesize % self.linesize == 0, (self.linesize, self.filesize, tmp)

        self.f.seek(0)

        self.chunksize = 10000 * self.linesize

        self._counter = 0
        self._line_number = 0

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
            line = tmp[:72]
            self._line_number += 1
            return line
        except IndexError:
            if self._new_data is not None:
                self._data = self._new_data
                self._new_data = None
                self._counter = 0
                self._line_number += 1

                return self._data[self._counter][:72]

            self._old_data = self._data

            self._old_data_read = self._data_read

            _data = self.f.read(self.chunksize)

            # print(_data)

            self._data_read += len(_data)

            if len(_data) == 0:
                # print('None 1')
                return None

            if b'\n\r' in _data:
                self._data = _data.split(b'\n\r')
            else:
                self._data = _data.split(b'\n')

            # self._data = _data.split(self.separator)

            self._data.pop()
            self._counter = 0

            # print(self.separator)
            # print(self._data)

            try:
                tmp = self._data[self._counter]
                self._counter += 1
                line = tmp[:72]
                self._line_number += 1
                return line
            except IndexError:
                # print('None 2')
                return None

    def previous_line(self):
        self._counter -= 1
        self._line_number -= 1

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
        return self._line_number
