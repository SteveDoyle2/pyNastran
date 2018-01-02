from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

from ._file_reader import FileReader
from ._table_data import PunchTableData


# python 2/3 compatibility for chr, is there a better way?
_test_chr = b'abcd'
try:
    chr(_test_chr[0])
except TypeError:
    def chr(x):
        return x


def _default_callback(table_data):
    print(table_data.header)


def _worker(data, lineno, d):
    try:
        tmp = PunchTableData(data)
        tmp.header.lineno = lineno
    except TypeError:
        d['error'] = True
        return

    d['error'] = False
    d['table_data'] = tmp.serialize()


class PunchReader(object):
    def __init__(self, filename):
        self.file = FileReader(filename)

        self._done_reading = False

        self._callback = _default_callback

    def register_callback(self, callback):
        assert callable(callback)

        self._callback = callback

    def close(self):
        self.file.close()

    def read_multiprocess(self):

        from multiprocessing import Process, Manager, cpu_count

        _count = 0

        _process = []
        _data = []
        _ln = []

        manager = Manager()

        max_count = cpu_count() - 1

        while not self._done_reading:
            table_data, line_number = self._read_table()

            # print(self.file.status())

            _ln.append(line_number)
            _data.append(manager.dict())

            d = _data[-1]

            p = Process(target=_worker, args=(table_data, line_number, d))
            p.start()

            _process.append(p)

            _count += 1

            breakout = False
            if _count == max_count:
                for i in range(len(_process)):
                    _process[i].join()
                    d = _data[i]

                    if d['error'] is True:
                        breakout = True
                        continue

                    td = PunchTableData()
                    td.load(d['table_data'])
                    self._callback(td)

                _process.clear()
                _data.clear()
                _count = 0

            if breakout:
                break

        for i in range(len(_process)):
            _process[i].join()
            d = _data[i]

            if d['error'] is True:
                continue

            td = PunchTableData()
            td.load(d['table_data'])
            self._callback(td)

    def read(self):
        while not self._done_reading:
            table_data, line_number = self._read_table()

            data = PunchTableData(table_data)
            data.header.lineno = line_number

            self._callback(data)

    def _read_table(self):

        table_data = []

        reading_data = False

        line_number = -1

        while True:
            next_line = self.file.next_line()

            if next_line is None:
                self._done_reading = True
                break

            if next_line.strip() == b'':
                break

            first_char = chr(next_line[0])

            if first_char == '$' and reading_data:
                self.file.previous_line()
                reading_data = False
                break

            if first_char != '$':
                reading_data = True

            if first_char == '-':
                try:
                    table_data[-1] += next_line[18:]
                except IndexError:
                    raise Exception('Error reading punch file %s!' % self.file.filename)

            else:
                if len(table_data) == 0:
                    line_number = self.file.line_number()

                table_data.append(next_line)

        if len(table_data) == 0:
            return None

        return table_data, line_number
