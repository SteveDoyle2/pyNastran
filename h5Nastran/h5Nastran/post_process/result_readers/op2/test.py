from typing import Union, List, Tuple
from collections import OrderedDict


class Lines(object):
    def __init__(self, lines):
        self.lines = lines
        self.i = -1

    def _next_line(self):
        self.i += 1
        try:
            return self.lines[self.i]
        except IndexError:
            return None

    def next_line(self):
        while True:
            line = self._next_line()
            if line is None or not line.startswith('$'):
                return line


class DataType(object):
    def __init__(self, data):

        def _split(_data):
            _tmp = _data.split('{')

            if len(_tmp) == 1:
                tmp = _data.split(',')
                if tmp[-1] in ('', ' '):
                    del tmp[-1]
                tmp.append('{}')
            else:
                tmp = _tmp[0].split(',')
                if tmp[-1] in ('', ' '):
                    del tmp[-1]
                tmp.append('{' + _tmp[1])

            if len(tmp) == 5 and tmp[2] == '':
                del tmp[2]

            return tmp

        tmp = _split(data)

        print(tmp)

        if '{' in tmp[-1]:
            try:
                int(tmp[-2])
            except ValueError:
                tmp.insert(-1, '1')

            data = ','.join(tmp)
        else:
            try:
                int(tmp[-1])
            except ValueError:
                tmp.append('1')

            tmp.append('{}')
            data = ','.join(tmp)

        tmp = _split(data)

        if tmp[0] == 'UNDEF':
            data = ','.join(['UNDEF', 'I', tmp[1], tmp[2]])

        if '(' in tmp[0]:
            _tmp = tmp[0].split('(')
            name = _tmp[0]
            _size = _tmp[1][:-1]

            try:
                _size = int(_size)
                data = [name, tmp[1], str(_size * int(tmp[2])), tmp[3]]
                data = ','.join(data)
            except ValueError:
                pass

        tmp = _split(data)
        tmp[1] = tmp[1].strip()

        if 'CHAR' in tmp[1]:
            print(tmp)
            count = tmp[2]
            if count.strip() == '':
                count = 1
            else:
                count = int(tmp[2])
            _type = 'S' + str(int(tmp[1][4:]) * count)
            tmp = [tmp[0], _type, '1', tmp[3]]

        elif tmp[1] == 'I':
            tmp[1] = '<i4'

        elif tmp[1] == 'RS':
            tmp[1] = '<f8'

        data = ','.join(tmp)

        self.data = data


    def print_lines(self, offset=''):
        return [offset + str(self.data)]

    def __repr__(self):
        return '\n'.join(self.print_lines())


class Entry(object):
    def __init__(self, lines=None):
        self.data = []  # type: List[DataType]
        self.count = None
        self.end_with = ''

        if lines is not None:
            self.read(lines)

    def read(self, lines):
        while True:
            line = lines.next_line()
            if line is None or line == 'EOF':
                break

            if line.startswith('ENDENTRY'):
                tmp = line.split(',')
                if 'COUNT' in tmp[1]:
                    self.count = tmp[1].split('=')[1]
                    self.end_with = ''
                else:
                    self.count = None
                    if '(' in line:
                        tmp = line.split('(')
                        tmp = '(' + tmp[1]
                        self.end_with = tmp
                    else:
                        try:
                            self.end_with = tmp[2]
                        except IndexError:
                            self.end_with = ''
                break

            self.data.append(get_data(line, lines))

    def print_lines(self, offset=''):
        result = [offset + 'ENTRY,']

        for data in self.data:
            result.extend(data.print_lines(offset + '    '))

        if self.count is not None:
            result.append(offset + 'ENDENTRY,COUNT=%s' % self.count)
        else:
            result.append(offset + 'ENDENTRY,WITH,%s' % self.end_with)

        return result

    def __repr__(self):
        return '\n'.join(self.print_lines())


class Either(object):
    def __init__(self, check, lines=None):
        self.check = check
        self.data = OrderedDict()  # type: OrderedDict[str, Tuple[str, List[Union[DataType, Entry, Either]]]]

        if lines is not None:
            self.read(lines)

    def read(self, lines):
        try:
            tmp = self.check.split('=')[1].split(',')
        except IndexError:
            tmp = self.check.split(',')[1:]
            raise IndexError('This is an error in the ddl file.  There should be an extra entry before the EITHER statement; the EITHER statement should use that value.')

        _check = tmp[0]
        _descr = tmp[1]

        if '{' in _check:
            _descr = _check
            _check = ''

        _data = []

        while True:
            line = lines.next_line()

            if line is None or line.startswith('ENDEITH'):
                break

            if line.startswith('OR,'):
                tmp = line.split(',')
                check = tmp[1]
                try:
                    descr = tmp[2]
                except IndexError:
                    descr = ''

                if '{' in check:
                    descr = check
                    check = ''

                if len(_data) > 0:
                    self.data[_check] = (_descr, list(_data))
                    del _data[:]

                _check = check
                _descr = descr

                continue

            _data.append(get_data(line, lines))

        if len(_data) > 0:
            self.data[_check] = (_descr, list(_data))
            del _data[:]

    def print_lines(self, offset=''):
        result = [offset + self.check]

        keys = list(self.data.keys())

        def _print_data(_data, _offset=''):
            _result = []
            if isinstance(_data, list):
                for i in range(len(_data)):
                    _result.extend(_data[i].print_lines(_offset))
            else:
                _result.extend(_data.print_lines(_offset))

            return _result

        result.extend(_print_data(self.data[keys[0]][1], offset + '    '))

        for key in keys[1:]:
            _descr, _data = self.data[key]
            result.append(offset + 'OR,%s,%s' % (key, _descr))
            result.extend(_print_data(_data, offset + '    '))

        result.append(offset + 'ENDEITH,')

        return result

    def __repr__(self):
        return '\n'.join(self.print_lines())


def get_data(line, lines):
    if line.startswith('ENTRY'):
        return Entry(lines)
    elif line.startswith('EITHER'):
        return Either(line, lines)
    elif line.startswith('RECORD'):
        name = line.split('=')[1]
        return Record(name, lines)
    else:
        return DataType(line)


class Record(object):
    def __init__(self, name, lines=None):
        self.name = name
        self.data = []  # type: List[Union(DataType, Entry, Either)]

        if lines is not None:
            self.read(lines)

    def read(self, lines):
        while True:
            line = lines.next_line()
            if line is None or line.startswith('EOR,'):
                break

            self.data.append(get_data(line, lines))

    def print_lines(self, offset=''):
        result = [offset + 'RECORD=%s' % self.name]

        for data in self.data:
            result.extend(data.print_lines(offset + '    '))

        result.append(offset + 'EOR,')

        return result

    def __repr__(self):
        return '\n'.join(self.print_lines())


class Table(object):
    def __init__(self):
        self.records = []  # type: List[Record]

    def read(self, filename):
        with open(filename, 'r') as f:
            lines = f.read().split('\n')

        for i in range(len(lines)):
            line = lines[i]
            lines[i] = line.strip().upper()

        lines = Lines(lines)

        while True:
            line = lines.next_line()
            if line is None or line.startswith('EOF'):
                break
            if line.startswith('RECORD'):
                name = line.split('=')[1]
                self.records.append(Record(name, lines))

        return self

    def __repr__(self):
        result = []

        for record in self.records:
            result.extend(record.print_lines())

        return '\n'.join(result)


if __name__ == '__main__':
    filename = r'C:\MSC.Software\MD_Nastran\20101\md20101\nast\del\oes.ddl'

    table = Table()
    table.read(filename)

    print(table)

