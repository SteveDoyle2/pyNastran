from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range


# python 2/3 compatibility for chr, is there a better way?
_test_chr = b'abcd'
try:
    chr(_test_chr[0])
except TypeError:
    def chr(x):
        return x


def convert_data(data):
    return data.strip()


class PunchHeaderData(object):
    def __init__(self):
        self.title = ''
        self.subtitle = ''
        self.label = ''
        self._subcase_id = ''
        self._random_id = ''
        self._results_type = ''
        self.real_output = False
        self.other = {}
        self.lineno = -1

    def clear(self):
        self.title = ''
        self.subtitle = ''
        self.label = ''
        self._subcase_id = ''
        self._results_type = ''
        self.real_output = False
        self.other.clear()

    def set_data(self, line):
        if line.startswith('$TITLE'):
            self.title = line.split('=')[1].strip()
        elif line.startswith('$SUBTITLE'):
            self.subtitle = line.split('=')[1].strip()
        elif line.startswith('$LABEL'):
            self.label = line.split('=')[1].strip()
        elif line.startswith('$SUBCASE ID'):
            self._subcase_id = line.split('=')[1].strip()
            self._random_id = None
        elif line.startswith('$RANDOM ID'):
            self._random_id = line.split('=')[1].strip()
            self._subcase_id = None
        elif line.startswith('$REAL OUTPUT'):
            self.real_output = True
        # elif line.startswith('$ELEMENT TYPE'):
        #     tmp = line[1:].split('=')
        #     right = tmp[1].strip().split()
        #     # TODO: needs to be done right
        #     if len(right) == 3:
        #         self.other['ELEMENT TYPE INFO'] = right[2]
        #     self.other['ELEMENT TYPE'] = ' '.join(right[:2])
        elif '=' in line:
            tmp = line[1:].split('=')
            left = tmp[0].strip()
            right = tmp[1].strip()
            self.other[left] = ' '.join(right.split())
        else:
            self._results_type = ' '.join(line[1:].strip().split())

    def __str__(self):
        return 'LINENO=%d, TITLE=%s, SUBTITLE=%s, LABEL=%s, SUBCASE ID=%s, Results Type=%s, REAL OUTPUT=%s, OTHER=%s' % (
            self.lineno, self.title, self.subtitle, self.label, self._subcase_id, self._results_type, str(self.real_output),
            str(self.other)
        )
    
    @property
    def options(self):
        results_type = set(self.results_type.split())
        no_options = set(self.results_type_no_options.split())
        return results_type - no_options

    @property
    def results_type(self):
        try:
            results_type = '%s %s' % (self._results_type, self.other['ELEMENT TYPE'])
        except KeyError:
            results_type = self._results_type

        if self.real_output:
            return '%s REAL' % results_type

        else:
            return '%s COMPLEX' % results_type

    @property
    def results_type_no_options(self):
        results_type = self.results_type
        try:
            elem_type = self.other['ELEMENT TYPE']
        except KeyError:
            return results_type

        tmp = elem_type.split()

        try:
            int(tmp[0])
            del tmp[0]
        except ValueError:
            pass

        del tmp[0]

        if len(tmp) == 0:
            return results_type

        results_type_ = results_type.split()
        for option in tmp:
            results_type_.remove(option)

        return ' '.join(results_type_)
    
    @property
    def results_type_basic(self):
        results_type = self.results_type_no_options
        if 'ELEMENT' in results_type:
            tmp = results_type.split()
            del tmp[-2]
            results_type = ' '.join(tmp)
            
        return results_type

    @property
    def subcase_id(self):
        try:
            return '%s Load Factor=%s' % (self._subcase_id, self.other['LOAD FACTOR'])
        except KeyError:
            return self._subcase_id

    def serialize(self):
        return self.title, self.subtitle, self.label, self._subcase_id, self._results_type, self.real_output, dict(self.other), self.lineno

    def load(self, data):
        self.title = data[0]
        self.subtitle = data[1]
        self.label = data[2]
        self._subcase_id = data[3]
        self._results_type = data[4]
        self.real_output = data[5]
        self.other = dict(data[6])
        self.lineno = data[7]


class PunchTableData(object):
    def __init__(self, table_data=None):
        self.header = PunchHeaderData()
        self.data = []

        if table_data is not None:
            self._load_data(table_data)

    def _load_data(self, table_data):
        self.header.clear()
        del self.data[:]

        for line in table_data:
            if chr(line[0]) == '$':
                self.header.set_data(line.decode())
                continue

            # _data = [convert_data(line[:10]), convert_data(line[10:18])]
            _data = [line[:10].strip(), line[10:18].strip()]

            i1 = 18
            i2 = 36

            line_len = len(line)

            while True:
                if i1 >= line_len:
                    break

                # _data.append(convert_data(line[i1:i2]))
                _data.append(line[i1:i2].strip())

                i1 = i2
                i2 += 18

            self.data.append(_data)

    def serialize(self):
        return list(self.data), self.header.serialize()

    def load(self, data):
        del self.data[:]
        self.data.extend(data[0])
        self.header.load(data[1])
