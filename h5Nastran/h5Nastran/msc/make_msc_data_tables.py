from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range


from collections import OrderedDict
import xml.etree.ElementTree as ET

import numpy as np



typedefs = {
    'integer': '<i8',
    'double': '<f8'
}


class Field(object):
    def __init__(self, name, data_type, shape, description):
        self.name = name
        self.data_type = data_type
        self.shape = shape
        self.description = description

    def __repr__(self):
        return self.name, self.data_type, self.shape, self.description

    def to_dtype(self):
        data_type = typedefs.get(self.data_type, self.data_type)
        try:
            data_type = data_type.to_dtype()
        except AttributeError:
            pass

        return self.name, data_type, self.shape


class Typedef(object):
    def __init__(self, name, fields):
        self.name = name
        self.fields = fields

    def __repr__(self):
        data = [self.name]
        data.append([field.__repr__() for field in self.fields])
        return str(data)

    def to_dtype(self):
        dtypes = []

        for field in self.fields:
            dtypes.append(field.to_dtype())

        return dtypes


def make_class(name, path, dtype, subtables=None, is_subtable=False, same_as=None):
    class_lines = [
        "class %s(object):" % name,
        "    name = '%s'" % name,
        "    path = '%s'" % path,
        "    dtype = %s" % str(dtype),
        "    is_subtable = %s" % str(is_subtable),
        "    same_as = '%s'" % str(same_as)
    ]

    if subtables is None:
        subtables = []

    class_lines.append('    subtables = %s' % str(subtables))

    return '\n'.join(class_lines)


class Group(object):
    def __init__(self, name, children):
        self.parent = None
        self.name = name
        self.children = children

        self.set_parent()

    def set_parent(self):
        for child in itervalues(self.children):
            child.parent = self
            child.set_parent()

    def path(self):
        try:
            parent_path = self.parent.path()
        except AttributeError:
            parent_path = []

        return parent_path + [self.name]

    def path_str(self):
        return '/'.join(self.path())

    def __getitem__(self, item):
        return self.children[item]

    def __repr__(self):
        return self.children.__repr__()

    def is_multitable(self):
        return 'IDENTITY' in self.children

    def make_class(self, prefix='', is_subtable=False):
        prefix = ''

        class_lines = []

        if self.is_multitable():
            prefix = prefix + self.name + '_'
            subtable_names = list(self.children.keys())
            subtable_names.remove('IDENTITY')
            children = list(itervalues(self.children))
            _subtables = [child.path_str() for child in children]
            subtables = []
            for sub in _subtables:
                if 'IDENTITY' not in sub:
                    subtables.append(sub)

            for subtable in subtable_names:
                _class_lines = self.children[subtable].make_class(prefix, is_subtable=True)
                if isinstance(_class_lines, str):
                    class_lines.append(_class_lines)
                else:
                    class_lines.extend(_class_lines)

            class_lines.append(self.children['IDENTITY'].make_class(prefix, subtables=subtables))

        else:
            for child in itervalues(self.children):
                _class_lines = child.make_class(prefix)
                if isinstance(_class_lines, str):
                    class_lines.append(_class_lines)
                else:
                    class_lines.extend(_class_lines)

        return class_lines


class Dataset(Typedef):
    def __init__(self, name, same_as, fields):
        super(Dataset, self).__init__(name, fields)
        self.same_as = same_as
        self.parent = None

    def set_parent(self):
        pass

    def path(self):
        try:
            parent_path = self.parent.path()
        except AttributeError:
            parent_path = []

        return parent_path + [self.name]

    def path_str(self):
        return '/'.join(self.path())

    def make_class(self, prefix='', subtables=None, is_subtable=False):
        prefix = ''
        if self.same_as is not None:
            same_as = self.same_as.split('/')
            path = self.path()
            print(path, same_as)
            try:
                index = path.index(same_as[0])
            except ValueError:
                index = len(path) - len(same_as)
            same_as_path = '/'.join(path[:index] + same_as)
            print(same_as_path)
        else:
            same_as_path = None

        return make_class(prefix + self.name, '/'.join(self.path()[:-1]), self.to_dtype(), subtables=subtables, is_subtable=is_subtable, same_as=same_as_path)


def get_field(data):
    items = OrderedDict(data.items())

    name = items['name']
    data_type = items['type']
    size = items.get('size', None)
    description = items.get('description', '')

    if size is None:
        shape = ()
    else:
        shape = (int(size),)

    if data_type == 'character':
        data_type = 'S%d' % int(size)
        shape = ()

    return Field(name, data_type, shape, description)


def get_typedef(data):
    items = OrderedDict(data.items())

    name = items['name']
    fields = []

    children = data.getchildren()

    for child in children:
        fields.append(get_field(child))

    return Typedef(name, fields)


def get_dataset(data):
    items = OrderedDict(data.items())

    name = items['name']
    same_as = items.get('sameAs', None)
    fields = []

    children = data.getchildren()

    for child in children:
        fields.append(get_field(child))

    return Dataset(name, same_as, fields)


def get_typedefs(data):
    children = data.getchildren()

    for child in children:
        typedef = get_typedef(child)
        typedefs[typedef.name] = typedef

    return typedefs


def get_group(parent):
    items = OrderedDict(parent.items())
    name = items.get('name', '')

    group_children = OrderedDict()

    children = parent.getchildren()

    for child in children:
        tag = child.tag
        if tag == 'group':
            child = get_group(child)
        elif tag == 'dataset':
            child = get_dataset(child)

        group_children[child.name] = child

    return Group(name, group_children)



tree = ET.parse('msc_datatype_2018.xml')
root = tree.getroot()


get_typedefs(root.getchildren()[0])


groups = get_group(root.getchildren()[1])

group = groups['NASTRAN']['INPUT']['PROPERTY']['PBARL']


class_lines = groups.make_class()

lines = [
    'from collections import OrderedDict',
    '',
    'data_tables = OrderedDict()',
    '',
    'def register_table(table):',
    '    table_id = table.path + "/" + table.name',
    '    data_tables[table_id] = table',
    '    return table',
    '',
    '',
    '@register_table'
]

lines.append('\n\n\n@register_table\n'.join(class_lines))


classes = '\n'.join(lines)

with open('msc_datatypes_classes.py', 'w') as f:
    f.write(classes)

