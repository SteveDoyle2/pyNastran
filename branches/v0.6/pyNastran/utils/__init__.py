## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
# -*- coding: utf-8 -*-
from types import MethodType
import os
from os.path import splitext, getsize
from os.path import join as pjoin
from numpy import ndarray
import io

def is_string(value):
    if isinstance(value, basestring):
        return True
    return False

def is_binary(filename):
    """
    Return true if the given filename is binary.

    :raises:  IOError if the file cannot be opened.
    :returns: True if filename is a binary file (contains null byte)
              and False otherwise.

    Based on the idea (..seealso:: http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text)
    that file is binary if it contains null.

    .. warning:: this may not work for unicode."""
    assert os.path.exists(filename), '%r does not exist' % filename
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b"\0" in chunk:  # found null byte
                return True
    return False


def obscure(num, debug=False):
    """
    Takes a large positive number and shrinks it down...similar to binary,
    but base 52.  A base 52 value takes up a fewer characters than a
    base 10 number, which helps to do Mat12345678 when there's an 8 character
    limit on variable names.

    >>> obscure(35214)
    'kbn'
    >>> de_obscure('kbn')
    35214

    :param num:   positive integer number
    :param debug: display additional information about conversion process
    :returns: shortened version of num"""
    vals = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    tmp = num
    pack = ['a'] if tmp == 0 else []
    while tmp > 0:
        if debug:
            print("num = %s\nfactor = %s" % (tmp, tmp % 52))
        pack.append(vals[tmp % 52])
        tmp //= 52

    if debug:
        print(pack,"\n\nsize chacnge %s > %s" % (len(str(num)), len(pack)))
    return "".join(pack)


def de_obscure(num, debug=False):
    """
    Unpacks an "obscured" positive number...similar to binary, but base 52.
    A base 52 value takes up a fewer characters than a base 10 number
    which helps to do Mat12345678 when there's an 8 character limit on variable names.

    >>> obscure(35214)
    'kbn'
    >>> de_obscure('kbn')
    35214

    :param debug: display additional information about conversion process
    :returns:     integer value of shortened version of a number
    .. seealso:: :func: `obscure`
    """
    dict_vals = dict(zip(list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'),xrange(52)))
    val = 0
    for i, letter in enumerate(list(num)):
        val += dict_vals[letter] * 52 ** i
        if debug:
            print("letter = ", letter, "\nfactor = ", dict_vals[letter] * 52 ** i)
    return val


def get_files_of_type(dirname, extension='.txt', maxSize=100.):
    """
    Gets the list of all the files with a given extension in the specified directory

    :param dirname:   the directory name
    :param extension: list of filetypes to get (default='.txt')
    :param maxSize:   size in MB for max file size
    :returns: list of all the files with a given extension in the specified directory
    """
    return [pjoin(dirname, f) for f in os.listdir(dirname) if extension in
            splitext(f)[1] and getsize(pjoin(dirname, f)) / (1048576.) <= maxSize]


def print_bad_path(path):
    """
    Prints information about the existence (access possibility) of the parts
    of the given path. Useful for debugging when the path to a given file
    is wrong.

    :param path: path to check
    :returns: string with informations whether access to parts of the path
              is possible
    """
    path = os.path.abspath(path)
    npath = os.path.dirname(path)
    res = [path]
    while path != npath:
        path, npath = npath, os.path.dirname(npath)
        res.append(path)
    msg = {True: "passed", False: "failed"}
    return "\n".join(["%s: %s" % (msg[os.path.exists(i)], i) for i in res])


def list_print(lst, float_fmt='%-4.2f'):
    """
    Prints a list, numpy array, or numpy matrix in an abbreviated format.
    Supported element types: None, string, numbers. Useful for debugging.

    :param lst: list, numpy array or numpy matrix
    :returns: the clean string representation of the object
    """
    def _print(val):
        if val is None or isinstance(val, basestring):
            return str(val)
        if isinstance(val, float):
            return float_fmt % val
        try:
            return '%g' % val
        except TypeError:
            print("parameter = %r" % val)
            raise

    try: # TODO: remove try block and fix bug in OP2 code or add a warning message
        if len(lst) == 0:
            return '[]'

        if isinstance(lst, ndarray) and lst.ndim == 2:
            r,c = lst.shape
            return ("["+",\n ".join(["["+",".join(['%-10g' % lst[i, j]
                                for j in range(c)])+"]" for i in range(r)])+"]")
        return "[" + ", ".join([_print(a) for a in lst]) + "]"
    except: # not a list
        return _print(lst)


def __object_attr(obj, mode, attr_type):
    """list object attributes of a given type"""
    test = {"public":  lambda k: not k.startswith('_'),
            "private": lambda k: (k.startswith('_') and not k.startswith('__')),
            "both": lambda k: not k.startswith('__'),
            "all":  lambda k: True
            }

    if not mode in test:
        print("Wrong mode! Accepted modes: public, private, both, all.")
        return None
    check = test[mode]

    out = []
    for k in dir(obj):
        if check(k) and attr_type(getattr(obj, k)):
            out.append(k)
    out.sort()
    return out
    #return sorted([k for k in dir(obj) if (check(k) and
    #                                           attr_type(getattr(obj, k)))])


def object_methods(obj, mode = "public"):
    """
    List the names of methods of a class as strings. Returns public methods
    as default.

    :param obj:  the object for checking
    :param mode: defines what kind of methods will be listed
       * "public" - names that do not begin with underscore
       * "private" - names that begin with single underscore
       * "both" - private and public
       * "all" - all methods that are defined for the object
    :returns:  sorted list of the names of methods of a given type
               or None if the mode is wrong
    """
    return __object_attr(obj, mode, lambda x: isinstance(x, MethodType))


def object_attributes(obj, mode = "public"):
    """
    List the names of attributes of a class as strings. Returns public attributes
    as default.

    :param obj:  the object for checking
    :param mode: defines what kind of attributes will be listed
       * "public" - names that do not begin with underscore
       * "private" - names that begin with single underscore
       * "both" - private and public
       * "all" - all attributes that are defined for the object
    :returns: sorted list of the names of attributes of a given type or None
              if the mode is wrong
    """
    return __object_attr(obj, mode, lambda x: not isinstance(x, MethodType))


def write_object_attributes(name, obj, nspaces=0, nbase=0, isClass=True, debug=False):
    """
    Writes a series of nested objects
    """
    spaces = (nbase+nspaces) * ' '
    msg = spaces
    xml = spaces
    if isClass:
        equals = '='
    else:
        equals = ':'

    if debug:
        print("attr=%s equals=|%s|" % (name, equals))
    # name
    if isinstance(obj, dict):
        if nspaces == 0:
            msg += '%s %s ' % (name, equals)
        else:
            if isinstance(name, tuple):
                msg += "%s %s " % (str(name), equals)
            else:
                msg += "'%s' %s " % (name, equals)
    elif isinstance(name, str):
        if isClass:
            key = '%s' % name
        else:
            key = "'%s'" % name
    elif isinstance(name, unicode):
        if isClass:
            key = u'%s' % name
        else:
            key = "u'%s'" % name
    elif isinstance(name, int) or isinstance(name, float) or isinstance(name, tuple) or name is None:
        key = "%s" % str(name)
    else:
        raise RuntimeError('key=%s is not a string.  Type=%s' % (name, type(name)))

    if debug:
        print("name=%s type=%s" % (name, type(obj)))

    # write the object
    if isinstance(obj, int) or isinstance(obj, float) or obj is None:
        xml += '<name=%s value=%s type=%s>' % (name, obj, type(obj))
        msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, isClass))
    elif is_string(obj):
        msg += "%s %s %s,\n" % (key, equals, write_value(obj, nspaces, nbase, isClass))

    elif isinstance(obj, dict):
        msg += write_dict(obj, nspaces, nbase, isClass) + ',\n'
    elif isinstance(obj, tuple) or isinstance(obj, list):
        msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, isClass))

    elif isinstance(obj, ndarray):
        starter = '%s%s %s' % (nspaces, key, equals)
        msg += '%s %s %s,\n' % (key, equals, write_array(obj, nspaces + 6 + len(starter)))
    else:  # generic class
        objectType = obj.__class__.__name__
        #raise RuntimeError('objectType=%s is not supported' % objectType)
        msg += "%s %s " % (key, equals)
        msg += write_class(name, obj, nspaces, nbase) + ',\n'  # comma for class
    if nspaces == 0:
        msg = msg[:-2]
    if debug:
        print("|%r|" % msg)
    return msg


def write_class(name, obj, nspaces=0, nbase=0):
    objectType = obj.__class__.__name__
    obj_attrs = object_attributes(obj, 'both')
    if not obj_attrs:
        return "%s()" % objectType

    spaces = ' ' * nspaces
    nspaces2 = nspaces + 4
    #spaces2 = nspaces2 * ' '
    msg = "%s(\n" % objectType
    for attr in obj_attrs[:-1]:
        value = getattr(obj, attr)
        #msg += '?'
        msg += write_object_attributes(attr, value, nspaces2, nbase, isClass=True)
    attr = obj_attrs[-1]
    value = getattr(obj, attr)
    msg += write_object_attributes(attr, value, nspaces2, nbase, isClass=True)
    msg += '%s)' % spaces

    #print "dir(obj) =", dir(obj)
    #print "obj_attrs =", obj_attrs
    return msg


def write_value(obj, nspaces, nbase=0, isClass=False):
    msg = ''
    if isinstance(obj, int) or isinstance(obj, float) or obj is None:
        msg += '%s' % (str(obj))
    elif isinstance(obj, str):
        msg += "'%s'" % obj
    elif isinstance(obj, unicode):
        msg += "u'%s'" % obj
    elif isinstance(obj, list):
        msg += write_list(obj, nspaces, nbase, isClass)
    elif isinstance(obj, tuple):
        msg += write_tuple(obj, nspaces, nbase, isClass)
    elif isinstance(obj, dict):
        msg += write_dict(obj, nspaces, nbase, isClass)
    else:
        objectType = type(obj)
        raise RuntimeError('objectType=%s is not supported; value=%s' % (objectType, obj))
    return msg


def write_dict(obj, nspaces, nbase, isClass):
    spaces = (nbase+nspaces) * ' '
    nspaces2 = nspaces + 4
    if len(obj) == 0:
        return '{}'

    msg = '{\n'
    for key, value in sorted(obj.iteritems()):
        #msg += '#'
        msg += write_object_attributes(key, value, nspaces2, nbase, isClass=False)
    msg += '%s}' % spaces
    return msg


def write_list(obj, nspaces, nbase, isClass):
    if len(obj) == 0:
        return '[]'

    spaces = ' ' * (nspaces + nbase)
    msg = '[\n%s    ' % spaces
    for value in obj[:-1]:
        msg += write_value(value, nspaces+4, nbase, isClass) + ', '
    msg += write_value(obj[-1], nspaces+4, nbase, isClass) + '\n%s]' % spaces
    return msg


def write_tuple(obj, nspaces, nbase, isClass):
    msg = '('
    for value in obj:
        msg += write_value(value, nspaces, nbase, isClass) + ', '
    msg += ')'
    return msg


def write_array(a, nspaces=0):
    return ' '*nspaces + '[???]'
    shape = a.shape
    dtype = a.dtype
    if len(shape) == 1:
        msg = 'array(['
        #print "a = ",a
        for ai in a[:-1]:
            #print "ai = ",ai
            if isinstance(ai, int) or isinstance(ai, float):
                msg += '%s, ' % ai
            elif isinstance(ai, str):
                msg += "'%s'," % ai
            else:
                objectType = type(ai)
                raise RuntimeError('objectType=%s is not supported; value=%s' % (objectType, ai))
                return "'array(.not supported type.)'"
            msg += '%s, ' % ai
        if len(a) > 0:
            if isinstance(a[-1], int) or isinstance(a[-1], float):
                msg += "%s], dtype='%s')" % (a[-1], dtype)
            elif isinstance(a[-1], str):
                msg += "'%s'], dtype='%s')" % (a[-1], dtype)
            else:
                objectType = type(ai)
                raise RuntimeError('objectType=%s is not supported; value=%s' % (objectType, ai))
                return "'array(.not supported type.)'"
        else:
            msg += '], dtype=%s)' % dtype
    elif len(shape) == 2:
        spaces = ' '*nspaces
        msg = 'array(['
        for i, ai in enumerate(a):
            if i > 0:
                msg += '%s[' % spaces
            for bi in ai[:-1]:
                msg += '%s, ' % bi
            msg += '%s' % ai[-1]

            if i+1 == len(a):
                msg += '], dtype=%s)' % dtype
            else:
                msg += '],\n'
    elif len(shape) == 3:
        return "'array(.not supported shape.)'"
    return msg


if __name__ == '__main__':
    from numpy import array, zeros
    class C(object):
        def __init__(self):
            pass
    class B(object):
        def __init__(self, x=None, e=None):
            self.x = 4
            self.e = C()

    class A(object):
        def __init__(self, a=None, b=None, c=None, d=None):
            self.a = a
            self.b = b
            self.c = c
            self.d = {'a' : 4,
                      'b' : [1,2,3],
                      'c' : {1:2},
                      'd' : B(),
                    (1,2) : 4,
            }

    z = zeros(2, dtype='float64')
    #print z
    #print z.dtype
    dictA = {
            'strString' : 'a string',
            'strFloat' : 1.0,
            'strInt': 2,
            'strTuple': (1,2),
            'strNone' : None,
            'strClass' : A('a', 'b', 'c'),
            'strList' : [1,2,3],
            'nullList' : [],
            'nullArray' : array([]),
            'stringArray' : array(['s']),
            'stringArray2' : array(['a', 'b']),
            'nullDict' : {},
            u'unicodStr' : u'',
            'ListOfLists' : [[[],[[],],2,{'a':3}]],
            1 : 1,
            None : 4,
            1.0 : 5,
            (1, 2) : 6,
            'strArray' : array([4,5,6]),
            'strArray2' : zeros((2,2)),
            'strArray3' : zeros((2,2,2)),
    }
    dictB = {
            'string2' : 'a string',
            'float2' : 1.0,
            'int2': 2,
            'dictA' : dictA,
    }


    dictC = {
        'dictA' : {
            None : 4,
            1 : 5,
            'strClass' : A(
                a = 'a',
                b = 'b',
                c = 'c',
            ),
            'strFloat' : 1.0,
            'strInt' : 2,
            'strNone' : None,
            'strString' : 'a string',
            'strTuple' : (1, 2),
            (1, 2) : 6,
        },
        'float2' : 1.0,
        'int2' : 2,
        'string2' : 'a string',
    }
    #assert sorted(dictB.items())==sorted(dictC.items())
    #print write_object_attributes('dictA', dictA, isClass=False)
    msg = write_object_attributes('dictB', dictB, nbase=0)
    print(msg)
    f = open('junk.py', 'wb')
    f.write(msg)
    f.close()

    import junk

    #dictB2 = eval(msg)
