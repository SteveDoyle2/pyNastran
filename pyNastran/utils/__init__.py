# -*- coding: utf-8 -*-
from __future__ import print_function
from types import MethodType
import os
import io
from codecs import open as codec_open
from itertools import count

from six import PY2, string_types, iteritems, StringIO

import numpy as np

if PY2:
    integer_types = (int, long, np.int32, np.int64)
    integer_float_types = (int, long, np.int32, np.int64, float)
else:
    integer_types = (int, np.int32, np.int64)
    integer_float_types = (int, np.int32, np.int64, float)


def is_file_obj(filename):
    """does this object behave like a file object?"""
    #if not (hasattr(out_filename, 'read') and hasattr(out_filename, 'write')) or
    #        isinstance(out_filename, file) or isinstance(out_filename, StringIO):
    return ((hasattr(filename, 'read') and hasattr(filename, 'write'))
            or isinstance(filename, file)
            or isinstance(filename, StringIO))

def b(string):
    """reimplementation of six.b(...) to work in Python 2"""
    return string.encode('latin-1')

def merge_dicts(dict_list, strict=True):
    """merges two or more dictionaries"""
    assert isinstance(dict_list, list), type(dict_list)
    dict_out = {}
    for adict in dict_list:
        assert isinstance(adict, dict), adict
        for key, value in iteritems(adict):
            if key not in dict_out:
                dict_out[key] = value
            elif strict:
                raise RuntimeError('key=%r exists in multiple dictionaries' % key)
            else:
                print('key=%r is dropped?' % key)
    return dict_out


def is_binary_file(filename):
    """
    Return true if the given filename is binary.

    Parameters
    ----------
    filename : str
        the filename to test

    Returns
    -------
    binary_flag : bool
        True if filename is a binary file (contains null byte)
        and False otherwise.

    :raises:  IOError if the file cannot be opened.

    Based on the idea (.. seealso:: http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text)
    that file is binary if it contains null.

    .. warning:: this may not work for unicode."""
    assert isinstance(filename, string_types), '%r is not a valid filename' % filename
    assert os.path.exists(filename), '%r does not exist\n%s' % (filename, print_bad_path(filename))
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b'\0' in chunk:  # found null byte
                return True
    return False


def print_bad_path(path):
    """
    Prints information about the existence (access possibility) of the parts
    of the given path. Useful for debugging when the path to a given file
    is wrong.

    Parameters
    ----------
    path : str
        path to check

    Returns
    -------
    msg : str
        string with informations whether access to parts of the path
        is possible
    """
    #raw_path = path
    if len(path) > 255:
        path = os.path.abspath(_filename(path))
        npath = os.path.dirname(path)
        res = [path]
        while path != npath:
            path, npath = npath, os.path.dirname(npath)
            res.append(path)
        msg = {True: 'passed', False: 'failed'}
        return '\n'.join(['%s: %s' % (msg[os.path.exists(i)], i[4:]) for i in res])
    else:
        path = os.path.abspath(path)
        npath = os.path.dirname(path)
        res = [path]
        while path != npath:
            path, npath = npath, os.path.dirname(npath)
            res.append(path)
        msg = {True: 'passed', False: 'failed'}
        return '\n'.join(['%s: %s' % (msg[os.path.exists(i)], i) for i in res])

def _filename(filename):
    """
    Prepends some magic data to a filename in order to have long filenames.

    .. warning:: This might be Windows specific.
    """
    if len(filename) > 255:
        return '\\\\?\\' + filename
    return filename

def __object_attr(obj, mode, keys_to_skip, attr_type):
    """list object attributes of a given type"""
    #print('keys_to_skip=%s' % keys_to_skip)
    keys_to_skip = [] if keys_to_skip is None else keys_to_skip
    test = {
        'public':  lambda k: (not k.startswith('_') and k not in keys_to_skip),
        'private': lambda k: (k.startswith('_') and not k.startswith('__') and k not in keys_to_skip),
        'both': lambda k: (not k.startswith('__') and k not in keys_to_skip),
        'all':  lambda k: (k not in keys_to_skip),
    }

    if not mode in test:
        print('Wrong mode! Accepted modes: public, private, both, all.')
        return None
    check = test[mode]

    out = []
    for k in dir(obj):
        if k in keys_to_skip:
            continue
        if check(k) and attr_type(getattr(obj, k)):
            out.append(k)
    out.sort()
    return out
    #return sorted([k for k in dir(obj) if (check(k) and
    #                                           attr_type(getattr(obj, k)))])


def object_methods(obj, mode='public', keys_to_skip=None):
    """
    List the names of methods of a class as strings. Returns public methods
    as default.

    Parameters
    ----------
    obj : instance
        the object for checking
    mode : str
        defines what kind of methods will be listed
        * "public" - names that do not begin with underscore
        * "private" - names that begin with single underscore
        * "both" - private and public
        * "all" - all methods that are defined for the object
    keys_to_skip : List[str]; default=None -> []
        names to not consider to avoid deprecation warnings

    Returns
    -------
    method : List[str]
        sorted list of the names of methods of a given type
        or None if the mode is wrong
    """
    return __object_attr(obj, mode, keys_to_skip, lambda x: isinstance(x, MethodType))


def object_attributes(obj, mode='public', keys_to_skip=None):
    """
    List the names of attributes of a class as strings. Returns public
    attributes as default.

    Parameters
    ----------
    obj : instance
        the object for checking
    mode : str
        defines what kind of attributes will be listed
        * 'public' - names that do not begin with underscore
        * 'private' - names that begin with single underscore
        * 'both' - private and public
        * 'all' - all attributes that are defined for the object
    keys_to_skip : List[str]; default=None -> []
        names to not consider to avoid deprecation warnings

    Returns
    -------
    attribute_names : List[str]
        sorted list of the names of attributes of a given type or None
        if the mode is wrong
    """
    return __object_attr(obj, mode, keys_to_skip, lambda x: not isinstance(x, MethodType))


#def write_object_attributes(name, obj, nspaces=0, nbase=0, is_class=True, debug=False):
    #"""
    #Writes a series of nested objects
    #"""
    #spaces = (nbase + nspaces) * ' '
    #msg = spaces
    #xml = spaces
    #if is_class:
        #equals = '='
    #else:
        #equals = ':'

    #if debug:
        #print('attr=%s equals=%r' % (name, equals))
    ## name
    #if isinstance(obj, dict):
        #if nspaces == 0:
            #msg += '%s %s ' % (name, equals)
        #else:
            #if isinstance(name, tuple):
                #msg += '%s %s ' % (str(name), equals)
            #else:
                #msg += "'%s' %s " % (name, equals)
    #elif isinstance(name, string_types):
        #if is_class:
            #key = '%s' % name
        #else:
            #key = "'%s'" % name
    ## elif isinstance(name, unicode):
        ## if is_class:
            ## key = u'%s' % name
        ## else:
            ## key = "u'%s'" % name
    #elif isinstance(name, (int, float, tuple)) or name is None:
        #key = '%s' % str(name)
    #else:
        #raise RuntimeError('key=%s is not a string.  Type=%s' % (name, type(name)))

    #if debug:
        #print('name=%s type=%s' % (name, type(obj)))

    ## write the object
    #if isinstance(obj, (int, float)) or obj is None:
        #xml += '<name=%s value=%s type=%s>' % (name, obj, type(obj))
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))
    #elif is_string(obj):
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))

    #elif isinstance(obj, dict):
        #msg += write_dict(obj, nspaces, nbase, is_class) + ',\n'
    #elif isinstance(obj, (tuple, list)):
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))

    #elif isinstance(obj, np.ndarray):
        #starter = '%s%s %s' % (nspaces, key, equals)
        #msg += '%s %s %s,\n' % (key, equals, write_array(obj, nspaces + 6 + len(starter)))
    #else:  # generic class
        #objectType = obj.__class__.__name__
        ##raise RuntimeError('objectType=%s is not supported' % objectType)
        #msg += '%s %s ' % (key, equals)
        #msg += write_class(name, obj, nspaces, nbase) + ',\n'  # comma for class
    #if nspaces == 0:
        #msg = msg[:-2]
    #if debug:
        #print('%r' % msg)
    #return msg
