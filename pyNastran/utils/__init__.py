# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import
from types import MethodType, FunctionType
import os
import io
import sys
from codecs import open
from itertools import count

from typing import List, Union, Optional
from six import PY2, string_types, StringIO

if PY2:
    FileNotFoundError = IOError
    unicode_type = unicode
else:
    unicode_type = str


#if PY2:
    #def ChainMap(*keys):
        #"""Python 2.7 hack to implement ChainMap"""
        #keys2 = []
        #for key in keys:
            #keys2 += list(key)
        #return keys2
#else:
    #from collections import ChainMap

def ipython_info():
    # type: () -> Optional[str]
    """determines if iPython/Jupyter notebook is running"""
    try:
        return get_ipython()
    except NameError:
        return None

def is_file_obj(filename):
    """does this object behave like a file object?"""
    if PY2:
        return (
            (hasattr(filename, 'read') and hasattr(filename, 'write'))
            or isinstance(filename, file)
            or isinstance(filename, StringIO)
        )
    return (
        (hasattr(filename, 'read') and hasattr(filename, 'write'))
        or isinstance(filename, io.IOBase)
        or isinstance(filename, StringIO)
    )

def b(string):
    # type: (str) -> bytes
    """reimplementation of six.b(...) to work in Python 2"""
    return string.encode('latin-1')

#def merge_dicts(dict_list, strict=True):
    #"""merges two or more dictionaries"""
    #assert isinstance(dict_list, list), type(dict_list)
    #dict_out = {}
    #for adict in dict_list:
        #assert isinstance(adict, dict), adict
        #for key, value in adict.items():
            #if key not in dict_out:
                #dict_out[key] = value
            #elif strict:
                #raise RuntimeError('key=%r exists in multiple dictionaries' % key)
            #else:
                #print('key=%r is dropped?' % key)
    #return dict_out


def is_binary_file(filename):
    # type: (str) -> bool
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
    check_path(filename)
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b'\0' in chunk:  # found null byte
                return True
    return False


def check_path(filename, name='file'):
    # type: (str, str) -> None
    try:
        exists = os.path.exists(filename)
    except TypeError:
        msg = 'cannot find %s=%r\n' % (name, filename)
        raise TypeError(msg)
    if not exists:
        msg = 'cannot find %s=%r\n%s' % (name, filename, print_bad_path(filename))
        raise FileNotFoundError(msg)

def print_bad_path(path):
    # type: (str) -> str
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
    # type: (str) -> str
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
        try:
            if check(k) and attr_type(getattr(obj, k)):
                out.append(k)
        except:
            pass
    out.sort()
    return out
    #return sorted([k for k in dir(obj) if (check(k) and
    #                                           attr_type(getattr(obj, k)))])


def object_methods(obj, mode='public', keys_to_skip=None):
    # type: (object, str, Optional[List[str]]) -> List[str]
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
    # type: (object, str, Optional[List[str]]) -> List[str]
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
    #if hasattr(obj, '__properties__'):
        #keys_to_skip += obj.__properties__()
    return __object_attr(
        obj, mode, keys_to_skip,
        lambda x: not isinstance(x, (MethodType, FunctionType))
    )


#def remove_files(*filenames):
    #"""delete a list of files"""
    #failed_list = []
    #for filename in filenames:
        #try:
            #os.remove(filename)
        #except OSError:  # OSError is the general version of WindowsError
            #failed_list.append(filename)
    #return failed_list


def int_version(name, version):
    """splits the version into a tuple of integers"""
    sversion = version.split('-')[0]
    #numpy
    #scipy
    #matplotlib
    #qtpy
    #vtk
    #cpylog
    #pyNastran
    if 'rc' not in name:
        # it's gotta be something...
        # matplotlib3.1rc1
        sversion = sversion.split('rc')[0]

    try:
        return [int(val) for val in sversion.split('.')]
    except ValueError:
        raise SyntaxError('cannot determine version for %s %s' % (name, sversion))
