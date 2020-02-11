"""
defines:
 - deprecated(old_name, new_name, deprecated_version, levels=None)
 - print_bad_path(path)
 - object_attributes(obj, mode='public', keys_to_skip=None)
 - object_methods(obj, mode='public', keys_to_skip=None)
"""
# -*- coding: utf-8 -*-
from types import MethodType, FunctionType
import os
import io
import sys
import getpass
import inspect
import warnings
from abc import abstractmethod
from typing import List, Optional, Any
import pyNastran


def ipython_info() -> Optional[str]:
    """determines if iPython/Jupyter notebook is running"""
    try:
        return get_ipython()
    except NameError:
        return None

def is_file_obj(filename: str) -> bool:
    """does this object behave like a file object?"""
    return (
        (hasattr(filename, 'read') and hasattr(filename, 'write'))
        or isinstance(filename, (io.IOBase, io.StringIO))
    )

def b(string: str) -> bytes:
    """reimplementation of six.b(...) to work in Python 2"""
    return string.encode('latin-1')

#def merge_dicts(dict_list, strict: bool=True):
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

def remove_files(filenames):
    """remvoes a series of files; quietly continues if the file can't be removed"""
    for filename in filenames:
        try:
            os.remove(filename)
        except OSError:
            pass

def is_binary_file(filename: str) -> bool:
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
    assert isinstance(filename, str), '%r is not a valid filename' % filename
    check_path(filename)
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b'\0' in chunk:  # found null byte
                return True
    return False


def check_path(filename: str, name: str='file') -> None:
    """checks that the file exists"""
    try:
        exists = os.path.exists(filename)
    except TypeError:
        msg = 'cannot find %s=%r\n' % (name, filename)
        raise TypeError(msg)
    if not exists:
        msg = 'cannot find %s=%r\n%s' % (name, filename, print_bad_path(filename))
        raise FileNotFoundError(msg)

def print_bad_path(path: str) -> str:
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

    path = os.path.abspath(path)
    npath = os.path.dirname(path)
    res = [path]
    while path != npath:
        path, npath = npath, os.path.dirname(npath)
        res.append(path)
    msg = {True: 'passed', False: 'failed'}
    return '\n'.join(['%s: %s' % (msg[os.path.exists(i)], i) for i in res])

def _filename(filename: str) -> str:
    """
    Prepends some magic data to a filename in order to have long filenames.

    .. warning:: This might be Windows specific.
    """
    if len(filename) > 255:
        return '\\\\?\\' + filename
    return filename

def __object_attr(obj, mode, keys_to_skip, attr_type, filter_properties: bool=False):
    """list object attributes of a given type"""
    #print('keys_to_skip=%s' % keys_to_skip)
    keys_to_skip = [] if keys_to_skip is None else keys_to_skip
    test = {
        'public':  lambda k: (not k.startswith('_') and k not in keys_to_skip),
        'private': lambda k: (k.startswith('_') and not k.startswith('__')
                              and k not in keys_to_skip),
        'both': lambda k: (not k.startswith('__') and k not in keys_to_skip),
        'all':  lambda k: (k not in keys_to_skip),
    }

    if not mode in test:
        raise ValueError('Wrong mode! Accepted modes: public, private, both, all.')
    check = test[mode]

    out = []
    obj_type = type(obj)
    for key in dir(obj):
        #if isinstance(key, abstractmethod):  # doesn't work...
            #print(key, '  abstractmethod')
        #if isinstance(key, FunctionType):
            #print(key, '  FunctionType')
        #if isinstance(key, MethodType):
            #print(key, '  MethodType')

        if key in keys_to_skip or not check(key):
            continue

        try:
            value = getattr(obj, key)
            save_value = attr_type(value)
            if not save_value:
                continue
            if filter_properties:
                if not isinstance(getattr(obj_type, key, None), property):
                    out.append(key)
            else:
                out.append(key)
        except:
            pass
    out.sort()
    return out
    #return sorted([k for k in dir(obj) if (check(k) and
    #                                           attr_type(getattr(obj, k)))])


def object_methods(obj: Any, mode: str='public',
                   keys_to_skip: Optional[List[str]]=None) -> List[str]:
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

def object_stats(obj: Any, mode: str='public',
                 keys_to_skip: Optional[List[str]]=None,
                 filter_properties: bool=False) -> str:
    """Prints out an easy to read summary of the object"""
    msg = '%s:\n' % obj.__class__.__name__
    attrs = object_attributes(
        obj, mode=mode, keys_to_skip=keys_to_skip,
        filter_properties=filter_properties)

    for name in sorted(attrs):
        #if short and '_ref' in name:
            #continue
        value = getattr(obj, name)
        msg += '  %-6s : %r\n' % (name, value)
    return msg

def object_attributes(obj: Any, mode: str='public',
                      keys_to_skip: Optional[List[str]]=None,
                      filter_properties: bool=False) -> List[str]:
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
    filter_properties: bool: default=False
        filters the @property objects

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
        lambda x: not isinstance(x, (MethodType, FunctionType)),
        filter_properties=filter_properties,
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


def int_version(name: str, version: str) -> List[int]:
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


def deprecated(old_name: str, new_name: str, deprecated_version: str,
               levels: Optional[List[int]]=None) -> None:
    """
    Throws a deprecation message and crashes if past a specific version.

    Parameters
    ----------
    old_name : str
        the old function name
    new_name : str
        the new function name
    deprecated_version : float
        the version the method was first deprecated in
    levels : List[int]
        the deprecation levels to show
        [1, 2, 3] shows 3 levels up from this function (good for classes)
        None : ???

    TODO: turn this into a decorator?

    """
    assert isinstance(deprecated_version, str), type(deprecated_version)
    assert isinstance(levels, list), type(levels)
    assert old_name != new_name, "'%s' and '%s' are the same..." % (old_name, new_name)

    version = pyNastran.__version__.split('_')[0]
    dep_ver_tuple = tuple([int(i) for i in deprecated_version.split('.')])
    ver_tuple = tuple([int(i) for i in version.split('.')[:2]])

    #new_line = ''
    msg = "'%s' was deprecated in v%s (current=%s)" % (
        old_name, deprecated_version, version)
    if new_name:
        msg += "; replace it with '%s'\n" % new_name

    for level in levels:
        # jump to get out of the inspection code
        frame = sys._getframe(3 + level)
        line_no = frame.f_lineno
        code = frame.f_code
        try:
            #filename = os.path.basename(frame.f_globals['__file__'])
            filename = os.path.basename(inspect.getfile(code))
        except:
            print(code)
            raise

        source_lines, line_no0 = inspect.getsourcelines(code)
        delta_nlines = line_no - line_no0
        try:
            line = source_lines[delta_nlines]
        except IndexError:
            break
        msg += '  %-25s:%-4s %s\n' % (filename, str(line_no) + ';', line.strip())

    user_name = getpass.getuser()
    if ver_tuple > dep_ver_tuple: # or 'id' in msg:
        # fail
        raise NotImplementedError(msg)
    elif user_name not in ['sdoyle', 'travis']:
        warnings.warn(msg, DeprecationWarning)
