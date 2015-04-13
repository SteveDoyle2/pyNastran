# -*- coding: utf-8 -*-
from six import string_types
from types import MethodType
import os
from numpy import ndarray
import io


def is_binary_file(filename):
    """
    Return true if the given filename is binary.

    :raises:  IOError if the file cannot be opened.
    :returns: True if filename is a binary file (contains null byte)
              and False otherwise.

    Based on the idea (.. seealso:: http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text)
    that file is binary if it contains null.

    .. warning:: this may not work for unicode."""
    assert isinstance(filename, string_types), '%r is not a valid filename' % filename
    assert os.path.exists(filename), '%r does not exist' % filename
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b"\0" in chunk:  # found null byte
                return True
    return False


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


def object_attributes(obj, mode="public"):
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
