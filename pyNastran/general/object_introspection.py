# -*- coding: utf-8 -*-
"""
gets information about objects
"""
from types import MethodType

def list_methods(obj):
    """
    Lists the names of the public methods of a class as strings.  A public
    method is defined as having an underscore at the beginning of it.
    """
    methods = []
    for key in dir(obj):
        val = getattr(obj, key)
        if not key.startswith('_'):
            if isinstance(val, MethodType):
                methods.append(key)
    return sorted(methods)

def list_private_methods(obj):
    """
    Lists the names of the private methods of a class as strings.  A private
    method is defined as having an underscore at the beginning of it.
    """
    methods = []
    for key in dir(obj):
        val = getattr(obj, key)
        if key.startswith('_') and not key.startswith('__'):
            if isinstance(val, MethodType):
                methods.append(key)
    return sorted(methods)


def list_attributes(obj):
    """
    Lists the names of public attributes of a class as strings.  A public
    attribute is defined as not having an underscore at the beginning of it.
    """
    data = []
    for key in dir(obj):
        val = getattr(obj, key)
        if('_' not in key and not isinstance(val, MethodType)):
            data.append(key)
    return sorted(data)


def list_private_attributes(obj):
    """
    Lists the names of private attributes of a class.  A private
    attribute is defined as having an underscore at the beginning of it.
    """
    data = []
    for key in dir(obj):
        val = getattr(obj, key)
        if not key.startswith('__') and key.startswith('_'):
            if not isinstance(val, MethodType):
                data.append(key)
    return sorted(data)