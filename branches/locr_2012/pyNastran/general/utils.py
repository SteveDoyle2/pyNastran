# -*- coding: utf-8 -*-

from types import MethodType
import os
from os.path import splitext
from os.path import join as pjoin
from os.path import getsize
from numpy import ndarray

def is_binary(filename):
    """
    Return true if the given filename is binary.
    Based on the idea (see http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text) that file is binary if it contains null.
    
    @raise IOError if the file cannot be opened.
    @retval True if filename is a binary file (contains null byte) and False otherwise.
    @warning this may not work for unicode.
    """
    with open(filename, 'rb') as fil:
        for chunk in iter(lambda:fil.read(1024), ''):
            if b'\0' in chunk:  # found null byte
                return True
    return False

def obscure(num, debug=False):
    """
    Takes a large positive number and shrinks it down...similar to binary, but base 52.
    A base 52 value takes up a fewer characters than a base 10 number
    which helps to do Mat12345678 when there's an 8 character limit on variable names.
    @code
    >>> obscure(35214)
    'kbn'
    >>> de_obscure('kbn')
    35214
    @endcode
    @param num positive integer number
    @param debug display additional information about conversion process
    @retval shortened version of num
    """
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

def de_obscure(num, debug = False):
    """
    Unpacks an "obscured" positive number...similar to binary, but base 52.
    A base 52 value takes up a fewer characters than a base 10 number
    which helps to do Mat12345678 when there's an 8 character limit on variable names.
    @code
    >>> obscure(35214)
    'kbn'
    >>> de_obscure('kbn')
    35214
    @endcode
    @param debug display additional information about conversion process
    @retval integer value of shortened version of a number, @see obscure
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
    @param dirname the directory name
    @param extension list of filetypes to get (default='.txt')
    @param maxSize size in MB for max file size
    @retval list of all the files with a given extension in the specified directory 
    """

    return [pjoin(dirname, f) for f in os.listdir(dirname) if extension in 
            splitext(f)[1] and getsize(pjoin(dirname, f)) / (1048576.) <= maxSize]

def print_bad_path(path):
    """
    Prints information about the existence (access possibility) of the parts of the given path.
    Useful for debugging when the path to a given file is wrong.
    @param path path to check
    @retval string with informations whether access to parts of the path is possible
    """
    path = os.path.abspath(path)
    npath = os.path.dirname(path)
    res = [path]
    while path != npath:
        path, npath = npath, os.path.dirname(npath)
        res.append(path)
    msg = {True: "passed", False: "failed"}
    return "\n".join(["%s: %s" % (msg[os.path.exists(i)], i) for i in res])

def list_print(lst):
    """
    Prints a list, numpy array, or numpy matrix in an abbreviated format.
    Supported element types: None, string, numbers. Useful for debugging.
    @param lst list, numpy array or numpy matrix
    @retval the clean string representation of the object
    """
    if len(lst) == 0:
        return '[]'

    if isinstance(lst, ndarray) and lst.ndim == 2:
        r,c = lst.shape
        return ("["+",\n ".join(["["+",".join(['%-10g' % lst[i, j] 
                            for j in range(c)])+"]" for i in range(r)])+"]")

    def _print(a):
        if a is None or isinstance(a, basestring):
            return str(a)
        if isinstance(a, float):
            return '%-4.2f' % (a)
        try:
            return '%g' % (a)
        except TypeError:
            print("parameter = |%s|" % (a))
            raise
    return "[" + ", ".join([_print(a) for a in lst]) + "]"

# list object attributes of a given type 
def __object_attr(obj, mode, attr_type):
    test = {"public":  lambda k: not k.startswith('_'),
            "private": lambda k: (k.startswith('_') and not k.startswith('__')),
            "both": lambda k: not k.startswith('__'),
            "all":  lambda k: True  
            }

    if not mode in test:
        print("Wrong mode! Accepted modes: public, private, both, all.")
        return None
    check = test[mode]
    
    return sorted([k for k in dir(obj) if (check(k) and  
                                               attr_type(getattr(obj, k)))])
    
def object_methods(obj, mode = "public"):
    """
    List the names of methods of a class as strings. Returns public methos
    as default.
    
    @param obj
      the object for checking
    @param mode
      defines what kind of methods will be listed
       
       * "public" - names that do not begin with underscore
       * "private" - names that begin with single underscore
       * "both" - private and public
       * "all" - all methods that are defined for the object 
    @retval 
      sorted list of the names of methods of a given type or None if the mode
      is wrong
    """
    return __object_attr(obj, mode, lambda x: isinstance(x, MethodType))
    

def object_attributes(obj, mode = "public"):
    """
    List the names of attributes of a class as strings. Returns public methos
    as default.
    
    @param obj
      the object for checking
    @param mode 
      defines what kind of attributes will be listed
       
       * "public" - names that do not begin with underscore
       * "private" - names that begin with single underscore
       * "both" - private and public
       * "all" - all attributes that are defined for the object 
    @retval 
      sorted list of the names of attributes of a given type or None if the mode
      is wrong
    """
    return __object_attr(obj, mode, lambda x: not isinstance(x, MethodType))
    
