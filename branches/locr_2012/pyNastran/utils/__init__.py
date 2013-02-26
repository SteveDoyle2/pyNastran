# -*- coding: utf-8 -*-
from types import MethodType
import os
from os.path import splitext, getsize
from os.path import join as pjoin
from numpy import ndarray
import io

def is_binary(filename):
    """
    Return true if the given filename is binary.
    Based on the idea (see http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text) that file is binary if it contains null.
    
    @raise IOError if the file cannot be opened.
    @retval True if filename is a binary file (contains null byte) and False otherwise.
    @warning this may not work for unicode.
    """
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b"\0" in chunk:  # found null byte
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
    def _print(val):
        if val is None or isinstance(val, basestring):
            return str(val)
        if isinstance(val, float):
            return '%-4.2f' % (val)
        try:
            return '%g' % (val)
        except TypeError:
            print("parameter = |%s|" % (val))
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
    

def write_object_attributes(name, obj, nspaces=0, nbase=0, isClass=False):
    """
    writes a series of nested objects
    """
    spaces = (nbase+nspaces) * ' '
    msg = spaces
    if isClass:
        equals = '='
    else:
        equals = ':'

    ## name
    if isinstance(obj, dict):
        if nspaces == 0:
            msg += '%s = {\n' % name
        else:
            if isinstance(name, tuple):
                msg += "%s : {\n" % str(name)
            else:
                msg += "'%s' : {\n" % name
    elif isinstance(name, str) or isinstance(name, unicode):
        if isClass:
            key = '%s' % name
        else:
            key = "'%s'" % name
    
    elif isinstance(name, int) or isinstance(name, float) or isinstance(name, tuple) or name is None:
        key = "%s" % str(name)
    else:
        raise RuntimeError('key=%s is not a string.  Type=%s' % (name, type(name)))
    
    
    ## write the object
    if isinstance(obj, int) or isinstance(obj, float) or obj is None:
        msg += '%s %s %s,\n' % (key, equals, str(obj))
    elif isinstance(obj, str) or isinstance(obj, unicode):
        msg += "%s %s '%s',\n" % (key, equals, obj)

    elif isinstance(obj, dict):
        nspaces2 = nspaces + 4
        msg += write_dict(obj, spaces, nspaces, nspaces2, nbase, isClass)
    elif isinstance(obj, tuple):
        msg += '%s : %s,\n' % (key, str(obj))
    elif isinstance(obj, list):
        msg += '%s : [' % (key)
        for value in obj:
            msg += write_value(value) + ', '
        msg += '],\n'
        
    elif isinstance(obj, ndarray):
        #msg += '%s %s %s,\n' % (key, equals, write_array(obj))
        msg += '%s %s %s,\n' % (key, equals, write_array(obj))
    else:  # generic class
        objectType = obj.__class__.__name__
        #raise RuntimeError('objectType=%s is not supported' % objectType)
        obj_attrs = object_attributes(obj, 'public')
        
        nspaces2 = nspaces + 4
        spaces2 = nspaces2 * ' '
        msg += "'%s' : %s(\n" % (name, objectType)
        for attr in obj_attrs:
            value = getattr(obj, attr)
            msg += write_object_attributes(attr, value, nspaces2, nbase, isClass=True)
        if nspaces == 0: # bottom level, no comma
            msg += '%s)\n'  % spaces
        else:  # embedded dictionary
            msg += '%s),\n'  % spaces

        #print "dir(obj) =", dir(obj)
        #print "obj_attrs =", obj_attrs
    return msg

def write_dict(obj, spaces, nspaces, nspaces2, nbase, isClass):
    msg = ''
    spaces2 = nspaces2 * ' '
    for key2, value in sorted(obj.iteritems()):
        msg += write_object_attributes(key2, value, nspaces2, nbase, isClass=False)
    if nspaces == 0: # bottom level, no comma
        msg += '%s}\n'  % spaces
    else:  # embedded dictionary
        msg += '%s},\n'  % spaces
    return msg

def write_value(obj, nspaces=0):
    msg = ''
    if isinstance(obj, int) or isinstance(obj, float) or obj is None:
        msg += '%s' % (str(obj))
    elif isinstance(obj, str) or isinstance(obj, unicode):
        msg += "'%s'" % (obj)
    elif isinstance(obj, list):
        msg += write_list(obj)
    else:
        raise RuntimeError('objectType=%s is not supported; value=%s' % objectType, obj)
    return msg    

def write_list(obj):
    msg = '['
    for value in obj:
        msg += write_value(value) + ', '
    msg += ']'
    return msg

def write_array(a):
    shape = a.shape
    if len(shape)==1:
        msg = 'array(['
        for ai in a[:-1]:
            if not(isinstance(ai, int) or isinstance(ai, int) or isinstance(ai, float)):
                return 'array(.not supported type.)'
            msg += '%s, ' % ai
        msg += '%s])' % a[-1]
    else:
        return 'array(.not supported shape.)'
    return msg

if __name__ == '__main__':
    from numpy import array, zeros
    class A(object):
        def __init__(self, a=None, b=None, c=None):
            self.a = a
            self.b = b
            self.c = c

    print zeros((2,2))
    dictA = {
            'strString' : 'a string',
            'strFloat' : 1.0,
            'strInt': 2,
            'strTuple': (1,2),
            'strNone' : None,
            #'strClass' : A('a', 'b', 'c'),
            'strList' : [1,2,3],
            1 : 1,
            None : 4,
            1.0 : 5,
            (1, 2) : 6,
            'strArray' : array([4,5,6]),
            'strArray2' : zeros((2,2)),
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
    #print write_object_attributes('dictA', dictA)
    msg = write_object_attributes('dictB', dictB, nbase=0)
    print msg
    
    
    #dictB2 = eval(msg)
    
