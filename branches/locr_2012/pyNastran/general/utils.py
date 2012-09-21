# -*- coding: utf-8 -*-

from types import MethodType

def is_binary(filename):
    """
    Return true if the given filename is binary.
    Based on the idea that file is binary if it contains null. See 
    http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text on 6/08/2010
    @raise IOError if the file cannot be opened.
    @warning this may not work for unicode...
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
    @param debug display additional information about conversion process
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
    @param debug display additional information about conversion process
    """
    
    dict_vals = dict(zip(list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'),xrange(52)))
    val = 0
    for i, letter in enumerate(list(num)):
        val += dict_vals[letter] * 52 ** i
        if debug:
            print("letter = ", letter, "\nfactor = ", dict_vals[letter] * 52 ** i)
    return val

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
    
