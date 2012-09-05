# -*- coding: utf-8 -*-

from types import MethodType

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
    
