# -*- coding: utf-8 -*-

import unittest
import sys


from pyNastran.general.utils import (object_methods, object_attributes)

class A(object):
    def __init__(self):
        self.a = 5
        self._a = self.a**2

    def getA(self):
        return self.a

    def _getA(self):
        return self.a

class B(A):
    c = 7
    def __init__(self, b):
        A.__init__(self)
        self.b = b
        self._b = b**2

    def getB(self):
        return self.b

    def _getB(self):
        return self.b

class TestObjectIntrospection(unittest.TestCase):
    
    def setUp(self):
        self.b = B(7)

    def test_object_methods_introspection(self):
        methods = object_methods(self.b)
        self.assertEqual(methods,  ['getA', 'getB'])
        
        methods = object_methods(self.b, "private")
        self.assertEqual(methods, ['_getA', '_getB'])
        
        methods = object_methods(self.b, "both")
        self.assertEqual(methods, ['_getA', '_getB', 'getA', 'getB'])
                         
        methods = object_methods(self.b, "all")
        self.assertEqual(methods, ['__init__', '_getA', '_getB', 'getA', 
                                    'getB'])
        
    def test_object_attributes_introspection(self):
        
        attributes = object_attributes(self.b)
        self.assertEqual(attributes, ['a', 'b', 'c'])
        
        
        attributes = object_attributes(self.b, "private")
        self.assertEqual(attributes, ['_a', '_b'])
        
        attributes = object_attributes(self.b, "both")
        self.assertEqual(attributes, ['_a', '_b', 'a', 'b', 'c'])
        
        
    @unittest.skipIf(sys.version_info >= (3,0), "est for Python 2.x")
    def test_object_attributes_introspection_2(self):
        attributes = object_attributes(self.b, "all")
        self.assertEqual(attributes, ['__class__', '__delattr__', '__dict__',
                '__doc__', '__format__', '__getattribute__', '__hash__', 
                '__module__', '__new__', '__reduce__', '__reduce_ex__', 
                '__repr__', '__setattr__', '__sizeof__', '__str__', 
                '__subclasshook__', '__weakref__', '_a', '_b', 'a', 'b', 'c'])
        
    @unittest.skipIf(sys.version_info < (3,0), "est for Python 3.x")
    def test_object_attributes_introspection_3(self):
        attributes = object_attributes(self.b, "all")
        self.assertEqual(attributes, ['__class__', '__delattr__', '__dict__', 
                '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__',
                '__gt__', '__hash__', '__le__', '__lt__', '__module__', 
                '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__',
                '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 
                '__weakref__', '_a', '_b', 'a', 'b', 'c'])
        
if __name__ == "__main__":
    unittest.main()
