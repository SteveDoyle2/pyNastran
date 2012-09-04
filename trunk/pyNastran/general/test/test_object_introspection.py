# -*- coding: utf-8 -*-

import unittest

from pyNastran.general.object_introspection import (list_methods,
                                                    list_private_methods,
                                                    list_attributes,
                                                    list_private_attributes)

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
        self.d = E()

    def getB(self):
        return self.b

    def _getB(self):
        return self.b

class E(object):
    def __init__(self):
        self.e = 3

    def getD(self):
        return self.d

class TestObjectIntrospection(unittest.TestCase):

    def test_object_introspection_1(self):
        b = B(7)
        methods = list_methods(b)
        _methods = list_private_methods(b)
        attributes = list_attributes(b)
        _attributes = list_private_attributes(b)

        self.assertEqual(methods,  ['getA', 'getB'])
        self.assertEqual(_methods, ['_getA', '_getB'])
        
        self.assertEqual(attributes,  ['a', 'b', 'c', 'd'])
        self.assertEqual(_attributes, ['_a', '_b'])
        
if __name__ == "__main__":
    unittest.main()
