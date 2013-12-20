# -*- coding: utf-8 -*-

import unittest
import sys
from numpy import matrix, array


from pyNastran.utils import (is_binary, obscure, de_obscure, list_print,
                                     object_methods, object_attributes)
from os.path import abspath

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


class TestUtils(unittest.TestCase):

    def setUp(self):
        self.b = B(7)

    def test_is_binary(self):
        self.assertTrue(is_binary(abspath("test_utils.pyc")))
        self.assertFalse(is_binary(abspath("test_utils.py")))

    def test_obscure(self):
        for num in [0,1,5,53,231123, 34567523, 1024, 65367, 14321324, 73123434,
                    1309872418439702897245, 955785585080806958106769948497824]:
            self.assertEqual(de_obscure(obscure(num)), num)

    def test_list_print(self):
        self.assertEqual(list_print(None), 'None')
        #self.assertRaises(TypeError, lambda: list_print(None))

        for a,b in [([],'[]'),(array([]), '[]'), (tuple(),'[]'), (matrix([]), '[[]]')]:
            self.assertEqual(list_print(a), b)
        r = ('[[1         ,2         ,3         ],\n [4         ,5         ,6'
             '         ],\n [7         ,8         ,9         ]]')
        self.assertEqual(list_print(array([(1,2,3),(4,5,6),(7,8,9)])), r)
        self.assertEqual(list_print(matrix([(1,2,3),(4,5,6),(7,8,9)])), r)
        self.assertEqual(list_print(array([(1.0,2,3.),(4.,5.,6),(7.0,8,9)])), r)
        self.assertEqual(list_print(matrix([(1,2,3.0),(4,5.0,6),(7.,8,9.0)])), r)

        r = "[[1.1       ,2.234     ,3.00001   ],\n [4.001     ,5         ,6.2       ]]"
        self.assertEqual(list_print(array([(1.1,2.234,3.00001),(4.001,5.0000005,6.2)])), r)
        self.assertEqual(list_print(matrix([(1.1,2.234,3.00001),(4.001,5.0000005,6.2)])), r)

        self.assertEqual(list_print(['a',None,11,'']), '[a, None, 11, ]')
        self.assertEqual(list_print(('a',None,11,'')), '[a, None, 11, ]')

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