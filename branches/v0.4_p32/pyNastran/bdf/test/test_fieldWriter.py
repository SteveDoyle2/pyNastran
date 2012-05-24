## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import random
import unittest

#from numpy import allclose

from pyNastran.bdf.fieldWriter import printField
from pyNastran.bdf.bdfInterface.bdf_cardMethods import getValue


class TestFieldWriter(unittest.TestCase):

    def test_1(self):
        field = printField(1e20);       assert '   1.+20' == field,'|%s|' %(field)
        field = printField(-.723476);   assert '-.723476' == field,'|%s|' %(field)
        field = printField(125000. );   assert ' 125000.' == field,'|%s|' %(field)
        field = printField(12500000.);  assert '  1.25+7' == field,'|%s|' %(field)
        field = printField(47.77267);   assert '47.77267' == field,'|%s|' %(field)
        field = printField(.001);       assert '    .001' == field,'|%s|' %(field)
        field = printField(.0000001);   assert '.0000001' == field,'|%s|' %(field)
        field = printField(-5.007e-3);  assert '-5.007-3' == field,'|%s|' %(field)
        field = printField(-0.0748662); assert '-.074866' == field,'|%s|' %(field)
    
    def a_test_2(self):
        for i in range(1000):
            a = random.uniform(-20,20)
            a2 = 10**a
            self.compare(a2)
            self.compare(-a2)

    def test_ints(self):
        self.assertEquals(printField( 1        ),'       1','a')
        self.assertEquals(printField( 12345678 ),'12345678','b')
        self.assertEquals(printField('12345678'),'12345678','c')
        #self.assertEquals(printField('1       '),'       1','|%s|' %(printField('1       ')))

    def test_floats_greater_than_1(self):
        self.assertEquals(printField( 1.2       ),'     1.2','a')
        self.assertEquals(printField( 1.23456789),'1.234568','b')
        self.assertEquals(printField( 12.234568 ),'12.23457','c')
        self.assertEquals(printField( 123.23457 ),'123.2346','d')
        self.assertEquals(printField( 1234.23468),'1234.235','e')
        self.assertEquals(printField( 12345.238 ),'12345.24','f')
        self.assertEquals(printField( 123456.28 ),'123456.3','g')
        self.assertEquals(printField( 1234567.25),'1234567.',printField( 1234567.25 ))  # 1.2346+6
        self.assertEquals(printField( 12345678. ),'1.2346+7','|%s|'%printField( 12345678. ))
        self.assertEquals(printField( 123456789.),'1.2346+8','|%s|'%printField( 12345678. ))

    def test_floats_small(self):
        self.assertEquals(printField( 0.1        ),'      .1','A|%s|'%printField( 0.1))
        self.assertEquals(printField( 0.0001     ),'   .0001','B|%s|'%printField( 0.0001))
        self.assertEquals(printField( 0.00001    ),'  .00001','C|%s|'%printField( 0.00001))
        self.assertEquals(printField( 0.000001   ),' .000001','D|%s|'%printField( 0.000001))
        self.assertEquals(printField( 0.0000001  ),'.0000001','E|%s|'%printField( 0.0000001))
        self.assertEquals(printField( 0.00000012 ),'   1.2-7','F|%s|'%printField( 0.00000012))
        self.assertEquals(printField( 0.000748519),'7.4852-4','G|%s|'%printField( 0.000748519))
        self.assertEquals(printField( 0.12345678 ),'.1234568','H|%s|'%printField( 0.12345678))
        self.assertEquals(printField( 0.00012349 ),'1.2349-4','I|%s|'%printField( 0.00012349))
        self.assertEquals(printField( 0.000012349),'1.2349-5','J|%s|'%printField( 0.000012349))
        self.assertEquals(printField( 1e-20      ),'   1.-20','K|%s|'%printField( 1e-20))

    def compare(self,valueIn):
        #print "a = |%s|" %(valueIn)
        field = printField(valueIn)
        print("a = |%s|" %(field))

        val = getValue(field)
        if val != 0:
            p = (val-valueIn)/val
            if p>0.01:
                raise Exception('val=%s valueIn=%s' %(val,valueIn))
        #assert allclose(val,valueIn),

if __name__ == "__main__":
    unittest.main()
