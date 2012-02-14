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
    
    def test_2(self):
        for i in range(100000):
            a = random.uniform(-20,20)
            a2 = 10**a
            self.compare(a2)
            self.compare(-a2)

    def compare(self,valueIn):
        #print "a = |%s|" %(valueIn)
        field = printField(valueIn)
        print "a = |%s|" %(field)

        val = getValue(field)
        if val != 0:
            p = (val-valueIn)/val
            if p>0.01:
                raise Exception('val=%s valueIn=%s' %(val,valueIn))
        #assert allclose(val,valueIn),

if __name__ == "__main__":
    unittest.main()
