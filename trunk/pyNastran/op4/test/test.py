import os
#from pyNastran.op4.cop4 import OP4
#print "f = ",op4.__file__
import unittest
from pyNastran.op4.op4 import OP4


def pass_test1():
    fh = OP4(os.path.abspath('mat_b_dn.op4'),'r')
    fh.print_header()
    #print fh.nmat = 9
    
    # crash
    a,b,c = fh.Load(nmat=3,skip=0)
    print a
    print b
    print c

def failed_test1():
    fh = OP4('mat_b_dn.op4','r')
    fh.print_header()
    #print fh.nmat = 9
    
    a,b,c = fh.Load(nmat=3,skip=0)
    print a
    print b
    print c

def pass_test2():
    fh = OP4(os.path.abspath('mat_b_dn.op4'),'r')
    #print fh.nmat = 9

    # crash with "unnamed is sparse, skipping for now"
    (a,b,c,d, f,g,h,i) = fh.Load(nmat=9,skip=0) 
    print a
    print b
    print c

def failed_test2():
    fh = OP4('mat_b_dn.op4','r')
    #print fh.nmat = 9

    # ValueError:  need more than 8 values to unpack
    (a,b,c,d, f,g,h,i) = fh.Load(nmat=9,skip=0) 
    print a
    print b
    print c

def testAscii():
    print "ascii"
    for fname in ['mat_t_dn.op4',
                  'mat_t_s1.op4',
                  'mat_t_s2.op4',
                  ]:
        op4 = OP4()
        matrices = op4.readOP4(fname)
        for name,matrix in sorted(matrices.items()):
            print "name = %s" %(name)
            print matrix

def testBinary():
    print "binary"
    for fname in ['mat_b_dn.op4' ,
                  'mat_t_s1.op4',
                  'mat_t_s2.op4',
                  ]:
        op4 = OP4()
        matrices = op4.readOP4(fname)
        for name,matrix in sorted(matrices.items()):
            print "name = %s" %(name)
            print matrix

class OP4_Test(unittest.TestCase):
    def op4_1(self):
        testBinary()
    def op4_2(self):
        testAscii()

print "-------------"
if __name__ == '__main__':
    #failed_test1()
    #print "*********"
    #failed_test2()
    
    #pass_test1()
    #pass_test2()
    testAscii()
    testBinary()
    unittest.main()