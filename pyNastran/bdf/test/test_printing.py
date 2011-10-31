import unittest

from test_helper import runBDF
from pyNastran.bdf.fieldWriter import printCard,printField

class Tester(unittest.TestCase):
    def test_ints(self):
        self.assertEquals(printField( 1        ),'       1','a')
        self.assertEquals(printField( 12345678 ),'12345678','b')
        self.assertEquals(printField('12345678'),'12345678','c')
        #self.assertEquals(printField('1       '),'       1','|%s|' %(printField('1       ')))

    def test_floats(self):
        self.assertEquals(printField( 1.2        ),'     1.2','a')
        self.assertEquals(printField( 1.234567890),'1.234568','b')
        self.assertEquals(printField( 12.234568  ),'12.23457','c')
        self.assertEquals(printField( 123.23457  ),'123.2346','d')
        self.assertEquals(printField( 1234.23468 ),'1234.235','e')
        self.assertEquals(printField( 12345.238  ),'12345.24','f')
        self.assertEquals(printField( 123456.28  ),'123456.3','g')
        self.assertEquals(printField( 1234567.25 ),'1234567.',printField( 1234567.25 ))  # 1.2346+6
        self.assertEquals(printField( 12345678. ), '1.2346+7','|%s|'%printField( 12345678. ))
        self.assertEquals(printField( 123456789.), '1.2346+8','|%s|'%printField( 12345678. ))
        
if __name__=='__main__':
    unittest.main()
    #print printField( 1234567. )
