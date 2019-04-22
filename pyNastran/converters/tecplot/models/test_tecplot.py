# coding: utf-8
from __future__ import  print_function
import os
import sys
import unittest

from pyNastran.converters.tecplot.tecplot import read_tecplot
from cpylog import get_logger2


#class TestTecplotAscii(unittest.TestCase):
    #def test_ascii(self):
        #log = get_logger2(debug=False)
        #dirname = 'ascii'
        #fnames = [os.path.join(dirname, fname) for fname in os.listdir(dirname)
                  #if not fname.endswith('.png')]
        #for fname in fnames:
            #try:
                #read_tecplot(fname, log=log)
                #log.info('read %r' % fname)
            #except Exception as error:
                #log.warning('failed reading %r' % fname)
                #log.error(error)
                #print('')

class TestTecplotBinary(unittest.TestCase):
    def test_binary(self):
        log = get_logger2(debug=True)
        fnames = [os.path.join('binary', fname) for fname in os.listdir('binary')
                  if not fname.endswith('.png')]
        print(fnames)
        for fname in fnames:
            print(fname)
            try:
                read_tecplot(fname, log=log)
                log.info('read %r' % fname)
            except Exception as error:
                sys.stdout.flush()
                raise
                #log.warning('failed reading %r' % fname)
                #log.error(error)
                #print('')

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

