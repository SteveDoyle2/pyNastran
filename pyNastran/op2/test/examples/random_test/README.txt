Attempting to parse PSDF, ATOC, RMS, and CRMS tables from a random analysis (modal frequency response).

1. run nastran on random_test.dat
2. in python do this:

from pyNastran.op2.op2 import OP2
op2 = OP2()
op2.read_op2('random_test.op2')

I'm getting the following error:

DEBUG:     fname=op2.py                    lineNo=375    combine=True
DEBUG:     fname=op2.py                    lineNo=376    -------- reading op2 with read_mode=1 --------
INFO:      fname=op2_scalar.py             lineNo=1177   op2_filename = 'random_test.op2'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='PVT0'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='CASECC'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='EQEXINS'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='OGPWG'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='FRL'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='KHH'
WARNING:   fname=op2_scalar.py             lineNo=1707   created KHH...verify the complex matrix
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='LAMA'
DEBUG:     fname=op2_scalar.py             lineNo=1368     table_name='OUGPSD1'

---------------------------------------------------------------------------
NotImplementedError                       Traceback (most recent call last)
<ipython-input-2-aaa8d394c671> in <module>()
      1 op2 = OP2()
----> 2 op2.read_op2('random_test.op2')

C:\Anaconda\lib\site-packages\pynastran-0.9.0+dev.no.checksum.error-py2.7.egg\pyNastran\op2\op2.pyc in read_op2(self, op2_filename, combine, build_dataframe, skip_undefined_matrices, encoding)
    379 
    380         # get GUI object names, build objects, but don't read data
--> 381         OP2_Scalar.read_op2(self, op2_filename=op2_filename)
    382 
    383         # TODO: stuff to figure out objects

C:\Anaconda\lib\site-packages\pynastran-0.9.0+dev.no.checksum.error-py2.7.egg\pyNastran\op2\op2_scalar.pyc in read_op2(self, op2_filename, combine)
   1299 
   1300         self._make_tables()
-> 1301         table_names = self._read_tables(table_name)
   1302         if self.is_debug_file:
   1303             self.binary_debug.write('-' * 80 + '\n')

C:\Anaconda\lib\site-packages\pynastran-0.9.0+dev.no.checksum.error-py2.7.egg\pyNastran\op2\op2_scalar.pyc in _read_tables(self, table_name)
   1412                     msg += 'If you have matrices that you want to read, see:\n'
   1413                     msg += '  model.set_additional_matrices(matrices)'
-> 1414                     raise NotImplementedError(msg)
   1415 
   1416             table_name = self._read_table_name(rewind=True, stop_on_failure=False)

NotImplementedError: geom/results split: 'OUGPSD1'

If you have matrices that you want to read, see:
  model.set_additional_matrices(matrices)


