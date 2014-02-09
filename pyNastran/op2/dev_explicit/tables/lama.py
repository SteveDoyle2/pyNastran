from struct import Struct

from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues)


class LAMA(object):
    def __init__(self):
        pass

    def _read_complex_eigenvalue_3(self, data):
        self.show_data(data)
        aaa

    def _read_complex_eigenvalue_4(self, data):
        self.show_data(data)
        bbb

    def _read_real_eigenvalue_3(self, data):
        #self.show_data(data)

        (three) = self.parse_approach_code(data)

        self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.add_data_parameter(data, 'fldFlag', 'i', 12, False)

        #print self.data_code
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code
        #self.add_data_parameter(data,'num_wide',     'i',10,False)  ## number of words per entry in record; .. note:: is this needed for this table ???

        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #self.print_block(data)
        self._read_title(data)


    def _read_real_eigenvalue_4(self, data):
        #self.show_data(data)
        nModes = len(data) // 28
        n = 0
        ntotal = 28
        lama = RealEigenvalues(self.isubcase)
        self.eigenvalues[self.isubcase] = lama
        s = Struct('ii5f')
        for i in xrange(nModes):
            edata = data[n:n+28]
            out = s.unpack(edata)
            if self.debug4():
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(iMode,order,eigen,omega,freq,mass,stiff) = out
            #(modeNum,extractOrder,eigenvalue,radian,cycle,genM,genK) = line
            #print out
            lama.addF06Line(out)
            n += ntotal
        return n
