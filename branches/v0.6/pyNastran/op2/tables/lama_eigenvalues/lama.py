import sys
from struct import unpack

from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues)
#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class LAMA(object):

    def readTable_LAMA(self):
        table_name = self.read_table_name(rewind=False)  # LAMA
        self.table_init(table_name)
        #print "tablename1 = |%r|" %(table_name)
        #print "tablename2 = |%r|" %(self.table_name)

        self.read_markers([-1, 7], 'LAMA')
        ints = self.read_int_block()
        #print "*ints = ",ints

        self.read_markers([-2, 1, 0], 'LAMA')
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words

        word = self.read_string_block()  # LAMA
        #print "word = |%s|" %(word)
        self.read_markers([-3, 1, 0], 'LAMA')

        #data = self.get_data(4*50)
        #print self.print_block(data)

        self.readTable_LAMA_3(-3)

        self.read_markers([-4, 1, 0], 'LAMA')
        self.readTable_LAMA_4(-4)

        self.read_markers([-5, 1, 0], 'LAMA')
        #data = self.get_data(4*30)
        #print self.print_block(data)
        return
        sys.exit('stopping in LAMA')
        if 0:
            iTable = -3
            #imax   = -244

            while buffer_words:  # read until buffer_words=0
                self.read_markers([iTable, 1, 0], 'LAMA')
                nOld = self.n
                buffer_words = self.get_marker()
                #print "buffer_words = ",buffer_words
                if buffer_words == 0:  # maybe read new buffer...
                    self.goto(nOld)
                    break
                data = self.read_block()
                #print "len(data) = ",len(data)
                self.readDesvar(data)
                iTable -= 1
            self.print_section(80)

        #self.op2Debug.write('buffer_words=%s\n' %(str(buffer_words)))
        #print "1-buffer_words = ",buffer_words,buffer_words*4

        #print self.print_section(300)
        #sys.exit('asdf')

    def readTable_LAMA_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2_debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

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
        self.read_title()

    def readTable_LAMA_4(self, iTable):  # iTable=-4
        buffer_words = self.get_marker()  # 70*4=280
        if self.make_op2_debug:
            self.op2_debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)  # dummy - 70*4=280
        #print self.print_block(data)
        #print "280/3 = ",280/4
        nModes = buffer_words // 7

        lama = RealEigenvalues(self.isubcase)
        self.eigenvalues[self.isubcase] = lama
        for i in xrange(nModes):
            data = self.get_data(28)  # 4*7
            out = unpack('iifffff', data)
            #(iMode,order,eigen,omega,freq,mass,stiff) = out
            #(modeNum,extractOrder,eigenvalue,radian,cycle,genM,genK) = line
            #print out
            lama.addF06Line(out)
            #print "mode=%s order=%s eigen=%s omega=%s freq=%s mass=%s stiff=%s" %(mode,order,eigen,omega,freq,mass,stiff)
        #print ""
        #print ''.join(msg)
        #print "self.isubcase = ",self.isubcase
        #print lama.write_f06([],'PAGE',1)[0]
        #sys.exit()
#                       '        1         1        8.232776E+06        2.869281E+03        4.566603E+02        8.719168E-03        7.178296E+04
#                       '        2         2        8.232776E+06        2.869281E+03        4.566603E+02        8.719168E-03        7.178296E+04

        data = self.get_data(4)
