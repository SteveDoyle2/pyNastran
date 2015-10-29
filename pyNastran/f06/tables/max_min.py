from pyNastran.f06.f06_classes import MaxDisplacement  # classes not in op2

class MAX_MIN(object):
    def __init__(self):
        pass

    def _get_max_spc_forces(self):  # .. todo:: not done
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_f06_table([int, float, float, float, float, float, float])
        #print "max SPC Forces   ",data
        #self.disp[isubcase] = DisplacementObject(isubcase,data)
        #print self.disp[isubcase]

    def _get_oload_resultant(self):  # .. todo:: not done
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_f06_table([int, float, float, float, float, float, float], debug=True)
        print(data)

    def _get_max_mpc_forces(self):  # .. todo:: not done
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_f06_table([int, float, float, float, float, float, float])
        #print "max SPC Forces   ", data
        #self.disp[isubcase] = DisplacementObject(isubcase, data)
        #print self.disp[isubcase]

    def _get_max_displacements(self):  # .. todo:: not done
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_f06_table([int, float, float, float, float, float, float])
        #print "max Displacements",data
        disp = MaxDisplacement(data)
        #print disp.write_f06()
        #self.disp[isubcase] = DisplacementObject(isubcase,data)
        #print self.disp[isubcase]

    def _get_max_applied_loads(self):  # .. todo:: not done
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_f06_table([int, float, float, float, float, float, float])
        #print "max Applied Loads",data
        #self.disp[isubcase] = DisplacementObject(isubcase,data)
        #print self.disp[isubcase]
