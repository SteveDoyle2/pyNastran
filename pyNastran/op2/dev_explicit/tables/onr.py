from struct import Struct

from pyNastran.op2.tables.oee_energy.oee_objects import StrainEnergyObject

class ONR(object):
    def __init__(self):
        pass

    def _read_onr1_3(self, data):  # TODO: this is wrong...
        self.read_opg1_3(data)

    def _read_onr1_4(self, data):
        if self.table_code == 18:  # element strain energy
            assert self.table_name in ['ONRGY1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_element_strain_energy(data)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_element_strain_energy(self, data):
        """
        table_code = 19
        """
        dt = self.nonlinear_factor
        n = 0
        if self.thermal == 0:
            result_name = 'strainEnergy'
            if self.num_wide == 4:
                self.create_transient_object(self.strainEnergy, StrainEnergyObject)
                s = Struct(b'i3f')

                ntotal = 16
                nnodes = len(data) // ntotal
                for i in xrange(nnodes):
                    eData = data[n:n+ntotal]

                    out = s.unpack(eData)
                    (eid_device, energy, percent, density) = out
                    eid = (eid_device - self.device_code) // 10
                    #print "eType=%s" %(eType)

                    dataIn = [eid, energy, percent, density]
                    #print "%s" %(self.get_element_type(self.element_type)),dataIn
                    self.obj.add(dt, dataIn)
                    n += ntotal
            elif self.num_wide == 5:
                self.create_transient_object(self.strainEnergy, StrainEnergyObject)  # why is this not different?
                s = Struct(b'i8s3f')
                ntotal = 20
                nnodes = len(data) // ntotal
                for i in xrange(nnodes):
                    eData = self.data[n:n+20]
                    out = unpack(format1, eData)
                    (word, energy, percent, density) = out
                    #print "out = ",out
                    word = word.strip()
                    #print "eType=%s" %(eType)

                    dataIn = [word, energy, percent, density]
                    #print "%s" %(self.get_element_type(self.element_type)),dataIn
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, dataIn)
                    n += ntotal
            else:
                raise NotImplementedError('num_wide = %s' % (self.num_wide))

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermalLoadVectors'
            #storage_obj = self.thermalLoadVectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)
        return n
