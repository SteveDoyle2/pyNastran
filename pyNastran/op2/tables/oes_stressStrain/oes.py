#pylint: disable=C0301,W0201,R0911,R0915,R0914
"""
Defines the Real/Complex Stresses/Strains created by:
    STRESS = ALL
    STRAIN = ALL
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import b
from six.moves import range
from struct import Struct
from numpy import fromstring, radians, sin, cos, vstack, repeat, array
import numpy as np

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import (RealBar10NodesStressArray, RealBar10NodesStrainArray)

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (RealBeamStressArray, RealBeamStrainArray,
                                                                  RealNonlinearBeamStressArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_compositePlates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import (RealPlateStressArray, RealPlateStrainArray,
                                                                   RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray)
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (RealSpringStressArray, RealSpringStrainArray,
                                                                    NonlinearSpringStress) # TODO: vectorize 1
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray


from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import ComplexBushStress, ComplexBushStrain    # TODO: vectorize 2
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexBush1DStress                   # TODO: vectorize 1
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray, ComplexPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStress, ComplexShearStrain    # TODO: vectorize 2
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexCelasStress, ComplexCelasStrain   # TODO: vectorize 2

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import (NonlinearRod, NonlinearQuad, HyperelasticQuad,  # TODO: vectorize 3
                                                                 RealNonlinearPlateArray)
ComplexShearStressArray = None
ComplexShearStrainArray = None


class OES(OP2Common):
    """
    Defines  the OES class that is used to read stress/strain data
    """
    def __init__(self):
        OP2Common.__init__(self)
        self.ntotal = 0

    #def _oes_cleanup():
    def _read_oes1_3(self, data, ndata):
        """
        reads OES1 subtable 3
        """
        self._data_factor = 1
        self.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', 'load_set'
            'format_code', 'num_wide', 's_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)  # 3

        ## element type
        self.element_type = self.add_data_parameter(data, 'element_type', 'i', 3, False)

        ## load set ID
        self.load_set = self.add_data_parameter(data, 'load_set', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## stress/strain codes
        self.s_code = self.add_data_parameter(data, 's_code', 'i', 11, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            #: mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            #: real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #: mode or cycle TODO confused on the type - F1 means float/int???
            self.mode2 = self.add_data_parameter(data, 'mode2', 'i', 7, False)
            self.cycle = self.add_data_parameter(data, 'cycle', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode2', 'cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.freq = self.add_data_parameter(data, 'freq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.dt = self.add_data_parameter(data, 'dt', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['dt'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        self.fix_format_code()
        self._parse_thermal_code()
        try:
            self.element_name = self.element_mapper[self.element_type]
        except KeyError:
            self.log.error(self.code_information())
            raise
        self.data_code['element_name'] = self.element_name
        if self.is_debug_file:
            self.binary_debug.write('  element_name   = %r\n' % self.element_name)
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)

        self._read_title(data)
        if self.element_type not in self.element_mapper:
            return self._not_implemented_or_skip(data, ndata, self.code_information())

        self._write_debug_bits()
        self._parse_stress_code()
        #print('self.nonlinear_factor =', self.nonlinear_factor)
        #assert self.num_wide != 146, self.code_information()

    def _parse_stress_code(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> isMaxShear=True       isVonMises=False
        stress_bits[0] = 1 -> isMaxShear=False      isVonMises=True

        stress_bits[1] = 0 -> is_stress=True        isStrain=False
        stress_bits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stress_bits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stress_bits[4] = 0 -> material coordinate system flag
        """
        bits = [0, 0, 0, 0, 0]

        s_code = self.s_code
        i = 4
        while s_code > 0:
            value = s_code % 2
            s_code = (s_code - value) // 2
            bits[i] = value
            i -= 1
        self.stress_bits = bits
        self.data_code['stress_bits'] = self.stress_bits

    def _read_oes1_4(self, data, ndata):
        """
        Reads the Stress Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert self.is_stress() == True, self.code_information()
        if self.is_sort1():
            n = self._read_oes1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_ostr1_4(self, data, ndata):
        """
        Reads the Strain Table 4
        """
        #assert self.isubtable == -4, self.isubtable
        #if self.is_debug_file:
            #self.binary_debug.write('  element_name = %r\n' % self.element_name)
        #print "element_name =", self.element_name
        assert self.is_strain() == True, self.code_information()
        if self.is_sort1():
            n = self._read_ostr1_4_sort1(data, ndata)
        else:
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def autojit3(func):
        """
        Debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            n = func(self, data)
            return n
        return new_func

    def print_obj_name_on_crash(func):
        """
        Debugging function to print the object name and an needed parameters
        """
        def new_func(self, data):
            """
            The actual function exec'd by the decorated function.
            """
            try:
                n = func(self, data)
            except NameError:
                raise
            except AttributeError:
                raise
            except:
                raise
                print("----------")
                print(self.obj)
                print(self.data_code)
                if self.obj is not None:
                    #from pyNastran.utils import object_attributes
                    #print object_attributes(self.obj)
                    print(self.obj.data_code)
                print("----------")
                raise
            return n
        return new_func

    #@print_obj_name_on_crash
    def _read_oes1_4_sort1(self, data, ndata):
        """
        Reads OES1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert self.is_sort1() == True
        if self.thermal == 0:
            n = self._read_oes_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    #@print_obj_name_on_crash
    def _read_ostr1_4_sort1(self, data, ndata):
        """
        Reads OSTR1 subtable 4
        """
        #if self.num_wide == 146:
            #assert self.num_wide != 146
            #assert ndata != 146, self.code_information()
        assert self.is_sort1() == True
        if self.thermal == 0:
            n = self._read_oes_loads(data, ndata)
        elif self.thermal == 1:
            n = self._read_oes_thermal(data, ndata)
        else:
            msg = 'thermal=%s' % self.thermal
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_oes_thermal(self, data, ndata):
        """
        Reads OES self.thermal=1 tables; uses a hackish method to just skip the table.
        """
        return ndata

    def _read_ostr_thermal(self, data, ndata):
        """
        Reads OSTR self.thermal=1 tables; uses a hackish method to just skip the table.
        """
        return ndata

    def _read_oes_loads(self, data, ndata):
        """
        Reads OES self.thermal=0 stress/strain
        """
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()
        dt = self.nonlinear_factor

        flag = 'element_id'
        if self.is_stress():
            result_name = 'stress'
            stress_name = 'STRESS'
        else:
            result_name = 'strain'
            stress_name = 'STRAIN'

        if self._results.is_not_saved(result_name):
            return ndata

        if self.element_type in [1, 3, 10]:  # rods
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
            if self.is_stress():
                obj_vector_real = RealRodStressArray
                obj_vector_complex = ComplexRodStressArray
                if self.element_type == 1: # CROD
                    result_name = 'crod_stress'
                elif self.element_type == 3:  # CTUBE
                    result_name = 'ctube_stress'
                elif self.element_type == 10:  # CONROD
                    result_name = 'conrod_stress'
                else:
                    msg = self.code_information()
                    return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealRodStrainArray
                obj_vector_complex = ComplexRodStrainArray
                if self.element_type == 1: # CROD
                    result_name = 'crod_strain'
                elif self.element_type == 3:  # CTUBE
                    result_name = 'ctube_strain'
                elif self.element_type == 10:  # CONROD
                    result_name = 'conrod_strain'
                else:
                    msg = self.code_information()
                    return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 5:  # real
                ntotal = 5 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[axial, torsion, SMa, SMt]
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i4f'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
                        obj.add_new_eid(dt, eid, axial, axial_margin, torsion, torsion_margin)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 5: # imag
                ntotal = 20
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(b(self._endian + 'i4f'))

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    if is_magnitude_phase:
                        print('mag_phase')
                        mag = floats[:, [1, 3]]
                        phase = floats[:, [2, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    else:
                        real = floats[:, [1, 3]]
                        imag = floats[:, [2, 4]]
                        real_imag = real + 1.j * imag
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag

                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = s.unpack(edata)
                        (eid_device, axial_real, axial_imag, torsion_real, torsion_imag) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)

                        if is_magnitude_phase:
                            axial = polar_to_real_imag(axial_real, axial_imag)
                            torsion = polar_to_real_imag(torsion_real, torsion_imag)
                        else:
                            axial = complex(axial_real, axial_imag)
                            torsion = complex(torsion_real, torsion_imag)

                        obj.add_sort1(dt, eid, axial, torsion)
                        n += ntotal
            #elif self.format_code in [2, 3] and self.num_wide == 8:  # is this imag ???
                #ntotal = 32
                #s = self.self.struct_i
                #nelements = ndata // ntotal
                #for i in range(nelements):
                    #edata = data[n:n + 4]
                    #eid_device, = s.unpack(edata)
                    #eid = (eid_device - self.device_code) // 10
                    #assert eid > 0, eid
                    #n += ntotal
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 2: # CBEAM
            # 2-CBEAM
            ## TODO: fix method to follow correct pattern...regarding???

            if self.is_stress():
                result_name = 'cbeam_stress'
            else:
                result_name = 'cbeam_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 111:  # real
                ntotal = 444 # 44 + 10*40  (11 nodes)

                if self.is_stress():
                    obj_vector_real = RealBeamStressArray
                else:
                    obj_vector_real = RealBeamStrainArray

                nelements = ndata // ntotal
                nlayers = nelements * 11
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 11
                    return nelements * self.num_wide * 4
                obj = self.obj
                s = self.struct_i
                nnodes = 10  # 11-1
                ntotal = self.num_wide * 4
                n1 = 44
                n2 = 40
                s1 = Struct(b(self._endian + 'ii9f'))
                s2 = Struct(b(self._endian + 'i9f'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+n1]
                    n += n1

                    out = s1.unpack(edata)
                    eid_device = out[0]
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

                    #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
                    obj.add_new_eid(dt, eid, out[1:])

                    for inode in range(nnodes):
                        edata = data[n:n+n2]
                        n += n2
                        out = s2.unpack(edata)
                        # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
                        obj.add(dt, eid, out)
            elif self.format_code in [2, 3] and self.num_wide == 111:  # imag and random?
                if self.read_mode == 1:
                    return ndata

                ntotal = 444 # 44 + 10*40  (11 nodes)
                #if self.is_stress():
                    #self.create_transient_object(self.beamStress, RealBeamStress)
                #else:
                    #self.create_transient_object(self.beamStrain, RealBeamStrain)

                nelements = ndata // ntotal
                s = self.struct_i

                nnodes = 10  # 11-1
                ntotal = self.num_wide * 4
                n1 = 44
                n2 = 40
                s1 = Struct(b(self._endian + 'ii9f'))
                s2 = Struct(b(self._endian + 'i9f'))

                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+n1]
                    n += n1

                    out = s1.unpack(edata)
                    eid_device = out[0]
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

                    #(grid, sd, ercr, exdr, exer, exfr,
                    #           exci, exdi, exei, exfi) = out
                    #obj.add_new_eid(dt, eid, out[1:])

                    for inode in range(nnodes):
                        edata = data[n:n+n2]
                        n += n2
                        out = s2.unpack(edata)
                        #obj.add(dt, eid, out)
                return ndata
            elif self.format_code == 1 and self.num_wide == 67: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 4: # CSHEAR
            # 4-CSHEAR
            if self.is_stress():
                #obj_real = RealShearStress
                obj_complex = ComplexShearStress

                obj_vector_real = RealShearStressArray
                obj_vector_complex = ComplexShearStressArray
                result_name = 'cshear_stress'
            else:
                #obj_real = RealShearStrain
                obj_complex = ComplexShearStrain

                obj_vector_real = RealShearStrainArray
                obj_vector_complex = ComplexShearStrainArray
                result_name = 'cshear_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 4:  # real
                ntotal = 16  # 4*4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                assert obj is not None
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 4)
                    itime = obj.itime
                    obj._times[itime] = dt
                    if itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 4)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #[max_strain, avg_strain, margin]
                    obj.data[itime, itotal:itotal2, :] = floats[:, 1:]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i3f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = s.unpack(edata)  # num_wide=5
                        if self.is_debug_file:
                            self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

                        (eid_device, max_strain, avg_strain, margin) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        obj.add_new_eid(dt, eid, max_strain, avg_strain, margin)
                        n += ntotal

            elif self.format_code in [2, 3] and self.num_wide == 5:  # imag
                ntotal = 20  # 4*5
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object3(nelements,
                                                       result_name, slot,
                                                       obj_complex, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, etmaxr, etmaxi, etavgr, etavgi)
                    if is_magnitude_phase:
                        real = floats[:, [1, 3]]
                        imag = floats[:, [2, 4]]
                        real_imag = real + 1.j * imag
                    else:
                        mag = floats[:, [1, 3]]
                        phase = floats[:, [2, 4]]
                        rtheta = radians(phase)
                        real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
                    obj.data[obj.itime, itotal:itotal2, :] = real_imag
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'i4f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = s.unpack(edata)  # num_wide=5
                        if self.is_debug_file:
                            self.binary_debug.write('CSHEAR-4 - %s\n' % str(out))
                        (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)

                        if is_magnitude_phase:
                            etmax = polar_to_real_imag(etmaxr, etmaxi)
                            etavg = polar_to_real_imag(etavgr, etavgi)
                        else:
                            etmax = complex(etmaxr, etmaxi)
                            etavg = complex(etavgr, etavgi)
                        obj.add_new_eid_sort1(dt, eid, (etmax, etavg))
                        n += ntotal
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [11, 12, 13, 14]:  # springs
            # 11-CELAS1
            # 12-CELAS2
            # 13-CELAS3
            # 14-CELAS4
            if self.is_stress():
                obj_real = RealSpringStressArray
                if self.element_type == 11:
                    result_name = 'celas1_stress'
                elif self.element_type == 12:
                    result_name = 'celas2_stress'
                elif self.element_type == 13:
                    result_name = 'celas3_stress'
                elif self.element_type == 14:
                    result_name = 'celas4_stress'
                else:
                    raise RuntimeError(self.element_type)
            else:
                obj_real = RealSpringStrainArray
                if self.element_type == 11:
                    result_name = 'celas1_strain'
                elif self.element_type == 12:
                    result_name = 'celas2_strain'
                elif self.element_type == 13:
                    result_name = 'celas3_strain'
                elif self.element_type == 14:
                    result_name = 'celas4_strain'
                else:
                    raise RuntimeError(self.element_type)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 2:  # real
                ntotal = 8 # 2 * 4
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    ielement2 = obj.itotal + nelements
                    itotal2 = ielement2

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 2)
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 2)
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()
                        obj.element[itotal:itotal2] = eids

                    #(eid_device, stress)
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1]
                    obj.itotal = itotal2
                    obj.ielement = ielement2
                else:
                    s = Struct(b(self._endian + 'if'))
                    for i in range(nelements):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (eid_device, ox) = out
                        #eid = (eid_device - self.device_code) // 10
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i result%i=[%i, %f]\n' % (eid, i, eid_device, ox))
                        obj.add_new_eid(dt, eid, ox)
                        n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 3:  # imag
                if self.read_mode == 1:
                    return ndata

                if self.is_stress():
                    # obj_vector_complex = ComplexCelasStressArray
                    self.create_transient_object(slot, ComplexCelasStress)
                else:
                    # obj_vector_complex = ComplexCelasStrainArray
                    self.create_transient_object(slot, ComplexCelasStrain)

                ntotal = 12
                s = Struct(b(self._endian + 'i2f'))
                nelements = ndata // ntotal
                obj = self.obj

                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    (eid_device, axial_real, axial_imag) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    if is_magnitude_phase:
                        axial = polar_to_real_imag(axial_real, axial_imag)
                    else:
                        axial = complex(axial_real, axial_imag)

                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i result%i=[%i, %f, %f]\n' % (eid, i, eid_device, axial_real, axial_imag))
                    obj.add_new_eid_sort1(dt, eid, axial)
                    n += ntotal
            elif self.format_code == 1 and self.num_wide == 3: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 34: # CBAR
            if self.is_stress():
                result_name = 'cbar_stress'
            else:
                result_name = 'cbar_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            slot_vector = slot

            if self.format_code == 1 and self.num_wide == 16:  # real
                if self.is_stress():
                    obj_vector_real = RealBarStressArray
                else:
                    obj_vector_real = RealBarStrainArray

                ntotal = 16 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                    self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                s = Struct(b(self._endian + 'i15f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device,
                     s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                     s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
                    n += ntotal
                    obj.add_new_eid(self.element_name, dt, eid,
                                         s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                                         s1b, s2b, s3b, s4b, smaxb, sminb, MSc)
            elif self.format_code in [2, 3] and self.num_wide == 19:  # imag
                if self.is_stress():
                    obj_vector_complex = ComplexBarStressArray
                else:
                    obj_vector_complex = ComplexBarStrainArray

                ntotal = 76
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
                    self.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                s = Struct(b(self._endian + 'i18f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    n += ntotal
                    out = s.unpack(edata)
                    (eid_device,
                     s1ar, s2ar, s3ar, s4ar, axialr,
                     s1ai, s2ai, s3ai, s4ai, axiali,
                     s1br, s2br, s3br, s4br,
                     s1bi, s2bi, s3bi, s4bi) = out

                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
                    if is_magnitude_phase:
                        s1a = polar_to_real_imag(s1ar, s1ai)
                        s1b = polar_to_real_imag(s1br, s1bi)
                        s2a = polar_to_real_imag(s2ar, s2ai)
                        s2b = polar_to_real_imag(s2br, s2bi)
                        s3a = polar_to_real_imag(s3ar, s3ai)
                        s3b = polar_to_real_imag(s3br, s3bi)
                        s4a = polar_to_real_imag(s4ar, s4ai)
                        s4b = polar_to_real_imag(s4br, s4bi)
                        axial = polar_to_real_imag(axialr, axiali)
                    else:
                        s1a = complex(s1ar, s1ai)
                        s1b = complex(s1br, s1bi)
                        s2a = complex(s2ar, s2ai)
                        s2b = complex(s2br, s2bi)
                        s3a = complex(s3ar, s3ai)
                        s3b = complex(s3br, s3bi)
                        s4a = complex(s4ar, s4ai)
                        s4b = complex(s4br, s4bi)
                        axial = complex(axialr, axiali)

                    obj.add_new_eid_sort1('CBAR', dt, eid,
                                               s1a, s2a, s3a, s4a, axial,
                                               s1b, s2b, s3b, s4b)
            elif self.format_code == 1 and self.num_wide == 19: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [39, 67, 68]: # solid stress
            # 39-CTETRA
            # 67-CHEXA
            # 68-CPENTA
            if self.is_stress():
                obj_vector_real = RealSolidStressArray
                obj_vector_complex = ComplexSolidStressArray
                if self.element_type == 39: # CTETRA
                    nnodes_expected = 5  # 1 centroid + 4 corner points
                    result_name = 'ctetra_stress'
                    element_name = 'CTETRA4'
                elif self.element_type == 67:  # CHEXA
                    nnodes_expected = 9
                    result_name = 'chexa_stress'
                    element_name = 'CHEXA8'
                elif self.element_type == 68:  # CPENTA
                    nnodes_expected = 7
                    result_name = 'cpenta_stress'
                    element_name = 'CPENTA6'
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealSolidStrainArray
                obj_vector_complex = ComplexSolidStrainArray

                if self.element_type == 39: # CTETRA
                    nnodes_expected = 5  # 1 centroid + 4 corner points
                    result_name = 'ctetra_strain'
                    element_name = 'CTETRA4'
                elif self.element_type == 67:  # CHEXA
                    nnodes_expected = 9
                    result_name = 'chexa_strain'
                    element_name = 'CHEXA8'
                elif self.element_type == 68:  # CPENTA
                    nnodes_expected = 7
                    result_name = 'cpenta_strain'
                    element_name = 'CPENTA6'
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            numwide_real = 4 + 21 * nnodes_expected
            numwide_imag = 4 + (17 - 4) * nnodes_expected
            numwide_random = 4 + (11 - 4) * nnodes_expected
            preline1 = '%s-%s' % (self.element_name, self.element_type)
            preline2 = ' ' * len(preline1)

            self._data_factor = nnodes_expected
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 16 + 84 * nnodes_expected
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    n = nelements * 4 * self.num_wide
                    itotal = obj.ielement
                    itotali = obj.itotal + nelements
                    itotalf = obj.itotal + nelements * nnodes_expected
                    obj._times[obj.itime] = dt
                    if obj.itime == 0:
                        # (eid_device, cid, abcd, nnodes)
                        ints = fromstring(data, dtype=self.idtype)
                        try:
                            ints1 = ints.reshape(nelements, numwide_real)
                        except ValueError:
                            msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                            msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (nelements, numwide_real,
                                                                                          nelements * numwide_real)
                            raise ValueError(msg)
                        eids = ints1[:, 0] // 10
                        cids = ints1[:, 1]
                        #nids = ints1[:, 4]
                        assert eids.min() > 0, eids.min()
                        obj.element_node[itotal:itotalf, 0] = repeat(eids, nnodes_expected)
                        ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                        #print('ints2 =', ints2)
                        grid_device = ints2[:, 0]#.reshape(nelements, nnodes_expected)

                        #print('%s-grid_device=%s' % (self.element_name, grid_device))
                        grid_device2 = repeat(grid_device, nnodes_expected)
                        try:
                            obj.element_node[itotal:itotalf, 1] = grid_device
                        except ValueError:
                            msg = '%s; nnodes=%s\n' % (self.element_name, nnodes_expected)
                            msg += 'itotal=%s itotalf=%s\n' % (itotal, itotalf)
                            msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                            #msg += 'nids=%s' % nids
                            raise ValueError(msg)
                        obj.element_cid[itotal:itotali, 0] = eids
                        obj.element_cid[itotal:itotali, 1] = cids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)[:, 4:]
                    # 1     9    15   2    10   16  3   11  17   8
                    #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
                    isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
                    #(grid_device,
                        #sxx, sxy, s1, a1, a2, a3, pressure, svm,
                        #syy, syz, s2, b1, b2, b3,
                        #szz, sxz, s3, c1, c2, c3)
                    floats1 = floats.reshape(nelements * nnodes_expected, 21)#[:, 1:] # drop grid_device

                    obj.data[obj.itime, itotal:itotalf, :] = floats1[:, isave]
                    obj.itotal = itotalf
                    obj.ielement = itotali
                else:
                    struct1 = Struct(b(self._endian + 'ii4si'))
                    struct2 = Struct(b(self._endian + 'i20f'))
                    if self.is_debug_file:
                        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
                            self.element_name, self.element_type, nelements, nnodes_expected)
                        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
                        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
                        self.binary_debug.write(msg)

                    for i in range(nelements):
                        edata = data[n:n+16]
                        out = struct1.unpack(edata)
                        (eid_device, cid, abcd, nnodes) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)

                        if self.is_debug_file:
                            self.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

                        assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

                        n += 16
                        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                            out = struct2.unpack(data[n:n + 84]) # 4*21 = 84
                            if self.is_debug_file:
                                self.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                            (grid_device,
                             sxx, sxy, s1, a1, a2, a3, pressure, svm,
                             syy, syz, s2, b1, b2, b3,
                             szz, sxz, s3, c1, c2, c3) = out

                            if self.is_debug_file:
                                self.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (eid, grid_device, ', '.join(['%r' % di for di in out])))

                            #if grid_device == 0:
                                #grid = 'CENTER'
                            #else:
                                ##grid = (grid_device - device_code) // 10
                                #grid = grid_device

                            grid = grid_device
                            aCos = [a1, a2, a3]
                            bCos = [b1, b2, b3]
                            cCos = [c1, c2, c3]
                            if inode == 0:
                                #element_name = self.element_name + str(nnodes) #  this is correct, but fails
                                obj.add_eid(element_name, cid, dt, eid, grid,
                                            sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                            aCos, bCos, cCos, pressure, svm)
                            else:
                                obj.add_node(dt, eid, inode, grid,
                                             sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                             aCos, bCos, cCos, pressure, svm)
                            n += 84

            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # complex
                ntotal = numwide_imag * 4
                nelements = ndata // ntotal
                self.ntotal += nelements * nnodes_expected
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj

                s1 = Struct(b(self._endian + '2i4si'))
                s2 = Struct(b(self._endian + 'i12f'))
                for i in range(nelements):
                    edata = data[n:n+16]
                    n += 16
                    out = s1.unpack(edata)
                    (eid_device, cid, ctype, nodef) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                    #element_name = self.element_name + str(nodef)  # this is correct, but has problems...
                    obj.add_eid_sort1(self.element_type, element_name, dt, eid, cid, ctype, nodef)
                    for inode in range(nnodes_expected):
                        edata = data[n:n+52]
                        n += 52
                        out = s2.unpack(edata)
                        (grid,
                         exr, eyr, ezr, etxyr, etyzr, etzxr,
                         exi, eyi, ezi, etxyi, etyzi, etzxi) = out
                        #if grid == 0:
                            #grid = 'CENTER'

                        if is_magnitude_phase:
                            ex = polar_to_real_imag(exr, exi)
                            ey = polar_to_real_imag(eyr, eyi)
                            ez = polar_to_real_imag(ezr, ezi)
                            etxy = polar_to_real_imag(etxyr, etxyi)
                            etyz = polar_to_real_imag(etyzr, etyzi)
                            etzx = polar_to_real_imag(etzxr, etzxi)
                        else:
                            ex = complex(exr, exi)
                            ey = complex(eyr, eyi)
                            ez = complex(ezr, ezi)
                            etxy = complex(etxyr, etxyi)
                            etyz = complex(etyzr, etyzi)
                            etzx = complex(etzxr, etzxi)

                        if self.is_debug_file:
                            self.binary_debug.write('       node%s=[%s]\n' % (grid, ', '.join(['%r' % di for di in out])))
                        obj.add_node_sort1(dt, eid, grid, inode,
                                                ex, ey, ez, etxy, etyz, etzx)
            elif self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        #=========================
        # plates
        elif self.element_type == 33:  # CQUAD4-centroidal
            # 33-QUAD4-centroidal
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                result_name = 'cquad4_stress'
            else:
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray
                result_name = 'cquad4_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 17:  # real
                ntotal = 68  # 4*17
                nelements = ndata // ntotal
                nlayers = nelements * 2

                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                s = Struct(b(self._endian + 'i16f'))
                cen = 0 # CEN/4
                obj = self.obj
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)

                    (eid_device,
                     fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
                     fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                    obj.add_new_eid('CQUAD4', dt, eid, cen, fd1, sx1, sy1,
                                         txy1, angle1, major1, minor1, max_shear1)
                    obj.add(dt, eid, cen, fd2, sx2, sy2, txy2,
                                 angle2, major2, minor2, max_shear2)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 15:  # imag
                nnodes = 0  # centroid + 4 corner points
                ntotal = 4 * (15 * (nnodes + 1))
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                obj = self.obj

                s1 = Struct(b(self._endian + 'i14f'))
                s2 = Struct(b(self._endian + 'i14f'))

                cen = 0 # 'CEN/4'
                for i in range(nelements):
                    edata = data[n:n+60]  # 4*15=60
                    n += 60
                    out = s1.unpack(edata)  # 15
                    (eid_device,
                     fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                     fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i C=%s\n' % (eid, str(out)))

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sx2 = complex(sx2r, sx2i)
                        sy1 = complex(sy1r, sy1i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)

                    obj.add_new_eid_sort1('CQUAD4', dt, eid, cen, fd1, sx1, sy1, txy1)
                    obj.add_sort1(dt, eid, cen, fd2, sx2, sy2, txy2)

                    for node_id in range(nnodes):  # nodes pts
                        edata = data[n:n+60]  # 4*15=60
                        n += 60
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('  %s\n' % str(out))
                        (grid,
                         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                        if is_magnitude_phase:
                            sx1 = polar_to_real_imag(sx1r, sx1i)
                            sx2 = polar_to_real_imag(sx2r, sx2i)
                            sy1 = polar_to_real_imag(sy1r, sy1i)
                            sy2 = polar_to_real_imag(sy2r, sy2i)
                            txy1 = polar_to_real_imag(txy1r, txy1i)
                            txy2 = polar_to_real_imag(txy2r, txy2i)
                        else:
                            sx1 = complex(sx1r, sx1i)
                            sx2 = complex(sx2r, sx2i)
                            sy1 = complex(sy1r, sy1i)
                            sy2 = complex(sy2r, sy2i)
                            txy1 = complex(txy1r, txy1i)
                            txy2 = complex(txy2r, txy2i)
                        obj.add_new_node_sort1(dt, eid, grid, fd1, sx1, sy1, txy1)
                        obj.add_sort1(dt, eid, grid, fd2, sx2, sy2, txy2)
            elif self.format_code == 1 and self.num_wide == 0: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type == 74:  # TRIA3
            if self.is_stress():
                #obj_real = RealPlateStress
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                result_name = 'ctria3_stress'
            else:
                #obj_real = RealPlateStrain
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray
                result_name = 'ctria3_strain'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            slot = getattr(self, result_name)
            if self.format_code == 1 and self.num_wide == 17:  # real
                ntotal = 68  # 4*17
                nelements = ndata // ntotal
                nlayers = nelements * 2  # 2 layers per node
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 2
                    return nelements * ntotal

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                    self.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                obj = self.obj
                if self.use_vector and is_vectorized:
                    nfields = 17 * nelements
                    nbytes = nfields * 4
                    istart = obj.itotal
                    iend = obj.itotal + nlayers

                    itime = obj.itime
                    if itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 17)
                        eids = ints[:, 0] // 10
                        ilayers = ints[:, 1]
                        ints2 = ints[:, 1:].reshape(nlayers, 8)
                        assert eids.min() > 0, eids
                        obj._times[obj.itime] = dt
                        obj.element_node[istart:iend:2, 0] = eids
                        obj.element_node[istart+1:iend+1:2, 0] = eids
                        #obj.element_node[istart:iend, 1] = 0
                        #print('obj.element_node\n', obj.element_node)

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 17)
                    floats1 = floats[:, 1:].reshape(nlayers, 8)
                    obj.data[obj.itime, istart:iend, :] = floats1
                    obj.itotal += nlayers
                    n = nbytes
                else:
                #if 1:
                    cen = 0 # 'CEN/3'
                    s = Struct(b(self._endian + 'i16f'))
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = s.unpack(edata)
                        (eid_device,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        if self.is_debug_file:
                            self.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        obj.add_new_eid('CTRIA3', dt, eid, cen, fd1, sx1, sy1,
                                             txy1, angle1, major1, minor1, vm1)
                        obj.add(dt, eid, cen, fd2, sx2, sy2, txy2,
                                     angle2, major2, minor2, vm2)
                        n += ntotal
                #ass
            elif self.format_code in [2, 3] and self.num_wide == 15:  # imag
                ntotal = 60  # 4*15
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj

                s = Struct(b(self._endian + 'i14f'))
                cen = 0 # CEN/3
                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    (eid_device,
                     fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                     fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    if self.is_debug_file:
                        self.binary_debug.write('  OESC CTRIA3-74 - eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sy1 = complex(sy1r, sy1i)
                        sx2 = complex(sx2r, sx2i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)
                    obj.add_new_eid_sort1('CTRIA3', dt, eid, cen, fd1, sx1, sy1, txy1)
                    obj.add_sort1(dt, eid, cen, fd2, sx2, sy2, txy2)
                    n += ntotal
            elif self.format_code == 1 and self.num_wide == 9: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [64, 70, 75, 82, 144]:  # bilinear plates
            # 64-CQUAD8
            # 70-CTRIAR
            # 75-CTRIA6
            # 82-CQUADR
            # 144-CQUAD4-bilinear
            if self.is_stress():
                obj_vector_real = RealPlateStressArray
                obj_vector_complex = ComplexPlateStressArray
                # 64-CQUAD8
                # 70-CTRIAR
                # 75-CTRIA6
                # 82-CQUADR
                # 144-CQUAD4-bilinear
                if self.element_type == 64: # CQUAD8
                    result_name = 'cquad8_stress'
                    #gridC = 'CEN/8'
                elif self.element_type == 70:  # CTRIAR
                    result_name = 'ctriar_stress'
                    #gridC = 'CEN/3'
                elif self.element_type == 75:  # CTRIA6
                    result_name = 'ctria6_stress'
                    #gridC = 'CEN/6'
                elif self.element_type == 82:  # CTRIA6
                    result_name = 'cquadr_stress'
                    #gridC = 'CEN/4'
                elif self.element_type == 144:  # CQUAD4-bilinear
                    # there's no nead to separate this with centroidal strain
                    # because you can only have one in a given OP2
                    result_name = 'cquad4_stress'
                    #gridC = 'CEN/4'
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    return self._not_implemented_or_skip(data, ndata, msg)
            else:
                obj_vector_real = RealPlateStrainArray
                obj_vector_complex = ComplexPlateStrainArray

                if self.element_type == 64: # CQUAD8
                    result_name = 'cquad8_strain'
                    #gridC = 'CEN/8'
                elif self.element_type == 70:  # CTRIAR
                    result_name = 'ctriar_strain'
                    #gridC = 'CEN/3'
                elif self.element_type == 75:  # CTRIA6
                    result_name = 'ctria6_strain'
                    #gridC = 'CEN/6'
                elif self.element_type == 82: # CQUADR
                    result_name = 'cquadr_strain'
                    #gridC = 'CEN/4'
                elif self.element_type == 144: # CQUAD4-bilinear
                    # there's no nead to separate this with centroidal strain
                    # because you can only have one in a given OP2
                    result_name = 'cquad4_strain'
                    #gridC = 'CEN/4'
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)

            if self.element_type in [64, 82, 144]:
                nnodes = 4 # + 1 centroid
            elif self.element_type in [70, 75]:
                nnodes = 3 # + 1 centroid
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

            gridC = 0
            slot = getattr(self, result_name)
            numwide_real = 2 + 17 * (nnodes + 1)
            numwide_imag = 2 + 15 * (nnodes + 1)
            numwide_random = 2 + 9 * (nnodes + 1)

            etype = self.element_name
            #gridC = 'CEN/%i' % nnodes
            if self.format_code == 1 and self.num_wide == numwide_real:  # real
                ntotal = 4 * (2 + 17 * (nnodes + 1))
                nelements = ndata // ntotal
                nlayers = 2 * nelements * (nnodes + 1)  # 2 layers per node
                nnodes_all = nnodes + 1

                self._data_factor = 10  # TODO: why is this 10?
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nlayers
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                        ints1 = ints[:, 2:].reshape(nlayers//2, 17)[:, 0].reshape(nelements, nnodes_all)
                        ints1[:, 0] = 0.
                        nids = ints1.ravel()

                        eids = ints[:, 0] // 10
                        eids2 = array([eids] * (nnodes_all * 2), dtype='int32').T.ravel()
                        #nids2 = array([nids, nids], dtype='int32').T.ravel()
                        #eids2 = vstack([eids] * (nnodes_all * 2)).T.ravel()
                        nids2 = vstack([nids, nids]).T.ravel()
                        #eids2 = repeat(eids, 2 * nnodes_all)
                        # nids2 = repeat(nids, 2)
                        obj.element_node[istart:iend, 0] = eids2
                        obj.element_node[istart:iend, 1] = nids2
                        #assert obj.element_node[:iend, 0].min() > 0, eids2
                        if obj.nonlinear_factor is not None:
                            float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                            float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                            obj.float_mask = float_mask1

                    if obj.nonlinear_factor is not None:
                        results = fromstring(data, dtype=self.fdtype)[obj.float_mask]
                    else:
                        floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                        floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                        results = floats1[:, 1:].reshape(nlayers, 8)

                    #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
                    obj.data[obj.itime, istart:iend, :] = results
                #if 1:
                else:
                    n = 0
                    center_format = b'i4si16f'
                    node_format = b'i16f'
                    cs = Struct(center_format)
                    ns = Struct(node_format)

                    if self.is_debug_file:
                        self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                        self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                        self.binary_debug.write('  #elementi = [centeri, node1i, node2i, node3i, node4i]\n')
                        self.binary_debug.write('  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                        self.binary_debug.write('  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                        self.binary_debug.write('  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
                        self.binary_debug.write('  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n')
                        self.binary_debug.write('  nelements=%i; nnodes=%i # +1 centroid\n' % (nelements, nnodes))

                    for i in range(nelements):
                        edata = data[n:n+76]

                        out = cs.unpack(edata)  # len=17*4
                        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
                        (eid_device, j,
                         grid,
                         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        #print(out[:3])
                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

                        obj.add_new_eid(etype, dt, eid, gridC, fd1, sx1, sy1,
                                             txy1, angle1, major1, minor1, vm1)
                        obj.add(dt, eid, gridC, fd2, sx2, sy2, txy2,
                                     angle2, major2, minor2, vm2)
                        n += 76
                        for inode in range(nnodes):
                            out = ns.unpack(data[n:n + 68])
                            (grid,
                             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

                            if self.is_debug_file:
                                d = tuple([grid,
                                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                                self.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
                            assert isinstance(grid, int), grid
                            assert grid > 0, grid

                            obj.addNewNode(dt, eid, grid, fd1, sx1, sy1,
                                                txy1, angle1, major1, minor1, vm1)
                            obj.add(dt, eid, grid, fd2, sx2, sy2,
                                         txy2, angle2, major2, minor2, vm2)
                            n += 68
            elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                ntotal = numwide_imag * 4
                assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide*4, ntotal)
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * ntotal
                obj = self.obj

                s1 = self.struct_2i  # 2
                s2 = Struct(b(self._endian + 'i14f')) # 15
                for i in range(nelements):
                    (eid_device, _) = s1.unpack(data[n:n+8])
                    n += 8

                    edata = data[n:n+60]  # 4*15
                    n += 60
                    out = s2.unpack(edata)  # len=15*4

                    eid = self._check_id(eid_device, 'element', stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('%s\n' % (str(out)))
                    (grid,
                     fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                     fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
                    #gridC = 'CEN/%i' % grid   # this is correct, but fails

                    if is_magnitude_phase:
                        sx1 = polar_to_real_imag(sx1r, sx1i)
                        sy1 = polar_to_real_imag(sy1r, sy1i)
                        sx2 = polar_to_real_imag(sx2r, sx2i)
                        sy2 = polar_to_real_imag(sy2r, sy2i)
                        txy1 = polar_to_real_imag(txy1r, txy1i)
                        txy2 = polar_to_real_imag(txy2r, txy2i)
                    else:
                        sx1 = complex(sx1r, sx1i)
                        sy1 = complex(sy1r, sy1i)
                        sx2 = complex(sx2r, sx2i)
                        sy2 = complex(sy2r, sy2i)
                        txy1 = complex(txy1r, txy1i)
                        txy2 = complex(txy2r, txy2i)

                    obj.add_new_eid_sort1(etype, dt, eid, gridC, fd1, sx1, sy1, txy1)
                    obj.add_sort1(dt, eid, gridC, fd2, sx2, sy2, txy2)

                    for node_id in range(nnodes):  # nodes pts
                        edata = data[n:n + 60]  # 4*15=60
                        n += 60
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s\n' % (str(out)))
                        (grid,
                         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
                         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

                        if is_magnitude_phase:
                            sx1 = polar_to_real_imag(sx1r, sx1i)
                            sx2 = polar_to_real_imag(sx2r, sx2i)
                            sy1 = polar_to_real_imag(sy1r, sy1i)
                            sy2 = polar_to_real_imag(sy2r, sy2i)
                            txy1 = polar_to_real_imag(txy1r, txy1i)
                            txy2 = polar_to_real_imag(txy2r, txy2i)
                        else:
                            sx1 = complex(sx1r, sx1i)
                            sx2 = complex(sx2r, sx2i)
                            sy1 = complex(sy1r, sy1i)
                            sy2 = complex(sy2r, sy2i)
                            txy1 = complex(txy1r, txy1i)
                            txy2 = complex(txy2r, txy2i)

                        obj.add_new_node_sort1(dt, eid, grid, fd1, sx1, sy1, txy1)
                        obj.add_sort1(dt, eid, grid, fd2, sx2, sy2, txy2)
            elif self.format_code == 1 and self.num_wide == numwide_random: # random
                msg = 'numwide_real=%s\n' % numwide_real
                msg = 'numwide_imag=%s\n' % numwide_imag
                msg = 'numwide_random=%s\n' % numwide_random
                msg += self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [88, 90]: # nonlinear shells
            # 88-CTRIA3NL
            # 90-CQUAD4NL
            obj_vector_real = RealNonlinearPlateArray
            if self.is_stress():
                if self.element_type == 88:
                    result_name = 'nonlinear_ctria3_stress'
                elif self.element_type == 90:
                    result_name = 'nonlinear_cquad4_stress'
                else:
                    raise RuntimeError(self.element_type)
            else:
                if self.element_type == 88:
                    result_name = 'nonlinear_ctria3_strain'
                elif self.element_type == 90:
                    result_name = 'nonlinear_cquad4_strain'
                else:
                    raise RuntimeError(self.element_type)

            slot = getattr(self, result_name)
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == 13 and self.element_type in [88, 90]:  # real
                # single layered hyperelastic (???) ctria3, cquad4
                ntotal = 52  # 4*13
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.use_vector and is_vectorized:
                # self.itime = 0
                # self.ielement = 0
                # self.itotal = 0
                #self.ntimes = 0
                #self.nelements = 0
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints1[:, 0] // 10
                    obj.element[istart:iend, 0] = eids

                results = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)

                #[fiber_distance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
                obj.data[obj.itime, istart:iend, :] = results[:, 1:]
                #else:
                if 0:
                    if self.is_stress():
                        self.create_transient_object(slot, NonlinearQuad)
                    else:
                        self.create_transient_object(slot, NonlinearQuad)

                    s = Struct(b(self._endian + 'i12f'))  # 1+12=13
                    for i in range(nelements):
                        edata = data[n:n + ntotal]
                        out = s.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))

                        (eid_device, fd1,
                         sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                         ex1, ey1, ez1, exy1) = out
                        eid = self._check_id(eid_device, flag, stress_name, out)
                        indata = (eid, fd1,
                                  sx1, sy1, sz1, txy1, es1, eps1, ecs1,
                                  ex1, ey1, ez1, exy1)
                        obj.add_new_eid(self.element_type, dt, indata)
                        #print("eid=%s axial=%s equivStress=%s totalStrain=%s effPlasticCreepStrain=%s effCreepStrain=%s linearTorsionalStresss=%s" % (
                            #eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss))
                    n += ntotal
            elif self.format_code == 1 and self.num_wide == 25 and self.element_type in [88, 90]:  # TODO: real?
                if self.read_mode == 1:
                    return ndata
                #cfjhfg
                if self.is_stress():
                    self.create_transient_object(slot, NonlinearQuad)
                else:
                    self.create_transient_object(slot, NonlinearQuad)

                ntotal = 100  # 4*25
                s = Struct(b(self._endian + 'i24f')) # 1+24=25
                nelements = ndata // ntotal
                obj = self.obj
                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('CQUADNL-90 - %s\n' % str(out))
                    (eid_device,
                     fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1,
                     fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    #in_data = (eid, fd1, sx1, sy1, txy1, es1, eps1, ecs1, ex1, ey1, etxy1,
                                    #fd2, sx2, sy2, txy2, es2, eps2, ecs2, ex2, ey2, etxy2)
                    obj.add_new_eid(self.element_type, dt, (
                        eid, fd1, sx1, sy1, undef1, txy1, es1, eps1, ecs1, ex1, ey1, undef2, etxy1))
                    obj.add(dt, (
                        eid, fd2, sx2, sy2, undef3, txy2, es2, eps2, ecs2, ex2, ey2, undef4, etxy2))
                    n += ntotal
            elif self.format_code == 1 and self.num_wide == 0: # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                raise RuntimeError(self.code_information())

        elif self.element_type in [95, 96, 97, 98]: # composite shell
            # 95 - CQUAD4
            # 96 - CQUAD8
            # 97 - CTRIA3
            # 98 - CTRIA6 (composite)
            if self.is_stress():
                # obj_complex = ComplexCompositePlateStress
                ComplexCompositePlateStressArray = None
                obj_vector_real = RealCompositePlateStressArray
                obj_vector_complex = ComplexCompositePlateStressArray
                if self.element_type == 95: # CQUAD4
                    result_name = 'cquad4_composite_stress'
                elif self.element_type == 96:  # CQUAD8
                    result_name = 'cquad8_composite_stress'
                elif self.element_type == 97:  # CTRIA3
                    result_name = 'ctria3_composite_stress'
                elif self.element_type == 98:  # CTRIA6
                    result_name = 'ctria6_composite_stress'
                #elif self.element_type == ???:  # CTRIA6
                    #result_name = 'ctriar_composite_stress'
                #elif self.element_type == 10:  # CTRIA6
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    raise RuntimeError(self.code_information())
                    #return self._not_implemented_or_skip(data, ndata, msg)
            else:
                ComplexCompositePlateStrainArray = None
                obj_vector_real = RealCompositePlateStrainArray
                obj_vector_complex = ComplexCompositePlateStrainArray
                if self.element_type == 95: # CQUAD4
                    result_name = 'cquad4_composite_strain'
                elif self.element_type == 96:  # CQUAD8
                    result_name = 'cquad8_composite_strain'
                elif self.element_type == 97:  # CTRIA3
                    result_name = 'ctria3_composite_strain'
                elif self.element_type == 98:  # CTRIA6
                    result_name = 'ctria6_composite_strain'
                #elif self.element_type == ???:  # CTRIA6
                    #result_name = 'ctriar_composite_strain'
                else:
                    msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                    raise RuntimeError(self.code_information())
                    #return self._not_implemented_or_skip(data, ndata, msg)

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            etype = self.element_name
            if self.format_code == 1 and self.num_wide == 11:  # real
                ntotal = 44
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

                if self.use_vector and is_vectorized: #  and self.element_type in [144]
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nelements * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nelements
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = fromstring(data, dtype=self.idtype).reshape(nelements, 11)
                        eids = ints[:, 0] // 10
                        nids = ints[:, 1]
                        obj.element_layer[istart:iend, 0] = eids
                        obj.element_layer[istart:iend, 1] = nids

                    floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 11)
                    #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 2:]
                else:
                    s = Struct(b(self._endian + 'ii9f')) # 11
                    eid_old = 0
                    if hasattr(self, 'eid_old'):
                        eid_old = self.eid_old

                    for i in range(nelements):
                        edata = data[n:n+44]  # 4*11
                        out = s.unpack(edata)
                        (eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
                        # TODO: this is a hack that I think is composite specific
                        eid = eid_device // 10

                        if self.is_debug_file:
                            self.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))

                        if eid != eid_old:  # originally initialized to None, the buffer doesnt reset it, so it is the old value
                            obj.add_new_eid(etype, dt, eid, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                        else:
                            obj.add(dt, eid, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)
                        eid_old = eid
                        n += 44
                    self.eid_old = eid_old
            elif self.format_code in [2, 3] and self.num_wide == 9:  # TODO: imag? - not done...
                aaa
                msg = self.code_information()
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_complex)
                if auto_return:
                    return nelements * self.num_wide * 4
                #return self._not_implemented_or_skip(data, ndata, msg)

                # TODO: this is an OEF result???
                #    furthermore the actual table is calle dout as
                #    'i8si4f4s', not 'i8si3fi4s'
                ntotal = 36
                nelements = ndata // ntotal
                s = self.struct_i
                s2 = Struct(b(self._endian + '8si3fi4s'))
                s3 = Struct(b(self._endian + '8si4f4s'))
                for i in range(nelements):
                    #out = s.unpack(data[n:n + ntotal])
                    eid_device, = s.unpack(data[n:n+4])
                    #t, = s.unpack(data[n:n+4])
                    #eid = (eid_device - self.device_code) // 10

                    if eid_device > 0:
                        out = s2.unpack(data[n+4:n+ntotal])
                    else:
                        out1 = s2.unpack(data[n+4:n+ntotal])
                        out = s3.unpack(data[n+4:n+ntotal])
                    (theory, lamid, fp, fm, fb, fmax, fflag) = out

                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - (%s) + %s\n' % (self.element_name, self.element_type, eid_device, str(out)))
                    obj.add_new_eid(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
                    n += ntotal
                raise NotImplementedError('this is a really weird case...')
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        #=========================
        elif self.element_type == 53: # axial plates - ctriax6
            # 53 - CTRIAX6
            if self.is_stress():
                result_name = 'ctriax_stress'
            else:
                result_name = 'ctriax_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 33: # real
                if self.is_stress():
                    obj_vector_real = RealTriaxStressArray
                else:
                    obj_vector_real = RealTriaxStrainArray

                ntotal = 132  # (1+8*4)*4 = 33*4 = 132
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 4
                    return nelements * self.num_wide * 4

                obj = self.obj

                s1 = Struct(b(self._endian + '2i7f'))  # 36
                s2 = Struct(b(self._endian + 'i7f'))
                for i in range(nelements):
                    out = s1.unpack(data[n:n + 36])
                    (eid_device, loc, rs, azs, As, ss, maxp, tmax, octs) = out
                    if self.is_debug_file:
                        self.binary_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                    n += 36
                    for i in range(3):
                        out = s2.unpack(data[n:n + 32])
                        (loc, rs, azs, As, ss, maxp, tmax, octs) = out
                        if self.is_debug_file:
                            self.binary_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
                        #print "eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" % (eid,loc,rs,azs,As,ss,maxp,tmax,octs)
                        obj.add_sort1(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
                        n += 32  # 4*8
            elif self.format_code in [2, 3] and self.num_wide == 37: # imag
                if self.read_mode == 1:
                    return ndata
                #msg = 'num_wide=%s' % self.num_wide
                #return self._not_implemented_or_skip(data, ndata, msg)

                #if self.is_stress():
                    #self.create_transient_object(self.ctriax_stress, ComplexTriaxStress)  # undefined
                    #raise NotImplementedError('ComplexTriaxStress')
                #else:
                    #self.create_transient_object(self.ctriax_strain, ComplexTriaxStrain)  # undefined
                    #raise NotImplementedError('ComplexTriaxStrain')
                s1 = Struct(b(self._endian + 'ii8f'))
                s2 = Struct(b(self._endian + 'i8f'))

                num_wide = 1 + 4 * 9
                ntotal = num_wide * 4
                assert num_wide == self.num_wide, num_wide
                nelements = ndata // ntotal  # (1+8*4)*4 = 33*4 = 132

                for i in range(nelements):
                    out = s1.unpack(data[n:n + 40])
                    (eid_device, loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))

                    if is_magnitude_phase:
                        rs = polar_to_real_imag(rsr, rsi)
                        azs = polar_to_real_imag(azsr, azsi)
                        As = polar_to_real_imag(Asr, Asi)
                        ss = polar_to_real_imag(ssr, ssi)
                    else:
                        rs = complex(rsr, rsi)
                        azs = complex(azsr, azsi)
                        As = complex(Asr, Asi)
                        ss = complex(ssr, ssi)
                    #obj.add_new_eid(dt, eid, loc, rs, azs, As, ss)

                    n += 40
                    for i in range(3):
                        out = s2.unpack(data[n:n + 36])
                        (loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
                        if self.is_debug_file:
                            self.binary_debug.write('    %s\n' % (str(out)))
                        #print("eid=%s loc=%s rs=%s azs=%s as=%s ss=%s" % (eid, loc, rs, azs, As, ss))

                        if is_magnitude_phase:
                            rs = polar_to_real_imag(rsr, rsi)
                            azs = polar_to_real_imag(azsr, azsi)
                            As = polar_to_real_imag(Asr, Asi)
                            ss = polar_to_real_imag(ssr, ssi)
                        else:
                            rs = complex(rsr, rsi)
                            azs = complex(azsr, azsi)
                            As = complex(Asr, Asi)
                            ss = complex(ssr, ssi)
                        #obj.add(dt, eid, loc, rs, azs, As, ss)
                        n += 36  # 4*8
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type == 102: # bush
            # 102-CBUSH
            if self.is_stress():
                result_name = 'cbush_stress'
            else:
                result_name = 'cbush_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 7:  # real
                if self.is_stress():
                    obj_vector_real = RealBushStressArray
                else:
                    obj_vector_real = RealBushStrainArray

                assert self.num_wide == 7, "num_wide=%s not 7" % self.num_wide
                ntotal = 28  # 4*7

                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4
                obj = self.obj

                s = Struct(b(self._endian + 'i6f'))
                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=7
                    if self.is_debug_file:
                        self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                    (eid_device, tx, ty, tz, rx, ry, rz) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    obj.add_sort1(dt, eid, tx, ty, tz, rx, ry, rz)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 13:  # imag
                if self.read_mode == 1:
                    return ndata
                if self.is_stress():
                    self.create_transient_object(self.cbush_stress, ComplexBushStress)
                else:
                    self.create_transient_object(self.cbush_strain, ComplexBushStrain)
                ntotal = 52  # 4*13

                nelements = ndata // ntotal
                obj = self.obj

                s = Struct(b(self._endian + 'i12f'))
                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=7
                    if self.is_debug_file:
                        self.binary_debug.write('CBUSH-102 - %s\n' % str(out))

                    (eid_device,
                     txr, tyr, tzr, rxr, ryr, rzr,
                     txi, tyi, tzi, rxi, ryi, rzi) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    if is_magnitude_phase:
                        tx = polar_to_real_imag(txr, txi)
                        ty = polar_to_real_imag(tyr, tyi)
                        tz = polar_to_real_imag(tzr, tzi)
                        rx = polar_to_real_imag(rxr, rxi)
                        ry = polar_to_real_imag(ryr, ryi)
                        rz = polar_to_real_imag(rzr, rzi)
                    else:
                        tx = complex(txr, txi)
                        ty = complex(tyr, tyi)
                        tz = complex(tzr, tzi)
                        rx = complex(rxr, rxi)
                        ry = complex(ryr, ryi)
                        rz = complex(rzr, rzi)
                    obj.add_new_eid(self.element_type, dt, eid, tx, ty, tz, rx, ry, rz)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 40:  # bush
            # 40-CBUSH1D
            if self.is_stress():
                result_name = 'cbush1d_stress_strain'
            else:
                result_name = 'cbush1d_stress_strain'
            self._results._found_result(result_name)
            slot = self.cbush1d_stress_strain

            if self.format_code == 1 and self.num_wide == 8:  # real
                if self.is_stress():
                    obj_vector_real = RealBush1DStressArray
                else:
                    #self.create_transient_object(self.cbush1d_stress_strain, Bush1DStrain)  # undefined
                    raise NotImplementedError('cbush1d_stress_strain; numwide=8')

                ntotal = 32  # 4*8
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(b(self._endian + 'i6fi'))
                nelements = ndata // ntotal
                obj = self.obj

                for i in range(nelements):
                    edata = data[n:n + ntotal]
                    out = s.unpack(edata)  # num_wide=25
                    if self.is_debug_file:
                        self.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
                    (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    # axial_force, axial_displacement, axial_velocity, axial_stress, axial_strain, plastic_strain, is_failed
                    obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
                    n += ntotal
            elif self.format_code in [2, 3] and self.num_wide == 9:  # imag
                if self.read_mode == 1:
                    return ndata
                if self.is_stress():
                    self.create_transient_object(self.cbush1d_stress_strain, ComplexBush1DStress)  # undefined
                else:
                    #self.create_transient_object(self.cbush1d_stress_strain, ComplexBush1DStress)  # undefined
                    raise NotImplementedError('self.cbush1d_stress_strain; complex strain')

                ntotal = 36  # 4*9
                s = Struct(b(self._endian + 'i8f'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)  # num_wide=25
                    (eid_device,
                     fer, uer, aor, aer,
                     fei, uei, aoi, aei) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    if is_magnitude_phase:
                        fe = polar_to_real_imag(fer, fei)
                        ue = polar_to_real_imag(uer, uei)
                        ao = polar_to_real_imag(aor, aoi)
                        ae = polar_to_real_imag(aer, aei)
                    else:
                        fe = complex(fer, fei)
                        ue = complex(uer, uei)
                        ao = complex(aor, aoi)
                        ae = complex(aer, aei)
                    obj.add_new_eid(self.element_type, dt, eid, fe, ue, ao, ae)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)


        elif self.element_type in [87, 89, 92]:  # nonlinear rods
            # 87-CTUBENL
            # 89-RODNL
            # 92-CONRODNL
            if self.read_mode == 1:
                return ndata

            if self.is_stress():
                if self.element_type == 87:
                    result_name = 'nonlinear_ctube_stress'
                    name = 'CTUBENL-87'
                elif self.element_type == 89:
                    result_name = 'nonlinear_crod_stress'
                    name = 'RODNL-89'
                elif self.element_type == 92:
                    result_name = 'nonlinear_conrod_stress'
                    name = 'CONRODNL-92'
                else:
                    return self._not_implemented_or_skip(data, ndata, self.code_information())
            else:
                if self.element_type == 87:
                    result_name = 'nonlinear_ctube_strain'
                    name = 'CTUBENL-87'
                elif self.element_type == 89:
                    result_name = 'nonlinear_crod_strain'
                    name = 'RODNL-89'
                elif self.element_type == 92:
                    result_name = 'nonlinear_conrod_strain'
                    name = 'CONRODNL-92'
                else:
                    return self._not_implemented_or_skip(data, ndata, self.code_information())
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 7:  # real
                self.create_transient_object(slot, NonlinearRod)
                s = Struct(b(self._endian + 'i6f'))  # 1+6=7
                ntotal = 28  #  7*4 = 28
                nelements = ndata // ntotal
                obj = self.obj

                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)

                    (eid_device, axial, equivStress, totalStrain, effPlasticCreepStrain,
                     effCreepStrain, linearTorsionalStresss) = out
                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('%s - %s\n' % (name, str(out)))
                    indata = (eid, axial, equivStress, totalStrain, effPlasticCreepStrain, effCreepStrain, linearTorsionalStresss)
                    obj.add(self.element_type, dt, indata)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [224, 225]: # nonlinear spring
            # 224-CELAS1
            # 225-CELAS3
            # nonlinearSpringStress
            numwide_real = 3
            if self.read_mode == 1:
                return ndata
            if self.is_stress():
                if self.element_type == 224:
                    result_name = 'nonlinear_celas1_stress'
                    slot = self.nonlinear_celas1_stress
                elif self.element_type == 225:
                    result_name = 'nonlinear_celas3_stress'
                    slot = self.nonlinear_celas3_stress
            else:
                raise NotImplementedError('NonlinearSpringStrain')

            self._results._found_result(result_name)
            if self.num_wide == numwide_real:
                if self.is_stress():
                    self.create_transient_object(slot, NonlinearSpringStress)
                else:
                    #self.create_transient_object(self.nonlinearSpringStrain, NonlinearSpringStrain)  # undefined
                    raise NotImplementedError('NonlinearSpringStrain')

                assert self.num_wide == 3, "num_wide=%s not 3" % self.num_wide
                ntotal = 12  # 4*3
                nelements = ndata // ntotal
                obj = self.obj

                s = Struct(b(self._endian + 'i2f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)  # num_wide=3
                    (eid_device, force, stress) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (self.element_name, self.element_type, str(out)))
                    obj.add_new_eid(self.element_name, dt, eid, force, stress)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 35:
            # 35-CON
            return ndata
        elif self.element_type in [60, 61]:
            # 60-DUM8
            # 61-DUM9
            return ndata
        elif self.element_type == 69:  # cbend
            # 69-CBEND
            return ndata
        elif self.element_type == 86:  # cgap
            # 86-GAPNL
            if self.is_stress():
                result_name = 'nonlinear_cgap_stress'
            else:
                result_name = 'nonlinear_cgap_strain'
            self._results._found_result(result_name)
            slot = getattr(self, result_name)
            #print('self.nonlinear_factor = ', self.nonlinear_factor)
            if self.format_code == 1 and self.num_wide == 11:  # real?
                if self.is_stress():
                    obj_vector_real = NonlinearGapStressArray
                else:
                    raise NotImplementedError('NonlinearGapStrain')

                ntotal = 44  # 4*11
                nelements = ndata // ntotal
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                s = Struct(b(self._endian + 'i8f4s4s'))
                nelements = ndata // ntotal
                obj = self.obj

                for i in range(nelements):
                    edata = data[n:n + ntotal]

                    out = s.unpack(edata)  # num_wide=25
                    (eid_device, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('CGAPNL-86 - %s\n' % str(out))
                    obj.add_sort1(dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 94:
            if self.read_mode == 0:
                return ndata
            # 94-BEAMNL
            numwide_real = 51
            numwide_random = 0

            if self.is_stress():
                result_name = 'nonlinear_cbeam_stress'
                slot = self.nonlinear_cbeam_stress
            else:
                result_name = 'nonlinear_cbeam_strain'
                slot = self.nonlinear_cbeam_strain
            self._results._found_result(result_name)

            if self.format_code == 1 and self.num_wide == numwide_real:
                msg = result_name
                if self.is_stress():
                    obj_vector_real = RealNonlinearBeamStressArray
                else:
                    raise NotImplementedError('Nonlinear CBEAM Strain...this should never happen')

                ntotal = numwide_real * 4
                nelements = ndata // ntotal

                nlayers = nelements * 8
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = 8
                    return ndata
                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    #self.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
                    #self.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
                    #self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)


                s = Struct(b(self._endian + '2i 4s5f 4s5f 4s5f 4s5f i 4s5f 4s5f 4s5f 4s5f'))  # 2 + 6*8 + 1 = 51
                for i in range(nelements):  # num_wide=51
                    edata = data[n:n + 204]
                    out = s.unpack(edata)

                    if self.is_debug_file:
                        self.binary_debug.write('BEAMNL-94 - %s\n' % str(out))

                    #gridA, CA, long_CA, eqS_CA, tE_CA, eps_CA, ecs_CA,
                    #       DA, long_DA, eqS_DA, tE_DA, eps_DA, ecs_DA,
                    #       EA, long_EA, eqS_EA, tE_EA, eps_EA, ecs_EA,
                    #       FA, long_FA, eqS_FA, tE_FA, eps_FA, ecs_FA,
                    #gridB, CB, long_CB, eqS_CB, tE_CB, eps_CB, ecs_CB,
                    #       DB, long_DB, eqS_DB, tE_DB, eps_DB, ecs_DB,
                    #       EB, long_EB, eqS_EB, tE_EB, eps_EB, ecs_EB,
                    #       FB, long_FB, eqS_FB, tE_FB, eps_FB, ecs_FB,
                    # A
                    assert out[3-1] == b'   C', out[3-1]
                    assert out[9-1] == b'   D', out[9-1]
                    assert out[15-1] == b'   E', out[15-1]
                    assert out[21-1] == b'   F', out[21-1]

                    # B
                    assert out[28-1] == b'   C', out[28-1]
                    assert out[34-1] == b'   D', out[34-1]
                    assert out[40-1] == b'   E', out[40-1]
                    assert out[46-1] == b'   F', out[46-1]

                    eid_device = out[0]
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    obj.add_new_eid_sort1(dt, eid, out)
                    n += 204

            elif self.format_code == 1 and self.num_wide == numwide_random:  # random
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [85, 91, 93]:
            if self.read_mode == 1:
                return ndata
            # 85-TETRANL
            # 91-PENTANL
            # 93-HEXANL

            #85: 2 + (18 - 2) * 5,  # Nonlinear CTETRA
            #91: 4 + (25 - 4) * 7,  # Nonlinear CPENTA
            #93: 4 + (25 - 4) * 9,  # Nonlinear CHEXA
            if self.element_type == 85:
                etype = 'CTETRANL'
                nnodes = 5
            elif self.element_type == 91:
                etype = 'CPENTANL'
                nnodes = 7
            elif self.element_type == 93:
                etype = 'CHEXANL'
                nnodes = 9
            else:
                return self._not_implemented_or_skip(data, ndata, self.code_information())

            numwide_real = 4 + (25 - 4) * nnodes  # real???
            numwide_random = 2 + (18 - 2) * nnodes  # imag???
            #print("numwide=%s numwide_real=%s numwide_random=%s" % (self.num_wide, numwide_real, numwide_random))

            #numwide_real = 0
            #numwide_imag = 2 + 16 * nnodes
            #ntotal = 8 + 64 * nnodes

            if self.format_code == 1 and self.num_wide == numwide_real:
                ntotal = numwide_real * 4
                #if self.is_stress():
                    #self.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
                #else:
                    #self.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)
                #self.handle_results_buffer(self.OES_CQUAD4NL_90, resultName, name)
                raise RuntimeError('OES_CQUAD4NL_90')
            elif self.format_code == 1 and self.num_wide == numwide_random:  # random
            #elif self.format_code in [2, 3] and self.num_wide == numwide_imag:  # imag
                ntotal = numwide_random * 4
                #if self.is_stress():
                    #self.create_transient_object(self.nonlinearPlateStress, NonlinearSolid)
                #else:
                    #self.create_transient_object(self.nonlinearPlateStrain, NonlinearSolid)

                n = 0
                s1 = Struct(b(self._endian + 'i4s'))
                s2 = Struct(b(self._endian + 'i15f'))
                nelements = ndata // ntotal
                for i in range(nelements):  # 2+16*9 = 146 -> 146*4 = 584
                    edata = data[n:n+8]
                    n += 8

                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))
                    (eid_device, ctype) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    for i in range(nnodes):
                        edata = data[n:n+64]
                        n += 64
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s-%sB - %s\n' % (etype, self.element_type, str(out)))

                        assert len(out) == 16
                        (grid,
                         sx, sy, sz, sxy, syz, sxz, se, eps, ecs,
                         ex, ey, ez, exy, eyz, exz) = out
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 100:  # bars
            # 100-BARS
            if self.read_mode == 0:
                return ndata
            if self.is_stress():
                result_name = 'cbar_stress_10nodes'
            else:
                result_name = 'cbar_strain_10nodes'

            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.format_code == 1 and self.num_wide == 10:  # real
                # if this isn't vectorized, it will be messed up...
                RealBar10NodesStress = RealBar10NodesStressArray
                RealBar10NodesStrain = RealBar10NodesStrainArray
                if self.is_stress():
                    obj_real = RealBar10NodesStress
                    obj_vector_real = RealBar10NodesStressArray
                else:
                    obj_real = RealBar10NodesStrain
                    obj_vector_real = RealBar10NodesStrainArray

                ntotal = 10 * 4
                nelements = ndata // ntotal

                auto_return, is_vectorized = self._create_oes_object3(nelements,
                                                       result_name, slot,
                                                       obj_real, obj_vector_real)
                if auto_return:
                    return ndata

                if self.is_debug_file:
                    self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                    #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                    self.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
                    self.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
                obj = self.obj

                s = Struct(b(self._endian + 'i9f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]
                    out = s.unpack(edata)
                    (eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS) = out
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
                    n += ntotal
                    # continue
                    obj.add_new_eid(self.element_name, dt, eid,
                                         sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)
            else:
                raise NotImplementedError(self.code_information())
        elif self.element_type == 101:
            # 101-AABSF
            return ndata
        elif self.element_type in [140, 201]:
            # 140-HEXA8FD, 201-QUAD4FD
            return ndata
        elif self.element_type in [145, 146, 147]:
            if self.read_mode == 1:
                return ndata
            # 145-VUHEXA  (8 nodes)
            # 146-VUPENTA (6 nodes)
            # 147-VUTETRA (4 nodes)
            if self.element_type == 147:
                etype = 'VUTETRA'
                nnodes = 4
                #numwide_a = 2 + (14 - 2) * nnodes  # 50
                #numwide_b = 2 + (9 - 2) * nnodes  # 30
                #numwide_c = 2 + 13 * nnodes  # 54
            elif self.element_type == 146:
                etype = 'VUPENTA'
                nnodes = 6
            elif self.element_type == 145:
                etype = 'VUHEXA'
                nnodes = 8
            else:
                return self._not_implemented_or_skip(data, ndata, self.code_information())

            #num_wideA = 2 + 12 * nnodes
            #ntotal = 8 + 48 * nnodes

            # assuming TETRA...
            numwide_a = 2 + (14 - 2) * nnodes  # 50
            numwide_b = 2 + (9 - 2) * nnodes  # 30
            numwide_c = 2 + 13 * nnodes  # 54
            if self.num_wide == numwide_a:
                ntotal = numwide_a * 4
                s1 = self.struct_2i
                s2 = Struct(b(self._endian + 'i11f'))
                nelements = ndata // ntotal  # 2+16*9 = 146 -> 146*4 = 584
                for i in range(nelements):
                    edata = data[n:n+8]
                    out = s1.unpack(edata)
                    (eid_device, parent_id) = out
                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)

                    for i in range(nnodes):
                        edata = data[n:n+48]
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))
                        assert len(out) == 12
                        (grid, xnorm, ynorm, znorm, txy, tyz, txz,
                         prin1, prin2, prin3, smean, vono_roct) = out
                return ndata
            elif self.num_wide == numwide_b:
                ntotal = numwide_b * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            elif self.num_wide == numwide_c:
                ntotal = numwide_c * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            else:
                msg = 'numwide=%s A=%s B=%s C=%s' % (self.num_wide, numwide_a, numwide_b, numwide_c)
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type == 139:
            # 139-QUAD4FD
            if self.read_mode == 1:
                return ndata
            if self.num_wide == 30:
                if self.is_stress():
                    self.create_transient_object(self.hyperelastic_cquad4_strain, HyperelasticQuad)
                    result_name = 'hyperelastic_cquad4_strain'
                else:
                    msg = 'HyperelasticQuad???'
                    return self._not_implemented_or_skip(data, ndata, msg)

                self._results._found_result(result_name)
                n = 0
                ntotal = 120  # 36+28*3
                s1 = Struct(b(self._endian + 'i4si6f'))  # 1 + 4+1+6 = 12
                s2 = Struct(b(self._endian + 'i6f'))
                nelements = ndata // ntotal
                obj = self.obj
                for i in range(nelements):
                    edata = data[n:n+36]  # 4*9
                    out = s1.unpack(edata)
                    if self.is_debug_file:
                        self.binary_debug.write('CQUAD4FD-139A- %s\n' % (str(out)))

                    (eid_device, Type, ID, sx, sy, sxy, angle, smj, smi) = out
                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    obj.add_new_eid(dt, [eid, Type, sx, sy, sxy, angle, smj, smi])
                    n += 36

                    for i in range(3):
                        edata = data[n:n + 28]  # 4*7
                        out = s2.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('               %s\n' % (str(out)))
                        (ID, sx, sy, sxy, angle, smj, smi) = out
                        obj.add(dt, eid, out)
                        n += 28
            else:
                msg = 'numwide=%s' % self.num_wide
                return self._not_implemented_or_skip(data, ndata, msg)
        elif self.element_type == 189:
            # 189-VUQUAD
            if self.element_type == 189:  # VQUAD
                if self.read_mode == 1:
                    return ndata
                #ntotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
                nnodes = 4    # 4 corner points + 1 centroid
                etype = 'VUQUAD4'
            #elif self.element_type == 144:  # CQUAD4
                #ntotal = 440  # 6+(33-7)*4 =  -> 110*4 = 440
                #nnodes = 4    # 4 corner points
                #etype = 'CQUAD4'
            #elif self.element_type == 64:  # CQUAD8
                #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
                #nnodes = 4    # centroid + 4 corner points
                #eType = 'CQUAD8'
            #elif self.element_type == 82:  # CQUADR
                #ntotal = 348  # 2+17*5 = 87 -> 87*4 = 348
                #nnodes = 4    # centroid + 4 corner points
                #eType = 'CQUAD4'  # TODO write the word CQUADR
            #elif self.element_type == 75:  # CTRIA6
                #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
                #nnodes = 3    # centroid + 3 corner points
                #eType = 'CTRIA6'
            #elif self.element_type == 70:  # CTRIAR
                #ntotal = 280  # 2+17*3 = 70 -> 70*4 = 280
                #nnodes = 3    # centroid + 3 corner points
                #eType = 'CTRIAR'  # TODO write the word CTRIAR
            else:
                return self._not_implemented_or_skip(data, ndata, self.code_information())

            numwide_real = 6 + (23 - 6) * nnodes  # real???
            numwide_imag = 6 + (33 - 6) * nnodes  # imag???

            if self.format_code == 1 and self.num_wide == numwide_real:  # real???
                ntotal = numwide_real * 4
                s2 = Struct(b(self._endian + '3i4s2i'))
                s3 = Struct(b(self._endian + 'i16f'))
                nelements = ndata // ntotal
                for i in range(nelements):
                    out = s2.unpack(data[n:n + 24])
                    (eid_device, parent, coord, icord, theta, itype) = out
                    n += 24
                    #eid = (eid_device - self.device_code) // 10
                    eid = self._check_id(eid_device, flag, stress_name, out)
                    edata = data[n:n + 68]
                    out = s3.unpack(edata)  # len=17*4
                    n += 68

                    if self.is_debug_file:
                        self.binary_debug.write('%s-%s - %s\n' % (etype, self.element_type, str(out)))

                    #obj.addNewNode(dt, eid, parent, coord, icord, theta, itype)
                    #obj.add_new_eid(eType, dt, eid, parent, coord, icord, theta, itype)
                    for node_id in range(nnodes - 1):  # nodes pts
                        edata = data[n:n + 68]
                        n += 68
                        out = s3.unpack(edata)
                        if self.is_debug_file:
                            self.binary_debug.write('              %s\n' % (str(out)))

                        (vuid, dummy, dummy2, msx, msy, mxy, dummy3, dummy4, dummy5,
                         bcx, bcy, bcxy, tyz, tzx, dummy6, dummy7, dummy8) = out
                        #obj.add(vuid, dummy, dummy2, msx, msy, mxy,
                                     #dummy3, dummy4, dummy5,
                                     #bcx, bcy, bcxy, tyz, tzx,
                                     #dummy6, dummy7, dummy8)
            elif self.num_wide == numwide_imag:
                ntotal = numwide_imag * 4
                nelements = ndata // ntotal
                n = nelements * ntotal
            else:
                msg = 'numwide=%s' % self.num_wide
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [47, 48, 189, 190]:
            # 47-AXIF2
            # 48-AXIF3
            # 190-VUTRIA
            return ndata
        elif self.element_type == 191:
            # 191-VUBEAM
            return ndata
        elif self.element_type in [50, 51, 203]:
            # 203-SLIF1D?
            # 50-SLOT3
            # 51-SLOT4
            return ndata
        elif self.element_type in [160, 161, 162, 163, 164, 165, 166, 167, 168,
                                   169, 170, 171, 172, 202,
                                   204, 218, 211, 213, 214,
                                   216, 217, 219, 220, 221, 222, 223,
                                   226, 232, 233, 235]:
            # 160-PENTA6FD
            # 161-TETRA4FD
            # 162-TRIA3FD
            # 163-HEXAFD
            # 164-QUADFD
            # 165-PENTAFD
            # 166-TETRAFD
            # 167-TRIAFD
            # 168-TRIAX3FD
            # 169-TRIAXFD
            # 170-QUADX4FD
            # 171-QUADXFD
            # 172-QUADRNL
            # 202-HEXA8FD
            # 204-PENTA6FD
            # 211-TRIAFD
            # 213-TRIAXFD
            # 214-QUADX4FD
            # 216-TETRA4FD
            # 217-TRIA3FD
            # 218-HEXAFD
            # 219-QUADFD
            # 220-PENTAFD
            # 221-TETRAFD
            # 223-QUADXFD
            # 222-TRIAX3FD
            # 226-BUSH
            # 232-QUADRLC
            # 233-TRIARLC
            # 235-CQUADR
            return ndata
        #elif self.element_type in [255]:
            #return ndata
        elif self.element_type in [271, 275]:
            if self.element_type == 271:
                # CPLSTS3
                result_name = 'cplstn3'
                nnodes = 1
                ntotal = 4 * 6
            elif self.element_type == 275:
                # CPLSTS3
                result_name = 'cplsts3'
                nnodes = 1
                ntotal = 4 * 6
            else:
                raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))
            if self.is_stress():
                obj_vector_real = RealCPLSTRNPlateStressArray
                result_name += '_stress'
            else:
                obj_vector_real = RealCPLSTRNPlateStrainArray
                result_name += '_strain'

            numwide_real = ntotal // 4
            if self.format_code == 1 and self.num_wide == numwide_real:
                #ntotal = 4 * (1 + 6 * (nnodes))
                nelements = ndata // ntotal

                #self._data_factor = 10  # TODO: why is this 10?
                if self.is_stress():
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_stress'
                else:
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_strain'
                slot = getattr(self, result_name)

                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, obj_vector_real)
                if auto_return:
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.use_vector and is_vectorized:
                n = nelements * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nelements
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    print(fromstring(data, dtype=self.idtype).size)
                    print('nelements=%s numwide=%s' % (nelements, numwide_real))
                    print('ndata=', ndata)
                    print('self.element_name=%s' % self.element_name)
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints[:, 0] // 10
                    obj.element[istart:iend] = eids

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                results = floats[:, 1:]
                print('results.shape', results.shape)

                #[oxx, oyy, ozz, txy, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)

        elif self.element_type in [276, 277, 278]:
            if self.element_type == 276:
                result_name = 'cplsts4'
                nnodes = 5  # 4 + 1
                ntotal = 4 * 32
            elif self.element_type == 277:
                result_name = 'cplsts6'
                nnodes = 4
                ntotal = 4 * 26
            elif self.element_type == 278:
                result_name = 'cplsts8'
                nnodes = 5
                ntotal = 4 * 32
            else:
                raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

            if self.is_stress():
                obj_vector_real = RealCPLSTRNPlateStressArray
                result_name += '_stress'
            else:
                obj_vector_real = RealCPLSTRNPlateStrainArray
                result_name += '_strain'

            numwide_real = 2 + 6 * (nnodes)
            assert ntotal // 4 == numwide_real, 'notal/4=%s numwide_real=%s\n%s' % (
                ntotal // 4, numwide_real, self.code_information())

            ntotal = numwide_real * 4
            if self.format_code == 1 and self.num_wide == numwide_real:
                nelements = ndata // ntotal

                #self._data_factor = 10  # TODO: why is this 10?
                if self.is_stress():
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_stress'
                else:
                    obj_vector_real = RealCPLSTRNPlateStressArray
                    #result_name = 'cplstn3_strain'
                slot = getattr(self, result_name)

                nlayers = nelements * nnodes
                auto_return, is_vectorized = self._create_oes_object4(
                    nlayers, result_name, slot, obj_vector_real)
                if auto_return:
                    self._data_factor = nnodes
                    return nelements * self.num_wide * 4

                obj = self.obj
                #if self.use_vector and is_vectorized:
                n = nlayers * self.num_wide * 4

                istart = obj.itotal
                iend = istart + nlayers
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    print(fromstring(data, dtype=self.idtype).size)
                    print('nelements=%s numwide=%s' % (nelements, numwide_real))
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, numwide_real)
                    eids = ints[:, 0] // 10
                    #obj.element[istart:iend] = eids

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, numwide_real)
                print('floats[:, 2:].shape', floats[:, 2:].shape)
                print('nnelements=%s nnodes=%s numwide//nodes=%s' % (nelements, nnodes, (numwide_real-2) / nnodes))
                results = floats[:, 2:].reshape(nelements, nnodes * 6)

                #[oxx, oyy, ozz, txy, ovm]
                obj.data[obj.itime, istart:iend, :] = results
            else:
                msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
                return self._not_implemented_or_skip(data, ndata, msg)
        else:
            msg = 'sort1 Type=%s num=%s' % (self.element_name, self.element_type)
            return self._not_implemented_or_skip(data, ndata, msg)

        assert ndata > 0, ndata
        assert nelements > 0, 'nelements=%r element_type=%s element_name=%r' % (nelements, self.element_type, self.element_name)
        #assert ndata % ntotal == 0, '%s n=%s nwide=%s len=%s ntotal=%s' % (self.element_name, ndata % ntotal, ndata % self.num_wide, ndata, ntotal)
        assert self.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (self.num_wide * 4, ntotal)
        assert self.thermal == 0, "thermal = %%s" % self.thermal
        assert n > 0, "n = %s result_name=%s" % (n, result_name)
        return n

    def OESRT_CQUAD4_95(self, data, ndata):
        """unsupported element"""
        assert self.num_wide == 9, "num_wide=%s not 9" % self.num_wide
        ntotal = 36  # 4*9

        n = 0
        s = Struct(b(self._endian + 'i8si3fi4s'))
        nelements = ndata // ntotal
        obj = self.obj
        for i in range(nelements):
            edata = data[n:n + ntotal]
            out = s.unpack(edata)  # num_wide=9
            if self.is_debug_file:
                self.binary_debug.write('CQUAD4-95 - %s\n' % str(out))
            #eid, failure, ply, failureIndexPly, failureIndexBonding, failureIndexMax, flag
            # 3,TSAIWU,1,8.5640,0.0,None

            (eid, failure, ply, strengthRatioPly, failureIndexBonding, strengthRatioBonding, flag, flag2) = out
            #strengthRatioPly
            #print("eid=%s failure=%r ply=%s failureIndexPly=%s  failureIndexBonding=%s strengthRatioBonding=%s flag=%s flag2=%s" % (eid, failure.strip(), ply, failureIndexPly, failureIndexBonding, strengthRatioBonding, flag, flag2))
            print("eid=%s strengthRatioPly=%g failureIndexBonding=%s strengthRatioBonding=%s" % (eid, strengthRatioPly, failureIndexBonding, strengthRatioBonding))
            #obj.add_new_eid(element_name, dt, eid, force, stress)
            n += ntotal

    def _create_oes_object3(self, nelements, result_name, slot, obj, obj_vector):
        """
        Creates the self.obj parameter based on if this is vectorized or not.

        Parameters
        ----------
        self : OES()
            the object pointer
        nelements :  int
            the number of elements to preallocate for vectorization
        result_name : str
            unused
        slot : dict[(int, int, str)=obj
            the self dictionary that will be filled with a
            non-vectorized result
        obj : OES
            a pointer to the non-vectorized class
        obj_vector : OESArray
            a pointer to the vectorized class

        Returns
        -------
        auto_return : bool
            a flag indicating a return n should be called
        is_vectorized : bool
            True/False

        Since that's confusing, let's say we have real CTETRA stress data.
        We're going to fill self.solidStress with the class
        RealSolidStress.  If it were vectorized, we'd fill
        self.ctetra_stress. with RealSolidStressArray.  So we call:

        if self._is_vectorized(RealSolidStressArray, self.ctetra_stress):
            if self._results.is_not_saved(result_vector_name):
                return ndata
        else:
            if self._results.is_not_saved(result_name):
                return ndata

        auto_return, is_vectorized = self._create_oes_object3(self, nelements,
                            'ctetra_stress', self.ctetra_stress,
                            RealSolidStress, RealSolidStressArray)
        if auto_return:
            return nelements * ntotal
        """
        auto_return = False
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print('is_vectorized=%s result_name=%r' % (is_vectorized, result_name))
        if is_vectorized:
            #print("vectorized...read_mode=%s...%s" % (self.read_mode, result_name))
            if self.read_mode == 1:
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nelements += nelements
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                self.obj.build()

        else:  # not vectorized
            self.code = self._get_code()
            self.result_names.add(result_name)
            #print("not vectorized...read_mode=%s...%s" % (self.read_mode, result_name))
            #self.log.info("code = %s" % str(self.code))
            if self.read_mode == 1:
                self.result_names.add(result_name)
                auto_return = True
            # pass = 0/2
            self.create_transient_object(slot, obj)

        if auto_return and self.read_mode == 2:
            raise RuntimeError('this should never happen...auto_return=True read_mode=2')
        return auto_return, is_vectorized

    def _create_oes_object4(self, nelements, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it doesn't support unvectorized objects"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nelements += nelements
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1...%s" % result_name
                    self.log.error(msg)
                    raise
                #obj.update_data_code(self.data_code)
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _create_nodes_object(self, nnodes, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the nnodes parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nnodes += nnodes
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _create_ntotal_object(self, ntotal, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it adds to the ntotal parameter"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.ntotal += ntotal
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _is_vectorized(self, obj_vector, slot_vector):
        """
        Checks to see if the data array has been vectorized

        Parameters
        ----------
        obj_vector:  the object to check
            (obj or None; None happens when vectorization hasn't been implemented)
        slot_vector: the dictionary to put the object in
            (dict or None; None happens when obj hasn't been implemented)

        Returns
        -------
        is_vectorized : bool
            should the data object be vectorized

        .. note :: the Vectorized column refers to the setting given by the user
        """
        is_vectorized = False
        #print('self.is_vectorized =', is_vectorized)
        if self.is_vectorized:
            if obj_vector is not None:
                is_vectorized = True
                #print("***vectorized...")
            #else:
                #print("***not vectorized...")
                #self.log.info('obj_vector=%s slot_vector=%s' % (obj_vector, slot_vector))
        return is_vectorized
