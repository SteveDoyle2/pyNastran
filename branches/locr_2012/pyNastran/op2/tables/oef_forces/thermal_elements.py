from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from struct import unpack

from .oef_thermalObjects import (HeatFlux_CHBDYx, HeatFlux_2D_3D, HeatFlux_1D,
                                 HeatFlux_VU, HeatFlux_VUBEAM, HeatFlux_VU_3D,
                                 HeatFlux_CONV)


class ThermalElements(object):

    def readOEF_Thermal(self):
        assert self.format_code == 1, self.code_information()
        #print "self.element_type = ",self.element_type
        if self.element_type in [107, 108, 109]:  # CHBDYE, CHBDYG, CHBDYP
            assert self.num_wide == 8, self.code_information()
            self.create_transient_object(self.thermalLoad_CHBDY, HeatFlux_CHBDYx)
            self.handle_results_buffer(self.OEF_CHBDYx,
                                       resultName='thermalLoad_CHBDY')
        elif self.element_type in [33, 39, 67, 68]:  # QUAD4,TETRA,HEXA,PENTA
            assert self.num_wide in [9, 10], self.code_information()
            self.create_transient_object(self.thermalLoad_2D_3D, HeatFlux_2D_3D)
            self.handle_results_buffer(self.OEF_2D_3D,
                                       resultName='thermalLoad_2D_3D')
        elif self.element_type in [53, 64, 74, 75]:  # TRIAX6,QUAD8,TRIA3,TRIA6
            assert self.num_wide == 9, self.code_information()
            self.create_transient_object(self.thermalLoad_2D_3D, HeatFlux_2D_3D)
            self.handle_results_buffer(self.OEF_2D_3D,
                                       resultName='thermalLoad_2D_3D')
        elif self.element_type in [1, 2, 3, 10, 34, 69]:  # ROD,BEAM,TUBE,CONROD,BAR,BEND
            assert self.num_wide == 9, self.code_information()
            self.create_transient_object(self.thermalLoad_1D, HeatFlux_1D)
            self.handle_results_buffer(self.OEF_1D, resultName='thermalLoad_1D')
        elif self.element_type in [189, 190]:  # VUQUAD,VUTRIA
            #assert self.num_wide==27,self.code_information()
            self.create_transient_object(self.thermalLoad_VU, HeatFlux_VU)
            self.handle_results_buffer(
                self.OEF_VU_Element, resultName='thermalLoad_VU')
        elif self.element_type in [191]:  # VUBEAM
            #assert self.num_wide==27,self.code_information()
            self.create_transient_object(
                self.thermalLoad_VUBeam, HeatFlux_VUBEAM)
            self.handle_results_buffer(self.OEF_VUBeam_Element,
                                      resultName='thermalLoad_VUBeam')
        elif self.element_type in [145, 146, 147]:  # VUHEXA,VUPENTA,VUTETRA
            self.create_transient_object(self.thermalLoad_VU_3D, HeatFlux_VU_3D)
            self.handle_results_buffer(self.OEF_VU_3D_Element,
                                      resultName='thermalLoad_VU_3D')
        elif self.element_type in [110]:
            self.create_transient_object(self.thermalLoad_CONV, HeatFlux_CONV)
            self.handle_results_buffer(self.OEF_CONV,
                                       resultName='thermalLoad_CONV')
        else:
            self.not_implemented_or_skip()

    def OEF_CHBDYx(self):  # [107,108,109]  CHBDYE, CHBDYG, CHBDYP
        if self.make_op2_debug:
            self.op2Debug.write('---OEF_CHBDYx---\n')
        #device_code = self.device_code

        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'i8s5f'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'f8s5f'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        format1 = bytes(format1)

        while len(self.data) >= 32:  # 8*4
            eData = self.data[0:32]
            self.data = self.data[32:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, eType, fApplied, freeConv, forceConv, fRad, fTotal) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, eType, fApplied, freeConv, forceConv, fRad, fTotal]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_CHBDYx")
        #print(self.thermalLoad_CHBDY)

    def OEF_CONV(self):  # [110]  CONV
        if self.make_op2_debug:
            self.op2Debug.write('---OEF_CONV---\n')
        #device_code = self.device_code

        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'ifif'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'ffif'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        format1 = bytes(format1)

        while len(self.data) >= 16:  # 4*4
            eData = self.data[0:16]
            self.data = self.data[16:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, cntlNode, freeConv, freeConvK) = out
            eid2 = extract(eid, dt)

            dataIn = [eid2, cntlNode, freeConv, freeConvK]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)

        if self.make_op2_debug:
            print("done with OEF_CONV")
        #print(self.thermalLoad_CHBDY)

    def OEF_VU_Element(self):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        #print "num_wide = ",self.num_wide

        if self.element_type in [189]:
            nNodes = 4
        elif self.element_type in [190]:
            nNodes = 3
        elif self.element_type in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.code_information())

        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'iii4sii'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'fii4sii'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor

        formatAll = 'iffffff'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 24 + 28 * nNodes
        while len(self.data) >= n:
            eData = self.data[0:24]  # 6*4
            self.data = self.data[24:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, parent, coord, icord, theta, null) = out

            eid2 = extract(eid, dt)
            dataIn = [eid2, parent, coord, icord, theta]

            gradFluxes = []
            for i in xrange(nNodes):
                eData = self.data[0:28]  # 7*4
                self.data = self.data[28:]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)

            #dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(nNodes, dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_1D")
        #print self.thermalLoad_VU

    def OEF_VUBeam_Element(self):  # 191-VUBEAM
        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        #print "num_wide = ",self.num_wide

        if self.element_type in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.code_information())

        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'iii4s'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'fii4s'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        formatAll = 'i6f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 16 + 28 * nNodes
        while len(self.data) >= n:
            eData = self.data[0:16]  # 4*4
            self.data = self.data[16:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, parent, coord, icord) = out

            eid2 = extract(eid, dt)
            dataIn = [eid2, parent, coord, icord]

            gradFluxes = []
            for i in xrange(nNodes):
                eData = self.data[0:28]  # 7*4
                self.data = self.data[28:]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)

            #dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(nNodes, dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_1D")
        #print self.thermalLoad_VUBeam

    def OEF_VU_3D_Element(self):  # 146-VUPENTA, 147-VUTETRA, 148-VUPENTA
        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        #print "num_wide = ",self.num_wide

        if self.element_type in [147]:
            nNodes = 4
        elif self.element_type in [146]:
            nNodes = 6
        elif self.element_type in [145]:
            nNodes = 8
        else:
            raise NotImplementedError(self.code_information())

        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'ii'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'fi'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        formatAll = 'i6f'
        format1 = bytes(format1)
        formatAll = bytes(formatAll)

        n = 8 + 28 * nNodes
        while len(self.data) >= n:
            eData = self.data[0:8]  # 2*4
            self.data = self.data[8:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, parent) = out

            eid2 = extract(eid, dt)
            dataIn = [eid2, parent]

            gradFluxes = []
            for i in xrange(nNodes):
                eData = self.data[0:7 * 4]
                self.data = self.data[7 * 4:]
                #print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)

            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            self.obj.add(nNodes, dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_VU_3D_Element")
        #print self.thermalLoad_VU_3D

    def OEF_1D(self):  # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        #print "num_wide = ",self.num_wide
        n = 36
        if is_sort1:
            #print "SORT1 - %s" %(self.get_element_type(self.element_type))
            format1 = 'i8sffffff'  # SORT1
            extract = self.extractSort1
            #dt = self.nonlinear_factor
        else:
            #print "SORT2 - %s" %(self.get_element_type(self.element_type))
            format1 = 'f8sffffff'  # SORT2
            extract = self.extractSort2
            #eid = self.nonlinear_factor
        format1 = bytes(format1)

        n = 36
        while len(self.data) >= n:  # 10*4
            eData = self.data[0:n]
            self.data = self.data[n:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_1D")
        #print self.thermalLoad_1D

    def OEF_2D_3D(self):  # 33-QUAD4, 39-TETRA, 53-TRIAX6,64-QUAD8, 67-HEXA, 68-PENTA, 74-TRIA3, 75-TRIA6
        """num_wide==10"""
        dt = self.nonlinear_factor
        is_sort1 = self.is_sort1()
        #print "num_wide = ",self.num_wide
        #print "dt = ",dt
        if self.element_type in [39, 67, 68]:  # HEXA,PENTA
            n = 40
            if is_sort1:
                #print "SORT1 - %s" %(self.get_element_type(self.element_type))
                format1 = 'i8s6fi'  # SORT1
                extract = self.extractSort1
                #dt = self.nonlinear_factor
            else:
                #print "SORT2 - %s" %(self.get_element_type(self.element_type))
                format1 = 'f8s6fi'  # SORT2
                extract = self.extractSort2
                #eid = self.nonlinear_factor
        elif self.element_type in [33, 53, 64, 74, 75]:  # no zed on this element for some reason...
            n = 36
            if is_sort1:
                #print "SORT1 - %s" %(self.get_element_type(self.element_type))
                format1 = 'i8s6f'  # SORT1
                extract = self.extractSort1
                #dt = self.nonlinear_factor
            else:
                #print "SORT2 - %s" %(self.get_element_type(self.element_type))
                format1 = 'f8s6f'  # SORT2
                extract = self.extractSort2
                #eid = self.nonlinear_factor
        else:
            raise NotImplementedError(self.code_information())
        format1 = bytes(format1)

        while len(self.data) >= n:  # 10*4
            eData = self.data[0:n]
            self.data = self.data[n:]

            out = unpack(format1, eData)
            #print("len(out)=",len(out))
            # len=8
            (eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux) = out[:8]
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux]
            #print "heatFlux %s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)
        if self.make_op2_debug:
            print("done with OEF_2D_3D")

        #print self.thermalLoad_2D_3D
