#VTK_TRIANGLE = 5
#VTK_QUADRATIC_TRIANGLE = 22

#VTK_QUAD = 9
#VTK_QUADRATIC_QUAD = 23

#VTK_TETRA = 10
#VTK_QUADRATIC_TETRA = 24

#VTK_WEDGE = 13
#VTK_QUADRATIC_WEDGE = 26

#VTK_HEXAHEDRON = 12
#VTK_QUADRATIC_HEXAHEDRON = 25



import vtk
from vtk import (vtkTriangle,vtkQuad,vtkTetra,vtkWedge,vtkHexahedron,
                 vtkQuadraticTriangle,vtkQuadraticQuad,vtkQuadraticTetra,vtkQuadraticWedge,vtkQuadraticHexahedron)

from pyNastran.bdf.bdf import *
from pyNastran.op2.op2 import OP2


class NastranIO(object):
    def __init__(self):
        pass

    def loadGeometry(self,bdfFileName,dirname):
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]
        
        if bdfFileName is None:
            self.grid       = vtk.vtkUnstructuredGrid()
            self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.grid2      = vtk.vtkUnstructuredGrid()
            self.scalarBar.VisibilityOff()
            return
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()
            self.gridResult = vtk.vtkFloatArray()
            self.gridResult.Reset()
            self.gridResult.Modified()
            self.eidMap = {}
            self.nidMap = {}

            self.resultCases = {}
            self.nCases = 0
            try:
                del self.caseKeys
                del self.iCase
                del self.iSubcaseNameMap
            except:
                print "cant delete geo"
                pass
            ###
            #print dir(self)
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        model = BDF()
        model.readBDF(bdfFileName,includeDir=dirname)

        nNodes    = model.nNodes()
        nElements = model.nElements()
        nCAeros   = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        #print "nNodes = ",nNodes
        print "nElements = ",nElements

        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.cardCount:
            nCONM2 = model.cardCount['CONM2']
        else:
            nCONM2 = 0
        self.grid.Allocate(nElements, 1000)
        self.gridResult.SetNumberOfComponents(nElements)
        self.grid2.Allocate(nCAeros +nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nNodes)
        self.gridResult.Allocate(nNodes,1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            i=0
            fraction = 1./nNodes # so you can color the nodes by ID
            for nid,node in sorted(model.nodes.iteritems()):
                #print "i = ",i
                point = node.Position()
                #print "point = ",point
                points.InsertPoint(i, *point)
                self.gridResult.InsertNextValue(i*fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

                self.nidMap[nid] = i
                i+=1
        if 1:
            i = 0
            for nid,node in sorted(model.nodes.iteritems()):
                point = node.Position()
                points.InsertPoint(i, *point)
                self.nidMap[nid] = i
                i+=1
            #print "nidMap = ",self.nidMap

        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros*4+nCONM2)
        for eid,element in sorted(model.caeros.iteritems()):
            if (isinstance(element,CAERO1) or isinstance(element,CAERO3) or
                isinstance(element,CAERO4)  or isinstance(element,CAERO5)):
                cpoints = element.Points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j+1)
                elem.GetPointIds().SetId(2, j+2)
                elem.GetPointIds().SetId(3, j+3)
                points2.InsertPoint(j,   *cpoints[0])
                points2.InsertPoint(j+1, *cpoints[1])
                points2.InsertPoint(j+2, *cpoints[2])
                points2.InsertPoint(j+3, *cpoints[3])
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=4
            #elif isinstance(element,CAERO2): # cylinder
                #pass
            else:
                print "skipping %s" %(element.type)
                
        self.mapElements(points,points2,self.nidMap,model,j)

    def mapElements(self,points,points2,nidMap,model,j):
        self.eidMap = {}
        i = 0
        for eid,element in sorted(model.elements.iteritems()):
            self.eidMap[eid] = i
            #print element.type
            if isinstance(element,CTRIA3) or isinstance(element,CTRIAR):
                #print "ctria3"
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTRIA6):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTRIAX6): # midside nodes are required, nodes out of order
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[1]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[4]])
                #a = [0,2,4]
                #msg = "CTRIAX6 %i %i %i" %(nidMap[nodeIDs[a[0]]], nidMap[nodeIDs[a[1]]], nidMap[nodeIDs[a[2]]] )
                #raise Exception(msg)
                #sys.stdout.flush()
                
                #elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                #elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                #elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element,CQUAD4) or isinstance(element,CSHEAR) or isinstance(element,CQUADR):
                nodeIDs = element.nodeIDs()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD8):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA10):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTetra()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                else:
                    elem = vtkTetra()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element,CPENTA15):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6,  nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7,  nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                else:
                    elem = vtkWedge()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA8):
                nodeIDs = element.nodeIDs()
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA20):
                nodeIDs = element.nodeIDs()
                #print "nodeIDs = ",nodeIDs
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                    elem.GetPointIds().SetId(15, nidMap[nodeIDs[15]])
                    elem.GetPointIds().SetId(16, nidMap[nodeIDs[16]])
                    elem.GetPointIds().SetId(17, nidMap[nodeIDs[17]])
                    elem.GetPointIds().SetId(18, nidMap[nodeIDs[18]])
                    elem.GetPointIds().SetId(19, nidMap[nodeIDs[19]])
                else:
                    elem = vtkHexahedron()

                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,LineElement) or isinstance(element,SpringElement):
                elem = vtk.vtkLine()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            ###
            elif isinstance(element,CONM2): # not perfectly located
                del self.eidMap[eid]; i-=1
                nid  = element.Nid()
                c    = element.Centroid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)

                points2.InsertPoint(j,*c)
                elem.GetPointIds().SetId(0, j)
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=1
            else:
                del self.eidMap[eid]; i-=1
                print "skipping %s" %(element.type)
            i+=1
        ###
        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print dir(self.grid) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        print "updated grid"


    def loadResults(self,op2FileName):
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        op2 = OP2(op2FileName,debug=True)
        op2.readOP2()
        #print op2.printResults()
        
        #case = op2.displacements[1]
        #print "case = ",case
        #for nodeID,translation in sorted(case.translations.iteritems()):
            #print "nodeID=%s t=%s" %(nodeID,translation)
        #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = op2.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = op2.iSubcaseNameMap
        
        nElements = len(self.eidMap)
        #print "nElements = ",nElements
        nidsSet = False # set to False to disable nodeIDs
        eidsSet = True
        for ID in subcaseIDs:
            if nidsSet:
                nids = zeros(self.nNodes,'d')
                for nid,nid2 in self.nidMap.iteritems():
                    nids[nid2] = nid
                cases[(ID,'Node_ID',1,'node','%.0f')] = nids
                nidsSet = True

            if eidsSet:
                eids = zeros(nElements,'d')
                for eid,eid2 in self.eidMap.iteritems():
                    eids[eid2] = eid
               
                eKey = (ID,'isElementOn',1,'centroid','%.0g')
                cases[(ID,'Element_ID',1,'centroid','%.0f')] = eids
                cases[eKey] = zeros(nElements) # is the element supported
                eidsSet = True
            
            if False:
                if ID in op2.displacements: # not correct?
                    case = op2.displacements[ID]
                    key = (ID,'DisplacementX',3,'node','%g')
                    #cases[key] = case.translations

                if ID in op2.temperatures:
                    case = op2.temperatures[ID]
                    #print case
                    temps = zeros(self.nNodes)
                    key = (ID,'Temperature',1,'node','%g')
                    for nid,T in case.temperatures.iteritems():
                        #print T
                        nid2 = self.nidMap[nid]
                        temps[nid2] = T
                    ###
                    #cases[key] = temps

            if self.isStress(op2,ID):
                cases = self.fillStressCase(cases,op2,ID,eKey,nElements)
            ###
        ###
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases)-1 # number of keys in dictionary
        self.cycleResults() # start at nCase=0

    def fillStressCase(self,cases,op2,ID,eKey,nElements):
        oxx = zeros(nElements)
        oyy = zeros(nElements)
        ozz = zeros(nElements)

        o1  = zeros(nElements)
        o2  = zeros(nElements)
        o3  = zeros(nElements)
        ovm = zeros(nElements)

        vmWord = 'N/A'
        if ID in op2.rodStress:
            case = op2.rodStress[ID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.
                #print "bar eid=%s" %(eid)

                axial   = case.axial[eid]
                torsion = case.torsion[eid]

                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2]  = max(axial,torsion)  # not really
                #oyy[eid2] = torsion
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = min(axial,torsion)

        if ID in op2.barStress:
            #self.s1    = {}
            #self.s2    = {}
            #self.s3    = {}
            #self.s4    = {}
            case = op2.barStress[ID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "bar eid=%s" %(eid)
                oxxi = case.axial[eid]
                o1i  = max(case.smax[eid])
                o3i  = min(case.smin[eid])
                
                oxx[eid2] = oxxi
                #oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                #ovm[eid2] = ovmi

        if ID in op2.plateStress:
            #self.txy    = {}
            case = op2.plateStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "plate eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                #self.oyy[eid][nid][iLayer]

                oxxi = max(case.oxx[eid]['C'])
                oyyi = max(case.oyy[eid]['C'])
                ozzi = min(case.oxx[eid]['C'],min(case.oyy[eid]['C']))

                o1i = max(case.majorP[eid]['C'])
                o2i = max(case.minorP[eid]['C'])
                o3i = min(case.majorP[eid]['C'],min(case.minorP[eid]['C']))
                ovmi = max(case.ovmShear[eid]['C'])


                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if ID in op2.solidStress:
            case = op2.solidStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "solid eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                oyyi = case.oyy[eid]['C']
                ozzi = case.ozz[eid]['C']
                o1i = case.o1[eid]['C']
                o2i = case.o2[eid]['C']
                o3i = case.o3[eid]['C']
                ovmi = case.ovmShear[eid]['C']
                #if ID==1:
                    #print "ovm[%s] = %s" %(eid,ovmi)
                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = o2i
                o3[eid2] = o3i
                ovm[eid2] = ovmi
            ###

        cases[(ID,'StressXX',1,'centroid','%.3f')] = oxx
        cases[(ID,'StressYY',1,'centroid','%.3f')] = oyy
        cases[(ID,'StressZZ',1,'centroid','%.3f')] = ozz

        cases[(ID,'Stress1',1,'centroid','%.3f')] = o1
        cases[(ID,'Stress2',1,'centroid','%.3f')] = o2
        cases[(ID,'Stress3',1,'centroid','%.3f')] = o3
        cases[(ID,vmWord   ,1,'centroid','%.3f')] = ovm
        return cases
    
