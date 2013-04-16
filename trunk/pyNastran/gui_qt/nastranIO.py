# pylint: disable=C0103,C0111,E1101


from numpy import zeros

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron)

from pyNastran.bdf.bdf import (BDF, CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                               CQUAD4, CQUAD8, CQUADR, CSHEAR,
                               CTRIA3, CTRIA6, CTRIAR, CTRIAX6,
                               CTETRA4, CTETRA10, CPENTA6, CPENTA15,
                               CHEXA8, CHEXA20,
                               CONM2,
                               LineElement, SpringElement)
from pyNastran.op2.op2 import OP2


class NastranIO(object):


    def load_nastran_geometry(self, bdfFileName, dirname, isNodal, isCentroidal):
        self.isNodal = isNodal
        self.isCentroidal = isCentroidal

        self.gridResult = vtk.vtkFloatArray()

        if bdfFileName is None:
            self.grid = vtk.vtkUnstructuredGrid()
            self.grid2 = vtk.vtkUnstructuredGrid()
            #self.scalarBar.VisibilityOff()
            return

        #self.TurnTextOff()
        self.grid.Reset()
        self.grid2.Reset()
        self.gridResult.Modified()
        self.eidMap = {}
        self.nidMap = {}

        self.resultCases = {}
        self.nCases = 0        

        for i in ('caseKeys', 'iCase', 'iSubcaseNameMap'):
            if hasattr(self, i):
                del i

        #self.scalarBar.VisibilityOff()
        #self.scalarBar.Modified()

        self.log_info("reading file " + bdfFileName)
        model = BDF(debug = True, log = self.log)
        self.modelType = model.modelType
        model.readBDF(bdfFileName, includeDir=dirname)

        nNodes = model.nNodes()
        nElements = model.nElements()
        nCAeros = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        self.log_info("nElements = %i" % (self.nElements))
        
        nCONM2 = model.cardCount['CONM2'] if 'CONM2' in model.cardCount else 0
        
        self.grid.Allocate(self.nElements, 1000)
        self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(nCAeros + nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        self.gridResult.Allocate(self.nNodes, 1000)
        self.nidMap = {}

        for i, (nid, node) in enumerate(sorted(model.nodes.iteritems())):
            point = node.Position()
            points.InsertPoint(i, *point)
            self.nidMap[nid] = i

        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros * 4 + nCONM2)
        for (eid, element) in sorted(model.caeros.iteritems()):
            if (isinstance(element, CAERO1) or isinstance(element, CAERO3) or
                isinstance(element, CAERO4) or isinstance(element, CAERO5)):
                cpoints = element.Points()
                elem = vtkQuad()
                for _i in range(4):
                    elem.GetPointIds().SetId(_i, j + _i)
                    points2.InsertPoint(j + _i, *cpoints[_i])
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 4
            else:
                self.log_info("skipping %s" % (element.type))

        self.mapElements(points, points2, self.nidMap, model, j)

    def mapElements(self, points, points2, nidMap, model, j): 
        _mk2 = lambda x: (x[0], x[1]) if isinstance(x, tuple) else (x, x)
        # parameter in a list can be given as one number or a tuple of two numbers
        # the idea is: lst = [1, (2,3)] --> map(_mk2, lst) = [(1,1), (2,3)]
        
        def _simple_set(obj, ids): # create new vtk object and set elements
            elem = obj()
            for i, j in map(_mk2, ids):
                elem.GetPointIds().SetId(i, nidMap[nodeIDs[j]])
            return elem
        # conditional new vtk object creating and initialisation
        def _complex_set(nids, obj1, obj2, ids, sec_ids): 
            if None not in nids:
                elem = obj1()
                for i, j in map(_mk2, ids):
                    elem.GetPointIds().SetId(i, nidMap[nodeIDs[j]])
            else:
                elem = obj2()
            for i, j in map(_mk2, sec_ids):
                elem.GetPointIds().SetId(i, nidMap[nodeIDs[j]])
            return elem
        
                
        self.eidMap = {}
        i = 0
        for (eid, element) in sorted(model.elements.iteritems()):
            self.eidMap[eid] = i
            nodeIDs = element.nodeIDs() if hasattr(element, 'nodeIDs') else None
                
            if isinstance(element, CTRIA3) or isinstance(element, CTRIAR):
                elem = _simple_set(vtkTriangle, [0, 1, 2])
            elif isinstance(element, CTRIA6):
                elem = _complex_set(nodeIDs, vtkQuadraticTriangle, vtkTriangle,
                                    [3, 4, 5], [0, 1, 2])
            elif isinstance(element, CTRIAX6):
                # midside nodes are required, nodes out of order
                elem = _complex_set(nodeIDs, vtkQuadraticTriangle, vtkTriangle,
                                    [(3, 1), (4, 3), 5], [0, (1, 2), (2, 4)])
            elif (isinstance(element, CQUAD4) or isinstance(element, CSHEAR) or
                  isinstance(element, CQUADR)):
                elem = _simple_set(vtkQuad, range(4))
            elif isinstance(element, CQUAD8):
                elem = _complex_set(nodeIDs, vtkQuadraticQuad, vtkQuad,
                                    range(4, 8), range(4))
            elif isinstance(element, CTETRA4):
                elem = _simple_set(vtkTetra, range(4))
            elif isinstance(element, CTETRA10):
                elem = _complex_set(nodeIDs, vtkQuadraticTetra, vtkTetra,
                                    range(4, 10), range(4))
            elif isinstance(element, CPENTA6):
                elem = _simple_set(vtkWedge, range(6))
            elif isinstance(element, CPENTA15):
                elem = _complex_set(nodeIDs, vtkQuadraticWedge, vtkWedge,
                                    range(6, 15), range(6))
            elif isinstance(element, CHEXA8):
                elem = _simple_set(vtkHexahedron, range(8))
            elif isinstance(element, CHEXA20):
                elem = _complex_set(nodeIDs, vtkQuadraticHexahedron, vtkHexahedron,
                                    range(8, 20), range(8))
            elif (isinstance(element, LineElement) or
                  isinstance(element, SpringElement)):
                elem = _simple_set(vtk.vtkLine, [0, 1])
            elif isinstance(element, CONM2):  # not perfectly located
                del self.eidMap[eid]
                i -= 1

                elem = vtk.vtkVertex()
                points2.InsertPoint(j, *element.Centroid() )
                elem.GetPointIds().SetId(0, j)
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
                i += 1
                continue # we already inserted cell to grid2 and do not want inserting to grid 
            else:
                del self.eidMap[eid]
                self.log_info("skipping %s" % (element.type))
                continue
            self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            i += 1

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
        self.log_info("updated grid")

    def loadNastranResults(self, op2FileName, dirname, isNodal, isCentroidal):
        #self.gridResult.SetNumberOfComponents(self.nElements)
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        op2 = OP2(op2FileName, debug=True)
        op2.readOP2()
        #print op2.print_results()

        #case = op2.displacements[1]
        #print "case = ",case
        #for nodeID,translation in sorted(case.translations.iteritems()):
            #print "nodeID=%s t=%s" %(nodeID,translation)
        #self.iSubcaseNameMap[self.isubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = op2.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = op2.iSubcaseNameMap

        nElements = len(self.eidMap)
        #print "nElements = ",nElements
        nidsSet = False  # set to False to disable nodeIDs
        eidsSet = True
        for subcaseID in subcaseIDs:
            if nidsSet:
                nids = zeros(self.nNodes, 'd')
                for (nid, nid2) in self.nidMap.iteritems():
                    nids[nid2] = nid
                cases[(subcaseID, 'Node_ID', 1, 'node', '%.0f')] = nids
                nidsSet = True

            if eidsSet:
                eids = zeros(nElements, 'd')
                for (eid, eid2) in self.eidMap.iteritems():
                    eids[eid2] = eid

                eKey = (subcaseID, 'isElementOn', 1, 'centroid', '%.0g')
                cases[(subcaseID, 'Element_ID', 1, 'centroid', '%.0f')] = eids
                cases[eKey] = zeros(nElements)  # is the element supported
                eidsSet = True

            if False:
                if subcaseID in op2.displacements:  # not correct?
                    case = op2.displacements[subcaseID]
                    key = (subcaseID, 'DisplacementX', 3, 'node', '%g')
                    #cases[key] = case.translations

                if subcaseID in op2.temperatures:
                    case = op2.temperatures[subcaseID]
                    #print case
                    temps = zeros(self.nNodes)
                    key = (subcaseID, 'Temperature', 1, 'node', '%g')
                    for (nid, T) in case.temperatures.iteritems():
                        #print T
                        nid2 = self.nidMap[nid]
                        temps[nid2] = T
                    #cases[key] = temps

            if self.isStress(op2, subcaseID):
                cases = self.fillStressCase(cases, op2, subcaseID,
                                            eKey, nElements)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
        self.cycleResults()  # start at nCase=0

    def fillStressCase(self, cases, op2, subcaseID, eKey, nElements):
        oxx = zeros(nElements)
        oyy = zeros(nElements)
        ozz = zeros(nElements)
        o1 = zeros(nElements)
        o2 = zeros(nElements)
        o3 = zeros(nElements)
        ovm = zeros(nElements)

        vmWord = 'N/A'
        if subcaseID in op2.rodStress:
            case = op2.rodStress[subcaseID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.
                oxx[eid2] = case.axial[eid]   # axial
                oyy[eid2] = case.torsion[eid] # torsion

                o1[eid2] = max(oxx[eid2], oyy[eid2])  # not really
                o3[eid2] = min(oxx[eid2], oyy[eid2])

        if subcaseID in op2.barStress:
            case = op2.barStress[subcaseID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.
                
                oxx[eid2] = case.axial[eid]
                o1[eid2] = max(case.smax[eid])
                o3[eid2] = min(case.smin[eid])

        if subcaseID in op2.plateStress:
            #self.txy    = {}
            case = op2.plateStress[subcaseID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                oxx[eid2] = max(case.oxx[eid]['C'])
                oyy[eid2] = max(case.oyy[eid]['C'])
                #ozz[eid2] = min(case.oxx[eid]['C'], min(case.oyy[eid]['C']))

                o1[eid2] = max(case.majorP[eid]['C'])
                o3[eid2] = min(case.majorP[eid]['C'], min(case.minorP[eid]['C']))
                ovm[eid2] = max(case.ovmShear[eid]['C'])

        if subcaseID in op2.solidStress:
            case = op2.solidStress[subcaseID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.


                oxx[eid2] = case.oxx[eid]['C']
                oyy[eid2] = case.oyy[eid]['C']
                ozz[eid2] = case.ozz[eid]['C']

                o1[eid2] = case.o1[eid]['C']
                o2[eid2] = case.o2[eid]['C']
                o3[eid2] = case.o3[eid]['C']
                ovm[eid2] = case.ovmShear[eid]['C']

        # subcaseID,resultType,vectorSize,location,dataFormat
        for (i,j) in [('StressXX', oxx), ('StressYY', oyy), ('StressZZ', ozz),
            ('Stress1', o1), ('Stress2', o2),('Stress3', o3), (vmWord, ovm)]:
            cases[(subcaseID, i, 1, 'centroid', '%.3f')] = j
        return cases
