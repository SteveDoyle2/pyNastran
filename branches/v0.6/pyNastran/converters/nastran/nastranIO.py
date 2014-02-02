# pylint: disable=C0103,C0111,E1101

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

import os
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
                               ShellElement, LineElement, SpringElement)
from pyNastran.op2.test.test_op2 import OP2
from pyNastran.f06.f06 import F06


class NastranIO(object):
    def __init__(self):
        pass

    def load_nastran_geometry(self, bdf_filename, dirname):
        self.eidMap = {}
        self.nidMap = {}
        #print('bdf_filename=%r' % bdf_filename)
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        #skipReading = self.removeOldGeometry(bdfFileName)
        #if skipReading:
            #return
        if bdf_filename is None or bdf_filename is '':
            #self.grid = vtk.vtkUnstructuredGrid()
            #self.gridResult = vtk.vtkFloatArray()
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            #self.grid2 = vtk.vtkUnstructuredGrid()
            #self.scalarBar.VisibilityOff()
            return
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()
            #self.gridResult = vtk.vtkFloatArray()
            #self.gridResult.Reset()
            #self.gridResult.Modified()
            #self.eidMap = {}
            #self.nidMap = {}

            self.resultCases = {}
            self.nCases = 0
        for i in ('caseKeys', 'iCase', 'iSubcaseNameMap'):
            if hasattr(self, i):
                del i

            #print dir(self)
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        fname_base, ext = os.path.splitext(bdf_filename)
        punch = False
        if ext.lower() in '.pch':
            punch = True

        if bdf_filename.lower().endswith('.op2'):  # read the OP2
            op2_filename = bdf_filename
            model = OP2(op2_filename, make_geom=True,
                       save_skipped_cards=False,
                       debug=True, log=self.log)
            model.read_op2(op2_filename)
            model.crossReference(xref=True)
        else:  # read the bdf/punch
            model = BDF(log=self.log, debug=True)
            self.modelType = model.modelType
            model.readBDF(bdf_filename, includeDir=dirname, punch=punch)

        nNodes = model.nNodes()
        assert nNodes > 0
        nElements = model.nElements()
        nCAeros = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        print "nNodes = ",self.nNodes
        self.log_info("nElements = %i" % self.nElements)
        msg = model.card_stats(return_type='list')
        #self.log_info(msg)
        for msgi in msg:
            model.log.debug(msgi)

        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.card_count:
            nCONM2 = model.card_count['CONM2']
        else:
            nCONM2 = 0
        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(nCAeros + nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            i = 0
            fraction = 1. / nNodes  # so you can color the nodes by ID
            for (nid, node) in sorted(model.nodes.iteritems()):
                point = node.Position()
                points.InsertPoint(i, *point)
                #self.gridResult.InsertNextValue(i * fraction)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(),
                #                              elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

                self.nidMap[nid] = i
                i += 1

        # add the nodes
        for i, (nid, node) in enumerate(sorted(model.nodes.iteritems())):
            point = node.Position()
            points.InsertPoint(i, *point)
            self.nidMap[nid] = i

        # add the CAERO/CONM2 elements
        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros * 4 + nCONM2)
        for (eid, element) in sorted(model.caeros.iteritems()):
            if (isinstance(element, CAERO1) or isinstance(element, CAERO3) or
                isinstance(element, CAERO4) or isinstance(element, CAERO5)):
                cpoints = element.Points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j + 1)
                elem.GetPointIds().SetId(2, j + 2)
                elem.GetPointIds().SetId(3, j + 3)
                points2.InsertPoint(j, *cpoints[0])
                points2.InsertPoint(j + 1, *cpoints[1])
                points2.InsertPoint(j + 2, *cpoints[2])
                points2.InsertPoint(j + 3, *cpoints[3])
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 4
            #elif isinstance(element,CAERO2): # cylinder
                #pass
            else:
                self.log_info("skipping %s" % element.type)
        self.mapElements(points, points2, self.nidMap, model, j)

    def mapElements(self, points, points2, nidMap, model, j):
        #self.eidMap = {}
        i = 0
        for (eid, element) in sorted(model.elements.iteritems()):
            self.eidMap[eid] = i
            #print element.type
            if isinstance(element, CTRIA3) or isinstance(element, CTRIAR):
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTRIA6):
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTRIAX6):
                # midside nodes are required, nodes out of order
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
                #msg = "CTRIAX6 %i %i %i" %(nidMap[nodeIDs[a[0]]],
                #                           nidMap[nodeIDs[a[1]]],
                #                           nidMap[nodeIDs[a[2]]] )
                #raise RuntimeError(msg)
                #sys.stdout.flush()

                #elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                #elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                #elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

            elif (isinstance(element, CQUAD4) or isinstance(element, CSHEAR) or
                  isinstance(element, CQUADR)):
                nodeIDs = element.nodeIDs()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CQUAD8):
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTETRA10):
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

            elif isinstance(element, CPENTA15):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CHEXA8):
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CHEXA20):
                nodeIDs = element.nodeIDs()
                #print "nodeIDs = ",nodeIDs
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
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
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif (isinstance(element, LineElement) or
                  isinstance(element, SpringElement) or
                  element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                      'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                      'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC', ]):

                    nodeIDs = element.nodeIDs()
                    if 0 not in nodeIDs:
                        elem = vtk.vtkLine()
                        elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                        elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                        self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CONM2):  # not perfectly located
                del self.eidMap[eid]
                i -= 1

                #nid = element.Nid()
                c = element.Centroid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)

                points2.InsertPoint(j, *c)
                elem.GetPointIds().SetId(0, j)
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                del self.eidMap[eid]
                self.log_info("skipping %s" % element.type)
                continue
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

        cases = {}
        nelements = len(model.elements)
        print "len(elements) =", nelements
        pids = [] # zeros(nelements, 'int32')
        nxs = []
        nys = []
        nzs = []
        i = 0
        for eid, element in sorted(model.elements.iteritems()):
            pids.append(element.Pid())
            if isinstance(element, ShellElement):
                (nx, ny, nz) = element.Normal()
            else:
                nx = ny = nz = 0.0
            nxs.append(nx)
            nys.append(ny)
            nzs.append(nz)

        self.iSubcaseNameMap = {1: ['Nastran', '']}

        # subcaseID, resultType, vectorSize, location, dataFormat
        if 1:
                cases[(0, 'Pid', 1, 'centroid', '%.0f')] = pids

            # if not a flat plate???
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
                # subcaseID, resultType, vectorSize, location, dataFormat
                cases[(0, 'Normal_x', 1, 'centroid', '%.1f')] = nxs
                cases[(0, 'Normal_y', 1, 'centroid', '%.1f')] = nys
                cases[(0, 'Normal_z', 1, 'centroid', '%.1f')] = nzs
        self.log.info(cases.keys())
        self.finish_io(cases)

    def load_nastran_results(self, op2FileName, dirname):
        #gridResult.SetNumberOfComponents(self.nElements)
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        print("tring to read...", op2FileName)
        if '.op2' in op2FileName:
            model = OP2(op2FileName, log=self.log, debug=True)
            model.read_op2()
        elif '.f06' in op2FileName:
            model = F06(op2FileName, log=self.log, debug=True)
            model.read_f06()
        else:
            print("error...")
            raise NotImplementedError(op2FileName)

        #print model.print_results()

        #case = model.displacements[1]
        #print "case = ",case
        #for nodeID,translation in sorted(case.translations.iteritems()):
            #print "nodeID=%s t=%s" %(nodeID,translation)
        #self.iSubcaseNameMap[self.isubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = model.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = model.iSubcaseNameMap

        nElements = len(self.eidMap)
        print "nElements = ", nElements

        # set to True to enable nodeIDs as an result
        nidsSet = False
        subcaseID = subcaseIDs[0]
        if nidsSet:
            nids = zeros(self.nNodes, 'd')
            for (nid, nid2) in self.nidMap.iteritems():
                nids[nid2] = nid
            cases[(0, 'Node_ID', 1, 'node', '%.0f')] = nids
            nidsSet = True

        # set to True to enable elementIDs as a result
        eidsSet = True
        if eidsSet and self.is_centroidal:
            eids = zeros(nElements, 'd')
            for (eid, eid2) in self.eidMap.iteritems():
                eids[eid2] = eid

            eKey = (subcaseID, 'isElementOn', 1, 'centroid', '%.0g')
            cases[(subcaseID, 'Element_ID', 1, 'centroid', '%.0f')] = eids
            cases[eKey] = zeros(nElements)  # is the element supported
            eidsSet = True
        else:
            eKey = None

        print('is_nodal =', self.is_nodal)

        print('model.displacements.keys()', model.displacements.keys())
        for subcaseID in subcaseIDs:
            if self.is_nodal: # nodal results don't work
                if subcaseID in model.displacements:  # not correct?
                    case = model.displacements[subcaseID]
                    nnodes = self.nNodes
                    displacements = zeros((nnodes, 3), dtype='float32')
                    for (nid, txyz) in case.translations.iteritems():
                        nid2 = self.nidMap[nid]
                        displacements[nid2] = txyz
                    print("disp =", displacements)
                    key = (subcaseID, 'DisplacementX', 3, 'node', '%g')
                    print('adding case=%s' % str(key))
                    cases[key] = displacements

            if 0:
                if subcaseID in model.temperatures:
                    case = model.temperatures[subcaseID]
                    #print case
                    temps = zeros(self.nNodes, dtype='float32')
                    for (nid, T) in case.temperatures.iteritems():
                        #print T
                        nid2 = self.nidMap[nid]
                        temperatures[nid2] = T

                    key = (subcaseID, 'Temperature', 1, 'node', '%g')
                    #cases[key] = temperatures

            cases = self.fill_stress_case(cases, model, subcaseID, eKey)

        #self.resultCases = cases
        self.finish_io(cases)
        #return
        #self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        #self.iCase = -1
        #self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
        #self.cycleResults()  # start at nCase=0

    def finish_io(self, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print("caseKeys = ",self.caseKeys)
        #print "type(caseKeys) = ",type(self.caseKeys)

        if len(self.resultCases) == 0:
            self.nCases = 1
            self.iCase = 0
        elif len(self.resultCases) == 1:
            self.nCases = 1
            self.iCase = 0
        else:
            self.nCases = len(self.resultCases) - 1  # number of keys in dictionary
            self.iCase = -1
        self.cycleResults()  # start at nCase=0

        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

    def _finish_io(self, cases):
        #self.finish()
        ncases = len(cases)
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        #print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        #self.nCases = (ncases - 1) if (ncases - 1) > 0 else 0  # number of keys in dictionary
        self.nCases = ncases #if self.ncases == 0
        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

        self.iCase = 0 if self.nCases == 0 else -1
        self.cycleResults()  # start at nCase=0
        self.log.info('end of finish io')

    def fill_stress_case(self, cases, model, subcaseID, eKey):
        print("fill_stress_case")
        if self.is_centroidal:
            cases = self._fill_stress_case_centroidal(cases, model, subcaseID, eKey)
        else:
            #cases = self._fill_stress_case_nodal(cases, model, subcaseID, eKey)
            pass
        return cases

    def _fill_stress_case_nodal(self, cases, model, subcaseID, eKey):
        nnodes = self.nNodes

    def _fill_stress_case_centroidal(self, cases, model, subcaseID, eKey):
        nElements = self.nElements
        oxx = zeros(nElements)
        oyy = zeros(nElements)
        ozz = zeros(nElements)

        o1 = zeros(nElements)
        o2 = zeros(nElements)
        o3 = zeros(nElements)
        ovm = zeros(nElements)

        vmWord = 'N/A'
        if subcaseID in model.rodStress:
            case = model.rodStress[subcaseID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.
                #print "bar eid=%s" %(eid)

                axial = case.axial[eid]
                torsion = case.torsion[eid]

                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2] = max(axial, torsion)  # not really
                #oyy[eid2] = torsion
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = min(axial, torsion)

        if subcaseID in model.barStress:
            #self.s1    = {}
            #self.s2    = {}
            #self.s3    = {}
            #self.s4    = {}
            case = model.barStress[subcaseID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "bar eid=%s" %(eid)
                oxxi = case.axial[eid]
                o1i = max(case.smax[eid])
                o3i = min(case.smin[eid])

                oxx[eid2] = oxxi
                #oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                #ovm[eid2] = ovmi

        if subcaseID in model.plateStress:
            #self.txy    = {}
            case = model.plateStress[subcaseID]
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
                ozzi = min(case.oxx[eid]['C'], min(case.oyy[eid]['C']))

                o1i = max(case.majorP[eid]['C'])
                o2i = max(case.minorP[eid]['C'])
                o3i = min(case.majorP[eid]['C'], min(case.minorP[eid]['C']))
                ovmi = max(case.ovmShear[eid]['C'])

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if subcaseID in model.solidStress:
            case = model.solidStress[subcaseID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                oxxi = case.oxx[eid]['C']
                oyyi = case.oyy[eid]['C']
                ozzi = case.ozz[eid]['C']

                o1i = case.o1[eid]['C']
                o2i = case.o2[eid]['C']
                o3i = case.o3[eid]['C']
                ovmi = case.ovmShear[eid]['C']

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = o2i
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        # subcaseID,resultType,vectorSize,location,dataFormat
        cases[(subcaseID, 'StressXX', 1, 'centroid', '%.3f')] = oxx
        cases[(subcaseID, 'StressYY', 1, 'centroid', '%.3f')] = oyy
        cases[(subcaseID, 'StressZZ', 1, 'centroid', '%.3f')] = ozz

        cases[(subcaseID, 'Stress1', 1, 'centroid', '%.3f')] = o1
        cases[(subcaseID, 'Stress2', 1, 'centroid', '%.3f')] = o2
        cases[(subcaseID, 'Stress3', 1, 'centroid', '%.3f')] = o3
        cases[(subcaseID, vmWord, 1, 'centroid', '%.3f')] = ovm
        return cases
