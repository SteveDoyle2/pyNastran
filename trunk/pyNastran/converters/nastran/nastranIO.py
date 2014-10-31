# pylint: disable=C0103,C0111,E1101
from __future__ import print_function
from six import iteritems
from six.moves import zip

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
from numpy import zeros, abs, mean, where, nan_to_num, amax, amin, array
from numpy import nan as NaN
from numpy.linalg import norm

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
                               ShellElement, LineElement, SpringElement,
                               LOAD)
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

            #print(dir(self))
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        fname_base, ext = os.path.splitext(bdf_filename)
        punch = False
        if ext.lower() in '.pch':
            punch = True

        #if bdf_filename.lower().endswith('.op2'):  # read the OP2; requires make_geom
            #op2_filename = bdf_filename
            #model = OP2(# make_geom=True,
                       #debug=True, log=self.log)
            #model.read_op2(op2_filename)
            #model.cross_reference(xref=True, xref_loads=False, xref_constraints=False)
        if 0:
            pass
        else:  # read the bdf/punch
            model = BDF(log=self.log, debug=True)
            self.modelType = model.modelType
            model.read_bdf(bdf_filename, include_dir=dirname, punch=punch, xref=True)

        nNodes = model.nNodes()
        assert nNodes > 0
        nElements = model.nElements()
        assert nElements > 0
        nCAeros = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        self.log_info("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        msg = model.get_bdf_stats(return_type='list')
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
            for (nid, node) in sorted(iteritems(model.nodes)):
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
        node0 = model.nodes.keys()[0]
        position0 = model.nodes[node0].Position()
        xmin = position0[0]
        xmax = position0[0]

        ymin = position0[1]
        ymax = position0[1]

        zmin = position0[2]
        zmax = position0[2]
        for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
            point = node.Position()
            xmin = min(xmin, point[0])
            xmax = max(xmax, point[0])

            ymin = min(ymin, point[1])
            ymax = max(ymax, point[1])

            zmin = min(zmin, point[2])
            zmax = max(zmax, point[2])
            points.InsertPoint(i, *point)
            self.nidMap[nid] = i

        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.update_axes_length(dim_max)

        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s xmax=%s dz=%s" % (zmin, zmax, zmax-zmin))

        # add the CAERO/CONM2 elements
        j = 0
        points2 = vtk.vtkPoints()

        nsprings = 0
        if 0:
            for (eid, element) in sorted(iteritems(model.elements)):
                if (isinstance(element, LineElement) or
                      isinstance(element, SpringElement) or
                      element.type in ['CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                          'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                          'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC', ]):

                        nodeIDs = element.nodeIDs()
                        if None in nodeIDs:
                            nsprings += 1

        points2.SetNumberOfPoints(nCAeros * 4 + nCONM2 + nsprings)
        for (eid, element) in sorted(iteritems(model.caeros)):
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
                #elem.GetProperty().SetLineWidth(5)
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #self.grid2.GetProperty().SetLineWidth(5)
                j += 4
            #elif isinstance(element, CAERO2): # cylinder
                #pass
            else:
                self.log_info("skipping %s" % element.type)

        sphere_size = self._get_sphere_size(dim_max)
        for (eid, element) in sorted(iteritems(model.masses)):
            if isinstance(element, CONM2):
                #del self.eidMap[eid]

                #print("element", element)
                #print("element.nid", element.nid)
                #print('nodeIDs', model.nodes.keys())
                xyz = element.nid.Position()
                c = element.Centroid()
                d = norm(xyz-c)
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()

                if d == 0.:
                    d = sphere_size
                #elem.SetRadius(d)
                #elem.SetCenter(points.GetPoint(self.nidMap[nid]))
                #print(str(element))

                points2.InsertPoint(j, *c)
                elem.GetPointIds().SetId(0, j)
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                self.log_info("skipping %s" % element.type)

        self.mapElements(points, points2, self.nidMap, model, j, dim_max)

    def _get_sphere_size(self, dim_max):
        return 0.05 * dim_max

    def mapElements(self, points, points2, nidMap, model, j, dim_max):
        sphere_size = self._get_sphere_size(dim_max)
        #self.eidMap = {}

        # :param i: the element id in grid
        # :param j: the element id in grid2
        i = 0

        #nids = self.eid_to_nid_map[eid]
        self.eid_to_nid_map = {}

        # the list of all pids
        pids = []

        # pid = pids_dict[eid]
        pids_dict = {}
        nelements = len(model.elements)
        pids = zeros(nelements, 'int32')
        for (eid, element) in sorted(iteritems(model.elementss)):
            self.eidMap[eid] = i
            #print(element.type)
            pid = 0
            if isinstance(element, CTRIA3) or isinstance(element, CTRIAR):
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTRIA6):
                nodeIDs = element.nodeIDs()
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:3]
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
                pid = element.Pid()
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
                self.eid_to_nid_map[eid] = [nodeIDs[0], nodeIDs[2], nodeIDs[4]]
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
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CQUAD8):
                nodeIDs = element.nodeIDs()
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:4]
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
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:4]
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
            elif isinstance(element, CTETRA10):
                nodeIDs = element.nodeIDs()
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:4]
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
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:6]
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
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:6]
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
                pid = element.Pid()
                self.eid_to_nid_map[eid] = nodeIDs[:8]
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
                pid = element.Pid()
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

                self.eid_to_nid_map[eid] = nodeIDs[:8]
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

                # TODO: verify
                # CBUSH, CBUSH1D, CFAST, CROD, CELAS1, CELAS3
                # CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CVISC
                if hasattr(element, 'pid'):
                    pid = element.Pid()
                else:
                    # CONROD
                    # CELAS2, CELAS4?
                    pid = 0
                nodeIDs = element.nodeIDs()
                if nodeIDs[0] is None and  nodeIDs[0] is None: # CELAS2
                    del self.eidMap[eid]
                    continue
                if None in nodeIDs:  # used to be 0...
                    if nodeIDs[0] is None:
                        slot = 1
                    elif nodeIDs[1] is None:
                        slot = 0
                    #print('nodeIDs=%s slot=%s' % (str(nodeIDs), slot))
                    self.eid_to_nid_map[eid] = nodeIDs[slot]
                    nid = nodeIDs[slot]
                    if nid not in nidMap:
                        # SPOINT
                        continue
                    c = nidMap[nid]
                    elem = vtk.vtkVertex()
                    elem.GetPointIds().SetId(0, j)

                    elem = vtk.vtkSphere()
                    if d == 0.:
                        d = sphere_size
                    elem.SetRadius(d)
                else:
                    # 2 points
                    #d = norm(element.nodes[0].Position() - element.nodes[1].Position())
                    self.eid_to_nid_map[eid] = nodeIDs
                    elem = vtk.vtkLine()
                    try:
                        elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                        elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    except KeyError:
                        print("nodeIDs =", nodeIDs)
                        print(str(element))
                        continue

                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                del self.eidMap[eid]
                self.log_info("skipping %s" % element.type)
                continue
            # what about MPCs, RBE2s (rigid elements)?
            #   are they plotted as elements?
            #   and thus do they need a property?

            if pid is None:
                # CONROD
                #print(element)
                #pids[i] = 0
                #pids_dict[eid] = 0
                pass
            else:
                pids[i] = pid
                pids_dict[eid] = pid
            i += 1

        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print(dir(self.grid)) #.SetNumberOfComponents(0)
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        self.log_info("updated grid")

        cases = {}
        nelements = len(model.elements)
        pids = array(pids, 'int32')
        if not len(pids) == len(self.eidMap):
            msg = 'ERROR:  len(pids)=%s len(eidMap)=%s\n' % (len(pids), len(self.eidMap))
            for eid, pid in sorted(iteritems(pids_dict)):
                if eid not in self.eidMap:
                    msg += 'eid=%s %s' % (eid, str(model.elements[eid]))
            raise RuntimeError(msg)
        del pids_dict

        #nxs = []
        #nys = []
        #nzs = []
        #i = 0

        #for eid, element in sorted(iteritems(model.elements)):
            #if isinstance(element, ShellElement):
                #(nx, ny, nz) = element.Normal()
            #else:
                #nx = ny = nz = 0.0
            #nxs.append(nx)
            #nys.append(ny)
            #nzs.append(nz)

        self.iSubcaseNameMap = {1: ['Nastran', '']}

        nElements = len(self.eidMap)
        #print("nElements = ", nElements)

        # set to True to enable nodeIDs as an result
        nidsSet = True
        if nidsSet and self.is_nodal:
            nids = zeros(self.nNodes, 'd')
            for (nid, nid2) in iteritems(self.nidMap):
                nids[nid2] = nid
            cases[(0, 'Node_ID', 1, 'node', '%.0f')] = nids
            nidsSet = True

        # set to True to enable elementIDs as a result
        eidsSet = True
        if eidsSet and self.is_centroidal:
            eids = zeros(nElements, dtype='int32')
            for (eid, eid2) in iteritems(self.eidMap):
                eids[eid2] = eid
            cases[(0, 'Element_ID', 1, 'centroid', '%.0f')] = eids
            eidsSet = True

        # subcaseID, resultType, vectorSize, location, dataFormat
        if len(model.properties) and self.is_centroidal:
            cases[(0, 'Pid', 1, 'centroid', '%.0f')] = pids

        #self._plot_pressures(model, cases)
        #self._plot_applied_loads(model, cases)

        if 0:
            # if not a flat plate???
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
                # subcaseID, resultType, vectorSize, location, dataFormat
                cases[(0, 'Normal_x', 1, 'centroid', '%.1f')] = nxs
                cases[(0, 'Normal_y', 1, 'centroid', '%.1f')] = nys
                cases[(0, 'Normal_z', 1, 'centroid', '%.1f')] = nzs
        self.log.info(cases.keys())
        #self.finish_io(cases)

        self.finish_nastran_io(cases)

    def finish_nastran_io(self, cases):  # same as Cart3d version
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print("ncases =", len(cases))
        print("caseKeys =", self.caseKeys)

        if len(self.caseKeys) > 1:
            print("finish_io case A")
            self.iCase = -1
            self.nCases = len(self.resultCases)  # number of keys in dictionary
        elif len(self.caseKeys) == 1:
            print("finish_io case B")
            self.iCase = -1
            self.nCases = 1
        else:
            print("finish_io case C")
            self.iCase = -1
            self.nCases = 0

        self.cycleResults_explicit()  # start at nCase=0
        if self.nCases:
            self.scalarBar.VisibilityOn()
            self.scalarBar.Modified()

    def _plot_pressures(self, model, cases):
        """
        does pressure act normal or antinormal?
        """
        if not self.is_centroidal:
            return
        sucaseIDs = model.caseControlDeck.get_subcase_list()
        for subcaseID in sucaseIDs:
            if subcaseID == 0:
                continue
            load_case_id, options = model.caseControlDeck.get_subcase_parameter(subcaseID, 'LOAD')
            loadCase = model.loads[load_case_id]

            # account for scale factors
            loads2 = []
            scale_factors2 = []
            for load in loadCase:
                if isinstance(load, LOAD):
                    scale_factors, loads = load.getReducedLoads()
                    scale_factors2 += scale_factors
                    loads2 += loads
                else:
                    scale_factors2.append(1.)
                    loads2.append(load)

            eids = sorted(model.elements.keys())
            pressures = zeros(len(model.elements), dtype='float32')

            # loop thru scaled loads and plot the pressure
            for load, scale in zip(loads2, scale_factors2):
                if load.type == 'PLOAD4':
                    elem = load.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                     'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        p = load.pressures[0] * scale

                        # single element per PLOAD
                        #eid = elem.eid
                        #pressures[eids.index(eid)] = p

                        # multiple elements
                        for el in load.eids:
                            pressures[eids.index(el.eid)] += p
                    #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        #A, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                        #r = centroid - p

            # if there is no applied pressure, don't make a plot
            if abs(pressures).max():
                # subcaseID, resultType, vectorSize, location, dataFormat
                cases[(0, 'Pressure Case=%i' % subcaseID, 1, 'centroid', '%.1f')] = pressures

    def _plot_applied_loads(self, model, cases):
        if not self.is_nodal:
            return
        sucaseIDs = model.caseControlDeck.get_subcase_list()
        for subcaseID in sucaseIDs:
            if subcaseID == 0:
                continue
            load_case_id, options = model.caseControlDeck.get_subcase_parameter(subcaseID, 'LOAD')
            loadCase = model.loads[load_case_id]

            # account for scale factors
            loads2 = []
            scale_factors2 = []
            for load in loadCase:
                if isinstance(load, LOAD):
                    scale_factors, loads = load.getReducedLoads()
                    scale_factors2 += scale_factors
                    loads2 += loads
                else:
                    scale_factors2.append(1.)
                    loads2.append(load)

            eids = sorted(model.elements.keys())
            nids = sorted(model.nodes.keys())
            loads = zeros((self.nNodes, 3), dtype='float32')

            # loop thru scaled loads and plot the pressure
            for load, scale in zip(loads2, scale_factors2):
                if load.type == 'FORCE':
                    scale2 = load.mag * scale  # does this need a magnitude?
                    nid = load.node
                    loads[nids.index(nid)] += load.xyz * scale2
                elif load.type == 'PLOAD4':  # centrodial, skipping
                    continue
                    elem = load.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        eid = elem.eid
                        node_ids = elem.nodeIDs()
                        k = load.pressures[0] * scale / 3.
                        for nid in node_ids[3:]:
                            pressures[eids.index(nid)] += k
                    if elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        eid = elem.eid
                        node_ids = elem.nodeIDs()
                        k = load.pressures[0] * scale / 4.
                        for nid in node_ids[4:]:
                            pressures[eids.index(nid)] += k
            if loads[:, 0].min() != loads[:, 0].max():
                cases[(subcaseID, 'LoadX Case=%i' % subcaseID, 1, 'node', '%.1f')] = loads[:,0]
            if loads[:, 1].min() != loads[:, 1].max():
                cases[(subcaseID, 'LoadY Case=%i' % subcaseID, 1, 'node', '%.1f')] = loads[:,1]
            if loads[:, 2].min() != loads[:, 2].max():
                cases[(subcaseID, 'LoadZ Case=%i' % subcaseID, 1, 'node', '%.1f')] = loads[:,2]

    def load_nastran_results(self, op2_filename, dirname):
        #gridResult.SetNumberOfComponents(self.nElements)
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        print("tring to read...", op2_filename)
        if '.op2' in op2_filename:  # TODO: do this based on lower & file extension
            model = OP2(log=self.log, debug=True)
            model._saved_results = set([])
            all_results = model.get_all_results()
            if self.is_nodal:
                desired_results = [
                    'displacements', 'velocities', 'accelerations', 'temperatures',
                    'constraint_forces', 'spcForces', 'mpcForces',

                    #'gridPointForces',
                    'stress', 'solidStress', 'plateStress', 'compositePlateStress', 'barStress', 'rodStress',
                    #'strain','solidStrain', 'plateStrain', 'compositePlateStrain', 'barStrain', 'rodStrain',

                    # untested
                    'loadVectors',
                    'appliedLoads',
                    'forceVectors',
                ]
            else:
                desired_results = [
                    'stress', 'solidStress', 'plateStress', 'compositePlateStress', 'barStress', 'rodStress',
                    #'strain','solidStrain', 'plateStrain', 'compositePlateStrain', 'barStrain', 'rodStrain',
                ]
            for result in desired_results:
                if result in all_results:
                    model._saved_results.add(result)
            model.read_op2(op2_filename)

        elif '.f06' in op2_filename:  # TODO: do this based on lower & file extension
            model = F06(log=self.log, debug=True)
            model.read_f06(op2_filename)
        else:
            print("error...")
            raise NotImplementedError(op2_filename)

        #print(model.print_results())

        #case = model.displacements[1]
        #print("case = ",case)
        #for nodeID,translation in sorted(iteritems(case.translations)):
            #print("nodeID=%s t=%s" %(nodeID,translation))
        #self.iSubcaseNameMap[self.isubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = model.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = model.iSubcaseNameMap

        for subcaseID in subcaseIDs:
            cases = self.fill_oug_oqg_case(cases, model, subcaseID)
            cases = self.fill_stress_case(cases, model, subcaseID)

        self.finish_nastran_io(cases)

    def fill_oug_oqg_case(self, cases, model, subcaseID):
        if self.is_nodal: # nodal results don't work with centroidal ones
            displacement_like = [
                [model.displacements, 'Displacement'],
                [model.velocities,    'Velocity'],
                [model.accelerations, 'Acceleration'],
                [model.spcForces,     'SPC Forces'],
                [model.mpcForces,     'MPC Forces'],

                # untested
                [model.loadVectors, 'loadVectors'],
                [model.appliedLoads, 'appliedLoads'],
                [model.forceVectors, 'forceVectors'],
                #[model.gridPointForces, 'GridPointForces'],  # TODO: this is buggy...
            ]
            temperature_like = [
                [model.temperatures,    'Temperature']
            ]
            nnodes = self.nNodes

            # size = 3
            for (result, name) in displacement_like:
                if subcaseID in result:
                    case = result[subcaseID]

                    if case.nonlinear_factor is not None: # transient
                        return
                    displacements = zeros((nnodes, 3), dtype='float32')
                    x_displacements = zeros(nnodes, dtype='float32')
                    y_displacements = zeros(nnodes, dtype='float32')
                    z_displacements = zeros(nnodes, dtype='float32')
                    xyz_displacements = zeros(nnodes, dtype='float32')

                    if hasattr(case, 'translations'):
                        word = 'translations'
                    elif hasattr(case, 'forces'):
                        word = 'forces'
                    else:
                        print('ERROR!!!!', case.__dict__.keys())
                        raise RuntimeError(case.__dict__.keys())

                    res = getattr(case, word)

                    print('case.type =', case.__class__.__name__)
                    for (nid, txyz) in iteritems(res):
                        nid2 = self.nidMap[nid]
                        displacements[nid2] = txyz
                        xyz_displacements[nid2] = norm(txyz)
                        x_displacements[nid2] = txyz[0]
                        y_displacements[nid2] = txyz[1]
                        z_displacements[nid2] = txyz[2]

                    #cases[(subcaseID, name + 'Vector', 3, 'node', '%g')] = displacements
                    cases[(subcaseID, name + 'X', 1, 'node', '%g')] = x_displacements
                    cases[(subcaseID, name + 'Y', 1, 'node', '%g')] = y_displacements
                    cases[(subcaseID, name + 'Z', 1, 'node', '%g')] = z_displacements
                    cases[(subcaseID, name + 'XYZ', 1, 'node', '%g')] = xyz_displacements
                    del res, word

            # size = 1
            for (result, name) in temperature_like:
                if subcaseID in result:
                    case = result[subcaseID]
                    if case.nonlinear_factor is not None: # transient
                        return
                    temperatures = zeros(nnodes, dtype='float32')
                    for (nid, txyz) in iteritems(case.translations):
                        nid2 = self.nidMap[nid]
                        displacements[nid2] = txyz
                        temperatures[nid2] = norm(txyz)

                        #cases[(subcaseID, name + 'Vector', 3, 'node', '%g')] = displacements
                        cases[(subcaseID, name, 1, 'node', '%g')] = temperatures
        return cases

    def finish_io(self, cases):
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print("caseKeys = ",self.caseKeys)

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

    def clear_nastran(self):
        self.eidMap = {}
        self.nidMap = {}
        self.eid_to_nid_map = {}

    def fill_stress_case(self, cases, model, subcaseID):
        return cases
        if self.is_centroidal:
            self._fill_stress_case_centroidal(cases, model, subcaseID)
        elif self.is_nodal:
            self._fill_stress_case_nodal(cases, model, subcaseID)
        else:
            raise RuntimeError('this shouldnt happen...')
        return cases

    def _fill_stress_case_nodal(self, cases, model, subcaseID):
        oxx_dict = {}
        oyy_dict = {}
        ozz_dict = {}
        o1_dict = {}
        o2_dict = {}
        o3_dict = {}
        ovm_dict = {}

        for nid in self.nidMap:
            oxx_dict[nid] = []
            oyy_dict[nid] = []
            ozz_dict[nid] = []
            o1_dict[nid] = []
            o2_dict[nid] = []
            o3_dict[nid] = []
            ovm_dict[nid] = []

        vmWord = None
        if subcaseID in model.rodStress:
            case = model.rodStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.axial:
                axial = case.axial[eid]
                torsion = case.torsion[eid]
                node_ids = self.eid_to_nid_map[eid]
                o1i = max(axial, torsion)  # not really
                o3i = min(axial, torsion)
                ovmi = max(abs(axial), abs(torsion))
                for nid in node_ids:
                    oxx_dict[nid].append(axial)
                    oyy_dict[nid].append(torsion)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcaseID in model.barStress:
            case = model.barStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.axial:
                node_ids = self.eid_to_nid_map[eid]
                oxxi = case.axial[eid]
                o1i = max(case.smax[eid])
                o3i = min(case.smin[eid])
                ovmi = max(    max(case.smax[eid]),
                           abs(min(case.smin[eid])))
                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcaseID in model.beamStress:
            case = model.beamStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.smax:
                node_ids = self.eid_to_nid_map[eid]
                oxxi = max(max(case.sxc[eid]),
                           max(case.sxd[eid]),
                           max(case.sxe[eid]),
                           max(case.sxf[eid]))
                o1i = max(case.smax[eid])
                o3i =min(case.smin[eid])
                ovmi =max(    max(case.smax[eid]),
                          abs(min(case.smin[eid])))
                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcaseID in model.plateStress:
            case = model.plateStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]

                eType = case.eType[eid]
                if eType in ['CQUAD4', 'CQUAD8']:
                    #cen = 'CEN/%s' % eType[-1]
                    for nid in node_ids:
                        oxxi = max(case.oxx[eid][nid])
                        oyyi = max(case.oyy[eid][nid])
                        ozzi = min(case.oxx[eid][nid], min(case.oyy[eid][nid]))
                        o1i = max(case.majorP[eid][nid])
                        o2i = max(case.minorP[eid][nid])
                        o3i = min(case.majorP[eid][nid], min(case.minorP[eid][nid]))
                        ovmi = max(case.ovmShear[eid][nid])

                        oxx_dict[nid].append(oxxi)
                        oyy_dict[nid].append(oyyi)
                        o1_dict[nid].append(o1i)
                        o3_dict[nid].append(o3i)
                        ovm_dict[nid].append(ovmi)

                elif eType in ['CTRIA3', 'CTRIA6']:
                    cen = 'CEN/%s' % eType[-1]
                    oxxi = case.oxx[eid][cen]

                    oxxi = max(case.oxx[eid][cen])
                    oyyi = max(case.oyy[eid][cen])
                    ozzi = min(case.oxx[eid][cen], min(case.oyy[eid][cen]))

                    o1i = max(case.majorP[eid][cen])
                    o2i = max(case.minorP[eid][cen])
                    o3i = min(case.majorP[eid][cen], min(case.minorP[eid][cen]))
                    ovmi = max(case.ovmShear[eid][cen])

                    for nid in node_ids:
                        oxx_dict[nid].append(oxxi)
                        oyy_dict[nid].append(oyyi)
                        o1_dict[nid].append(o1i)
                        o3_dict[nid].append(o3i)
                        ovm_dict[nid].append(ovmi)

        if subcaseID in model.compositePlateStress:
            case = model.compositePlateStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'

            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]

                oxxi = max(case.o11[eid])
                oyyi = max(case.o22[eid])
                o1i = max(case.majorP[eid])
                o3i = min(case.minorP[eid])
                ovmi = max(case.ovmShear[eid])

                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    oyy_dict[nid].append(oyyi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcaseID in model.solidStress:
            case = model.solidStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                node_ids = self.eid_to_nid_map[eid]
                for nid in node_ids:
                    oxxi = case.oxx[eid][nid]
                    oyyi = case.oyy[eid][nid]
                    ozzi = case.ozz[eid][nid]
                    o1i = case.o1[eid][nid]
                    o2i = case.o2[eid][nid]
                    o3i = case.o3[eid][nid]
                    ovmi = case.ovmShear[eid][nid]

                    oxx_dict[nid].append(oxxi)
                    oyy_dict[nid].append(oyyi)
                    ozz_dict[nid].append(ozzi)
                    o1_dict[nid].append(o1i)
                    o2_dict[nid].append(o2i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        nnodes = self.nNodes
        oxx = zeros(nnodes, dtype='float32')
        oyy = zeros(nnodes, dtype='float32')
        ozz = zeros(nnodes, dtype='float32')
        o1 = zeros(nnodes, dtype='float32')
        o2 = zeros(nnodes, dtype='float32')
        o3 = zeros(nnodes, dtype='float32')
        ovm = zeros(nnodes, dtype='float32')
        for i, nid in enumerate(sorted(self.nidMap)):
            oxx[i] = mean(oxx_dict[nid])
            oyy[i] = mean(oyy_dict[nid])
            ozz[i] = mean(ozz_dict[nid])
            o1[i] = mean(o1_dict[nid])
            o2[i] = mean(o2_dict[nid])
            o3[i] = mean(o3_dict[nid])
            ovm[i] = mean(ovm_dict[nid])

        # do this to prevent screwy stresses at points that have no stress
        oxx = nan_to_num(oxx)
        oyy = nan_to_num(oyy)
        ozz = nan_to_num(ozz)
        o1 = nan_to_num(o1)
        o2 = nan_to_num(o2)
        o3 = nan_to_num(o3)
        ovm = nan_to_num(ovm)
        if oxx.min() != oxx.max():
            cases[(subcaseID, 'StressXX', 1, 'node', '%.3f')] = oxx
        if oyy.min() != oyy.max():
            cases[(subcaseID, 'StressYY', 1, 'node', '%.3f')] = oyy
        if ozz.min() != ozz.max():
            cases[(subcaseID, 'StressZZ', 1, 'node', '%.3f')] = ozz

        if o1.min() != o1.max():
            cases[(subcaseID, 'Stress1', 1, 'node', '%.3f')] = o1
        if o2.min() != o2.max():
            cases[(subcaseID, 'Stress2', 1, 'node', '%.3f')] = o2
        if o3.min() != o3.max():
            cases[(subcaseID, 'Stress3', 1, 'node', '%.3f')] = o3
        if vmWord is not None:
            cases[(subcaseID, vmWord, 1, 'node', '%.3f')] = ovm

    def _fill_stress_case_centroidal(self, cases, model, subcaseID):
        nElements = self.nElements
        isElementOn = zeros(nElements)  # is the element supported

        oxx = zeros(nElements, dtype='float32')
        oyy = zeros(nElements, dtype='float32')
        ozz = zeros(nElements, dtype='float32')

        o1 = zeros(nElements, dtype='float32')
        o2 = zeros(nElements, dtype='float32')
        o3 = zeros(nElements, dtype='float32')
        ovm = zeros(nElements, dtype='float32')

        vmWord = None
        if subcaseID in model.rodStress:
            case = model.rodStress[subcaseID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                axial = case.axial[eid]
                torsion = case.torsion[eid]

                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2] = max(axial, torsion)  # not really
                o3[eid2] = min(axial, torsion)
                ovm[eid2] = max(abs(axial), abs(torsion))

        if subcaseID in model.barStress:
            case = model.barStress[subcaseID]
            if case.nonlinear_factor is not None:
                return
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                oxx[eid2] = case.axial[eid]

                o1[eid2] = max(case.smax[eid])
                o3[eid2] = min(case.smin[eid])
                ovm[eid2] = max(    max(case.smax[eid]),
                                abs(min(case.smin[eid])))

        if subcaseID in model.beamStress:
            case = model.beamStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.smax:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                oxx[eid2] = max(max(case.sxc[eid]),
                                max(case.sxd[eid]),
                                max(case.sxe[eid]),
                                max(case.sxf[eid]))

                o1[eid2] = max(case.smax[eid])
                o3[eid2] = min(case.smin[eid])
                ovm[eid2] = max(    max(case.smax[eid]),
                                abs(min(case.smin[eid])))

        if subcaseID in model.plateStress:
            case = model.plateStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                eType = case.eType[eid]
                if eType in ['CQUAD4', 'CQUAD8', 'CTRIA3', 'CTRIA6']:
                    cen = 'CEN/%s' % eType[-1]
                oxxi = case.oxx[eid][cen]

                oxxi = max(case.oxx[eid][cen])
                oyyi = max(case.oyy[eid][cen])
                ozzi = min(case.oxx[eid][cen], min(case.oyy[eid][cen]))

                o1i = max(case.majorP[eid][cen])
                o2i = max(case.minorP[eid][cen])
                o3i = min(case.majorP[eid][cen], min(case.minorP[eid][cen]))
                ovmi = max(case.ovmShear[eid][cen])

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi

                o1[eid2] = o1i
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if subcaseID in model.compositePlateStress:
            case = model.compositePlateStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                eType = case.eType[eid]
                oxxi = max(case.o11[eid])
                oyyi = max(case.o22[eid])

                o1i = max(case.majorP[eid])
                o3i = min(case.minorP[eid])
                ovmi = max(case.ovmShear[eid])

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi

                o1[eid2] = o1i
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if subcaseID in model.solidStress:
            case = model.solidStress[subcaseID]
            if case.nonlinear_factor is not None: # transient
            #if case.isTransient():
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                isElementOn[eid2] = 1.

                #print('case.oxx[%i].keys() = %s' % (eid, case.oxx[eid].keys()))
                oxx[eid2] = case.oxx[eid]['CENTER']
                oyy[eid2] = case.oyy[eid]['CENTER']
                ozz[eid2] = case.ozz[eid]['CENTER']

                o1[eid2] = case.o1[eid]['CENTER']
                o2[eid2] = case.o2[eid]['CENTER']
                o3[eid2] = case.o3[eid]['CENTER']
                ovm[eid2] = case.ovmShear[eid]['CENTER']

        # subcaseID,resultType,vectorSize,location,dataFormat
        if isElementOn.min() != isElementOn.max():
            cases[(1, 'isElementOn', 1, 'centroid', '%.0f')] = isElementOn
        if oxx.min() != oxx.max():
            cases[(subcaseID, 'StressXX', 1, 'centroid', '%.3f')] = oxx
        if oyy.min() != oyy.max():
            cases[(subcaseID, 'StressYY', 1, 'centroid', '%.3f')] = oyy
        if ozz.min() != ozz.max():
            cases[(subcaseID, 'StressZZ', 1, 'centroid', '%.3f')] = ozz

        if o1.min() != o1.max():
            cases[(subcaseID, 'Stress1', 1, 'centroid', '%.3f')] = o1
        if o2.min() != o2.max():
            cases[(subcaseID, 'Stress2', 1, 'centroid', '%.3f')] = o2
        if o3.min() != o3.max():
            cases[(subcaseID, 'Stress3', 1, 'centroid', '%.3f')] = o3
        if vmWord is not None:
            cases[(subcaseID, vmWord, 1, 'centroid', '%.3f')] = ovm
        return cases

def main():
    from pyNastran.gui.testing_methods import add_dummy_gui_functions
    test = NastranIO()
    test.is_nodal = False
    test.is_centroidal = True

    add_dummy_gui_functions(test)

    #test.load_panair_geometry('SWB.INP','')
    test.load_nastran_geometry('bottle_shell_w_holes_pmc.bdf', '')
    test.load_nastran_results('bottle_shell_w_holes_pmc.op2', '')

    keys = test.resultCases.keys()
    assert (1, 'Stress1', 1, 'centroid', '%.3f') in keys, keys

if __name__ == '__main__':
    main()