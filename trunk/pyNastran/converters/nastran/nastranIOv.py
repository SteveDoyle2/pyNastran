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
from numpy import nan as NaN, searchsorted, sqrt, pi, arange, unique
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
from pyNastran.op2.op2_vectorized import OP2_Vectorized
from pyNastran.f06.f06 import F06


class NastranIO(object):
    def __init__(self):
        #: flips the nastran CAERO subpaneling
        #:   False -> borders of CAEROs can be seen
        #:   True  -> individual subpanels can be seen
        self.is_sub_panels = False

        #: coordinate systems can be messy, so this is the
        #: list of coords to show
        self.show_cids = []
        self.save_data = False
        self.show_alt_actor = True  # show the caero mesh

    def show_caero_mesh(self, is_shown=None):
        """
        :param is_shown: should the mesh be shown/hidden (default=None -> flip between shown/not shown)
        """
        if is_shown is None:
            is_shown = not(self.show_alt_actor)

        self.show_alt_actor = is_shown
        if is_shown:
            if not self.show_alt_actor:
                return
            self.alt_geometry_actor.VisibilityOn()
        else:
            if self.show_alt_actor:
                return
            self.alt_geometry_actor.VisibilityOff()


    def load_nastran_geometry(self, bdf_filename, dirname, plot=True):
        #print('plot=%s' % plot)
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
        if plot:
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
            model.read_bdf(bdf_filename, include_dir=dirname, punch=punch, xref=False)
            model.cross_reference(xref=True, xref_loads=False,
                                  xref_constraints=False)

        nNodes = model.nNodes()
        assert nNodes > 0
        nElements = model.nElements()
        assert nElements > 0

        if self.is_sub_panels:
            nsub_elements_caeros = 0
            nsub_points_caeros = 0
            for key, caero in iteritems(model.caeros):
                if hasattr(caero, 'panel_points_elements'):
                    npoints, nelements = caero.get_npanel_points_elements()
                    nsub_elements_caeros += npoints
                    nsub_points_caeros += nelements
                else:
                    print('%r doesnt support panel_points_elements' % caero.type)
            nCAeros = nsub_elements_caeros
            nCAerosPoints = nsub_points_caeros
        else:
            nCAeros = model.nCAeros()
            nCAerosPoints = nCAeros * 4

        self.nNodes = nNodes
        self.nElements = nElements  # approximate...

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
        self.show_caero_mesh(is_shown=self.show_alt_actor)

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
        if self.save_data:
            self.model = model
            n = len(model.nodes)
            xyz_cid0 = zeros((n, 3), dtype='float64')
            for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
                xyz = node.Position()
                xyz_cid0[i, :] = xyz
            self.xyz_cid0 = xyz_cid0

        cid_types = {
            'R' : 'xyz',
            'C' : 'Rtz',
            'S' : 'Rtp',
        }
        for cid, coord in sorted(iteritems(model.coords)):
            if cid == 0:
                continue
            if cid in self.show_cids:
                origin = coord.origin
                beta = coord.beta()
                Type = cid_types[coord.Type]
                self.create_coordinate_system(label=cid, origin=origin, matrix_3x3=beta, Type=Type)
            else:
                print('skipping cid=%s; use a script and set self.show_cids=[%s] to view' % (cid, cid))

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

        points2.SetNumberOfPoints(nCAerosPoints * 4 + nCONM2 + nsprings)
        for (eid, element) in sorted(iteritems(model.caeros)):
            if (isinstance(element, CAERO1) or isinstance(element, CAERO3) or
                isinstance(element, CAERO4) or isinstance(element, CAERO5)):
                if self.is_sub_panels:
                    pointsi, elementsi = element.panel_points_elements()
                    for ipoint, pointii in enumerate(pointsi):
                        points2.InsertPoint(j + ipoint, *pointii)

                    elem = vtkQuad()
                    eType = elem.GetCellType()
                    for ielement, elementsi in enumerate(elementsi):
                        elem = vtkQuad()
                        elem.GetPointIds().SetId(0, j + elementsi[0])
                        elem.GetPointIds().SetId(1, j + elementsi[1])
                        elem.GetPointIds().SetId(2, j + elementsi[2])
                        elem.GetPointIds().SetId(3, j + elementsi[3])
                        self.grid2.InsertNextCell(eType, elem.GetPointIds())
                    j += ipoint + 1
                    #isubpanel += ipoint
                else:
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

        self.mapElements(points, points2, self.nidMap, model, j, dim_max, plot=plot)

    def _get_sphere_size(self, dim_max):
        return 0.05 * dim_max

    def mapElements(self, points, points2, nidMap, model, j, dim_max, plot=True):
        sphere_size = self._get_sphere_size(dim_max)
        #self.eidMap = {}

        # :param i: the element id in grid
        # :param j: the element id in grid2
        i = 0

        #nids = self.eid_to_nid_map[eid]
        self.eid_to_nid_map = {}

        # the list of all pids
        #pids = []

        # pid = pids_dict[eid]
        pids_dict = {}
        nelements = len(model.elements)
        pids = zeros(nelements, 'int32')
        for (eid, element) in sorted(iteritems(model.elements)):
            self.eidMap[eid] = i
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
                    print('removing CELASx eid=%i -> no node %i' % (eid, nodeIDs[0]))
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
                        print('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
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
                print('removing eid=%s' % eid)
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

        nelements = i
        self.nElements = nelements
        #print('nelements=%s pids=%s' % (nelements, list(pids)))
        pids = pids[:nelements]
        #print('len(pids) = ', len(pids))
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
        #pids = array(pids, 'int32')
        if 0:
            if not len(pids) == len(self.eidMap):
                msg = 'ERROR:  len(pids)=%s len(eidMap)=%s\n' % (len(pids), len(self.eidMap))
                for eid, pid in sorted(iteritems(pids_dict)):
                    #self.eidMap[eid] = i
                    #pids_dict[eid] = pid
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

        #nelements = len(self.eidMap)
        icase = 0
        form = ['Geometry', None, []]
        form0 = form[2]

        # set to True to enable nodeIDs as an result
        nidsSet = True
        if nidsSet and self.is_nodal:
            nids = zeros(self.nNodes, dtype='int32')
            for (nid, nid2) in iteritems(self.nidMap):
                nids[nid2] = nid
            cases[(0, icase, 'Node_ID', 1, 'node', '%i')] = nids
            form0.append(('Node_ID', icase, []))
            icase += 1
            self.node_ids = nids
            nidsSet = True

        # set to True to enable elementIDs as a result
        eidsSet = True
        if eidsSet and self.is_centroidal:
            eids = zeros(nelements, dtype='int32')
            for (eid, eid2) in iteritems(self.eidMap):
                eids[eid2] = eid
            cases[(0, icase, 'Element_ID', 1, 'centroid', '%i')] = eids
            form0.append(('Element_ID', icase, []))
            icase += 1
            self.element_ids = eids
            eidsSet = True

        # subcaseID, resultType, vectorSize, location, dataFormat
        if len(model.properties) and self.is_centroidal:
            cases[(0, icase, 'Property_ID', 1, 'centroid', '%i')] = pids
            form0.append(('Property_ID', icase, []))
            icase += 1

        #self._plot_pressures(model, cases)
        #self._plot_applied_loads(model, cases)

        if 0:
            # if not a flat plate???
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
                # subcaseID, resultType, vectorSize, location, dataFormat
                cases[(0, icase, 'Normal_x', 1, 'centroid', '%.1f')] = nxs
                form0.append(('Normal_x', icase, []))
                icase += 1

                cases[(0, icase, 'Normal_y', 1, 'centroid', '%.1f')] = nys
                form0.append(('Normal_y', icase, []))
                icase += 1

                cases[(0, icase, 'Normal_z', 1, 'centroid', '%.1f')] = nzs
                form0.append(('Normal_z', icase, []))
                icase += 1

        self.log.info(cases.keys())
        #self.finish_io(cases)

        if plot:
            #self._finish_results_io(cases)
            self._finish_results_io2([form], cases)

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
        #self.show_caero_mesh()

        print("tring to read...", op2_filename)
        if '.op2' in op2_filename:  # TODO: do this based on lower & file extension
            model = OP2_Vectorized(log=self.log, debug=True)

            if 0:
                model._saved_results = set([])
                all_results = model.get_all_results()
                if self.is_nodal:
                    desired_results = [
                        'displacements', 'velocities', 'accelerations', 'temperatures',
                        'constraint_forces', 'spcForces', 'mpcForces', 'eigenvectors',

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

        self.log.info(model.get_op2_stats())
        print(model.get_op2_stats())

        #print(model.print_results())

        #self.iSubcaseNameMap[self.isubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = model.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = model.iSubcaseNameMap

        form = []
        icase = 0
        for subcaseID in subcaseIDs:
            subcase_name = 'Subcase %i' % subcaseID
            form0 = (subcase_name, None, [])
            formi = form0[2]
            icase = self.fill_oug_oqg_case(cases, model, subcaseID, formi, icase)
            icase = self.fill_stress_case(cases, model, subcaseID, formi, icase)
            if len(formi):
                form.append(form0)

        self._finish_results_io2(form, cases)

    def fill_oug_oqg_case(self, cases, model, subcaseID, formi, icase):
        nnodes = self.nNodes
        if self.is_nodal: # nodal results don't work with centroidal ones
            displacement_like = [
                (model.displacements, 'Displacement'),
                (model.velocities,    'Velocity'),
                (model.accelerations, 'Acceleration'),
                (model.eigenvectors,  'Eigenvectors'),
                (model.spcForces,     'SPC Forces'),
                (model.mpcForces,     'MPC Forces'),

                # untested
                (model.loadVectors, 'loadVectors'),
                (model.appliedLoads, 'appliedLoads'),
                (model.forceVectors, 'forceVectors'),
                #[model.gridPointForces, 'GridPointForces'],  # TODO: this is buggy...
            ]
            temperature_like = [
                (model.temperatures, 'Temperature'),
            ]
            nids = self.node_ids

            for (result, name) in displacement_like:
                if subcaseID in result:
                    case = result[subcaseID]
                    if not case.is_real():
                        return icase
                    if case.nonlinear_factor is not None: # transient
                        ntimes = len(case._times)
                        code_name = case.data_code['name']
                        has_cycle = hasattr(case, 'mode_cycle')

                        itime0 = 0
                        t1 = case.data[itime0, :, 0]
                        ndata, = t1.shape
                        if nnodes != ndata:
                            nidsi = case.node_gridtype[:, 0]
                            j = searchsorted(nids, nidsi)

                        for itime in range(case.ntimes):
                            dt = case._times[itime]
                            t1 = case.data[itime, :, 0]
                            t2 = case.data[itime, :, 1]
                            t3 = case.data[itime, :, 2]

                            if(t1.min() == t1.max() and t2.min() == t2.max() and
                               t3.min() == t3.max()):
                                continue
                            if nnodes != ndata:
                                t1i = zeros(nnodes, dtype='float32')
                                t2i = zeros(nnodes, dtype='float32')
                                t3i = zeros(nnodes, dtype='float32')
                                t1i[j] = t1
                                t2i[j] = t2
                                t3i[j] = t3
                                t1 = t1i
                                t2 = t2i
                                t3 = t3i

                            if isinstance(dt, float):
                                header = ' %s = %.4E' % (code_name, dt)
                            else:
                                header = ' %s = %i' % (code_name, dt)

                            if has_cycle:
                                freq = case.eigrs[itime]
                                #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
                                cycle = sqrt(abs(freq))/(2. * pi)
                                header += '; freq=%g' % cycle

                            form0 = (header, None, [])
                            formi2 = form0[2]

                            cases[(subcaseID, icase, name + 'X', 1, 'node', '%g', header)] = t1
                            formi2.append((name + 'X', icase, []))
                            icase += 1

                            cases[(subcaseID, icase, name + 'Y', 1, 'node', '%g', header)] = t2
                            formi2.append((name + 'Y', icase, []))
                            icase += 1

                            cases[(subcaseID, icase, name + 'Z', 1, 'node', '%g', header)] = t3
                            formi2.append((name + 'Z', icase, []))
                            icase += 1

                            formi.append(form0)
                    else:
                        t1 = case.data[0, :, 0]
                        t2 = case.data[0, :, 1]
                        t3 = case.data[0, :, 2]
                        t123 = norm(case.data[0, :, :3], axis=1)


                        if(t1.min() == t1.max() and t2.min() == t2.max() and
                           t3.min() == t3.max() and t123.min() == t123.max()):
                            continue
                        #if t1.min() != t1.max():
                        cases[(subcaseID, icase, name + 'X', 1, 'node', '%g')] = t1
                        formi.append((name + 'X', icase, []))
                        icase += 1

                        #if t2.min() != t2.max():
                        cases[(subcaseID, icase, name + 'Y', 1, 'node', '%g')] = t2
                        formi.append((name + 'Y', icase, []))
                        icase += 1

                        #if t3.min() != t3.max():
                        cases[(subcaseID, icase, name + 'Z', 1, 'node', '%g')] = t3
                        formi.append((name + 'Z', icase, []))
                        icase += 1

                        #if t123.min() != t123.max():
                        cases[(subcaseID, icase, name + 'XYZ', 1, 'node', '%g')] = t123
                        formi.append((name + 'XYZ', icase, []))
                        icase += 1

            for (result, name) in temperature_like:
                if subcaseID in result:
                    case = result[subcaseID]
                    if case.nonlinear_factor is not None: # transient
                        continue
                    temperatures = case.data[0, :, 0]
                    cases[(subcaseID, name, 1, 'node', '%g')] = temperatures
                    formi.append((name, icase, []))
                    icase += 1
        else: # centroidal
            pass
        return icase

    def clear_nastran(self):
        self.eidMap = {}
        self.nidMap = {}
        self.eid_to_nid_map = {}
        self.element_ids = None
        self.node_ids = None

    def fill_stress_case(self, cases, model, subcaseID, formi, icase):
        if self.is_centroidal:
            icase = self._fill_stress_case_centroidal(cases, model, subcaseID, formi, icase)
        elif self.is_nodal:
            icase = self._fill_stress_case_nodal(cases, model, subcaseID, formi, icase)
        else:
            raise RuntimeError('this shouldnt happen...')
        return icase

    def _fill_stress_case_nodal(self, cases, model, subcaseID, formi, icase):
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
        return icase
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

        if is_stress:
            bars = model.barStress
        else:
            bars = model.barStrain

        if subcaseID in bars:
            case = bars[subcaseID]
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
            cases[(subcaseID, icase, word + 'XX', 1, 'node', '%.3f')] = oxx
            icase += 1
        if oyy.min() != oyy.max():
            cases[(subcaseID, icase, word + 'YY', 1, 'node', '%.3f')] = oyy
            icase += 1
        if ozz.min() != ozz.max():
            cases[(subcaseID, icase, word + 'ZZ', 1, 'node', '%.3f')] = ozz
            icase += 1

        if o1.min() != o1.max():
            cases[(subcaseID, icase, word + '1', 1, 'node', '%.3f')] = o1
            icase += 1
        if o2.min() != o2.max():
            cases[(subcaseID, icase, word + '2', 1, 'node', '%.3f')] = o2
            icase += 1
        if o3.min() != o3.max():
            cases[(subcaseID, icase, word + '3', 1, 'node', '%.3f')] = o3
            icase += 1
        if vmWord is not None:
            cases[(subcaseID, icase, vmWord, 1, 'node', '%.3f')] = ovm
            icase += 1
        return icase

    def _is_nonlinear(self, model, isubcase):
        table_types = model.get_table_types()
        msg = []
        for table_type in table_types:
            table = getattr(model, table_type)
            if isubcase in table:
                case = table[isubcase]
                if case.nonlinear_factor:
                    return True
                else:
                    return False
        raise RuntimeError('self._is_nonlinear(...) failed')

    def _get_stress_times(self, model, isubcase):
        table_types = self._get_stress_table_types()
        msg = []
        is_real = True
        is_data = False
        for table_type in table_types:
            if not hasattr(model, table_type):
                print('no table_type=%s' % table_type)
                continue
            table = getattr(model, table_type)
            if isubcase in table:
                is_data = True
                print('table_type=%s' % table_type)
                case = table[isubcase]
                if case.nonlinear_factor is not None:
                    return is_data, False, case.is_real(), case._times
                else:
                    return is_data, True, case.is_real(), zeros(1, dtype='int32')

        return is_data, False, is_real, None

    def _get_stress_table_types(self):
        table_types = [
            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress',  # vectorized
            'celas2_stress',
            'celas3_stress',
            'celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain',  # vectorized
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'crod_stress',  # vectorized
            'conrod_stress',
            'ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'crod_strain',  # vectorized
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress
            'barStress',
            # OES - isotropic CBAR strain
            'barStrain',
            # OES - isotropic CBEAM stress
            'beamStress',
            # OES - isotropic CBEAM strain
            'beamStrain',

            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress',
            'cquad4_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain',
            'cquad4_strain',

            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress',  # vectorized
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain',  # vectorized
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress
            'shearStress',
            # OES - CSHEAR strain
            'shearStrain',
            # OES - CEALS1 224, CELAS3 225
            'nonlinearSpringStress',
            # OES - GAPNL 86
            'nonlinearGapStress',
            # OES - CBUSH 226
            'nolinearBushStress',
        ]

        table_types += [
            # OES - CTRIAX6
            'ctriaxStress',
            'ctriaxStrain',

            'bushStress',
            'bushStrain',
            'bush1dStressStrain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinearRodStress',
            'nonlinearRodStrain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinearPlateStress',
            'nonlinearPlateStrain',
            'hyperelasticPlateStress',
            'hyperelasticPlateStrain',

            # OES - composite CTRIA3/CQUAD4 stress
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',

            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',

            # OGS1 - grid point stresses
            'gridPointStresses',        # tCode=26
            'gridPointVolumeStresses',  # tCode=27
        ]
        return table_types

    def _fill_stress_case_centroidal(self, cases, model, subcaseID, form, icase):
        assert isinstance(subcaseID, int), type(subcaseID)
        assert isinstance(icase, int), type(icase)
        is_data, is_static, is_real, times = self._get_stress_times(model, subcaseID)

        if not is_data:
            times = []

        #print('times = %s' % times)
        for itime, dt in enumerate(times):
            if is_static:
                formi = form
            else:
                header = 'dummy'
                form_time = [header, None, []]
                formi = form_time[2]
            is_form_time = False

            # stress
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_stress(cases, model, subcaseID, form, icase, itime, dt,
                                                             is_stress=True, is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            # strain
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_stress(cases, model, subcaseID, form, icase, itime, dt,
                                                             is_stress=False, is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            # ese
            icase, ncase, case, header, form0 = self._get_nastran_time_centroidal_strain_energy(cases, model, subcaseID, form, icase, itime, dt,
                                                            is_real=is_real, is_static=is_static)
            if ncase:
                assert ncase > 0, ncase
                if is_static:
                    formi.append(form0)
                else:
                    #form_time[0] = header
                    form_time[2].append(form0)
                    is_form_time = True

            #--------------------------
            if is_form_time:
                form.append(form_time)

        return icase

    def _get_nastran_header(self, case, dt, itime):
        if case is None:
            return None
        code_name = case.data_code['name']
        has_cycle = hasattr(case, 'mode_cycle')

        if isinstance(dt, float):
            header = ' %s = %.4E' % (code_name, dt)
        else:
            header = ' %s = %i' % (code_name, dt)

        if hasattr(case, 'mode_cycle'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = sqrt(abs(freq)) / (2. * pi)
            header += '; freq=%g' % cycle
        elif hasattr(case, 'eigrs'):
            freq = case.eigrs[itime]
            #msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))
            cycle = sqrt(abs(freq)) / (2. * pi)
            header += '; freq=%g' % cycle
        return header

    def _get_nastran_time_centroidal_strain_energy(self, cases, model, subcaseID, form, icase, itime, dt,
                                                   is_real=True, is_static=False):
        oxx = zeros(self.nElements, dtype='float32')
        oyy = zeros(self.nElements, dtype='float32')
        ozz = zeros(self.nElements, dtype='float32')
        fmt = '%g'
        header = ''
        ncase = 0
        form0 = ('Element Strain Energy', None, [])

        #op2.strainEnergy[1]
          #type=StrainEnergyObject ntimes=3 nelements=16
          #energy, percent, density
          #modes = [1, 2, 3]

        ese = model.strainEnergy[subcaseID]
        times = sorted(ese.energy.keys())  # TODO: not vectorized
        assert times[itime] == dt, 'actual=%s expected=%s' % (times[itime], dt)

        if is_static:
            percent = ese.percent
            energy = ese.energy
            density = ese.density
        else:
            percent = ese.percent[dt]
            energy = ese.energy[dt]
            density = ese.density[dt]
        for eid, p in sorted(iteritems(percent)):
            if eid not in self.eidMap:
                continue
            i = self.eidMap[eid]
            oxx[i] = energy[eid]
            oyy[i] = p
            ozz[i] = density[eid]

        case = ese
        fmt = '%.4f'
        header = self._get_nastran_header(case, dt, itime)
        cases[(subcaseID, icase, 'StrainEnergy', 1, 'centroid', fmt, header)] = oxx
        form0[2].append(('StrainEnergy', icase, []))
        icase += 1
        ncase += 1

        cases[(subcaseID, icase, 'PercentOfTotal', 1, 'centroid', fmt, header)] = oyy
        form0[2].append(('PercentOfTotal', icase, []))
        icase += 1
        ncase += 1

        cases[(subcaseID, icase, 'Density', 1, 'centroid', fmt, header)] = ozz
        form0[2].append(('Density', icase, []))
        icase += 1
        ncase += 1

        return icase, ncase, case, header, form0

    def _get_nastran_time_centroidal_stress(self, cases, model, subcaseID, form, icase, itime, dt,
                                            is_stress=True, is_real=True, is_static=False):
        ncase = 0
        case = None
        assert isinstance(subcaseID, int), type(subcaseID)
        assert isinstance(icase, int), icase
        assert isinstance(itime, int), type(itime)
        assert is_real in [True, False], is_real
        assert is_stress in [True, False], is_stress
        assert is_static in [True, False], is_static
        eids = self.element_ids
        nElements = self.nElements

        isElementOn = zeros(nElements, dtype='int8')  # is the element supported
        oxx = zeros(nElements, dtype='float32')
        oyy = zeros(nElements, dtype='float32')
        ozz = zeros(nElements, dtype='float32')

        txy = zeros(nElements, dtype='float32')
        tyz = zeros(nElements, dtype='float32')
        txz = zeros(nElements, dtype='float32')

        o1 = zeros(nElements, dtype='float32')  # max
        o2 = zeros(nElements, dtype='float32')  # mid
        o3 = zeros(nElements, dtype='float32')  # min
        ovm = zeros(nElements, dtype='float32')

        vmWord = None
        if is_stress:
            rods = [model.crod_stress, model.conrod_stress, model.ctube_stress,]
        else:
            rods = [model.crod_strain, model.conrod_strain, model.ctube_strain,]

        for result in rods:
            if subcaseID not in result:
                continue

            case = result[subcaseID]
            eidsi = case.element
            i = searchsorted(eids, eidsi)
            if len(i) != len(unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                print('eids = %s\n' % str(list(eids)))
                print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)

            #print('irod = %s' % str(i))
            isElementOn[i] = 1

            # data=[1, nnodes, 4] where 4=[axial, SMa, torsion, SMt]
            oxx[i] = case.data[itime, :, 0]
            txy[i] = case.data[itime, :, 2]
            ovm[i] = sqrt(oxx[i]**2 + 3*txy[i]**2)
            o1[i] = sqrt(oxx[i]**2 + txy[i]**2)
            o3[i] = -o1[i]

        if is_stress:
            bars = model.bar_stress
        else:
            bars = model.bar_strain

        if subcaseID in bars:  # vectorized....
            case = bars[subcaseID]
            print('class_name = %s' % str(case))

            #s1a = case.data[itime, :, 0]
            #s2a = case.data[itime, :, 1]
            #s3a = case.data[itime, :, 2]
            #s4a = case.data[itime, :, 3]

            axial = case.data[itime, :, 4]
            smaxa = case.data[itime, :, 5]
            smina = case.data[itime, :, 6]
            #MSt = case.data[itime, :, 7]

            #s1b = case.data[itime, :, 8]
            #s2b = case.data[itime, :, 9]
            #s3b = case.data[itime, :, 10]
            #s4b = case.data[itime, :, 11]

            smaxb = case.data[itime, :, 12]
            sminb = case.data[itime, :, 13]
            #MSc   = case.data[itime, :, 14]

            eidsi = case.element_node # [:, 0]

            i = searchsorted(eids, eidsi)
            #print('ibar = %s' % i)
            assert len(i) == len(unique(i)), 'ibar=%s is not unique' % str(i)

            isElementOn[i] = 1.
            oxx[i] = axial

            amax = max(smaxa, smaxb)
            amin = min(smina, sminb)
            avm = max(abs([smina, sminb,
                           smaxa, smaxb, axial]))

            #print('avm.shape', avm.shape)

            o1[i] = amax
            o3[i] = amin
            ovm[i] = avm

        if 0 and subcaseID in model.beamStress:
            case = model.beamStress[subcaseID]
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


        if is_stress:
            plates = [model.ctria3_stress, model.cquad4_stress,
                      model.ctria6_stress, model.cquad8_stress]
        else:
            plates = [model.ctria3_strain, model.cquad4_strain,
                      model.ctria6_strain, model.cquad8_strain]

        for result in plates:
            ## is tria6, quad8, bilinear quad handled?
            if subcaseID not in result:
                continue

            case = result[subcaseID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'

            nnodes_per_element = case.nnodes * 2  # *2 for every other layer
            eidsi = case.element_node[::nnodes_per_element, 0]  # ::2 is for layer skipping
            #print('eidsi =', case.element_node[:, 0])

            i = searchsorted(eids, eidsi)
            assert len(i) == len(unique(i)), 'iplate=%s is not unique' % str(i)
            #print('iplate = %s' % i)
            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            isElementOn[i] = 1.
            ndt, ntotal, nresults = case.data.shape
            if nnodes_per_element == 1:
                j = None
            else:
                j = arange(ntotal)[::nnodes_per_element]

            oxx[i] = case.data[itime, j, 1]
            oyy[i] = case.data[itime, j, 2]
            txy[i] = case.data[itime, j, 3]
            #txy[i] = case.data[itime, ::nnodes_per_element, 3]
            #ozz[i] = 0
            o1[i]  = case.data[itime, j, 5]
            #o2[i]  = case.data[itime, :, ???]
            o3[i]  = case.data[itime, j, 6]
            ovm[i] = case.data[itime, j, 7]

        if subcaseID in model.compositePlateStress:
            case = model.compositePlateStress[subcaseID]
            asdf
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


        if is_stress:
            solids = [(5, model.ctetra_stress),
                      (7, model.cpenta_stress),
                      (11, model.chexa_stress),]
        else:
            solids = [(5, model.ctetra_strain),
                      (7, model.cpenta_strain),
                      (11, model.chexa_strain),]

        for nnodes_per_element, result in solids:
            if subcaseID not in result:
                continue

            case = result[subcaseID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'

            eidsi = case.element_cid[:, 0]
            i = searchsorted(eids, eidsi)
            assert len(i) == len(unique(i)), 'isolid=%s is not unique' % str(i)
            #print('isolid = %s' % str(i))
            isElementOn[i] = 1
            #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
            #                                         txy, tyz, txz,
            #                                         o1, o2, o3, ovm]
            oxx[i] = case.data[itime, ::nnodes_per_element, 0]
            oyy[i] = case.data[itime, ::nnodes_per_element, 1]
            ozz[i] = case.data[itime, ::nnodes_per_element, 2]
            txy[i] = case.data[itime, ::nnodes_per_element, 3]
            tyz[i] = case.data[itime, ::nnodes_per_element, 4]
            txz[i] = case.data[itime, ::nnodes_per_element, 5]
            o1[i]  = case.data[itime, ::nnodes_per_element, 6]
            o2[i]  = case.data[itime, ::nnodes_per_element, 7]
            o3[i]  = case.data[itime, ::nnodes_per_element, 8]
            ovm[i] = case.data[itime, ::nnodes_per_element, 9]

        if is_stress:
            word = 'Stress'
            fmt = '%.3f'
        else:
            word = 'Strain'
            fmt = '%.4e'

        # a form is the table of output...
        # Subcase 1         <--- formi  - form_isubcase
        #    Time 1
        #        Stress     <--- form0  - the root level
        #            oxx    <--- formis - form_itime_stress
        #            oyy
        #            ozz

        header = ''
        if not is_static:
            #print('is_static = %s' % is_static)
            if case is None:
                return icase, ncase, case, formis
            header = self._get_nastran_header(case, dt, itime)
            #form_time[0] = header

        form0 = (word, None, [])
        formis = form0[2]
        # subcaseID,resultType,vectorSize,location,dataFormat
        if is_stress and itime == 0:
            if isElementOn.min() != isElementOn.max():
                if isElementOn.min() == 0:
                    ioff = where(isElementOn == 0)[0]
                    print('eids_off = %s' % self.element_ids[ioff])
                    self.log_error('eids_off = %s' % self.element_ids[ioff])
                cases[(1, icase, 'isElementOn', 1, 'centroid', '%i')] = isElementOn
                form.append(('IsElementOn', icase, []))
                icase += 1
                ncase += 1

        if oxx.min() != oxx.max():
            cases[(subcaseID, icase, word + 'XX', 1, 'centroid', fmt, header)] = oxx
            formis.append((word + 'XX', icase, []))
            icase += 1
            ncase += 1
        if oyy.min() != oyy.max():
            cases[(subcaseID, icase, word + 'YY', 1, 'centroid', fmt, header)] = oyy
            formis.append((word + 'YY', icase, []))
            icase += 1
            ncase += 1
        if ozz.min() != ozz.max():
            cases[(subcaseID, icase, word + 'ZZ', 1, 'centroid', fmt, header)] = ozz
            formis.append((word + 'ZZ', icase, []))
            ncase += 1
            icase += 1

        if txy.min() != txy.max():
            cases[(subcaseID, icase, word + 'XY', 1, 'centroid', fmt, header)] = txy
            formis.append((word + 'XY', icase, []))
            icase += 1
            ncase += 1
        if tyz.min() != tyz.max():
            cases[(subcaseID, icase, word + 'YZ', 1, 'centroid', fmt, header)] = tyz
            formis.append((word + 'YZ', icase, []))
            icase += 1
            ncase += 1
        if txz.min() != txz.max():
            cases[(subcaseID, icase, word + 'XZ', 1, 'centroid', fmt, header)] = txz
            formis.append((word + 'XZ', icase, []))
            icase += 1
            ncase += 1

        if o1.min() != o1.max():
            cases[(subcaseID, icase, word + '1', 1, 'centroid', fmt, header)] = o1
            formis.append((word + '1', icase, []))
            icase += 1
            ncase += 1
        if o2.min() != o2.max():
            cases[(subcaseID, icase, word + '2', 1, 'centroid', fmt, header)] = o2
            formis.append((word + '2', icase, []))
            icase += 1
            ncase += 1
        if o3.min() != o3.max():
            cases[(subcaseID, icase, word + '3', 1, 'centroid', fmt, header)] = o3
            formis.append((word + '3', icase, []))
            icase += 1
            ncase += 1
        if vmWord is not None:
            if not is_stress:
                max_min = max(ovm.max(), abs(ovm.min()))
                if max_min > 100:
                    raise RuntimeError('vm strain = %s' % ovm)
            cases[(subcaseID, icase, vmWord, 1, 'centroid', fmt, header)] = ovm
            formis.append((vmWord, icase, []))
            icase += 1
            ncase += 1
        #print('form0 = %s' % str(form0))
        return icase, ncase, case, header, form0

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

if __name__ == '__main__':  # pragma: no cover
    main()
