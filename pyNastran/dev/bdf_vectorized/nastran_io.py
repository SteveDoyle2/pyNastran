# pylint: disable=C0103,C0111,E1101
import os

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

from numpy import zeros, abs, mean, nan_to_num, amax, amin, array
from numpy import nan as NaN
from numpy.linalg import norm  # type: ignore

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron)

from pyNastran.dev.bdf_vectorized.bdf import BDF
    #CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    #CQUAD4, CQUAD8, CQUADR, CSHEAR,
    #CTRIA3, CTRIA6, CTRIAR, CTRIAX6,
    #CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    #CHEXA8, CHEXA20,
    #CONM2,
    #ShellElement, LineElement, SpringElement,
    #LOAD)
from pyNastran.op2.test.test_op2 import OP2
from pyNastran.f06.f06 import F06
from pyNastran.converters.nastran.nastranIO import NastranIO as NastranIO_xref
from pyNastran.utils.numpy_utils import integer_types

class NastranIO(NastranIO_xref):
    def __init__(self):
        self.is_sub_panels = False
        self.save_data = False

    def load_nastran_geometry(self, bdf_filename: str, dirname: str,
                              name: str='main') -> None:
        self.eid_map = {}
        self.nid_map = {}
        if bdf_filename is None or bdf_filename == '':
            return
        else:
            self.turn_text_off()
            self.grid.Reset()

            self.result_cases = {}
            self.ncases = 0
        for i in ('case_keys', 'icase', 'isubcase_name_map'):
            if hasattr(self, i):
                del i

        self.scalar_bar_actor.VisibilityOff()
        self.scalar_bar_actor.Modified()

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
            self.model_type = model.model_type
            model.read_bdf(bdf_filename, include_dir=dirname, punch=punch, xref=True)

        nNodes = model.grid.n
        nElements = model.elements.ne

        #nNodes = model.nNodes()
        assert nNodes > 0
        #nElements = model.nElements()
        assert nElements > 0

        if 0:
            if self.is_sub_panels:
                nsub_elements_caeros = 0
                nsub_points_caeros = 0
                for key, caero in model.caeros.items():
                    if hasattr(caero, 'panel_points_elements'):
                        npoints, nelements = caero.get_npanel_points_elements()
                        nsub_elements_caeros += npoints
                        nsub_points_caeros += nelements
                    else:
                        print('%r doesnt support panel_points_elements' % caero.type)
                ncaeros = nsub_elements_caeros
                ncaeros_points = nsub_points_caeros
            else:
                ncaeros = model.nCAeros()
                ncaeros_points = ncaeros * 4
        else:
            ncaeros = 0
            ncaeros_points = 0

        # old
        self.nNodes = nNodes
        self.nElements = nElements

        # new
        #self.nNodes = model.grid.n
        #self.nElements = model.elements.n

        self.log_info("nNodes=%i nElements=%i" % (self.nNodes, self.nElements))
        msg = model.get_bdf_stats(return_type='list')
        #self.log_info(msg)
        for msgi in msg:
            model.log.debug(msgi)

        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.card_count:
            nconm2 = model.card_count['CONM2']
        else:
            nconm2 = 0
        self.grid.Allocate(self.nElements, 1000)
        self.grid2.Allocate(ncaeros + nconm2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)

        # add the nodes
        node_ids = model.grid.node_id
        xyz_cid0 = model.grid.get_position_by_node_index()
        xmax, ymax, zmax = amax(xyz_cid0, axis=0)
        xmin, ymin, zmin = amin(xyz_cid0, axis=0)
        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s xmax=%s dz=%s" % (zmin, zmax, zmax-zmin))

        if self.save_data:
            self.xyz_cid0 = xyz_cid0
        for i, (xyz, node_id) in enumerate(zip(xyz_cid0, node_ids)):
            points.InsertPoint(i, *xyz)
            self.nid_map[node_id] = i

        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)
        self.gui.create_global_axes(dim_max)


        # add the CAERO/CONM2 elements
        j = 0
        points2 = vtk.vtkPoints()

        nsprings = 0
        #if 0:
            #elements_no_mass = [
                #'CBUSH', 'CBUSH1D', 'CFAST', 'CROD', 'CONROD',
                #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC', ]
            #for (eid, element) in sorted(model.elements.items()):
                #if (isinstance(element, (LineElement, SpringElement)) or
                    #element.type in elements_no_mass):
                        #node_ids = element.node_ids
                        #if None in node_ids:
                            #nsprings += 1

        points2.SetNumberOfPoints(ncaeros_points * 4 + nconm2 + nsprings)
        for (eid, element) in sorted(model.caeros.items()):
            if isinstance(element, (CAERO1, CAERO3, CAERO4, CAERO5)):
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
        #if 0:
            #for (eid, element) in sorted(model.elements.mass.items()):
                #if isinstance(element, CONM2):
                    ##del self.eid_map[eid]

                    ##print("element", element)
                    ##print("element.nid", element.nid)
                    ##print('node_ids', model.nodes.keys())
                    #xyz = element.nid.get_position()
                    #c = element.Centroid()
                    #d = norm(xyz-c)
                    #elem = vtk.vtkVertex()
                    ##elem = vtk.vtkSphere()

                    #if d == 0.:
                        #d = sphere_size
                    ##elem.SetRadius(d)
                    ##elem.SetCenter(points.GetPoint(self.nid_map[nid]))
                    ##print(str(element))

                    #points2.InsertPoint(j, *c)
                    #elem.GetPointIds().SetId(0, j)
                    #self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                    #j += 1
                #else:
                    #self.log_info("skipping %s" % element.type)

        self.map_elements(points, points2, self.nid_map, model, j, dim_max)

    #def _get_sphere_size(self, dim_max):
        #return 0.05 * dim_max

    def map_elements(self, points, points2, nidMap, model, j, dim_max):
        """
        Creates the elements
        """
        sphere_size = self._get_sphere_size(dim_max)
        #self.eid_map = {}

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
        elements = model.elements

        ie = 0
        eidMap = {}
        self.eid_map = eidMap

        #======================================================================
        rods = [
            elements.crod,
            elements.conrod,
        ]
        for rod in rods:
            if rod.n:
                self.model.log.debug('n%s = %s' % (rod.type, rod.n))
                i1, i2 = rod.get_node_indicies()
                eids = rod.element_id
                for eid, ii1, ii2 in zip(eids, i1, i2):
                    eidMap[eid] = ie
                    elem = vtk.vtkLine()
                    elem.GetPointIds().SetId(0, ii1)
                    elem.GetPointIds().SetId(1, ii2)
                    self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                    ie += 1

        tria3s = [
            elements.elements_shell.ctria3,
            #elements.elements_shell.ctriar,
            #elements.elements_shell.ctriax6,
        ]
        for tria3 in tria3s:
            if tria3.n:
                self.model.log.debug('n%s = %s' % (tria3.type, tria3.n))
                i1, i2, i3 = tria3.get_node_indicies()
                eids = tria3.element_id
                for eid, ii1, ii2, ii3 in zip(eids, i1, i2, i3):
                    eidMap[eid] = ie
                    elem = vtkTriangle()
                    elem.GetPointIds().SetId(0, ii1)
                    elem.GetPointIds().SetId(1, ii2)
                    elem.GetPointIds().SetId(2, ii3)
                    ie += 1

        tria6s = [
            elements.elements_shell.ctria6,
        ]
        for tria6 in tria6s:
            if tria6.n:
                i1, i2, i3 = tria3.get_node_indicies()
                eids = tria6.element_id
                for eid, ii1, ii2, ii3 in zip(eids, i1, i2, i3):
                    eidMap[eid] = ie
                    elem = vtkTriangle()
                    elem.GetPointIds().SetId(0, ii1)
                    elem.GetPointIds().SetId(1, ii2)
                    elem.GetPointIds().SetId(2, ii3)
                    ie += 1

        quad4s = [
            elements.elements_shell.cquad4,
            #elements.elements_shell.cshear,
            #elements.elements_shell.cquad,
            #elements.elements_shell.cquadr,
        ]
        for quad4 in quad4s:
            if quad4.n:
                i1, i2, i3, i4 = quad4.get_node_indicies()
                eids = quad4.element_id
                for eid, ii1, ii2, ii3, ii4 in zip(eids, i1, i2, i3, i4):
                    eidMap[eid] = ie
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, ii1)
                    elem.GetPointIds().SetId(1, ii2)
                    elem.GetPointIds().SetId(2, ii3)
                    elem.GetPointIds().SetId(3, ii4)
                    ie += 1

        quad8s = [
            elements.elements_shell.cquad8,
        ]
        for quad8 in quad8s:
            if quad8.n:
                i1, i2, i3, i4 = quad8.get_node_indicies()
                eids = quad8.element_id
                for eid, ii1, ii2, ii3, ii4 in zip(eids, i1, i2, i3, i4):
                    eidMap[eid] = ie
                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, ii1)
                    elem.GetPointIds().SetId(1, ii2)
                    elem.GetPointIds().SetId(2, ii3)
                    elem.GetPointIds().SetId(3, ii4)
                    ie += 1


        if elements.elements_solid.ctetra4.n:
            i1, i2, i3, i4 = elements.elements_solid.ctetra4.get_node_indicies()
            eids = elements.elements_solid.ctetra4.element_id
            for eid, ii1, ii2, ii3, ii4 in zip(eids, i1, i2, i3, i4):
                eidMap[eid] = ie
                elem = vtkTetra()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1
        if elements.elements_solid.ctetra10.n:
            i1, i2, i3, i4 = elements.elements_solid.ctetra10.get_node_indicies()
            eids = elements.elements_solid.ctetra10.element_id
            for eid, ii1, ii2, ii3, ii4 in zip(eids, i1, i2, i3, i4):
                eidMap[eid] = ie
                elem = vtkTetra()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1


        if elements.elements_solid.cpenta6.n:
            i1, i2, i3, i4, i5, i6 = elements.elements_solid.cpenta6.get_node_indicies()
            eids = elements.elements_solid.cpenta6.element_id
            for eid, ii1, ii2, ii3, ii4, ii5, ii6 in zip(eids, i1, i2, i3, i4, i5, i6):
                eidMap[eid] = ie
                elem = vtkWedge()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                elem.GetPointIds().SetId(4, ii5)
                elem.GetPointIds().SetId(5, ii6)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1
        if elements.elements_solid.cpenta15.n:
            i1, i2, i3, i4, i5, i6 = elements.elements_solid.cpenta15.get_node_indicies()
            eids = elements.elements_solid.cpenta15.element_id
            for eid, ii1, ii2, ii3, ii4, ii5, ii6 in zip(eids, i1, i2, i3, i4, i5, i6):
                eidMap[eid] = ie
                elem = vtkWedge()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                elem.GetPointIds().SetId(4, ii5)
                elem.GetPointIds().SetId(5, ii6)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1


        if elements.elements_solid.chexa8.n:
            i1, i2, i3, i4, i5, i6, i7, i8 = elements.elements_solid.chexa8.get_node_indicies()
            eids = elements.elements_solid.chexa8.element_id
            for eid, ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8 in zip(eids, i1, i2, i3, i4, i5, i6, i7, i8):
                eidMap[eid] = ie
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                elem.GetPointIds().SetId(4, ii5)
                elem.GetPointIds().SetId(5, ii6)
                elem.GetPointIds().SetId(6, ii7)
                elem.GetPointIds().SetId(7, ii8)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1
        if elements.elements_solid.chexa20.n:
            i1, i2, i3, i4, i5, i6, i7, i8 = elements.elements_solid.chexa20.get_node_indicies()
            eids = elements.elements_solid.chexa20.element_id
            for eid, ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8 in zip(eids, i1, i2, i3, i4, i5, i6, i7, i8):
                eidMap[eid] = ie
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, ii1)
                elem.GetPointIds().SetId(1, ii2)
                elem.GetPointIds().SetId(2, ii3)
                elem.GetPointIds().SetId(3, ii4)
                elem.GetPointIds().SetId(4, ii5)
                elem.GetPointIds().SetId(5, ii6)
                elem.GetPointIds().SetId(6, ii7)
                elem.GetPointIds().SetId(7, ii8)
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
                ie += 1

        if 0:
            for (eid, element) in sorted(model.elements.items()):
                self.eid_map[eid] = i
                #print(element.type)
                pid = 0
                if isinstance(element, (CTRIA3, CTRIAR)):
                    elem = vtkTriangle()
                    nodeIDs = element.node_ids
                    pid = element.Pid()
                    self.eid_to_nid_map[eid] = nodeIDs
                    elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                    elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                    self.grid.InsertNextCell(elem.GetCellType(),
                                             elem.GetPointIds())
                elif isinstance(element, CTRIA6):
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                    #                           nidMap[nodeIDs[a[2]]])
                    #raise RuntimeError(msg)
                    #sys.stdout.flush()

                    #elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                    #elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    #elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                    self.grid.InsertNextCell(elem.GetCellType(),
                                             elem.GetPointIds())

                elif isinstance(element, (CQUAD4, CSHEAR, CQUADR)):
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
                    pid = element.Pid()
                    self.eid_to_nid_map[eid] = nodeIDs[:4]
                    elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                    elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    self.grid.InsertNextCell(elem.GetCellType(),
                                             elem.GetPointIds())
                elif isinstance(element, CTETRA10):
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                    nodeIDs = element.node_ids
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
                elif (isinstance(element, (LineElement, SpringElement)) or
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
                    nodeIDs = element.node_ids
                    if nodeIDs[0] is None and  nodeIDs[0] is None: # CELAS2
                        del self.eid_map[eid]
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

                        #elem = vtk.vtkSphere()
                        #if d == 0.:
                            #d = sphere_size
                        #elem.SetRadius(d)
                    else:
                        # 2 points
                        #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
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
                    del self.eid_map[eid]
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
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        self.log_info("updated grid")

        cases = {}

        if 0:
            nelements = len(model.elements)
            pids = array(pids, 'int32')
            if not len(pids) == len(self.eid_map):
                msg = 'ERROR:  len(pids)=%s len(eidMap)=%s\n' % (len(pids), len(self.eid_map))
                for eid, pid in sorted(pids_dict.items()):
                    if eid not in self.eid_map:
                        msg += 'eid=%s %s' % (eid, str(model.elements[eid]))
                raise RuntimeError(msg)
            del pids_dict

        #nxs = []
        #nys = []
        #nzs = []
        #i = 0

        #for eid, element in sorted(model.elements.items()):
            #if isinstance(element, ShellElement):
                #(nx, ny, nz) = element.Normal()
            #else:
                #nx = ny = nz = 0.0
            #nxs.append(nx)
            #nys.append(ny)
            #nzs.append(nz)

        self.isubcase_name_map = {1: ['Nastran', '']}

        nElements = len(self.eid_map)
        #print("nElements = ", nElements)

        # set to True to enable nodeIDs as an result
        nidsSet = True
        if nidsSet:
            nids = model.grid.node_id
            #for (nid, nid2) in self.nid_map.items():
            #    nids[nid2] = nid
            cases[(0, 'Node_ID', 1, 'node', '%i')] = nids
            nidsSet = True

        # set to True to enable elementIDs as a result
        eidsSet = True
        Types, eids, pids = model.elements.get_element_properties()
        if eidsSet:
            #eids = zeros(nElements, dtype='int32')
            #for (eid, eid2) in self.eid_map.items():
               #eids[eid2] = eid
            #eids = model.elements.element_id
            cases[(0, 'Element_ID', 1, 'centroid', '%i')] = eids
            eidsSet = True

        # subcase_id, resultType, vectorSize, location, dataFormat
        if len(model.properties):
            #pids = model.elements.property_id
            cases[(0, 'Property_ID', 1, 'centroid', '%i')] = pids

        #self._plot_pressures(model, cases)
        #self._plot_applied_loads(model, cases)

        #if 0:
            #if not a flat plate???
            #if min(nxs) == max(nxs) and min(nxs) != 0.0:
                #subcase_id, resultType, vectorSize, location, dataFormat
                #cases[(0, 'Normal_x', 1, 'centroid', '%.1f')] = nxs
                #cases[(0, 'Normal_y', 1, 'centroid', '%.1f')] = nys
                #cases[(0, 'Normal_z', 1, 'centroid', '%.1f')] = nzs
        self.log.info(cases.keys())

        self.finish_nastran_io(cases)

    def finish_nastran_io(self, cases):  # same as Cart3d version
        self.result_cases = cases
        self.case_keys = sorted(cases.keys())
        #print("ncases = %s" % len(cases))
        #print("case_keys = %s" % self.case_keys)

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0

        self.cycle_results_explicit()  # start at ncase=0
        if self.ncases:
            self.scalar_bar_actor.VisibilityOn()
            self.scalar_bar_actor.Modified()

    def _plot_pressures(self, model, cases):
        """
        pressure act normal to a shell (as opposed to anti-normal to a solid face)
        """
        sucaseIDs = model.case_control_deck.get_subcase_list()
        for subcase_id in sucaseIDs:
            if subcase_id == 0:
                continue
            load_case_id = model.case_control_deck.get_subcase_parameter(subcase_id, 'LOAD')[0]
            load_case = model.loads[load_case_id]

            # account for scale factors
            loads2 = []
            scale_factors2 = []
            for load in load_case:
                if isinstance(load, LOAD):
                    scale_factors, loads = load.get_reduced_loads(
                        resolve_load_card=False, filter_zero_scale_factors=False)
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
                        #A, centroid, normal = elem.get_face_area_centroid_normal(
                            #load.g34_ref.nid, load.g1_ref.nid)
                        #r = centroid - p

            # if there is no applied pressure, don't make a plot
            if abs(pressures).max():
                # subcase_id, resultType, vectorSize, location, dataFormat
                cases[(0, 'Pressure Case=%i' % subcase_id, 1, 'centroid', '%.1f')] = pressures

    def _plot_applied_loads(self, model, cases):
        sucaseIDs = model.case_control_deck.get_subcase_list()
        for subcase_id in sucaseIDs:
            if subcase_id == 0:
                continue
            load_case_id, options = model.case_control_deck.get_subcase_parameter(subcase_id, 'LOAD')
            if not isinstance(load_case_id, integer_types):
                msg = 'subcase_id LOAD=%r type=%s' % (
                    subcase_id, load_case_id, type(load_case_id))
                raise TypeError(msg)
            loadCase = model.loads[load_case_id]

            # account for scale factors
            loads2 = []
            scale_factors2 = []
            for load in loadCase:
                if isinstance(load, LOAD):
                    scale_factors, loads = load.get_reduced_loads()
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
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 3.
                        for nid in node_ids[3:]:
                            pressures[eids.index(nid)] += k
                    if elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        eid = elem.eid
                        node_ids = elem.node_ids
                        k = load.pressures[0] * scale / 4.
                        for nid in node_ids[4:]:
                            pressures[eids.index(nid)] += k
            if loads[:, 0].min() != loads[:, 0].max():
                cases[(subcase_id, 'LoadX Case=%i' % subcase_id, 1, 'node', '%.1f')] = loads[:, 0]
            if loads[:, 1].min() != loads[:, 1].max():
                cases[(subcase_id, 'LoadY Case=%i' % subcase_id, 1, 'node', '%.1f')] = loads[:, 1]
            if loads[:, 2].min() != loads[:, 2].max():
                cases[(subcase_id, 'LoadZ Case=%i' % subcase_id, 1, 'node', '%.1f')] = loads[:, 2]

    def load_nastran_results(self, op2_filename, dirname):
        """
        Loads the Nastran results into the GUI
        """
        self.scalar_bar_actor.VisibilityOn()
        self.scalar_bar_actor.Modified()

        print('trying to read...%s' % op2_filename)
        if '.op2' in op2_filename:  # TODO: do this based on lower & file extension
            model = OP2(log=self.log, debug=True)
            model._results.saved = set()
            all_results = model.get_all_results()
            desired_results = [
                # nodal
                'displacements', 'velocities', 'accelerations', 'temperatures',
                'constraint_forces', 'spc_forces', 'mpc_forces',

                #'gridPointForces',
                'stress', 'solidStress', 'plateStress', 'compositePlateStress',
                'barStress', 'rodStress',
                #'strain','solidStrain', 'plateStrain', 'compositePlateStrain',
                #'barStrain', 'rodStrain',

                # untested
                'Load Vectors',
                'Applied Loads',
                'Force Vectors',

                # centroidal
                    'stress', 'solidStress', 'plateStress', 'compositePlateStress',
                    'barStress', 'rodStress',
                    #'strain','solidStrain', 'plateStrain', 'compositePlateStrain',
                    #'barStrain', 'rodStrain',
            ]
            for result in desired_results:
                if result in all_results:
                    model._results.saved.add(result)
            model.read_op2(op2_filename)

        elif '.f06' in op2_filename:  # TODO: do this based on lower & file extension
            model = F06(log=self.log, debug=True)
            model.read_f06(op2_filename)
        else:
            print("error...")
            raise NotImplementedError(op2_filename)

        #print(model.print_results())

        #case = model.displacements[1]
        #print("case = %s" % case)
        #for nodeID,translation in sorted(case.translations.items()):
            #print("nodeID=%s t=%s" % (nodeID, translation))
        #self.isubcase_name_map[self.isubcase] = [Subtitle, Label]

        cases = {}
        subcase_ids = model.isubcase_name_map.keys()
        self.isubcase_name_map = model.isubcase_name_map

        for subcase_id in subcase_ids:
            cases = self.fill_oug_oqg_case(cases, model, subcase_id)
            cases = self.fill_stress_case(cases, model, subcase_id)

        self.finish_nastran_io(cases)

    def fill_oug_oqg_case(self, cases, model, subcase_id):
        displacement_like = [
            [model.displacements, 'Displacement'],
            [model.velocities, 'Velocity'],
            [model.accelerations, 'Acceleration'],
            [model.spc_forces, 'SPC Forces'],
            [model.mpc_forces, 'MPC Forces'],

            # untested
            [model.load_vectors, 'Load Vectors'],
            #[model.applied_loads, 'Applied Loads'],
            [model.force_vectors, 'Force Vectors'],
            #[model.gridPointForces, 'GridPointForces'],  # TODO: this is buggy...
        ]
        temperature_like = [
            [model.temperatures, 'Temperature']
        ]
        nnodes = self.nNodes

        # size = 3
        for (result, name) in displacement_like:
            if subcase_id in result:
                case = result[subcase_id]

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
                    self.log.error('ERROR!!!!', case.__dict__.keys())
                    raise RuntimeError(case.__dict__.keys())

                res = getattr(case, word)

                self.log.debug('case.type =', case.__class__.__name__)
                for (nid, txyz) in res.items():
                    nid2 = self.nid_map[nid]
                    displacements[nid2] = txyz
                    xyz_displacements[nid2] = norm(txyz)
                    x_displacements[nid2] = txyz[0]
                    y_displacements[nid2] = txyz[1]
                    z_displacements[nid2] = txyz[2]

                #cases[(subcase_id, name + 'Vector', 3, 'node', '%g')] = displacements
                cases[(subcase_id, name + 'X', 1, 'node', '%g')] = x_displacements
                cases[(subcase_id, name + 'Y', 1, 'node', '%g')] = y_displacements
                cases[(subcase_id, name + 'Z', 1, 'node', '%g')] = z_displacements
                cases[(subcase_id, name + 'XYZ', 1, 'node', '%g')] = xyz_displacements
                del res, word

        # size = 1
        for (result, name) in temperature_like:
            if subcase_id in result:
                case = result[subcase_id]
                if case.nonlinear_factor is not None: # transient
                    return
                temperatures = zeros(nnodes, dtype='float32')
                for (nid, txyz) in case.translations.items():
                    nid2 = self.nid_map[nid]
                    displacements[nid2] = txyz
                    temperatures[nid2] = norm(txyz)

                    #cases[(subcase_id, name + 'Vector', 3, 'node', '%g')] = displacements
                    cases[(subcase_id, name, 1, 'node', '%g')] = temperatures
        return cases

    def clear_nastran(self):
        self.eid_map = {}
        self.nid_map = {}
        self.eid_to_nid_map = {}

    def fill_stress_case(self, cases, model, subcase_id):
        #return cases
        self._fill_stress_centroidal(cases, model, subcase_id)
        self._fill_stress_case_nodal(cases, model, subcase_id)
        return cases

    def _fill_stress_case_nodal(self, cases, model, subcase_id):
        oxx_dict = {}
        oyy_dict = {}
        ozz_dict = {}
        o1_dict = {}
        o2_dict = {}
        o3_dict = {}
        ovm_dict = {}

        for nid in self.nid_map:
            oxx_dict[nid] = []
            oyy_dict[nid] = []
            ozz_dict[nid] = []
            o1_dict[nid] = []
            o2_dict[nid] = []
            o3_dict[nid] = []
            ovm_dict[nid] = []

        vmWord = None
        if subcase_id in model.rodStress:
            case = model.rodStress[subcase_id]
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

        if subcase_id in model.barStress:
            case = model.barStress[subcase_id]
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

        if subcase_id in model.beamStress:
            case = model.beamStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.smax:
                node_ids = self.eid_to_nid_map[eid]
                oxxi = max(max(case.sxc[eid]),
                           max(case.sxd[eid]),
                           max(case.sxe[eid]),
                           max(case.sxf[eid]))
                o1i = max(case.smax[eid])
                o3i = min(case.smin[eid])
                ovmi = max(max(case.smax[eid]),
                           abs(min(case.smin[eid])))
                for nid in node_ids:
                    oxx_dict[nid].append(oxxi)
                    o1_dict[nid].append(o1i)
                    o3_dict[nid].append(o3i)
                    ovm_dict[nid].append(ovmi)

        if subcase_id in model.plateStress:
            case = model.plateStress[subcase_id]
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

        if subcase_id in model.compositePlateStress:
            case = model.compositePlateStress[subcase_id]
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

        if subcase_id in model.solidStress:
            case = model.solidStress[subcase_id]
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
        for i, nid in enumerate(sorted(self.nid_map)):
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
            cases[(subcase_id, 'StressXX', 1, 'node', '%.3f')] = oxx
        if oyy.min() != oyy.max():
            cases[(subcase_id, 'StressYY', 1, 'node', '%.3f')] = oyy
        if ozz.min() != ozz.max():
            cases[(subcase_id, 'StressZZ', 1, 'node', '%.3f')] = ozz

        if o1.min() != o1.max():
            cases[(subcase_id, 'Stress1', 1, 'node', '%.3f')] = o1
        if o2.min() != o2.max():
            cases[(subcase_id, 'Stress2', 1, 'node', '%.3f')] = o2
        if o3.min() != o3.max():
            cases[(subcase_id, 'Stress3', 1, 'node', '%.3f')] = o3
        if vmWord is not None:
            cases[(subcase_id, vmWord, 1, 'node', '%.3f')] = ovm

    def _fill_stress_centroidal(self, cases, model, subcase_id):
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
        if subcase_id in model.rodStress:
            case = model.rodStress[subcase_id]
            for eid in case.axial:
                eid2 = self.eid_map[eid]
                isElementOn[eid2] = 1.

                axial = case.axial[eid]
                torsion = case.torsion[eid]

                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2] = max(axial, torsion)  # not really
                o3[eid2] = min(axial, torsion)
                ovm[eid2] = max(abs(axial), abs(torsion))

        if subcase_id in model.barStress:
            case = model.barStress[subcase_id]
            if case.nonlinear_factor is not None:
                return
            for eid in case.axial:
                eid2 = self.eid_map[eid]
                isElementOn[eid2] = 1.

                oxx[eid2] = case.axial[eid]

                o1[eid2] = max(case.smax[eid])
                o3[eid2] = min(case.smin[eid])
                ovm[eid2] = max(
                    max(case.smax[eid]),
                    abs(min(case.smin[eid]))
                )

        if subcase_id in model.beamStress:
            case = model.beamStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            for eid in case.smax:
                eid2 = self.eid_map[eid]
                isElementOn[eid2] = 1.

                oxx[eid2] = max(
                    max(case.sxc[eid]),
                    max(case.sxd[eid]),
                    max(case.sxe[eid]),
                    max(case.sxf[eid])
                )

                o1[eid2] = max(case.smax[eid])
                o3[eid2] = min(case.smin[eid])
                ovm[eid2] = max(
                    max(case.smax[eid]),
                    abs(min(case.smin[eid]))
                )

        if subcase_id in model.plateStress:
            case = model.plateStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eid_map[eid]
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

        if subcase_id in model.compositePlateStress:
            case = model.compositePlateStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eid_map[eid]
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

        if subcase_id in model.solidStress:
            case = model.solidStress[subcase_id]
            if case.nonlinear_factor is not None: # transient
            #if case.isTransient():
                return
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eid_map[eid]
                isElementOn[eid2] = 1.

                #print('case.oxx[%i].keys() = %s' % (eid, case.oxx[eid].keys()))
                oxx[eid2] = case.oxx[eid]['CENTER']
                oyy[eid2] = case.oyy[eid]['CENTER']
                ozz[eid2] = case.ozz[eid]['CENTER']

                o1[eid2] = case.o1[eid]['CENTER']
                o2[eid2] = case.o2[eid]['CENTER']
                o3[eid2] = case.o3[eid]['CENTER']
                ovm[eid2] = case.ovmShear[eid]['CENTER']

        # subcase_id,resultType,vectorSize,location,dataFormat
        if isElementOn.min() != isElementOn.max():
            cases[(1, 'isElementOn', 1, 'centroid', '%i')] = isElementOn
        if oxx.min() != oxx.max():
            cases[(subcase_id, 'StressXX', 1, 'centroid', '%.3f')] = oxx
        if oyy.min() != oyy.max():
            cases[(subcase_id, 'StressYY', 1, 'centroid', '%.3f')] = oyy
        if ozz.min() != ozz.max():
            cases[(subcase_id, 'StressZZ', 1, 'centroid', '%.3f')] = ozz

        if o1.min() != o1.max():
            cases[(subcase_id, 'Stress1', 1, 'centroid', '%.3f')] = o1
        if o2.min() != o2.max():
            cases[(subcase_id, 'Stress2', 1, 'centroid', '%.3f')] = o2
        if o3.min() != o3.max():
            cases[(subcase_id, 'Stress3', 1, 'centroid', '%.3f')] = o3
        if vmWord is not None:
            cases[(subcase_id, vmWord, 1, 'centroid', '%.3f')] = ovm
        return cases

def main():
    from pyNastran.gui.testing_methods import add_dummy_gui_functions
    import pyNastran
    from numpy import array_equal
    pkg_path = pyNastran.__path__[0]

    bdf_op2_filenames = [
        (os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'),
         os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.op2'),),

        (os.path.join(pkg_path, '..', 'models', 'sol_101_elements', 'static_solid_shell_bar.bdf'),
         os.path.join(pkg_path, '..', 'models', 'sol_101_elements', 'static_solid_shell_bar.op2')),
    ]

    #bdf_filename = 'bottle_shell_w_holes_pmc.bdf'
    #op2_filename = 'bottle_shell_w_holes_pmc.op2'

    for bdf_filename, op2_filename in bdf_op2_filenames:
        test = NastranIO()
        test.save_data = True
        add_dummy_gui_functions(test)
        test.load_nastran_geometry(bdf_filename, '')
        test.load_nastran_results(op2_filename, '')

        keys = test.result_cases.keys()
        assert (1, 'Stress1', 1, 'centroid', '%.3f') in keys, keys

        test2 = NastranIO_xref()
        add_dummy_gui_functions(test2)
        test2.load_nastran_geometry(bdf_filename, '')
        test2.load_nastran_results(op2_filename, '')

        assert test.eidMap == test2.eidMap
        assert test.nidMap == test2.nidMap
        assert array_equal(test.xyz_cid0, test2.xyz_cid0)
        print(test.nidMap)
        print(test2.nidMap)
        del test
        del test2


if __name__ == '__main__':  # pragma: no cover
    main()
