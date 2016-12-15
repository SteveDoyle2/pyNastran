#VTK_TRIANGLE = 5
from six import iteritems
from six.moves import range
import os
from numpy import zeros, arange, mean, amax, amin, array, where, unique, cross
from numpy.linalg import norm

import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron)

from pyNastran.converters.openfoam.blockMesh import BlockMesh, Boundary
from pyNastran.utils import print_bad_path


class OpenFoamIO(object):
    def __init__(self):
        pass

    def get_openfoam_hex_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Hex - BlockMeshDict',
            'OpenFOAM Hex (*)', self.load_openfoam_geometry_hex,
             None, None)
        return data

    def get_openfoam_shell_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Shell - BlockMeshDict',
            'OpenFOAM Shell (*)', self.load_openfoam_geometry_shell,
             None, None)
        return data

    def get_openfoam_faces_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Face - BlockMeshDict',
            'OpenFOAM Face (*)', self.load_openfoam_geometry_faces,
             None, None)
        return data

    def remove_old_openfoam_geometry(self, filename):
        self.eidMap = {}
        self.nidMap = {}
        if filename is None:
            self.scalarBar.VisibilityOff()
            skipReading = True
        else:
            self.TurnTextOff()
            self.grid.Reset()
            self.grid2.Reset()

            self.resultCases = {}
            self.nCases = 0
            try:
                del self.caseKeys
                del self.iCase
                del self.iSubcaseNameMap
            except:
                print("cant delete geo")
                pass
            skipReading = False
        self.scalarBar.Modified()
        return skipReading

    def load_openfoam_geometry_hex(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'hex')

    def load_openfoam_geometry_shell(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'shell')

    def load_openfoam_geometry_faces(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'faces')

    def load_openfoam_geometry(self, openfoam_filename, dirname, mesh_3d, plot=True):
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skipReading = self.remove_old_openfoam_geometry(openfoam_filename)
        if skipReading:
            return
        #print('self.modelType=%s' % self.modelType)
        print('mesh_3d = %s' % mesh_3d)
        if mesh_3d in ['hex', 'shell']:
            model = BlockMesh(log=self.log, debug=False) # log=self.log, debug=False
        elif mesh_3d == 'faces':
            model = BlockMesh(log=self.log, debug=False) # log=self.log, debug=False
            boundary = Boundary(log=self.log, debug=False)


        self.modelType = 'openfoam'
        #self.modelType = model.modelType
        print('openfoam_filename = %s' % openfoam_filename)

        is_face_mesh = False
        if mesh_3d == 'hex':
            is_3d_blockmesh = True
            is_surface_blockmesh = False
            (nodes, hexas, quads, names, patches) = model.read_openfoam(openfoam_filename)
        elif mesh_3d == 'shell':
            is_3d_blockmesh = False
            is_surface_blockmesh = True
            (nodes, hexas, quads, names, patches) = model.read_openfoam(openfoam_filename)
        elif mesh_3d == 'faces':
            is_3d_blockmesh = False
            is_surface_blockmesh = False
            is_face_mesh = True
            #(nodes, hexas, quads, names, patches) = model.read_openfoam(openfoam_filename)
        else:
            raise RuntimeError(mesh_3d)

        tris = []


        if mesh_3d == 'hex':
            self.nElements = len(hexas)
        elif mesh_3d == 'shell':
            self.nElements = len(quads)
        elif mesh_3d == 'faces':
            point_filename = os.path.join(dirname, 'points')
            face_filename = os.path.join(dirname, 'faces')
            boundary_filename = os.path.join(dirname, 'boundary')
            assert os.path.exists(face_filename), print_bad_path(face_filename)
            assert os.path.exists(point_filename), print_bad_path(point_filename)
            assert os.path.exists(boundary_filename), print_bad_path(boundary_filename)

            hexas = None
            patches = None
            nodes, quads, names = boundary.read_openfoam(point_filename, face_filename, boundary_filename)
            self.nElements = len(quads) + len(tris)
        else:
            raise RuntimeError(mesh_3d)

        self.nNodes = len(nodes)
        print("nNodes = %s" % self.nNodes)
        print("nElements = %s" % self.nElements)


        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes, three = nodes.shape

        nid = 0
        #print("nnodes=%s" % nnodes)
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.update_axes_length(dim_max)

        f = open('points.bdf', 'wb')
        f.write('CEND\n')
        f.write('BEGIN BULK\n')

        unames = unique(names)
        for pid in unames:
            f.write('PSHELL,%i,1,0.1\n' % pid)
        f.write('MAT1,1,1.0e7,,0.3\n')

        if is_face_mesh:
            unodes = unique(quads)
            unodes.sort()
            # should stop plotting duplicate nodes
            for inode, node in enumerate(nodes):
                if inode in unodes:
                    f.write('GRID,%i,,%s,%s,%s\n' % (inode + 1, node[0], node[1], node[2], ))
                points.InsertPoint(inode, node)
        else:
            #print(nodes)
            for inode, node in enumerate(nodes):
                points.InsertPoint(inode, node)

        #elements -= 1
        normals = None
        if is_3d_blockmesh:
            nelements, three = hexas.shape
            for eid, element in enumerate(hexas):
                #print(element)
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                elem.GetPointIds().SetId(3, element[3])
                elem.GetPointIds().SetId(4, element[4])
                elem.GetPointIds().SetId(5, element[5])
                elem.GetPointIds().SetId(6, element[6])
                elem.GetPointIds().SetId(7, element[7])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

                #elem = vtkTriangle()
                #node_ids = elements[eid, :]
                #elem.GetPointIds().SetId(0, node_ids[0])
                #elem.GetPointIds().SetId(1, node_ids[1])
                #elem.GetPointIds().SetId(2, node_ids[2])
                #self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        elif is_surface_blockmesh:
            nelements, four = quads.shape
            for eid, element in enumerate(quads):
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                elem.GetPointIds().SetId(3, element[3])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
        elif is_face_mesh:
            elems = quads
            nelements = quads.shape[0]
            nnames = len(names)
            normals = zeros((nelements, 3), dtype='float32')
            assert nnames == nelements, 'nnames=%s nelements=%s names.max=%s names.min=%s' % (nnames, nelements, names.max(), names.min())
            for eid, element in enumerate(elems):
                #print('element = %s' % element)
                ineg = where(element == -1)[0]

                nnodes = 4
                if ineg:
                    nnodes = ineg.max()

                #pid = 1
                pid = names[eid]
                if nnodes == 3: # triangle!
                    f.write('CTRIA3,%i,%i,%i,%i,%i\n' % (eid+1, pid, element[0]+1, element[1]+1, element[2]+1))
                    elem = vtkTriangle()
                    a = nodes[element[1], :] - nodes[element[0], :]
                    b = nodes[element[2], :] - nodes[element[0], :]
                    n = cross(a, b)
                    normals[eid, :] = n / norm(n)

                    elem.GetPointIds().SetId(0, element[0])
                    elem.GetPointIds().SetId(1, element[1])
                    elem.GetPointIds().SetId(2, element[2])
                    self.grid.InsertNextCell(elem.GetCellType(),
                                             elem.GetPointIds())
                elif nnodes == 4:
                    f.write('CQUAD4,%i,%i,%i,%i,%i,%i\n' % (eid+1, pid, element[0]+1, element[1]+1, element[2]+1, element[3]+1))
                    a = nodes[element[2], :] - nodes[element[0], :]
                    b = nodes[element[3], :] - nodes[element[1], :]
                    n = cross(a, b)
                    normals[eid, :] = n / norm(n)

                    elem = vtkQuad()
                    elem.GetPointIds().SetId(0, element[0])
                    elem.GetPointIds().SetId(1, element[1])
                    elem.GetPointIds().SetId(2, element[2])
                    elem.GetPointIds().SetId(3, element[3])
                    self.grid.InsertNextCell(elem.GetCellType(),
                                             elem.GetPointIds())
                else:
                    asdf
        else:
            asdf

        self.nElements = nelements
        self.grid.SetPoints(points)
        #self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #print(dir(self.grid) #.SetNumberOfComponents(0))
        #self.grid.GetCellData().SetNumberOfTuples(1);
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        #self.grid2.Modified()
        self.grid.Update()
        #self.grid2.Update()
        print("updated grid")

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()


        #assert loads is not None
        #if 'Mach' in loads:
            #avgMach = mean(loads['Mach'])
            #note = ':  avg(Mach)=%g' % avgMach
        #else:
            #note = ''
        #self.iSubcaseNameMap = {1: ['Cart3d%s' % note, '']}
        self.iSubcaseNameMap = {1: ['OpenFoam BlockMeshDict', '']}
        cases = {}
        ID = 1

        #print("nElements = ",nElements)
        f.write('ENDDATA\n')
        f.close()

        if mesh_3d == 'hex':
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements, patches, names, normals, is_surface_blockmesh)
        elif mesh_3d == 'shell':
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements, patches, names, normals, is_surface_blockmesh)
        elif mesh_3d == 'faces':
            if len(names) == nelements:
                is_surface_blockmesh = True
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements, patches, names, normals, is_surface_blockmesh)
        else:
            raise RuntimeError(mesh_3d)

        if plot:
            self._finish_results_io2(form, cases)

    def clear_openfoam(self):
        pass

    def _load_openfoam_results(self, openfoam_filename, dirname):
        model = Cart3DReader(log=self.log, debug=False)
        #self.modelType = model.modelType
        #(nodes, elements, regions, loads) = model.read_cart3d(cart3dFileName)

        model.infilename = cart3d_filename
        if is_binary(infilename):
            model.infile = open(cart3d_filename, 'rb')
            (model.nPoints, model.nElements) = self.read_header_binary()
            points = model.read_points_binary(self.nPoints)
            elements = model.read_elements_binary(self.nElements)
            regions = model.read_regions_binary(self.nElements)
            #loads = {}
        else:
            model.infile = open(cart3d_filename, 'r')
            model.read_header_ascii()
            points = model.read_points_ascii(bypass=True)
            elements = model.read_elements_ascii(bypass=True)
            regions = model.read_regions_ascii(bypass=True)
            loads = model.read_results_ascii(0, model.infile, result_names=result_names)
        self.load_cart3d_geometry(cart3d_filename, dirname)


    def _fill_openfoam_case(self, cases, ID, nodes, nelements, patches, names, normals, is_surface_blockmesh):
        #result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                        #'rhoU', 'rhoV', 'rhoW', 'rhoE']
        #nelements, three = elements.shape
        nnodes, three = nodes.shape

        new = False
        results_form = []
        geometry_form = [
            #('Region', 0, []),
            ('ElementID', 0, []),
            ('NodeID', 1, []),
        ]

        eids = arange(nelements) + 1
        nids = arange(0, nnodes)
        if new:
            cases_new = [None]
            cases_new[0] = (ID, eids, 'ElementID', 'centroid', '%.0f')
            cases_new[1] = (ID, nids, 'NodeID', 'node', '%.0f')
            #cases_new[2] = (ID, regions, 'Region', 'centroid', '%.0f')
        else:
            cases[(ID, 0, 'ElementID', 1, 'centroid', '%.0f')] = eids
            cases[(ID, 1, 'NodeID', 1, 'node', '%.0f')] = nids
            icase = 2
            if is_surface_blockmesh:
                if patches is not None:
                    cases[(ID, icase, 'Patch', 1, 'centroid', '%.0f')] = patches
                    gf = ('PatchType', icase, [])
                    geometry_form.append(gf)
                    icase += 1

                if names is not None:
                    cases[(ID, icase, 'Name', 1, 'centroid', '%.0f')] = names
                    gf = ('Names', icase, [])
                    geometry_form.append(gf)
                    icase += 1
                else:
                    asdf

            if normals is not None:
                cases[(ID, icase, 'NormalX', 1, 'centroid', '%.0f')] = normals[:, 0]
                cases[(ID, icase + 1, 'NormalY', 1, 'centroid', '%.0f')] = normals[:, 1]
                cases[(ID, icase + 2, 'NormalZ', 1, 'centroid', '%.0f')] = normals[:, 2]
                geometry_form.append( ('NormalX', icase, []) )
                geometry_form.append( ('NormalY', icase, []) )
                geometry_form.append( ('NormalZ', icase, []) )
                icase += 3

        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
