from __future__ import print_function
import os
from six import iteritems
from six.moves import range
from numpy import zeros, arange, mean, amax, amin, array, where, unique, cross
from numpy.linalg import norm

import vtk
#VTK_TRIANGLE = 5
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
                 vtkQuadraticWedge, vtkQuadraticHexahedron)

from pyNastran.converters.openfoam.block_mesh import BlockMesh, Boundary
from pyNastran.gui.gui_objects.gui_result import GuiResult
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
            skip_reading = True
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
            skip_reading = False
        self.scalarBar.Modified()
        return skip_reading

    def load_openfoam_geometry_hex(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'hex')

    def load_openfoam_geometry_shell(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'shell')

    def load_openfoam_geometry_faces(self, openfoam_filename, dirname, plot=True):
        self.load_openfoam_geometry(openfoam_filename, dirname, 'faces')

    def load_openfoam_geometry(self, openfoam_filename, dirname, mesh_3d, plot=True):
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        if skip_reading:
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
            nodes, quads, names = boundary.read_openfoam(
                point_filename, face_filename, boundary_filename)
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

        with open('points.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')

            unames = unique(names)
            for pid in unames:
                bdf_file.write('PSHELL,%i,1,0.1\n' % pid)
            bdf_file.write('MAT1,1,1.0e7,,0.3\n')

            if is_face_mesh:
                unodes = unique(quads)
                unodes.sort()
                # should stop plotting duplicate nodes
                for inode, node in enumerate(nodes):
                    if inode in unodes:
                        bdf_file.write('GRID,%i,,%s,%s,%s\n' % (
                            inode + 1, node[0], node[1], node[2], ))
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

                    #elem.GetCellType() = 5  # vtkTriangle
                    #self.grid.InsertNextCell(5, elem.GetPointIds())
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
                if nnames != nelements:
                    msg = 'nnames=%s nelements=%s names.max=%s names.min=%s' % (
                        nnames, nelements, names.max(), names.min())
                    raise RuntimeError(msg)
                for eid, element in enumerate(elems):
                    #print('element = %s' % element)
                    ineg = where(element == -1)[0]

                    nnodes = 4
                    if ineg:
                        nnodes = ineg.max()

                    #pid = 1
                    pid = names[eid]
                    if nnodes == 3: # triangle!
                        bdf_file.write('CTRIA3,%i,%i,%i,%i,%i\n' % (
                            eid+1, pid, element[0]+1, element[1]+1, element[2]+1))
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
                        bdf_file.write('CQUAD4,%i,%i,%i,%i,%i,%i\n' % (
                            eid+1, pid, element[0]+1, element[1]+1, element[2]+1, element[3]+1))
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
                        raise RuntimeError('nnodes=%s' % nnodes)
            else:
                msg = 'is_surface_blockmesh=%s is_face_mesh=%s; pick one' % (
                    is_surface_blockmesh, is_face_mesh)
                raise RuntimeError(msg)
            bdf_file.write('ENDDATA\n')

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
        self.iSubcaseNameMap = {0: ['OpenFoam BlockMeshDict', '']}
        cases = {}
        ID = 1

        #print("nElements = ",nElements)
        if mesh_3d == 'hex':
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements,
                                                   patches, names, normals, is_surface_blockmesh)
        elif mesh_3d == 'shell':
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements,
                                                   patches, names, normals, is_surface_blockmesh)
        elif mesh_3d == 'faces':
            if len(names) == nelements:
                is_surface_blockmesh = True
            form, cases = self._fill_openfoam_case(cases, ID, nodes, nelements, patches,
                                                   names, normals, is_surface_blockmesh)
        else:
            raise RuntimeError(mesh_3d)

        if plot:
            self._finish_results_io2(form, cases)

    def clear_openfoam(self):
        pass

    def _load_openfoam_results(self, openfoam_filename, dirname):
        raise NotImplementedError()

    def _fill_openfoam_case(self, cases, ID, nodes, nelements, patches, names, normals,
                            is_surface_blockmesh):
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
        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nids)

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))

        if is_surface_blockmesh:
            if patches is not None:
                patch_res = GuiResult(0, header='Patch', title='Patch',
                                      location='centroid', scalar=patches)
                cases[icase] = (patch_res, (0, 'Patch'))
                formi = ('PatchType', icase, [])
                geometry_form.append(formi)
                icase += 1

            if names is not None:
                name_res = GuiResult(0, header='Name', title='Name',
                                     location='centroid', scalar=names)
                cases[icase] = (name_res, (0, 'Name'))
                formi = ('Names', icase, [])
                geometry_form.append(formi)
                icase += 1
            else:
                raise RuntimeError('names is None...')

        if normals is not None:
            nx_res = GuiResult(0, header='NormalX', title='NormalX',
                               location='node', data_format='%.1f',
                               scalar=normals[:, 0])
            ny_res = GuiResult(0, header='NormalY', title='NormalY',
                               location='node', data_format='%.1f',
                               scalar=normals[:, 1])
            nz_res = GuiResult(0, header='NormalZ', title='NormalZ',
                               location='node', data_format='%.1f',
                               scalar=normals[:, 2])
            geometry_form.append(('NormalX', icase, []))
            geometry_form.append(('NormalY', icase, []))
            geometry_form.append(('NormalZ', icase, []))
            cases[icase] = (nx_res, (0, 'NormalX'))
            cases[icase + 1] = (ny_res, (0, 'NormalY'))
            cases[icase + 2] = (nz_res, (0, 'NormalZ'))
            icase += 3

        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
