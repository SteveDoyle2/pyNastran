"""
defines:
 - OpenFoamIO

"""
import os
from collections import OrderedDict

import numpy as np
from numpy import zeros, arange, where, unique, cross
from numpy.linalg import norm  # type: ignore

import vtk
#VTK_TRIANGLE = 5
from vtk import vtkTriangle, vtkQuad, vtkHexahedron

from pyNastran.converters.openfoam.block_mesh import BlockMesh
from pyNastran.converters.openfoam.boundary_file import Boundary
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.utils import check_path
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class OpenFoamIO:
    def __init__(self, gui):
        """creates OpenFoamIO"""
        self.gui = gui

    def get_openfoam_hex_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Hex - BlockMeshDict',
            'OpenFOAM Hex (*)', self.load_openfoam_hex_geometry,
            None, None)
        return data

    def get_openfoam_shell_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Shell - BlockMeshDict',
            'OpenFOAM Shell (*)', self.load_openfoam_shell_geometry,
            None, None)
        return data

    def get_openfoam_faces_wildcard_geometry_results_functions(self):
        data = (
            'OpenFOAM Face - BlockMeshDict',
            'OpenFOAM Face (*)', self.load_openfoam_faces_geometry,
            None, None)
        return data

    #def remove_old_openfoam_geometry(self, filename):
        #self.eid_map = {}
        #self.nid_map = {}
        #if filename is None:
            #self.scalar_bar_actor.VisibilityOff()
            #skip_reading = True
        #else:
            ##self.TurnTextOff()
            #self.grid.Reset()

            #self.resultCases = {}
            #self.nCases = 0
            #try:
                #del self.caseKeys
                #del self.iCase
                #del self.gui.isubcase_name_map
            #except:
                #print("cant delete geo")
            #skip_reading = False
        #self.scalar_bar_actor.Modified()
        #return skip_reading

    def load_openfoam_hex_geometry(self, openfoam_filename, name='main', plot=True, **kwargs):
        self.load_openfoam_geometry(openfoam_filename, 'hex')

    def load_openfoam_shell_geometry(self, openfoam_filename, name='main', plot=True, **kwargs):
        self.load_openfoam_geometry(openfoam_filename, 'shell')

    def load_openfoam_faces_geometry(self, openfoam_filename, name='main', plot=True, **kwargs):
        self.load_openfoam_geometry(openfoam_filename, 'faces')

    def load_openfoam_geometry(self, openfoam_filename, mesh_3d, name='main', plot=True, **kwargs):
        model_name = name
        #key = self.caseKeys[self.iCase]
        #case = self.resultCases[key]

        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        skip_reading = self.gui._remove_old_geometry(openfoam_filename)
        if skip_reading:
            return
        log = self.gui.log
        reset_labels = True
        #self.log.info('self.modelType=%s' % self.modelType)
        log.info('mesh_3d = %s' % mesh_3d)
        if mesh_3d in ['hex', 'shell']:
            model = BlockMesh(log=log, debug=False) # log=self.log, debug=False
        elif mesh_3d == 'faces':
            model = BlockMesh(log=log, debug=False) # log=self.log, debug=False
            boundary = Boundary(log=log, debug=False)


        self.gui.modelType = 'openfoam'
        #self.modelType = model.modelType
        log.info('openfoam_filename = %s' % openfoam_filename)

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
            self.gui.nelements = len(hexas)
        elif mesh_3d == 'shell':
            self.gui.nelements = len(quads)
        elif mesh_3d == 'faces':
            dirname = os.path.dirname(openfoam_filename)
            point_filename = os.path.join(dirname, 'points')
            face_filename = os.path.join(dirname, 'faces')
            boundary_filename = os.path.join(dirname, 'boundary')
            check_path(face_filename, 'face_filename')
            check_path(point_filename, 'point_filename')
            check_path(boundary_filename, 'boundary_filename')

            hexas = None
            patches = None
            nodes, quads, names = boundary.read_openfoam(
                point_filename, face_filename, boundary_filename)
            self.gui.nelements = len(quads) + len(tris)
        else:
            raise RuntimeError(mesh_3d)

        self.gui.nnodes = len(nodes)
        log = self.gui.log
        log.debug("nnodes = %s" % self.gui.nnodes)
        log.debug("nelements = %s" % self.gui.nelements)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        self.gui.nid_map = {}

        assert nodes is not None
        nnodes = nodes.shape[0]

        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        nodes -= np.array([xmin, ymin, zmin])
        log.info('xmax=%s xmin=%s' % (xmax, xmin))
        log.info('ymax=%s ymin=%s' % (ymax, ymin))
        log.info('zmax=%s zmin=%s' % (zmax, zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        #dim_max = (mmax - mmin).max()
        assert dim_max > 0


        # breaks the model without subracting off the delta
        #self.update_axes_length(dim_max)
        self.gui.create_global_axes(dim_max)

        #print('is_face_mesh=%s is_3d_blockmesh=%s is_surface_blockmesh=%s' % (
            #is_face_mesh, is_3d_blockmesh, is_surface_blockmesh))
        with open('points.bdf', 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')

            unames = unique(names)
            for pid in unames:
                bdf_file.write('PSHELL,%i,1,0.1\n' % pid)
            bdf_file.write('MAT1,1,1.0e7,,0.3\n')

            if is_face_mesh:
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(self.gui.nnodes)

                unodes = unique(quads)
                unodes.sort()
                # should stop plotting duplicate nodes
                for inode, node in enumerate(nodes):
                    if inode in unodes:
                        bdf_file.write('GRID,%i,,%s,%s,%s\n' % (
                            inode + 1, node[0], node[1], node[2], ))
                    points.InsertPoint(inode, node)
            else:
                points = numpy_to_vtk_points(nodes)

            #elements -= 1
            normals = None
            if is_3d_blockmesh:
                nelements = hexas.shape[0]
                cell_type_hexa8 = vtkHexahedron().GetCellType()
                create_vtk_cells_of_constant_element_type(grid, hexas, cell_type_hexa8)

            elif is_surface_blockmesh:
                nelements = quads.shape[0]
                cell_type_quad4 = vtkQuad().GetCellType()
                create_vtk_cells_of_constant_element_type(grid, quads, cell_type_quad4)

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
                        grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
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
                        grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                    else:
                        raise RuntimeError('nnodes=%s' % nnodes)
            else:
                msg = 'is_surface_blockmesh=%s is_face_mesh=%s; pick one' % (
                    is_surface_blockmesh, is_face_mesh)
                raise RuntimeError(msg)
            bdf_file.write('ENDDATA\n')

        self.gui.nelements = nelements
        grid.SetPoints(points)
        grid.Modified()

        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {0: ['OpenFoam BlockMeshDict', '']}
        cases = OrderedDict()
        ID = 1

        #print("nElements = ",nElements)
        if mesh_3d == 'hex':
            form, cases, node_ids, element_ids = self._fill_openfoam_case(
                cases, ID, nodes, nelements,
                patches, names, normals, is_surface_blockmesh)

        elif mesh_3d == 'shell':
            form, cases, node_ids, element_ids = self._fill_openfoam_case(
                cases, ID, nodes, nelements,
                patches, names, normals, is_surface_blockmesh)

        elif mesh_3d == 'faces':
            if len(names) == nelements:
                is_surface_blockmesh = True
            form, cases, node_ids, element_ids = self._fill_openfoam_case(
                cases, ID, nodes, nelements, patches,
                names, normals, is_surface_blockmesh)
        else:
            raise RuntimeError(mesh_3d)

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        if plot:
            self.gui._finish_results_io2(model_name, form, cases, reset_labels=reset_labels)
        else:
            self.gui._set_results(form, cases)


    def clear_openfoam(self):
        pass

    #def _load_openfoam_results(self, openfoam_filename):
        #raise NotImplementedError()

    def _fill_openfoam_case(self, cases, unused_ID, nodes, nelements, patches, names, normals,
                            is_surface_blockmesh):
        #nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        #new = False
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
        icase += 2

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
        return form, cases, nids, eids
