"""
defines:
 - MarkActions
"""
from __future__ import print_function
from six import string_types
import numpy as np
import vtk

from pyNastran.utils.numpy_utils import integer_types


class MarkActions(object):
    """defines MarkActions"""
    def __init__(self, gui):
        """creates MarkActions"""
        self.gui = gui

    def create_annotation(self, text, slot, x, y, z):
        """
        Creates the actual annotation and appends it to slot

        Parameters
        ----------
        text : str
            the text to display
        label_actors[icase] : List[annotation]
            where to place the annotation
            icase : int
                the key in label_actors to slot the result into
            annotation : vtkBillboardTextActor3D
                the annotation object
        x, y, z : float
            the position of the label
        """
        if not isinstance(slot, list):
            msg = 'slot=%r type=%s' % (slot, type(slot))
            raise TypeError(msg)
        # http://nullege.com/codes/show/src%40p%40y%40pymatgen-2.9.6%40pymatgen%40vis%40structure_vtk.py/395/vtk.vtkVectorText/python

        #self.convert_units(icase, result_value, x, y, z)

        text_actor = vtk.vtkBillboardTextActor3D()
        text_actor.SetPosition(x, y, z)
        text_actor.SetInput(text)
        text_actor.PickableOff()
        text_actor.DragableOff()
        #text_actor.SetPickable(False)

        #text_actor.SetPosition(actor.GetPosition())
        text_prop = text_actor.GetTextProperty()
        text_prop.SetFontSize(self.gui.settings.annotation_size)
        text_prop.SetFontFamilyToArial()
        text_prop.BoldOn()
        text_prop.ShadowOn()
        text_prop.SetColor(self.gui.settings.annotation_color)
        text_prop.SetJustificationToCentered()

        # finish adding the actor
        self.gui.rend.AddActor(text_actor)

        #self.label_actors[icase].append(text_actor)
        slot.append(text_actor)

        #print('added label actor %r; icase=%s' % (text, icase))
        #print(self.label_actors)

        #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
        #camera.GetPosition()
        #camera.GetClippingRange()
        #camera.GetFocalPoint()

    def get_result_by_xyz_cell_id(self, node_xyz, cell_id):
        """won't handle multiple cell_ids/node_xyz"""
        case_key = self.gui.case_keys[self.gui.icase_fringe]
        result_name = self.gui.result_name

        cell = self.gui.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()

        #node_xyz = array(node_xyz, dtype='float32')
        #point0 = array(points.GetPoint(0), dtype='float32')
        #dist_min = norm(point0 - node_xyz)
        point0 = points.GetPoint(0)
        dist_min = vtk.vtkMath.Distance2BetweenPoints(point0, node_xyz)

        point_min = point0
        imin = 0
        for ipoint in range(1, nnodes):
            #point = array(points.GetPoint(ipoint), dtype='float32')
            #dist = norm(point - node_xyz)
            point = points.GetPoint(ipoint)
            dist = vtk.vtkMath.Distance2BetweenPoints(point, node_xyz)
            if dist < dist_min:
                dist_min = dist
                imin = ipoint
                point_min = point

        node_id = cell.GetPointId(imin)
        xyz = np.array(point_min, dtype='float32')
        case = self.gui.result_cases[case_key]
        assert isinstance(case_key, integer_types), case_key
        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        case = obj.get_result(i, res_name)
        result_values = case[node_id]
        assert not isinstance(xyz, int), xyz
        return result_name, result_values, node_id, xyz

    def mark_elements_by_different_case(self, eids, icase_result, icase_to_apply):
        """
        Marks a series of elements with custom text labels

        Parameters
        ----------
        eids : int, List[int]
            the elements to apply a message to
        icase_result : int
            the case to draw the result from
        icase_to_apply : int
            the key in label_actors to slot the result into

        TODO: fix the following
        correct   : applies to the icase_to_apply
        incorrect : applies to the icase_result

        Examples
        --------
        .. code-block::

          eids = [16563, 16564, 8916703, 16499, 16500, 8916699,
                  16565, 16566, 8916706, 16502, 16503, 8916701]
          icase_result = 22
          icase_to_apply = 25
          self.mark_elements_by_different_case(eids, icase_result, icase_to_apply)
        """
        if icase_result not in self.gui.label_actors:
            msg = 'icase_result=%r not in label_actors=[%s]' % (
                icase_result, ', '.join(self.gui.label_actors))
            self.gui.log_error(msg)
            return
        if icase_to_apply not in self.gui.label_actors:
            msg = 'icase_to_apply=%r not in label_actors=[%s]' % (
                icase_to_apply, ', '.join(self.gui.label_actors))
            self.gui.log_error(msg)
            return

        eids = np.unique(eids)
        unused_neids = len(eids)
        #centroids = np.zeros((neids, 3), dtype='float32')
        ieids = np.searchsorted(self.gui.element_ids, eids)
        #print('ieids = ', ieids)

        for cell_id in ieids:
            centroid = self.gui.cell_centroid(cell_id)
            unused_result_name, result_values, unused_xyz = self.get_result_by_cell_id(
                cell_id, centroid, icase_result)
            texti = '%s' % result_values
            xi, yi, zi = centroid
            self.create_annotation(texti, self.gui.label_actors[icase_to_apply], xi, yi, zi)
        self.gui.log_command('mark_elements_by_different_case(%s, %s, %s)' % (
            eids, icase_result, icase_to_apply))
        self.gui.vtk_interactor.Render()

    def get_result_by_cell_id(self, cell_id, world_position, icase=None):
        """should handle multiple cell_ids"""
        if icase is None:
            icase = self.gui.icase_fringe
        case_key = self.gui.case_keys[icase] # int for object
        case = self.gui.result_cases[case_key]

        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        case = obj.get_result(i, res_name)

        try:
            result_values = case[cell_id]
        except IndexError:
            msg = ('case[cell_id] is out of bounds; length=%s\n'
                   'result_name=%r cell_id=%r case_key=%r\n' % (
                       len(case), res_name, cell_id, case_key))
            raise IndexError(msg)

        cell = self.gui.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        cell_type = cell.GetCellType()
        if cell_type in [5, 9, 22, 23, 28]:  # CTRIA3, CQUAD4, CTRIA6, CQUAD8, CQUAD
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [10, 12, 13, 14]: # CTETRA4, CHEXA8, CPENTA6, CPYRAM5
            # TODO: No idea how to get the center of the face
            #       vs. a point on a face that's not exposed
            #faces = cell.GetFaces()
            #nfaces = cell.GetNumberOfFaces()
            #for iface in range(nfaces):
                #face = cell.GetFace(iface)
                #points = face.GetPoints()
            #faces
            xyz = world_position
        elif cell_type in [24, 25, 26, 27]: # CTETRA10, CHEXA20, CPENTA15, CPYRAM13
            xyz = world_position
        elif cell_type in [3]: # CBAR, CBEAM, CELASx, CDAMPx, CBUSHx
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in [21]: # CBEND
            # 21-QuadraticEdge
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        else:
            #self.log.error(msg)
            msg = 'cell_type=%s nnodes=%s; icase=%s result_values=%s' % (
                cell_type, nnodes, icase, result_values)
            self.gui.log.error(msg)
            #VTK_LINE = 3

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

            #VTK_PYRAMID = 14
            #VTK_QUADRATIC_PYRAMID = 27
            raise NotImplementedError(msg)
        return res_name, result_values, xyz

    def mark_nodes(self, nids, icase, text):
        """
        Marks a series of nodes with custom text labels

        Parameters
        ----------
        nids : int, List[int]
            the nodes to apply a message to
        icase : int
            the key in label_actors to slot the result into
        text : str, List[str]
            the text to display

        0 corresponds to the NodeID result
        self.mark_nodes(1, 0, 'max')
        self.mark_nodes(6, 0, 'min')
        self.mark_nodes([1, 6], 0, 'max')
        self.mark_nodes([1, 6], 0, ['max', 'min'])
        """
        if icase not in self.gui.label_actors:
            msg = 'icase=%r not in label_actors=[%s]' % (
                icase, ', '.join(self.gui.label_actors))
            self.gui.log_error(msg)
            return
        i = np.searchsorted(self.gui.node_ids, nids)
        if isinstance(text, string_types):
            text = [text] * len(i)
        else:
            assert len(text) == len(i)

        xyz = self.gui.xyz_cid0[i, :]
        for (xi, yi, zi), texti in zip(xyz, text):
            self.create_annotation(texti, self.gui.label_actors[icase], xi, yi, zi)
        self.gui.vtk_interactor.Render()

    #def __mark_nodes_by_result(self, nids, icases):
        #"""
        ## mark the node 1 with the NodeID (0) result
        #self.mark_nodes_by_result_case(1, 0)

        ## mark the nodes 1 and 2 with the NodeID (0) result
        #self.mark_nodes_by_result_case([1, 2], 0)

        ## mark the nodes with the NodeID (0) and ElementID (1) result
        #self.mark_nodes_by_result_case([1, 2], [0, 1])
        #"""
        #i = np.searchsorted(self.gui.node_ids, nids)
        #if isinstance(icases, int):
            #icases = [icases]

        #for icase in icases:
            #if icase not in self.gui.label_actors:
                #msg = 'icase=%r not in label_actors=[%s]' % (
                    #icase, ', '.join(self.gui.label_actors))
                #self.gui.log_error(msg)
                #continue

            #for node_id in i:
                #jnid = np.where(node_id == self.gui.node_ids)[0]
                #world_position = self.gui.xyz_cid0[jnid, :]
                #out = self.get_result_by_xyz_node_id(world_position, node_id)
                #_result_name, unused_result_value, node_id, node_xyz = out
                #xi, yi, zi = node_xyz
                #texti = 'test'
                #self.create_annotation(texti, self.gui.label_actors[icase], xi, yi, zi)
        #self.gui.vtk_interactor.Render()
