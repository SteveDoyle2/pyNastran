"""
defines:
 - MarkActions

"""
from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np
from vtkmodules.vtkCommonDataModel import vtkSelection
from vtkmodules.vtkFiltersExtraction import vtkExtractSelection

from pyNastran.gui.vtk_rendering_core import vtkBillboardTextActor3D

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.vtk_common_core import vtkMath
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points, create_unstructured_point_grid
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow


class MarkActions:
    """defines MarkActions"""
    def __init__(self, gui):
        """creates MarkActions"""
        self.gui = gui

    @property
    def log(self):
        return self.gui.log

    def create_annotation(self, text: str, x, y, z):
        """
        Creates the actual annotation and appends it to slot

        Parameters
        ----------
        text : str
            the text to display
        x, y, z : float
            the position of the label

        Returns
        -------
        annotation : vtkBillboardTextActor3D
            the annotation object

        """
        annotation = create_annotation(self.gui, text, x, y, z)
        return annotation

    def get_result_by_xyz_cell_id(self, node_xyz: np.ndarray,
                                  cell_id: int,
                                  icase: Optional[int]=None):
        """won't handle multiple cell_ids/node_xyz"""
        if icase is None:
            icase = self.gui.icase_fringe
        case_key = self.gui.case_keys[icase]
        result_name = self.gui.result_name

        grid = self.gui.grid_selected
        cell = grid.GetCell(cell_id)
        if cell is None:
            return
        #if cell is None:
            #print('cell_id = %s' % cell)
            #self.log.error('cell_id = %s' % cell)
            #return
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()

        #node_xyz = array(node_xyz, dtype='float32')
        #point0 = array(points.GetPoint(0), dtype='float32')
        #dist_min = norm(point0 - node_xyz)

        try:
            point0 = points.GetPoint(0)
        except ValueError:
            #ValueError: expects 0 <= id && id < GetNumberOfPoints()
            return None
        dist_min = vtkMath.Distance2BetweenPoints(point0, node_xyz)

        point_min = point0
        imin = 0
        for ipoint in range(1, nnodes):
            #point = array(points.GetPoint(ipoint), dtype='float32')
            #dist = norm(point - node_xyz)
            point = points.GetPoint(ipoint)
            dist = vtkMath.Distance2BetweenPoints(point, node_xyz)
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
        #method = obj.get_methods(i, res_name)[0]

        is_value, values = get_plottable_values(obj, i, res_name)
        if not is_value:
            return None
        result_values = values[node_id]
        assert not isinstance(xyz, int), xyz
        return result_name, result_values, node_id, xyz

    def mark_elements_by_different_case(self, eids, icase_result: int, icase_to_apply: int,
                                        stop_on_failure: bool=False, show_command: bool=True):
        """
        Marks a series of elements with custom text labels

        Parameters
        ----------
        eids : int, list[int]
            the elements to apply a message to
        icase_result : int
            the case to draw the result from
        icase_to_apply : int
            the key in label_actors to slot the result into
        show_command : bool; default=True
            should the command be shown

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
            actorsi = list(self.gui.label_actors.keys())
            msg = 'icase_result=%r not in label_actors=%s' % (icase_result, actorsi)
            if stop_on_failure:
                raise RuntimeError(msg)
            self.gui.log_error(msg)
            return
        if icase_to_apply not in self.gui.label_actors:
            actorsi = list(self.gui.label_actors.keys())
            msg = 'icase_to_apply=%r not in label_actors=%s' % (icase_to_apply, actorsi)
            if stop_on_failure:
                raise RuntimeError(msg)
            self.gui.log.error(msg)
            return
        label_actors = self.gui.label_actors[icase_to_apply]

        eids = np.unique(eids)
        unused_neids = len(eids)
        ieids = np.searchsorted(self.gui.element_ids, eids)
        #print('ieids = ', ieids)

        cell_ids = ieids
        #cell_centroids
        unused_result_name, result_values = self.get_result_by_cell_ids(
            cell_ids, icase=icase_result)
        for cell_id, result_value in zip(cell_ids, result_values):
            centroid = self.gui.cell_centroid(cell_id)
            #unused_result_name, result_values, unused_xyz = self.get_result_by_cell_id(
                #cell_id, centroid, icase_result)
            texti = '%s' % result_value
            xi, yi, zi = centroid
            label_actors.append(self.create_annotation(texti, xi, yi, zi))
        if show_command:
            self.gui.log_command('self.mark_elements_by_different_case(%s, %s, %s)' % (
                eids, icase_result, icase_to_apply))
        self.gui.vtk_interactor.Render()

    def get_result_by_cell_id(self, cell_id: int,
                              world_position: np.ndarray,
                              icase: Optional[int]=None) -> tuple[str, Any, np.ndarray]:
        """should handle multiple cell_ids"""
        if icase is None:
            icase = self.gui.icase_fringe
        case_key = self.gui.case_keys[icase] # int for object
        case = self.gui.result_cases[case_key]

        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        #method = obj.get_methods(i, res_name)[0]
        is_value, values = get_plottable_values(obj, i, res_name)
        if not is_value:
            return None

        try:
            result_values = values[cell_id]
        except IndexError:
            msg = ('values[cell_id] is out of bounds; length=%s\n'
                   'result_name=%r cell_id=%r case_key=%r\n' % (
                       len(values), res_name, cell_id, case_key))
            raise IndexError(msg)

        grid = self.gui.grid_selected
        cell = grid.GetCell(cell_id)
        if cell is None:
            return
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        cell_type = cell.GetCellType()
        if cell_type in {5, 9, 22, 23, 28}:  # CTRIA3, CQUAD4, CTRIA6, CQUAD8, CQUAD
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type in {10, 12, 13, 14}: # CTETRA4, CHEXA8, CPENTA6, CPYRAM5
            # TODO: No idea how to get the center of the face
            #       vs. a point on a face that's not exposed
            #faces = cell.GetFaces()
            #nfaces = cell.GetNumberOfFaces()
            #for iface in range(nfaces):
                #face = cell.GetFace(iface)
                #points = face.GetPoints()
            #faces
            xyz = world_position
        elif cell_type in {24, 25, 26, 27}: # CTETRA10, CHEXA20, CPENTA15, CPYRAM13
            xyz = world_position
        elif cell_type == 3: # CBAR, CBEAM, CELASx, CDAMPx, CBUSHx
            node_xyz = np.zeros((nnodes, 3), dtype='float32')
            for ipoint in range(nnodes):
                point = points.GetPoint(ipoint)
                node_xyz[ipoint, :] = point
            xyz = node_xyz.mean(axis=0)
        elif cell_type == 21: # CBEND
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

    def get_result_by_cell_ids(self, cell_ids, icase: Optional[int]=None):
        """should handle multiple cell_ids"""
        if icase is None:
            icase = self.gui.icase_fringe
        case_key = self.gui.case_keys[icase] # int for object
        case = self.gui.result_cases[case_key]

        (obj, (i, res_name)) = case
        unused_subcase_id = obj.subcase_id
        #method = obj.get_methods(i, res_name)
        is_value, values = get_plottable_values(obj, i, res_name)
        if not is_value:
            return None

        try:
            result_values = values[cell_ids]
        except IndexError:
            msg = (
                f'values[cell_id] is out of bounds; length={len(values)}\n'
                f'result_name={res_name!r} cell_id={cell_ids} case_key={case_key}\n')
            raise IndexError(msg)

        return res_name, result_values

    def highlight_nodes_elements(self, nids=None, eids=None,
                                 representation: str='wire',
                                 model_name=None,
                                 #callback=self.on_focus_callback,
                                 force: bool=True):
        """
        Highlights a series of nodes/elements

        Parameters
        ----------
        nids / eids : int, list[int]
            the nodes/elements to highlight
        representation : str; default='wire'
            representation for elements
            allowed = {'wire', 'surface'}
        model_name : str, list[str]
            the name of the actor to highlight the nodes for

        """
        actors = []
        if eids and representation in {'wire', 'surface'}:
            actor = self.highlight_elements(self, eids, model_name=model_name)
            actors.append(actor)
        if nids and representation == 'points':
            actor = self.highlight_nodes(self, nids, model_name=model_name)
            actors.append(actor)
        return actors

    def highlight_nodes(self, nids, model_name: str='', add_actor: bool=True):
        """
        Highlights a series of nodes

        Parameters
        ----------
        nids : int, list[int]
            the nodes to apply a message to
        model_name : str, list[str]
            the name of the actor to highlight the nodes for

        """
        npoints = len(nids)
        all_nids = self.gui.get_node_ids(model_name=model_name, ids=None)
        all_xyz = self.gui.get_xyz_cid0(model_name=model_name)
        i = np.searchsorted(all_nids, nids)
        xyz = all_xyz[i, :]

        #ugrid = vtkUnstrucuturedGrid()
        points = numpy_to_vtk_points(xyz, points=None, dtype='<f', deep=1)
        ugrid = create_unstructured_point_grid(points, npoints)
        actor = self.gui.mouse_actions.create_highlighted_actor(
            ugrid, representation='points', add_actor=add_actor)
        self.gui.vtk_interactor.Render()
        return actor

    def highlight_elements(self, eids, model_name: str='', add_actor: bool=True):
        """
        Highlights a series of elements

        Parameters
        ----------
        eids : int, list[int]
            the elements to apply a message to
        model_name : str, list[str]
            the name of the actor to highlight the elements for

        """
        raise NotImplementedError('need to define active_ugrid')
        #grid = active_ugrid
        #all_eids = self.gui.get_element_ids(model_name=model_name, ids=None)
        #cell_ids = np.searchsorted(all_eids, eids)

        #ids = vtkIdTypeArray()
        #ids.SetNumberOfComponents(1)
        #for cell_id in cell_ids:
            #ids.InsertNextValue(cell_id)

        #selection_node = vtkSelectionNode()
        #selection_node.SetFieldType(vtkSelectionNode.CELL)
        #selection_node.SetContentType(vtkSelectionNode.INDICES)
        #selection_node.SetSelectionList(ids)

        #actor = self._highlight_picker_by_selection_node(
            #grid, selection_node, representation='surface', add_actor=add_actor)
        #return actor

    def _highlight_picker_by_selection_node(self, grid, selection_node,
                                            representation: str='surface',
                                            add_actor: bool=True):
        """
        helper method for:
            - _highlight_picker_cell
            #- _highlight_picker_node

        """
        selection = vtkSelection()
        selection.AddNode(selection_node)

        extract_selection = vtkExtractSelection()
        extract_selection.SetInputData(0, grid)
        extract_selection.SetInputData(1, selection)
        extract_selection.Update()

        ugrid = extract_selection.GetOutput()
        actor = self.gui.mouse_actions.create_highlighted_actor(
            ugrid, representation=representation, add_actor=add_actor)
        return actor

    def mark_nodes(self, nids: list[int], icase: int, text: str) -> None:
        """
        Marks a series of nodes with custom text labels

        Parameters
        ----------
        nids : int, list[int]
            the nodes to apply a message to
        icase : int
            the key in label_actors to slot the result into
        text : str, list[str]
            the text to display

        0 corresponds to the NodeID result
        self.mark_nodes(1, 0, 'max')
        self.mark_nodes(6, 0, 'min')
        self.mark_nodes([1, 6], 0, 'max')
        self.mark_nodes([1, 6], 0, ['max', 'min'])

        """
        if icase not in self.gui.label_actors:
            names = ', '.join(self.gui.label_actors)
            msg = f'icase={icase} not in label_actors=[{names}]'
            self.gui.log_error(msg)
            return
        slot = self.gui.label_actors[icase]
        slot.extend(
            create_marked_node_actors(self.gui, self.gui.node_ids, nids, text,
                                      self.gui.xyz_cid0))
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
                #self.gui.label_actors[icase].append(self.create_annotation(texti, xi, yi, zi))
        #self.gui.vtk_interactor.Render()

def create_marked_node_actors(gui, node_ids, nids, text, xyz_cid0):
    """
    Marks a series of nodes with custom text labels

    Parameters
    ----------
    gui : MainWindow
        ???
    node_ids : (nnodes, ) int ndarray
        all the nodes
    nids : int, list[int]
        the nodes to apply a message to
    text : str, list[str]
        the text to display
    xyz_cid0 : (nnodes, 3) float ndarray
        the xyz locations of node_ids

    0 corresponds to the NodeID result
    self.mark_nodes(1, 0, 'max')
    self.mark_nodes(6, 0, 'min')
    self.mark_nodes([1, 6], 0, 'max')
    self.mark_nodes([1, 6], 0, ['max', 'min'])

    """
    i = np.searchsorted(node_ids, nids)
    if isinstance(text, str):
        text = [text] * len(i)
    else:
        assert len(text) == len(i)

    xyz = xyz_cid0[i, :]
    slot = []
    for (xi, yi, zi), texti in zip(xyz, text):
        slot.append(create_annotation(gui, texti, xi, yi, zi))
    return slot

def create_annotation(gui: MainWindow,
                      text: str,
                      x: float, y: float, z: float) -> vtkBillboardTextActor3D:
    """
    Creates the actual annotation and appends it to slot

    Parameters
    ----------
    gui : MainWindow
       has access to the settings
    text : str
        the text to display
    x, y, z : float
        the position of the label

    Returns
    -------
    annotation : vtkBillboardTextActor3D
        the annotation object

    """
    # http://nullege.com/codes/show/src%40p%40y%40pymatgen-2.9.6%40pymatgen%40vis%40structure_vtk.py/395/vtk.vtkVectorText/python

    #self.convert_units(icase, result_value, x, y, z)

    settings = gui.settings
    text_actor = vtkBillboardTextActor3D()
    text_actor.SetPosition(x, y, z)
    text_actor.SetInput(str(text))
    text_actor.PickableOff()
    text_actor.DragableOff()
    #text_actor.SetPickable(False)

    #text_actor.SetPosition(actor.GetPosition())
    text_prop = text_actor.GetTextProperty()
    text_prop.SetFontSize(settings.annotation_size)
    text_prop.SetFontFamilyToArial()
    text_prop.BoldOn()
    text_prop.ShadowOn()
    text_prop.SetColor(settings.annotation_color)
    text_prop.SetJustificationToCentered()

    # finish adding the actor
    gui.rend.AddActor(text_actor)

    #self.label_actors[icase].append(text_actor)

    #print('added label actor %r; icase=%s' % (text, icase))
    #print(self.label_actors)

    #self.picker_textMapper.SetInput("(%.6f, %.6f, %.6f)"% pickPos)
    #camera.GetPosition()
    #camera.GetClippingRange()
    #camera.GetFocalPoint()
    return text_actor

def get_plottable_values(obj, i, res_name) -> tuple[bool, Optional[np.ndarray]]:
    """returns a vector -> fringe depending on what's available"""
    is_value = False
    fringe, vector = obj.get_fringe_vector_result(i, res_name)

    if vector is None:
        if fringe is not None:
            vector = fringe
        else:
            return is_value, None
    is_values = True
    return is_values, vector
