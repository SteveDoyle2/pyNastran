"""
defines the AreaPick class

http://www.paraview.org/Wiki/Selection_Implementation_in_VTK_and_ParaView_III
http://ruby-vtk.rubyforge.org/svn/trunk/VTK/Rendering/Testing/Cxx/TestAreaSelections.cxx
http://vtk.1045678.n5.nabble.com/picking-objects-from-a-subset-of-a-grid-td5143877.html
http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/HighlightSelectedPoints
http://www.vtk.org/doc/nightly/html/classvtkExtractSelectedFrustum.html
http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGridAlgorithm.html
http://www.vtk.org/doc/nightly/html/classvtkExtractCells.html
http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/HighlightSelection
http://public.kitware.com/pipermail/vtkusers/2012-January/072046.html
http://vtk.1045678.n5.nabble.com/Getting-the-original-cell-id-s-from-vtkExtractUnstructuredGrid-td1239667.html
"""
from __future__ import print_function, division
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from pyNastran.bdf.utils import write_patran_syntax_dict

#left_button_down=self._zoom_picker,
#left_button_up=self._zoom_picker,
#right_button_down=self._zoom_reset,

#class AreaPickStyle(vtk.vtkInteractorStyleRubberBandPick):
    #"""Custom Rubber Band Picker"""
    #def __init__(self, parent=None):
        #"""creates the AreaPickStyle instance"""
        #pass
        ##super(AreaPickStyle, self).__init__()

        ## for vtk.vtkInteractorStyleRubberBandZoom
        #self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        #self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        ##self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        #self.parent = parent
        #self.area_pick_button = self.parent.actions['area_pick']
        #self.picker_points = []
        #self.parent.area_picker.SetRenderer(self.parent.rend)

class AreaPickStyle(vtk.vtkInteractorStyleRubberBandZoom):
    def __init__(self, parent=None):
        """creates the AreaPickStyle instance"""
        pass
        # for vtk.vtkInteractorStyleRubberBandZoom
        self.AddObserver("LeftButtonPressEvent", self._left_button_press_event)
        self.AddObserver("LeftButtonReleaseEvent", self._left_button_release_event)
        self.AddObserver("RightButtonPressEvent", self.right_button_press_event)
        self.parent = parent
        self.area_pick_button = self.parent.actions['area_pick']
        self.picker_points = []
        self.parent.area_picker.SetRenderer(self.parent.rend)

    #def leftButtonPressEvent(self, obj, event):
        #pass

    def _left_button_press_event(self, obj, event):
        """
        gets the first point
        """
        print('area_picker - left_button_press_event')
        self.OnLeftButtonDown()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        self.picker_points.append((pixel_x, pixel_y))

    def _left_button_release_event(self, obj, event):
        """
        gets the second point and zooms

        TODO: doesn't handle panning of the camera to center the image
              with respect to the selected limits
        """
        #self.OnLeftButtonUp()
        pixel_x, pixel_y = self.parent.vtk_interactor.GetEventPosition()
        #selector = vtk.vtkVisibleCellSelector()
        #selector


        self.picker_points.append((pixel_x, pixel_y))

        area_picker = self.parent.area_picker
        area_picker.Pick()

        print(self.picker_points)
        if len(self.picker_points) == 2:
            p1x, p1y = self.picker_points[0]
            p2x, p2y = self.picker_points[1]
            self.picker_points = []
            dx = abs(p1x - p2x)
            dy = abs(p1y - p2y)
            xmin = min(p1x, p2x)
            ymin = min(p1y, p2y)
            xmax = max(p1x, p2x)
            ymax = max(p1y, p2y)
            #self.area_picker.SetPickCoords(xmin, ymin, xmax, ymax)
            area_picker.AreaPick(xmin, ymin, xmax, ymax, self.parent.rend)
            print(self.picker_points)
            #print('_rotation_center_cell_picker', cell_id)

            self.picker_points = []
            if dx > 0 and dy > 0:
                if 0:
                    # works very inefficiently
                    cell_ids = set([])
                    picker = self.cell_picker
                    nx = max((xmax - xmin) // 100, 1)
                    ny = max((ymax - ymin) // 100, 1)
                    #print('nx=%s ny=%s' % (nx, ny))
                    for xi in range(xmin, xmax, nx):
                        for yi in range(ymin, ymax, ny):
                            picker.Pick(xi, yi, 0, self.rend)
                            cell_id = picker.GetCellId()
                            if cell_id >= 0:
                                cell_ids.add(int(cell_id))
                else:
                    data = area_picker.GetDataSet()
                    pick_list = area_picker.GetPickList()

                    props = area_picker.GetProp3Ds()
                    #print('props=', props)
                    #print('data=', data)

                    cells = data.GetCells()  # vtkCellArray
                    #cells.GetCell(int, vtkIdList)

                    cell_data = data.GetCellData() # vtkCellData

                    #data.GetIdsOfCellsOfType(int, vtkIdTypeArray)
                    for i in range(data.GetNumberOfCells()):
                        c = data.GetCell(i)  # vtkTetra
                        #cell_ids
                        pass
                    frustrum = area_picker.GetFrustum() # vtkPlanes
                    grid = self.parent.grid

                    selection_cells = vtk.vtkSelectionNode()
                    selection_cells.SetFieldType(vtk.vtkSelectionNode.CELL)
                    selection_cells.SetContentType(vtk.vtkSelectionNode.INDICES)


                    #extract_cells = vtk.vtkExtractCells()

                    extract_ids = vtk.vtkExtractSelectedIds()
                    extract_ids.AddInputData(grid)



                    idsname = "Ids"
                    ids = vtk.vtkIdFilter()
                    #ids.SetInputConnection(grid.GetOutputPort()) #  fails
                    ids.SetInputData(grid)
                    ids.PointIdsOn()
                    ids.CellIdsOn()
                    ids.FieldDataOn()
                    ids.SetIdsArrayName(idsname)

                    # instead of extract surface, put your operation here ...
                    #surface = vtkDataSetSurfaceFilter()
                    #surface.SetInputConnection(ids.GetOutputPort())
                    #surface.Update()

                    #pointIdArray = surface.GetOutput().GetPointData().GetArray(idsname)
                    #cellIdArray = surface.GetOutput().GetCellData().GetArray(idsname)


                    selected_frustrum = vtk.vtkExtractSelectedFrustum()
                    selected_frustrum.SetFrustum(frustrum)
                    #selected_frustrum.SetInputData(grid)  $ was on
                    selected_frustrum.PreserveTopologyOff()
                    #maskPts.SetInputData(src)
                    #selected_frustrum.SetSelectionConnection(selection_cells.GetOutputPort())
                    selected_frustrum.Update()
                    #selected_frustrum.RequestData()
                    #frustrum_cells = selected_frustrum.GetContainingCells()
                    #sdata = vtk.vtkDataObject()
                    #extractor.GetOutput()
                    selected_frustrum.SetInputConnection(ids.GetOutputPort())  # was grid?
                    selected_frustrum.Update()

                    ugrid = selected_frustrum.GetOutput()
                    ugrid_cell_data = ugrid.GetCellData()
                    cell_data_array = ugrid_cell_data.GetArray('Ids')
                    cell_ids = vtk_to_numpy(cell_data_array)
                    #print('eids = ', eids)

                    if 0:
                        #cells = ugrid.GetCells()
                        idtype = cells.GetData()
                        algorithm = selected_frustrum.GetOutputPort()
                        extract_ids.SetSelectionConnection(algorithm)

                        #idtype.GetValue(8)
                        numpy_array = vtk_to_numpy(idtype)
                        unumpy_array = np.unique(numpy_array)
                        eids = self.parent.element_ids[unumpy_array]

                        if 1:
                            extract_geometry = vtk.vtkExtractGeometry()
                            extract_geometry.ExtractInsideOn()
                            extract_geometry.SetImplicitFunction(frustrum)

                            if vtk.VTK_MAJOR_VERSION <= 5:
                                extract_geometry.SetInput(grid)
                            else:
                                extract_geometry.SetInputData(grid)
                                extract_geometry.Update()

                            data = vtk.vtkDataObject()
                            #data.SetCells(grid.GetCells())
                            extract_geometry.SetOutput(data)
                            #extract_geometry.SetInputConnection(0, grid.GetProducerPort())
                            #extract_geometry.SetInputConnection(1, selection.GetProducerPort())

                        extract_grid = vtk.vtkExtractUnstructuredGrid()
                        extract_grid.SetInputData(grid)
                        extract_grid.SetInputConnection(connect2.GetOutputPort())
                        extract_grid.CellClippingOn()

                    #print('pick_list=', pick_list)
                    if 0:
                        # used to be active and didn't crash
                        for i in range(props.GetNumberOfItems()):
                            prop = props.GetNextProp3D()
                            print("Picked prop: ", prop)

                eids = self.parent.element_ids[cell_ids]
                eids_str = write_patran_syntax_dict({'Elem' : eids})
                self.parent.log_info('cell_ids = %s' % eids_str)

                #area_pick_button = self.actions['area_pick']
                self.area_pick_button.setChecked(False)
                self.parent.setup_mouse_buttons(mode='default')
            self.parent.vtk_interactor.Render()
        self.picker_points = []

    def right_button_press_event(self, obj, event):
        """cancels the button"""
        self.area_pick_button.setChecked(False)
        self.parent.setup_mouse_buttons(mode='default')
        self.parent.vtk_interactor.Render()
