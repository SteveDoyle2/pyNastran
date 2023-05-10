import h5py
import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, vtkMultiBlockDataSet


class vtkH5NastranReader(VTKPythonAlgorithmBase):
    """https://www.paraview.org/Wiki/Python_Programmable_Filter"""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=0,
            nOutputPorts=1,
            outputType='vtkMultiBlockDataSet')
        self.__FileName = ""

    def SetFileName(self, fname):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname
            assert fname.lower().endswith('.h5'), fname

    def GetFileName(self):
        return self.__FileName

    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtk.vtkMultiBlockDataSet.GetData(outInfo))
        hdf5_filename = self.__FileName

        geom_model, vtk_ugrid, root, alt_grids = get_paraview_nastran_ugrid(hdf5_filename)
        #ugrid = vtkUnstructuredGrid()
        idx = 0
        output.SetBlock(idx, ugrid)
        ug = dsa.WrapDataObject(ugrid)
        return
        f = h5py.File(self.__FileName, 'r')
        idx = 0
        for grp_name in f:
            ug = vtkUnstructuredGrid()
            output.SetBlock(idx, ug)
            idx += 1
            ug = dsa.WrapDataObject(ug)
            fill_paraview_vtk_unstructured_grid(geom_model, vtk_ugrid,
                                                add_property=True, add_material=True)
            #grp = f[grp_name]
            #cells = grp['cells'][:]
            #locations = grp['cell_locations'][:]
            #types = grp['cell_types'][:]
            #ug.SetCells(types, locations, cells)
            pts = grp['points'][:]
            ug.Points = pts
            pt_arrays = grp['point_data']
            for pt_array in pt_arrays:
                array = pt_arrays[pt_array][:]
                ug.PointData.append(array, pt_array)
        return 1
