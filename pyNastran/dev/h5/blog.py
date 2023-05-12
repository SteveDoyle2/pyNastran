"""
a transcription of:
 - https://blog.kitware.com/developing-hdf5-readers-using-vtkpythonalgorithm/
 - https://blog.kitware.com/vtkpythonalgorithm-is-great/
"""
import os
import h5py
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk import *
from pyNastran.gui.vtk_renering_core import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor,
    vtkActor,
    vtkPolyDataMapper,
)

def write_file(data, xfreq, dirname: str):
    wdata = dsa.WrapDataObject(data)
    array = wdata.PointData['RTData']
    # Note that we flip the dimensions here because
    # VTK's order is Fortran whereas h5py writes in
    # C order. We don't want to do deep copies so we write
    # with dimensions flipped and pretend the array is
    # C order.
    array = array.reshape(wdata.GetDimensions()[::-1])
    hdf5_filename = os.path.join(dirname, 'data%d.h5' % xfreq)
    f = h5py.File(hdf5_filename, 'w')
    f.create_dataset("RTData", data=array)


def setup(dirname):
    rt = vtk.vtkRTAnalyticSource()

    for xfreq in range(60, 80):
        rt.SetXFreq(xfreq)
        rt.Update()
        write_file(rt.GetOutput(), xfreq)


def run_source_class():
    asdf
class MyPointSource(vtk.vtkProgrammableSource):
    def __init__(self):
        vtk.vtkProgrammableSource(self)
        def genPoints():
            points = vtk.vtkPoints()
            output = self.GetPolyDataOutput()
            output.SetPoints(points)
            # Left the points in a *. txt file , and read them from here
            n = 5
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        points.InsertNextPoint(i, j, k)
        self.SetExecuteMethod(genPoints)

def get_point_source():
    #return MyPointSource()
    pointSource = vtk.vtkProgrammableSource() # Source
    def genPoints():
        points = vtk.vtkPoints()
        output = pointSource.GetPolyDataOutput()
        output.SetPoints(points)
        # Left the points in a *. txt file , and read them from here
        n = 3
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    points.InsertNextPoint(i, j, k)
    pointSource.SetExecuteMethod(genPoints)
    return pointSource

def run_source():
    """https://www.creatis.insa-lyon.fr/~davila/bbtk/Software/new/doc/VTK_Documentation/report-about-vtk-2007-02.pdf"""
    pointSource = get_point_source()
    popSplatter = vtk.vtkGaussianSplatter() # Splatter :
    popSplatter.SetSampleDimensions(300, 300, 300) # points -> gauss sphere

    popSplatter.SetInputConnection(pointSource.GetOutputPort())
    # f(x) = scale * exp ( expF *(( r/ Radius )^2) )
    popSplatter.SetRadius(0.1)
    popSplatter.SetExponentFactor(1.0)
    popSplatter.SetScaleFactor(1.0)

    cast = vtk.vtkImageCast() # uchar conversion
    cast.SetInputConnection(popSplatter.GetOutputPort()) # ( required by VRayCast )
    cast.SetOutputScalarTypeToUnsignedChar()

    alphaTF = vtk.vtkPiecewiseFunction() # Opacity (A-TF)
    alphaTF.AddPoint(0, 0)
    alphaTF.AddPoint(1, 0.01)
    alphaTF.AddPoint(255, 0.5)
    colorTF = vtk.vtkColorTransferFunction() # Color (RGB -TF)
    colorTF.AddRGBPoint(0, 1, 0.4 , 0)
    colorTF.AddRGBPoint(255, 1, 0.0 , 0)

    volumeProperty = vtk.vtkVolumeProperty() # Property
    volumeProperty.SetColor(colorTF)
    volumeProperty.SetScalarOpacity(alphaTF)
    volumeProperty.SetInterpolationTypeToLinear()
    #compositeFunction = vtkVolumeRayCastCompositeFunction()# Mapper

    volumeMapper = vtk.vtkFixedPointVolumeRayCastMapper() # ( Software )
    #volumeMapper.SetVolumeRayCastFunction(compositeFunction)
    volumeMapper.SetBlendModeToComposite()
    volumeMapper.SetInputConnection(cast.GetOutputPort())

    volume = vtk.vtkVolume() # Volume = Mapper + Property
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    ren = vtkRenderer() # Renderer
    ren.AddVolume(volume)
    ren.ResetCamera()
    ren.GetActiveCamera().ParallelProjectionOn()

    renwin = vtkRenderWindow() # Window
    renwin.AddRenderer (ren)
    renwin.SetSize(300, 300)
    renwin.SetDesiredUpdateRate(1.0)

    iren = vtkRenderWindowInteractor() # Interactor
    iren.SetRenderWindow(renwin)
    iren.Start()

class HDF5Source(VTKPythonAlgorithmBase):
    def __init__(self):
        """
        1 input port + 1 output port == your common VTK filter.
        0 input port + 1 output port == your common VTK source.
        1 input port + 0 output port == your common VTK sink (writer, mapper etc.).
        """
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=0,
            nOutputPorts=1,
            outputType='vtkImageData')
        self.__FileName = ""

    def RequestData(self, request, inInfo, outInfo):
        f = h5py.File(self.__FileName, 'r')
        data = f['RTData'][:]
        output = dsa.WrapDataObject(vtk.vtkImageData.GetData(outInfo))
        # Note that we flip the dimensions here because
        # VTK's order is Fortran whereas h5py writes in
        # C order.
        output.SetDimensions(data.shape[::-1])
        output.PointData.append(data.ravel(), 'RTData')
        output.PointData.SetActiveScalars('RTData')
        return 1

    def SetFileName(self, fname):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname

    def GetFileName(self):
        return self.__FileName

    def RequestInformation(self, request, inInfo, outInfo):
        """
        As I discussed previously, RequestInformation provides meta-data downstream.
        This meta-data is most of the time lightweight. In this example, we used
        f[‘RTData’].shape to read extent meta-data from the HDF5 file. This does not
        read any heavyweight data. Later on, we will see other examples of meta-data
        that is provided during RequestInformation.
        """
        f = h5py.File(self.__FileName, 'r')
        # Note that we flip the shape because VTK is Fortran order
        # whereas h5py reads in C order. When writing we pretend that the
        # data was C order so we have to flip the extents/dimensions.
        dims = f['RTData'].shape[::-1]
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
            (0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1), 6)
        return 1

class RequestSubset(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=1, inputType='vtkImageData',
            nOutputPorts=1, outputType='vtkImageData')
        self.__UpdateExtent = None

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(), \
            self.__UpdateExtent, 6)
        return 1

    def RequestUpdateExtent(self, request, inInfo, outInfo):
        if self.__UpdateExtent is not None:
            info = inInfo[0].GetInformationObject(0)
            info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_EXTENT(), \
                self.__UpdateExtent, 6)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        inp = vtk.vtkImageData.GetData(inInfo[0])
        opt = vtk.vtkImageData.GetData(outInfo)
        opt.ShallowCopy(inp)
        return 1

    def SetUpdateExtent(self, ue):
        if ue != self.__UpdateExtent:
            self.Modified()
            self.__UpdateExtent = ue

    def GetUpdateExtent(self):
        return self.__UpdateExtent

def run_vtk(dirname: str):
    alg = HDF5Source()
    hdf5_filename = os.path.join(dirname, 'data60.h5')
    alg.SetFileName(hdf5_filename)

    cf = vtk.vtkContourFilter()
    cf.SetInputConnection(alg.GetOutputPort())
    cf.SetValue(0, 200)

    m = vtkPolyDataMapper()
    m.SetInputConnection(cf.GetOutputPort())

    a = vtkActor()
    a.SetMapper(m)

    ren = vtkRenderer()
    ren.AddActor(a)

    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(600, 600)

    for xfreq in range(60, 80):
        hdf5_filename = os.path.join(dirname, 'data%d.h5' % xfreq)
        alg.SetFileName(hdf5_filename)
        renWin.Render()
        import time
        time.sleep(0.1)

def run_vtk2(dirname: str):
    alg = HDF5Source()
    hdf5_filename = os.path.join(dirname, 'data60.h5')
    alg.SetFileName(hdf5_filename)

    rs = RequestSubset()
    rs.SetInputConnection(alg.GetOutputPort())
    rs.SetUpdateExtent((5, 10, 5, 10, 0, 20))

    cf = vtk.vtkContourFilter()
    cf.SetInputConnection(rs.GetOutputPort())
    cf.SetValue(0, 200)

    m = vtkPolyDataMapper()
    m.SetInputConnection(cf.GetOutputPort())

    a = vtkActor()
    a.SetMapper(m)

    ren = vtkRenderer()
    ren.AddActor(a)

    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(300, 300)
    renWin.Render()

    for xfreq in range(60, 80):
        hdf5_filename = os.path.join(dirname, 'data%d.h5' % xfreq)
        alg.SetFileName(hdf5_filename)
        renWin.Render()
    asdf

def main():
    dirname = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\utils\hdf5'
    #setup(dirname)
    #run_vtk(dirname)
    run_vtk2(dirname)


if __name__ == '__main__':
    #main()
    run_source()
