import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase

class vtkRequestSubset(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=1, inputType='vtkImageData',
            nOutputPorts=1, outputType='vtkImageData')
        self.__UpdateExtent = None

    def RequestInformation(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
                 self.__UpdateExtent, 6)
        return 1

    def RequestUpdateExtent(self, request, inInfo, outInfo):
        if self.__UpdateExtent is not None:
            info = inInfo[0].GetInformationObject(0)
            info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_EXTENT(),
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
