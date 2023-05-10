try:
    from vtkmodules.vtkCommonDataModel import (
        vtkLine,
        vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode)
except ImportError:
    from vtk import (
        vtkLine,
        vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, )
