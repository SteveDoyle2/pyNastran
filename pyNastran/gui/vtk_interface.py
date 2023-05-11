try:
    from vtkmodules.vtkCommonDataModel import (
        vtkLine,
        vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode,
        VTK_POLYHEDRON)
except ImportError:
    from vtk import (
        vtkLine,
        vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode,
        VTK_POLYHEDRON)
