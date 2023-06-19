try:
    from vtkmodules.vtkCommonDataModel import (
        vtkLine,
        vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
        vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
        vtkQuadraticWedge, vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

except ImportError:
    from vtk import (
        vtkLine,
        vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
        vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
        vtkQuadraticWedge, vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

