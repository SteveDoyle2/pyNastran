try:
    from vtkmodules.vtkCommonDataModel import (
        vtkVertex, vtkLine,
        vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
        vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
        vtkQuadraticWedge, vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkQuadraticEdge, vtkBiQuadraticQuad,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

except ImportError:
    from vtk import (
        vtkVertex, vtkLine,
        vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
        vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
        vtkQuadraticWedge, vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkQuadraticEdge, vtkBiQuadraticQuad,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

