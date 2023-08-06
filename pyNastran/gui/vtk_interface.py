try:
    from vtkmodules.vtkCommonDataModel import (
        vtkVertex,
        vtkLine, vtkQuadraticEdge,
        vtkTriangle, vtkQuadraticTriangle,
        vtkQuad, vtkQuadraticQuad, vtkBiQuadraticQuad,
        vtkTetra, vtkQuadraticTetra,
        vtkWedge, vtkQuadraticWedge,
        vtkHexahedron,vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

except ImportError:
    print('error vtk_interface')
    from vtk import (
        vtkVertex,
        vtkLine, vtkQuadraticEdge,
        vtkTriangle, vtkQuadraticTriangle,
        vtkQuad, vtkQuadraticQuad, vtkBiQuadraticQuad,
        vtkTetra, vtkQuadraticTetra,
        vtkWedge, vtkQuadraticWedge,
        vtkHexahedron,vtkQuadraticHexahedron,
        vtkPyramid, vtkQuadraticPyramid,
        vtkUnstructuredGrid, vtkMultiBlockDataSet,
        vtkSelectionNode, vtkCellArray,
        vtkPolyData,
        VTK_POLYHEDRON,
    )

