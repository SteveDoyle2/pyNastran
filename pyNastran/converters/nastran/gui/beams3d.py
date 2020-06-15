"""creates 3d beams"""
from typing import Any
import vtk

def faces_to_element_facelist(faces: Any, node0: int) -> vtk.vtkIdList:
    """creates a series of faces for the custom elements"""
    face_idlist = vtk.vtkIdList()

    nfaces = len(faces)
    face_idlist.InsertNextId(nfaces) # Number faces that make up the cell.
    for face in faces: # Loop over all the faces
        #print(face)
        face_idlist.InsertNextId(len(face)) # Number of points in face

        # Insert the pointIds for the face
        #for i in face:
            #face_idlist.InsertNextId(i + node0)
        [face_idlist.InsertNextId(i + node0) for i in face]
    return face_idlist

def get_bar_type(ptype: str, pid_ref):
    """helper method for _get_bar_yz_arrays"""
    if ptype in ['PBAR', 'PBEAM']:
        bar_type = 'bar'
    #if ptype == 'PBAR':
        #bar_type = 'pbar'
    #elif ptype == 'PBEAM':
        #bar_type = 'pbeam'
    elif ptype in ['PBARL', 'PBEAML']:
        bar_type = pid_ref.Type
    elif ptype == 'PBCOMP':
        bar_type = 'pbcomp'
    else:  # pragma: no cover
        raise NotImplementedError(pid_ref)
    return bar_type
