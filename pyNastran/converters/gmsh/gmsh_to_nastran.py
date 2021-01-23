from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.bdf.bdf import BDF
if TYPE_CHECKING:  # pragma: no cover
    from pygmsh.geo import Geometry

ELEMENT_TYPE_MAPPER = {
    'line' : ['CROD', 'CTUBE', 'CBAR'],
    'triangle' : ['CTRIA3', 'CTRIAR'],
    'quad' : ['CQUAD4', 'CQUADR'],
    'tetra' : ['CTETRA'],
    'hexahedron' : ['CHEXA'],
    'vertex' : ['CONM2', 'CELAS1'],
}

def _load_grids(mesh: Geometry, model: BDF):
    nid = 1
    for xyz in mesh.points:
        model.add_grid(nid, xyz, cp=0, cd=0, ps='', seid=0, comment='')
        nid += 1

def write_nastran(mesh: Geometry, nastran_filename: str,
                  element_type_mapper=None) -> BDF:
    """converts a meshio Geometry/mesh object to a pyNastran BDF"""
    if element_type_mapper is None:
        element_type_mapper = {
            'line' : 'CBAR',
            'triangle' : 'CTRIA3',
            'quad' : 'CQUAD4',
            'tetra' : 'CTETRA',
            'hexahedron' : 'CHEXA',
            'vertex' : 'CONM2',
        }

    model = BDF()
    eid = 1
    pid = 1
    _load_grids(mesh, model)

    for cells in mesh.cells:
        eid_old = eid
        cell_type = cells.type
        element_type = element_type_mapper[cell_type]
        allowed_element_types = ELEMENT_TYPE_MAPPER[cell_type]
        print(cell_type, element_type)
        if element_type not in allowed_element_types:
            raise RuntimeError(f'element_type={element_type} not in '
                               f'{allowed_element_types}')
        #cells_2d = {"triangle", "quad"}
        #cells_3d = {
            #"tetra",
            #"hexahedron",
            #"wedge",
            #"pyramid",
            #"penta_prism",
            #"hexa_prism",
        #}
        if cell_type == 'line':
            if element_type == 'CROD':
                for nids in cells.data + 1:
                    model.add_crod(eid, pid, nids, comment='')
                    eid += 1
            elif element_type == 'CTUBE':
                for nids in cells.data + 1:
                    model.add_ctube(eid, pid, nids, comment='')
                    eid += 1
            elif element_type == 'CBAR':
                x_vector = [0., 0., 1.]
                g0 = None
                for nids in cells.data + 1:
                    model.add_cbar(eid, pid, nids, x_vector, g0, offt='GGG',
                                   pa=0, pb=0, wa=None, wb=None, comment='')
                    eid += 1
                del x_vector, g0

        # shells
        elif cell_type == 'triangle':
            if element_type == 'CTRIA3':
                for nids in cells.data + 1:
                    model.add_ctria3(eid, pid, nids, zoffset=0., theta_mcid=0.0, tflag=0,
                                     T1=None, T2=None, T3=None, comment='')
                    eid += 1
            elif element_type == 'CTRIAR':
                for nids in cells.data + 1:
                    model.add_ctriar(eid, pid, nids, theta_mcid=0.0, zoffset=0.0, tflag=0,
                                     T1=None, T2=None, T3=None, comment='')
                    eid += 1

        elif cell_type == 'quad':
            if element_type == 'CQUAD4':
                for nids in cells.data + 1:
                    model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                     T1=None, T2=None, T3=None, T4=None, comment='')
                    eid += 1
            elif element_type == 'CQUADR':
                for nids in cells.data + 1:
                    model.add_cquadr(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                     T1=None, T2=None, T3=None, T4=None, comment='')
                    eid += 1

        # solids
        elif cell_type == 'tetra':
            for nids in cells.data + 1:
                model.add_ctetra(eid, pid, nids, comment='')
                eid += 1
        elif cell_type == 'hexahedron':
            for nids in cells.data + 1:
                model.add_chexa(eid, pid, nids, comment='')
                eid += 1
        #elif cell_type == 'wedge':
            #for nids in cells.data + 1:
                #model.add_cpenta(eid, pid, nids, comment='')
                #eid += 1
        #elif cell_type == 'pyramid':
            #for nids in cells.data + 1:
                #model.add_cpyram(eid, pid, nids, comment='')
                #eid += 1

        elif cell_type == 'vertex':
            if element_type == 'CONM2':
                mass = 0.0
                for nids in cells.data + 1:
                    nid = nids[0]
                    model.add_conm2(eid, nid, mass, cid=0, X=None, I=None, comment='')
                    eid += 1
                del mass
            elif element_type == 'CELAS1':
                for nids in cells.data + 1:
                    nid = nids[0]
                    nids = [nid]
                    model.add_celas1(eid, pid, nids, c1=0, c2=0, comment='')
                    eid += 1
            del nid


        else:
            raise NotImplementedError(cells.type)
        assert eid != eid_old, f'element_type={element_type} is not supported'
        pid += 1
    assert eid != 1, eid
    model.write_bdf(nastran_filename)

def main():
    import sys
    gmsh_filename = sys.argv[1]
    nastran_filename = sys.argv[2]
    print(gmsh_filename)
    make_geometry(gmsh_filename, nastran_filename)

def make_geometry(mesh_filename: str, nastran_filename: str):
    import meshio
    geom = meshio.read(mesh_filename, file_format=None)
    geom.write('spike.bdf')
    model = write_nastran(geom, nastran_filename)

if __name__ == '__main__':
    main()
