import os
from pyNastran.bdf.bdf import read_bdf
from cpylog import SimpleLogger

def nastran_to_fluent(nastran_filename: str, fluent_filename: str, log: SimpleLogger=None):
    model = read_bdf(nastran_filename, log=log)
    vrt_filename = os.path.splitext(fluent_filename)[0] + '.vrt'
    cel_filename = os.path.splitext(fluent_filename)[0] + '.cel'
    daten_filename = os.path.splitext(fluent_filename)[0] + '.daten'
    with open(vrt_filename, 'w') as vrt_file:
        vrt_file.write('PROSTAR_VERTEX\n')
        vrt_file.write('asdf\n')
        for nid, node in sorted(model.nodes.items()):
            x, y, z = node.get_position()
            vrt_file.write(f'{nid} {x} {y} {z}\n')

    eids = []
    thicknesses = []
    with open(cel_filename, 'w') as cel_file:
        cel_file.write('PROSTAR_CELL\n')
        cel_file.write('asdf\n')
        for eid, elem in sorted(model.elements.items()):
            if elem.type == 'CTRIA3':
                nnodes = 3
                pid = elem.pid
                n1, n2, n3 = elem.nodes
                cel_file.write(f'{eid}          3          {nnodes}         {pid}         4\n'
                               f'{eid}          {n1}          {n2}          {n3}\n')

            elif elem.type == 'CQUAD4':
                nnodes = 4
                pid = elem.pid
                n1, n2, n3, n4 = elem.nodes
                cel_file.write(f'{eid}          3          {nnodes}         {pid}         4\n'
                               f'{eid}          {n1}          {n2}          {n3}          {n4}\n')
            else:
                continue
            eids.append(eid)
            thicknesses.append(elem.Thickness())

    with open(daten_filename, 'w') as daten_file:
        daten_file.write('# Shell Id, Thickness\n')
        for eid, thickness in zip(eids, thicknesses):
            daten_file.write(f'{eid} {thickness}\n')
