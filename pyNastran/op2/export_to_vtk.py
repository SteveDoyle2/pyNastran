from __future__ import print_function, unicode_literals
import struct
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from numpy import zeros


def export_to_vtk(bdf_filename, op2_filename, vtk_filename):
    with open(vtk_filename, 'w') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write('num_scheme_test_steve\n')
        f.write('BINARY\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')

        model = BDF()
        model.read_bdf(bdf_filename)
        #op2 = OP2()
        #op2.read_op2(op2_filename)

        out = bdf.get_card_ids_by_card_types()
        grids = sorted(out['GRID'])
        spoints = sorted(out['SPOINT'])
        epoints = sorted(out['EPOINT'])
        ngrid = len(grids)
        nspoint = len(spoints)
        nepoint = len(epoints)
        nnodes = ngrid + nspoint + nepoint

        ncrod = len(out['CROD'])
        nconrod = len(out['CONROD'])
        nctube = len(out['CTUBE'])
        nlines = ncrod + nconrod + nctube

        nctria3 = len(out['CTRIA3'])
        ncquad4 = len(out['CQUAD4'])
        nshells = nctria3 + ncquad4
        nproperties = nelements

        #nsolids = nctetra4 + nctetra10

        # SPOINT & EPOINT are implicitly defined
        xyz_cid0 = zeros((nnodes, 3), dtype='float32')
        for nid in grids:
            xyz_cid[i, :] = model.nodes[nid].get_position()
        nids[:i] = grids
        if nspoints:
            nids[i:i+nspoint] = spoints
        if nepoints:
            nids[i+nspoint:] = epoints

        nelements = nlines + nshells + nsolids

        f.write('POINTS %i float\n' % nnodes)
        f.write(struct.pack(fmt, nids.ravel()))
        f.write('CELLS %i %i\n' % (nnodes, nelements))

        f.write('NodeID %i float\n' % nnodes)
        f.write(struct.pack(fmt, nids.ravel()))

        fmt = b'%si' % nelements
        if nelements:
            f.write('ElementID %i float\n' % nelements)
            f.write(struct.pack(fmt, eids.ravel()))
        if nproperties:
            f.write('PropertyID %i float\n' % nproperties)
            f.write(struct.pack(fmt, pids.ravel()))
        if nmaterials:
            f.write('MaterialID %i float\n' % nmaterials)
            f.write(struct.pack(fmt, mids.ravel()))

        nodal_cases = [op2.eigenvectors, op2.displacements, op2.velocities, op2.accelerations]
        for cases in nodal_cases:
            for isubcase, case in iteritems(sorted(cases)):
                pass
        #CELLS 217 1039

