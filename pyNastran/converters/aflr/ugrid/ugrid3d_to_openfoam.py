import os
from copy import deepcopy
from collections import defaultdict

from numpy import zeros, unique, where, argsort, searchsorted, allclose, array

from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.aflr.surf.surf_reader import TagReader

def write_foam(ugrid, foam_filename, tag_filename):
    """writes an OpenFOAM file"""
    dirname = os.path.dirname(foam_filename)
    unused_base = os.path.splitext(foam_filename)[0]

    points_filename = os.path.join(dirname, 'points')
    boundary_filename = os.path.join(dirname, 'boundary')
    faces_filename = os.path.join(dirname, 'faces')
    #neighbor_filename = os.path.join(dirname, 'neighbor')
    #owner_filename = os.path.join(dirname, 'owner')

    # boundary
    # 1. get array of unique properties
    # 2. loop over in sorted order
    # 3.   find where values are equal to the property id in ctrias
    #      and cquads (surface elements)
    # 4.   find the min/max

    # points
    # 1. loop over points and write

    # faces
    # 1. loop over tets/penta5s/penta/chexas:  DO CPENTA5s last!
    # 2.    element.faces
    # 3.    only write it once!  (Sort the data for checking, but
    #       write out with proper connectivity)

    #f.write('CEND\n')
    #f.write('BEGIN BULK\n')
    #f.write('PSHELL,1,1, 0.1\n')
    #f.write('MAT1, 1, 1.0e7,, 0.3\n')

    #pids = self.pids
    #mid = 1
    #points_filename = foam_filename  # remove...
    _write_points(ugrid, points_filename)
    _write_boundary(ugrid, boundary_filename + '2', tag_filename)
    _write_faces(ugrid, faces_filename + '2')

def _write_points(ugrid, points_filename):
    """writes an OpenFOAM points file"""
    with open(points_filename, 'wb') as points_file:
        nnodes = ugrid.nodes.shape[0]

        points_file.write(
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\\\    /   O peration     | Version:  1.7.1                                 |\n'
            '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n'
            '|    \\\\/     M anipulation  |                                                 |\n'
            '\\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     2.0;\n'
            '    format      ascii;\n'
            '    class       vectorField;\n'
            '    location    "constant/polyMesh";\n'
            '    object      points;\n'
            '}\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /\n'
        )


        points_file.write('\n\n')
        points_file.write('%i\n' % (nnodes))
        points_file.write('(\n')
        for node in enumerate(ugrid.nodes):
            points_file.write(('    (%-12s %-12s %-12s)\n') % (node[0], node[1], node[2]))
        points_file.write(')\n')

def _write_boundary(ugrid, boundary_filename, tag_filename):
    """writes an OpenFOAM boundary file"""
    with open(boundary_filename, 'wb') as boundary_file:
        boundary_file.write('\n\n')
        #f.write('%i\n' % (nnodes))
        boundary_file.write('(\n')

        uboundaries = unique(ugrid.pids)
        nboundaries = len(uboundaries)
        boundary_file.write('%i\n' % nboundaries)
        boundary_file.write('(\n')

        tagger = TagReader()
        tag_data = tagger.read_tag_filename(tag_filename)


        isort = argsort(ugrid.pids)
        ugrid.pids.sort()
        #print(isort)
        pids = ugrid.pids
        for iboundary in uboundaries:
            data = tag_data[iboundary]
            #name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
            #is_trans, is_delete, bl_spacing, bl_thickness, nlayers = data
            name = data[0]

            i = where(iboundary == pids)[0]
            nfaces = i.max() - i.min() + 1
            startface = i.min()

            if len(i) != nfaces:
                msg = 'The data is unsorted...len(i)=%s nfaces=%s' % (len(i), nfaces)
                raise RuntimeError(msg)
            boundary_file.write('    %s\n' % name)
            boundary_file.write('    {\n')
            boundary_file.write('        type patch;\n')
            #f.write('        type 1(wall);\n')  # create a new group
            boundary_file.write('        nFaces %i;\n' % nfaces)
            boundary_file.write('        startFace %i;\n' % startface)
            boundary_file.write('    }\n')
        boundary_file.write(')\n')
    ugrid.isort = isort

def _write_faces(ugrid, faces_filename):
    """writes an OpenFOAM faces file"""
    nhexas = ugrid.hexas.shape[0]
    npenta6s = ugrid.penta6s.shape[0]
    npenta5s = ugrid.penta5s.shape[0]
    ntets = ugrid.tets.shape[0]

    nquad_faces = nhexas * 6 + npenta5s + npenta6s * 3
    ntri_faces = ntets * 4 + npenta5s * 4 + npenta6s * 2
    nfaces = ntri_faces + nquad_faces
    assert nfaces > 0, nfaces

    #tri_face_to_eids = ones((nt, 2), dtype='int32')
    tri_face_to_eids = defaultdict(list)

    #quad_face_to_eids = ones((nq, 2), dtype='int32')
    quad_face_to_eids = defaultdict(list)

    tri_faces = zeros((ntri_faces, 3), dtype='int32')
    quad_faces = zeros((nquad_faces, 4), dtype='int32')

    with open(faces_filename, 'wb') as faces_file:
        faces_file.write('\n\n')
        #faces_file.write('%i\n' % (nnodes))
        faces_file.write('(\n')

        it_start = {}
        iq_start = {}
        min_eids = {}
        it = 0
        iq = 0
        eid = 1
        it_start[1] = it
        iq_start[1] = iq
        min_eids[eid] = ugrid.tets
        for element in ugrid.tets - 1:
            (n1, n2, n3, n4) = element
            face1 = [n3, n2, n1]
            face2 = [n1, n2, n4]
            face3 = [n4, n3, n1]
            face4 = [n2, n3, n4]

            tri_faces[it, :] = face1
            tri_faces[it+1, :] = face2
            tri_faces[it+2, :] = face3
            tri_faces[it+3, :] = face4

            face1.sort()
            face2.sort()
            face3.sort()
            face4.sort()
            tri_face_to_eids[tuple(face1)].append(eid)
            tri_face_to_eids[tuple(face2)].append(eid)
            tri_face_to_eids[tuple(face3)].append(eid)
            tri_face_to_eids[tuple(face4)].append(eid)
            it += 4
            eid += 1

        it_start[2] = it
        iq_start[2] = iq
        min_eids[eid] = ugrid.hexas
        ugrid.log.debug('HEXA it=%s iq=%s' % (it, iq))
        for element in ugrid.hexas-1:
            (n1, n2, n3, n4, n5, n6, n7, n8) = element

            face1 = [n1, n2, n3, n4]
            face2 = [n2, n6, n7, n3]
            face3 = [n6, n5, n8, n7]
            face4 = [n5, n1, n4, n8]
            face5 = [n4, n3, n7, n8]
            face6 = [n5, n6, n2, n1]

            quad_faces[iq, :] = face1
            quad_faces[iq+1, :] = face2
            quad_faces[iq+2, :] = face3
            quad_faces[iq+3, :] = face4
            quad_faces[iq+4, :] = face5
            quad_faces[iq+5, :] = face6

            face1.sort()
            face2.sort()
            face3.sort()
            face4.sort()
            face5.sort()
            face6.sort()

            quad_face_to_eids[tuple(face1)].append(eid)
            quad_face_to_eids[tuple(face2)].append(eid)
            quad_face_to_eids[tuple(face3)].append(eid)
            quad_face_to_eids[tuple(face4)].append(eid)
            quad_face_to_eids[tuple(face5)].append(eid)
            quad_face_to_eids[tuple(face6)].append(eid)
            iq += 6
            eid += 1

        it_start[3] = it
        iq_start[3] = iq
        min_eids[eid] = ugrid.penta5s
        ugrid.log.debug('PENTA5 it=%s iq=%s' % (it, iq))
        for element in ugrid.penta5s-1:
            (n1, n2, n3, n4, n5) = element

            face1 = [n2, n3, n5]
            face2 = [n1, n2, n5]
            face3 = [n4, n1, n5]
            face4 = [n5, n3, n4]
            face5 = [n4, n3, n2, n1]

            tri_faces[it, :] = face1
            tri_faces[it+1, :] = face2
            tri_faces[it+2, :] = face3
            tri_faces[it+3, :] = face4
            quad_faces[iq, :] = face5

            face1.sort()
            face2.sort()
            face3.sort()
            face4.sort()
            face5.sort()

            tri_face_to_eids[tuple(face1)].append(eid)
            tri_face_to_eids[tuple(face2)].append(eid)
            tri_face_to_eids[tuple(face3)].append(eid)
            tri_face_to_eids[tuple(face4)].append(eid)
            quad_face_to_eids[tuple(face5)].append(eid)

            it += 4
            iq += 1
            eid += 1

        it_start[4] = it
        iq_start[4] = iq
        min_eids[eid] = ugrid.penta6s
        ugrid.log.debug('PENTA6 it=%s iq=%s' % (it, iq))
        for element in ugrid.penta6s-1:
            (n1, n2, n3, n4, n5, n6) = element

            face1 = [n1, n2, n3]
            face2 = [n5, n4, n6]
            face3 = [n2, n5, n6, n3]
            face4 = [n4, n1, n3, n6]
            face5 = [n4, n5, n2, n1]

            tri_faces[it, :] = face1
            tri_faces[it+1, :] = face2
            quad_faces[iq, :] = face3
            quad_faces[iq+1, :] = face4
            quad_faces[iq+2, :] = face5

            face1.sort()
            face2.sort()
            face3.sort()
            face4.sort()
            face5.sort()

            tri_face_to_eids[tuple(face1)].append(eid)
            tri_face_to_eids[tuple(face2)].append(eid)
            quad_face_to_eids[tuple(face3)].append(eid)
            quad_face_to_eids[tuple(face4)].append(eid)
            quad_face_to_eids[tuple(face5)].append(eid)
            it += 2
            iq += 3
            eid += 1

        # find the unique faces
        tri_faces_sort = deepcopy(tri_faces)
        quad_faces_sort = deepcopy(quad_faces)
        #print('t0', tri_faces_sort[0, :])
        #print('t1', tri_faces_sort[1, :])

        ugrid.log.debug('nt=%s nq=%s' % (ntri_faces, nquad_faces))
        tri_faces_sort.sort(axis=1)
        #for i, tri in enumerate(tri_faces):
            #assert tri[2] > tri[0], 'i=%s tri=%s' % (i, tri)
        #print('*t0', tri_faces_sort[0, :])
        #print('*t1', tri_faces_sort[1, :])

        quad_faces_sort.sort(axis=1)
        #for i, quad in enumerate(quad_faces):
            #assert quad[3] > quad[0], 'i=%s quad=%s' % (i, quad)


        #iq_start_keys = iq_start.keys()
        #it_start_keys = it_start.keys()
        #iq_start_keys.sort()
        #it_start_keys.sort()

        unused_face_to_eid = []

        eid_keys = min_eids.keys()
        eid_keys.sort()

        type_mapper = {
            1 : 'tets',
            2 : 'hexas',
            3 : 'penta5s',
            4 : 'penta6s',
        }
        ugrid.log.info("eid_keys = %s" % eid_keys)
        for face, eids in tri_face_to_eids.items():
            if len(eids) == 1:
                #if it's a boundary face, wer're fine, otherwise, error...
                #print('*face=%s eids=%s' % (face, eids))
                #pid = lookup from quads/tris
                eid = eids[0]
                unused_owner = eid
                unused_neighbor = -1
                continue
                #raise RuntimeError()

            e1, e2 = eids
            i1 = searchsorted(eid_keys, e1)
            i2 = searchsorted(eid_keys, e2)

            if i1 == 1: # tet
                it1 = (e1-1) * 4
                it2 = (e1-1) * 4 + 4
                faces1_sort = tri_faces_sort[it1:it2, :]
                faces1_unsorted = tri_faces[it1:it2, :]

                #print("faces1 = \n", faces1_sort, '\n')

                # figure out irow; 3 for the test case
                face = array(face, dtype='int32')

                #print('face  = %s' % face)
                #print('face3 = %s' % faces1_sort[3, :])

                if allclose(face, faces1_sort[0, :]):
                    n1 = 0
                elif allclose(face, faces1_sort[1, :]):
                    n1 = 1
                elif allclose(face, faces1_sort[2, :]):
                    n1 = 2
                elif allclose(face, faces1_sort[3, :]):
                    n1 = 3
                else:
                    raise RuntimeError('cant find face=%s in faces for eid1=%s' % (face, e1))

                if allclose(face, faces1_unsorted[n1, :]):
                    unused_owner = e1
                    unused_neighbor = e2
                else:
                    unused_owner = e2
                    unused_neighbor = e1
                unused_face_new = faces1_unsorted[n1, :]

            elif i1 == 2:  # CHEXA
                iq1 = iq_start[2]
                unused_iq2 = iq1 + 6

            elif i1 == 3:  # CPENTA5
                #e1_new = e1 - eid_keys[2]
                iq1 = iq_start[3]
                unused_iq2 = iq1 + 1
                it1 = it_start[3]
                it2 = it1 + 4
            elif i1 == 4:  # CPENTA6
                iq1 = iq_start[4]
                unused_iq2 = iq1 + 3
                it1 = it_start[4]
                it2 = it1 + 2
            else:
                raise NotImplementedError('This is a %s and is not supported' % type_mapper[i1])

            # do we need to check this???
            if 0:
                if i2 == 1: # tet

                    it1 = it_start_keys[i1]
                    it2 = it1 + 4
                    unused_faces2 = tri_faces_sort[it1:it2, :]
                    #print('face=%s eids=%s' % (face, eids))
                    #print("faces2 = \n", unused_faces2)
                    # spits out 3
                else:
                    asdf
            #type1 = type_mapper[i1]
            #type2 = type_mapper[i2]
            #if type1:
        faces_file.write(')\n')
    return


def main():  # pragma: no cover
    """Tests UGrid"""
    ugrid_filename = 'bay_steve_recon1_fixed0.b8.ugrid'
    #bdf_filename = 'bay_steve_recon1_fixed0.b8.bdf'
    foam_filename = 'bay_steve_recon1_fixed0.b8.foam'
    tag_filename = 'bay_steve.tags'
    assert os.path.exists(tag_filename)
    ugrid_model = read_ugrid(ugrid_filename)
    #ugrid_model.write_bdf(bdf_filename)
    write_foam(ugrid_model, foam_filename, tag_filename)


if __name__ == '__main__':  # pragma: no cover
    main()
