"""
defines:
 - get_element_faces(model, element_ids=None)
 - get_solid_skin_faces(model)
 - write_skin_solid_faces(model, skin_filename,
                          write_solids=False, write_shells=True,
                          size=8, is_double=False, encoding=None)
"""
from __future__ import print_function
import sys
from codecs import open
from copy import deepcopy
from collections import defaultdict

from six import string_types

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf import read_bdf

def get_element_faces(model, element_ids=None):
    """
    Gets the elements and faces that are skinned from solid elements.
    This includes internal faces.

    Parameters
    ----------
    model : BDF()
        the BDF object
    element_ids : List[int] / None
        skin a subset of element faces
        default=None -> all elements

    Returns
    -------
    eid_faces : (int, List[(int, int, ...)])
       value1 : element id
       value2 : face

    """
    if element_ids is None:
        element_ids = model.element_ids

    eid_faces = []
    for eid in element_ids:
        elem = model.elements[eid]
        if elem.type in ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM']:
            faces = elem.faces
            #print(faces)
            for face_id, face in faces.items():
                if None in face:
                    msg = 'There is a None in the face.\n'
                    msg = 'face_id=%s face=%s\n%s' % (face_id, str(face), str(elem))
                    raise RuntimeError(msg)
                eid_faces.append((eid, face))
    return eid_faces


def get_solid_skin_faces(model):
    # type: (BDF) -> Any
    """
    Gets the elements and faces that are skinned from solid elements
    This doesn't include internal faces.

    Parameters
    ----------
    model : BDF()
        the BDF object

    Returns
    -------
    eid_set : Dict[tuple(int, int, ...)] = List[int]
       key : sorted face
       value : list of element ids with that face
    face_map : Dict[tuple(int, int, ...)] = List[int]
       key : sorted face
       value : unsorted face

    """
    eid_faces = get_element_faces(model)
    face_set = defaultdict(int)
    eid_set = defaultdict(list)
    face_map = {}
    for eid, face in eid_faces:
        #print(eid, face)
        raw_face = deepcopy(face)
        try:
            face.sort()
        except:
            print('face = %s' % str(face))
            raise
        tface = tuple(face)
        #print(tface)
        face_set[tface] += 1
        eid_set[tface].append(eid)
        face_map[tface] = raw_face

    #print('eid_set:')
    #for tface, eidset in eid_set.items():
        #print(tface, eidset)

    #print('face_set:')
    #for tface, faceset in face_set.items():
        #print(tface, faceset)

    #print('face_map:')
    #for tface, facemap in face_map.items():
        #print(tface, facemap)

    del_faces = []
    for face, face_count in face_set.items():
        if face_count == 2:
            del_faces.append(face)

    for face in del_faces:
        del face_set[face]
        del eid_set[face]
    return eid_set, face_map


def write_skin_solid_faces(model, skin_filename,
                           write_solids=False, write_shells=True,
                           size=8, is_double=False, encoding=None,
                           punch=False, log=None):
    """
    Writes the skinned elements

    Parameters
    ----------
    model : BDF() or str
        BDF : the BDF object
        str : bdf_filename and the read_bdf method is called
    skin_filename : str
        the file to write
    write_solids : bool; default=False
        write solid elements that have skinned faces
    write_shells : bool; default=False
        write shell elements
    size : int; default=8
        the field width
    is_double : bool; default=False
        double precision flag
    encoding : str; default=None -> system default
        the string encoding
    log : logger; default=None
        a python logging object
    punch : bool; default=False
        is this a punch file; should  be used by the read_bdf if model is a string
        unused

    """
    if isinstance(model, string_types):
        model = read_bdf(model, log=log)
    if(len(model.element_ids) == 0 or len(model.material_ids) == 0 or
       len(model.property_ids) == 0):
        msg = 'returning due to no elements/materials/properties\n'
        msg += '  nelements=%s nmaterials=%s nproperties=%s' % (
            len(model.element_ids), len(model.material_ids),
            len(model.property_ids)
        )
        model.log.warning(msg)
        return
    eid_set, face_map = get_solid_skin_faces(model)
    if len(eid_set) == 0:
        model.log.info('returning due to no elements in set')
        return

    eid_set_to_write = set()
    nid_set_to_write = set()
    mid_set_to_write = set()
    if write_solids:
        for face, eids in eid_set.items():
            eid_set_to_write.update(eids)
            for eid in eids:
                elem = model.elements[eid]
                pid = elem.Pid()
                prop = model.properties[pid] # PSOLID
                mid = prop.Mid()
                #print(prop)
                nid_set_to_write.update(elem.node_ids)
                mid_set_to_write.add(mid)
                #print('added_mid (a) =', mid)
    elif write_shells:
        for face, eids in eid_set.items():
            eid_set_to_write.update(eids)
            nid_set_to_write.update(face)
            for eid in eids:
                elem = model.elements[eid]
                pid = elem.Pid()
                prop = model.properties[pid] # PSOLID
                #print(prop)
                try:
                    #print(prop.mid)
                    mid = prop.Mid()
                    mid_set_to_write.add(mid)
                    #print('added eid=%s pid=%s mid=%s (b)' % (eid, pid, mid))
                except AttributeError:
                    continue
    else:
        raise RuntimeError('write_solids=False write_shells=False')

    eids_to_write = list(eid_set_to_write)
    nids_to_write = list(nid_set_to_write)
    mids_to_write = list(mid_set_to_write)

    #element_ids_to_delete = set(model.element_ids) - eids_to_write

    eid_shell = max(model.elements) + 1
    pid_shell = max(model.properties) + 1
    mid_shell = max(model.materials) + 1
    _write_skin_solid_faces(model, skin_filename, face_map,
                            nids_to_write, eids_to_write, mids_to_write, eid_set,
                            eid_shell, pid_shell, mid_shell,
                            write_solids=write_solids, write_shells=write_shells,
                            size=size, is_double=is_double, encoding=encoding)

def _write_skin_solid_faces(model, skin_filename, face_map,
                            nids_to_write, eids_to_write, mids_to_write, eid_set,
                            eid_shell, pid_shell, mid_shell,
                            write_solids=False, write_shells=True,
                            size=8, is_double=False, encoding=None):
    """
    helper method for ``write_skin_solid_faces``

    Parameters
    ----------
    model : BDF()
        the BDF object
    skin_filename : str
        the file to write
    face_map : dict[sorted_face] : face
        sorted_face : List[int, int, int] / List[int, int, int, int]
        face : List[int, int, int] / List[int, int, int, int]
    nids_to_write : List[int, int, ...]
        list of node ids to write
    eids_to_write : List[int, int, ...]
        list of element ids to write
    mids_to_write : List[int, int, ...]
        list of material ids to write
    eid_set : Set[int]
        is the type right???
    eid_shell : int
        the next id to use for the shell id
    pid_shell : int
        the next id to use for the shell property
    mid_shell : int
        the next id to use for the shell material
    write_solids : bool; default=False
        write solid elements that have skinned faces
    write_shells : bool; default=True
        write shell elements
    size : int; default=/8
        the field width
    is_double : bool; default=False
        double precision flag
    encoding : str; default=None -> system default
        the string encoding

    """
    #encoding = model.get_encoding(encoding)
    with open(skin_filename, 'w') as bdf_file:
        bdf_file.write('$ pyNastran: punch=True\n')
        for nid in sorted(nids_to_write):
            if nid is None:
                continue
            node = model.nodes[nid]
            bdf_file.write(node.write_card(size=size, is_double=is_double))

        for cid, coord in model.coords.items():
            if cid == 0:
                continue
            bdf_file.write(coord.write_card(size=size, is_double=is_double))

        if write_solids:
            for eid in sorted(eids_to_write):
                elem = model.elements[eid]
                bdf_file.write(elem.write_card(size=size))
            for pid, prop in model.properties.items():
                bdf_file.write(prop.write_card(size=size, is_double=is_double))
            for mid in sorted(mids_to_write):
                material = model.materials[mid]
                bdf_file.write(material.write_card(size=size, is_double=is_double))
            del eid, pid, mid

        if write_shells:
            eid_shell = _write_shells(bdf_file, model, eid_set, face_map, eid_shell, pid_shell, mid_shell, mids_to_write)
        bdf_file.write('ENDDATA\n')
    #if 0:
        #model = model.__class__.__init__()
        #model.read_bdf(skin_filename)


def _write_shells(bdf_file, model, eid_set, face_map, eid_shell, pid_shell, mid_shell, mids_to_write):
    """helper method for ``_write_skin_solid_faces``"""
    mids_to_write.sort()
    for imid, mid in enumerate(mids_to_write):
        card = ['PSHELL', pid_shell + imid, mid_shell + imid, 0.1]
        bdf_file.write(print_card_8(card))

        card = ['MAT1', mid_shell + imid, 3.e7, None, 0.3]
        #bdf_file.write(model.materials[mid].comment)
        bdf_file.write(print_card_8(card))

    for face, eids in eid_set.items():
        face_raw = face_map[face]
        nface = len(face)
        #assert len(eids) == 1, eids
        #for eid in sorted(eids):
            #elem = model.elements[eid]
            #print(elem)
            #break

        #elem = next(itervalues(model.elements)) # old
        elem = model.elements[eids[0]]
        #pid = next(iterkeys(model.properties))
        pid = elem.Pid()
        prop = model.properties[pid]
        try:
            mid = prop.Mid()
        except AttributeError:
            continue
        #print('mids_to_write = %s' % mids_to_write)
        #print('mids = ', model.materials.keys())
        imid = mids_to_write.index(mid)

        if nface == 3:
            card = ['CTRIA3', eid_shell, pid_shell + imid] + list(face_raw)
        elif nface == 4:
            card = ['CQUAD4', eid_shell, pid_shell + imid] + list(face_raw)
        elif nface == 4:
            card = ['CQUAD4', eid_shell, pid_shell + imid] + list(face_raw)
        elif nface == 6:
            card = ['CTRIA6', eid_shell, pid_shell + imid] + list(face_raw)
        elif nface == 8:
            card = ['CQUAD8', eid_shell, pid_shell + imid] + list(face_raw)
        else:
            raise NotImplementedError('face=%s len(face)=%s' % (face, nface))
        bdf_file.write(print_card_8(card))
        eid_shell += 1

        #elem = model.elements[eid]
        #bdf_file.write(elem.write_card(size=size))
    #for pid, prop in model.properties.items():
        #bdf_file.write(prop.write_card(size=size, is_double=is_double))
    return eid_shell


def main():  # pragma: no cover
    """docopt interface"""
    encoding = sys.getdefaultencoding()
    import pyNastran
    from docopt import docopt
    msg = (
        'Usage:\n'
        '  free_faces [-d] [-f] [--encoding ENCODE] BDF_FILENAME SKIN_FILENAME\n'
        '  free_faces [-l] [-f] [--encoding ENCODE] BDF_FILENAME SKIN_FILENAME\n'

        '  free_faces -h | --help\n'
        '  free_faces -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  BDF_FILENAME   path to input BDF/DAT/NAS file\n"
        "  SKIN_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        '  -l, --large        writes the BDF in large field, single precision format (default=False)\n'
        '  -d, --double       writes the BDF in large field, double precision format (default=False)\n'
        '  --encoding ENCODE  the encoding method (default=None -> %r)\n'
        '\n'
        'Developer:\n'
        '  -f, --profile    Profiles the code (default=False)\n'
        '\n'
        "Info:\n"
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n" % encoding
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    #data['--xref'] = not data['--xref']
    if not data['--encoding']:
        data['--encoding'] = None

    for key, value in sorted(data.items()):
        print("%-12s = %r" % (key.strip('--'), value))

    import time
    time0 = time.time()

    #is_double = False
    if data['--double']:
        size = 16
        #is_double = True
    elif data['--large']:
        size = 16
    else:
        size = 8
    bdf_filename = data['BDF_FILENAME']
    skin_filename = data['SKIN_FILENAME']


    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    tol = 1e-005
    bdf_filename_merged = 'merged.bdf'
    bdf_equivalence_nodes(bdf_filename, bdf_filename_merged, tol,
                          renumber_nodes=False, neq_max=10, xref=True,
                          node_set=None,
                          size=8, is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False, debug=True)
    print('done with equivalencing')
    write_skin_solid_faces(
        bdf_filename_merged, skin_filename,
        write_solids=False, write_shells=True,
        size=size, is_double=False, encoding=None,
    )
    print("total time:  %.2f sec" % (time.time() - time0))

if __name__ == '__main__':  # pragma: no cover
    main()
