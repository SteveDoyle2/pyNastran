"""
defines:
 - get_solid_skin_faces(model)

"""
from collections import defaultdict
from copy import deepcopy

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


def write_skin_solid_faces(model, skin_filename,
                           write_solids=False, write_shells=True,
                           size=8, is_double=False, encoding=None):
    """
    Writes the skinned elements

    Parameters
    ----------
    model : BDF()
        the BDF object
    skin_filename : str
        the file to write
    write_solids : bool; default=False
        write solid elements that have skinned faces
    write_shells : bool; default=False
        write newly created shell elements
        if there are shells in the model, doesn't write these
    size : int; default=8
        the field width
    is_double : bool; default=False
        double precision flag
    encoding : str; default=None -> system default
        the string encoding

    """
    if(len(model.element_ids) == 0 or len(model.material_ids) == 0 or
       len(model.property_ids) == 0):
        return
    eid_set, face_map = get_solid_skin_faces(model)
    if len(eid_set) == 0:
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
                if prop.type in ['PSOLID', 'PLSOLID']:
                    mid = prop.Mid()
                elif prop.type in ['PCOMPS', 'PCOMP', 'PCOMPG']:
                    mid = prop.mids[0]
                else:
                    raise NotImplementedError(prop)
                #except TypeError:
                    #model.log.warning('TypeError: skipping:%s' % prop)
                    #raise
                #except AttributeError:
                    #model.log.warning('skipping:%s' % prop)
                    #continue
                mid_set_to_write.add(mid)
                #print('added eid=%s pid=%s mid=%s (b)' % (eid, pid, mid))
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


def get_solid_skin_faces(model):
    """
    Gets the elements and faces that are skinned from solid elements.
    This doesn't include internal faces or existing shells.

    Parameters
    ----------
    model : BDF()
        the BDF object

    Returns
    -------
    eid_set : Dict[sorted_face] = eids
       sorted_face : tuple(int, int, ...)
           the face nids in sorted order
       eids : List[int]
           list of element ids with that face
    face_map : Dict[sorted_face] = face
       sorted_face : tuple(int, int, ...)
           the face nids in sorted order
       face : List(int, int, ...)
           the face nids

    """
    eid_faces = model.get_element_faces()
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
    eid_set : dict[face] : eids
        ???
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
    size : int; default=8
        the field width
    is_double : bool; default=False
        double precision flag
    encoding : str; default=None -> system default
        the string encoding

    """
    encoding = model.get_encoding(encoding)
    with open(skin_filename, 'w', encoding=encoding) as bdf_file:
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
            mids_to_write.sort()
            for imid, mid in enumerate(mids_to_write):
                card = ['PSHELL', pid_shell + imid, mid_shell + imid, 0.1]
                try:
                    msg = print_card_8(card)
                except RuntimeError:
                    msg = print_card_16(card)
                bdf_file.write(msg)

                card = ['MAT1', mid_shell + imid, 3.e7, None, 0.3]
                #bdf_file.write(model.materials[mid].comment)
                try:
                    msg = print_card_8(card)
                except RuntimeError:
                    msg = print_card_16(card)
                bdf_file.write(msg)

            for face, eids in eid_set.items():
                face_raw = face_map[face]
                nface = len(face)
                #print("eids =", eids)
                #assert len(eids) == 1, eids
                #for eid in sorted(eids):
                    #elem = model.elements[eid]
                    #print(elem)
                    #break

                assert len(eids) == 1, eids
                elem = model.elements[eids[0]]
                #pid = next(model.properties.keys())
                pid = elem.Pid()
                prop = model.properties[pid]
                if prop.type in ['PSOLID']: # 'PSHELL',
                    mid = prop.Mid()
                #elif prop.type in ['PCOMP', 'PCOMPG']:
                    #mid = prop.mids[0]
                else:
                    raise NotImplementedError(prop)

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
                try:
                    msg = print_card_8(card)
                except RuntimeError:
                    msg = print_card_16(card)
                bdf_file.write(msg)
                eid_shell += 1

                #elem = model.elements[eid]
                #bdf_file.write(elem.write_card(size=size))
            #for pid, prop in model.properties.items():
                #bdf_file.write(prop.write_card(size=size, is_double=is_double))
        bdf_file.write('ENDDATA\n')
    #if 0:
        #model = model.__class__.__init__()
        #model.read_bdf(skin_filename)
