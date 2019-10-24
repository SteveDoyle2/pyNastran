import numpy as np
from pyparsing import ParseResults

def todicti(data, log):
    i = 0
    dicti = {}
    log.debug('------')
    while i < len(data):
        #print(dicti)
        key = data[i]
        #print('%r' % modeli)
        i += 1
        value = data[i]
        if isinstance(value, str):
            dicti[key] = value
        elif isinstance(value, ParseResults):
            value2 = value.asList()
            log.debug(f'  key={key} value={value2}')
            dicti[key] = value2
        else:
            print(f'  key={key} value={value}')
            print(dicti)
            raise NotImplementedError(f'value={value!r} type={type(value)}')
        i += 1
    #print(data)
    #print(dicti)
    return dicti

def read_transforms(pmodel):
    #print('**************************')
    transforms = []
    i = 0
    i = 0
    while i < len(pmodel):
        modeli = pmodel[i]
        if modeli == 'Shape':
            i += 1
            data = pmodel[i]
            shape = toshape(data)
            #print('shape =', shape)
            i += 1
            transforms.append(shape)
        else:
            raise NotImplementedError(f'modeli={modeli!r} type={type(modeli)}')
    return transforms

def tolist(data):
    data2 = []
    for datai in data:
        if isinstance(datai, ParseResults):
            datai2 = tolist(datai)
            data2.append(datai2)
        elif isinstance(datai, (int, float, str)):
            data2.append(datai)
        else:
            raise NotImplementedError(f'datai={datai!r} type={type(datai)}')
    return data2

def to_quads_tris(coord_indexs):
    coord_indexs = np.array(coord_indexs, dtype='int32')
    #print(coord_indexs)
    iminus1 = np.where(coord_indexs == -1)[0]
    #print(iminus1)
    i0 = 0
    quads = []
    tris = []
    for i1 in iminus1:
        datai = coord_indexs[i0:i1]
        #print(datai)
        ndata = len(datai)
        if ndata == 3:
            tris.append(datai)
        elif ndata == 4:
            quads.append(datai)
        else:
            raise NotImplementedError(datai)
        i0 = i1 + 1

    if quads:
        quads = np.array(quads)
    if tris:
        tris = np.array(tris)
    return coord_indexs, quads, tris

def read_indexed_face_set(data):
    indexed_face_set = {}
    #print(data)
    i = 0
    while i < len(data):
        datai = data[i]
        if datai == 'coord':
            i += 1
            datai = data[i]
            assert datai == 'Coordinate', datai
            i += 1
            datai = data[i]

            assert len(datai) == 2, datai
            point = datai[1]
            points = np.array(tolist(point), dtype='float64')
            indexed_face_set['coord'] = {
                'point' : points,
            }

            i += 1
        elif datai == 'coordIndex':
            i += 1
            coord_index = data[i]
            coord_indexs = tolist(coord_index)
            coord_indexs, quads, tris = to_quads_tris(coord_indexs)
            #indexed_face_set['coord_index'] = coord_index
            if len(quads):
                indexed_face_set['quads'] = quads
            if len(tris):
                indexed_face_set['tris'] = tris
            i += 1

        elif datai == 'normal':
            i += 1
            datai = data[i]
            assert datai == 'Normal', datai
            i += 1
            vector_normals = data[i]
            #print(vector_normals)
            assert len(vector_normals) == 2, vector_normals
            normals = np.array(tolist(vector_normals[1]), dtype='float64')
            indexed_face_set['normals'] = normals
            i += 1
        elif datai == 'texCoord':
            i += 1
            datai = data[i]
            assert datai == 'TextureCoordinate', datai

            i += 1
            datai = data[i]
            assert datai == 'point', datai

            i += 1
            point = data[i]
            #print('point!', point)

            points = np.array(tolist(point), dtype='float64')
            indexed_face_set['tex_coord'] = {
                'point' : points,
            }
            i += 1
        elif datai == 'creaseAngle':
            i += 1
            datai = data[i]
            indexed_face_set['crease_angle'] = datai
            i += 1
        else:
            raise NotImplementedError(f'datai={datai!r} type={type(datai)}')
    #print(type(indexed_face_set))
    return indexed_face_set

def toshape(data):
    shape = {}
    i = 0
    while i < len(data):
        datai = data[i]
        if datai == 'appearance':
            #key = datai
            i += 1
            datai = data[i]
            assert data[i] == 'Appearance', data[i]
            i += 1
            #datai = data[i]
            #i += 1
            #datai = data[i]
            #print(datai)
            #aa
            shape['appearance'] = None
            datai = data[i]

            if isinstance(datai, ParseResults):
                datai = datai.asList()
                if datai[0] == 'texture':
                    i += 1
            datai = data[i]

            #return shape
            if isinstance(datai, ParseResults):
                datai = datai.asList()

                #['material', 'Material',
                #   ['ambientIntensity', 0.21,
                #    'diffuseColor', [0.49, 0.49, 0.49],
                #    'specularColor', [0.5, 0.5, 0.5],
                #    'transparency', 0.0,
                #    'shininess', 0.6]]
                i += 1
                #j = 0
                #while j < len(datai):
            datai = data[i]
            if datai == 'geometry':
                continue

            raise NotImplementedError(f'datai={datai!r} type={type(datai)}')

        elif datai == 'geometry':
            i += 1
            datai = data[i]
            #print(datai)
            if datai == 'Sphere':
                shape['geometry'] = 'sphere'
            elif datai == 'IndexedFaceSet':
                i += 1
                datai = data[i]
                indexed_face_set = read_indexed_face_set(datai)
                shape['indexed_face_set'] = indexed_face_set
            else:
                raise NotImplementedError(f'datai={datai!r} type={type(datai)}')
            i += 1
        else:
            raise NotImplementedError(f'datai={datai!r} type={type(datai)}')
    return shape


def todict(pmodel, log):
    out_model = {
        'transforms' : [],
    }
    i = 0
    while i < len(pmodel):
        modeli = pmodel[i]
        if modeli == 'WorldInfo':
            i += 1
            data = pmodel[i]
            out_model[modeli] = todicti(data, log)
            i += 1
        elif modeli == 'Background':
            i += 1
            data = pmodel[i]
            out_model[modeli] = todicti(data, log)
            i += 1
        elif modeli == 'NavigationInfo':
            i += 1
            data = pmodel[i]
            out_model[modeli] = todicti(data, log)
            i += 1
        elif modeli == 'Shape':
            i += 1
            data = pmodel[i]
            out_model[modeli] = toshape(data)
            i += 1
        elif modeli == 'Transform':
            i += 1
            data = pmodel[i] # children
            assert data == 'children', data
            #out_model[modeli] = totransform(data)
            i += 1
            data = pmodel[i]
            #print(data)
            transforms = read_transforms(data)
            #assert 'transforms' not in out_model
            out_model['transforms'].append(transforms)
            i += 1
        elif isinstance(modeli, ParseResults):
            modeli = modeli.asList()
            #print(modeli)

            #['DirectionalLight',
            #   ['direction', [0.577, -0.577, -0.577],
            #    'color', [1.0, 1.0, 1.0],
            #    'intensity', 0.45,
            #    'ambientIntensity', 1.0]]
            j = 0
            while j < len(modeli):
                datai = modeli[j]
                if datai == 'DirectionalLight':
                    j += 2
                else:
                    raise NotImplementedError(f'datai={datai!r} type={type(datai)}')
            i += 1
        elif modeli == 'DEF': # DEF Body__115117 Transform
            i += 1
            modeli = pmodel[i]
            log.debug('%r' % modeli)
            i += 1
            #modeli = pmodel[i]
            #print('%r' % modeli)
            #i += 1
        else:
            raise NotImplementedError(f'modeli={modeli!r} type={type(modeli)}')

    #print('----------------------------')
    return out_model
