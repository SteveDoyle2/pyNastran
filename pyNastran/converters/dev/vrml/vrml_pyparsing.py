import time
from pyparsing import (
    Group, Word,
    OneOrMore, Forward, Literal,
    delimitedList,
)
from pyNastran.converters.dev.vrml.parsing_help import (
    pfloat, cvt_int, pint, pword, name_str, dict_open, dict_close,
    xyz, pword_num_underscore,
    directional_light, world_info, background, navigation_info, shape, geometry)

def remove_comments(lines):
    lines2 = []
    for line in lines:
        line = line.split('#')[0].rstrip()
        if not line:
            continue
        lines2.append(line)
    return lines2


def get_vrml_format_old():
    name_float = pword + pfloat
    name_float3 = pword + Group(pfloat * 3)

    pminus1 = Word('-1')
    list_open = Literal('[').suppress()
    list_close = Literal(']').suppress()
    int_list = list_open + delimitedList(pint.setParseAction(cvt_int) | pminus1.setParseAction(cvt_int)) + list_close
    float_list = list_open + delimitedList(pfloat) + list_close

    name_int_list = pword + Group(int_list)
    name_float_list = pword + Group(float_list)
    data_value = name_int_list | name_float_list | name_float + name_float3

    name_class_dict = pword + pword + dict_open + data_value + dict_close
    name_dict = Forward()
    name_dict <<= pword + dict_open + data_value + dict_close
    #list_item = real | integer # | quotedString
    #pfloat = Combine(Optional(oneOf("+-")) + Word(nums) + '.' + Word(nums))
    #triple = nums + White + nums + White + nums
    #pfloat = pfloat_neg
    #names_dict = Optional(pword) + name_dict

    # name_dict.parse_string(simple_shape)
    # name_class_dict.parse_string(simple_shape)

    msg = """
      coordIndex [
       0, 1, 2, -1, 0, 3, 4, -1, 4, 1, 0, -1,
       0, 5, 6, -1, 0, 7, 8, -1, 0, 9, 10, -1,
       0, 11, 3, -1, 0, 12, 13, -1, 0, 14, 5, -1,
       8, 15, 0, -1, 0, 16, 9, -1, 0, 17, 18, -1,
       0, 19, 20, -1, 13, 11, 0, -1, 0, 21, 22, -1,
       15, 14, 0, -1, 18, 16, 0, -1, 2, 19, 0, -1,
       0, 23, 24, -1, 0, 25, 26, -1, 0, 27, 28, -1,
       22, 29, 0, -1, 0, 30, 31, -1, 6, 32, 0, -1,
       20, 17, 0, -1, 24, 33, 0, -1, 0, 34, 35, -1,
       0, 36, 25, -1, 29, 27, 0, -1, 0, 37, 38, -1,
       32, 30, 0, -1, 0, 39, 12, -1, 35, 23, 0, -1,
       28, 36, 0, -1, 31, 37, 0, -1, 33, 39, 0, -1,
       38, 21, 0, -1, 26, 34, 0, -1,
       ]
    """
    # print(name_int_list.parse_string(msg).asList())
    # print(name_int_list.parse_string(msg).asDict())


def get_vrml_format():
    #sample = """
    ##VRML V2.0 utf8

    #DirectionalLight {
     #direction 0.577 -0.577 -0.577
     #color    1.000 1.000 1.000
     #intensity 0.450
     #ambientIntensity 1.0
    #}

    #DirectionalLight {
     #direction -0.577 -0.577 -0.577
     #color    1.000 1.000 1.000
     #intensity 0.680
     #ambientIntensity 1.0
    #}
    #"""

    # 1.0
    # -1.0
    # -1.
    #pminus1 = Word('-1')
    list_open = Literal('[').suppress()
    list_close = Literal(']').suppress()

    name_float = pword + pfloat
    name_float3 = pword + Group(pfloat * 3)
    name_float4 = pword + Group(pfloat * 4)

    #pfloat_pos = Combine(Optional('+') + pint + '.' + Optional(pint))
    #pfloat_neg = Combine(Optional('-') + pint + '.' + Optional(pint))
    #pfloat_pos2 = Combine(Optional('+') + Optional(pint) + '.' + pint)
    #pfloat_neg2 = Combine(Optional('-') + Optional(pint) + '.' + pint)
    #pfloat = pfloat_neg | pfloat_neg2 | pfloat_pos | pfloat_pos2

    #name_name_float_list = pword + pword + float_list

    #pword_triple = pword + pfloat * 3
    #pword_float = pword + pfloat


    #print(pword_float.parse_string('1.'))

    # color    1.000 1.000 1.000
    # direction -0.577 -0.577 -0.577
    # out = pword_triple.parse_string('color    1.000 1.000 1.000', parseAll=False)  # works
    unused_out = name_float3.parse_string('direction -0.577 -0.577 -0.577', parseAll=False)  # works

    # ambientIntensity 1.0
    unused_out = name_float.parse_string('ambientIntensity 1.0', parseAll=False)
    #print(out, dir(type(out)))  # ParsingResults
    #print(out.asDict())
    #print(out.asList())

    #----------------------------------------
    #simple_shape = """
    #Shape {
        #appearance Appearance{
             #texture DEF PICBAND PixelTexture {
                 #image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                              #0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                              #0x7700FF77 0x444444FF
             #}
        #}
        #geometry Sphere{}
    #}
    #"""

    #hexi = Word("047xF")
    #simple_shape = """
    #appearance Appearance{
         #texture DEF PICBAND PixelTexture {
             #image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                          #0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                          #0x7700FF77 0x444444FF
         #}
    #}
    #"""

    xyz_vector = list_open + OneOrMore(Group(xyz)) + list_close
    name_xyz_vector = pword + xyz_vector

    data_value = name_float3 | name_str | name_float | name_xyz_vector
    data_values = OneOrMore(data_value)
    #name_dict = pword + Group(dict_open + data_value + dict_close)
    name_dict = pword + Group(dict_open.suppress() + data_values + dict_close.suppress())
    name_float3.parse_string('skyColor 0.1 0.3 1')

    name_dict.parse_string("""
    Background {
        skyColor 0.1 0.3 1
    }
    """)

    name_str.parse_string('''
        title "Texture-mapped pyramid"
    ''')

    name_dict.parse_string("""
    WorldInfo {
        title "Texture-mapped pyramid"
    }
    """)

    name_dict.parse_string("""
    WorldInfo {
        title "Texture-mapped pyramid"
        info  "Gravity: on"
    }
    """)

    name_dict.parse_string("""
    DirectionalLight {
     direction 0.577 -0.577 -0.577
     color    1.000 1.000 1.000
     intensity 0.450
     ambientIntensity 1.0
    }
    """)
    xyz.parse_string("""
      0 0 -1,
    """)

    xyz_vector.parse_string("""
    [
      0 0 -1,  0 0 -1,  0 0 -1,
      0 0 -1,  0 0 -1,  0 0 -1,
    ]
    """)

    name_xyz_vector.parse_string("""
    vector [
      0 0 -1,  0 0 -1,  0 0 -1,
      0 0 -1,  0 0 -1,  0 0 -1,
    ]
    """)

    #print(names_dict.parse_string("""
    #normal Normal {
     #vector [
      #0 0 -1,  0 0 -1,  0 0 -1,
      #0 0 -1,  0 0 -1,  0 0 -1,
     #]
    #}
    #"""))

    shape.parse_string("""
    Shape {
        appearance Appearance{
             texture DEF PICBAND PixelTexture {
                 image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                              0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                              0x7700FF77 0x444444FF
             }
        }
        geometry Sphere{}
    }
    """)
    #txt = read_vrml('pyramid_sphere.wrl')
    #----------------------------------------
    txt = """
    WorldInfo {
            title "Texture-mapped pyramid"
            info  "Gravity: on"
    }
    Background {
            skyColor 0.1 0.3 1
    }
    NavigationInfo {
        type "EXAMINE"
            headlight TRUE
    }
    Shape {
        appearance Appearance{
             texture DEF PICBAND PixelTexture {
                 image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                              0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                              0x7700FF77 0x444444FF
             }
        }
        geometry Sphere{}
    }
    """
    child = OneOrMore(shape)
    children = Literal('children') + list_open + Group(child) + list_close
    translation = Literal('translation') + name_float3
    rotation = Literal('rotation') + name_float4
    transform_values = children | translation | rotation
    transform1 = Literal('Transform') + dict_open + transform_values + dict_close
    transform2 = Literal('DEF') + pword_num_underscore + transform1
    transform = transform1 | transform2
    geometry_str = """
    geometry IndexedFaceSet {
        coord Coordinate{
            point[
               -2 -2 0,
                2 -2 0,
                2  2 0,
                -2  2 0,
                0  0 5,
            ]
        }
        coordIndex [
            0, 1, 4, -1,
            1, 2, 4, -1,
            2, 3, 4, -1,
            3, 0, 4, -1,
            3, 2, 1, 0, -1,
        ]
        texCoord TextureCoordinate {
            point [
                0 0,
                0 .3,
                0 .5,
                0 .7,
                0 1,
            ]
        }
    }
    """
    #print('geometry...')
    geometry.parse_string(geometry_str)

    shape_str = """
    Shape{
        appearance Appearance{
            texture DEF PICBAND ImageTexture {
                url "http://www.rt.cs.boeing.com/people/davidk/wrl/geo/colors.jpg"
                repeatS FALSE
                repeatT FALSE
            }
        }
        geometry IndexedFaceSet {
            coord Coordinate{
                point[
                   -2 -2 0,
                    2 -2 0,
                    2  2 0,
                    -2  2 0,
                    0  0 5,
                ]
            }
            coordIndex [
                0, 1, 4, -1,
                1, 2, 4, -1,
                2, 3, 4, -1,
                3, 0, 4, -1,
                3, 2, 1, 0, -1,
            ]
            texCoord TextureCoordinate {
                point [
                    0 0,
                    0 .3,
                    0 .5,
                    0 .7,
                    0 1,
                ]
            }
        }
    }
    """
    transform_str = """
    Transform{
        children[
            Shape{
                appearance Appearance{
                    texture DEF PICBAND ImageTexture {
                        url "http://www.rt.cs.boeing.com/people/davidk/wrl/geo/colors.jpg"
                        repeatS FALSE
                        repeatT FALSE
                    }
                }
                geometry IndexedFaceSet {
                   coord Coordinate{
                       point[
                           -2 -2 0,
                            2 -2 0,
                            2  2 0,
                           -2  2 0,
                            0  0 5,
                       ]
                   }
                   coordIndex [
                       0, 1, 4, -1,
                       1, 2, 4, -1,
                       2, 3, 4, -1,
                       3, 0, 4, -1,
                       3, 2, 1, 0, -1,
                   ]
                   texCoord TextureCoordinate {
                       point [
                            0 0,
                            0 .3,
                            0 .5,
                            0 .7,
                            0 1,
                       ]
                    }
                }
            }
        ]
    }
    """
    #print(txt)

    vrml_format = OneOrMore(
        Group(directional_light) | world_info | background | navigation_info | shape | transform
    )
    vrml_format.parse_string(txt)
    #print('shape...')
    shape.parse_string(shape_str)
    #print('transform...')
    transform.parse_string(transform_str)


    if 0:
        print('ready for gbu')
        #txt = read_vrml('gbu.wrl')

        # t_no_float_regex = 63 sec
        # t_float_regex = 31 sec
        t0 = time.time()
        vrml_format.parse_string(txt, parseAll=True)
        print(time.time() - t0)

        #for datai in data:
            #print('  ', datai)
        #print('done!!!!!!!!')

        #import json
        #with open('gbu.json', 'w') as fp:
            #json.dump(data, fp, indent=4)

    return vrml_format
