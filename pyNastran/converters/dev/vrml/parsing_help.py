import pyparsing as pp
#from pyparsing import (
    #Suppress, Group, Optional, Word, ZeroOrMore, White, Combine,
    #Dict, Literal, OneOrMore, Regex,
    #alphas, alphanums, nums, oneOf, delimitedList, quotedString
#)

pword = pp.Word(pp.alphas).setName('word')
pword_underscore = pp.Word(pp.alphas + '_').setName('word_underscore')
pword_num_underscore = pp.Word(pp.alphas + pp.nums + '_').setName('word_num_underscore')
pint = pp.Word(pp.nums).setName('integer')
pint_sign = pp.Combine(pp.Optional(pp.oneOf("+ -")) + pp.Word(pp.nums)).setName('signed_integer')
pminus1 = pp.Word('-1')

#pint = Regex('/^[-+]?\d+$/')  # all integers
boolean = (pp.Literal('TRUE') | pp.Literal('FALSE')).setName('boolean')

# fourth pass, add parsing of dicts
cvt_int = lambda s, l, toks: int(toks[0])
cvt_float = lambda s, l, toks: float(toks[0])
#cvtDict = lambda s, l, toks: dict(toks[0])

# super fast at 31 sec
#float_regex = r"\d+\.\d+([Ee][+-]?\d+)?"
#float_regex = '[-+]?[0-9]*\.?[0-9]*'
float_regex = '[+-]?([0-9]*[.])?[0-9]+'

#  31.8 sec -> 28.2 sec if we drop casting...
pfloat = pp.Regex(float_regex).setName('real').setParseAction(cvt_float)
#pfloat_lazy = (pfloat1 | pfloat2 | pint_sign).setName('real').setParseAction(cvt_float)

pfloat.parseString('1.0')
pfloat.parseString('+1.0')
pfloat.parseString('-1.0')
pfloat.parseString('1.')
pfloat.parseString('+1.')
pfloat.parseString('-1.')
pfloat.parseString('-1')
pfloat.parseString('0')
pfloat.parseString('3')
#------------------------------------------------------

name_str = pword + pp.quotedString


comma = pp.Word(',').setName('comma')
xyz = pp.Group(pfloat * 3 + pp.Optional(comma.suppress())).setName('xyz')
xy = pp.Group(pfloat * 2 + pp.Optional(comma.suppress())).setName('xy')

# 0xFFFFFF77
hexa = pp.Word('0123456789ABCDEFx', min=10, max=10).setName('hex')

hexa.parseString("0xFFFFFF77")
hexa.parseString("0xFF0000FF")
hexa.parseString("0xFFCC0077")
hexa.parseString("0xFFFF00FF")
hexa.parseString("0x77FF00FF")
hexa.parseString("0x00FF00FF")
hexa.parseString("0x00FFFFFF")
hexa.parseString("0x7700FF77")
hexa.parseString("0x444444FF")

#0xFF0000FF 0xFFCC0077 0xFFFF00FF
#0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
#0x7700FF77 0x444444FF


#-----------------------------------------
list_open = pp.Literal('[').suppress()
list_close = pp.Literal(']').suppress()
dict_open = pp.Literal('{').suppress()
dict_close = pp.Literal('}').suppress()

color_datai = pp.Group(pfloat * 3)
#-----------------------------------------------
ambient_intensity = pp.Literal('ambientIntensity') + pfloat
diffuse_color = pp.Literal('diffuseColor') + color_datai
specular_color = pp.Literal('specularColor') + color_datai
transparency = pp.Literal('transparency') + pfloat
shininess = pp.Literal('shininess') + pfloat
material_values = pp.Group(pp.OneOrMore(
    ambient_intensity | diffuse_color | specular_color |
    transparency | shininess))
material = (
    pword + pp.Literal('Material') +
    dict_open +
    material_values +
    dict_close)
material_values.parseString("""
    ambientIntensity 0.210
    diffuseColor 0.596 0.667 0.686
    specularColor 0.500 0.500 0.500
    transparency 0.000
    shininess 0.600
""")

material.parseString("""
material Material {
    ambientIntensity 0.210
    diffuseColor 0.596 0.667 0.686
    specularColor 0.500 0.500 0.500
    transparency 0.000
    shininess 0.600
}
""")

# --------------------------------------
direction = pp.Literal('direction') + xyz
color = pp.Literal('color') + color_datai
#specular_color = pp.Literal('specularColor') + color_datai
intensity = pp.Literal('intensity') + pfloat
#shininess = pp.Literal('shininess') + pfloat
directional_light_values = pp.Group(pp.OneOrMore(
    direction | color | intensity | ambient_intensity))
directional_light = pp.Literal('DirectionalLight') + dict_open + directional_light_values + dict_close
directional_light.parseString("""
DirectionalLight {
 direction 0.577 -0.577 -0.577
 color    1.000 1.000 1.000
 intensity 0.450
 ambientIntensity 1.0
}
""")
# --------------------------------------
title = pp.Literal('title') + pp.quotedString
info = pp.Literal('info') + pp.quotedString
world_info_values = pp.OneOrMore(title | info)
world_info = (
    pp.Literal('WorldInfo') +
    dict_open +
    pp.Group(world_info_values) +
    dict_close)
world_info.parseString("""
WorldInfo {
    title "Texture-mapped pyramid"
    info  "Gravity: on"
}
""")
# --------------------------------------
sky_color = (pp.Literal('skyColor') + color_datai).setName('sky_color')
background_values = pp.OneOrMore(sky_color)
background = (
    pp.Literal('Background') +
    dict_open +
    pp.Group(background_values) +
    dict_close).setName('background')
background.parseString("""
Background {
    skyColor 0.1 0.3 1
}
""")
# --------------------------------------
typei = (pp.Literal('type') + pp.quotedString).setName('type')
headlight = (pp.Literal('headlight') + boolean).setName('headlight')
navigation_info_values = pp.OneOrMore(typei | headlight)
navigation_info = (
    pp.Literal('NavigationInfo') +
    dict_open +
    pp.Group(navigation_info_values) +
    dict_close).setName('navigation_info')
navigation_info.parseString("""
NavigationInfo {
 type "EXAMINE"
    headlight TRUE
}
""")
#-----------------------------------------------------
image = pp.Group(
    pp.Literal('image') + pp.Group(pint * 3) +
    pp.Group(pp.OneOrMore(hexa))
).setName('image')
image.parseString("""
image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
             0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
             0x7700FF77 0x444444FF
""")

pixel_texturei = (
    pp.Literal('PixelTexture')  +
    dict_open +
    image +
    dict_close).setName('pixel_texture')

pixel_texturei.parseString("""
PixelTexture {
    image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                 0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                 0x7700FF77 0x444444FF
}
""")


# url "http://www.rt.cs.boeing.com/people/davidk/wrl/geo/colors.jpg"
# repeatS FALSE
# repeatT FALSE
url = (pp.Literal('url') + pp.quotedString).setName('url')
repeat_s = (pp.Literal('repeatS') + boolean).setName('repeat_s')
repeat_t = (pp.Literal('repeatT') + boolean).setName('repeat_t')

image_texture_values = pp.Group(pp.OneOrMore(url | repeat_s | repeat_t))
image_texturei = (
    pp.Literal('ImageTexture') +
    dict_open +
    image_texture_values +
    dict_close).setName('image_texture')


texture_types = pixel_texturei | image_texturei
texture = (
    pp.Literal('texture') + pp.Literal('DEF').suppress() + pword +
    texture_types).setName('texture')

texture.parseString("""
texture DEF PICBAND ImageTexture {
   url "http://www.rt.cs.boeing.com/people/davidk/wrl/geo/colors.jpg"
   repeatS FALSE
   repeatT FALSE
}
""")


texture.parseString("""
texture DEF PICBAND PixelTexture {
    image 1 10 4 0xFFFFFF77 0xFF0000FF 0xFFCC0077 0xFFFF00FF
                 0x77FF00FF 0x00FF00FF 0x00FFFFFF 0x0000FFFF
                 0x7700FF77 0x444444FF
}
""")
#-----------------------------------------
point3d = (pp.Literal('point') + list_open + pp.Group(pp.OneOrMore(xyz)) + list_close).setName('point')
point2d = (pp.Literal('point') + list_open + pp.Group(pp.OneOrMore(xy)) + list_close).setName('point')
coord_values = point3d
coord = (pp.Literal('coord') + pp.Literal('Coordinate') + dict_open + pp.Group(coord_values) + dict_close).setName('coord')
coord.parseString("""
coord Coordinate {
    point [
        3.303 -6.738 -16.931,  3.275 -6.738 -16.932,  3.285 -6.821 -17.012,
        3.641 -6.636 -16.832,  3.642 -6.624 -16.82,  3.509 -6.622 -16.819,
        2.885 -7.116 -17.299,  3.019 -7.116 -17.299
    ]
}
""")

vector = (pp.Literal('vector') + list_open + pp.Group(pp.OneOrMore(xyz)) + list_close).setName('vector')
normal_values = vector
normal = (pp.Literal('normal') + pp.Literal('Normal') + dict_open + pp.Group(normal_values) + dict_close).setName('normal')

vector.parseString("""
vector [
    0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
    0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
]
""")

normal.parseString("""
normal Normal {
    vector [
        0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
        0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
    ]
}
""")


# this shoud be (1, 2, 3, -1), (1, 2, 3, -1), etc.
def cast_to_ints(args):
    ints = np.array(list(args), dtype='int32')
    return ints

if 1:
    # 54.56 sec, 52.45 sec, 49, 35 by moving numpy import
    coord_indicies = pp.OneOrMore(comma.suppress() | pint.setParseAction(cvt_int) | pminus1.setParseAction(cvt_int)) # works and parses
    coord_index = (pp.Literal('coordIndex') + list_open + pp.Group(coord_indicies) + list_close).setName('coord_index')  # works wiht A
elif 0:  # pragma: no cover
    # back to 51.2 sec
    coord_indicies = pp.OneOrMore(comma.suppress() | pp.Word(pp.nums + '-').setParseAction(cvt_int)) # works and parses
    coord_index = (pp.Literal('coordIndex') + list_open + pp.Group(coord_indicies) + list_close).setName('coord_index')  # works wiht A
elif 0:  # pragma: no cover
    # has issues with the big problem
    coord_indicies = pp.delimitedList(pint.setParseAction(cvt_int) | pminus1.setParseAction(cvt_int)) # good
    coord_index = (pp.Literal('coordIndex') + list_open + coord_indicies + list_close).setName('coord_index')
else:  # pragma: no cover
    # has issues with the big problem
    # probably will be benefical in other cases
    import numpy as np
    coord_indicies = pp.delimitedList(pp.Word(pp.nums + '-')).setParseAction(cast_to_ints)  # single numpy array cast
    #coord_indicies = pp.pyparsing_common.comma_separated_list # bad...
    #coord_indicies = OneOrMore(comma.suppress() | pint.setParseAction(cvt_int)) + pminus1.setParseAction(cvt_int)))

coord_indicies.parseString("0, 1, 2, -1, 3, 4, 5, -1, 1, 6, 2, -1")
coord_index.parseString("coordIndex [0, 1, 2, -1, 3, 4, 5, -1, 1, 6, 2, -1]")
#aaab

crease_angle = (pp.Literal('creaseAngle') + pfloat).setName('crease_angle')
tex_coord = (pp.Literal('texCoord') + pp.Literal('TextureCoordinate') + dict_open + point2d + dict_close).setName('tex_coord')

index_face_set_values = pp.OneOrMore(crease_angle | coord | normal | coord_index | tex_coord)
index_face_set = (
    pp.Literal('IndexedFaceSet') + dict_open +
    pp.Group(index_face_set_values) + dict_close).setName('indexed_face_set')
#-----------------------------------------
appearance_values = texture | material
appearance = (
    pp.Literal('appearance') + pp.Literal('Appearance')
    + dict_open + pp.Group(appearance_values) + dict_close).setName('appearance')
sphere = (pp.Literal('Sphere') + dict_open + dict_close).setName('sphere')

geometry_values = sphere | index_face_set
geometry = (pp.Literal('geometry') + geometry_values).setName('geometry')
shape_values = pp.Group(pp.OneOrMore(appearance | geometry))
shape = (pp.Literal('Shape') + dict_open + shape_values + dict_close).setName('shape')

#print(geometry.parseString("""
#geometry IndexedFaceSet {
    #creaseAngle 0.1
    #coord Coordinate {
         #[
            #3.303 -6.738 -16.931,  3.275 -6.738 -16.932,  3.285 -6.821 -17.012,
            #3.641 -6.636 -16.832,  3.642 -6.624 -16.82,  3.509 -6.622 -16.819,
            #2.885 -7.116 -17.299,  3.019 -7.116 -17.299
        #]
    #}
    #normal Normal {
        #vector [
            #0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
            #0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
        #]
    #}
    #coordIndex [
        #0, 1, 2, -1, 3, 4, 5, -1, 1, 6, 2, -1,
    #]
#}
#"""))

#print('geometry...')
# crease_angle + coord + normal + coord_index
geometry.parseString("""
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
""")


shape.parseString("""
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

# wow...super dead link; google has info from him back in 1994...
shape.parseString("""
Shape{
    appearance Appearance{
         texture DEF PICBAND ImageTexture {
             url "http://www.rt.cs.boeing.com/people/davidk/wrl/geo/colors.jpg"
             repeatS FALSE
             repeatT FALSE
         }
    }
}
""")

appearance.parseString("""
    appearance Appearance {
        material Material {
            ambientIntensity 0.210
            diffuseColor 1.000 0.000 0.000
            specularColor 0.500 0.500 0.500
            transparency 0.000
            shininess 0.600
        }
    }
""")
shape.parseString("""
Shape {
    appearance Appearance {
        material Material {
            ambientIntensity 0.210
            diffuseColor 1.000 0.000 0.000
            specularColor 0.500 0.500 0.500
            transparency 0.000
            shininess 0.600
        }
    }
}
""")

#---------------------------------------
shape.parseString("""
Shape {
    appearance Appearance {
        material Material {
            ambientIntensity 0.210
            diffuseColor 1.000 0.000 0.000
            specularColor 0.500 0.500 0.500
            transparency 0.000
            shininess 0.600
        }
    }
    geometry IndexedFaceSet {
        creaseAngle 0.1
        coord Coordinate {
            point [
                3.303 -6.738 -16.931,  3.275 -6.738 -16.932,  3.285 -6.821 -17.012,
                3.641 -6.636 -16.832,  3.642 -6.624 -16.82,  3.509 -6.622 -16.819,
                2.885 -7.116 -17.299,  3.019 -7.116 -17.299
            ]
        }
        normal Normal {
            vector [
                0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
                0 -0.697 0.717,  0 -0.697 0.717,  0 -0.697 0.717,
            ]
        }
        coordIndex [
            0, 1, 2, -1, 3, 4, 5, -1, 1, 6, 2, -1,
            ]
        }
    }
}
""")
#print('done with pre-parsing!')
