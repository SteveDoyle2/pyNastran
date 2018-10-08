from __future__ import print_function
import sys

def read_avl(avl_filename):
    avl = AVL()
    avl.read_avl(avl_filename)

class AVL(object):
    def __init__(self):
        pass

    def read_avl(self, avl_filename):
        with open(avl_filename, 'r') as avl_file:
            lines = avl_file.readlines()

        # remove comments and whitespace lines
        lines = [line.rstrip() for line in lines if line.split('#')[0].split('!')[0].strip()]

        name = lines[0].strip()
        mach = float(lines[1])
        sym_iy, sym_iz, symz = lines[2].split()
        sref, cref, bcref = lines[3].split()
        xref, yref, zref = lines[4].split()

        i = 5
        surfaces = []
        surface = {}
        sections = []
        while i < len(lines):
            line = lines[i]
            if line == 'SURFACE':
                if surface:
                    surfaces.append(surface)
                    surface = {}
                    sections = []
                #print('---------------')
                #print
                name = lines[i+1]
                #print('name = %r' % name)
                #print(lines[i+2])
                sline2 = lines[i+2].split()
                if len(sline2) == 4:
                    nchord, chord_spacing, nspan, span_spacing = sline2
                    print('nchord=%s chord_spacing=%s nspan=%s span_spacing=%s' % (nchord, chord_spacing, nspan, span_spacing))
                    nspan = int(nspan)
                    span_spacing = float(span_spacing)
                    surface['span'] = (nspan, span_spacing)
                else:
                    nchord, chord_spacing = sline2
                    print('nchord=%s chord_spacing=%s' % (nchord, chord_spacing))
                nchord = int(nchord)
                chord_spacing = float(chord_spacing)
                surface['name'] = name
                surface['chord'] = (nchord, chord_spacing)
                surface['sections'] = sections

                i += 3
                print()
                print('---------------')
                continue
            elif line == 'COMPONENT':
                ncomponents = int(lines[i+1])
                assert ncomponents == 1, ncomponents
                assert 'component' not in surface
                surface['component'] = 1
                i += 1

            elif line == 'YDUPLICATE':
                yduplicate = float(lines[i+1])
                assert 'yduplicate' not in surface
                surface['yduplicate'] = yduplicate
                i += 1

            elif line == 'ANGLE':
                angle = float(lines[i+1])
                assert 'angle' not in surface
                surface['angle'] = angle
                i += 1

            elif line == 'NOWAKE':
                assert 'nowake' not in surface
                surface['nowake'] = True

            elif line == 'SCALE':
                xscale, yscale, zscale = lines[i+1].split()
                xscale = float(xscale)
                yscale = float(yscale)
                zscale = float(zscale)
                assert 'scale' not in surface
                surface['scale'] = (xscale, yscale, zscale)
                i += 1

            elif line == 'TRANSLATE':
                xtranslate, ytranslate, ztranslate = lines[i+1].split()
                xtranslate = float(xtranslate)
                ytranslate = float(ytranslate)
                ztranslate = float(ztranslate)
                assert 'translate' not in surface
                surface['translate'] = (xtranslate, ytranslate, ztranslate)
                i += 1
            elif line == 'SECTION':
                section_data = {
                    'afile' : [],
                    'control' : [],
                }
                sline = lines[i+1].split()
                if len(sline) == 6:
                    # #Xle Yle Zle Chord    Nspanwise Sspace
                    xle, yle, zle, chord, nspan, span_spacing = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    chord = float(chord)
                    nspan = float(nspan)
                    span_spacing = float(span_spacing)
                    assert 'station' not in surface
                    section = [xle, yle, zle, chord, nspan, span_spacing]
                elif len(sline) == 5:
                    # #Xle Yle Zle Chord Ainc
                    xle, yle, zle, chord, ainc = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    chord = float(chord)
                    ainc = float(ainc)
                    section = [xle, yle, zle, chord, ainc]
                elif len(sline) == 7:
                    #Xle   Yle    Zle      Chord   Ainc  Nspanwise  Sspace
                    xle, yle, zle, chord, ainc, nspan, span_spacing = sline
                    xle = float(xle)
                    yle = float(yle)
                    zle = float(zle)
                    ainc = float(ainc)
                    chord = float(chord)
                    nspan = float(nspan)
                    span_spacing = float(span_spacing)
                    assert 'station' not in surface
                    section = [xle, yle, zle, chord, ainc, nspan, span_spacing]

                else:
                    print(sline, len(sline))
                    asdf
                section_data['section'] = section
                sections.append(section_data)
                i += 1
            elif line == 'AFILE':
                section_data['afile'].append(lines[i+1])
                i += 1
            elif line == 'CONTROL':
                sline = lines[i+1].split()
                #print(sline)
                name, afloat, bfloat, cfloat, dfloat, efloat, fint = sline
                section_data['control'].append([name, afloat, bfloat, cfloat, dfloat, efloat, fint])
                i += 1

            #elif line == 'CONTROL':


            else:
                print(line)
            i += 1

        for surface in surfaces:
            sections = surface['sections']
            del surface['sections']
            print(surface)
            for section in sections:
                if not section['afile']:
                    del section['afile']
                if not section['control']:
                    del section['control']
                print('  ', section)
            print('')

def main():
    avl_filename = sys.argv[1]
    read_avl(avl_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
