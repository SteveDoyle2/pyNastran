import sys
from numpy import float32


def read_lsdyna(key_filename):
    model = Lsdyna()
    model.read_lsdyna(key_filename)
    return model

class Lsdyna:
    def __init__(self, log=None, debug=False):
        self.log = log
        self.nodes = None

    def read_lsdyna(self, key_filename):
        with open(key_filename) as key_file:
            lines = key_file.readlines()

        lines = remove_comments(lines)
        istop = 0
        iline = 0
        cards = [
            '*CONSTRAINED_EXTRA_NODES_NODE', '*CONSTRAINED_EXTRA_NODES_SET',
            '*CONSTRAINED_JOINT_LOCKING_LOCAL', '*CONSTRAINED_JOINT_REVOLUTE',
            '*CONSTRAINED_JOINT_STIFFNESS_GENERALIZED', '*CONSTRAINED_JOINT_TRANSLATIONAL_LOCAL',
            '*CONSTRAINED_RIGID_BODIES', '*CONTACT_AUTOMATIC_SURFACE_TO_SURFACE',
            '*CONTROL_CONTACT', '*CONTROL_ENERGY', '*CONTROL_OUTPUT', '*CONTROL_PARALLEL',
            '*CONTROL_SHELL', '*CONTROL_TERMINATION', '*CONTROL_TIMESTEP',

            '*DATABASE_BINARY_D3PLOT', '*DATABASE_DEFORC', '*DATABASE_ELOUT',
            '*DATABASE_EXTENT_BINARY', '*DATABASE_GLSTAT', '*DATABASE_HISTORY_NODE',
            '*DATABASE_HISTORY_NODE_LOCAL', '*DATABASE_JNTFORC', '*DATABASE_MATSUM',
            '*DATABASE_NODOUT', '*DATABASE_RCFORC', '*DATABASE_SBTOUT',
            '*DATABASE_SECFORC', '*DATABASE_SLEOUT',

            '*DEFINE_COORDINATE_NODES', '*DEFINE_CURVE', '*DEFINE_SD_ORIENTATION',
            '*ELEMENT_DISCRETE', '*ELEMENT_MASS', '*ELEMENT_SEATBELT',
            '*ELEMENT_SEATBELT_SLIPRING', '*ELEMENT_SOLID', '*HOURGLASS',
            '*MAT_ELASTIC', '*MAT_RIGID', '*MAT_SEATBELT', '*MAT_SPRING_MAXWELL',
            '*MAT_VISCOELASTIC', '*NODE', '*PART', '*SECTION_DISCRETE',
            '*SECTION_SEATBELT', '*SECTION_SOLID', '*SET_NODE_LIST',
            '*SET_SEGMENT', '*TITLE',
        ]
        card_map = {
            '*NODE' : self.parse_node,
        }

        while iline < len(lines):
            line = lines[iline].strip()
            #if line.startswith('$'):
            #    continue
            keyword = line
            if line in ['*KEYWORD']:
                pass
            elif line in cards:
                iline, line, data = get_block(lines, iline)
                if keyword in card_map:
                    card_map[keyword](data)
            elif line == '*END':
                return
            elif line[0] == '*':
                raise NotImplementedError(line)
            else:
                asdf
                print(line)
                istop += 1
            if istop == 10:
                sys.exit()
            iline += 1

    def parse_node(self, lines):
        nnodes = len(lines)
        nids = []
        nodes = []
        for iline, line in enumerate(lines):
            nid, x, y, z = line[:7], line[7:22], line[22:38], line[38:]
            #print(line)
            #print('nid=%r x=%r y=%r z=%r' % (nid, x, y, z))
            #assert len(nid) == 7, len(nid)
            #assert len(x) == 15, len(x)
            #assert len(y) == 16, len(y)
            assert len(z) <= 16, len(z)
            nid = int(nid)
            try:
                x, y, z = double(x, 1, 'x'), double(y, 2, 'y'), double(z, 3, 'z')
            except:
                print(iline)
                print(line)
                raise
            nids.append(nid)
            nodes.append([x, y, z])
        self.nids = nids
        self.nodes = nodes


def double(svalue, ifield=-1, fieldname='NA'):
    """
    Converts a field into a double

    Parameters
    ----------
    value : str
        the value to parse
    fieldname : str
        name of field

    Returns
    -------
    value : float
        the value from the desired field
    """
    if isinstance(svalue, (float, float32)):
        return svalue
    elif isinstance(svalue, int):
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n' % (
            fieldname, svalue, ifield, dtype))
    elif svalue is None or len(svalue) == 0:  ## None
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n'  % (
            fieldname, svalue, ifield, dtype))

    if svalue.isdigit():
        # if only int
        dtype = _get_dtype(int(svalue))
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n' % (
            fieldname, svalue, ifield, dtype))

    #svalue = svalue.strip()
    try:
        # 1.0, 1.0E+3, 1.0E-3
        value = float(svalue)
    except TypeError:
        dtype = _get_dtype(svalue)
        raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n' % (
            fieldname, svalue, ifield, dtype))
    except ValueError:
        # 1D+3, 1D-3, 1-3
        try:
            svalue = svalue.upper()
            if 'D' in svalue:
                # 1.0D+3, 1.0D-3
                svalue2 = svalue.replace('D', 'E')
                return float(svalue2)

            # 1.0+3, 1.0-3
            sign = ''
            if svalue[0] in ('+', '-'):
                sign = svalue[0]
                svalue = svalue[1:]
            if '+' in svalue:
                svalue = sign + svalue.replace('+', 'E+')
            elif '-' in svalue:
                svalue = sign + svalue.replace('-', 'E-')

            value = float(svalue)
        except ValueError:
            dtype = _get_dtype(svalue)
            raise SyntaxError('%s = %r (field #%s) on card must be a float (not %s).\n' % (
                fieldname, svalue, ifield, dtype))
    return value

def _get_dtype(svalue):
    return type(svalue)

def get_block(lines, iline):
    iline += 1
    line = lines[iline]
    data = []
    while '*' not in line:
        data.append(line.rstrip())
        line = lines[iline]
        iline += 1
    iline -= 2
    line = lines[iline]
    assert '*' in lines[iline+1], line
    return iline, line, data

def remove_comments(lines):
    return [line.strip().split('$')[0] for line in lines if line.strip().split('$')[0]]
