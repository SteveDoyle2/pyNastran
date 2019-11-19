from itertools import chain

class Var:
    def __init__(self, name, var, vartype='lineedit', pulldown_objs=None,
                 pulldown_type_limit=None, enabled=True, pulldown_allow_zero=False, required=True):
        self.name = name
        self.var = var
        self.vartype = vartype
        self.pulldown_objs = pulldown_objs
        self.pulldown_type_limit = pulldown_type_limit
        self.pulldown_allow_zero = pulldown_allow_zero
        self.enabled = enabled
        self.required = required
        assert vartype in ['lineedit', 'lineedit_table', 'pulldown', 'spinner'], vartype
    def __repr__(self):
        return (f'Var(name={self.name}, var={self.var}, vartype={self.vartype}, '
                f'pulldown_objs={self.pulldown_objs}, enabled={self.enabled}, required={self.required})')

class TransposedVars:
    def __init__(self, variables):
        self.variables = variables

shell_ptypes = ['PSHELL', 'PCOMP']
bar_ptypes = ['PBAR', 'PBARL']
beam_ptypes = ['PBEAM', 'PBEAML', 'PBCOMP']
ELEMENTS_MAP = {
    'CQUAD4' : [
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='properties',
            pulldown_type_limit=shell_ptypes),
        Var('Nodes', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        #Var('Node 1', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        #Var('Node 2', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        #Var('Node 3', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        #Var('Node 4', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        Var('Theta/Material Coord', 'theta_mcid', vartype='lineedit'),
        Var('Z Offset', 'zoffset', vartype='lineedit'),
    ],
    'CTRIA3' : [
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='properties',
            pulldown_type_limit=shell_ptypes),
        Var('Nodes', 'nodes', vartype='pulldown', pulldown_objs='nodes'),
        Var('Theta/Material Coord', 'theta_mcid', vartype='lineedit'),
        Var('Z Offset', 'zoffset', vartype='lineedit'),
    ],
    'CBAR' : [
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='properties',
            pulldown_type_limit=bar_ptypes, ),
        Var('Node 1', 'ga', vartype='pulldown', pulldown_objs='nodes'),
        Var('Node 2', 'gb', vartype='pulldown', pulldown_objs='nodes'),
        Var('OFFT', 'offt', vartype='pulldown',
            pulldown_objs=['GGG', 'GOO', 'BOO']),
        Var('g0', 'g0', vartype='pulldown', pulldown_objs='nodes', required=False),
        Var('x', 'x', vartype='lineedit'),
        Var('wa', 'wa', vartype='lineedit'),
        Var('wb', 'wa', vartype='lineedit'),
    ],
    'CBEAM' : [
        #bit    : None
        #is_bit : False
        #is_offt : True
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='properties',
            pulldown_type_limit=beam_ptypes, ),
        Var('Node 1', 'ga', vartype='pulldown', pulldown_objs='nodes'),
        Var('Node 2', 'gb', vartype='pulldown', pulldown_objs='nodes'),
        Var('OFFT', 'offt', vartype='pulldown',
            pulldown_objs=['GGG', 'GOO', 'BOO']),
        Var('g0', 'g0', vartype='pulldown', pulldown_objs='nodes', required=False),
        Var('x', 'x', vartype='lineedit'),
        Var('Pin A, pa', 'pa', vartype='lineedit'),
        Var('Pin B, pa', 'pa', vartype='lineedit'),
        Var('Warping A, sa', 'sa', vartype='lineedit'),
        Var('Warping B, sb', 'sb', vartype='lineedit'),
        Var('wa', 'wa', vartype='lineedit'),
        Var('wb', 'wa', vartype='lineedit'),
    ],
}

MASSES_MAP = {
    'CONM2' : [
        Var('Element ID', 'eid', enabled=False),
        Var('Node', 'nid', vartype='pulldown', pulldown_objs='nodes'),
        Var('Coord', 'cid', vartype='pulldown', pulldown_objs='coords'),
        Var('Offset', 'X', vartype='lineedit'),

        Var('Mass', 'mass', vartype='lineedit'),
        Var('I, Inertia', 'I', vartype='lineedit'),
    ],
}

bar_types = [
    'ROD', 'TUBE', 'TUBE2', 'I', 'CHAN', 'T', 'BOX', 'BAR', 'CROSS', 'H', 'T1',
    'I1', 'CHAN1', 'Z', 'CHAN2', 'T2', 'BOX1', 'HEXA', 'HAT', 'HAT1', 'DBOX', 'L']
#shell_mtypes = ['MAT1', 'MAT2', 'MAT8']
PROPERTIES_MAP = {
    'PSHELL' : [
        Var('Property ID', 'pid', enabled=False),
        [
            Var('Material ID 1', 'mid1', pulldown_objs='materials', required=False),
            Var('Thickness', 't', vartype='lineedit', required=False),
        ],
        [
            Var('Material ID 2', 'mid2', pulldown_objs='materials', required=False),
            Var('12I/t^3', 'twelveIt3', vartype='lineedit', required=False),
        ],
        [
            Var('Material ID 3', 'mid3', pulldown_objs='materials', required=False),
            Var('ts/t', 'tst', vartype='lineedit', required=False),
        ],
        Var('Material ID 4', 'mid4', pulldown_objs='materials', required=False),
        Var('z1', 'z1', vartype='lineedit', required=False),
        #Var('z2', 'z2', vartype='lineedit', required=False),
        Var('nsm', 'nsm', vartype='lineedit', required=False),
    ],
    'PCOMP' : [
        Var('Property ID', 'pid', enabled=False),
        TransposedVars([
            Var('Material ID', 'mids', vartype='pulldown', pulldown_objs='materials'),
            Var('Material Angle, Theta', 'thetas'),
            Var('Thickness', 'thicknesses'),
            Var('SOUT', 'souts', vartype='pulldown', pulldown_objs=['YES', 'NO']),
        ]),
        Var('TRef', 'tref', vartype='lineedit'),
        Var('lam', 'lam', vartype='lineedit', required=False),
        Var('ft', 'ft', vartype='lineedit', required=False),
        Var('sb', 'sb', vartype='lineedit', required=False),
        Var('z0', 'z0', vartype='lineedit'),
        Var('Damping', 'ge', vartype='lineedit'),
        Var('nsm', 'nsm', vartype='lineedit', required=False),
    ],
    'PBARL' : [
        Var('Property ID', 'pid', enabled=False),
        Var('Material ID 1', 'mid', vartype='pulldown', pulldown_objs='materials',
            pulldown_type_limit=['MAT1'], ),
        Var('Bar Type', 'beam_type', vartype='pulldown', pulldown_objs=bar_types, enabled=False),
        Var('Dimensions', 'dim', vartype='lineedit', required=False),
        Var('Group', 'group', vartype='lineedit', required=False),
        Var('nsm', 'nsm', vartype='lineedit', required=False),
    ],
    'PBAR' : [
        Var('Property ID', 'pid', enabled=False),
        Var('Material ID', 'mid', vartype='pulldown', pulldown_objs='materials',
            pulldown_type_limit=['MAT1'], ),
        [
            Var('Area', 'A', vartype='lineedit'),
            Var('I1', 'i1', vartype='lineedit'),
            Var('I2', 'i2', vartype='lineedit'),
            Var('I12', 'i12', vartype='lineedit'),
            Var('J', 'j', vartype='lineedit'),
        ],
        [
            Var('C1', 'c1', vartype='lineedit'),
            Var('D1', 'd1', vartype='lineedit'),
            Var('E1', 'e1', vartype='lineedit'),
            Var('F1', 'f1', vartype='lineedit'),
        ],
        [
            Var('C2', 'c2', vartype='lineedit'),
            Var('D2', 'd2', vartype='lineedit'),
            Var('E2', 'e2', vartype='lineedit'),
            Var('F2', 'f2', vartype='lineedit'),
        ],
        [
            Var('K1', 'k1', vartype='lineedit'),
            Var('K2', 'k2', vartype='lineedit'),
        ],
        Var('nsm', 'nsm', vartype='lineedit', required=False),
    ],
    'PBEAML' : [
        Var('Property ID', 'pid', enabled=False),
        Var('Material ID', 'mid', vartype='pulldown', pulldown_objs='materials',
            pulldown_type_limit=['MAT1'], ),
        Var('Beam Type', 'beam_type', vartype='pulldown', pulldown_objs=bar_types, enabled=False),
        Var('Group', 'group', vartype='lineedit', required=False),

        #TransposedVars([
        Var('SOUT', 'so', vartype='pulldown', pulldown_objs=['YES', 'NO']),
        Var('x/xb', 'xxb', vartype='lineedit'),
        Var('nsm', 'nsm', vartype='lineedit', required=False),
        Var('Dimensions', 'dim', vartype='lineedit_table'),
    ],
}


MATERIALS_MAP = {
    'MAT1' : [
        Var('Material ID', 'mid', enabled=False),
        Var("E, Young's Modulus", 'e', vartype='lineedit', required=False),
        Var("G, Shear Modulus", 'g', vartype='lineedit', required=False),
        Var("nu, Poisson's Ratio", 'nu', vartype='lineedit', required=False),
        Var('TRef', 'tref', vartype='lineedit'),
        Var('Damping', 'ge', vartype='lineedit'),
        Var('Density', 'rho', vartype='lineedit'),
        Var('Material Coord', 'mcsid', vartype='lineedit'),
    ],
    'MAT8' : [
        Var('Material ID', 'mid', enabled=False),
        Var("E11", 'e11', vartype='lineedit'),
        Var("E22", 'e22', vartype='lineedit'),
        [
            Var("G12", 'g12', vartype='lineedit'),
            Var("G1z", 'g1z', vartype='lineedit'),
            Var("G2z", 'g2z', vartype='lineedit'),
        ],
        Var("nu12", 'nu12', vartype='lineedit'),

        Var('Density', 'rho', vartype='lineedit'),
        Var('Damping', 'ge', vartype='lineedit'),

        Var('TRef', 'tref', vartype='lineedit'),
        Var('A1', 'a1', vartype='lineedit'),
        Var('A2', 'a2', vartype='lineedit'),
        [
            Var('Xt', 'Xt', vartype='lineedit'),
            Var('Xc', 'Xc', vartype='lineedit'),
            Var('Yt', 'Yt', vartype='lineedit'),
            Var('Yc', 'Yc', vartype='lineedit'),
            Var('S', 'S', vartype='lineedit'),
            Var('F12', 'F12', vartype='lineedit'),
        ],
        Var('strn', 'strn', vartype='lineedit'),
    ],
}

CAERO_MAP = {
    'CAERO1': [
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='paeros'),
        Var('iGroup', 'igroup'),
        [
            Var('nSpan Boxes', 'nspan', vartype='spinner'),
            Var('AEFACT Span', 'lspan', vartype='pulldown', pulldown_objs='aefacts', pulldown_allow_zero=True),
        ],
        [
            Var('nChord Boxes', 'nchord', vartype='spinner'),
            Var('AEFACT Chord', 'lchord', vartype='pulldown', pulldown_objs='aefacts', pulldown_allow_zero=True),
        ],
        Var('Point 1', 'p1'),
        Var('Distance 12', 'x12'),
        Var('Point 4', 'p4'),
        Var('Distance 43', 'x43'),
    ],
}

MODIFY_MAP = dict(chain(
    ELEMENTS_MAP.items(),
    MASSES_MAP.items(),
    PROPERTIES_MAP.items(),
    MATERIALS_MAP.items(),
    CAERO_MAP.items(),
))

UPDATE_MAP = {
    #'GRID' : 'update_grid',
    #'CONROD' : 'update_element',
    'CAERO1' : 'update_caeros',
}
