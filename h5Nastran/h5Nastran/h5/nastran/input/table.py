from __future__ import print_function, absolute_import

from six.moves import range

from h5Nastran.h5nastrannode import H5NastranNode
from .input_table import InputTable, TableDef


class Table(H5NastranNode):
    def __init__(self, h5n, input):
        self._h5n = h5n
        self._input = input

        self.mkaero1 = MKAERO1(self._h5n, self)
        self.tabdmp1 = TABDMP1(self._h5n, self)
        self.tabled1 = TABLED1(self._h5n, self)
        self.tabled2 = TABLED2(self._h5n, self)
        self.tabled3 = TABLED3(self._h5n, self)
        self.tabled4 = TABLED4(self._h5n, self)
        self.tablem1 = TABLEM1(self._h5n, self)
        self.tablem2 = TABLEM2(self._h5n, self)
        self.tablem3 = TABLEM3(self._h5n, self)
        self.tablem4 = TABLEM4(self._h5n, self)

    def path(self):
        return self._input.path() + ['TABLE']


########################################################################################################################


class MKAERO1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/MKAERO1')

    def from_bdf(self, cards):
        result = {
            'IDENTITY': {
                'M1': [], 'K1': [], 'U2': [], 'M2': [], 'K2': [], 'U3': [], 'M3': [], 'K3': [],
                'U4': [], 'M4': [], 'K4': [], 'U5': [], 'M5': [], 'K5': [], 'U6': [], 'M6': [], 'K6': [],
                'U7': [], 'M7': [], 'K7': [], 'U8': [], 'M8': [], 'K8': []
            }
        }

        data = result['IDENTITY']

        m1 = data['M1']
        k1 = data['K1']
        u2 = data['U2']
        m2 = data['M2']
        k2 = data['K2']
        u3 = data['U3']
        m3 = data['M3']
        k3 = data['K3']
        u4 = data['U4']
        m4 = data['M4']
        k4 = data['K4']
        u5 = data['U5']
        m5 = data['M5']
        k5 = data['K5']
        u6 = data['U6']
        m6 = data['M6']
        k6 = data['K6']
        u7 = data['U7']
        m7 = data['M7']
        k7 = data['K7']
        u8 = data['U8']
        m8 = data['M8']
        k8 = data['K8']

        m = [m1, m2, m3, m4, m5, m6, m7, m8]
        k = [k1, k2, k3, k4, k5, k6, k7, k8]
        u = [None, u2, u3, u4, u5, u6, u7, u8]

        for card in cards:
            m_ = card.machs.tolist()
            m_len = len(m_)
            if m_len < 8:
                m_ += [0.] * (8 - m_len)

            k_ = card.reduced_freqs.tolist()
            k_len = len(k_)
            if k_len < 8:
                k_ += [0.] * (8 - k_len)

            min_i = min(m_len, k_len) - 1

            for i in range(8):
                m[i].append(m_[i])
                k[i].append(k_[i])

                # TODO: MKAERO1 - verify u's are correct, m and k are not always the same length
                #       but the msc spec implies that they should be

                if i > 0:
                    if i <= min_i:
                        u[i].append(1)
                    else:
                        u[i].append(0)

        return result


########################################################################################################################


# TODO: TABDMP1 is missing Type in msc spec... why?

class TABDMP1_DAMP(object):
    name = 'DAMP'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('F', '<f8', ()), ('G', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


class TABDMP1_SPEC(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S4', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = [TABDMP1_DAMP]


class TABDMP1(InputTable):
    table_def = TableDef.create(TABDMP1_SPEC, rename={'DAMP_POS': 'POS', 'DAMP_LEN': 'LEN'})

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        damp = {'IDENTITY': {'F': [], 'G': []}}

        result = {'IDENTITY': {'ID': [], 'TYPE': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'DAMP': damp,
                  '_subtables': ['DAMP']}

        f = damp['IDENTITY']['F']
        g = damp['IDENTITY']['G']
        identity = result['IDENTITY']
        id_ = identity['ID']
        type_ = identity['TYPE']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            type_.append(card.Type)
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            f += list(card.x)
            g += list(card.y)

        return result


########################################################################################################################


class TABLED1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED1/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        xy = {'IDENTITY': {'X': [], 'Y': []}}

        result = {'IDENTITY': {'ID': [], 'CODEX': [], 'CODEY': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'XY': xy,
                  '_subtables': ['XY']}

        x = xy['IDENTITY']['X']
        y = xy['IDENTITY']['Y']
        identity = result['IDENTITY']
        id_ = identity['ID']
        codex = identity['CODEX']
        codey = identity['CODEY']
        pos = identity['POS']
        len_ = identity['LEN']

        # TODO: TABLED1 confirm codes are correct
        _xy_code = {'LINEAR': 1, 'LOG': 2, '': 1, None: 1}

        _pos = 0

        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            codex.append(_xy_code[card.xaxis])
            codey.append(_xy_code[card.yaxis])
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            _pos += _len
            x.extend(list(card.x))
            y.extend(list(card.y))

        return result

########################################################################################################################


class TABLED2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED2/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        xy = {'IDENTITY': {'X': [], 'Y': []}}

        result = {'IDENTITY': {'ID': [], 'X1': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'XY': xy,
                  '_subtables': ['XY']}

        x = xy['IDENTITY']['X']
        y = xy['IDENTITY']['Y']
        identity = result['IDENTITY']
        id_ = identity['ID']
        x1 = identity['X1']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            x1.append(card.x1)
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            _pos += _len
            x.extend(list(card.x))
            y.extend(list(card.y))

        return result

########################################################################################################################


class TABLED3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED3/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        xy = {'IDENTITY': {'X': [], 'Y': []}}

        result = {'IDENTITY': {'ID': [], 'X1': [], 'X2': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'XY': xy,
                  '_subtables': ['XY']}

        x = xy['IDENTITY']['X']
        y = xy['IDENTITY']['Y']
        identity = result['IDENTITY']
        id_ = identity['ID']
        x1 = identity['X1']
        x2 = identity['X2']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            x1.append(card.x1)
            x2.append(card.x2)
            pos.append(_pos)
            _len = len(card.x)
            len_.append(_len)
            _pos += _len
            x.extend(list(card.x))
            y.extend(list(card.y))

        return result

########################################################################################################################


class TABLED4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLED4/IDENTITY', rename={'COEF_POS': 'POS', 'COEF_LEN': 'LEN'})

    def from_bdf(self, cards):
        card_ids = sorted(cards.keys())

        coef = {'IDENTITY': {'A': []}}

        result = {'IDENTITY': {'ID': [], 'X1': [], 'X2': [], 'X3': [], 'X4': [], 'POS': [], 'LEN': [], 'DOMAIN_ID': []},
                  'COEF': coef,
                  '_subtables': ['COEF']}

        a = coef['IDENTITY']['A']
        identity = result['IDENTITY']
        id_ = identity['ID']
        x1 = identity['X1']
        x2 = identity['X2']
        x3 = identity['X3']
        x4 = identity['X4']
        pos = identity['POS']
        len_ = identity['LEN']

        _pos = 0
        for card_id in card_ids:
            card = cards[card_id]

            id_.append(card.tid)
            x1.append(card.x1)
            x2.append(card.x2)
            x3.append(card.x3)
            x4.append(card.x4)
            pos.append(_pos)
            _len = len(card.a)
            len_.append(_len)
            _pos += _len
            a.extend(list(card.a))

        return result

########################################################################################################################


class TABLEM1(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLEM1/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})
    from_bdf = TABLED1.__dict__['from_bdf']


########################################################################################################################


class TABLEM2(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLEM2/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})
    from_bdf = TABLED2.__dict__['from_bdf']

########################################################################################################################


class TABLEM3(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLEM3/IDENTITY', rename={'XY_POS': 'POS', 'XY_LEN': 'LEN'})
    from_bdf = TABLED3.__dict__['from_bdf']


########################################################################################################################


class TABLEM4(InputTable):
    table_def = TableDef.create('/NASTRAN/INPUT/TABLE/TABLEM4/IDENTITY', rename={'COEF_POS': 'POS', 'COEF_LEN': 'LEN'})
    from_bdf = TABLED4.__dict__['from_bdf']

