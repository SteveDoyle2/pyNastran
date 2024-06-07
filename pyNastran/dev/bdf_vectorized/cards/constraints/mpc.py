from io import StringIO
from pyNastran.bdf.bdf_interface.assign_type import (integer, double_or_blank, components_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

def get_mpc_constraint(card):
    #: Set identification number. (Integer > 0)
    constraint_id = integer(card, 1, 'conid')
    #: Identification number of grid or scalar point. (Integer > 0)
    #gids = []
    #: Component number. (Any one of the Integers 1 through 6 for grid
    #: points; blank or zero for scalar points.)
    #constraints = []
    #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
    #enforced = []
    constraints = []

    fields = card.fields(0)
    nfields = len(fields) - 1
    for iField in range(2, nfields, 8):
        grid = integer(card, iField, 'gid')
        component = components_or_blank(card, iField + 1, 'constraint', '0')  # scalar point
        value = double_or_blank(card, iField + 2, 'enforced', 0.0)

        constraints.append([grid, int(component), value])
        #self.gids.append(grid)
        #self.constraints.append(component)
        #self.enforced.append(value)

        if iField + 3 > nfields:
            break
        grid = integer(card, iField + 3, 'gid')
        component = components_or_blank(card, iField + 4, 'constraint', '0')  # scalar point
        value = double_or_blank(card, iField + 5, 'enforced')

        constraints.append([grid, int(component), value])
        #self.gids.append(grid)
        #self.constraints.append(component)
        #self.enforced.append(value)

    return constraint_id, constraints


class MPC:
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).::

    +-----+-----+----+----+------+----+----+----+-----+
    |  1  |  2  |  3 |  4 |  5   |  6 |  7 |  8 |  9  |
    +=====+=====+====+====+======+====+====+====+=====+
    | MPC | SID | G1 | C1 |  A1  | G2 | C2 | A2 |     |
    +-----+-----+----+----+------+----+----+----+-----+
    |     |  G3 | C3 | A3 | ...  |    |    |    |     |
    +-----+-----+----+----+------+----+----+----+-----+

    +-----+-----+----+----+------+----+----+----+-----+
    | MPC | 2   | 32 | 3  | -2.6 |  5 |    |    |     |
    +-----+-----+----+----+------+----+----+----+-----+
     """
    type = 'MPC'

    def __init__(self, model):
        self.model = model
        self.n = 0

        self._comments = []
        self.constraint_id = None
        self.constraints = []

    def add(self, constraint_id, constraint, comment):
        self._comments.append(comment)
        self.constraints.append(constraint)

        if self.constraint_id is None:
            self.constraint_id = constraint_id
        elif self.constraint_id != constraint_id:
            raise RuntimeError('self.constraint_id == constraint_id; '
                               f'constraint_id={self.constraint_id} expected; found={constraint_id}')

    def build(self):
        #float_fmt = self.model.float_fmt
        self.n = len(self.constraints)

    def write_card(self, bdf_file, size=8):
        assert self.n > 0, self.n
        if self.n:
            for constraint in self.constraints:
                card = ['MPC', self.constraint_id]
                for (G, C, A) in constraint:
                    card += [G, C, A]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue()
