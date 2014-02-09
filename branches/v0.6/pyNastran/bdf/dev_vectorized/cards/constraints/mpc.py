import StringIO
from pyNastran.bdf2.bdf_interface.assign_type import (integer, double_or_blank, components_or_blank)
from pyNastran.bdf.fieldWriter import print_card

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
    nFields = len(fields) - 1
    for iField in xrange(2, nFields, 8):
        grid = integer(card, iField, 'gid')
        component = components_or_blank(card, iField + 1, 'constraint', 0)  # scalar point
        value = double_or_blank(card, iField + 2, 'enforced', 0.0)

        constraints.append([grid, int(component), value])
        #self.gids.append(grid)
        #self.constraints.append(component)
        #self.enforced.append(value)

        if iField + 3 > nFields:
            break
        grid = integer(card, iField + 3, 'gid')
        component = components_or_blank(card, iField + 4, 'constraint', 0)  # scalar point
        value = double_or_blank(card, iField + 5, 'enforced')

        constraints.append([grid, int(component), value])
        #self.gids.append(grid)
        #self.constraints.append(component)
        #self.enforced.append(value)

    return constraint_id, constraints


class MPC(object):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).::

      MPC SID G1 C1 D1   G2 C2 D2
      MPC 2   32 3  -2.6  5
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
            raise RuntimeError('self.constraint_id == constraint_id; constraint_id=%r expected; found=%r' % (self.constraint_id. constraint_id))

    def build(self):
        float_fmt = self.model.float
        self.n = len(self.constraints)

    def write_bdf(self, f, size=8):
        if self.n:
            for constraint in self.constraints:
                card = ['MPC', self.constraint_id]
                for (G, C, A) in constraint:
                    card += [G, C, A]
                f.write(print_card(card))

    def __repr__(self):
        f = StringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue()