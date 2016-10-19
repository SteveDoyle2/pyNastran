from six.moves import zip
from numpy import array


class SPCD(object):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).::

      SPCD SID G1 C1 D1   G2 C2 D2
      SPCD 2   32 3  -2.6  5
    """
    type = 'SPC'

    def __init__(self, model):
        self.model = model
        self.n = 0

        self._comments = []
        self.constraint_id = None
        self.grid_id = []
        self.components = []

    def add(self, i, constraint_id, node_id, dofs, enforced_motion, comment):
        if i == 0:
            n = 2
        elif i == 1:
            n = 5
        else:
            raise RuntimeError('i =', i)

        assert enforced_motion != 0.0

        self._comments.append(comment)
        self.grid_id.append(node_id)
        self.components.append(dofs)
        self.enforced_motion.append(enforced_motion)

        if self.constraint_id is None:
            self.constraint_id = constraint_id
        else:
            msg = ('self.constraint_id == constraint_id; constraint_id=%r '
                   'expected; found=%r' % (self.constraint_id. constraint_id))
            raise RuntimeError(msg)
            #assert self.constraint_id == constraint_id, ''

    def build(self):
        #float_fmt = self.model.float_fmt
        self.n = len(self.self.constraint_id)

        self.grid_id = array(self.grid_id)
        self.components = array(self.components)
        self.enforced_motion = array(self.enforced_motion)

    def write_card(self, f, size=8):
        if self.n:
            fields = ['SPCD', self.constraint_id]
            for (nid, constraint, enforced) in zip(self.gids, self.components, self.enforced_motion):
                fields += [nid, constraint, enforced]
            f.write(fields)
