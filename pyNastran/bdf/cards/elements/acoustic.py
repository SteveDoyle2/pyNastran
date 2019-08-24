from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.field_writer_8 import print_card_8


class CHACAB(BaseCard):
    """
    Acoustic Absorber Element Connection
    Defines the acoustic absorber element in coupled fluid-structural analysis.

    | CHACAB | EID | PID | G1  | G2  | G3  | G4  | G5 | G6 |
    |        | G7  | G8  | G9  | G10 | G11 | G12 |    |    |
    |        |     |     | G17 | G18 | G19 | G20 |    |    |
    """
    type = 'CHACAB'
    def __init__(self, eid, pid, nodes):
        """
        EID Element identification number. (0 < Integer < 100,000,000)
        PID Property identification number of a PACABS entry. (Integer > 0)
        Gi Grid point identification numbers of connection points. (Integer > 0 or blank)
        """
        self.eid = eid
        self.pid = pid
        self.nodes = nodes

    def raw_fields(self):
        nodes1 = self.nodes[:12]
        nodes2 = self.nodes[16:]
        assert len(nodes1) == 12, len(nodes1)
        assert len(nodes2) == 4, len(nodes2)
        return ['CHACAB', self.eid, self.pid, ] + nodes1 + [None, None, None, None] + nodes2

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class CHACBR(BaseCard):
    """

    | CHACBR | EID | PID | G1  | G2  | G3  | G4  | G5 | G6 |
    |        | G7  | G8  | G9  | G10 | G11 | G12 |    |    |
    |        |     |     | G17 | G18 | G19 | G20 |    |    |
    """
    type = 'CHACBR'
    def __init__(self, eid, pid, nodes):
        self.eid = eid
        self.pid = pid
        self.nodes = nodes

    def raw_fields(self):
        nodes1 = self.nodes[:12]
        nodes2 = self.nodes[16:]
        assert len(nodes1) == 12, len(nodes1)
        assert len(nodes2) == 4, len(nodes2)
        return ['CHACBR', self.eid, self.pid, ] + nodes1 + [None, None, None, None] + nodes2

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class CAABSF(BaseCard):
    """
    Frequency-Dependent Acoustic Absorber Element

    Defines a frequency-dependent acoustic absorber element in coupled fluid-structural
    analysis.
    CAABSF EID PID G1 G2 G3 G4
    """
    type = 'CAABSF'
    def __init__(self, eid, pid, nodes):
        self.eid = eid
        self.pid = pid
        self.nodes = nodes

    def raw_fields(self):
        nodes1 = self.nodes[:12]
        nodes2 = self.nodes[16:]
        assert len(nodes1) == 12, len(nodes1)
        assert len(nodes2) == 4, len(nodes2)
        return ['CHACBR', self.eid, self.pid, ] + self.nodes

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)
