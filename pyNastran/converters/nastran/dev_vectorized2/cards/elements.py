from __future__ import print_function

class Elements(object):
    """stores all the elements"""
    def __init__(self, model):
        self.model = model

        self.springs = model.springs
        self.dampers = model.dampers
        self.bushes = model.bushes
        self.rods = model.rods

        self.bars = model.bars
        self.beams = model.beams

        self.shells = model.shells
        self.shears = model.shears
        self.solids = model.solids
        self.eids = []

    def make_current(self):
        elems = [self.springs, self.dampers, self.bushes, self.rods,
                 self.bars, self.beams,
                 self.shells, self.shells, self.shears, self.solids]
        self.eids = []
        for elem in elems:
            elem.make_current()

    def repr_indent(self, indent=''):
        indent = ''
        msg = '%s<Elements> : nelements=%s\n' % (indent, len(self))

        msg += '%s  springs:  %s\n' % (indent, len(self.springs))
        msg += self.springs.repr_indent('    ')

        msg += '%s  dampers:  %s\n' % (indent, len(self.dampers))
        msg += self.dampers.repr_indent('    ')

        msg += '%s  bushes:  %s\n' % (indent, len(self.bushes))
        msg += self.bushes.repr_indent('    ')

        msg += '%s  rods:  %s\n' % (indent, len(self.rods))
        msg += self.rods.repr_indent('    ')

        msg += '%s  bars:  %s\n' % (indent, len(self.bars))
        msg += self.bars.repr_indent('    ')
        msg += '%s  beams:  %s\n' % (indent, len(self.beams))
        msg += self.beams.repr_indent('    ')

        msg += '%s  shells:  %s\n' % (indent, len(self.shells))
        msg += self.shells.repr_indent('    ')

        msg += '%s  shears:  %s\n' % (indent, len(self.shears))
        msg += self.shears.repr_indent('    ')

        msg += '%s  solids:  %s\n' % (indent, len(self.solids))
        msg += self.solids.repr_indent('    ')
        return msg

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.springs.write_card(size, is_double, bdf_file)  # celas
        self.dampers.write_card(size, is_double, bdf_file)  # cdamp
        self.bushes.write_card(size, is_double, bdf_file)  # cbush, cbush1d, cbush2d
        self.rods.write_card(size, is_double, bdf_file)   # crod, conrod, ctube
        self.bars.write_card(size, is_double, bdf_file)   # cbar
        self.beams.write_card(size, is_double, bdf_file)  # cbeam
        self.shears.write_card(size, is_double, bdf_file) # cshear
        self.shells.write_card(size, is_double, bdf_file) # cquad4, ctria3, cquad8, ctria6, cquad
        self.solids.write_card(size, is_double, bdf_file) # ctetra, cpenta, chexa, cpyram

    def get_objects_and_lengths(self):
        """I don't love this as is.  The idea is that it's easier to loop over elements."""
        return (
            (len(self.springs.celas1), self.springs.celas1),
            (len(self.springs.celas2), self.springs.celas2),
            (len(self.springs.celas3), self.springs.celas3),
            (len(self.springs.celas4), self.springs.celas4),
            (len(self.dampers.cdamp1), self.dampers.cdamp1),
            (len(self.dampers.cdamp2), self.dampers.cdamp2),
            (len(self.dampers.cdamp3), self.dampers.cdamp3),
            (len(self.dampers.cdamp4), self.dampers.cdamp4),
            #(len(self.dampers.cdamp5), self.dampers.cdamp5),
            (len(self.bushes.cbush), self.bushes.cbush),
            #(len(self.bushes.cbush1d), (self.bushes.cbush1d),
            #(len(self.bushes.cbush2d), (self.bushes.cbush2d),
            (len(self.rods.conrod), self.rods.conrod),
            (len(self.rods.crod), self.rods.crod),
            (len(self.rods.crod), self.rods.crod),
            (len(self.bars.cbar), self.bars.cbar),
            (len(self.beams.cbeam), self.beams.cbeam),
            (len(self.shells.ctria3), self.shells.ctria3),
            (len(self.shells.cquad4), self.shells.cquad4),
            (len(self.shells.ctria6), self.shells.ctria6),
            (len(self.shells.cquad8), self.shells.cquad8),
            #(len(self.shells.cquad), self.shells.cquad),
            (len(self.shears.cshear), self.shears.cshear),
            (len(self.solids.ctetra4), self.solids.ctetra4),
            (len(self.solids.ctetra10), self.solids.ctetra10),
            (len(self.solids.cpenta6), self.solids.cpenta6),
            (len(self.solids.cpenta15), self.solids.cpenta15),
            (len(self.solids.chexa8), self.solids.chexa8),
            (len(self.solids.chexa20), self.solids.chexa20),
        )

    def __len__(self):
        return (len(self.springs) + len(self.dampers) + len(self.bushes) +
                len(self.rods) + len(self.bars) + len(self.beams) +
                len(self.shells) + len(self.shears) +
                len(self.solids))

    def __repr__(self):
        return self.repr_indent('')
