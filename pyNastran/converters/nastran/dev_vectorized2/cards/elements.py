from __future__ import print_function

class Elements(object):
    """stores all the elements"""
    def __init__(self, model):
        self.model = model

        self.springs = model.springs
        self.dampers = model.dampers
        self.rods = model.rods
        self.bars = model.bars
        self.beams = model.beams
        self.shells = model.shells
        self.shears = model.shears
        self.solids = model.solids
        self.eids = []

    def make_current(self):
        elems = [self.springs, self.dampers, self.rods,
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
        self.rods.write_card(size, is_double, bdf_file)   # crod, conrod, ctube
        self.bars.write_card(size, is_double, bdf_file)   # cbar
        self.beams.write_card(size, is_double, bdf_file)  # cbeam
        self.shears.write_card(size, is_double, bdf_file) # cshear
        self.shells.write_card(size, is_double, bdf_file) # cquad4, ctria3, cquad8, ctria6, cquad
        self.solids.write_card(size, is_double, bdf_file) # ctetra, cpenta, chexa, cpyram

    def __len__(self):
        return (len(self.springs) +  + len(self.dampers) +
                len(self.rods) + len(self.bars) + len(self.beams) +
                len(self.shells) + len(self.shears) +
                len(self.solids))

    def __repr__(self):
        return self.repr_indent('')
