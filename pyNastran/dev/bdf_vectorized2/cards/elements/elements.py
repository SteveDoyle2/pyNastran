class Elements:
    """stores all the elements"""
    def __init__(self, model):
        """creates the Elements object that stores all the elements"""
        self.model = model

        self.springs = model.springs
        self.dampers = model.dampers
        self.bushes = model.bushes
        self.masses = model.masses2
        self.rods = model.rods

        self.bars = model.bars
        self.beams = model.beams

        self.shells = model.shells
        self.shears = model.shears
        self.solids = model.solids
        self.eids = []

    def make_current(self):
        """calls make_current() for each group"""
        self.eids = []
        for group in self.groups:
            group.make_current()

    @property
    def groups(self):
        """gets the sub-element groups"""
        groups = [
            self.springs, self.dampers, self.bushes, self.rods,
            self.bars, self.beams,
            self.shells, self.shears, self.solids]
        return groups

    @property
    def elements(self):
        """gets all the elements in the expected sorted order"""
        elems = []
        for group in self.groups:
            elems += group.elements
        return elems

    def repr_indent(self, indent=''):
        indent = ''
        msg = '%s<Elements> : nelements=%s\n' % (indent, len(self))

        msg += '%s  masses:  %s\n' % (indent, len(self.masses))
        msg += self.masses.repr_indent('    ')

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
        msg += '%s  shears:  %s\n' % (indent, len(self.shears))
        msg += self.shears.repr_indent('    ')

        msg += '%s  shells:  %s\n' % (indent, len(self.shells))
        msg += self.shells.repr_indent('    ')
        msg += '%s  solids:  %s\n' % (indent, len(self.solids))
        msg += self.solids.repr_indent('    ')
        return msg

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        #for group in self.groups:
            #group.write_card(size, is_double, bdf_file)
        self.masses.write_card(size, is_double, bdf_file)  # conm1, conm2, cmass1, cmass2, cmass3, cmass4
        self.springs.write_card(size, is_double, bdf_file)  # celas
        self.dampers.write_card(size, is_double, bdf_file)  # cdamp
        self.bushes.write_card(size, is_double, bdf_file)  # cbush, cbush1d, cbush2d
        self.rods.write_card(size, is_double, bdf_file)   # crod, conrod, ctube
        self.bars.write_card(size, is_double, bdf_file)   # cbar
        self.beams.write_card(size, is_double, bdf_file)  # cbeam
        self.shears.write_card(size, is_double, bdf_file) # cshear
        self.shells.write_card(size, is_double, bdf_file) # cquad4, ctria3, cquad8, ctria6, cquad
        self.solids.write_card(size, is_double, bdf_file) # ctetra, cpenta, chexa, cpyram

    def __len__(self):
        """gets the number of elements"""
        return sum([len(group) for group in self.groups])

    def __repr__(self):
        """simplified printout"""
        return self.repr_indent('')
