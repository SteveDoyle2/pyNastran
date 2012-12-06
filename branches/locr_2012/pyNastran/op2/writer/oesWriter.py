from struct import pack


class Oes1Writer(object):
    def writeOES1(self):
        """
        writes isotropic/composite stress/strain
        @todo assumes s_code=0 (stress) or 10 (strain)
        """
        msg += self.print_header('OES1X1', 8)
        # OES1X1
        stress = [
            self.rodStress,
            self.barStress,
            self.beamStress,
            self.plateStress,
            self.solidStress,
        ]
        cstress = [self.compositePlateStress, ]  # OES1C

        # approach_code=1, table_code=1
        for data in stress:
            for isubcase in data:
                msg += self.write_markers([self.iTable, 1, 0])
                msg += self.writeOES(isubcase, data)
                self.iTable -= 1
        msg += self.write_markers([self.iTable, 1, 0])

        # ----------------
        msg += self.print_header('OSTR1X', 8)
        # OSTR1X
        strain = [
            self.rodStrain,
            self.barStrain,
            self.beamStrain,
            self.plateStrain,
            self.solidStrain,
            self.compositePlateStrain,
        ]
        for data in strain:
            for isubcase in data:
                msg += self.write_markers([self.iTable, 1, 0])
                msg += self.writeOES(isubcase, data)
                self.iTable -= 1
        msg += self.write_markers([self.iTable, 1, 0])

        # ----------------
        msg = self.print_header('OSTRIC', 8)
        cstrain = [self.compositePlateStrain, ]  # OSTR1C
        for data in cstrain:
            for isubcase in data:
                msg += self.write_markers([self.iTable, 1, 0])
                msg += self.writeOES(isubcase, data)
                self.iTable -= 1
        msg += self.write_markers([self.iTable, 1, 0])
