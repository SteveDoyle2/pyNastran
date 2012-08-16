from struct import pack


class Oes1Writer(object):
    def writeOES1(self):
        """
        writes isotropic/composite stress/strain
        @todo assumes sCode=0 (stress) or 10 (strain)
        """
        msg += self.printHeader('OES1X1', 8)
        # OES1X1
        stress = [
            self.rodStress,
            self.barStress,
            self.beamStress,
            self.plateStress,
            self.solidStress,
        ]
        cstress = [self.compositePlateStress, ]  # OES1C

        # approachCode=1, tableCode=1
        for data in stress:
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable, 1, 0])
                msg += self.writeOES(iSubcase, data)
                self.iTable -= 1
        msg += self.writeMarkers([self.iTable, 1, 0])

        # ----------------
        msg += self.printHeader('OSTR1X', 8)
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
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable, 1, 0])
                msg += self.writeOES(iSubcase, data)
                self.iTable -= 1
        msg += self.writeMarkers([self.iTable, 1, 0])

        # ----------------
        msg = self.printHeader('OSTRIC', 8)
        cstrain = [self.compositePlateStrain, ]  # OSTR1C
        for data in cstrain:
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable, 1, 0])
                msg += self.writeOES(iSubcase, data)
                self.iTable -= 1
        msg += self.writeMarkers([self.iTable, 1, 0])
