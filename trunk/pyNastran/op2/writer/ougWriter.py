from struct import pack


class Ougv1Writer(object):
    def writeOUGV1(self):
        msg = ''
        self.deviceCode = 1  # print the OP2...

        self.writeStringBlock('OUGV1', 8)

        msg += self.writeMarkers([-1, 7])
        out = [101, 0, 4136, 0, 0, 0, 1]  ## @todo what this is - DMAP -> "no def or month,year,one,one"...huh???
        msg += pack('iiiiiii', *out)

        msg += self.writeMarkers([-2, 1, 0])

        # approachCode=1, tableCode=1
        self.iTable = -3
        data = self.displacements
        for iSubcase in data:
            msg += self.writeMarkers([self.iTable, 1, 0])
            msg += self.writeOUG_displacements(iSubcase, data)
            self.iTable -= 1

        data = self.temperatures
        for iSubcase in data:
            msg += self.writeMarkers([self.iTable, 1, 0])
            msg += self.writeOUG_displacements(iSubcase, data, thermal=1)
            self.iTable -= 1

        # approachCode=3, tableCode=1
        #for iSubcase in self.fluxes:
        #    msg += writeMarkers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(iSubcase,iTable)
        #    iTable-=1

        # approachCode=6, tableCode=1
        #data = self.nonlinearDisplacements
        #for iSubcase in self.nonlinearDisplacements:
        #    msg += writeMarkers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(iSubcase,iTable)
        #    iTable-=1

        #data = self.nonlinearTemperatures
        #for iSubcase in self.nonlinearTemperatures:
        #    msg += writeMarkers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(iSubcase,iTable)
        #    iTable-=1

        msg += self.writeMarkers([self.iTable - 1, 1, 0])
        #self.displacements = {}
        #self.temperatures  = {}

    def writeOUG_displacements(self, iSubcase, data, thermal=0):
        """
        this function writes table 3 for the OUGV1 Table
        @todo add the buffer and block caps
        """
        msg = ''
        disp = data[iSubcase]
        self.writeMarkers([146, 584])

        lsdvmn = iSubcase  ## @todo is this correct???

        if disp.dt is None:
            approachCode = 1  # statics
            five = lsdvmn
            FiveSixSeven = [lsdvmn, 0, 0]  # fields five,six,seven
        else:
            approachCode = 6  # transient
            FiveSixSeven = [disp.dt, 0, 0]  # fields five,six,seven
        ###

        tableCode = 1  # statics
        sortCode = 0
        randomCode = 0  # 8 @todo no idea...
        formatCode = 1  # 9 - Real numbers
        numWide = 7  # 10
        #thermal = 0 # 23
        (aCode, tCode) = self.aCode_tCode(approachCode, tableCode, sortCode)

        #12345
        zero = pack('i', 0)
        msg += pack('iiii', aCode, tCode, 0, iSubcase)  # 1,2,3,4
        msg += pack('iii', *FiveSixSeven)  # 5,6,7
        msg += pack('iii', randomCode, formatCode, numWide)  # 8,9,10

        #22-11 = 11
        msg += zero * 11

        msg += pack('i', thermal)

        #51-23 = 28
        msg += zero * 28
        msg += self.packTitle(iSubcase)

        msg += data[iSubcase].writeOp2('', self.deviceCode)
        #msg += pack('i',4) # data length...
        return msg
