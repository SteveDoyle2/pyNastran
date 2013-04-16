from struct import pack


class Ougv1Writer(object):
    def writeOUGV1(self):
        msg = ''
        self.device_code = 1  # print the OP2...

        self.write_string_block('OUGV1', 8)

        msg += self.write_markers([-1, 7])
        out = [101, 0, 4136, 0, 0, 0, 1]  # TODO what this is - DMAP -> "no def or month,year,one,one"...huh???
        msg += pack('iiiiiii', *out)

        msg += self.write_markers([-2, 1, 0])

        # approach_code=1, table_code=1
        self.iTable = -3
        data = self.displacements
        for isubcase in data:
            msg += self.write_markers([self.iTable, 1, 0])
            msg += self.writeOUG_displacements(isubcase, data)
            self.iTable -= 1

        data = self.temperatures
        for isubcase in data:
            msg += self.write_markers([self.iTable, 1, 0])
            msg += self.writeOUG_displacements(isubcase, data, thermal=1)
            self.iTable -= 1

        # approach_code=3, table_code=1
        #for isubcase in self.fluxes:
        #    msg += write_markers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(isubcase,iTable)
        #    iTable-=1

        # approach_code=6, table_code=1
        #data = self.nonlinearDisplacements
        #for isubcase in self.nonlinearDisplacements:
        #    msg += write_markers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(isubcase,iTable)
        #    iTable-=1

        #data = self.nonlinearTemperatures
        #for isubcase in self.nonlinearTemperatures:
        #    msg += write_markers([iTable,1,0])
        #    msg += self.writeOUG_temperatures(isubcase,iTable)
        #    iTable-=1

        msg += self.write_markers([self.iTable - 1, 1, 0])
        #self.displacements = {}
        #self.temperatures  = {}

    def writeOUG_displacements(self, isubcase, data, thermal=0):
        """
        this function writes table 3 for the OUGV1 Table
        @todo add the buffer and block caps
        """
        msg = ''
        disp = data[isubcase]
        self.write_markers([146, 584])

        lsdvmn = isubcase  # TODO is this correct???

        if disp.dt is None:
            approach_code = 1  # statics
            five = lsdvmn
            FiveSixSeven = [lsdvmn, 0, 0]  # fields five,six,seven
        else:
            approach_code = 6  # transient
            FiveSixSeven = [disp.dt, 0, 0]  # fields five,six,seven

        table_code = 1  # statics
        sort_code = 0
        randomCode = 0  # 8 @todo no idea...
        format_code = 1  # 9 - Real numbers
        num_wide = 7  # 10
        #thermal = 0 # 23
        (aCode, tCode) = self.aCode_tCode(approach_code, table_code, sort_code)

        #12345
        zero = pack('i', 0)
        msg += pack('iiii', aCode, tCode, 0, isubcase)  # 1,2,3,4
        msg += pack('iii', *FiveSixSeven)  # 5,6,7
        msg += pack('iii', randomCode, format_code, num_wide)  # 8,9,10

        #22-11 = 11
        msg += zero * 11

        msg += pack('i', thermal)

        #51-23 = 28
        msg += zero * 28
        msg += self.packTitle(isubcase)

        msg += data[isubcase].writeOp2('', self.device_code)
        #msg += pack('i',4) # data length...
        return msg
