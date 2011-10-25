from pyNastran.bdf.fieldWriter import printCard

class writeMesh(object):
    def echoBDF(self,infileName):
        self.cardsToRead = set([])
        return self.read(infileName)

    def writeElementsAsCTRIA3(self):
        eids = self.elementIDs()
        #print "eids = ",eids
        nextEID = max(eids)+1  # set the new ID
        msg = '$ELEMENTS\n'
        for key,element in sorted(self.elements.items()):
            if element.Is('CQUAD4'):
                msg += element.writeAsCTRIA3(nextEID)
                nextEID+=1
            else:
                msg += str(element)
            ###
        ###
        return msg

    def writeAsPatran(self,outfilename='fem.out.bdf',debug=False):
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()
        msg += self.writeElementsProperties()
        msg += self.writeMaterials()
        msg += self.writeLoads()
        msg += self.writeAero()
        msg += self.writeConstraints()
        msg += self.writeRejects()
        msg += self.writeCoords()
        msg += 'ENDDATA\n'

        self.log().info("***writing %s" %(outfilename))
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

    def write(self,outfilename='fem.out.bdf',debug=False):
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()

        msg += self.writeElements()
        msg += self.writeProperties()

        msg += self.writeMaterials()
        msg += self.writeLoads()
        msg += self.writeAero()
        msg += self.writeConstraints()
        msg += self.writeRejects()
        msg += self.writeCoords()
        msg += 'ENDDATA\n'

        self.log().info("***writing %s" %(outfilename))
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

    def writeAsCTRIA3(self,outfilename='fem.out.bdf',debug=False):
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()
        msg += self.writeElementsAsCTRIA3()
        msg += self.writeProperties()
        msg += self.writeMaterials()
        msg += self.writeLoads()
        msg += self.writeConstraints()
        msg += self.writeRejects()
        msg += self.writeCoords()
        msg += 'ENDDATA\n'

        self.log().info("***writing %s" %(outfilename))
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

    def writeHeader(self):
        msg = '$EXECUTIVE CONTROL DECK\n'
        for line in self.executiveControlLines:
            msg += line

        if self.caseControlDeck:
            msg += '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
        #for line in self.caseControlLines:
        #    msg += line
        return msg

    def writeParams(self):
        msg = ''
        if self.params:
            msg += '$PARAMS\n'
        #print "self.nodes = ",self.nodes
        for key,param in sorted(self.params.items()):
            #print "param = ",param
            msg += str(param)
        return msg

    def writeNodes(self):
        msg = ''
        if self.nodes:
            msg += '$NODES\n'
            #print "nNodes = ",len(self.nodes)
        #print "self.nodes = ",self.nodes
        for key,node in sorted(self.nodes.items()):
            #print "node = ",node
            msg += str(node)
        return msg

    def writeElements(self):
        msg = ''
        if self.elements:
            msg += '$ELEMENTS\n'
        for eid,element in sorted(self.elements.items()):
            msg += str(element)
        return msg

    def writeProperties(self):
        msg = ''
        if self.properties:
            msg += '$PROPERTIES\n'
        for pid,prop in sorted(self.properties.items()):
            msg += str(prop)
        return msg

    def writeElementsProperties(self):
        msg = ''
        if self.properties:
            msg += '$ELEMENTS_WITH_PROPERTIES\n'
        for pid,prop in sorted(self.properties.items()):
            msg += str(prop)
            eids = self.getElementIDsWithPID(pid)

            for eid in eids:
                element = self.Element(eid)
                msg += str(element)
            ###
        ###
        msg += '$ELEMENTS_WITH_NO_PROPERTIES\n'
        eids = self.getElementIDsWithPID(0)
        for eid in eids:
            element = self.Element(eid)
            msg += str(element)
        ###
        return msg

    def writeMaterials(self):
        msg = ''
        if self.materials:
            msg += '$MATERIALS\n'
        for key,material in sorted(self.materials.items()):
            msg += str(material)
        return msg

    def writeConstraints(self):
        msg = ''
        if self.constraints:
            msg += '$CONSTRAINTS\n'
        for key,loadcase in sorted(self.constraints.items()):
            for constraint in loadcase:
                msg += str(constraint)
        return msg

    def writeLoads(self):
        msg = ''
        if self.loads:
            msg += '$LOADS\n'
        for key,loadcase in sorted(self.loads.items()):
            for load in loadcase:
                msg += str(load)
        return msg

    def writeAero(self):
        #print "output aero cards..."
        msg = ''
        if self.flfacts:  msg = '$AERO\n'
        flfactKeys = self.flfacts.keys()
        #self.log().info("flfactKeys = %s" %(flfactKeys))
        for ID,flfact in sorted(self.flfacts.items()):
            #if ID!=0:
            msg += str(flfact)
        for ID,aero in sorted(self.aeros.items()):
            msg += str(aero)
        for ID,gust in sorted(self.gusts.items()):
            msg += str(gust)
        return msg

    def writeCoords(self):
        #print "output coords..."
        msg = ''
        #if self.coords:
        msg += '$COORDS\n'
        coordKeys = self.coords.keys()
        self.log().info("coordKeys = %s" %(coordKeys))
        for ID,coord in sorted(self.coords.items()):
            if ID!=0:
                msg += str(coord)
        return msg

    def writeRejects(self):
        msg = '$REJECTS\n'
        for rejectCard in self.rejectCards:
            #print "rejectCard = ",rejectCard
            #print ""
            msg += printCard(rejectCard)
        for rejectLines in self.rejects:
            if rejectLines[0][0]==' ':
                continue
            else:
                #print "rejectLines = ",rejectLines
                for reject in rejectLines:
                    reject2 = reject.rstrip()
                    if reject2:
                        msg += str(reject2)+'\n'
                    ###
                ###
            ###
        ###
        return msg

