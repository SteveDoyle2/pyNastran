import os
from pyNastran.bdf.fieldWriter import printCard

class writeMesh(object):
    def __init__(self):
        pass

    def echoBDF(self,infileName):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.
        @todo maybe add the write method
        """
        self.cardsToRead = set([])
        return self.read(infileName)

    def autoRejectBDF(self,infileName):
        """
        This method parses supported cards, but does not group them into nodes,
        elements, properties, etc.
        @todo maybe add the write method
        """
        self.autoReject = True
        return self.read(infileName)

    def writeElementsAsCTRIA3(self):
        """
        takes the cquad4 elements and splits them
        @retval msg  string representation of the elements
        """
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
        """
        Writes a bdf with properties & elements interspersed like how
        Patran writes the bdf.  This takes longer than the write method
        but makes it easier to compare to a Patran-formatted bdf.
        """
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()

        msg += self.writeElementsProperties()

        msg += self.writeMaterials()
        msg += self.writeLoads()
        msg += self.writeAero()
        msg += self.writeThermal()
        msg += self.writeConstraints()
        msg += self.writeRejects()
        msg += self.writeCoords()
        msg += 'ENDDATA\n'

        self.log().info("***writing %s" %(os.path.relpath(outfilename)))
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

    def write(self,outfilename='fem.out.bdf',debug=False):
        """
        Writes the bdf.  It groups the various sections together to make it
        easy to find cards.
        """
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()

        msg += self.writeElements()
        msg += self.writeProperties()

        msg += self.writeMaterials()
        msg += self.writeLoads()
        msg += self.writeAero()
        msg += self.writeThermal()
        msg += self.writeConstraints()
        msg += self.writeRejects()
        msg += self.writeCoords()
        msg += 'ENDDATA\n'

        self.log().info("***writing %s" %(outfilename))
        outfile = open(outfilename,'wb')
        outfile.write(msg)
        outfile.close()

    def writeAsCTRIA3(self,outfilename='fem.out.bdf',debug=False):
        """
        Writes a series of CQUAD4s as CTRIA3s.  All other cards are echoed.
        """
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
        """
        Writes the executive and case control decks.
        """
        msg  = self.writeExecutiveControlDeck()
        msg += self.writeCaseControlDeck()
        return msg

    def writeExecutiveControlDeck(self):
        msg = '$EXECUTIVE CONTROL DECK\n'
        
        if self.sol==600:
            newSol = 'SOL 600,%s\n' %(self.solMethod)
        else:
            newSol = 'SOL %s\n' %(self.sol)
        self.executiveControlLines[self.iSolLine] = newSol

        for line in self.executiveControlLines:
            msg += line
        return msg

    def writeCaseControlDeck(self):
        msg = ''
        if self.caseControlDeck:
            msg += '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
        #for line in self.caseControlLines:
        #    msg += line

        #if 'BEGIN BULK' not in msg:
        #    msg += 'BEGIN BULK\n'
        return msg

    def writeParams(self):
        """writes the PARAM cards"""
        msg = ''
        if self.params:
            msg += '$PARAMS\n'
        #print "self.nodes = ",self.nodes
        for key,param in sorted(self.params.items()):
            #print "param = ",param
            msg += str(param)
        return msg

    def writeNodes(self):
        """writes the NODE-type cards"""
        msg = ''
        if self.nodes:
            msg += '$NODES\n'
            #print "nNodes = ",len(self.nodes)
        #print "self.nodes = ",self.nodes
        if self.gridSet:
            msg += str(self.gridSet)
        for key,node in sorted(self.nodes.items()):
            #print "node = ",node
            msg += str(node)
        return msg

    def writeElements(self):
        """writes the elements in a sorted order"""
        msg = ''
        if self.elements:
            msg += '$ELEMENTS\n'
        for eid,element in sorted(self.elements.items()):
            msg += str(element)
        return msg

    def writeProperties(self):
        """writes the properties in a sorted order"""
        msg = ''
        if self.properties:
            msg += '$PROPERTIES\n'
        for pid,prop in sorted(self.properties.items()):
            msg += str(prop)
        return msg

    def writeElementsProperties(self):
        """
        writes the elements and properties in and interspersed order
        """
        msg = ''
        missingProperties = []
        if self.properties:
            msg += '$ELEMENTS_WITH_PROPERTIES\n'
        for pid,prop in sorted(self.properties.items()):
            #print "pid = ",pid
            eids = self.getElementIDsWithPID(pid)
            #print "   eids = ",eids
            if eids:
                msg += str(prop)
                for eid in eids:
                    element = self.Element(eid)
                    #print "e.type = ",element.type
                    msg += str(element)
                ###
            else:
                #print "*MISSING",prop
                missingProperties.append(str(prop))
            ###
        ###

        eids = self.getElementIDsWithPID(0)
        if eids:
            msg += '$ELEMENTS_WITH_NO_PROPERTIES (PID=0)\n'
            for eid in sorted(eids):
                element = self.Element(eid)
                msg += str(element)
            ###

        if missingProperties:
            msg += '$UNASSOCIATED_PROPERTIES\n'
            for missingProperty in missingProperties:
                #print 'missingProperty = ',missingProperty
                #print missingProperty
                msg += str(missingProperty)
        return msg

    def writeMaterials(self):
        """writes the materials in a sorted order"""
        msg = ''
        if self.materials:
            msg += '$MATERIALS\n'
        for key,material in sorted(self.materials.items()):
            msg += str(material)
        return msg

    def writeConstraints(self):
        """writes the constraint cards sorted by ID"""
        msg = ''
        #msg += '$ where are my constraints...\n'
        if self.constraints:
            msg += '$CONSTRAINTS\n'
        for key,loadcase in sorted(self.constraints.items()):
            for constraint in loadcase:
                msg += str(constraint)

        if self.spcObject:
            msg += '$SPCs\n'
            msg += str(self.spcObject)

        if self.mpcObject:
            msg += '$MPCs\n'
            msg += str(self.mpcObject)
        return msg

    def writeLoads(self):
        """writes the load cards sorted by ID"""
        msg = ''
        if self.loads:
            msg += '$LOADS\n'
        for key,loadcase in sorted(self.loads.items()):
            for load in loadcase:
                msg += str(load)
        return msg

    def writeAero(self):
        """writes the aero cards"""
        #print "output aero cards..."
        msg = ''
        if self.flfacts or self.aeros or self.gusts or self.flutters:  msg = '$AERO\n'
        flfactKeys = self.flfacts.keys()
        #self.log().info("flfactKeys = %s" %(flfactKeys))
        for ID,flfact in sorted(self.flfacts.items()):
            #if ID!=0:
            msg += str(flfact)
        for ID,spline in sorted(self.splines.items()):
            msg += str(spline)
        for ID,caero in sorted(self.caeros.items()):
            msg += str(caero)
        for ID,aero in sorted(self.aeros.items()):
            msg += str(aero)
        for ID,gust in sorted(self.gusts.items()):
            msg += str(gust)
        for ID,flutter in sorted(self.flutters.items()):
            msg += str(flutter)

        return msg

    def writeThermal(self):
        """writes the thermal cards"""
        msg = ''
        # PHBDY
        if self.phbdys or self.convectionProperties or self.bcs:  # self.thermalProperties or
            msg = '$THERMAL\n'

        for key,phbdy in sorted(self.phbdys.items()):
            msg += str(phbdy)

        #for key,prop in sorted(self.thermalProperties.items()):
        #    msg += str(prop)
        for key,prop in sorted(self.convectionProperties.items()):
            msg += str(prop)

        # BCs
        for key,bcs in sorted(self.bcs.items()):
            for bc in bcs: # list
                msg += str(bc)
            ###
        ###
        return msg

    def writeCoords(self):
        """writes the coordinate cards in a sorted order"""
        #print "output coords..."
        msg = ''
        if len(self.coords)>1:
            msg += '$COORDS\n'
        coordKeys = self.coords.keys()
        #self.log().info("coordKeys = %s" %(coordKeys))
        for ID,coord in sorted(self.coords.items()):
            if ID!=0:
                msg += str(coord)
        return msg

    def writeRejects(self):
        """
        writes the rejected (processed) cards and the rejected unprocessed
        cardLines
        """
        msg = ''
        if self.rejectCards:
            msg += '$REJECTS\n'
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

