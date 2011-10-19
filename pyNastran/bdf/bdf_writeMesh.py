class writeMesh(object):
    def writeHeader(self):
        msg = '$EXECUTIVE CONTROL DECK\n'
        for line in self.executiveControlLines:
            msg += line

        msg += '$CASE CONTROL DECK\n'
        msg += str(self.caseControlDeck)
        #for line in self.caseControlLines:
        #    msg += line
        return msg

    def writeParams(self):
        msg = '$PARAMS\n'
        #print "self.nodes = ",self.nodes
        for key,param in sorted(self.params.items()):
            #print "param = ",param
            msg += str(param)
        return msg

    def writeNodes(self):
        msg = '$NODES\n'
        #print "self.nodes = ",self.nodes
        print "nNodes = ",len(self.nodes)
        for key,node in sorted(self.nodes.items()):
            #print "node = ",node
            msg += str(node)
        return msg

    def writeElements(self):
        msg = '$ELEMENTS\n'
        for key,element in sorted(self.elements.items()):
            msg += str(element)
        return msg

    def writeProperties(self):
        msg = '$PROPERTIES\n'
        for key,prop in sorted(self.properties.items()):
            msg += str(prop)
        return msg

    def writeMaterials(self):
        msg = '$MATERIALS\n'
        for key,material in sorted(self.materials.items()):
            msg += str(material)
        return msg

    def writeConstraints(self):
        msg = '$CONSTRAINTS\n'
        for key,loadcase in sorted(self.constraints.items()):
            for constraint in loadcase:
                msg += str(constraint)
        return msg

    def writeLoads(self):
        msg = '$LOADS\n'
        for key,loadcase in sorted(self.loads.items()):
            for load in loadcase:
                msg += str(load)
        return msg

    def writeCoords(self):
        #print "output coords..."
        msg = '$COORDS\n'
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

