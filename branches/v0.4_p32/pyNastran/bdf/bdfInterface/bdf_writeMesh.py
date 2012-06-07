## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os
from pyNastran.bdf.fieldWriter import printCard
from pyNastran.bdf.errors import *

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
        return self.readBDF(infileName)

    def autoRejectBDF(self,infileName):
        """
        This method parses supported cards, but does not group them into nodes,
        elements, properties, etc.
        @todo maybe add the write method
        """
        self.autoReject = True
        return self.readBDF(infileName)

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

    def writeCommon(self):
        """
        method to write the common outputs so none get missed...
        @param self the object pointer
        @retval msg part of the bdf
        """
        msg = ''
        msg += self.writeRigidElements()
        msg += self.writeLoads()
        msg += self.writeDynamic()
        msg += self.writeAero()
        msg += self.writeAeroControl()
        msg += self.writeFlutter()
        msg += self.writeThermal()
        msg += self.writeThermalMaterials()

        msg += self.writeConstraints()
        msg += self.writeOptimization()
        msg += self.writeTables()
        msg += self.writeSets()
        msg += self.writeRejects()
        msg += self.writeCoords()
        return msg

    def writeBDFAsPatran(self,outFileName='fem.out.bdf',debug=False):
        """
        Writes a bdf with properties & elements interspersed like how
        Patran writes the bdf.  This takes longer than the write method
        but makes it easier to compare to a Patran-formatted bdf.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()

        msg += self.writeElementsProperties()
        msg += self.writeMaterials()

        msg += self.writeCommon()
        msg += 'ENDDATA\n'

        fname = self.printFileName(outFileName)
        self.log.debug("***writing %s" %(fname))
        
        outfile = open(outFileName,'w',encoding='ascii')
        outfile.write(msg)
        outfile.close()
    
    def writeBDF(self,outFileName='fem.out.bdf',debug=False):
        """
        Writes the bdf.  It groups the various sections together to make it
        easy to find cards.  This method is slightly more stable than 
        writeAsPatran due to the properties sometimes being a little funny.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()

        msg += self.writeElements()
        msg += self.writeProperties()
        msg += self.writeMaterials()

        msg += self.writeCommon()
        msg += 'ENDDATA\n'

        fname = self.printFileName(outFileName)
        self.log.debug("***writing %s" %(fname))

        outfile = open(outFileName,'w',encoding='ascii')
        outfile.write(msg)
        outfile.close()

    def writeAsCTRIA3(self,outFileName='fem.out.bdf',debug=False):
        """
        Writes a series of CQUAD4s as CTRIA3s.  All other cards are echoed.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        @warning not tested in a long time
        """
        msg  = self.writeHeader()
        msg += self.writeParams()
        msg += self.writeNodes()
        msg += self.writeElementsAsCTRIA3()
        msg += self.writeProperties()
        msg += self.writeMaterials()

        msg += self.writeCommon()
        msg += 'ENDDATA\n'

        fname = self.printFileName(outFileName)
        self.log.debug("***writing %s" %(fname))

        outfile = open(outFileName,'w',encoding='ascii')
        outfile.write(msg)
        outfile.close()

    def writeHeader(self):
        """
        Writes the executive and case control decks.
        @param self the object pointer
        """
        msg  = self.writeExecutiveControlDeck()
        msg += self.writeCaseControlDeck()
        return msg

    def writeExecutiveControlDeck(self):
        """
        Writes the executive control deck.
        @param self the object pointer
        """
        msg = '$EXECUTIVE CONTROL DECK\n'
        
        if self.sol==600:
            newSol = 'SOL 600,%s' %(self.solMethod)
        else:
            newSol = 'SOL %s' %(self.sol)
        
        if self.iSolLine is not None:
            self.executiveControlLines[self.iSolLine] = newSol

        for line in self.executiveControlLines:
            msg += line+'\n'
        return msg

    def writeCaseControlDeck(self):
        """
        Writes the case control deck.
        @param self the object pointer
        """
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

        if self.spoints:
            msg += '$SPOINTS\n'
            msg += str(self.spoints)
        ###
        return msg

    def writeElements(self):
        """writes the elements in a sorted order"""
        msg = ''
        if self.elements:
            msg += '$ELEMENTS\n'
            for eid,element in sorted(self.elements.items()):
                try:
                    msg += str(element)
                except:
                    print('failed printing element...type=%s eid=%s' %(element.type,eid))
                    raise
                ###
        return msg

    def writeRigidElements(self):
        """writes the rigid elements in a sorted order"""
        msg = ''
        if self.rigidElements:
            msg += '$RIGID ELEMENTS\n'
            for eid,element in sorted(self.rigidElements.items()):
                try:
                    msg += str(element)
                except:
                    print('failed printing element...type=%s eid=%s' %(element.type,eid))
                    raise
                ###
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
        """writes the elements and properties in and interspersed order"""
        msg = ''
        missingProperties = []
        if self.properties:
            msg += '$ELEMENTS_WITH_PROPERTIES\n'

        eidsWritten = []
        for pid,prop in sorted(self.properties.items()):
            #print "pid = ",pid
            eids = self.getElementIDsWithPID(pid)
            eids.sort()
            #print "   eids = ",eids
            if eids:
                msg += str(prop)
                for eid in eids:
                    element = self.Element(eid)
                    #print "e.type = ",element.type
                    try:
                        msg += str(element)
                    except:
                        print('failed printing element...type=%s eid=%s' %(element.type,eid))
                        raise
                    ###
                ###
                eidsWritten+=eids
            else:
                #print "*MISSING",prop
                missingProperties.append(str(prop))
            ###
        ###

        eidsMissing = set(self.elements.keys()).difference(set(eidsWritten))
        if eidsMissing:
            msg += '$ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)\n'
            for eid in sorted(eidsMissing):
                element = self.Element(eid)
                try:
                    msg += str(element)
                except:
                    print('failed printing element...type=%s eid=%s' %(element.type,eid))
                    raise
                ###
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
            for key,material in sorted(self.creepMaterials.items()):
                msg += str(material)
            for key,material in sorted(self.materialDeps.items()):
                msg += str(material)
        return msg

    def writeThermalMaterials(self):
        """writes the thermal materials in a sorted order"""
        msg = ''
        if self.thermalMaterials:
            msg += '$THERMAL MATERIALS\n'
            for key,material in sorted(self.thermalMaterials.items()):
                msg += str(material)
        return msg

    def writeConstraints(self):
        """writes the constraint cards sorted by ID"""
        msg = ''
        if self.suports:
            msg += '$CONSTRAINTS\n'
            for suport in self.suports:
                msg += str(suport)
            ###

        if self.spcs or self.spcadds:
            msg += '$SPCs\n'
            strSPC = str(self.spcObject2)
            if strSPC:
                msg += strSPC
            else:
                for spcID, spcadd in sorted(self.spcadds.items()):
                    msg += str(spcadd)
                for spcID, spcs in sorted(self.spcs.items()):
                    for spc in spcs:
                        msg += str(spc)
                    ###
                ###
            ###
        ###
        
        if self.mpcs or self.mpcadds:
            msg += '$MPCs\n'
            strMPC = str(self.mpcObject2)
            if strMPC:
                msg += strMPC
            else:
                for mpcID,mpcadd in sorted(self.mpcadds.items()):
                    msg += str(mpcadd)
                for mpcID,mpcs in sorted(self.mpcs.items()):
                    for mpc in mpcs:
                        msg += str(mpc)
                    ###
                ###
            ###
        ###
        return msg
        #msg += '$ where are my constraints...\n'
        if self.constraints or self.suports:
            msg += '$CONSTRAINTS\n'
            for key,loadcase in sorted(self.constraints.items()):
                for constraint in loadcase:
                    msg += str(constraint)
            for suport in self.suports:
                msg += str(suport)
            ###
        ###

    def writeLoads(self):
        """writes the load cards sorted by ID"""
        msg = ''
        if self.loads or self.gravs:
            msg += '$LOADS\n'
            for key,loadcase in sorted(self.loads.items()):
                for load in loadcase:
                    try:
                        msg += str(load)
                    except:
                        print('failed printing load...type=%s key=%s' %(load.type,key))
                        raise
                    ###
            for ID,grav in sorted(self.gravs.items()):
                msg += str(grav)
            ###
        ###
        return msg

    def writeOptimization(self):
        msg = ''
        if (self.dconstrs or self.desvars or self.ddvals or self.dresps 
          or self.dvprels or self.dvmrels or self.doptprm or self.dlinks):
            msg += '$OPTIMIZATION\n'
            for ID,dconstr in sorted(self.dconstrs.items()):
                msg += str(dconstr)
            for ID,desvar in sorted(self.desvars.items()):
                msg += str(desvar)
            for ID,ddval in sorted(self.ddvals.items()):
                msg += str(ddval)
            for ID,dlink in sorted(self.dlinks.items()):
                msg += str(dlink)
            for ID,dresp in sorted(self.dresps.items()):
                msg += str(dresp)
            for ID,dvmrel in sorted(self.dvmrels.items()):
                msg += str(dvmrel)
            for ID,dvprel in sorted(self.dvprels.items()):
                msg += str(dvprel)
            for ID,equation in sorted(self.dequations.items()):
                msg += str(equation)
            if self.doptprm is not None:
                msg += str(self.doptprm)
        return msg

    def writeTables(self):
        msg = ''
        if self.tables:
            msg += '$TABLES\n'
            for ID,table in sorted(self.tables.items()):
                msg += str(table)
        return msg

    def writeSets(self):
        msg = ''
        if self.sets or self.setsSuper:
            msg += '$SETS\n'
            for ID,setObj in sorted(self.sets.items()):
                msg += str(setObj)
            for ID,setObj in sorted(self.setsSuper.items()):
                msg += str(setObj)
        return msg

    def writeDynamic(self):
        msg = ''
        if self.dareas or self.nlparms or self.frequencies or self.methods or self.tsteps:
            msg += '$DYNAMIC\n'
            for ID,method in sorted(self.methods.items()):
                msg += str(method)
            for ID,darea in sorted(self.dareas.items()):
                msg += str(darea)
            for ID,nlparm in sorted(self.nlparms.items()):
                msg += str(nlparm)
            for ID,tstep in sorted(self.tsteps.items()):
                msg += str(tstep)
            for ID,freq in sorted(self.frequencies.items()):
                msg += str(freq)
        return msg
        
    def writeAero(self):
        """writes the aero cards"""
        #print "output aero cards..."
        msg = ''
        if (self.aero or self.aeros or self.gusts or self.caeros 
        or self.paeros or self.trims):
            msg = '$AERO\n'
            for ID,caero in sorted(self.caeros.items()):
                msg += str(caero)
            for ID,paero in sorted(self.paeros.items()):
                msg += str(paero)
            for ID,spline in sorted(self.splines.items()):
                msg += str(spline)
            for ID,trim in sorted(self.trims.items()):
                msg += str(trim)

            for ID,aero in sorted(self.aero.items()):
                msg += str(aero)
            for ID,aero in sorted(self.aeros.items()):
                msg += str(aero)

            for ID,gust in sorted(self.gusts.items()):
                msg += str(gust)
        return msg

    def writeAeroControl(self):
        """writes the aero control surface cards"""
        msg = ''
        if (self.aefacts or self.aeparams or self.aelinks or self.aelists or
            self.aestats or self.aesurfs):
            msg = '$AERO CONTROL SURFACES\n'
            for ID,aelinks in sorted(self.aelinks.items()):
                for aelink in aelinks:
                    msg += str(aelink)
            for ID,aeparam in sorted(self.aeparams.items()):
                msg += str(aeparam)
            for ID,aestat in sorted(self.aestats.items()):
                msg += str(aestat)

            for ID,aelist in sorted(self.aelists.items()):
                msg += str(aelist)
            for ID,aesurf in sorted(self.aesurfs.items()):
                msg += str(aesurf)
            for ID,aefact in sorted(self.aefacts.items()):
                msg += str(aefact)
        return msg

    def writeFlutter(self):
        """writes the flutter cards"""
        msg = ''
        if (self.flfacts or self.flutters or self.mkaeros):
            msg = '$FLUTTER\n'
            for ID,flfact in sorted(self.flfacts.items()):
                #if ID!=0:
                msg += str(flfact)
            for ID,flutter in sorted(self.flutters.items()):
                msg += str(flutter)
            for mkaero in self.mkaeros:
                msg += str(mkaero)
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
        #self.log.info("coordKeys = %s" %(coordKeys))
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
                try:
                    msg += printCard(rejectCard)
                except:
                    
                    for field in rejectCard:
                        if field is not None and '=' in field:
                            msg = 'cannot reject equal signed cards\ncard=%s\n' %(rejectCard)
                            raise BDF_SyntaxError(msg)
                        ###
                    ###
                    raise

        if self.rejects:
            msg += '$REJECT_LINES\n'
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

