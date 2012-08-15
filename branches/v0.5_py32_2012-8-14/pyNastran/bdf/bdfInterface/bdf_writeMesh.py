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
# pylint: disable=C0103,C0111,W0612,R0912,R0914,R0904,W0613,E1101



from pyNastran.bdf.fieldWriter import printCard

class WriteMesh(object):
    def __init__(self):
        pass

    def echoBDF(self, infileName):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.
        @todo maybe add the write method
        """
        self.cardsToRead = set([])
        return self.readBDF(infileName)

    def autoRejectBDF(self, infileName):
        """
        This method parses supported cards, but does not group them into
        nodes, elements, properties, etc.
        @todo maybe add the write method
        """
        self.autoReject = True
        return self.readBDF(infileName)

    def write_elements_as_CTRIA3(self):
        """
        takes the cquad4 elements and splits them
        @retval msg  string representation of the elements
        """
        eids = self.elementIDs()
        #print "eids = ",eids
        nextEID = max(eids)+1  # set the new ID
        msg = '$ELEMENTS\n'
        for eid, element in sorted(self.elements.items()):
            if element.Is('CQUAD4'):
                msg += element.writeAsCTRIA3(nextEID)
                nextEID += 1
            else:
                msg += str(element)
            ###
        ###
        return msg

    def write_DMIGs(self):
        msg = ''
        for (name, dmig) in sorted(self.dmigs.items()):
            msg += str(dmig)
        for (name, dmi) in sorted(self.dmis.items()):
            msg += str(dmi)
        for (name, dmij) in sorted(self.dmijs.items()):
            msg += str(dmij)
        for (name, dmiji) in sorted(self.dmijis.items()):
            msg += str(dmiji)
        for (name, dmik) in sorted(self.dmiks.items()):
            msg += str(dmik)
        return msg

    def write_common(self):
        """
        method to write the common outputs so none get missed...
        @param self the object pointer
        @retval msg part of the bdf
        """
        msg = ''
        msg += self.write_rigid_elements()
        msg += self.write_DMIGs()
        msg += self.write_loads()
        msg += self.write_dynamic()
        msg += self.write_aero()
        msg += self.write_aero_control()
        msg += self.write_flutter()
        msg += self.write_thermal()
        msg += self.write_thermal_materials()

        msg += self.write_constraints()
        msg += self.write_optimization()
        msg += self.write_tables()
        msg += self.write_sets()
        msg += self.write_rejects()
        msg += self.write_coords()
        return msg

    def writeBDFAsPatran(self, outFileName='fem.out.bdf', debug=False):
        """
        Writes a bdf with properties & elements interspersed like how
        Patran writes the bdf.  This takes longer than the write method
        but makes it easier to compare to a Patran-formatted bdf.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        msg  = self.write_header()
        msg += self.write_params()
        msg += self.write_nodes()

        msg += self.write_elements_properties()
        msg += self.write_materials()

        msg += self.write_common()
        msg += 'ENDDATA\n'

        fname = self.print_filename(outFileName)
        self.log.debug("***writing %s" %(fname))
        
        outfile = open(outFileName,'w',encoding='utf-8')
        outfile.write(msg)
        outfile.close()
    
    def writeBDF(self, outFileName='fem.out.bdf', debug=False):
        """
        Writes the bdf.  It groups the various sections together to make it
        easy to find cards.  This method is slightly more stable than 
        writeAsPatran due to the properties sometimes being a little funny.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        msg  = self.write_header()
        msg += self.write_params()
        msg += self.write_nodes()

        msg += self.write_elements()
        msg += self.write_properties()
        msg += self.write_materials()

        msg += self.write_common()
        msg += 'ENDDATA\n'

        fname = self.print_filename(outFileName)
        self.log.debug("***writing %s" %(fname))

        outfile = open(outFileName, 'w',encoding='utf-8')
        outfile.write(msg)
        outfile.close()

    def write_as_CTRIA3(self, outFileName='fem.out.bdf', debug=False):
        """
        Writes a series of CQUAD4s as CTRIA3s.  All other cards are echoed.
        @param self the object pointer
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        @warning not tested in a long time
        """
        msg  = self.write_header()
        msg += self.write_params()
        msg += self.write_nodes()
        msg += self.write_elements_as_CTRIA3()
        msg += self.write_properties()
        msg += self.write_materials()

        msg += self.write_common()
        msg += 'ENDDATA\n'

        fname = self.print_filename(outFileName)
        self.log.debug("***writing %s" % (fname))

        outfile = open(outFileName, 'w',encoding='utf-8')
        outfile.write(msg)
        outfile.close()

    def write_header(self):
        """
        Writes the executive and case control decks.
        @param self the object pointer
        """
        msg  = self.write_executive_control_deck()
        msg += self.write_case_control_deck()
        return msg

    def write_executive_control_deck(self):
        """
        Writes the executive control deck.
        @param self the object pointer
        """
        msg = '$EXECUTIVE CONTROL DECK\n'
        
        if self.sol == 600:
            newSol = 'SOL 600,%s' % (self.solMethod)
        else:
            newSol = 'SOL %s' % (self.sol)
        
        if self.iSolLine is not None:
            self.executive_control_lines[self.iSolLine] = newSol

        for line in self.executive_control_lines:
            msg += line+'\n'
        return msg

    def write_case_control_deck(self):
        """
        Writes the Case Control Deck.
        @param self the object pointer
        """
        msg = ''
        if self.caseControlDeck:
            msg += '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
        assert 'BEGIN BULK' in msg, msg

        return msg

    def write_params(self):
        """writes the PARAM cards"""
        msg = ''
        if self.params:
            msg += '$PARAMS\n'
        #print "self.nodes = ",self.nodes
        for (key, param) in sorted(self.params.items()):
            #print "param = ",param
            msg += str(param)
        return msg

    def write_nodes(self):
        """writes the NODE-type cards"""
        msg = []
        if self.nodes:
            msg = ['$NODES\n']
            if self.gridSet:
                msg.append(str(self.gridSet))
            for (nid, node) in sorted(self.nodes.items()):
                msg.append(str(node))

        if 0:
            self.write_nodes_associated()

        if self.spoints:
            msg.append('$SPOINTS\n')
            msg.append(str(self.spoints))
        ###
        return ''.join(msg)

    def write_nodes_associated(self):
        """
        Writes the NODE-type in associated and unassociated groups.
        @warning Sometimes crashes, probably on invalid BDFs.
        """
        msg = ''

        associatedNodes = set([])
        for (eid, element) in self.elements.items():
            print(element)
            associatedNodes = associatedNodes.union(set(element.nodeIDs()))

        allNodes = set(self.nodes.keys())
        unassociatedNodes = list(allNodes.difference(associatedNodes))
        #missingNodes = allNodes.difference(
        associatedNodes = list(associatedNodes)


        if associatedNodes:
            msg += '$ASSOCIATED NODES\n'
            if self.gridSet:
                msg += str(self.gridSet)
            for key in sorted(associatedNodes):
                msg += str(self.nodes[key])

        if unassociatedNodes:
            msg += '$UNASSOCIATED NODES\n'
            if self.gridSet and not associatedNodes:
                msg += str(self.gridSet)
            for key in sorted(unassociatedNodes):
                if key in self.nodes:
                    msg += str(self.nodes[key])
                else:
                    msg += '$ Missing NodeID=%s' % (key)
        return msg

    def write_elements(self):
        """writes the elements in a sorted order"""
        msg = []
        if self.elements:
            msg = ['$ELEMENTS\n']
            for (eid, element) in sorted(self.elements.items()):
                try:
                    msg.append(str(element))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' %(element.type,eid))
                    raise
                ###
        return ''.join(msg)

    def write_rigid_elements(self):
        """writes the rigid elements in a sorted order"""
        msg = []
        if self.rigidElements:
            msg += '$RIGID ELEMENTS\n'
            for (eid, element) in sorted(self.rigidElements.items()):
                try:
                    msg.append(str(element))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' %(element.type,eid))
                    raise
                ###
        return ''.join(msg)

    def write_properties(self):
        """writes the properties in a sorted order"""
        msg = ''
        if self.properties:
            msg += '$PROPERTIES\n'
            for (pid, prop) in sorted(self.properties.items()):
                msg += str(prop)
        return msg

    def write_elements_properties(self):
        """writes the elements and properties in and interspersed order"""
        msg = []
        missingProperties = []
        if self.properties:
            msg = ['$ELEMENTS_WITH_PROPERTIES\n']

        eidsWritten = []
        for (pid, prop) in sorted(self.properties.items()):
            eids = self.getElementIDsWithPID(pid)

            if eids:
                msg.append(str(prop))
                eids.sort()
                for eid in eids:
                    element = self.Element(eid)
                    try:
                        msg.append(str(element))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' %(element.type,eid))
                        raise
                ###
                eidsWritten += eids
            else:
                missingProperties.append(str(prop))
            ###
        ###

        eidsMissing = set(self.elements.keys()).difference(set(eidsWritten))
        if eidsMissing:
            msg.append('$ELEMENTS_WITH_NO_PROPERTIES '
                       '(PID=0 and unanalyzed properties)\n')
            for eid in sorted(eidsMissing):
                element = self.Element(eid)
                try:
                    msg.append(str(element))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' %(element.type, eid))
                    raise
                ###
            ###

        if missingProperties:
            msg.append('$UNASSOCIATED_PROPERTIES\n')
            for missingProperty in missingProperties:
                msg.append(str(missingProperty))
        return ''.join(msg)

    def write_materials(self):
        """writes the materials in a sorted order"""
        msg = ''
        if self.materials:
            msg += '$MATERIALS\n'
            for (mid, material) in sorted(self.materials.items()):
                msg += str(material)
            for (mid, material) in sorted(self.creepMaterials.items()):
                msg += str(material)
            for (mid, material) in sorted(self.materialDeps.items()):
                msg += str(material)
        return msg

    def write_thermal_materials(self):
        """writes the thermal materials in a sorted order"""
        msg = ''
        if self.thermalMaterials:
            msg += '$THERMAL MATERIALS\n'
            for (mid, material) in sorted(self.thermalMaterials.items()):
                msg += str(material)
        return msg

    def write_constraints(self):
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
                for (spcID, spcadd) in sorted(self.spcadds.items()):
                    msg += str(spcadd)
                for (spcID, spcs) in sorted(self.spcs.items()):
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
                for (mpcID, mpcadd) in sorted(self.mpcadds.items()):
                    msg += str(mpcadd)
                for (mpcID, mpcs) in sorted(self.mpcs.items()):
                    for mpc in mpcs:
                        msg += str(mpc)
        return msg

    def write_loads(self):
        """writes the load cards sorted by ID"""
        msg = ''
        if self.loads:
            msg += '$LOADS\n'
            for (key, loadcase) in sorted(self.loads.items()):
                for load in loadcase:
                    try:
                        msg += str(load)
                    except:
                        print('failed printing load...type=%s key=%s'
                            %(load.type, key))
                        raise
        return msg

    def write_optimization(self):
        """writes the optimization cards sorted by ID"""
        msg = ''
        if (self.dconstrs or self.desvars or self.ddvals or self.dresps 
          or self.dvprels or self.dvmrels or self.doptprm or self.dlinks):
            msg += '$OPTIMIZATION\n'
            for (ID, dconstr) in sorted(self.dconstrs.items()):
                msg += str(dconstr)
            for (ID, desvar) in sorted(self.desvars.items()):
                msg += str(desvar)
            for (ID, ddval) in sorted(self.ddvals.items()):
                msg += str(ddval)
            for (ID, dlink) in sorted(self.dlinks.items()):
                msg += str(dlink)
            for (ID, dresp) in sorted(self.dresps.items()):
                msg += str(dresp)
            for (ID, dvmrel) in sorted(self.dvmrels.items()):
                msg += str(dvmrel)
            for (ID, dvprel) in sorted(self.dvprels.items()):
                msg += str(dvprel)
            for (ID, equation) in sorted(self.dequations.items()):
                msg += str(equation)
            if self.doptprm is not None:
                msg += str(self.doptprm)
        return msg

    def write_tables(self):
        """writes the TABLEx cards sorted by ID"""
        msg = ''
        if self.tables:
            msg += '$TABLES\n'
            for (ID, table) in sorted(self.tables.items()):
                msg += str(table)
        if self.randomTables:
            msg += '$RANDOM TABLES\n'
            for (ID, table) in sorted(self.randomTables.items()):
                msg += str(table)
        return msg

    def write_sets(self):
        """writes the SETx cards sorted by ID"""
        msg = ''
        if (self.sets or self.setsSuper or self.asets or self.bsets or
            self.csets or self.qsets):
            msg += '$SETS\n'
            for (ID, setObj) in sorted(self.sets.items()):  # dict
                msg += str(setObj)
            for setObj in self.asets:  # list
                msg += str(setObj)
            for setObj in self.bsets:  # list
                msg += str(setObj)
            for setObj in self.csets:  # list
                msg += str(setObj)
            for setObj in self.qsets:  # list
                msg += str(setObj)
            for (ID, setObj) in sorted(self.setsSuper.items()):  # dict
                msg += str(setObj)
        return msg

    def write_dynamic(self):
        """writes the dynamic cards sorted by ID"""
        msg = ''
        if (self.dareas or self.nlparms or self.frequencies or self.methods or
            self.cMethods or self.tsteps or self.tstepnls):
            msg += '$DYNAMIC\n'
            for (ID, method) in sorted(self.methods.items()):
                msg += str(method)
            for (ID, cMethod) in sorted(self.cMethods.items()):
                msg += str(cMethod)
            for (ID, darea) in sorted(self.dareas.items()):
                msg += str(darea)
            for (ID, nlparm) in sorted(self.nlparms.items()):
                msg += str(nlparm)
            for (ID, tstep) in sorted(self.tsteps.items()):
                msg += str(tstep)
            for (ID, tstepnl) in sorted(self.tstepnls.items()):
                msg += str(tstepnl)
            for (ID, freq) in sorted(self.frequencies.items()):
                msg += str(freq)
        return msg
        
    def write_aero(self):
        """writes the aero cards"""
        msg = ''
        if (self.aero or self.aeros or self.gusts or self.caeros 
        or self.paeros or self.trims):
            msg = '$AERO\n'
            for (ID, caero) in sorted(self.caeros.items()):
                msg += str(caero)
            for (ID, paero) in sorted(self.paeros.items()):
                msg += str(paero)
            for (ID, spline) in sorted(self.splines.items()):
                msg += str(spline)
            for (ID, trim) in sorted(self.trims.items()):
                msg += str(trim)

            for (ID, aero) in sorted(self.aero.items()):
                msg += str(aero)
            for (ID, aero) in sorted(self.aeros.items()):
                msg += str(aero)

            for (ID, gust) in sorted(self.gusts.items()):
                msg += str(gust)
        return msg

    def write_aero_control(self):
        """writes the aero control surface cards"""
        msg = ''
        if (self.aefacts or self.aeparams or self.aelinks or self.aelists or
            self.aestats or self.aesurfs):
            msg = '$AERO CONTROL SURFACES\n'
            for (ID, aelinks) in sorted(self.aelinks.items()):
                for aelink in aelinks:
                    msg += str(aelink)
            for (ID, aeparam) in sorted(self.aeparams.items()):
                msg += str(aeparam)
            for (ID, aestat) in sorted(self.aestats.items()):
                msg += str(aestat)

            for (ID, aelist) in sorted(self.aelists.items()):
                msg += str(aelist)
            for (ID, aesurf) in sorted(self.aesurfs.items()):
                msg += str(aesurf)
            for (ID, aefact) in sorted(self.aefacts.items()):
                msg += str(aefact)
        return msg

    def write_flutter(self):
        """writes the flutter cards"""
        msg = []
        if (self.flfacts or self.flutters or self.mkaeros):
            msg = ['$FLUTTER\n']
            for (ID, flfact) in sorted(self.flfacts.items()):
                #if ID!=0:
                msg.append(str(flfact))
            for (ID, flutter) in sorted(self.flutters.items()):
                msg.append(str(flutter))
            for mkaero in self.mkaeros:
                msg.append(str(mkaero))
        return ''.join(msg)

    def write_thermal(self):
        """writes the thermal cards"""
        msg = ''
        # PHBDY
        
        if self.phbdys or self.convectionProperties or self.bcs:
            # self.thermalProperties or
            msg = '$THERMAL\n'

            for (key, phbdy) in sorted(self.phbdys.items()):
                msg += str(phbdy)

            #for key,prop in sorted(self.thermalProperties.iteritems()):
            #    msg += str(prop)
            for (key, prop) in sorted(self.convectionProperties.items()):
                msg += str(prop)

            # BCs
            for (key, bcs) in sorted(self.bcs.items()):
                for bc in bcs: # list
                    msg += str(bc)
                ###
            ###
        return msg

    def write_coords(self):
        """writes the coordinate cards in a sorted order"""
        msg = ''
        if len(self.coords) > 1:
            msg += '$COORDS\n'
        for (ID, coord) in sorted(self.coords.items()):
            if ID != 0:
                msg += str(coord)
        return msg

    def write_rejects(self):
        """
        writes the rejected (processed) cards and the rejected unprocessed
        cardLines
        """
        msg = ''
        if self.reject_cards:
            msg += '$REJECTS\n'
            for reject_card in self.reject_cards:
                try:
                    msg += printCard(reject_card)
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                          'cards\ncard=%s\n' %(reject_card))
                    raise

        if self.rejects:
            msg += '$REJECT_LINES\n'
        for reject_lines in self.rejects:
            if reject_lines[0][0] == ' ':
                continue
            else:
                for reject in reject_lines:
                    reject2 = reject.rstrip()
                    if reject2:
                        msg += str(reject2)+'\n'
        return msg

