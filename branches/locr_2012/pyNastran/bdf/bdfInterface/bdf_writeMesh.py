# pylint: disable=C0103,C0111,W0612,R0912,R0914,R0904,W0613,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import warnings
from pyNastran.bdf.fieldWriter import print_card


class WriteMesh(object):
    def __init__(self):
        pass

    def writeBDF(self, outFileName='fem.out.bdf', size=8, debug=False):
        """
        @see write_bdf
        @warning will be removed after v0.7 in favor of write_bdf
        """
        warnings.warn('writeBDF has been deprecated; use '
                      'write_bdf', DeprecationWarning, stacklevel=2)
        self.write_bdf(outFileName, size, debug)

    def writeBDFAsPatran(self, outFileName='fem.out.bdf', size=8, debug=False):
        """
        @see write_bdf_as_patran
        @warning will be removed after v0.7 in favor of write_bdf_as_patran
        """
        warnings.warn('writeBDFAsPatran has been deprecated; use '
                      'write_bdf_as_patran', DeprecationWarning, stacklevel=2)
        self.write_bdf_as_patran(outFileName, size, debug)

    def echoBDF(self, infileName):
        """
        @see echo_bdf
        @warning will be removed after v0.7 in favor of echo_bdf
        """
        warnings.warn('echoBDF has been deprecated; use '
                      'echo_bdf', DeprecationWarning, stacklevel=2)
        self.echo_bdf(infileName)

    def echo_bdf(self, infileName):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.
        @todo maybe add the write method
        """
        self.cardsToRead = set([])
        return self.read_bdf(infileName)

    def auto_reject_bdf(self, infileName):
        """
        This method parses supported cards, but does not group them into
        nodes, elements, properties, etc.
        @todo maybe add the write method
        """
        self._auto_reject = True
        return self.read_bdf(infileName)

    def _write_elements_as_CTRIA3(self, size):
        """
        takes the cquad4 elements and splits them
        @retval msg  string representation of the elements
        """
        eids = self.elementIDs()
        #print "eids = ",eids
        nextEID = max(eids) + 1  # set the new ID
        msg = '$ELEMENTS\n'
        for eid, element in sorted(self.elements.iteritems()):
            if element.Is('CQUAD4'):
                msg += element.writeAsCTRIA3(nextEID)
                nextEID += 1
            else:
                msg += str(element)
        return msg

    def _write_dmigs(self, size):
        """
        @param self the BDF object
        @param size: large field (16) or small field (8)
        @return: msg part of the BDF
        """
        msg = []
        for (name, dmig) in sorted(self.dmigs.iteritems()):
            msg.append(str(dmig))
        for (name, dmi) in sorted(self.dmis.iteritems()):
            msg.append(str(dmi))
        for (name, dmij) in sorted(self.dmijs.iteritems()):
            msg.append(str(dmij))
        for (name, dmiji) in sorted(self.dmijis.iteritems()):
            msg.append(str(dmiji))
        for (name, dmik) in sorted(self.dmiks.iteritems()):
            msg.append(str(dmik))
        return ''.join(msg)

    def _write_common(self, size):
        """
        method to write the common outputs so none get missed...
        @param self the BDF object
        @retval msg part of the bdf
        """
        msg = ''
        msg += self._write_rigid_elements(size)
        msg += self._write_dmigs(size)
        msg += self._write_loads(size)
        msg += self._write_dynamic(size)
        msg += self._write_aero(size)
        msg += self._write_aero_control(size)
        msg += self._write_flutter(size)
        msg += self._write_thermal(size)
        msg += self._write_thermal_materials(size)

        msg += self._write_constraints(size)
        msg += self._write_optimization(size)
        msg += self._write_tables(size)
        msg += self._write_sets(size)
        msg += self._write_rejects(size)
        msg += self._write_coords(size)
        return msg

    def write_bdf_as_patran(self, outFileName='fem.out.bdf', size=8,
                            debug=False):
        """
        Writes a bdf with properties & elements interspersed like how
        Patran writes the bdf.  This takes longer than the write method
        but makes it easier to compare to a Patran-formatted bdf.
        @param self the BDF object
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        assert size in [8, 16]
        #size = 16
        fname = self.print_filename(outFileName)
        self.log.debug("***writing %s" % fname)

        outfile = open(outFileName, 'wb')

        msg = self._write_header()
        msg += self._write_params(size)
        outfile.write(msg)

        msg = self._write_nodes(size)
        outfile.write(msg)

        msg = self._write_elements_properties(size)
        outfile.write(msg)

        msg += self._write_materials(size)
        msg += self._write_common(size)
        msg += 'ENDDATA\n'
        outfile.close()

    def write_bdf(self, out_filename='fem.out.bdf', size=8, debug=False):
        """
        Writes the bdf.  It groups the various sections together to make it
        easy to find cards.  This method is slightly more stable than
        writeAsPatran due to the properties sometimes being a little funny.
        @param self the BDF object
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        """
        assert size in [8, 16]
        #size = 16
        fname = self.print_filename(out_filename)
        self.log.debug("***writing %s" % fname)

        outfile = open(out_filename, 'wb')
        msg = self._write_header()
        msg += self._write_params(size)
        outfile.write(msg)

        msg = self._write_nodes(size)
        outfile.write(msg)

        msg = self._write_elements(size)
        outfile.write(msg)

        msg = self._write_properties(size)
        msg += self._write_materials(size)
        msg += self._write_common(size)
        msg += 'ENDDATA\n'
        outfile.write(msg)
        outfile.close()

    def write_as_CTRIA3(self, outFileName='fem.out.bdf', size=8, debug=False):
        """
        Writes a series of CQUAD4s as CTRIA3s.  All other cards are echoed.
        @param self the BDF object
        @param outFileName the name to call the output bdf
        @param debug developer debug (unused)
        @warning not tested in a long time
        """
        fname = self.print_filename(outFileName)
        self.log.debug("***writing %s" % fname)

        outfile = open(outFileName, 'wb')

        msg = self._write_header()
        msg += self._write_params(size)
        outfile.write(msg)

        msg = self._write_nodes(size)
        outfile.write(msg)
        msg = self._write_elements_as_CTRIA3(size)
        outfile.write(msg)

        msg = self._write_properties(size)
        msg += self._write_materials(size)
        msg += self._write_common(size)
        msg += 'ENDDATA\n'
        outfile.write(msg)
        outfile.close()

    def _write_header(self):
        """
        Writes the executive and case control decks.
        @param self the BDF object
        """
        msg = self._write_executive_control_deck()
        msg += self._write_case_control_deck()
        return msg

    def _write_executive_control_deck(self):
        """
        Writes the executive control deck.
        @param self the BDF object
        """
        msg = '$EXECUTIVE CONTROL DECK\n'
        if self.sol == 600:
            newSol = 'SOL 600,%s' % (self.solMethod)
        else:
            newSol = 'SOL %s' % (self.sol)

        if self.iSolLine is not None:
            self.executive_control_lines[self.iSolLine] = newSol

        for line in self.executive_control_lines:
            msg += line + '\n'
        return msg

    def _write_case_control_deck(self):
        """
        Writes the Case Control Deck.
        @param self the BDF object
        """
        msg = ''
        if self.caseControlDeck:
            msg += '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
        assert 'BEGIN BULK' in msg, msg
        return msg

    def _write_params(self, size):
        """
        Writes the PARAM cards
        @param self the BDF object
        """
        msg = []
        if self.params:
            msg = ['$PARAMS\n']
            for (key, param) in sorted(self.params.iteritems()):
                msg.append(param.print_card(size))
        return ''.join(msg)

    def _write_nodes(self, size):
        """
        writes the NODE-type cards
        @param self the BDF object
        """
        msg = []
        if self.nodes:
            msg = ['$NODES\n']
            if self.gridSet:
                msg.append(str(self.gridSet))
            for (nid, node) in sorted(self.nodes.iteritems()):
                msg.append(node.print_card(size))

        if 0:
            self._write_nodes_associated(size)

        if self.spoints:
            msg.append('$SPOINTS\n')
            msg.append(str(self.spoints))
        return ''.join(msg)

    def _write_nodes_associated(self, size):
        """
        Writes the NODE-type in associated and unassociated groups.
        @param self the BDF object
        @warning Sometimes crashes, probably on invalid BDFs.
        """
        msg = []

        associatedNodes = set([])
        for (eid, element) in self.elements.iteritems():
            print(element)
            associatedNodes = associatedNodes.union(set(element.nodeIDs()))

        allNodes = set(self.nodes.keys())
        unassociatedNodes = list(allNodes.difference(associatedNodes))
        #missingNodes = allNodes.difference(
        associatedNodes = list(associatedNodes)

        if associatedNodes:
            msg += ['$ASSOCIATED NODES\n']
            if self.gridSet:
                msg.append(str(self.gridSet))
            for key, node in sorted(associatedNodes.iteritems()):
                msg.append(node.print_card(size))

        if unassociatedNodes:
            msg.append('$UNASSOCIATED NODES\n')
            if self.gridSet and not associatedNodes:
                msg.append(str(self.gridSet))
            for key, node in sorted(unassociatedNodes.iteritems()):
                if key in self.nodes:
                    msg.append(node.print_card(size))
                else:
                    msg.append('$ Missing NodeID=%s' % key)
        return ''.join(msg)

    def _write_elements(self, size):
        """
        Writes the elements in a sorted order
        @param self the BDF object
        """
        msg = []
        if self.elements:
            msg = ['$ELEMENTS\n']
            for (eid, element) in sorted(self.elements.iteritems()):
                try:
                    msg.append(str(element))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
        return ''.join(msg)

    def _write_rigid_elements(self, size):
        """writes the rigid elements in a sorted order"""
        msg = []
        if self.rigidElements:
            msg = ['$RIGID ELEMENTS\n']
            for (eid, element) in sorted(self.rigidElements.iteritems()):
                try:
                    msg.append(str(element))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
        return ''.join(msg)

    def _write_properties(self, size):
        """writes the properties in a sorted order"""
        msg = []
        if self.properties:
            msg += ['$PROPERTIES\n']
            for (pid, prop) in sorted(self.properties.iteritems()):
                msg.append(prop.print_card(size))
        return ''.join(msg)

    def _write_elements_properties(self, size):
        """writes the elements and properties in and interspersed order"""
        msg = []
        missingProperties = []
        if self.properties:
            msg.append('$ELEMENTS_WITH_PROPERTIES\n')

        eidsWritten = []
        for (pid, prop) in sorted(self.properties.iteritems()):
            eids = self.getElementIDsWithPID(pid)

            if eids:
                msg.append(prop.print_card(size))
                eids.sort()
                for eid in eids:
                    element = self.Element(eid)
                    try:
                        msg.append(str(element))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
                eidsWritten += eids
            else:
                missingProperties.append(str(prop))

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
                          'type=%s eid=%s' % (element.type, eid))
                    raise

        if missingProperties or self.pdampt or self.pbusht or self.pelast:
            msg.append('$UNASSOCIATED_PROPERTIES\n')
            for pbusht in sorted(self.pbusht.itervalues()):
                msg.append(str(pbusht))
            for pdampt in sorted(self.pdampt.itervalues()):
                msg.append(str(pdampt))
            for pelast in sorted(self.pelast.itervalues()):
                msg.append(str(pelast))
            for missingProperty in missingProperties:
                #print("missingProperty = ",missingProperty)
                #msg.append(missingProperty.print_card(size))
                msg.append(missingProperty)
        return ''.join(msg)

    def _write_materials(self, size):
        """writes the materials in a sorted order"""
        msg = []
        if self.materials:
            msg.append('$MATERIALS\n')
            for (mid, material) in sorted(self.materials.iteritems()):
                msg.append(material.print_card(size))
            for (mid, material) in sorted(self.creepMaterials.iteritems()):
                msg.append(material.print_card(size))
            for (mid, material) in sorted(self.materialDeps.iteritems()):
                msg.append(material.print_card(size))
        return ''.join(msg)

    def _write_thermal_materials(self, size):
        """writes the thermal materials in a sorted order"""
        msg = []
        if self.thermalMaterials:
            msg.append('$THERMAL MATERIALS\n')
            for (mid, material) in sorted(self.thermalMaterials.iteritems()):
                msg.append(material.print_card(size))
        return ''.join(msg)

    def _write_constraints(self, size):
        """writes the constraint cards sorted by ID"""
        msg = []
        if self.suports:
            msg.append('$CONSTRAINTS\n')
            for suport in self.suports:
                msg.append(str(suport))

        if self.spcs or self.spcadds:
            msg.append('$SPCs\n')
            strSPC = str(self.spcObject2)
            if strSPC:
                msg.append(strSPC)
            else:
                for (spcID, spcadd) in sorted(self.spcadds.iteritems()):
                    msg.append(str(spcadd))
                for (spcID, spcs) in sorted(self.spcs.iteritems()):
                    for spc in spcs:
                        msg.append(str(spc))

        if self.mpcs or self.mpcadds:
            msg.append('$MPCs\n')
            strMPC = str(self.mpcObject2)
            if strMPC:
                msg.append(strMPC)
            else:
                for (mpcID, mpcadd) in sorted(self.mpcadds.iteritems()):
                    msg.append(str(mpcadd))
                for (mpcID, mpcs) in sorted(self.mpcs.iteritems()):
                    for mpc in mpcs:
                        msg.append(str(mpc))
        return ''.join(msg)

    def _write_loads(self, size):
        """writes the load cards sorted by ID"""
        msg = []
        if self.loads:
            msg.append('$LOADS\n')
            for (key, loadcase) in sorted(self.loads.iteritems()):
                for load in loadcase:
                    try:
                        msg.append(load.print_card(size))
                    except:
                        print('failed printing load...type=%s key=%s'
                              % (load.type, key))
                        raise
        return ''.join(msg)

    def _write_optimization(self, size):
        """writes the optimization cards sorted by ID"""
        msg = []
        if (self.dconstrs or self.desvars or self.ddvals or self.dresps
            or self.dvprels or self.dvmrels or self.doptprm or self.dlinks
            or self.ddvals):
            msg.append('$OPTIMIZATION\n')
            for (ID, dconstr) in sorted(self.dconstrs.iteritems()):
                msg.append(dconstr.print_card(size))
            for (ID, desvar) in sorted(self.desvars.iteritems()):
                msg.append(desvar.print_card(size))
            for (ID, ddval) in sorted(self.ddvals.iteritems()):
                msg.append(ddval.print_card(size))
            for (ID, dlink) in sorted(self.dlinks.iteritems()):
                msg.append(dlink.print_card(size))
            for (ID, dresp) in sorted(self.dresps.iteritems()):
                msg.append(dresp.print_card(size))
            for (ID, dvmrel) in sorted(self.dvmrels.iteritems()):
                msg.append(dvmrel.print_card(size))
            for (ID, dvprel) in sorted(self.dvprels.iteritems()):
                msg.append(dvprel.print_card(size))
            for (ID, equation) in sorted(self.dequations.iteritems()):
                msg.append(str(equation))
            if self.doptprm is not None:
                msg.append(self.doptprm.print_card(size))
        return ''.join(msg)

    def _write_tables(self, size):
        """writes the TABLEx cards sorted by ID"""
        msg = []
        if self.tables:
            msg.append('$TABLES\n')
            for (ID, table) in sorted(self.tables.iteritems()):
                msg.append(table.print_card(size))
        if self.randomTables:
            msg.append('$RANDOM TABLES\n')
            for (ID, table) in sorted(self.randomTables.iteritems()):
                msg.append(table.print_card(size))
        return ''.join(msg)

    def _write_sets(self, size):
        """writes the SETx cards sorted by ID"""
        msg = []
        if (self.sets or self.setsSuper or self.asets or self.bsets or
            self.csets or self.qsets):
            msg.append('$SETS\n')
            for (ID, setObj) in sorted(self.sets.iteritems()):  # dict
                msg.append(str(setObj))
            for setObj in self.asets:  # list
                msg.append(str(setObj))
            for setObj in self.bsets:  # list
                msg.append(str(setObj))
            for setObj in self.csets:  # list
                msg.append(str(setObj))
            for setObj in self.qsets:  # list
                msg.append(str(setObj))
            for (ID, setObj) in sorted(self.setsSuper.iteritems()):  # dict
                msg.append(str(setObj))
        return ''.join(msg)

    def _write_dynamic(self, size):
        """writes the dynamic cards sorted by ID"""
        msg = []
        if (self.dareas or self.nlparms or self.frequencies or self.methods or
            self.cMethods or self.tsteps or self.tstepnls):
            msg.append('$DYNAMIC\n')
            for (ID, method) in sorted(self.methods.iteritems()):
                msg.append(method.print_card(size))
            for (ID, cMethod) in sorted(self.cMethods.iteritems()):
                msg.append(cMethod.print_card(size))
            for (ID, darea) in sorted(self.dareas.iteritems()):
                msg.append(darea.print_card(size))
            for (ID, nlparm) in sorted(self.nlparms.iteritems()):
                msg.append(nlparm.print_card(size))
            for (ID, tstep) in sorted(self.tsteps.iteritems()):
                msg.append(tstep.print_card(size))
            for (ID, tstepnl) in sorted(self.tstepnls.iteritems()):
                msg.append(tstepnl.print_card(size))
            for (ID, freq) in sorted(self.frequencies.iteritems()):
                msg.append(freq.print_card(size))
        return ''.join(msg)

    def _write_aero(self, size):
        """writes the aero cards"""
        msg = []
        if (self.aero or self.aeros or self.gusts or self.caeros
        or self.paeros or self.trims):
            msg.append('$AERO\n')
            for (ID, caero) in sorted(self.caeros.iteritems()):
                msg.append(caero.print_card(size))
            for (ID, paero) in sorted(self.paeros.iteritems()):
                msg.append(paero.print_card(size))
            for (ID, spline) in sorted(self.splines.iteritems()):
                msg.append(spline.print_card(size))
            for (ID, trim) in sorted(self.trims.iteritems()):
                msg.append(trim.print_card(size))

            for (ID, aero) in sorted(self.aero.iteritems()):
                msg.append(aero.print_card(size))
            for (ID, aero) in sorted(self.aeros.iteritems()):
                msg.append(aero.print_card(size))

            for (ID, gust) in sorted(self.gusts.iteritems()):
                msg.append(gust.print_card(size))
        return ''.join(msg)

    def _write_aero_control(self, size):
        """writes the aero control surface cards"""
        msg = []
        if (self.aefacts or self.aeparams or self.aelinks or self.aelists or
            self.aestats or self.aesurfs):
            msg.append('$AERO CONTROL SURFACES\n')
            for (ID, aelinks) in sorted(self.aelinks.iteritems()):
                for aelink in aelinks:
                    msg.append(aelink.print_card(size))
            for (ID, aeparam) in sorted(self.aeparams.iteritems()):
                msg.append(aeparam.print_card(size))
            for (ID, aestat) in sorted(self.aestats.iteritems()):
                msg.append(aestat.print_card(size))

            for (ID, aelist) in sorted(self.aelists.iteritems()):
                msg.append(aelist.print_card(size))
            for (ID, aesurf) in sorted(self.aesurfs.iteritems()):
                msg.append(aesurf.print_card(size))
            for (ID, aefact) in sorted(self.aefacts.iteritems()):
                msg.append(aefact.print_card(size))
        return ''.join(msg)

    def _write_flutter(self, size):
        """writes the flutter cards"""
        msg = []
        if (self.flfacts or self.flutters or self.mkaeros):
            msg.append('$FLUTTER\n')
            for (ID, flfact) in sorted(self.flfacts.iteritems()):
                #if ID!=0:
                msg.append(flfact.print_card(size))
            for (ID, flutter) in sorted(self.flutters.iteritems()):
                msg.append(flutter.print_card(size))
            for mkaero in self.mkaeros:
                msg.append(mkaero.print_card(size))
        return ''.join(msg)

    def _write_thermal(self, size):
        """writes the thermal cards"""
        msg = []
        # PHBDY
        if self.phbdys or self.convectionProperties or self.bcs:
            # self.thermalProperties or
            msg.append('$THERMAL\n')

            for (key, phbdy) in sorted(self.phbdys.iteritems()):
                msg.append(phbdy.print_card(size))

            #for key,prop in sorted(self.thermalProperties.iteritems()):
            #    msg.append(str(prop))
            for (key, prop) in sorted(self.convectionProperties.iteritems()):
                msg.append(prop.print_card(size))

            # BCs
            for (key, bcs) in sorted(self.bcs.iteritems()):
                for bc in bcs:  # list
                    msg.append(bc.print_card(size))
        return ''.join(msg)

    def _write_coords(self, size):
        """writes the coordinate cards in a sorted order"""
        msg = []
        if len(self.coords) > 1:
            msg.append('$COORDS\n')
        for (ID, coord) in sorted(self.coords.iteritems()):
            if ID != 0:
                msg.append(coord.print_card(size))
        return ''.join(msg)

    def _write_rejects(self, size):
        """
        writes the rejected (processed) cards and the rejected unprocessed
        cardLines
        """
        msg = []
        if self.reject_cards:
            msg.append('$REJECTS\n')
            for reject_card in self.reject_cards:
                try:
                    msg.append(print_card(reject_card))
                except RuntimeError:
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                          'cards\ncard=%s\n' % reject_card)
                    raise

        if self.rejects:
            msg.append('$REJECT_LINES\n')
        for reject_lines in self.rejects:
            if reject_lines[0][0] == ' ':
                continue
            else:
                for reject in reject_lines:
                    reject2 = reject.rstrip()
                    if reject2:
                        msg.append(str(reject2) + '\n')
        return ''.join(msg)