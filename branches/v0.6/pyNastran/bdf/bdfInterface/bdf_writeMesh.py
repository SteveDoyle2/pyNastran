# pylint: disable=C0103,C0111,W0612,R0912,R0914,R0904,W0613,E1101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import warnings
from pyNastran.bdf.fieldWriter import print_card, print_card_8, print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
#from pyNastran.bdf.bdfInterface.bdf_Reader import print_filename


class WriteMeshDeprecated(object):
    def writeBDF(self, outFileName='fem.out.bdf', size=8, debug=False):
        """
        .. seealso:: write_bdf
        .. deprecated:: will be replaced in version 0.7 with write_bdf with interspersed=False
        """
        warnings.warn('writeBDF has been deprecated; use '
                      'write_bdf', DeprecationWarning, stacklevel=2)
        self.write_bdf(out_filename=outFileName, interspersed=False, size=size, debug=debug)

    def writeBDFAsPatran(self, outFileName='fem.out.bdf', size=8, debug=False):
        """
        .. seealso:: write_bdf
        .. deprecated:: will be replaced in version 0.7 with write_bdf with an interspersed=True
        """
        warnings.warn('writeBDFAsPatran has been deprecated; use '
                      'write_bdf_as_patran', DeprecationWarning, stacklevel=2)
        self.write_bdf(outFileName, interspersed=True, size=size, debug=debug)

    def echoBDF(self, infileName):
        """
        .. seealso:: echo_bdf
        .. deprecated:: will be replaced in version 0.7 with echo_bdf
        """
        warnings.warn('echoBDF has been deprecated; use '
                      'echo_bdf', DeprecationWarning, stacklevel=2)
        self.echo_bdf(infileName)


class WriteMesh(WriteMeshDeprecated):
    def __init__(self):
        pass

    def echo_bdf(self, infile_name):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.

        .. todo:: maybe add the write method
        """
        self.cardsToRead = set([])
        return self.read_bdf(infile_name)

    def auto_reject_bdf(self, infile_name):
        """
        This method parses supported cards, but does not group them into
        nodes, elements, properties, etc.

        .. todo:: maybe add the write method
        """
        self._auto_reject = True
        return self.read_bdf(infile_name)

    def _write_elements_as_CTRIA3(self, size):
        """
        Takes the cquad4 elements and splits them

        :returns msg:  string representation of the elements
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
                msg += element.print_card(size)
        return msg

    def _write_dmigs(self, size, card_writer):
        """
        :param self:  the BDF object
        :param size:  large field (16) or small field (8)
        :returns msg: string representation of the DMIGs
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

    def _write_common(self, size, card_writer):
        """
        method to write the common outputs so none get missed...
        :param self: the BDF object
        :returns msg: part of the bdf
        """
        msg = ''
        msg += self._write_rigid_elements(size, card_writer)
        msg += self._write_dmigs(size, card_writer)
        msg += self._write_loads(size, card_writer)
        msg += self._write_dynamic(size, card_writer)
        msg += self._write_aero(size, card_writer)
        msg += self._write_aero_control(size, card_writer)
        msg += self._write_flutter(size, card_writer)
        msg += self._write_thermal(size, card_writer)
        msg += self._write_thermal_materials(size, card_writer)

        msg += self._write_constraints(size, card_writer)
        msg += self._write_optimization(size, card_writer)
        msg += self._write_tables(size, card_writer)
        msg += self._write_sets(size, card_writer)
        msg += self._write_contact(size, card_writer)
        msg += self._write_rejects(size, card_writer)
        msg += self._write_coords(size, card_writer)
        return msg

    def write_bdf(self, out_filename='fem.out.bdf', interspersed=True,
                  size=8, precision='single', debug=False):
        """
        Writes the BDF.

        :param self:         the BDF object
        :param out_filename: the name to call the output bdf
        :param debug:        developer debug (unused)
        :param interspersed: Writes a bdf with properties & elements
              interspersed like how Patran writes the bdf.  This takes
              slightly longer than if interspersed=False, but makes it
              much easier to compare to a Patran-formatted bdf and is
              more clear. (default=True)
        :param size:  the field size (8 is recommended)
        :param precision:  'single', 'double'
        :param debug: developer debug
        """
        if size == 8:
            assert precision == 'single', 'precision=%r' % precision
            card_writer = print_card_8
        elif size == 16:
            assert precision in ['single', 'double'], 'precision=%r' % precision
            if precision == 'single':
                card_writer = print_card_16
            else:
                card_writer = print_card_double
        else:
            assert size in [8, 16]

        assert isinstance(interspersed, bool)
        #size = 16
        fname = self.print_filename(out_filename)
        self.log.debug("***writing %s" % fname)

        outfile = open(out_filename, 'wb')
        msg = self._write_header()
        msg += self._write_params(size, card_writer)
        outfile.write(msg)

        msg = self._write_nodes(size, card_writer)
        outfile.write(msg)

        if interspersed:
            msg = self._write_elements_properties(size, card_writer)
        else:
            msg = self._write_elements(size, card_writer)
            outfile.write(msg)
            msg = self._write_properties(size, card_writer)

        outfile.write(msg)

        msg = self._write_materials(size, card_writer)
        msg += self._write_common(size, card_writer)
        msg += 'ENDDATA\n'
        outfile.write(msg)
        outfile.close()

    def write_as_CTRIA3(self, out_filename='fem.out.bdf', size=8, debug=False):
        """
        Writes a series of CQUAD4s as CTRIA3s.  All other cards are echoed.
        :param self:         the BDF object
        :param out_filename: the name to call the output bdf
        :param debug:        developer debug (unused)
        .. warning:: not tested in a long time
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
        :param self: the BDF object
        """
        msg = self._write_executive_control_deck()
        msg += self._write_case_control_deck()
        return msg

    def _write_executive_control_deck(self):
        """
        Writes the executive control deck.
        :param self: the BDF object
        """
        msg = ''
        if self.executive_control_lines:
            msg = '$EXECUTIVE CONTROL DECK\n'
            if self.sol == 600:
                newSol = 'SOL 600,%s' % self.solMethod
            else:
                newSol = 'SOL %s' % self.sol

            if self.iSolLine is not None:
                self.executive_control_lines[self.iSolLine] = newSol

            for line in self.executive_control_lines:
                msg += line + '\n'
        return msg

    def _write_case_control_deck(self):
        """
        Writes the Case Control Deck.
        :param self: the BDF object
        """
        msg = ''
        if self.caseControlDeck:
            msg += '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
            assert 'BEGIN BULK' in msg, msg
        return msg

    def _write_params(self, size, card_writer):
        """
        Writes the PARAM cards
        :param self: the BDF object
        """
        msg = []
        if self.params:
            msg = ['$PARAMS\n']
            for (key, param) in sorted(self.params.iteritems()):
                #msg.append(param.print_card(size))
                msg.append(param.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_nodes(self, size, card_writer):
        """
        Writes the NODE-type cards
        :param self: the BDF object
        """
        msg = []
        if self.spoints:
            msg.append('$SPOINTS\n')
            msg.append(str(self.spoints))

        if self.nodes:
            msg.append('$NODES\n')
            if self.gridSet:
                msg.append(self.gridSet.print_card(size))
            for (nid, node) in sorted(self.nodes.iteritems()):
                #msg.append(node.print_card(size))
                msg.append(node.write_bdf(size, card_writer))
        if 0:
            self._write_nodes_associated(size)

        return ''.join(msg)

    def _write_nodes_associated(self, size, card_writer):
        """
        Writes the NODE-type in associated and unassociated groups.
        :param self: the BDF object
        .. warning:: Sometimes crashes, probably on invalid BDFs.
        """
        msg = []
        associated_nodes = set([])
        for (eid, element) in self.elements.iteritems():
            #print(element)
            associated_nodes = associated_nodes.union(set(element.nodeIDs()))

        all_nodes = set(self.nodes.keys())
        unassociated_nodes = list(all_nodes.difference(associated_nodes))
        #missing_nodes = all_nodes.difference(
        associated_nodes = list(associated_nodes)

        if associated_nodes:
            msg += ['$ASSOCIATED NODES\n']
            if self.gridSet:
                msg.append(str(self.gridSet))
            for key, node in sorted(associated_nodes.iteritems()):
                msg.append(node.print_card(size))

        if unassociated_nodes:
            msg.append('$UNASSOCIATED NODES\n')
            if self.gridSet and not associated_nodes:
                msg.append(str(self.gridSet))
            for key, node in sorted(unassociated_nodes.iteritems()):
                if key in self.nodes:
                    msg.append(node.print_card(size))
                else:
                    msg.append('$ Missing NodeID=%s' % key)
        return ''.join(msg)

    def _write_elements(self, size, card_writer):
        """
        Writes the elements in a sorted order
        :param self: the BDF object
        """
        msg = []
        if self.elements:
            msg = ['$ELEMENTS\n']
            for (eid, element) in sorted(self.elements.iteritems()):
                try:
                    msg.append(element.print_card(size))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
        return ''.join(msg)

    def _write_rigid_elements(self, size, card_writer):
        """Writes the rigid elements in a sorted order"""
        msg = []
        if self.rigidElements:
            msg = ['$RIGID ELEMENTS\n']
            for (eid, element) in sorted(self.rigidElements.iteritems()):
                try:
                    msg.append(element.write_bdf(size, card_writer))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
        return ''.join(msg)

    def _write_properties(self, size, card_writer):
        """Writes the properties in a sorted order"""
        msg = []
        if self.properties:
            msg += ['$PROPERTIES\n']
            for (pid, prop) in sorted(self.properties.iteritems()):
                msg.append(prop.print_card(size))
        return ''.join(msg)

    def _write_elements_properties(self, size, card_writer):
        """
        Writes the elements and properties in and interspersed order
        """
        msg = []
        missing_properties = []
        if self.properties:
            msg.append('$ELEMENTS_WITH_PROPERTIES\n')

        eids_written = []
        pids = sorted(self.properties.keys())
        eids2 = self.getElementIDsWithPIDs(pids, mode='dict')

        for (pid, eids) in sorted(eids2.iteritems()):
            prop = self.properties[pid]
            if eids:
                #msg.append(prop.print_card(size))
                msg.append(prop.write_bdf(size, card_writer))
                eids.sort()
                for eid in eids:
                    element = self.Element(eid)
                    try:
                        #msg.append(element.print_card(size))
                        msg.append(element.write_bdf(size, card_writer))
                    except:
                        print('failed printing element...' 'type=%r eid=%s' % (element.type, eid))
                        raise
                eids_written += eids
            else:
                missing_properties.append(str(prop))

        eids_missing = set(self.elements.keys()).difference(set(eids_written))

        if eids_missing:
            msg.append('$ELEMENTS_WITH_NO_PROPERTIES '
                       '(PID=0 and unanalyzed properties)\n')
            for eid in sorted(eids_missing):
                element = self.Element(eid)
                try:
                    msg.append(element.write_bdf(size, card_writer))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise

        if missing_properties or self.pdampt or self.pbusht or self.pelast:
            msg.append('$UNASSOCIATED_PROPERTIES\n')
            for pbusht in sorted(self.pbusht.itervalues()):
                msg.append(str(pbusht))
            for pdampt in sorted(self.pdampt.itervalues()):
                msg.append(str(pdampt))
            for pelast in sorted(self.pelast.itervalues()):
                msg.append(str(pelast))
            for missing_property in missing_properties:
                #print("missing_property = ",missing_property)
                #msg.append(missing_property.print_card(size))
                msg.append(missing_property)
        return ''.join(msg)

    def _write_materials(self, size, card_writer):
        """Writes the materials in a sorted order"""
        msg = []
        if self.materials:
            msg.append('$MATERIALS\n')
            for (mid, material) in sorted(self.materials.iteritems()):
                #msg.append(material.print_card(size))
                msg.append(material.write_bdf(size, card_writer))
            for (mid, material) in sorted(self.creepMaterials.iteritems()):
                msg.append(material.print_card(size))
            for (mid, material) in sorted(self.materialDeps.iteritems()):
                msg.append(material.print_card(size))
        return ''.join(msg)

    def _write_thermal_materials(self, size, card_writer):
        """Writes the thermal materials in a sorted order"""
        msg = []
        if self.thermalMaterials:
            msg.append('$THERMAL MATERIALS\n')
            for (mid, material) in sorted(self.thermalMaterials.iteritems()):
                msg.append(material.print_card(size))
        return ''.join(msg)

    def _write_constraints(self, size, card_writer):
        """Writes the constraint cards sorted by ID"""
        msg = []
        if self.suports:
            msg.append('$CONSTRAINTS\n')
            for suport in self.suports:
                msg.append(str(suport))

        if self.spcs or self.spcadds:
            msg.append('$SPCs\n')
            strSPC = str(self.spcObject)
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
            strMPC = str(self.mpcObject)
            if strMPC:
                msg.append(strMPC)
            else:
                for (mpcID, mpcadd) in sorted(self.mpcadds.iteritems()):
                    msg.append(str(mpcadd))
                for (mpcID, mpcs) in sorted(self.mpcs.iteritems()):
                    for mpc in mpcs:
                        msg.append(str(mpc))
        return ''.join(msg)

    def _write_loads(self, size, card_writer):
        """Writes the load cards sorted by ID"""
        msg = []
        if self.loads:
            msg.append('$LOADS\n')
            for (key, loadcase) in sorted(self.loads.iteritems()):
                for load in loadcase:
                    try:
                        #msg.append(load.print_card(size))
                        msg.append(load.write_bdf(size, card_writer))
                    except:
                        print('failed printing load...type=%s key=%r' % (load.type, key))
                        raise
        return ''.join(msg)

    def _write_contact(self, size, card_writer):
        """Writes the contact cards sorted by ID"""
        msg = []
        if (self.bcrparas or self.bctadds or self.bctparas or self.bctsets
            or self.bsurf or self.bsurfs):
            msg.append('$CONTACT\n')
            for (ID, bcrpara) in sorted(self.bcrparas.iteritems()):
                msg.append(bcrpara.print_card(size))
            for (ID, bctadds) in sorted(self.bctadds.iteritems()):
                msg.append(bctadds.print_card(size))
            for (ID, bctpara) in sorted(self.bctparas.iteritems()):
                msg.append(bctpara.print_card(size))

            for (ID, bctset) in sorted(self.bctsets.iteritems()):
                msg.append(bctset.write_bdf(size, card_writer))
            for (ID, bsurfi) in sorted(self.bsurf.iteritems()):
                msg.append(bsurfi.write_bdf(size, card_writer))
            for (ID, bsurfsi) in sorted(self.bsurfs.iteritems()):
                #msg.append(bsurfsi.print_card(size))
                msg.append(bsurfsi.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_optimization(self, size, card_writer):
        """Writes the optimization cards sorted by ID"""
        msg = []
        if (self.dconstrs or self.desvars or self.ddvals or self.dresps
            or self.dvprels or self.dvmrels or self.doptprm or self.dlinks
            or self.ddvals):
            msg.append('$OPTIMIZATION\n')
            for (ID, dconstr) in sorted(self.dconstrs.iteritems()):
                msg.append(dconstr.write_bdf(size, card_writer))
            for (ID, desvar) in sorted(self.desvars.iteritems()):
                msg.append(desvar.write_bdf(size, card_writer))
            for (ID, ddval) in sorted(self.ddvals.iteritems()):
                msg.append(ddval.write_bdf(size, card_writer))
            for (ID, dlink) in sorted(self.dlinks.iteritems()):
                msg.append(dlink.write_bdf(size, card_writer))
            for (ID, dresp) in sorted(self.dresps.iteritems()):
                msg.append(dresp.write_bdf(size, card_writer))
            for (ID, dvmrel) in sorted(self.dvmrels.iteritems()):
                msg.append(dvmrel.write_bdf(size, card_writer))
            for (ID, dvprel) in sorted(self.dvprels.iteritems()):
                msg.append(dvprel.write_bdf(size, card_writer))
            for (ID, equation) in sorted(self.dequations.iteritems()):
                msg.append(str(equation))
            if self.doptprm is not None:
                msg.append(self.doptprm.write_bdf(size, card_writer))

        return ''.join(msg)

    def _write_tables(self, size, card_writer):
        """Writes the TABLEx cards sorted by ID"""
        msg = []
        if self.tables:
            msg.append('$TABLES\n')
            for (ID, table) in sorted(self.tables.iteritems()):
                msg.append(table.write_bdf(size, card_writer))

        if self.randomTables:
            msg.append('$RANDOM TABLES\n')
            for (ID, table) in sorted(self.randomTables.iteritems()):
                msg.append(table.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_sets(self, size, card_writer):
        """Writes the SETx cards sorted by ID"""
        msg = []
        if (self.sets or self.setsSuper or self.asets or self.bsets or
            self.csets or self.qsets):
            msg.append('$SETS\n')
            for (ID, setObj) in sorted(self.sets.iteritems()):  # dict
                msg.append(setObj.write_bdf(size, card_writer))
            for setObj in self.asets:  # list
                msg.append(setObj.write_bdf(size, card_writer))
            for setObj in self.bsets:  # list
                msg.append(setObj.write_bdf(size, card_writer))
            for setObj in self.csets:  # list
                msg.append(setObj.write_bdf(size, card_writer))
            for setObj in self.qsets:  # list
                msg.append(setObj.write_bdf(size, card_writer))
            for (ID, setObj) in sorted(self.setsSuper.iteritems()):  # dict
                msg.append(setObj.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_dynamic(self, size, card_writer):
        """Writes the dynamic cards sorted by ID"""
        msg = []
        if (self.dareas or self.nlparms or self.frequencies or self.methods or
            self.cMethods or self.tsteps or self.tstepnls):
            msg.append('$DYNAMIC\n')
            for (ID, method) in sorted(self.methods.iteritems()):
                #msg.append(method.print_card(size))
                msg.append(method.write_bdf(size, card_writer))
            for (ID, cMethod) in sorted(self.cMethods.iteritems()):
                #msg.append(cMethod.print_card(size))
                msg.append(cMethod.write_bdf(size, card_writer))
            for (ID, darea) in sorted(self.dareas.iteritems()):
                #msg.append(darea.print_card(size))
                msg.append(darea.write_bdf(size, card_writer))
            for (ID, nlparm) in sorted(self.nlparms.iteritems()):
                #msg.append(nlparm.print_card(size))
                msg.append(nlparm.write_bdf(size, card_writer))
            for (ID, nlpci) in sorted(self.nlpcis.iteritems()):
                msg.append(nlpci.print_card(size))
            for (ID, tstep) in sorted(self.tsteps.iteritems()):
                msg.append(tstep.print_card(size))
            for (ID, tstepnl) in sorted(self.tstepnls.iteritems()):
                msg.append(tstepnl.print_card(size))
            for (ID, freq) in sorted(self.frequencies.iteritems()):
                msg.append(freq.print_card(size))
        return ''.join(msg)

    def _write_aero(self, size, card_writer):
        """Writes the aero cards"""
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
                #msg.append(trim.print_card(size))
                msg.append(trim.write_bdf(size, card_writer))

            for (ID, aero) in sorted(self.aero.iteritems()):
                msg.append(aero.print_card(size))
            for (ID, aero) in sorted(self.aeros.iteritems()):
                msg.append(aero.print_card(size))

            for (ID, gust) in sorted(self.gusts.iteritems()):
                msg.append(gust.print_card(size))
        return ''.join(msg)

    def _write_aero_control(self, size, card_writer):
        """Writes the aero control surface cards"""
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
                msg.append(aelist.write_bdf(size, card_writer))
            for (ID, aesurf) in sorted(self.aesurfs.iteritems()):
                msg.append(aesurf.print_card(size))
            for (ID, aefact) in sorted(self.aefacts.iteritems()):
                msg.append(aefact.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_flutter(self, size, card_writer):
        """Writes the flutter cards"""
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

    def _write_thermal(self, size, card_writer):
        """Writes the thermal cards"""
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

    def _write_coords(self, size, card_writer):
        """Writes the coordinate cards in a sorted order"""
        msg = []
        if len(self.coords) > 1:
            msg.append('$COORDS\n')
        for (ID, coord) in sorted(self.coords.iteritems()):
            if ID != 0:
                msg.append(coord.write_bdf(size, card_writer))
        return ''.join(msg)

    def _write_rejects(self, size, card_writer):
        """
        Writes the rejected (processed) cards and the rejected unprocessed
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