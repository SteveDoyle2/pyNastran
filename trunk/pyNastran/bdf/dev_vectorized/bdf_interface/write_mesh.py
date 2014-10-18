"""
This file defines:
  - WriteMeshDeprecated (NA for major releases)
  - WriteMesh
"""
from pyNastran.bdf.utils import print_filename
from numpy import array, unique, concatenate, intersect1d, where

class WriteMesh(object):
    """
    Major methods:
      - model.write_bdf(...)
      - model.echo_bdf(...)
      - model.auto_reject_bdf(...)
    """
    def __init__(self):
        pass

    def echo_bdf(self, infile_name):
        """
        This method removes all comment lines from the bdf
        A write method is stil required.

        .. todo:: maybe add the write method
        """
        self.cardsToRead = set([])
        pass

    def auto_reject_bdf(self, infile_name):
        """
        This method parses supported cards, but does not group them into
        nodes, elements, properties, etc.

        .. todo:: maybe add the write method
        """
        self._auto_reject = True
        pass

    def _output_helper(self, out_filename, interspersed, size, precision):
        """
        Performs type checking on the write_bdf inputs
        """
        if out_filename is None:
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Save BDF/DAT/PCH'
            out_filename = save_file_dialog(title, wildcard_wx, wildcard_qt)
            assert out_filename is not None, out_filename
        if not isinstance(out_filename, basestring):
            raise TypeError('out_filename=%r must be a string' % out_filename)

        if size == 8:
            assert precision == 'single', 'precision=%r' % precision
        elif size == 16:
            assert precision in ['single', 'double'], 'precision=%r' % precision
        else:
            assert size in [8, 16], size

        assert isinstance(interspersed, bool)
        fname = print_filename(out_filename, self._relpath)
        self.log.debug("***writing %s" % fname)
        return out_filename

    def write_bdf(self, out_filename=None, interspersed=True,
                  size=8, precision='single', enddata=None):
        """
        Writes the BDF.

        :param self:         the BDF object
        :param out_filename: the name to call the output bdf
                             (default=None; pops a dialog)
        :param interspersed: Writes a bdf with properties & elements
              interspersed like how Patran writes the bdf.  This takes
              slightly longer than if interspersed=False, but makes it
              much easier to compare to a Patran-formatted bdf and is
              more clear. (default=True)
        :param size:  the field size (8 is recommended)
        :param precision:  'single', 'double'
        :param enddata:  Flag to enable/disable writing ENDDATA
                         (default=None -> depends on input BDF)
        """
        out_filename = self._output_helper(out_filename,
                                           interspersed, size, precision)

        if precision == 'double':
            double = True
        elif precision == 'single':
            double = False
        is_double = False

        outfile = open(out_filename, 'wb')
        self._write_header(outfile)
        self._write_params(outfile, size, is_double)

        #==========================
        self._write_nodes(outfile, size, is_double)

        self.write_elements_properties(outfile, size, is_double, interspersed)
        self.materials.write_bdf(outfile, size, is_double)
        self._write_common(outfile, size, is_double)

        #self._write_constraints(outfile, size)
        self._write_nonlinear(outfile, size)
        self._write_rejects(outfile, size)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            outfile.write('ENDDATA\n')
        outfile.close()


    def _write_header(self, outfile):
        """
        Writes the executive and case control decks.

        :param self: the BDF object
        """
        self._write_executive_control_deck(outfile)
        self._write_case_control_deck(outfile)

    def _write_executive_control_deck(self, outfile):
        """
        Writes the executive control deck.

        :param self: the BDF object
        """
        msg = ''
        if self.executive_control_lines:
            msg = '$EXECUTIVE CONTROL DECK\n'
            if self.sol == 600:
                new_sol = 'SOL 600,%s' % self.solMethod
            else:
                new_sol = 'SOL %s' % self.sol

            if self.iSolLine is not None:
                self.executive_control_lines[self.iSolLine] = new_sol

            for line in self.executive_control_lines:
                msg += line + '\n'
            outfile.write(''.join(msg))

    def _write_case_control_deck(self, outfile):
        """
        Writes the Case Control Deck.

        :param self: the BDF object
        """
        if self.caseControlDeck:
            msg = '$CASE CONTROL DECK\n'
            msg += str(self.caseControlDeck)
            assert 'BEGIN BULK' in msg, msg
            outfile.write(''.join(msg))

    def _write_elements_properties(self, f, size):
        """
        Writes the elements and properties in and interspersed order
        """
        #return self._write_elements_properties2(f, size)
        msg = []
        missing_properties = []

        eids_written = []
        #pids = sorted(self.properties.keys())

        ptypes = [self.properties_shell.pshell,
                  self.properties_shell.pcomp,
                  self.pshear,
                  self.prod,
                  self.properties_solid.psolid,

                  #self.properties_bar.pbar,
                  #self.properties_bar.pbarl,
                  #self.properties_beam.pbeam,
                  #self.properties_beam.pbeaml,
                  ]

        n = 0
        pids_all = None  # the actual properties
        for t in ptypes:
            if t.n and n == 0:
                pids_all = t.property_id
                n = 1
            elif t.n:
                print pids_all
                print t.property_id
                try:
                    pids_all = concatenate(pids_all, t.property_id)
                except ValueError:
                    pids_all = array(list(pids_all) + list(t.property_id))

        etypes = (self.elements_shell._get_types() +
                  self.elements_solid._get_types() +
                  [self.crod, self.cshear])

        #pids_set = None
        if pids_all is None:
            f.write('$MISSING_ELEMENTS because there are no properties\n')
            for t in etypes:
                #print "t.type =", t.type
                t.write_bdf(f, size=size)
            return

        # there are properties
        pids_set = set(list(pids_all))

        n = 0
        pids = None
        for t in etypes:
            #print "t.type =", t.type
            if t.n and n == 0:
                eids = t.element_id
                pids = t.property_id
                n = 1
            elif t.n:
                try:
                    eids = concatenate(eids, t.element_id)
                #except AttributeError:
                    #eids = array(list(eids) + list(t.element_id))
                except TypeError:
                    #print eids
                    #print t.element_id
                    eids = array(list(eids) + list(t.element_id))
                except ValueError:
                    #print eids
                    #print t.element_id
                    eids = array(list(eids) + list(t.element_id))
                    #asdf

                try:
                    pids = concatenate(pids, t.property_id)
                except AttributeError:
                    pids = array(list(pids) + list(t.property_id))
                except TypeError:
                    pids = array(list(pids) + list(t.property_id))
                except ValueError:
                    pids = array(list(pids) + list(t.property_id))
            #else:
                #print t.type

        elements_by_pid = {}
        if pids is not None:
            pids_unique = unique(pids)
            print "pids_unique =", pids_unique
            pids_unique.sort()
            if len(pids_unique) > 0:
                f.write('$ELEMENTS_WITH_PROPERTIES\n')

            for pid in pids_all:
                i = where(pid==pids)[0]
                eids2 = eids[i]

                for t in ptypes:
                    if t.n and pid in t.property_id:
                        print("prop.type = %s" % t.type)
                        t.write_bdf(f, size=size, property_ids=[pid])
                        pids_set.remove(pid)
                n = 0
                for t in etypes:
                    if not t.n:
                        continue
                    eids3 = intersect1d(t.element_id, eids2, assume_unique=False)
                    #print "eids3[pid=%s]" %(pid), eids3
                    if n == 0 and len(eids3):
                        elements_by_pid[pid] = eids3
                        n = 1
                    elif len(eids3):
                        try:
                            c = concatenate(elements_by_pid[pid], eids3)
                        except TypeError:
                            c = array(list(elements_by_pid[pid]) + list(eids3))
                        except ValueError:
                            c = array(list(elements_by_pid[pid]) + list(eids3))
                        elements_by_pid[pid] = c
                    else:
                        continue
                    try:
                        t.write_bdf(f, size=size, element_ids=eids3)
                    except TypeError:
                        print("t.type = %s" % t.type)
                        raise
                    del eids3
            #for pid, elements in elements_by_pid.items():
                #print "pid=%s n=%s" % (pid, len(elements))
            #print elements_by_pid

        # missing properties
        if pids_set:
            pids_list = list(pids_set)
            f.write('$UNASSOCIATED_PROPERTIES\n')
            for pid in pids_list:
                for t in ptypes:
                    if t.n and pid in t.property_id:
                        t.write_bdf(f, size=size, property_ids=[pid])

        #..todo:: finish...
        f.write('$UNASSOCIATED_ELEMENTS\n')
        # missing elements...

    def write_elements_properties(self, f, size, is_double, interspersed):
        self.elements.write_bdf(f, size=size, is_double=is_double, include_properties=True, interspersed=interspersed)

    def _write_aero(self, outfile, size, card_writer):
        """Writes the aero cards"""
        if (self.aero or self.aeros or self.gusts or self.caeros
        or self.paeros or self.trims):
            msg = ['$AERO\n']
            for (unused_id, caero) in sorted(self.caeros.iteritems()):
                msg.append(caero.write_bdf(size, card_writer))
            for (unused_id, paero) in sorted(self.paeros.iteritems()):
                msg.append(paero.write_bdf(size, card_writer))
            for (unused_id, spline) in sorted(self.splines.iteritems()):
                msg.append(spline.write_bdf(size, card_writer))
            for (unused_id, trim) in sorted(self.trims.iteritems()):
                msg.append(trim.write_bdf(size, card_writer))

            for (unused_id, aero) in sorted(self.aero.iteritems()):
                msg.append(aero.write_bdf(size, card_writer))
            for (unused_id, aero) in sorted(self.aeros.iteritems()):
                msg.append(aero.write_bdf(size, card_writer))

            for (unused_id, gust) in sorted(self.gusts.iteritems()):
                msg.append(gust.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    #def _write_aero(self, f, size):
        #self.paero.write_bdf(f, size)
        #self.caero.write_bdf(f, size)
        #for key, trim in sorted(self.trim.iteritems()):
            #trim.write_bdf(f, size)
        ##self.aero.write_bdf(f, size)
        ##self.aeros.write_bdf(f, size)

    def _write_aero_control(self, outfile, size, card_writer):
        """Writes the aero control surface cards"""
        if (self.aefacts or self.aeparams or self.aelinks or self.aelists or
            self.aestats or self.aesurfs):
            msg = ['$AERO CONTROL SURFACES\n']
            for (unused_id, aelinks) in sorted(self.aelinks.iteritems()):
                for aelink in aelinks:
                    msg.append(aelink.write_bdf(size, card_writer))
            for (unused_id, aeparam) in sorted(self.aeparams.iteritems()):
                msg.append(aeparam.write_bdf(size, card_writer))
            for (unused_id, aestat) in sorted(self.aestats.iteritems()):
                msg.append(aestat.write_bdf(size, card_writer))

            for (unused_id, aelist) in sorted(self.aelists.iteritems()):
                msg.append(aelist.write_bdf(size, card_writer))
            for (unused_id, aesurf) in sorted(self.aesurfs.iteritems()):
                msg.append(aesurf.write_bdf(size, card_writer))
            for (unused_id, aefact) in sorted(self.aefacts.iteritems()):
                msg.append(aefact.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_common(self, outfile, size, card_writer):
        """
        Write the common outputs so none get missed...

        :param self: the BDF object
        :returns msg: part of the bdf
        """
        is_double = False
        self._write_rigid_elements(outfile, size, card_writer)
        self._write_dmigs(outfile, size, card_writer)
        self.loads.write_bdf(outfile, size, is_double, sort_data=False)
        self._write_dynamic(outfile, size, card_writer)
        self._write_aero(outfile, size, card_writer)
        self._write_aero_control(outfile, size, card_writer)
        self._write_flutter(outfile, size, card_writer)
        self._write_thermal(outfile, size, card_writer)
        self._write_thermal_materials(outfile, size, card_writer)

        self._write_constraints(outfile, size, card_writer)
        self._write_optimization(outfile, size, card_writer)
        self._write_tables(outfile, size, card_writer)
        self._write_sets(outfile, size, card_writer)
        self._write_contact(outfile, size, card_writer)
        self._write_rejects(outfile, size)
        self._write_coords(outfile, size, card_writer)

    def _write_constraints(self, f, size, card_writer):
        spcs = [self.spcadd, self.spc, self.spcd, self.spc1]
        mpcs = [self.mpcadd, self.mpc]
        self._write_constraints_spc_mpc(f, size, spcs)
        self._write_constraints_spc_mpc(f, size, mpcs)

    def _write_constraints_spc_mpc(self, f, size, types):
        interspersed = False
        if interspersed:
            raise NotImplementedError()
        else:
            ids = []
            for t in types:
                ids += t.iterkeys()
            ids = unique(ids)
            ids.sort()
            if len(ids) > 0:
                f.write('$CONSTRAINTS\n')
                for ID in ids:
                    for t in types:
                        for constraint_id, constraint in sorted(t.iteritems()):
                            if ID == constraint_id:
                                constraint.write_bdf(f, size=size)

    def _write_contact(self, outfile, size, card_writer):
        """Writes the contact cards sorted by ID"""
        if (self.bcrparas or self.bctadds or self.bctparas or self.bctsets
            or self.bsurf or self.bsurfs):
            msg = ['$CONTACT\n']
            for (unused_id, bcrpara) in sorted(self.bcrparas.iteritems()):
                msg.append(bcrpara.write_bdf(size, card_writer))
            for (unused_id, bctadds) in sorted(self.bctadds.iteritems()):
                msg.append(bctadds.write_bdf(size, card_writer))
            for (unused_id, bctpara) in sorted(self.bctparas.iteritems()):
                msg.append(bctpara.write_bdf(size, card_writer))

            for (unused_id, bctset) in sorted(self.bctsets.iteritems()):
                msg.append(bctset.write_bdf(size, card_writer))
            for (unused_id, bsurfi) in sorted(self.bsurf.iteritems()):
                msg.append(bsurfi.write_bdf(size, card_writer))
            for (unused_id, bsurfsi) in sorted(self.bsurfs.iteritems()):
                msg.append(bsurfsi.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_coords(self, outfile, size, card_writer):
        """Writes the coordinate cards in a sorted order"""
        msg = []
        if len(self.coords) > 1:
            msg.append('$COORDS\n')
        for (unused_id, coord) in sorted(self.coords.iteritems()):
            if unused_id != 0:
                msg.append(coord.write_bdf(size, card_writer))
        outfile.write(''.join(msg))

    def _write_dmigs(self, outfile, size, card_writer):
        """
        Writes the DMIG cards

        :param self:  the BDF object
        :param size:  large field (16) or small field (8)
        :returns msg: string representation of the DMIGs
        """
        msg = []
        for (unused_name, dmig) in sorted(self.dmigs.iteritems()):
            msg.append(dmig.write_bdf(size, card_writer))
        for (unused_name, dmi) in sorted(self.dmis.iteritems()):
            msg.append(dmi.write_bdf(size, card_writer))
        for (unused_name, dmij) in sorted(self.dmijs.iteritems()):
            msg.append(dmij.write_bdf(size, card_writer))
        for (unused_name, dmiji) in sorted(self.dmijis.iteritems()):
            msg.append(dmiji.write_bdf(size, card_writer))
        for (unused_name, dmik) in sorted(self.dmiks.iteritems()):
            msg.append(dmik.write_bdf(size, card_writer))
        outfile.write(''.join(msg))

    def _write_dynamic(self, outfile, size, card_writer):
        """Writes the dynamic cards sorted by ID"""
        if (self.dareas or self.nlparms or self.frequencies or self.methods or
            self.cMethods or self.tsteps or self.tstepnls):
            msg = ['$DYNAMIC\n']
            for (unused_id, method) in sorted(self.methods.iteritems()):
                msg.append(method.write_bdf(size, card_writer))
            for (unused_id, cMethod) in sorted(self.cMethods.iteritems()):
                msg.append(cMethod.write_bdf(size, card_writer))
            for (unused_id, darea) in sorted(self.dareas.iteritems()):
                msg.append(darea.write_bdf(size, card_writer))
            for (unused_id, nlparm) in sorted(self.nlparms.iteritems()):
                msg.append(nlparm.write_bdf(size, card_writer))
            for (unused_id, nlpci) in sorted(self.nlpcis.iteritems()):
                msg.append(nlpci.write_bdf(size, card_writer))
            for (unused_id, tstep) in sorted(self.tsteps.iteritems()):
                msg.append(tstep.write_bdf(size, card_writer))
            for (unused_id, tstepnl) in sorted(self.tstepnls.iteritems()):
                msg.append(tstepnl.write_bdf(size, card_writer))
            for (unused_id, freq) in sorted(self.frequencies.iteritems()):
                msg.append(freq.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_flutter(self, outfile, size, card_writer):
        """Writes the flutter cards"""
        if self.flfacts or self.flutters or self.mkaeros:
            msg = ['$FLUTTER\n']
            for (unused_id, flfact) in sorted(self.flfacts.iteritems()):
                #if unused_id != 0:
                msg.append(flfact.write_bdf(size, card_writer))
            for (unused_id, flutter) in sorted(self.flutters.iteritems()):
                msg.append(flutter.write_bdf(size, card_writer))
            for mkaero in self.mkaeros:
                msg.append(mkaero.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_masses(self, outfile, size, card_writer):
        pass

    def _write_nodes(self, f, size, is_double):
        self.grdset.write_bdf(f, size)
        self.grid.write_bdf(f, size)
        self.point.write_bdf(f, size)
        self.epoint.write_bdf(f, size)
        self.spoint.write_bdf(f, size)
        self.pointax.write_bdf(f, size)

    def _write_nonlinear(self, f, size):
        for key, card in sorted(self.nlparm.iteritems()):
            card.write_bdf(f, size)
        for key, card in sorted(self.nlpci.iteritems()):
            card.write_bdf(f, size)
        #self.tables1.write_bdf(f, size)

    def _write_optimization(self, outfile, size, card_writer):
        """Writes the optimization cards sorted by ID"""
        if (self.dconstrs or self.desvars or self.ddvals or self.dresps
            or self.dvprels or self.dvmrels or self.doptprm or self.dlinks
            or self.ddvals):
            msg = ['$OPTIMIZATION\n']
            for (unused_id, dconstr) in sorted(self.dconstrs.iteritems()):
                msg.append(dconstr.write_bdf(size, card_writer))
            for (unused_id, desvar) in sorted(self.desvars.iteritems()):
                msg.append(desvar.write_bdf(size, card_writer))
            for (unused_id, ddval) in sorted(self.ddvals.iteritems()):
                msg.append(ddval.write_bdf(size, card_writer))
            for (unused_id, dlink) in sorted(self.dlinks.iteritems()):
                msg.append(dlink.write_bdf(size, card_writer))
            for (unused_id, dresp) in sorted(self.dresps.iteritems()):
                msg.append(dresp.write_bdf(size, card_writer))
            for (unused_id, dvmrel) in sorted(self.dvmrels.iteritems()):
                msg.append(dvmrel.write_bdf(size, card_writer))
            for (unused_id, dvprel) in sorted(self.dvprels.iteritems()):
                msg.append(dvprel.write_bdf(size, card_writer))
            for (unused_id, equation) in sorted(self.dequations.iteritems()):
                msg.append(str(equation))
            if self.doptprm is not None:
                msg.append(self.doptprm.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_params(self, outfile, size, is_double):
        """
        Writes the PARAM cards

        :param self: the BDF object
        """
        if self.params:
            msg = ['$PARAMS\n']
            for (unused_key, param) in sorted(self.params.iteritems()):
                msg.append(param.write_bdf(size))
            outfile.write(''.join(msg))

    def _write_rigid_elements(self, outfile, size, card_writer):
        """Writes the rigid elements in a sorted order"""
        if self.rigidElements:
            msg = ['$RIGID ELEMENTS\n']
            for (eid, element) in sorted(self.rigidElements.iteritems()):
                try:
                    msg.append(element.write_bdf(size, card_writer))
                except:
                    print('failed printing element...'
                          'type=%s eid=%s' % (element.type, eid))
                    raise
            outfile.write(''.join(msg))

    def _write_rejects(self, outfile, size):
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
        outfile.write(''.join(msg))

    def _write_sets(self, outfile, size, card_writer):
        """Writes the SETx cards sorted by ID"""
        if (self.sets or self.setsSuper or self.asets or self.bsets or
            self.csets or self.qsets):
            msg = ['$SETS\n']
            for (unused_id, set_obj) in sorted(self.sets.iteritems()):  # dict
                msg.append(set_obj.write_bdf(size, card_writer))
            for set_obj in self.asets:  # list
                msg.append(set_obj.write_bdf(size, card_writer))
            for set_obj in self.bsets:  # list
                msg.append(set_obj.write_bdf(size, card_writer))
            for set_obj in self.csets:  # list
                msg.append(set_obj.write_bdf(size, card_writer))
            for set_obj in self.qsets:  # list
                msg.append(set_obj.write_bdf(size, card_writer))
            for (set_id, set_obj) in sorted(self.setsSuper.iteritems()):  # dict
                msg.append(set_obj.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_tables(self, outfile, size, card_writer):
        """Writes the TABLEx cards sorted by ID"""
        if self.tables:
            msg = ['$TABLES\n']
            for (unused_id, table) in sorted(self.tables.iteritems()):
                msg.append(table.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

        if self.randomTables:
            msg = ['$RANDOM TABLES\n']
            for (unused_id, table) in sorted(self.randomTables.iteritems()):
                msg.append(table.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_thermal(self, outfile, size, card_writer):
        """Writes the thermal cards"""
        # PHBDY
        if self.phbdys or self.convectionProperties or self.bcs:
            # self.thermalProperties or
            msg = ['$THERMAL\n']

            for (unused_key, phbdy) in sorted(self.phbdys.iteritems()):
                msg.append(phbdy.write_bdf(size, card_writer))

            #for unused_key, prop in sorted(self.thermalProperties.iteritems()):
            #    msg.append(str(prop))
            for (unused_key, prop) in sorted(self.convectionProperties.iteritems()):
                msg.append(prop.write_bdf(size, card_writer))

            # BCs
            for (unused_key, bcs) in sorted(self.bcs.iteritems()):
                for bc in bcs:  # list
                    msg.append(bc.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

    def _write_thermal_materials(self, outfile, size, card_writer):
        """Writes the thermal materials in a sorted order"""
        if self.thermalMaterials:
            msg = ['$THERMAL MATERIALS\n']
            for (mid, material) in sorted(self.thermalMaterials.iteritems()):
                msg.append(material.write_bdf(size, card_writer))
            outfile.write(''.join(msg))

