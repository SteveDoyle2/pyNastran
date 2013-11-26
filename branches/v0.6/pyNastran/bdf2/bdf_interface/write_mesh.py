from numpy import array, unique, concatenate, intersect1d, where

class WriteMesh(object):
    def __init__(self):
        pass
    
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

    def _write_params(self, size):
        """
        Writes the PARAM cards
        :param self: the BDF object
        """
        msg = []
        if self.params:
            msg = ['$PARAMS\n']
            for (key, param) in sorted(self.params.iteritems()):
                msg.append(param.print_card(size))
        return ''.join(msg)

    def _write_elements_properties2(self, f, size):
        ptypes = [self.properties_shell.pshell,
                  self.properties_shell.pcomp,
                  self.properties_shell.pshear,
                  self.prod,
                  self.properties_solid.psolid,]
        
        #pid_to_prop = defaultdict(list)
        #for proptype in ptypes:
            #for pid, prop in proptype.iteritems():
                #pid_to_prop[pid]

    def _write_elements_properties(self, f, size):
        """
        Writes the elements and properties in and interspersed order
        """
        return self._write_elements_properties2(f, size)
        msg = []
        missing_properties = []

        eids_written = []
        #pids = sorted(self.properties.keys())
        
        ptypes = [self.properties_shell.pshell,
                  self.properties_shell.pcomp,
                  self.properties_shell.pshear,
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
                  [self.crod,])

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
            pids_unique.sort()
            if pids_unique:
                f.write('$ELEMENTS_WITH_PROPERTIES\n')

            for pid in pids_all:
                i = where(pid==pids)[0]
                eids2 = eids[i]

                for t in ptypes:
                    if t.n and pid in t.property_id:
                        t.write_bdf(f, size=size, pids=[pid])
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
                        t.write_bdf(f, size=size, eids=eids3)
                    except TypeError:
                        print "t.type =", t.type
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
                        t.write_bdf(f, size=size, pids=[pid])
            
        #..todo:: finish...
        f.write('$UNASSOCIATED_ELEMENTS\n')
        aaa
        # missing elements...

    def write_bdf(self, bdf_filename, interspersed=True, size=8):
        f = open(bdf_filename, 'wb')

        #==========================
        f.write(self._write_header())
        f.write(self._write_params(size))

        #==========================
        self.grid.write_bdf(f, size)

        interspersed = False
        if interspersed:
            #raise NotImplementedError('interspersed=False')
            self._write_elements_properties(f, size)
        else:
            #self.properties_springs.write_bdf(f)
            self.elements_springs.write_bdf(f)
            self.pelas.write_bdf(f, size)

            #self.properties_rods.write_bdf(f)
            #self.elements_rods.write_bdf(f)

            #self.properties_bars.write_bdf(f)
            #self.elements_bars.write_bdf(f)

            self.properties_shell.write_bdf(f, size)
            self.elements_shell.write_bdf(f)

            self.properties_solids.write_bdf(f, size)
            self.elements_solids.write_bdf(f)
        self.conrod.write_bdf(f, size)
        self._write_loads(f, size)
        self.materials.write_bdf(f, size)
        self._write_constraints(f, size)
        f.write(self._write_rejects(size))
        f.write('ENDDATA\n')

    def _write_constraints(self, f, size):
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


    def _write_loads(self, f, size, interspersed=False):
        #self.loadcase.write_bdf(f, size)
        for load_id, loads in sorted(self.load.iteritems()):
            for load in loads:
                load.write_bdf(f, size)

        for load_id, loads in sorted(self.dload.iteritems()):
            for load in loads:
                load.write_bdf(f, size)

        #self.loadset.write_bdf(f, size)
        self.force.write_bdf(f, size)
        #self.force1.write_bdf(f, size)
        #self.force2.write_bdf(f, size)
        self.moment.write_bdf(f, size)
        #self.moment1.write_bdf(f, size)
        #self.moment2.write_bdf(f, size)

    def _write_rejects(self, size):
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