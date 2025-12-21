from pyNastran.dev.bdf_vectorized.bdf_interface2.attributes import BDFAttributes


class AddCard(BDFAttributes):
    """defines methods to add card objects to the BDF"""
    def __init__(self):
        BDFAttributes.__init__(self)

    def _add_param_object(self, param, allow_overwrites=False):
        """adds a PARAM object"""
        key = param.key
        if key in self.params and not allow_overwrites:
            if not param == self.params[key]:
                # if param.key in self.params:
                #     msg = 'key=%s param=%s old_param=%s' % (key, param, self.params[key])
                #     raise KeyError(msg)
                self.log.warning('key=%s param=%s old_param=%s' %
                                 (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
            self._type_to_id_map[param.type].append(key)

    def _add_plotel_object(self, elem, allow_overwrites=False):
        """adds an PLOTEL object"""
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if not allow_overwrites:
            # if key in self.elements:
            #     if elem ==self.elements[key]:
            #         self._duplicate['elements'].append(elem)
            #         if self._stop_on_duplicate_error:
            #             self.pop_parse_errors()
            if key in self.plotels:
                if not elem == self.plotels[key]:
                    assert elem.eid not in self.plotels, 'eid=%s\nold_element=\n%snew_element=\n%s' % (elem.eid, self.plotels[elem.eid], elem)
        self.plotels[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_doptprm_object(self, doptprm, comment=''):
        """adds a DOPTPRM"""
        self.doptprm = doptprm

    def _add_rigid_element_object(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in self.rigid_elements and not allow_overwrites:
            assert elem.eid not in self.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigid_elements[elem.eid], elem)
        self.rigid_elements[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_set_object(self, set_obj):
        """adds an SET1/SET3 object"""
        key = set_obj.sid
        assert key >= 0
        if key in self.sets:
            self.sets[key].add_set(set_obj)
        else:
            self.sets[key] = set_obj
            self._type_to_id_map[set_obj.type].append(key)

    def _add_aset_object(self, set_obj):
        """adds an ASET/ASET1 object"""
        self.asets.append(set_obj)

    def _add_omit_object(self, set_obj):
        """adds an OMIT/OMIT1 object"""
        self.omits.append(set_obj)

    def _add_bset_object(self, set_obj):
        """adds an BSET/BSET1 object"""
        self.bsets.append(set_obj)

    def _add_cset_object(self, set_obj):
        """adds an CSET/USET1 object"""
        self.csets.append(set_obj)

    def _add_qset_object(self, set_obj):
        """adds an QSET/QSET1 object"""
        self.qsets.append(set_obj)

    def _add_uset_object(self, set_obj):
        """adds an USET/USET1 object"""
        key = set_obj.name
        if key in self.usets:
            self.usets[key].append(set_obj)
        else:
            self.usets[key] = [set_obj]

    def _add_sebset_object(self, set_obj):
        """adds an SEBSET/SEBSET1 object"""
        self.se_bsets.append(set_obj)

    def _add_secset_object(self, set_obj):
        """adds an SECSET/SECSTE1 object"""
        self.se_csets.append(set_obj)

    def _add_seqset_object(self, set_obj):
        """adds an SEQSET/SEQSET1 object"""
        self.se_qsets.append(set_obj)

    def _add_seuset_object(self, set_obj):
        """adds an SEUSET/SEUSET1 object"""
        key = set_obj.name
        if key in self.se_usets:
            self.se_usets[key].append(set_obj)
        else:
            self.se_usets[key] = [set_obj]

    def _add_seset_object(self, set_obj):
        """adds an SESET object"""
        key = set_obj.seid
        assert key >= 0
        if key in self.se_sets:
            old_set = self.se_sets[key]
            set_obj.add_seset(old_set)
        self.se_sets[key] = set_obj
        self._type_to_id_map[set_obj.type].append(key)

    def _add_method_object(self, method, allow_overwrites=False):
        """adds a EIGR/EIGRL object"""
        key = method.sid
        if key in self.methods and not allow_overwrites:
            if not method == self.methods[key]:
                assert key not in self.methods, 'sid=%s\nold_method=\n%snew_method=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method
            self._type_to_id_map[method.type].append(key)
