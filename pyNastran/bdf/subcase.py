"""Subcase creation/extraction class"""
from typing import List, Dict, Tuple, Any
from numpy import ndarray

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils import object_attributes, deprecated

from pyNastran.bdf.bdf_interface.case_control_cards import CLASS_MAP
from pyNastran.bdf.bdf_interface.subcase_utils import (
    write_stress_type, write_set, expand_thru_case_control)


INT_CARDS = (
    # these are cards that look like:
    #    LOAD = 6
    'SPC', 'MPC', 'TRIM', 'FMETHOD', 'METHOD', 'LOAD',
    'SUPORT', 'SUPORT1', 'TEMPERATURE(INITIAL)', 'TEMPERATURE(LOAD)',
    'DLOAD', 'MFLUID', 'CLOAD', 'NLPARM', 'CMETHOD',
    'FREQUENCY', 'TSTEP', 'TSTEPNL', 'SDAMPING', 'DESOBJ',
    'TEMPERATURE(INIT)', 'RANDOM', 'DESSUB', 'ADAPT', 'MAXLINES',
    'TFL', 'DESGLB', 'SMETHOD', 'DYNRED', 'GUST', 'TEMPERATURE(MATE)',
    'OTIME', 'NONLINEAR', 'AUXM', 'IC', 'BC', 'OUTRCV', 'DIVERG',
    'DATAREC', 'TEMPERATURE(BOTH)', 'DEFORM', 'MODES', 'CASE',
    'SEDR', 'SELG', 'SEFINAL', 'SEKR', 'TEMPERATURE(ESTIMATE)',
    'GPSDCON', 'AUXMODEL',
    'MODTRAK', 'OFREQ', 'DRSPAN', 'OMODES', 'ADACT', 'SERESP', 'STATSUB',
    'CURVESYM', 'ELSDCON', 'CSSCHD', 'NSM', 'TSTRU', 'RANDVAR',
    'RGYRO', 'SELR', 'TEMPERATURE(ESTI)', 'RCROSS', 'SERE', 'SEMR',
)

PLOTTABLE_TYPES = (
    # these are types that look like:
    #    STRESS(PLOT,PRINT,PUNCH,SORT1) = ALL
    # they all support PLOT
    'STRESS', 'STRAIN', 'SPCFORCES', 'DISPLACEMENT', 'MPCFORCES', 'SVECTOR',
    'VELOCITY', 'ACCELERATION', 'FORCE', 'ESE', 'OLOAD', 'SEALL', 'GPFORCE',
    'GPSTRESS', 'GPSTRAIN', 'FLUX', 'AEROF', 'THERMAL', 'STRFIELD',
    'NOUTPUT', 'SEDV', 'APRES', 'HTFLOW', 'NLSTRESS', 'GPKE',
    'SACCELERATION', 'SDISPLACEMENT', 'SEMG', 'HARMONICS', 'PRESSURE', 'VUGRID',
    'ELSUM', 'SVELOCITY', 'STRFIELD REAL', 'SENSITY', 'MONITOR',
    'NLLOAD', 'GPSDCON', 'BOUTPUT',
)


class Subcase:
    """
    Subcase creation/extraction class
    """
    allowed_param_types = [
        'SET-type', 'CSV-type', 'SUBCASE-type', 'KEY-type', 'STRESS-type', 'STRING-type',
        'OBJ-type',
    ]
    solCodeMap = {
        1: 101,
        21: 101,
        24: 101,
        26: 101,
        61: 101,
        64: 106,  # correct
        66: 106,  # correct
        68: 106,  # correct
        76: 101,
        99: 129,  # correct
        144: 101,  # correct
        187: 101,
    }

    def __init__(self, id=0):
        self.id = id
        self.params = {}
        self.sol = None
        self.log = None
        #print("\n***adding subcase %s***" % self.id)

    def load_hdf5_file(self, hdf5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast

        keys = list(hdf5_file.keys())
        for key in keys:
            #print(key)
            group = hdf5_file[key]
            if key in ['id']: # scalars
                value = _cast(group)
                setattr(self, key, value)
            elif key == 'params':
                #print(group)
                #print(group.keys())
                for group_key in group.keys():
                    #self.log.debug('%s %s' % (group_key, key))
                    value, options, param_type = _load_hdf5_param(group, group_key, encoding)

                    #self.log.debug('%s (%s, %s, %s)' % (key, value, options, param_type))
                    if isinstance(options, list):
                        options = [
                            option.decode(encoding) if isinstance(option, bytes) else option
                            for option in options]

                    self.params[group_key] = (value, options, param_type)
                    str(self)
            else:  # pragma: no cover
                raise RuntimeError(f'failed exporting Subcase/{key}')

    def export_to_hdf5(self, hdf5_file, encoding):
        keys_to_skip = ['log', 'solCodeMap', 'allowed_param_types']
        h5attrs = object_attributes(self, mode='both', keys_to_skip=keys_to_skip)
        #print('Subcase %i' % self.id)
        for h5attr in h5attrs:
            value = getattr(self, h5attr)
            if h5attr in ['id']: # scalars
                # simple export
                hdf5_file.create_dataset(h5attr, data=value)
            elif h5attr in ['sol']: # scalars/None
                if value is None:
                    continue
                hdf5_file.create_dataset(h5attr, data=value)
            elif h5attr in ['params']:
                if len(value) == 0:
                    continue
                keys = list(self.params.keys())
                params_group = hdf5_file.create_group('params')
                #print('keys =', keys)
                unused_keys_bytes = [key.encode(encoding) for key in keys]
                #params_group.create_dataset('keys', data=keys_bytes)
                for key, (value, options, param_type) in self.params.items():
                    #print('  %-14s: %-8r %r %r' % (key, value, options, param_type))
                    if key == '':
                        sub_group = params_group.create_group('blank')
                        sub_group.create_dataset('value', data=value)

                    else:
                        #print('key = %r' % key)
                        sub_group = params_group.create_group(key)
                        if value is not None:
                            if isinstance(value, list):
                                value_bytes = [
                                    valuei.encode(encoding) if isinstance(valuei, str) else valuei
                                    for valuei in value]
                                sub_group.create_dataset('value', data=value_bytes)
                            elif isinstance(value, (integer_types, float, str)):
                                sub_group.create_dataset('value', data=value)
                            elif hasattr(value, 'export_to_hdf5'):
                                sub_groupi = sub_group.create_group('object')
                                sub_groupi.attrs['type'] = key
                                value.export_to_hdf5(sub_groupi, encoding)
                            else:
                                print(f'value = {value!r}')
                                raise NotImplementedError(value)

                    if param_type is not None:
                        sub_group.create_dataset('param_type', data=param_type)

                    if options is not None:
                        if isinstance(options, list):
                            options_bytes = [
                                option.encode(encoding) if isinstance(option, str) else option
                                for option in options]
                            sub_group.create_dataset('options', data=options_bytes)
                        else:
                            sub_group.create_dataset('options', data=options)

            #if h5attr in ['_begin_count', 'debug', 'write_begin_bulk']: # scalars
                ## do nothing on purpose
                #hdf5_file.create_dataset(h5attr, data=value)
            #elif h5attr in ['reject_lines', 'begin_bulk', 'lines', 'output_lines']:
                # lists of strings
                #if len(value) == 0:
                    #continue
                #value_bytes = [line.encode(encoding) for line in value]
                ##print(value_bytes)
                #hdf5_file.create_dataset(h5attr, data=value_bytes)
            #elif h5attr == 'subcases':
                #keys = list(self.subcases.keys())
                #subcase_group = hdf5_file.create_group('subcases')
                #subcase_group.create_dataset('keys', data=keys)
                #for key, subcase in self.subcases.items():
                    #sub_group = subcase_group.create_group(str(key))
                    #subcase.export_to_hdf5(subcase_group, encoding)
            else:  # pragma: no cover
                print(key, value)
                raise RuntimeError(f'cant export to hdf5 Subcase/{h5attr}')


    def __deepcopy__(self, memo):
        """
        Custom method for copy.deepcopy to improve speed by more than
        2x (tested with timeit).

        Default method is a bit slow for a list of lists and can take a
        long time to read a bdf with many subcases this method removes
        some of the overhead if the subcase is the default subcase, it
        is shallow copied instead this greatly improves bdf read speed,
        since it avoids deepcopying large sets defined in the default
        subcase for every subcase that is defined.
        """
        _copy = self.__copy__()
        memo[id(self)] = _copy

        if _copy.id == 0:
            return _copy

        def _deepcopy(lst):
            """
            Deep copies objects that aren't lists; references lists.
            """
            _cpy = list(lst)

            for i in range(len(_cpy)):
                cpyi = _cpy[i]
                if isinstance(cpyi, list):
                    _cpy[i] = _deepcopy(cpyi)
            return _cpy

        params = _copy.params
        for key, val in self.params.items():
            if isinstance(val, list):
                val = _deepcopy(val)
            params[key] = val
        return _copy

    def __copy__(self):
        _copy = self.__class__()
        _copy.id = self.id
        _copy.sol = self.sol
        _copy.log = self.log
        _copy.params.update(self.params)
        return _copy

    def deprecated(self, old_name: str, new_name: str, deprecated_version: str) -> None:
        """
        Throws a deprecation message and crashes if past a specific version.

        Parameters
        ----------
        old_name : str
            the old function name
        new_name : str
            the new function name
        deprecated_version : float
            the version the method was first deprecated in
        """
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def add_op2_data(self, data_code, msg, log):
        """
        >>> self.subcase.add_op2_data(self.data_code, 'VECTOR')
        """
        assert log is not None, log
        #subtable_name = data_code['subtable_name']
        table_name = data_code['table_name']
        if not isinstance(table_name, str):
            # table_name is a byte string
            table_name = table_name.decode('latin1')
        else:
            raise NotImplementedError(f'table_name={table_name!r}')

        table_code = data_code['table_code']
        unused_sort_code = data_code['sort_code']
        unused_device_code = data_code['device_code']

        sort_str = 'SORT1'
        if table_name in {'OESNLBR2'}:
            sort_str = 'SORT1'

        #print(data_code)
        #print('table_name=%r table_code=%s sort_code=%r device_code=%r' % (
            #table_name, table_code, sort_code, device_code))
        table_name = table_name.strip()
        #if 'TITLE' in
        #print(data_code)
        #print(f'table_name={table_name!r} type={type(table_name)}')
        options = []
        if data_code['title']:
            self.add('TITLE', data_code['title'], [], 'STRING-type')
        if data_code['subtitle']:
            self.add('SUBTITLE', data_code['subtitle'], [], 'STRING-type')
        if data_code['label']:
            self.add('LABEL', data_code['label'], [], 'STRING-type')

        #OUGF1, # OUG1F
        if table_name in ['OUGV1', 'BOUGV1', 'OUGV2', 'OUG1', 'OUGV1PAT', 'OUGMC1', 'OUGMC2',
                          'OUGF1', 'OUGF2', 'BOUGF1', ]:
            # OUG1F - acoustic displacements
            if table_code == 1:
                thermal = data_code['thermal']
                if thermal == 0:
                    self.add('DISPLACEMENT', 'ALL', options, 'STRESS-type')
                elif thermal == 1:
                    self.add('ANALYSIS', 'HEAT', [], 'KEY-type')
                else:
                    self._write_op2_error_msg(log, self.log, msg, data_code)
            elif table_code == 7:
                self.add('VECTOR', 'ALL', options, 'STRESS-type')
            elif table_code == 10:
                self.add('VELOCITY', 'ALL', options, 'STRESS-type')
            elif table_code == 11:
                self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
            elif table_code == 44:  # OUGMC1
                thermal = data_code['thermal']
                assert thermal == 0, thermal
                self.add('DISPLACEMENT', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)

        elif table_name == 'OAG1':
            if table_code == 11:
                thermal = data_code['thermal']
                assert thermal == 0, data_code
                self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)

        # OAG-random
        elif table_name in ['OAGPSD1', 'OAGPSD2']:
            options = ['PSDF']
            self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
        elif table_name in ['OAGCRM1', 'OAGCRM2']:
            options = ['CRM']
            self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
        elif table_name in ['OAGRMS1', 'OAGRMS2']:
            options = ['RMS']
            self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
        elif table_name in ['OAGNO1', 'OAGNO2']:
            options = ['NO']
            self.add('ACCELERATION', 'ALL', options, 'STRESS-type')

        elif table_name == 'TOUGV1':
            thermal = data_code['thermal']
            if thermal == 1:
                self.add('ANALYSIS', 'HEAT', options, 'KEY-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['ROUGV1', 'ROUGV2']:
            thermal = data_code['thermal']
            if thermal == 0:
                self.add('DISPLACEMENT', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OPHIG', 'BOPHIG', 'BOPHIGF']:
            if table_code == 7:
                self.add('ANALYSIS', 'HEAT', [], 'KEY-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name == 'OUPV1':
            if table_code == 1:
                self.add('SDISPLACEMENT', 'ALL', options, 'STRESS-type')
            elif table_code == 10:
                self.add('SVELOCITY', 'ALL', options, 'STRESS-type')
            elif table_code == 11:
                self.add('SACCELERATION', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)

        elif table_name == 'OPHSA':
            if table_code == 14:
                self.add('SVECTOR', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OUXY1', 'OUXY2']:
            if table_code == 15:
                self.add('SDISPLACEMENT', 'ALL', options, 'STRESS-type')
            elif table_code == 16:
                self.add('SVELOCITY', 'ALL', options, 'STRESS-type')
            elif table_code == 17:
                self.add('SACCELERATION', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)


        elif table_name in ['OQG1', 'OQG2', 'OQGV1']:
            if table_code == 3:
                self.add('SPCFORCES', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OEF1X', 'OEF1', 'RAFCONS', 'RAFEATC']:
            if table_code in [4]:
                self.add('FORCE', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OEF2']:
            options.append('SORT2')
            if table_code == 4:
                self.add('FORCE', 'ALL', options, 'STRESS-type')
            else:
                self._write_op2_error_msg(log, self.log, msg, data_code)

        elif table_name in ['OEFATO1', 'OEFATO2']:
            options.append('PSDF')
            self.add('FORCE', 'ALL', options, 'STRESS-type')
        elif table_name in ['OEFCRM1', 'OEFCRM2']:
            options.append('CRM')
            self.add('FORCE', 'ALL', options, 'STRESS-type')
        elif table_name in ['OEFRMS1', 'OEFRMS2']:
            options.append('RMS')
            self.add('FORCE', 'ALL', options, 'STRESS-type')
        elif table_name in ['OEFNO1', 'OEFNO2']:
            options.append('NO')
            self.add('FORCE', 'ALL', options, 'STRESS-type')
        elif table_name in ['OEFPSD1', 'OEFPSD2']:
            options.append('PSDF')
            self.add('FORCE', 'ALL', options, 'STRESS-type')

        elif table_name in ['OESATO1', 'OESATO2']:
            options.append('PSDF')
            self.add('STRESS', 'ALL', options, 'STRESS-type')
        elif table_name in ['OESCRM1', 'OESCRM2']:
            options.append('CRM')
            self.add('STRESS', 'ALL', options, 'STRESS-type')
        elif table_name in ['OESRMS1', 'OESRMS2', 'OESXRMS1']:
            options.append('RMS')
            self.add('STRESS', 'ALL', options, 'STRESS-type')
        elif table_name in ['OESNO1', 'OESNO2', 'OESXNO1']:
            options.append('NO')
            self.add('STRESS', 'ALL', options, 'STRESS-type')
        elif table_name in ['OESPSD1', 'OESPSD2', 'OESPSD2C']:
            options.append('PSDF')
            self.add('STRESS', 'ALL', options, 'STRESS-type')

        elif table_name in ['OSTRATO1', 'OSTRATO2']:
            options.append('PSDF')
            self.add('STRAIN', 'ALL', options, 'STRESS-type')
        elif table_name in ['OSTRCRM1', 'OSTRCRM2']:
            options.append('CRM')
            self.add('STRAIN', 'ALL', options, 'STRESS-type')
        elif table_name in ['OSTRRMS1', 'OSTRRMS2']:
            options.append('RMS')
            self.add('STRAIN', 'ALL', options, 'STRESS-type')
        elif table_name in ['OSTRNO1', 'OSTRNO2']:
            options.append('NO')
            self.add('STRAIN', 'ALL', options, 'STRESS-type')
        elif table_name in ['OSTRPSD1', 'OSTRPSD2',
                            'OSTPSD2C']:
            options.append('PSDF')
            self.add('STRAIN', 'ALL', options, 'STRESS-type')

        elif table_name in ['OEFIT', 'OEFITSTN']:
            if table_code in [25]:
                self.add('FORCE', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OQMG1', 'OQMG2']:
            if table_code in [3, 39]:
                self.add('MPCFORCES', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OGPFB1', 'RAGCONS', 'RAGEATC']:
            if table_code == 19:
                self.add('GPFORCE', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)

        # stress
        elif table_name in ['OES1', 'OES1X', 'OES1X1', 'OES1C', 'OESCP',
                            'OESNL2', 'OESNLXD', 'OESNLXR', 'OESNLBR', 'OESTRCP',
                            'OESVM1', 'OESVM1C', 'OESNL1X',
                            'OESNLXR2', 'RASCONS', 'RASEATC', 'OESNLBR2']:
            #assert data_code['is_stress_flag'] == True, data_code
            options.append(sort_str)
            if table_code == 5:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OES2', 'OES2C', 'OESVM2', ]:
            options.append('SORT2')
            if table_code == 5:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OESXRM1C']:
            if table_code == 805:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OESXNO1C']:
            if table_code == 905:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)

        elif table_name in ['OSTR2', 'OSTR2C']:
            options.append('SORT2')
            if table_code == 5:
                self.add('STRAIN', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)

        elif table_name in ['OESRT']:
            #assert data_code['is_stress_flag'] == True, data_code
            if table_code in [25, 89]:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OCRUG']:
            if table_code in [1]:
                self.add('DISP', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)

        # strain
        elif table_name in ['OSTR1X', 'OSTR1C', 'OSTR1', 'RAECONS', 'RAEEATC']:
            assert data_code['is_strain_flag'] is True, data_code
            if table_code == 5:
                self.add('STRAIN', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)
        elif table_name in ['OSTRVM1', 'OSTRVM1C', 'OSTRVM2']:
            #assert data_code['is_stress_flag'] == True, data_code
            if table_code == 5:
                self.add('STRAIN', 'ALL', options, 'STRESS-type')
            else:  # pragma: no cover
                self._write_op2_error_msg(log, self.log, msg, data_code)

        # special tables
        elif table_name in ['RADCONS', 'RADEFFM', 'RADEATC', 'RAPEATC', 'RAQEATC', 'RADCONS',
                            'RASEATC', 'RAFEATC', 'RAEEATC', 'RANEATC', 'RAGEATC', 'RAQCONS',
                            'RAPCONS']:
            pass
        elif table_name in ['OUGPSD2',
                            'OSTRNO1', 'OSTNO1C', 'OSTRNO1C',
                            'OSTRMS1C', 'OSTRRMS1', 'OSTRRMS1C',
                            'OQMPSD2']:
            pass
        else:  # pragma: no cover
            self._write_op2_error_msg(log, self.log, msg, data_code)
        #print(self)

    def _write_op2_error_msg(self, log, log_error, msg, data_code):
        if log is not None:
            log.error(msg)
            log.error(data_code)
        elif log_error is not None:
            log_error.error(msg)
            log_error.error(data_code)
        else:  # pragma: no cover
            # log_error is None
            print('Error calling subcase.add_op2_data...')
            print(msg)
            print(data_code)
            raise RuntimeError(data_code)
        raise RuntimeError(data_code)

    def __contains__(self, param_name: str) -> bool:
        """
        Checks to see if a parameter name is in the subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to check for

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           case_control = model.case_control_deck
           subcase1 = case_control.subcases[1]
           if 'LOAD' in subcase1:
               print('found LOAD for subcase 1')

        """
        if param_name in self.params:
            return True
        return False

    def has_parameter(self, *param_names) -> List[bool]:
        """
        Checks to see if one or more parameter names are in the subcase.

        Parameters
        ----------
        param_names : str; List[str]
            the case control parameters to check for

        Returns
        -------
        exists : List[bool]
            do the parameters exist

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           case_control = model.case_control_deck
           subcase1 = case_control.subcases[1]
           if any(subcase1.has_parameter('LOAD', 'TEMPERATURE(LOAD)')):
               print('found LOAD for subcase 1')

        """
        exists = [param_name.upper() in self.params
                  for param_name in param_names]
        return exists

    def __getitem__(self, param_name):
        """
        Gets the [value, options] for a subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to get

        Returns
        -------
        value : varies
            the value of the parameter
            'ALL' in STRESS(PLOT,PRINT) = ALL
        options : List[varies]
            the values in parentheses
            ['PLOT', 'PRINT'] in STRESS(PLOT,PRINT) = ALL

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           case_control = model.case_control_deck
           subcase1 = case_control.subcases[1]
           value, options = subcase1['LOAD']

        """
        return self.get_parameter(param_name)

    def suppress_output(self, suppress_to='PLOT'):
        """
        Replaces F06 printing with OP2 printing

        Converts:
            STRESS(PRINT,SORT1,REAL)
            FORCE(PRINT,PLOT,SORT1,REAL)

        to:
            STRESS(PLOT,SORT1,REAL)
            FORCE(PLOT,SORT1,REAL)

        .. warning:: needs more validation
        """
        for key, param in self.params.items():
            (unused_value, options, unused_param_type) = param
            if key in INT_CARDS or key in ('SUBTITLE', 'LABEL', 'TITLE', 'ECHO'):
                pass
            elif key in PLOTTABLE_TYPES:
                if suppress_to not in options:
                    param[1].append(suppress_to)
                if 'PRINT' in options:
                    param[1].remove('PRINT')
            else:
                raise NotImplementedError(key)

    def get_parameter(self, param_name, msg='', obj=False):
        """
        Gets the [value, options] for a subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to get
        obj : bool; default=False
            should the object be returned

        Returns
        -------
        value : varies
            the value of the parameter
            'ALL' in STRESS(PLOT,PRINT) = ALL
        options : List[varies]
            the values in parentheses
            ['PLOT', 'PRINT'] in STRESS(PLOT,PRINT) = ALL

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           case_control = model.case_control_deck
           subcase1 = case_control.subcases[1]
           value, options = subcase1['LOAD']

        """
        param_name = update_param_name(param_name)
        if param_name not in self.params:
            raise KeyError(f'{param_name} doesnt exist in subcase={self.id} in the case '
                           f'control deck{msg}.')
        value, options, param_type = self.params[param_name]
        #print('param_name=%r value=%s options=%s param_type=%r' % (
            #param_name, value, options, param_type))
        if param_type == 'OBJ-type' and not obj:
            return value.value, options
        return value, options

    def _validate_param_type(self, param_type) -> None:
        """checks to see if a valid parmater type is selected"""
        if param_type not in self.allowed_param_types:
            msg = (
                f'param_type={param_type!r} is not supported\n'
                f'   allowed_types={self.allowed_param_types}\n'
                '  - SET-type:     SET 5 = 1,2,3,4\n'
                '  - CSV-type:     PARAM,FIXEDB,-1\n'
                '  - KEY-type:     ANALYSIS = HEAT\n'
                '  - STRESS-type:  LOAD = 5\n'
                '  - STRESS-type:  STRESS = ALL\n'
                '  - STRESS-type:  STRESS(PLOT) = ALL\n'
                '  - STRESS-type:  DISP(PLOT) = ALL\n'
                '  - STRING-type:  TITLE = SOME TITLE\n'
                '  - OBJ-type:     ???\n'
            )
            raise TypeError(msg)

    def add(self, key: str, value: Any, options: List[Any], param_type: str):
        self._validate_param_type(param_type)
        self._add_data(key, value, options, param_type)

    def update(self, key, value, options, param_type):
        self._validate_param_type(param_type)
        assert key in self.params, f'key={key!r} is not in isubcase={self.id}'
        self._add_data(key, value, options, param_type)

    def add_set_from_values(self, set_id: int, values: List[int]):
        """
        Simple way to add SET cards

        Example
        -------
        >>> set_id = 42
        >>> values = [1, 2, 3]
        >>> subcase.add_set_from_values(set_id, values)
        >>> subcase

        SUBCASE 1
            SET 42 = 1, 2, 3

        """
        key = f'SET {set_id:d}'
        param_type = 'SET-type'
        assert isinstance(values, list), values
        self.params[key] = [values, set_id, param_type]

    def _add_data(self, key: str, value: Any, options: List[Any], param_type: str):
        """
        Adds the data to the subcase  in the KEY(OPTIONS)=VALUE style.

        Parameters
        ----------
        key : str
            the name

        >>> subcase._add_data(key, value, options, param_type)

        """
        key = update_param_name(key)
        if key == 'ANALYSIS' and value == 'FLUT':
            value = 'FLUTTER'

        #print("adding isubcase=%s key=%r value=%r options=%r "
              #"param_type=%r" %(self.id, key, value, options, param_type))
        if isinstance(value, str) and value.isdigit():
            value = int(value)

        if param_type == 'OBJ-type':
            self.params[key] = value
        else:
            (key, value, options) = self._simplify_data(key, value, options, param_type)
            self.params[key] = [value, options, param_type]

    def _simplify_data(self, key: str, value: Any, options: List[Any], param_type: str):
        if param_type == 'SET-type':
            #print("adding isubcase=%s key=%r value=%r options=%r "
                  #"param_type=%r" % (self.id, key, value, options, param_type))

            #print("adding isubcase=%s key=%r value=%r options=%r "
                  #"param_type=%r" % (self.id, key, value, options, param_type))
            values2 = expand_thru_case_control(value)

            assert isinstance(values2, list), type(values2)
            if isinstance(options, list):
                msg = f'invalid type for options={options!r} value; expected an integer; got a list'
                raise TypeError(msg)
            options = int(options)
            return (key, values2, options)

        elif param_type == 'CSV-type':
            #print(f'adding isubcase={self.id} key={key!r} value={value!r} options={options!r} '
                  #'param_type={param_type!r}')
            if value.isdigit():  # PARAM,DBFIXED,-1
                value = value
        #elif param_type == 'OBJ-type':
            ##self.params[value.type] = value
            #return value.type,
        elif param_type == 'OBJ-type':
            raise RuntimeError('this function should never be called with an OBJ-type...')
        else:
            #if 0:
                #a = 'key=%r' % key
                #b = 'value=%r' % value
                #c = 'options=%r' % options
                #d = 'param_type=%r' % param_type
                #print("_adding isubcase=%s %-18s %-12s %-12s %-12s" %(self.id, a, b, c, d))
            if isinstance(value, integer_types) or value is None:
                pass
            elif isinstance(value, (list, ndarray)):  # new???
                msg = f'invalid type for key={key} value; expected an integer; got a list'
                raise TypeError(msg)
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return key, value, options

    def get_op2_data(self, sol, solmap_to_value):
        self.sol = sol
        label = 'SUBCASE %s' % (self.id)
        op2_params = {
            'isubcase': None, 'tables': [], 'analysis_codes': [],
            'device_codes': [], 'sort_codes': [], 'table_codes': [],
            'label': label, 'subtitle': None, 'title': None,
            'format_codes': [], 'stress_codes': [], 'thermal': None}

        results = ['DISPLACEMENT', 'EKE', 'EDE', 'ELSDCON', 'ENTHALPY',
                   'EQUILIBRIUM', 'ESE', 'FLUX', 'FORCE', 'GPFORCE', 'GPKE',
                   'GPSDCON', 'GPSTRAIN', 'GPSTRESS', 'HOUTPUT', 'MODALKE',
                   'MODALSE', 'MPCFORCES', 'NLSTRESS', 'NOUTPUT', 'OLOAD',
                   'PFGRID', 'PFMODE', 'PFPANEL', 'RCROSS', 'RESVEC',
                   'SACCELERATION', 'SDISPACEMENT', 'SPCFORCES', 'STRAIN',
                   'STRESS', 'SVECTOR', 'SVELOCITY', 'THERMAL', 'VECTOR',
                   'VELOCITY', 'VUGRID', 'WEIGHTCHECK']

        # converts from solution 200 to solution 144
        if self.sol == 200 or 'ANALYSIS' in self.params:
            param = self.params['ANALYSIS']
            (value, options, param_type) = param

            sol = solmap_to_value[value.upper()]
            #print("***value=%s sol=%s" % (value, sol))
        else:  # leaves SOL the same
            sol = self.sol

        if sol in self.solCodeMap:  # reduces SOL 144 to SOL 101
            sol = self.solCodeMap[sol]

        for (key, param) in self.params.items():
            key = key.upper()
            (value, options, param_type) = param
            #msg = ("  -key=|%s| value=|%s| options=%s param_type=|%s|"
            #    % (key, value, options, param_type))

        thermal = 0
        for (key, param) in self.params.items():
            key = key.upper()
            (value, options, param_type) = param
            #msg = ("  *key=|%s| value=|%s| options=%s param_type=|%s|"
            #    % (key, value, options, param_type)
            #print(msg)
            #msg += self.printParam(key, param)
            if param_type == 'SUBCASE-type':
                op2_params['isubcase'].append(value)
            elif key in ['BEGIN', 'ECHO', 'ANALYSIS'] or 'SET' in key:
                pass
            elif key == 'TEMPERATURE':
                thermal = 1
            elif key in results:
                sort_code = get_sort_code(options, value)
                device_code = get_device_code(options, value)

                if key in ['STRESS', 'STRAIN']:
                    stress_code = get_stress_code(key, options, value)
                    op2_params['stress_codes'].append(stress_code)
                else:
                    op2_params['stress_codes'].append(0)

                format_code = get_format_code(options, value)
                table_code = get_table_code(sol, key, options)
                analysis_code = get_analysis_code(sol)

                approach_code = analysis_code * 10 + device_code
                tcode = table_code * 1000 + sort_code
                op2_params['tables'].append(key)

                op2_params['analysis_codes'].append(analysis_code)
                op2_params['approach_codes'].append(approach_code)
                op2_params['device_codes'].append(device_code)
                op2_params['format_codes'].append(format_code)
                op2_params['sort_codes'].append(sort_code)
                op2_params['table_codes'].append(table_code)
                op2_params['tcodes'].append(tcode)
                #analysisMethod = value

            #elif key in ['ADACT', 'ADAPT', 'AERCONFIG', 'TITLE', 'SUBTITLE',
            #             'LABEL', 'LOAD', 'SUPORT', 'SUPORT1', 'MPC', 'SPC',
            #            'TSTEPNL', 'NLPARM', 'TRIM', 'GUST', 'METHOD',
            #            'DESOBJ', 'DESSUB', 'FMETHOD', 'SEALL']:
            else:
                op2_params[key.lower()] = value

            #else:
            #    raise NotImplementedErrror('unsupported entry...\n%s' %(msg))

        op2_params['thermal'] = thermal

        #print("\nThe estimated results...")
        #for (key, value) in sorted(op2_params.items()):
            #if value is not None:
                #print("   key=%r value=%r" % (key, value))

    def print_param(self, key: str, param: Tuple[Any, Any, str]) -> str:
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        """
        msg = ''
        #msg += 'id=%s   ' %(self.id)
        (value, options, param_type) = param

        spaces = ''
        if self.id > 0:
            spaces = '    '

        #print('key=%s param=%s param_type=%s' % (key, param, param_type))
        if param_type == 'SUBCASE-type':
            if self.id > 0:
                msgi = 'SUBCASE %s\n' % (self.id)
                assert len(msgi) < 72, f'len(msg)={len(msgi)}; msg=\n{msgi}'
                msg += msgi
            #else:  global subcase ID=0 and is not printed
            #    pass
        elif param_type == 'KEY-type':
            #print (KEY-TYPE:  %r" % value)
            assert value is not None, param
            if ',' in value:
                sline = value.split(',')
                two_spaces = ',\n' + 2 * spaces
                msgi = spaces + two_spaces.join(sline) + '\n'
                assert len(msgi) < 68, f'len(msg)={len(msgi)}; msg=\n{msgi}'
                msg += msgi
            else:
                msgi = spaces + f'{value}\n'
                #assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                msg += msgi
        elif param_type == 'STRING-type':
            msgi = spaces + '%s = %s\n' % (key, value)
            if key not in ['TITLE', 'LABEL', 'SUBTITLE']:
                assert len(msgi) < 68, f'len(msg)={len(msgi)}; msg=\n{msgi}'
            msg += msgi
        elif param_type == 'CSV-type':
            msgi = spaces + '%s,%s,%s\n' % (key, value, options)
            assert len(msgi) < 68, f'len(msg)={len(msgi)}; msg=\n{msgi}'
            msg += msgi
        elif param_type == 'STRESS-type':
            msg += write_stress_type(key, options, value, spaces)

        elif param_type == 'SET-type':
            #: .. todo:: collapse data...not written yet
            msg += write_set(options, value, spaces)
        elif param_type == 'OBJ-type':
            msg += value.write(spaces)
        else:
            # SET-type is not supported yet...
            msg = ('\nkey=%r param=%r\n'
                   'allowed_params=[SET-type, STRESS-type, STRING-type, SUBCASE-type, KEY-type]\n'
                   'CSV-type     -> PARAM,FIXEDB,-1\n'
                   'KEY-type     -> ???\n' # the catch all
                   'SET-type     -> SET 99 = 1234\n'
                   'SUBCASE-type -> ???\n'
                   'STRESS-type  -> DISP(PLOT, SORT1)=ALL\n'
                   '                STRESS(PLOT, SORT1)=ALL\n'
                   'STRING-type  -> LABEL = A label\n'
                   '                TITLE = A title\n'
                   '                SUBTITLE = A subtitle\n'
                   ''% (key, param))
            raise NotImplementedError(msg)

        #print("msg = %r" % (msg))
        return msg

    #def cross_reference(self, model: BDF) -> None:
        #"""
        #Method crossReference:

        #Parameters
        #----------
        #model : BDF()
            #the BDF object

        #.. note:: this is not integrated and probably never will be as it's
          #not really that necessary.  it's only really useful when running an
          #analysis.
        #"""
        #raise NotImplementedError()
        #print("keys = %s" % (sorted(self.params.keys())))
        #if 'LOAD' in self.params:
            #load_id = self.params['LOAD'][0]
            #load_obj = model.loads[load_id]
            #load_obj.cross_reference(model)
        #if 'SUPORT' in self.params:
            #pass
        #if 'MPC' in self.params:
            ##mpc_id = self.params['MPC'][0]
            ##mpc_obj = model.mpcs[mpc_id]
            ##mpc_obj.cross_reference(model)
            #pass
        #if 'SPC' in self.params:
            ##spcID = self.params['SPC'][0]
            ##print "SPC ID = %r" % spcID
            ##spcObj = model.spcObject
            ##spcObj.cross_reference(spcID, model)
            #pass
        #if 'TSTEPNL' in self.params:
            #tstepnl_id = self.params['TSTEPNL'][0]
            #tstepnl_obj = model.tstepnl[tstepnl_id]
            #tstepnl_obj.cross_reference(model)
        #if 'NLPARM' in self.params:
            #nlparm_id = self.params['NLPARM'][0]
            #nlparm_obj = model.nlparms[nlparm_id]
            #nlparm_obj.cross_reference(model)
        #if 'TRIM' in self.params:
            #trim_id = self.params['TRIM'][0]
            #trim_obj = model.trims[trim_id]
            #trim_obj.cross_reference(model)
        #if 'GUST' in self.params:
            #gust_id = self.params['GUST'][0]
            #gust_obj = model.gusts[gust_id]
            #gust_obj.cross_reference(model)
        #if 'DLOAD' in self.params:  # ???
            #pass

    def finish_subcase(self):
        """
        Removes the subcase parameter from the subcase to avoid printing
        it in a funny spot
        """
        if 'SUBCASE' in self.params:
            del self.params['SUBCASE']
        #print "self.params %s = %s" %(self.id,self.params)

    def write_subcase(self, subcase0):
        """
        Internal method to print a subcase

        Parameters
        ----------
        subcase0 : Subcase()
            the global Subcase object

        Returns
        -------
        msg : str
            the string of the current Subcase
        """
        if self.id == 0:
            msg = str(self)
        else:
            msg = f'SUBCASE {self.id:d}\n'
            nparams = 0
            for (key, param) in self.subcase_sorted(self.params.items()):
                if key in subcase0.params and subcase0.params[key] == param:
                    pass  # dont write global subcase parameters
                else:
                    #print("key=%s param=%s" %(key, param))
                    #(unused_value, unused_options, unused_param_type) = param
                    #print("  *key=%r value=%r options=%s "
                          #"param_type=%r" % (key, value, options, param_type))
                    msg += self.print_param(key, param)
                    nparams += 1
                    #print ""
            if nparams == 0:
                for (key, param) in self.subcase_sorted(self.params.items()):
                    #print("key=%s param=%s" %(key, param))
                    #(unused_value, unused_options, unused_param_type) = param
                    #print("  *key=%r value=%r options=%s "
                          #"param_type=%r" % (key, value, options, param_type))
                    msg += self.print_param(key, param)
                    nparams += 1
                assert nparams > 0, f'No subcase parameters are defined for isubcase={self.id:d}...'

        return msg

    def subcase_sorted(self, lst):
        """
        Does a "smart" sort on the keys such that SET cards increment in
        numerical order.  Also puts the sets first.

        Parameters
        ----------
        lst : List[str]
            the list of subcase list objects (list_a)

        Returns
        -------
        list_b : List[str]
            the sorted version of list_a
        """
        # presort the list to put all the SET cards next to each other
        # instead of list_a.sort() as it allows lst to be any iterable
        lst = sorted(lst)

        i = 0  # needed in case the loop doesn't execute
        iset = None  # index of first SET card in the deck
        set_dict = {}
        list_before = []
        set_keys = []
        for (i, entry) in enumerate(lst):
            key = entry[0]  # type: str
            if 'SET' in key[0:3]:
                if key == 'SET':  # handles "SET = ALL"
                    key = 0
                else:  # handles "SET 100 = 1,2,3"
                    sline = key.split(' ')
                    try:
                        key = int(sline[1])
                    except:
                        msg = f'error caclulating key; sline={sline}'
                        raise RuntimeError(msg)

                # store the integer ID and the SET-type list
                set_dict[key] = entry
                set_keys.append(key)
            else:
                # only store the entries before the SET cards
                list_before.append(entry)
                if iset:
                    break

        # grab the other entries
        list_after = lst[i + 1:]

        # write the SET cards in a sorted order
        set_list = []
        for key in sorted(set_keys):
            set_list.append(set_dict[key])

        # combine all the cards
        list_b = set_list + list_before + list_after
        return list_b

    def __repr__(self):
        """
        Prints out EVERY entry in the subcase.  Skips parameters already in
        the global subcase.

        .. note:: this function is only used for debugging.
        """
        #msg = "-------SUBCASE %s-------\n" %(self.id)
        msg = ''
        if self.id > 0:
            msg += f'SUBCASE {self.id:d}\n'

        nparams = 0
        for key, param in self.subcase_sorted(self.params.items()):
            #(unused_value, unused_options, unused_param_type) = param
            #print('key=%r value=%s options=%s' % (key, value, options))
            msg += self.print_param(key, param)
            nparams += 1
        if self.id > 0:
            assert nparams > 0, f'No subcase parameters are defined for isubcase={self.id:d}...'
        return msg

def _load_hdf5_param(group, key, encoding):
    import h5py
    from pyNastran.utils.dict_to_h5py import _cast
    #print('-----------------------------------------')
    #print(type(key), key)
    sub_group = group[key]
    keys = list(sub_group.keys())
    #print('subgroup.keys() =', sub_group.keys())

    if key == 'blank':
        key = ''

    if 'options' in sub_group:
        keys.remove('options')
        options = _cast(sub_group['options'])
        if isinstance(options, (integer_types, str)):
            pass
        else:
            options = options.tolist()
            options = [
                option.decode(encoding) if isinstance(option, bytes) else option
                for option in options]
    else:
        options = None

    param_type = None
    if 'param_type' in sub_group:
        param_type = _cast(sub_group['param_type'])
        keys.remove('param_type')

    #print('param_type ', param_type)
    value = None
    if 'value' in sub_group:
        keys.remove('value')
        value = _cast(sub_group['value'])
        if isinstance(value, bytes):
            value = value.decode(encoding)
        elif isinstance(value, (integer_types, str)):
            pass
        else:
            value = value.tolist()

    elif 'object' in sub_group:
        value, options = _load_hdf5_object(key, keys, sub_group, encoding)


    if len(keys) > 0:
        #keyi = _cast(sub_group['key'])
        #print('keyi = %r' % keyi)
        raise RuntimeError(f'keys = {keys}')

    #print(value, options, param_type)
    return value, options, param_type

def _load_hdf5_object(key, keys, sub_group, encoding):
    import h5py
    from pyNastran.utils.dict_to_h5py import _cast
    keys.remove('object')
    sub_groupi = sub_group['object']

    Type = sub_groupi.attrs['type']
    cls = CLASS_MAP[Type]

    if hasattr(cls, 'load_hdf5'):
        class_obj, options = cls.load_hdf5(sub_groupi, 'utf8')
        value = class_obj
        return value, options

    use_data = True
    if 'options' in sub_groupi:
        options2 = _cast(sub_groupi['options']).tolist()
        value = _cast(sub_groupi['value'])
        #print('sub_keys =', sub_groupi, sub_groupi.keys())

        options_str = [
            option.decode(encoding) if isinstance(option, bytes) else option
            for option in options2]
        use_data = False

    data_group = sub_groupi['data']
    keys2 = _cast(data_group['keys']).tolist()

    h5_values = data_group['values']
    if isinstance(h5_values, h5py._hl.group.Group):
        values2 = [None] * len(keys2)
        for ih5 in h5_values.keys():
            ih5_int = int(ih5)
            h5_value = _cast(h5_values[ih5])
            values2[ih5_int] = h5_value
    else:
        values2 = _cast(h5_values).tolist()
    #print('data_keys =', data_group, data_group.keys())

    unused_keys_str = [
        keyi.decode(encoding) if isinstance(keyi, bytes) else keyi
        for keyi in keys2]
    unused_values_str = [
        valuei.decode(encoding) if isinstance(valuei, bytes) else valuei
        for valuei in values2]

    #print('keys2 =', keys2)
    #print('values2 =', values2)
    #print('options2 =', options2)

    if use_data:
        #print('keys2 =', keys2)
        #print('values2 =', values2)
        data = []
        for keyi, valuei in zip(keys2, values2):
            data.append((keyi, valuei))
        class_obj = cls(data)
        assert options is None, options
    else:
        class_obj = cls(key, value, options_str)
        options = options_str
    value = class_obj
    #print(class_obj)
    #class_obj.load_hdf5_file(hdf5_file, encoding)
    return value, options

def update_param_name(param_name):
    """
    Takes an abbreviated name and expands it so the user can type DISP or
    DISPLACEMENT and get the same answer

    Parameters
    ----------
    param_name : str
        the parameter name to be standardized (e.g. DISP vs. DIPLACEMENT)

    .. todo:: not a complete list
    """
    param_name = param_name.strip().upper()
    if param_name.startswith('ACCE'):
        param_name = 'ACCELERATION'
    elif param_name.startswith('DESO'):
        param_name = 'DESOBJ'
    elif param_name.startswith('DESS'):
        param_name = 'DESSUB'
    elif param_name.startswith('DISP'):
        param_name = 'DISPLACEMENT'
    elif param_name.startswith('EXPO'):
        param_name = 'EXPORTLID'
    elif param_name.startswith('ELFO'):
        param_name = 'FORCE'
    elif param_name.startswith('ELST'):
        param_name = 'STRESS' # or ELSTRESS
    elif param_name.startswith('FORC'):
        param_name = 'FORCE'
    elif param_name.startswith('FREQ'):
        param_name = 'FREQUENCY'
    elif param_name.startswith('GPFO'):
        param_name = 'GPFORCE'
    elif param_name == 'GPST':
        raise SyntaxError('invalid GPST stress/strain')
    elif param_name.startswith('HARMONIC'):
        param_name = 'HARMONICS'
    elif param_name.startswith('METH'):
        param_name = 'METHOD'
    elif param_name.startswith('MPCF'):
        param_name = 'MPCFORCES'
    elif param_name.startswith('NLPAR'):
        param_name = 'NLPARM'
    elif param_name.startswith('OLOA'):
        param_name = 'OLOAD'
    elif param_name.startswith('PRES'):
        param_name = 'PRESSURE'
    elif param_name.startswith('SDAMP'):
        param_name = 'SDAMPING'
    elif param_name.startswith('SDISP'):
        param_name = 'SDISPLACEMENT'
    elif param_name.startswith('SMETH'):
        param_name = 'SMETHOD'
    elif param_name.startswith('SPCF'):
        param_name = 'SPCFORCES'
    elif param_name.startswith('STRA'):
        param_name = 'STRAIN'
    elif param_name.startswith('STRE'):
        param_name = 'STRESS'
    elif param_name.startswith('SUBT'):
        param_name = 'SUBTITLE'
    elif param_name.startswith('SUPO'):
        param_name = 'SUPORT1'
    elif param_name.startswith('SVEC'):
        param_name = 'SVECTOR'
    elif param_name.startswith('SVELO'):
        param_name = 'SVELOCITY'
    elif param_name.startswith('THER'):
        param_name = 'THERMAL'
    elif param_name.startswith('VECT'):
        #param_name = 'PRESSURE' # or VECTOR
        param_name = 'DISPLACEMENT' # or VECTOR
    elif param_name.startswith('VELO'):
        param_name = 'VELOCITY'
    elif param_name.startswith('TITL'):
        param_name = 'TITLE'
    elif param_name.startswith('MAXLINE'):
        param_name = 'MAXLINES'
    elif param_name.startswith('LINE'):
        param_name = 'LINE'
    elif param_name.startswith('AXISYM'):
        param_name = 'AXISYMMETRIC'
    elif param_name.startswith('SUBSE'):
        param_name = 'SUBSEQ'
    elif param_name.startswith('XTIT'):
        param_name = 'XTITLE'
    elif param_name.startswith('YTIT'):
        param_name = 'YTITLE'
    elif param_name.startswith('SACCE'):
        param_name = 'SACCELERATION'
    elif param_name.startswith('GPSTRE'):
        param_name = 'GPSTRESS'
    elif param_name.startswith('GPSTR'):
        param_name = 'GPSTRAIN'
    elif param_name in ['DEFO', 'DEFOR']:
        param_name = 'DEFORM'
    elif param_name == 'TEMPERATURE(INIT)':
        param_name = 'TEMPERATURE(INITIAL)'

    #elif param_name.startswith('DFRE'):  param_name = 'D'

    # handled in caseControlDeck.py
    #elif param_name.startswith('TEMP'):  param_name = 'TEMPERATURE'
    #print '*param_name = ',param_name
    return param_name

def get_analysis_code(sol):
    """
    Maps the solution number to the OP2 analysis code.

     * 8 - post-buckling (maybe 7 depending on NLPARM???)

    # not important
     * 3/4 - differential stiffness (obsolete)
     * 11  - old geometric nonlinear statics
     * 12  - contran (???)

    .. todo:: verify
    """
    codes = {
        101 : 1,  # staics
        103 : 2,  # modes
        105 : 7,  # pre-buckling
        106 : 10, # nonlinear statics
        107 : 9,  # complex eigenvalues
        108 : 5,  # frequency
        111 : 5,
        112 : 6,
        114 : 1,
        115 : 2,
        116 : 7,
        118 : 5,
        129 : 6,  # nonlinear
        144 : 1,  # static aero
        145 : 1,
        146 : 1,  # flutter
        153 : 10,
        159 : 6,  # transient thermal
    }
    #print("sol=%s" % sol)
    approach_code = codes[sol]
    #print('approach_code = %s' % approach_code)
    return approach_code

def get_device_code(options: Any, unused_value: Any) -> int:
    """
    Gets the device code of a given set of options and value

    Parameters
    ----------
    options : list[int/float/str]
        the options for a parameter
    unused_value : int/float/str
        the value of the parameter

    Returns
    -------
    device_code : int
       The OP2 device code

       0 - No output
       1 - PRINT
       2 - PLOT
       3 - PRINT, PLOT
       4 - PUNCH
       5 - PRINT, PUNCH
       6 - PRINT, PLOT, PUNCH
    """
    device_code = 0
    if 'PRINT' in options:
        device_code += 1
    if 'PLOT' in options:
        device_code += 2
    if 'PUNCH' in options:
        device_code += 4
    device_code = max(device_code, 1)
    #if device_code==0:
    #    device_code=1  # PRINT
    return device_code

def get_table_code(sol, table_name, unused_options):
    """
    Gets the table code of a given parameter.  For example, the
    DISPLACMENT(PLOT,POST)=ALL makes an OUGV1 table and stores the
    displacement.  This has an OP2 table code of 1, unless you're running a
    modal solution, in which case it makes an OUGV1 table of eigenvectors
    and has a table code of 7.

    Parameters
    ----------
    options : list[int/float/str]
        the options for a parameter
    value : int/float/str
        the value of the parameter

    Returns
    -------
    table_code : int
       the OP2 table_code
    """
    if table_name in ['VECTOR', 'PRESSURE']:
        table_name = 'DISPLACEMENT'  # equivalent tables...

    key = (sol, table_name)
    tables = {
        # SOL, table_name      table_code
        (101, 'ACCELERATION'): 11,
        (103, 'ACCELERATION'): 11,
        (106, 'ACCELERATION'): 11,
        (107, 'ACCELERATION'): 11,
        (108, 'ACCELERATION'): 11,
        (129, 'ACCELERATION'): 11,
        #(144, 'ACCELERATION'): 11,
        (145, 'ACCELERATION'): 11,
        (146, 'ACCELERATION'): 11,

        (101, 'DISPLACEMENT'): 1,
        (103, 'DISPLACEMENT'): 7,  # VECTOR
        (105, 'DISPLACEMENT'): 7,
        (106, 'DISPLACEMENT'): 1,
        (107, 'DISPLACEMENT'): 7,
        (108, 'DISPLACEMENT'): 1,
        (109, 'DISPLACEMENT'): 1,
        (111, 'DISPLACEMENT'): 7,
        (112, 'DISPLACEMENT'): 1,
        (129, 'DISPLACEMENT'): 7,
        #(144, 'DISPLACEMENT'): 1,
        (145, 'DISPLACEMENT'): 1,
        (146, 'DISPLACEMENT'): 1,

        (101, 'ESE'): 18,  # energy
        (103, 'ESE'): 18,  # energy
        (105, 'ESE'): 18,  # energy
        (106, 'ESE'): 18,  # energy
        (107, 'ESE'): 18,  # energy
        (108, 'ESE'): 18,  # energy
        (109, 'ESE'): 18,  # energy
        (110, 'ESE'): 18,  # energy
        (111, 'ESE'): 18,  # energy
        (112, 'ESE'): 18,  # energy
        (145, 'ESE'): 18,  # energy
        (146, 'ESE'): 18,  # energy

        (101, 'FORCE'): 3,  # ???
        (103, 'FORCE'): 3,  # ???
        (105, 'FORCE'): 3,  # ???
        (106, 'FORCE'): 3,  # ???
        (107, 'FORCE'): 4,  # ???
        (108, 'FORCE'): 3,  # ???
        (111, 'FORCE'): 3,  # ???
        (112, 'FORCE'): 3,  # ???
        (129, 'FORCE'): 3,  # ???
        (145, 'FORCE'): 3,  # ???
        (146, 'FORCE'): 3,  # ???

        (101, 'GPFORCE'): 19,
        (105, 'GPFORCE'): 19,
        (106, 'GPFORCE'): 19,
        (107, 'GPFORCE'): 19,
        (108, 'GPFORCE'): 19,
        (111, 'GPFORCE'): 19,
        (112, 'GPFORCE'): 19,
        (129, 'GPFORCE'): 19,
        (145, 'GPFORCE'): 19,
        (146, 'GPFORCE'): 19,

        (101, 'GPSTRESS'): 20,
        (105, 'GPSTRESS'): 20,
        (106, 'GPSTRESS'): 20,
        (107, 'GPSTRESS'): 20,
        (108, 'GPSTRESS'): 20,
        (111, 'GPSTRESS'): 20,
        (112, 'GPSTRESS'): 20,
        (129, 'GPSTRESS'): 20,
        (145, 'GPSTRESS'): 20,
        (146, 'GPSTRESS'): 20,

        (101, 'GPSTRAIN'): 21,
        (105, 'GPSTRAIN'): 21,
        (106, 'GPSTRAIN'): 21,
        (107, 'GPSTRAIN'): 21,
        (108, 'GPSTRAIN'): 21,
        (111, 'GPSTRAIN'): 21,
        (112, 'GPSTRAIN'): 21,
        (129, 'GPSTRAIN'): 21,
        (145, 'GPSTRAIN'): 21,
        (146, 'GPSTRAIN'): 21,

        (101, 'MPCFORCES'): 3,
        (103, 'MPCFORCES'): 3,
        (106, 'MPCFORCES'): 3,
        (108, 'MPCFORCES'): 3,
        (112, 'MPCFORCES'): 3,
        (129, 'MPCFORCES'): 3,
        #(144, 'MPCFORCES'): 3,
        (145, 'MPCFORCES'): 3,
        (146, 'MPCFORCES'): 3,

        (101, 'OLOAD'): 2,
        (103, 'OLOAD'): 2,
        (105, 'OLOAD'): 2,
        (106, 'OLOAD'): 2,
        (107, 'OLOAD'): 2,
        (108, 'OLOAD'): 2,
        (111, 'OLOAD'): 2,
        (112, 'OLOAD'): 2,
        (129, 'OLOAD'): 2,
        #(144, 'OLOAD'): 2,
        (145, 'OLOAD'): 2,
        (146, 'OLOAD'): 2,

        (101, 'SPCFORCES'): 3,
        (103, 'SPCFORCES'): 3,
        (105, 'SPCFORCES'): 3,
        (106, 'SPCFORCES'): 3,
        (107, 'SPCFORCES'): 3,
        (108, 'SPCFORCES'): 3,
        (110, 'SPCFORCES'): 3,
        (111, 'SPCFORCES'): 3,
        (112, 'SPCFORCES'): 3,
        (129, 'SPCFORCES'): 3,
        #(144, 'SPCFORCES'): 3,
        (145, 'SPCFORCES'): 3,
        (146, 'SPCFORCES'): 3,

        (101, 'STRAIN'): 5,  # 5/20/21 ???
        (105, 'STRAIN'): 5,
        (106, 'STRAIN'): 5,
        (107, 'STRAIN'): 5,
        (108, 'STRAIN'): 5,
        (110, 'STRAIN'): 5,
        (111, 'STRAIN'): 5,
        (112, 'STRAIN'): 5,
        (129, 'STRAIN'): 5,
        (145, 'STRAIN'): 5,
        (146, 'STRAIN'): 5,

        (101, 'STRESS'): 5,  # 5/20/21 ???
        (103, 'STRESS'): 5,
        (105, 'STRESS'): 5,
        (106, 'STRESS'): 5,
        (107, 'STRESS'): 5,
        (108, 'STRESS'): 5,
        (111, 'STRESS'): 5,
        (112, 'STRESS'): 5,
        (129, 'STRESS'): 5,
        (145, 'STRESS'): 5,
        (146, 'STRESS'): 5,

        (145, 'SVECTOR'): 14,

        (101, 'FLUX'): 4,
        (103, 'FLUX'): 4,
        (106, 'FLUX'): 4,
        (112, 'FLUX'): 4,
        (108, 'FLUX'): 4,
        (153, 'FLUX'): 4,
        (159, 'FLUX'): 4,

        (101, 'THERMAL'): 3,  # 3/4 ???
        (159, 'THERMAL'): 3,  # 3/4 ???

        (101, 'VELOCITY'): 10,
        (103, 'VELOCITY'): 10,
        (106, 'VELOCITY'): 10,
        (107, 'VELOCITY'): 10,
        (108, 'VELOCITY'): 10,
        (111, 'VELOCITY'): 10,
        (112, 'VELOCITY'): 10,
        (129, 'VELOCITY'): 10,
        #(144, 'VELOCITY'): 10,
        (145, 'VELOCITY'): 10,
        (146, 'VELOCITY'): 10,

        (101, 'VUGRID'): 10,
    }
    #print("key=%s" % str(key))
    if key not in tables:
        raise KeyError(key)
    table_code = tables[key]
    return table_code

def get_sort_code(options, unused_value):
    """
    Gets the sort code of a given set of options and value

    Parameters
    ----------
    options : List[int/str]
        the options for a parameter
    unused_value : int; str
        the value of the parameter
    """
    sort_code = 0
    if 'COMPLEX' in options:
        sort_code += 1
    if 'SORT2' in options:
        sort_code += 2
    if 'RANDOM' in options:
        sort_code += 4
    return sort_code

def get_format_code(options: Any, unused_value: Any) -> int:
    """
    Gets the format code that will be used by the op2 based on
    the options.

    Parameters
    ----------
    options : list[int/float/str]
        the options for a parameter
    unused_value : int/float/str
        the value of the parameter

    .. todo::  not done...only supports REAL, IMAG, PHASE, not RANDOM
    """
    format_code = 0
    if 'REAL' in options:
        format_code += 1
    if 'IMAG' in options:
        format_code += 2
    if 'PHASE' in options:
        format_code += 4
    format_code = max(format_code, 1)
    return format_code

def get_stress_code(key: str, options: Dict[str, Any], unused_value: Any) -> int:
    """
    Method get_stress_code:

    .. note:: the individual element must take the stress_code and reduce
    it to what the element can return.  For example, for an isotropic
    CQUAD4 the fiber field doesnt mean anything.

    BAR - no von mises/fiber
    ISOTROPIC - no fiber

    .. todo:: how does the MATERIAL bit get turned on?  I'm assuming it's
              element dependent...
    """
    stress_code = 0
    if 'VONMISES' in options:
        stress_code += 1
    if key == 'STRAIN':
        stress_code += 10  # 2+8=10 - fields 2 and 4
    if 'FIBER' in options:
        stress_code += 4
    #if 'MATERIAL' in options:
    #    stress_code += 16  material coord (1) vs element (0)
    return stress_code
