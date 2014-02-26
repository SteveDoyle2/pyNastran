def _clean_lines(case_control, lines):
    """
    Removes comment characters defined by a *$*.

    :param case_control:  the CaseControlDeck object
    :param lines: the lines to clean.
    """
    lines2 = []
    for line in lines:
        line = line.strip(' \n\r').split('$')[0].rstrip()
        if line:
            lines2.append(line)
    return lines2

class Subcase(object):
    def __init__(self, id):
        self.id = id
        self.params = []

    def __repr__(self):
        spaces = '    ' if self.id != 0 else ''
        msg = ''
        if self.id > 0:
            msg += 'SUBCASE %i\n' % self.id

        end = ''
        for param in self.params:
            if param.__class__.__name__ == 'SET':
                msg += param.write(spaces) + '\n'
            else:
                end += param.write(spaces) + '\n'
        return msg + end

    def add(self, obj):
        self.params.append(obj)

    def has_parameter(self, param_name):
        names = [param.__class__.__name__ for param in self.params]
        #print names
        return True if param_name.upper() in names else False

    def get_parameter(self, param_name):
        names = [param.__class__.__name__ for param in self.params]
        i = names.index(param_name.upper())
        return self.params[i]

class CaseControlDeck(object):
    """
    CaseControlDeck parsing and extraction class
    """
    def __init__(self, lines, log=None):
        self.lines = lines
        self.subcases = {0: Subcase(id=0)}

        self._map = {
        }

        lines = _clean_lines(self, lines)
        self._read(lines)

    def _parse_new_subcase(self, line):
        try:
            sline = line.strip().split()
        except:
            raise RuntimeError(line)
        assert len(sline) == 2, sline
        subcase_id = int(sline[1])
        #print "Subcase ID=%s" % subcase_id
        subcase_id = self.create_new_subcase(subcase_id)
        return subcase_id

    def create_new_subcase(self, subcase_id):
        #print "making subcase = %r" % subcase_id
        assert isinstance(subcase_id, int), 'subcase_id=%r' % subcase_id
        assert subcase_id not in self.subcases, 'ids=%s' % self.subcases.keys()
        self.subcases[subcase_id] = Subcase(subcase_id)
        return subcase_id

    def has_parameter(self, subcase_id, param_name):
        assert isinstance(subcase_id, int), 'subcase_id=%r' % subcase_id
        assert subcase_id in self.subcases, 'ids=%s' % self.subcases.keys()
        subcase = self.subcases[subcase_id]
        return subcase.has_parameter(param_name)

    def get_subcase_parameter(self, subcase_id, param_name):
        assert isinstance(subcase_id, int), 'subcase_id=%r' % subcase_id
        assert subcase_id in self.subcases, 'ids=%s' % self.subcases.keys()
        subcase = self.subcases[subcase_id]
        return subcase.get_parameter(param_name)

    def _add_obj_to_subcase(self, subcase_id, obj):
        assert isinstance(subcase_id, int), 'subcase_id=%r' % subcase_id
        assert subcase_id in self.subcases, 'ids=%s' % self.subcases.keys()
        subcase = self.subcases[subcase_id]
        subcase.add(obj)

    def add_parameter_to_local_subcase(self, subcase_id, strobj):
        self._add_line(subcase_id, strobj)

    def add_parameter_to_global_subcase(self, strobj):
        subcase_id = 0
        self._add_line(subcase_id, strobj)

    def _read(self, lines):
        subcase_id = 0
        for line in lines:
            line = line.strip().upper()
            subcase_id = self._add_line(subcase_id, line)

    def _add_line(self, subcase_id, line):
        #print "%%", line
        if '=' in line:
            if line.startswith('SET'):
                #print "SET.."
                obj = SET(line)

            #======================
            # DISP(PLOT, PUNCH) = 5
            elif line.startswith('DISP') or line.startswith('VECTOR') or line.startswith('PRESSURE'):  # DISP / DISPLACEMENT
                obj = DISPLACEMENT(line)
            elif line.startswith('VELO'):  # VELOCITY
                obj = VELOCITY(line)
            elif line.startswith('ACCE'):  # ACCELERATION
                obj = ACCELERATION(line)

            elif line.startswith('FORC') or line.startswith('ELFO'): # FORCE / ELFORCE
                obj = FORCE(line)
            elif line.startswith('STRESS') or line.startswith('ELST'):  # ELSTRESS
                obj = STRESS(line)
            elif line.startswith('STRAIN'):
                obj = STRAIN(line)

            elif line.startswith('GPFO'):  # GPFORCE
                obj = GPFORCE(line)
            elif line.startswith('GPSTRESS'):
                obj = GPSTRESS(line)
            elif line.startswith('GPSTRAIN'):
                obj = GPSTRAIN(line)

            elif line.startswith('LOAD'):
                obj = LOAD(line)
            elif line.startswith('DLOAD'):
                obj = DLOAD(line)
            elif line.startswith('OLOAD'):
                obj = OLOAD(line)
            elif line.startswith('NLSTRESS'):
                obj = NLSTRESS(line)

            elif line.startswith('BCRESULTS'):
                obj = BCRESULTS(line)
            elif line.startswith('BGRESULTS'):
                obj = BGRESULTS(line)
            elif line.startswith('GKRESULTS'):
                obj = GKRESULTS(line)

            elif line.startswith('SHELLTHK'):
                obj = SHELLTHK(line)

            #========================
            # MPC = 5
            elif line.startswith('BC'):
                obj = BC(line)
            elif line.startswith('LINE'):
                obj = LINE(line)
            elif line.startswith('MAXLINES'):
                obj = MAXLINES(line)

            elif line.startswith('BC'):
                obj = BC(line)

            elif line.startswith('SPCFORCE'):
                obj = SPCFORCES(line)
            elif line.startswith('MPCFORCE'):
                obj = MPCFORCES(line)

            elif line.startswith('METH'):  # METHOD
                obj = METHOD(line)

            elif line.startswith('CMETHOD'):
                obj = CMETHOD(line)
            elif line.startswith('MPC'):
                obj = MPC(line)
            elif line.startswith('SPC'):
                obj = SPC(line)

            elif line.startswith('FREQ'): # FREQUENCY
                obj = FREQUENCY(line)
            elif line.startswith('TSTEP'):
                obj = TSTEP(line)
            elif line.startswith('TSTEPNL'):
                obj = TSTEPNL(line)
            elif line.startswith('NLPARM'):
                obj = NLPARM(line)
            elif line.startswith('MPRES'):
                obj = MPRES(line)

            #========================
            # OMODES = ALL
            # OMODES = 5
            elif line.startswith('OMODES'):
                obj = OMODES(line)
            elif line.startswith('OTIME'):
                obj = OTIME(line)
            elif line.startswith('STRFIELD'):
                obj = STRFIELD(line)

            #========================

            elif line.startswith('DESS'):  # DESSUB, DESS
                obj = DESSUB(line)
            elif line.startswith('DESO'): # DESOBJ / DESO
                obj = DESOBJ(line)

            #================================
            elif line.startswith('A2GG'):
                obj = M2GG(line)
            elif line.startswith('B2GG'):
                obj = B2GG(line)
            elif line.startswith('K2GG'):
                obj = K2GG(line)
            elif line.startswith('M2GG'):
                obj = M2GG(line)

            elif line.startswith('B2PP'):
                obj = B2PP(line)
            elif line.startswith('K2PP'):
                obj = K2PP(line)
            elif line.startswith('M2PP'):
                obj = M2PP(line)

            #================================

            elif line.startswith('TITLE'):
                obj = TITLE(line)
            elif line.startswith('SUBTITLE'):
                obj = SUBTITLE(line)
            elif line.startswith('LABEL'):
                obj = LABEL(line)

            elif line.startswith('DSAPRT'):
                obj = DSAPRT(line)
            elif line.startswith('EDE'):
                obj = EDE(line)
            elif line.startswith('EKE'):
                obj = EKE(line)
            elif line.startswith('ESE'):
                obj = ESE(line)
            elif line.startswith('EXTSEOUT'):
                obj = EXTSEOUT(line)
            elif line.startswith('GPKE'):
                obj = GPKE(line)
            elif line.startswith('GROUNDCHECK'):
                obj = GROUNDCHECK(line)


            elif line.startswith('ECHO'):
                obj = ECHO(line)

            elif line.startswith('ANALYSIS'):
                obj = ANALYSIS(line)
            else:
                print "*line =", line
                raise NotImplementedError(line)
            self._add_obj_to_subcase(subcase_id, obj)
        else:
            if line.startswith('SUBCASE'):
                subcase_id = self._parse_new_subcase(line)
                assert isinstance(subcase_id, int)
            elif line.startswith('SUPER'):
                obj = SUPER(line)
                self._add_obj_to_subcase(subcase_id, obj)
            elif line.startswith('BEGIN'):
                obj = BEGIN(line)
                self._add_obj_to_subcase(subcase_id, obj)
            elif line.startswith('OUTPUT'):
                obj = OUTPUT(line)
                self._add_obj_to_subcase(subcase_id, obj)
            else:
                print "***line =", line
                #raise NotImplementedError(line)
        return subcase_id

    def __repr__(self):
        msg = ''
        for subcase_id, subcase in sorted(self.subcases.iteritems()):
            msg += str(subcase)
        #msg += 'BEGIN BULK\n'
        return msg

#========================================================
class SUPER(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class OUTPUT(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

#========================================================
class BEGIN(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return self.line

#========================================================
class LABEL(object):
    def __init__(self, line):
        i = line.index('=')
        self.name = line[i+1:].strip()
    def write(self, spaces):
        return spaces + 'LABEL = %s' % self.name

class SUBTITLE(object):
    def __init__(self, line):
        i = line.index('=')
        self.name = line[i+1:].strip()
    def write(self, spaces):
        return spaces + 'SUBTITLE = %s' % self.name

class TITLE(object):
    def __init__(self, line):
        i = line.index('=')
        self.name = line[i+1:].strip()
    def write(self, spaces):
        return spaces + 'TITLE = %s' % self.name

#========================================================
class ECHO(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return self.line

class ANALYSIS(object):
    def __init__(self, line):
        self.line = line
        sline = line.split('=')
        assert len(sline) == 2, sline
        self.value = sline[1].strip()
        assert self.value in [
            # SOL 200
            'STATIC',
            'MODES',
            'BUCK',
            'DFREQ',
            'MFREQ',
            'MTRAN',
            'DCEIG',
            'MCEIG',
            'SAERO',
            'DIVERGE',
            'FLUTTER',

            # SOL 153/159
            'HEAT',
            'STRUCTURE',

            # SOL 400
            'NLSTATICS',
            'LNSTATICS',
            ], 'value=%r' % self.value
    def write(self, spaces):
        return spaces + self.line

#========================================================
class KeyParenKeyiEqValParenEqVal(object):
    def __init__(self, line):
        self.line = line
        print "line =", line

        if ')' in line:
            # KEY(OPTION) = 5
            left, right = line.split(')')
            key, options = left.split('(')
            right = right.lstrip(' ')
            soptions = options.split(',')
        else:
            # KEY = 5
            left, right = line.split('=')
            left = left.strip(' ')
            right = right.strip(' ')
            soptions = []


        if right:
            # KEY(OPTION) = 5
            # KEY = 5
            self.value = right
        else:
            # KEY(OPTION)
            # KEY
            self.value = None

        self.options = soptions
        #print "key=%r options=%r" % (key, options)

    def write(self, spaces):
        name = self.__class__.__name__
        if len(self.options) > 0:
            options = '(%s)' % ','.join(self.options)
        else:
            options = ''
        if self.value is not None:
            msg = spaces + '%s%s = %s' % (name, options, self.value)
        else:
            msg = spaces + '%s%s' % (name, options)
        return msg

class DSAPRT(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class EDE(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class EKE(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class ESE(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class EXTSEOUT(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class GPKE(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)

class GROUNDCHECK(KeyParenKeyiEqValParenEqVal):
    def __init__(self, line):
        KeyParenKeyiEqValParenEqVal.__init__(self, line)



#========================================================
class KeyOptionsEqValue(object):
    def __init__(self, line):
        sline = line.split('=')
        #print "sline =", sline
        assert len(sline) == 2, sline
        key_options = sline[0].strip().strip()
        self.value = sline[1].strip()
        #print "key_options = %r" % key_options

        if '(' in key_options:
            key, options = key_options.rstrip(')').split('(')
            soptions = options.split(',')
            #print "key=%r options=%r" % (key, options)
        else:
            key = key_options
            soptions = None

        self.key = key
        self.options = soptions
        self.line = line

    def write(self, spaces):
        name = self.__class__.__name__
        if self.options:
            msg = '%s%s(%s) = %s' % (spaces, name, ','.join(self.options), self.value)
        else:
            msg = '%s%s = %s' % (spaces, name, self.value)
        return msg

class SPCFORCES(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class MPCFORCES(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class VELOCITY(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class DISPLACEMENT(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class ACCELERATION(KeyOptionsEqValue):
    """
    ACCELERATION(SORT1/SORT2, PRINT/PUNCH/PLOT, REAL/IMAG/PHASE, PSD/ATOC/CRMS/RALL,
                 RPRINT/NORPRINT, RPUNCH, CID) = ALL, N, NONE
    """
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
        options2 = []
        self.sort1 = None
        self.sort2 = None
        self.plot = None
        self.pprint = None
        self.punch = None
        self.real = None
        self.phase = None
        for option in self.options:
            if option == 'SORT1':
                assert self.sort1 is None
                self.sort1 = True
                self.sort2 = False
            elif option == 'SORT2':
                assert self.sort1 is None
                self.sort1 = False
                self.sort2 = True

            elif option == 'PLOT':
                assert self.plot is None
                self.plot = True
                #self.print = False
                #self.punch = False
            elif option == 'PRINT':
                assert self.pprint is None
                #self.plot = False
                self.pprint = True
                #self.punch = False
            elif option == 'PUNCH':
                assert self.punch is None
                #self.plot = False
                #self.print = False
                self.punch = True
            elif option in ['REAL', 'IMAG']:
                assert self.real is None
                self.real = True
                self.phase = False
            elif option == 'PHASE':
                assert self.real is None
                self.real = False
                self.phase = True
            else:
                raise NotImplementedError(option)

class OLOAD(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class FORCE(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class STRESS(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

class STRAIN(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)


class BCRESULTS(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class BGRESULTS(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class ELSDCON(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class ENTHALPY(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class EQUILIBRIUM(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class FLUX(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class GPSDCON(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class GKRESULTS(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class GPFORCE(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class GPSTRAIN(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class GPSTRESS(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class HDOT(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class HOUTPUT(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)
class SHELLTHK(KeyOptionsEqValue):
    def __init__(self, line):
        KeyOptionsEqValue.__init__(self, line)

#========================================================
class KeyOptionsEqInt(object):
    def __init__(self, line):
        sline = line.split('=')
        assert len(sline) == 2, sline
        key_options = sline[0].strip().strip()
        self.value = int(sline[1].strip())

        if '(' in key_options:
            key, options = key_options.rstrip(')').split('(')
            soptions = options.split(',')
        else:
            key = key_options
            soptions = None

        self.key = key
        self.options = soptions
        self.line = line

    def write(self, spaces):
        name = self.__class__.__name__
        if self.options:
            msg = '%s%s(%s) = %i' % (spaces, name, ','.join(self.options), self.value)
        else:
            msg = '%s%s = %i' % (spaces, name, self.value)
        return msg

class IC(KeyOptionsEqInt):
    def __init__(self, line):
        KeyOptionsEqInt.__init__(self, line)
class METHOD(KeyOptionsEqInt):
    def __init__(self, line):
        KeyOptionsEqInt.__init__(self, line)

#========================================================
# DLOAD = 5

class KeyEqInt(object):
    def __init__(self, line):
        try:
            self.ID = int(line.split('=')[1])
        except:
            raise RuntimeError(line)

    def write(self, spaces):
        return spaces + '%s = %i' % (self.__class__.__name__, self.ID)

class ADAPT(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class AUXMODEL(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class BC(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class CSSCHD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class LOAD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class CLOAD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class DLOAD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class LOADSET(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class SPC(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class MPC(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class FREQUENCY(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class TSTEP(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class TSTEPNL(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)


class DEFORM(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class DESGLB(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class DESSUB(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class LINE(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class MAXLINES(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class CMETHOD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class FMETHOD(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

class DIVERG(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class DRSPAN(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class GUST(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class MFLUID(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class MODES(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)
class NLPARM(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)


#========================================================
#DESOBJ = 5
#DESOBJ(MIN/MAX) = 5
class DESOBJ(KeyEqInt):
    def __init__(self, line):
        KeyEqInt.__init__(self, line)

#========================================================
# SET = ALL
# SET = 5
# SET 3 = 6, 10, 12, 100 THRU 200
class SET(object):
    def __init__(self, line):
        sline = line[3:].split('=')
        ID = sline[0].strip()
        if ID:
            self.ID = int(sline[0])
        else:
            self.ID = None
        self.values = [value.strip() for value in sline[1].strip().split(',') ]

    def write(self, spaces):
        msg = ''
        if self.ID is None:
            msg = spaces + 'SET = %s' % ','.join(self.values)
        else:
            starter = 'SET %s = ' % self.ID
            msg2 = spaces + starter
            nChars = len(msg2)

            i = 0
            value = self.values
            while i < len(value):
                newString = '%s, ' % value[i]
                if len(msg2 + newString) > 70:
                    msg += msg2 + '\n'
                    msg2 = ' ' * nChars + newString
                else:
                    msg2 += newString
                i += 1

            msg += msg2.rstrip(' \n,')

            #msg = spaces + 'SET %i = %s' % (self.ID, ','.join(self.values))
        return msg

#========================================================
class A2GG(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class B2GG(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class K2GG(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class M2GG(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class B2PP(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class K2PP(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class M2PP(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

#========================================================
# STRFIELD = ALL
# STRFIELD = N

class OTIME(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class STRFIELD(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class AEROF(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class APRESSURE(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class DESVAR(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line

class GPRSORT(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line
class HTFLOW(object):
    def __init__(self, line):
        self.line = line
    def write(self, spaces):
        return spaces + self.line
#========================================================
# ADACT = ALL
# ADACT = NONE
# ADACT = 5

class ADACT(object):
    def __init__(self, line):
        sline = line.split('=')
        assert len(sline) == 2
        self.value = sline[1]
    def write(self, spaces):
        return spaces + 'ADACT = %s' % self.value

class HARMONICS(object):
    def __init__(self, line):
        sline = line.split('=')
        assert len(sline) == 2
        self.value = sline[1]
    def write(self, spaces):
        return spaces + 'HARMONICS = %s' % self.value

#========================================================
class AEUXREF(object):
    """
    AEUXREF = 5
    AEUXREF = TRIM
    """
    def __init__(self, line):
        sline = line.split('=')
        assert len(sline) == 2
        self.value = sline[1]
    def write(self, spaces):
        return spaces + 'AEUXREF = %s' % self.value

#========================================================
# skip these
# ADAMSMNF
# AECONFIG
# AESYMXY
# AESYMXZ
# AUTOSPC
# AUXCASE
# AXISYMMETRIC
# BCONTACT
# BOUTPUT
# CMSENERGY
# DATAREC
# DSYM
# ELSUM
# ENDTIME
# FLSFSEL
# FLSPOUT
# FLSTCNT
# MASTER
# MAXMIN
# MCFRACTION
# MEFFMASS
# MODALKE
# MODALSE

#========================================================

if __name__ == '__main__':  # pragma: no cover
    lines = [
        'SUBCASE 1',
        '    ACCELERATION(PLOT,PRINT,PHASE) = ALL',
        '    DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
        '    DLOAD = 32',
        '    M2GG = 111',
        '    SET 88  = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, 100 THRU 300',
        '    SET 99  = 1 THRU 10',
        '    SET 105 = 1.009, 10.2, 13.4, 14.0, 15.0',
        '    SET 111 = MAAX1,MAAX2',
        '    SET 1001 = 101/T1, 501/T3, 991/R3',
        '    SET = ALL',
        '    SPC = 42',
        '    TSTEPNL = 22',
        '    VELOCITY(PLOT,PRINT,PHASE) = ALL',
        'BEGIN BULK',
        ]
    deck = CaseControlDeck(lines)
