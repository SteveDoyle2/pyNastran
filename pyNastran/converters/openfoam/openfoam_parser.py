from collections import OrderedDict


def remove_c_comments(lines):
    # remove // comments
    lines2 = []
    for line in lines:
        line2 = line.rstrip().split('//')[0]
        if line2.strip():
            lines2.append(line2)

    lines3 = []
    for line in lines2:
        line3 = line.split('|')[0].split('/*')[0].split('\\*')[0]
        if line3.strip():
            lines3.append(line3)
    return lines3

def convert_to_dict(self, lines, debug=True):
    """
    should use regex...

    FoamFile
    {
        version     0.0508;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }

    -> data = {
        'FoamFile' : {
            'version' : '0.0508',
            'format'  : 'ascii',
            'class'   : 'dictionary',
            'object'  : 'blockMeshDict',
        }
    }
    """
    #lines = ' '.join(lines).replace('\n', ' ')
    #lines = lines.replace('  ', ' ')
    i = 0
    data = OrderedDict()

    active_key = None
    active_data = data
    active_keys = []
    step = ''

    while i < len(lines):
        line = lines[i].strip()
        sline = line.split()
        if debug:
            self.log.debug("line = %r" % line)
        #if active_keys:
            #assert active_keys[0] in ['FoamFile', 'convertToMeters', 'vertices',
                                      #'blocks', 'edges', 'boundary', 'mergeMatchPairs'], active_keys[0]
        #for key in data.keys():
            #assert key in ['FoamFile', 'convertToMeters', 'vertices',
                           #'blocks', 'edges', 'boundary', 'mergeMatchPairs'], data.keys()

        if line == '{':
            if debug:
                self.log.debug('%s*A1a %r' % (step, line))
            #active_keys.append(active_key)
            active_data[active_key] = OrderedDict()
            active_data = active_data[active_key]
            if debug:
                self.log.debug('   %sstepping; active_key=%r' % (step, active_key))
            step += '  '
        elif line == '(':
            if debug:
                self.log.debug('%s*A1b %r' % (step, line))
            #active_keys.append(active_key)
            active_data[active_key] = OrderedDict()
            active_data = active_data[active_key]
            if debug:
                self.log.debug('   %sstepping; active_key=%r' % (step, active_key))
            step += '  '
            inode = 0
        elif '}' in line[0] or ');' in line:
            if debug:
                self.log.debug('%s*A2a %r' % (step, line))
                self.log.debug('   %sunstepping; active_key=%r' % (step, active_key))
                self.log.debug('   %sactive_keys1=%r' % (step, active_keys))
            key_list = active_keys[:-1]
            if debug:
                self.log.debug('key_list = %s' % key_list)

            active_data = data
            if debug:
                self.log.debug('active_data.keys() = %s' % active_data.keys())
            for key in key_list:
                if debug:
                    self.log.debug('  key=%r active_data.keys()=%s' % (key, active_data.keys()))
                active_data = active_data[key]

            active_keys.pop()
            if len(active_keys):
                active_key = active_keys[-1]
            else:
                active_key = None
            if debug:
                self.log.debug('   %sactive_keys2=%r' % (step, active_keys))
            step = step[:-2]
        elif line == ')':
            step = step[:-2]
            if debug:
                self.log.debug('%s*A2b %r' % (step, line))
                self.log.debug('   %sun); active_key=%r' % (step, active_key))
        else:
            if len(sline) == 1:  # FoamFile
                if debug:
                    self.log.debug('%s*B %r' % (step, line))
                active_key = sline[0]
                if debug:
                    self.log.debug(active_data.keys())
                    self.log.debug(active_key)
                if isinstance(active_data, list):
                    active_data.append()
                active_data[active_key] = None
                active_keys.append(active_key)
                #active_data = active_data[active_key]
                step += '  '
                if debug:
                    self.log.debug('%ssingle - adding active_key=%r' % (step, active_key))
                    self.log.debug('%sactive_keys=%s' % (step, active_keys))
                assert active_key is not None
                assert active_key != '};', active_key
                assert active_key != ');', active_key
                assert ';' not in active_key, active_key
                #active_data =
            else:
                #if debug:
                    #print('   %sdataline = %r  active_keys=%s active_key=%s' % (
                        #step, line, active_keys, active_key))
                if line.endswith(';'): # single-value
                    if debug:
                        self.log.debug('*C1 %s' % line)
                    #assert active_key is not None
                    key, value = sline
                    if active_key is None:
                        active_data[key] = value[:-1]
                    else:
                        active_data[key] = value[:-1]
                    if debug:
                        self.log.debug('  %sadding %r -> %s' % (step, key, value))
                else: # multi-value

                    #if debug:
                        #print '%s*C2 %r; active_keys=%s; active_key=%s' % (
                            #step, line, active_keys, active_key, )
                    value = line
                    #if debug:
                        #print('active_data = %s' % active_data)
                        #self.log.debug('active_keys = %s' % active_keys)
                    #if isinstance(active_data, dict):
                        #active_data = []
                        #inode = 0
                    active_data[inode] = value
                    inode += 1
        i += 1
    return data


def write_dict(openfoam_dict, nbase=0, baseword='name'):
    msg = ''
    space1 = ' ' * nbase
    space2 = ' ' * (nbase + 4)

    msg += '%s{\n' % baseword
    msg += '%s#keys (%s)=%s\n' % (space2, type(openfoam_dict), openfoam_dict.keys())
    for key, value in sorted(openfoam_dict.items()):
        if isinstance(value, basestring):
            value = value.strip()

        if isinstance(key, basestring):
            msg += '%s%s : ' % (space2, key)
        else:
            msg += '%s%s (%s) : ' % (space2, key, type(value))
        if isinstance(value, dict):
            value = write_dict(value, nbase + 4, baseword='')
        msg += '%s\n' % (str(value).strip())
    msg += '%s{\n' % space1
    return msg

class FoamFile:
    def __init__(self, filename, log=None):
        self.filename = filename
        self.log = log

    def read_foam_file(self):
        with open(self.filename, 'r') as foam_file:
            lines = foam_file.readlines()

        lines = remove_c_comments(lines)
        #for line in lines:
            #print line
        return lines
