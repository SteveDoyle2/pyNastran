from __future__ import absolute_import, print_function
from functools import reduce
import re
import sys
from applyLicense import get_folders_files
   

def fix_files(fils):
    
    def fix_string(code_string):
        """
        Function that transforms given Python 2 program code_string into (hopefully)
        correct Python 3 program.
        """
        _multisub = lambda s, lst: reduce(lambda _s, l: _s.replace(*l), lst, s)
        
        res = []
        for line in code_string.splitlines():
            line = _multisub(line, (('from StringIO import StringIO','from io import StringIO'),
                ('from itertools import izip, count','from itertools import count'),
                ('from itertools import izip', ''), ('izip(','zip('), 
                ('xrange(','range('), ('.iteritems()','.items()'), 
                ('.itervalues()','.values()'), ('.iterkeys()','.keys()'),
                ('unicode)', 'str)'),("unpack(b'","unpack('")))
        
            # look for .keys() and values() method calls that are not inside
            # set() or sorted(), list() 
            keyval = re.search('[\w\.\[\]]+\.(keys|values)\(\)',line)
            if keyval:
                kv = keyval.group()
                if reduce(lambda r, s: r and ("%s(%s)" % (s, kv) not in line), 
                              ("set", "sorted", "list"), True):
                    line = line.replace(kv, "list(%s)" % (kv))
            # remove var = bytes(var) lines
            if re.search("([\w\.\[\]]+) *= *bytes\(\1\)", line):
                line= ''
                
            #file handling - add utf8 coding option, do not change io.open -> would brake is_binary in utils
            line = re.sub(r'((?<!\.)open\(.*,.*?)([\'"])(r|w)b{0,1}\2.*?(\))', 
                          r"\1'\3', encoding='utf-8'\4",line) 
        
            if "data = b''" in line:
                pass
            elif " b'" in line and 'data' not in line:
                line = line.replace(" b'"," '")
                
            #remove unicode mark 'u' at th beginig of strings
            line = re.sub(r'u((["\']).*?(?<!\\)(\\\\)*\2)', r'\1', line) 
            res.append(line)
        return "\n".join(res)

    # Fix every python file in the project source directory
    for fil in fils:
        with open(fil, 'r') as f:
            txt = f.read()
        with open(fil, 'w') as f:
            f.write(fix_string(txt))
   
if __name__ == "__main__":
    files = [i for i in  get_folders_files('../pyNastran')[1] if i.endswith('.py')]
    print("Applying changess to ", len(files), " files... ", end = "")
    sys.stdout.flush()
    fix_files(files)
    print("done")
