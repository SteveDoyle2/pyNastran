import os
from pyNastran.utils import print_bad_path
#versions = ['26', '27', '31', '32', '33']
#versions = ['27', '32']
versions = ['27']

for version in versions:
    v0 = version[0]
    egginfo_path = r'C:\Python%s\Scripts\egginfo.exe' % version
    python_path = r'C:\Python%s\python.exe' % version
    
    if v0 == '2':
        pkg_path = r'D:\work\pyNastran_py2x\pyNastran'
    elif v0 == '3':
        pkg_path = r'D:\work\pyNastran_py3x\pyNastran'
    else:
        raise NotImplementedError(version)

    test_dir = pkg_path + r'\pyNastran'
    egginfo = '%s pyNastran' % egginfo_path
    #print "egginfo =",egginfo
    test_path = r'%s\all_tests.py' % test_dir

    #print print_bad_path(egginfo_path)
    print print_bad_path(test_path)

    tests = r'%s %s' % (python_path, test_path)
    print "tests = %r" % tests

    code1 = os.system(egginfo)
    assert code1 == 0
    code2 = os.system(tests)
    assert code2 == 0

