from __future__ import print_function
from six.moves import urllib

import pyNastran


def check_for_newer_version():
    """
    Checks to see if a newer version of pyNastran has been released.
    Only checks this for the GUI.

    Looks for:
        ## pyNastran v0.7.2 has been Released (4/25/2015)

    Specifically, it finds, 'has been released'
       then takes the the part that:
         - starts with 'v',
         - strips the 'v'
         - makes a version tuple:
           - (0,7,2)
       and compares that to the current version
    """
    is_newer = False
    version_current = pyNastran.__version__
    target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/README.md'
    try:
        # it's a file like object and works just like a file
        data = urllib.request.urlopen(target_url)
    except (urllib.error.HTTPError, urllib.error.URLError):
        #print(help(urllib))
        #raise
        return None, None, False
    for btye_line in data: # files are iterable
        line_lower = btye_line.lower().decode('utf-8')
        if 'has been released' in line_lower:
            sline = line_lower.split()
            version_latest = [slot for slot in sline if slot.startswith('v')][0][1:]
            break

    is_dev = False
    if 'dev' in version_current:
        is_dev = True

    try:
        major, minor, rev = version_current.split('+')[0].split('.')
    except ValueError:
        print('sline = %s' % sline)
        print('version_current = %s' % version_current)
        raise
    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_current_version = (major, minor, rev)

    try:
        major, minor, rev = version_latest.split('_')[0].split('.')
    except ValueError:
        print('sline = %s' % sline)
        print('version_latest = %s' % version_latest)
        raise

    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_latest_version = (major, minor, rev)
    #print('tuple_latest_version = %s' % str(tuple_latest_version))  # (0,7,2)
    #print('tuple_current_version = %s' % str(tuple_current_version))  # (0,8,0)

    #is_newer = True
    if (tuple_current_version < tuple_latest_version or
            (is_dev and tuple_current_version <= tuple_latest_version)):
        print('pyNastran %s is now availible; current=%s' % (version_latest, version_current))
        is_newer = True
    #print('*pyNastran %s is now availible; current=%s' % (version_latest, version_current))
    return version_latest, version_current, is_newer

if __name__ == '__main__':  # pragma: no cover
    check_for_newer_version()
