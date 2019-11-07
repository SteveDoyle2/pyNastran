"""
defines:
 - version_latest, version_current, is_newer = check_for_newer_version(version_current=None)
"""
import urllib.request
import urllib.error
from typing import List, Tuple, Optional

import pyNastran


def get_data_from_website(target_url: str) -> Tuple[Optional[List[str]], bool]:
    """downloads the byte data from the website"""
    is_failed = True
    data = ''
    try:
        # it's a file like object and works just like a file
        data_bytes = urllib.request.urlopen(target_url)
        is_failed = False
    except (urllib.error.HTTPError, urllib.error.URLError):  #  forbidden, page not found
        pass
    #except: #  urllib2.URLError  # e.g., timeout
        #print(help(urllib))
        #raise

    lines = [] # type: List[str]
    if not is_failed:
        data = data_bytes.read().decode('utf-8')
        lines = data.split('\n')  # basically data.readlines()
    return lines, is_failed

def split_version(version: str, sline: List[str], msg: str) -> Tuple[int, int, int]:
    try:
        major, minor, rev = version.split('+')[0].split('.')
    except ValueError:
        print('sline = %s' % sline)
        print(f'version_{msg} = {version}')
        raise
    imajor = int(major)
    iminor = int(minor)
    irev = int(rev)
    tuple_version = (imajor, iminor, irev)
    return tuple_version

def get_latest_version_from_data(data: str) -> Tuple[Optional[str], List[str]]:
    """finds the latest released version"""
    version_latest = None
    for line in data: # files are iterable
        line_lower = line.lower()
        #print(line_lower.rstrip())
        if 'has been released' in line_lower:
            sline = line_lower.split()
            version_latest = [slot for slot in sline if slot.startswith('v')][0][1:]
            break
    if version_latest is None:
        return version_latest, []
    return version_latest, sline

def check_for_newer_version(version_current: Optional[str]=None,
                            quiet=False) -> Tuple[Optional[str], Optional[str], bool]:
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
    if version_current is None:
        version_current = pyNastran.__version__
    target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/README.md'
    data, is_failed = get_data_from_website(target_url)
    if is_failed:
        return None, None, is_newer

    version_latest, sline = get_latest_version_from_data(data)  # type: ignore
    if version_latest is None:
        raise RuntimeError("can't parse website")
        #return None, None, is_newer

    is_dev = 'dev' in version_current
    tuple_current_version = split_version(version_current, sline, 'current')
    tuple_latest_version = split_version(version_latest, sline, 'latest')

    #print('tuple_latest_version = %s' % str(tuple_latest_version))  # (0,7,2)
    #print('tuple_current_version = %s' % str(tuple_current_version))  # (0,8,0)

    if (tuple_current_version < tuple_latest_version or
            (is_dev and tuple_current_version == tuple_latest_version)):
        msg = 'pyNastran %s is now availible; current=%s' % (version_latest, version_current)

        if not quiet:  # pragma: no cover
            print(msg)
        is_newer = True
    return version_latest, version_current, is_newer

if __name__ == '__main__':  # pragma: no cover
    check_for_newer_version()
