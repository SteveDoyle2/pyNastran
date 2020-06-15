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

def split_version(version: str, msg: str) -> Tuple[int, int, int]:
    try:
        major, minor, rev = version.split('+')[0].split('.')
    except ValueError:
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
        version_latest = line.lower()
        break
    if version_latest is None:
        return version_latest
    return version_latest

def check_for_newer_version(version_current: Optional[str]=None,
                            version_latest: Optional[str]=None,
                            quiet=False) -> Tuple[Optional[str], Optional[str], bool]:
    """
    Checks to see if a newer version of pyNastran has been released.
    Only checks this for the GUI.

    Looks in latest.txt for:
        0.7.2
    and compares that to the current version.

    Returns
    -------
    version_latest : str
        the latest version found
        '1.3.2'
    version_current : str
        the current version (may be a dev version)
        '1.3.2'
        '1.4.0+dev.5378fd363'
    is_newer : bool
        is there a newer version

    """
    is_newer = False
    if version_current is None:
        version_current = pyNastran.__version__

    if version_latest is None:
        target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/latest.txt'
        data, is_failed = get_data_from_website(target_url)
        if is_failed:
            return None, None, is_newer
        version_latest = get_latest_version_from_data(data)  # type: ignore
        if version_latest is None:
            raise RuntimeError("can't parse website")
            #return None, None, is_newer

    is_dev = 'dev' in version_current
    tuple_current_version = split_version(version_current, 'current')
    tuple_latest_version = split_version(version_latest, 'latest')

    #print('tuple_latest_version = %s' % str(tuple_latest_version))  # (0,7,2)
    #print('tuple_current_version = %s' % str(tuple_current_version))  # (0,8,0)

    # (1, 0, 0) (1, 3, 2) True
    is_self_newer = tuple_current_version > tuple_latest_version
    is_newer_release_version = tuple_current_version < tuple_latest_version
    is_newer_dev_version = is_dev and (tuple_current_version <= tuple_latest_version)

    if is_self_newer:
        pass
    elif is_newer_release_version or is_newer_dev_version:
        msg = 'pyNastran %s is now availible; current=%s' % (version_latest, version_current)

        if not quiet:  # pragma: no cover
            print(msg)
        is_newer = True
    return version_latest, version_current, is_newer

if __name__ == '__main__':  # pragma: no cover
    check_for_newer_version()
