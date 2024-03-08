"""
defines:
 - version_latest, version_current, is_newer = check_for_newer_version(version_current=None)
"""
import urllib.request
import urllib.error
from typing import Optional

import pyNastran


def get_data_from_website(target_url: str) -> tuple[Optional[list[str]], bool]:
    """downloads the byte data from the website"""
    is_failed = True
    data = ''
    try:
        # it's a file like object and works just like a file
        data_bytes = urllib.request.urlopen(target_url)
        is_failed = False
    except (urllib.error.HTTPError, urllib.error.URLError):  #  forbidden, page not found
        #print(f'error...target_url={target_url}')
        #import os
        #print(os.environ.keys())
        #raise
        pass
    #except Exception: #  urllib2.URLError  # e.g., timeout
        #print(help(urllib))
        #raise

    lines: list[str] = []
    if not is_failed:
        data = data_bytes.read().decode('utf-8')
        lines = data.split('\n')  # basically data.readlines()
    return lines, is_failed

def split_version(version: str, msg: str) -> tuple[int, int, int]:
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

def get_latest_version_from_data(data: str) -> Optional[str]:
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
                            quiet: bool=False) -> tuple[Optional[str], Optional[str], bool]:
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
    assert version_current is not None, version_current

    is_failed, version_latest = _get_latest_version(version_latest)
    if is_failed:
        return None, None, is_newer
    is_newer = _check_if_version_is_newer(version_latest, version_current, quiet)
    return version_latest, version_current, is_newer

def _get_latest_version(version_latest: Optional[str]) -> tuple[bool, str]:
    is_failed = False
    if version_latest is None:
        target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/latest.txt'
        data, is_failed = get_data_from_website(target_url)
        if is_failed:
            return is_failed, version_latest
        version_latest = get_latest_version_from_data(data)  # type: ignore
        if version_latest is None:
            raise RuntimeError("can't parse website")
    return is_failed, version_latest

def _check_if_version_is_newer(version_latest: str,
                               version_current: str,
                               quiet: bool) -> bool:
    """
    Not 100% on this list, but this is the general idea...

    Current   Release   Dev?   -> is_newer   Result         Description
    ========= ========= ====      ========   ========       =================
    (1, 0, 0) (1, 3, 2) True   -> True       download       typical user
    (1, 4, 0) (1, 3, 2) True   -> False      don't download typical dev
    (1, 4, 0) (1, 4, 0) False  -> True       download       intermediate gui release
    (1, 4, 0) (1, 3, 2) False  -> False      dont download  release candidate

    """
    is_newer = False
    is_dev = 'dev' in version_current
    tuple_current_version = split_version(version_current, 'current')
    tuple_latest_version = split_version(version_latest, 'latest')

    #print('tuple_latest_version = %s' % str(tuple_latest_version))  # (0,7,2)
    #print('tuple_current_version = %s' % str(tuple_current_version))  # (0,8,0)

    is_self_newer = tuple_current_version > tuple_latest_version
    is_newer_release_version = tuple_current_version < tuple_latest_version
    is_newer_dev_version = (not is_dev) and (tuple_current_version <= tuple_latest_version)

    #print(f"tuple_current_version={tuple_current_version} tuple_latest_version={tuple_latest_version}")
    #print(f"is_dev={is_dev} is_self_newer={is_self_newer} is_newer_release_version={is_newer_release_version} is_newer_dev_version={is_newer_dev_version}")
    if is_self_newer:
        pass
    elif is_newer_release_version or is_newer_dev_version:
        msg = f'pyNastran {version_latest} is now available; current={version_current}'

        if not quiet:  # pragma: no cover
            print(msg)
        is_newer = True
    return is_newer

if __name__ == '__main__':  # pragma: no cover
    check_for_newer_version()
