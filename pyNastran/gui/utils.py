import urllib2
import pyNastran

def check_for_newer_version(window=None, pop_msg=False):
    current_version = pyNastran.__version__
    target_url = 'https://raw.githubusercontent.com/SteveDoyle2/pyNastran/master/README.md'
    try:
        # it's a file like object and works just like a file
        data = urllib2.urlopen(target_url)
    except urllib2.URLError:
        return
    for btye_line in data: # files are iterable
        line_lower = btye_line.lower().decode('utf-8')
        if 'has been released' in line_lower:
            sline = line_lower.split()
            version = [slot for slot in sline if slot.startswith('v')][0][1:]
            break

    is_dev = False
    if 'dev' in current_version:
        is_dev = True

    major, minor, rev = current_version.split('_')[0].split('.')
    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_current_version = (major, minor, rev)

    major, minor, rev = version.split('_')[0].split('.')
    major = int(major)
    minor = int(minor)
    rev = int(rev)
    tuple_latest_version = (major, minor, rev)

    if tuple_current_version < tuple_latest_version or (is_dev and tuple_current_version <= tuple_latest_version):
        print('pyNastran %s is now availible; current=%s' % (version, pyNastran.__version__))
