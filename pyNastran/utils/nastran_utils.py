import subprocess
from six import iteritems

def run_nastran(fname, keywords=None):
    """
    Call a nastran subprocess with the given filename

    Parameters
    -----------
    fname : string
        Filename of the Nastran .bdf file
    keywords : dict/list of strings, optional
        Default keywords are `'scr=yes'`, `'bat=no'`, `'old=no'`, and `'news=no'`
    """
    if keywords is None:
        keywords_list = ['scr=yes', 'bat=no', 'old=no','news=no'] # 'mem=1024mb',
    else:
        if isinstance(keywords, (list, tuple)):
            keywords_list = keywords
        else:
            keywords_list = []
            for keyword, value in iteritems(keywords):
                if value is None:
                    continue
                keywords_list.append('%s=%s' % (keyword, value))

    call_args = ['nastran', fname] + keywords_list
    return subprocess.call(call_args)

