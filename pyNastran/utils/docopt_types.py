from docopt import docopt

def docopt_types(doc, argv=None, help=True, version=None, options_first=False, type_defaults=None):
    """
    `docopt` creates your command-line interface based on its
    description that you pass as `doc`. Such description can contain
    --options, <positional-argument>, commands, which could be
    [optional], (required), (mutually | exclusive) or repeated...

    Parameters
    ----------
    doc : str
        Description of your command-line interface.
    argv : list of str, optional
        Argument vector to be parsed. sys.argv[1:] is used if not
        provided.
    help : bool (default: True)
        Set to False to disable automatic help on -h or --help
        options.
    version : any object
        If passed, the object will be printed if --version is in
        `argv`.
    options_first : bool (default: False)
        Set to True to require options preceed positional arguments,
        i.e. to forbid options and positional arguments intermix.
    type_defaults : dict (default:None)
        key : name of argument
        value : type, default

    Returns
    -------
    args : dict
        A dictionary, where keys are names of command-line elements
        such as e.g. "--verbose" and "<path>", and values are the
        parsed values of those elements.

    Example
    -------
    >>> from docopt import docopt
    >>> doc = '''
    Usage:
        my_program tcp <host> <port> [--timeout=<seconds>]
        my_program serial <port> [--baud=<n>] [--timeout=<seconds>]
        my_program (-h | --help | --version)

    Options:
        -h, --help  Show this screen and exit.
        --baud=<n>  Baudrate [default: 9600]
    '''
    >>> argv = ['tcp', '127.0.0.1', '80', '--timeout', '30']
    >>> docopt(doc, argv)
    {'--baud': '9600',
     '--help': False,
     '--timeout': '30',
     '--version': False,
     '<host>': '127.0.0.1',
     '<port>': '80',
     'serial': False,
     'tcp': True}

    See also
    --------
    * For video introduction see http://docopt.org
    * Full documentation is available in README.rst as well as online
      at https://github.com/docopt/docopt#readme

    """
    data = docopt(doc, argv=argv, help=help, version=version, options_first=options_first)
    if type_defaults:
        data2 = data
        for key, (type_func, default) in sorted(type_defaults.items()):
            try:
                value_in = data[key]
            except KeyError:
                msg = 'Key not found : key=%r\ndata.keys()=%s' % (key, list(data.keys()))
                raise KeyError(msg)
            if value_in is None:
                data[key] = default
                continue
            data[key] = type_func(value_in)
    return data