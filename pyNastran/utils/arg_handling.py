"""interface to argparse"""
import sys
import pyNastran

def swap_key(mydict, key_orig, key_new, flip=False):
    """replaces a key in a dictionary"""
    if flip:
        mydict[key_new] = not mydict[key_orig]
    else:
        mydict[key_new] = mydict[key_orig]
    del mydict[key_orig]

def argparse_to_dict(args):
    """converts the argparse output into a dictionary"""
    argdict = {}
    for name, value in args._get_args():
        argdict[name] = value
    for name, value in args._get_kwargs():
        argdict[name] = value
    return argdict

def update_message(parser, usage, arg_msg, examples):
    """overwrites the default argparse help message"""
    def _print_message(message, file=None):
        """overwrites the argparse print to get a better help message"""
        mymsg = usage + arg_msg + examples
        if message:
            if file is None:
                file = sys.stderr

            #file.write(message)
            #return
            is_error = (
                'unrecognized arguments' in message or         # pyNastranGUI --at
                'ignored explicit argument' in message or      # ???
                'ambiguous option' in message or               # pyNastranGUI --q  (--qt, --quiet)
                'following arguments are required' in message) # test_bdf <no_args>
            if is_error:
                mymsg = message
            elif message.strip() == pyNastran.__version__:
                mymsg = message
            #else:
                #mymsg = message
            file.write(mymsg)

    parser._print_message = _print_message
