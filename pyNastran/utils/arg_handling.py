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
            if 'unrecognized arguments' in message or 'ignored explicit argument' in message:
                mymsg = message
            elif message.strip() == pyNastran.__version__:
                mymsg = message
            file.write(mymsg)

    parser._print_message = _print_message
