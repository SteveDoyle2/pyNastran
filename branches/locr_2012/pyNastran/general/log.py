import sys
import platform
import os
import inspect
import logging
from functools import partial

def make_log(display=False):
    """
    Creates 'pyNastran.log' file with information about working environment,
    such as Python version, platform, architecture, etc. Useful for debugging.
    @param do not only create file but also print log information
    """
    smsg = [("sys.version", sys.version), ("sys.version_info", sys.version_info)]
    pmsg = ["machine", "platform", "processor", "architecture","python_branch", 
           "python_revision", "win32_ver", "version", "uname", "system",
           "python_build", "python_compiler", "python_implementation", "system",
           "mac_ver", "linux_distribution", "libc_ver"]
    
    fmt = "%-{0}s = %s\n".format(max(map(len, pmsg + [j[0] for j in smsg])))
    msg = "".join([fmt % (i, str(j).replace("\n", "; ")) for (i,j) in smsg])
    msg += "".join([fmt % (i, str(getattr(platform, i)())) for i in pmsg])
    if display:
        print(msg)
    
    with open('pyNastran.log', 'w') as fil:
        fil.write(msg)



class debugLogger(object):

    def properties(self):
        """Return tuple: line number and filename"""
        f =  inspect.currentframe(3)  # jump 2 levels down to get out of the logger code
        return (f.f_lineno, os.path.basename(f.f_globals['__file__']))

    def debug(self, msg):
        lines = str(msg).split('\n')
        self.msg_typ('DEBUG',''.join([lines[0]] + [' ' * 54 + line + '\n' 
                                                   for line in lines[1:]]))
        
    def msg_typ(self,typ,msg):
        n, fn = self.properties()
        sys.stdout.write('%s:    fname=%-25s lineNo=%-4s   %s\n' % (typ, fn, n, msg))

    def info(self, msg):
        self.msg_typ("INFO", msg)
    
    def warning(self, msg):
        self.msg_typ("WARNING", msg)
    
    def error(self, msg):
        self.msg_typ("ERROR", msg)
    
    def critical(self, msg):
        self.msg_typ("CRITICAL", msg)
    

class infoLogger(debugLogger):
    def debug(self, msg):
        pass


class dummyLogger(object):
    
    def startLog(self, level):
        if level == 'debug':
            return debugLogger()
        elif level == 'info':
            return infoLogger()
        raise RuntimeError("invalid logger:  debug, info ONLY!")


def get_logger(log=None, level='debug'):
    """
    This function is useful as it will instantiate a dummy logger object if log=None
    log may be a logger object or None
    pyFile is the python file that the code was called from (unused)
    """
    return dummyLogger().startLog(level) if log is None else log
    


def buildDummyLogger2(level='debug'):
    try:
        version = sys.version_info
        fname = 'pyNastran.py%s%s.log' % (version.major, version.minor)
        logPath = os.path.join(os.getcwd(), fname)
        if os.path.exists(logPath):
            os.remove(logPath)

        # create logger with 'pyNastran'
        logger = logging.getLogger('pyNastran')
        logger.setLevel(logging.DEBUG)
        # create file handler which logs even debug messages
        fh = logging.FileHandler(logPath)
        if level == 'debug':
            fh.setLevel(logging.DEBUG)
        elif level == 'info':
            fh.setLevel(logging.INFO)
        else:
            raise RuntimeError("invalid logger:  debug, info ONLY!")
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(levelname)-8s:  %(filename)-25s linenum=%(lineno)-4s %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)
        logger.addHandler(ch)

        logger.info('logger is initialized---')
    except Exception:
        logger = logging.getLogger('pyNastran')

    return logger

if __name__ == '__main__':
    # how to use a dummy logger
    logger = dummyLogger()
    log = logger.startLog('debug')  # or info
    log.debug('test message')
    log.warning('asf')
    log.error('errorss')