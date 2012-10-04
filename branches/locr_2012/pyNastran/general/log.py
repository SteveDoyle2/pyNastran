import sys
import platform
import os

def make_log(display=False):
    """
    Creates 'pyNastran.log' file with information about working environment,
    such as Python version, platform, architecture, etc. Useful for debugging.
    @param display do not only create file but also print log information
    """
    smsg = [("sys.version", sys.version), ("sys.version_info", sys.version_info)]
    pmsg = ["machine", "platform", "processor", "architecture", "python_branch", 
           "python_revision", "win32_ver", "version", "uname", "system",
           "python_build", "python_compiler", "python_implementation", "system",
           "mac_ver", "linux_distribution", "libc_ver"]
    
    fmt = "%-{0}s = %s\n".format(max(map(len, pmsg + [j[0] for j in smsg])))
    msg = "".join([fmt % (i, str(j).replace("\n", "; ")) for (i, j) in smsg])
    msg += "".join([fmt % (i, str(getattr(platform, i)())) for i in pmsg])
    if display:
        print(msg)
    
    with open('pyNastran.log', 'w') as fil:
        fil.write(msg)


class simpleLogger(object):
    """
    Simple logger object. In future might be changed to use Python logging module.
    Two levels are supported: 'debug' and 'info'. Info level discards debug
    messages, 'debug' level displays all messages.
    """
    
    def __init__(self, level='debug'):
        """
        @param level level of logging: 'info' or 'debug'
        """
        assert level in ('info','debug')
        self.level = level

    def properties(self):
        """Return tuple: line number and filename"""
        _fr =  sys._getframe(3)  # jump to get out of the logger code
        return (_fr.f_lineno, os.path.basename(_fr.f_globals['__file__']))

    def debug(self, msg):
        """
        Log DEBUG message
        @param msg message to be logged
        """
        if self.level != 'debug':
            return
        lines = str(msg).split('\n')
        self.msg_typ('DEBUG', ''.join([lines[0]] + [' ' * 54 + line + '\n' 
                                                   for line in lines[1:]]))
        
    def msg_typ(self, typ, msg):
        """
        Log message of a given type
        @param typ type of a message (e.g. INFO)
        @param msg message to be logged
        """
        n, fn = self.properties()
        sys.stdout.write('%s:    fname=%-25s lineNo=%-4s   %s\n' % (typ, fn, n, msg))

    def info(self, msg):
        """
        Log INFO message
        @param msg message to be logged
        """
        self.msg_typ("INFO", msg)
    
    def warning(self, msg):
        """
        Log WARNING message
        @param msg message to be logged
        """
        self.msg_typ("WARNING", msg)
    
    def error(self, msg):
        """
        Log ERROR message
        @param msg message to be logged
        """
        self.msg_typ("ERROR", msg)
    
    def critical(self, msg):
        """
        Log CRITICAL message
        @param msg message to be logged
        """
        self.msg_typ("CRITICAL", msg)


def get_logger(log = None, level = 'debug'):
    """
    This function is useful as it will instantiate a simpleLogger object if log=None.
    @param log a logger object or None
    @param level level of logging: 'info' or 'debug'
    """
    return simpleLogger(level) if log is None else log

if __name__ == '__main__':
    # how to use a simple logger
    for nam in ["debug", "info"]:      
        print('--- %s logger ---' % (nam))
        test_log = simpleLogger(nam)
        test_log.debug('test message')
        test_log.warning('asf')
        test_log.error('errorss')
    
