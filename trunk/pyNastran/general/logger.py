import os
import sys
import inspect

def frame(n):
    return inspect.currentframe(n+4) # jump 4 levels down to get out of the logger code

def lineNum():
    return frame().f_lineno

def fileNamE():
    return os.path.basename(frame().f_globals['__file__'])

def fileNameLineNum(n=0):
    f = frame(n)
    lineNum = f.f_lineno
    fname = os.path.basename(f.f_globals['__file__'])
    return(fname,lineNum)
    
class debugLogger(object):
    def __init__(self):
        pass

    #def frame(self): # good...
    #    return inspect.currentframe()

    def frame(self):
       #print help(inspect)
       return inspect.currentframe(4) # jump 4 levels down to get out of the logger code

    def lineno(self):
        """Returns the current line number in our program."""
        return self.frame().f_lineno

    def fname(self):
        """Returns the current file in our program."""
        return os.path.basename(self.frame().f_globals['__file__'])

    def funcName(self):
        """Returns the current function name in our program."""
        #print inspect.currentframe().f_back.f_back.f_locals['self']
        #print self.frame().f_globals.keys()
        #print self.frame().f_locals.keys()

        #return os.path.basename(self.frame().f_globals['__name__'])
        #return self.frame().__file__
        pass

    def properties(self):
        n = self.lineno()
        fname = self.fname()
        #funcName = self.funcName()
        funcName = ''
        return (n,fname,funcName)

    def fixMessage(self,msg,n=54):
        msg = str(msg)
        #print "msg = |%s| type=%s" %(msg,type(msg))
        lines = msg.split('\n')
        lines2 = [lines[0]]

        spaces = ' '*n
        #print "lines = ",lines
        for line in lines[1:]:
            lines2.append(spaces+line+'\n')
        msg2 = ''.join(lines2)
        return msg2

    def debug(self,msg):
        (n,fname,funcName) = self.properties()
        msg = self.fixMessage(msg)
        sys.stdout.write('DEBUG:   fname=%-25s lineNo=%-4s   %s\n' %(fname,n,msg))

    def info(self,msg):
        (n,fname,funcName) = self.properties()
        sys.stdout.write('INFO:    fname=%-25s lineNo=%-4s   %s\n' %(fname,n,msg))

    def warning(self,msg):
        (n,fname,funcName) = self.properties()
        sys.stderr.write('WARNING: fname=%-25s lineNo=%-4s   %s\n' %(fname,n,msg))

    def error(self,msg):
        (n,fname,funcName) = self.properties()
        print 'ERROR:   fname=%-25s lineNo=%-4s   %s\n' %(fname,n,msg)

    def critical(self,msg):
        (n,fname,funcName) = self.properties()
        sys.stderr.write('CRITICAL fname=%-25s lineNo=%-4s   %s\n' %(fname,n,msg))
    ###

class infoLogger(debugLogger):
    def debug(self,msg):
        pass

class dummyLogger(object):
    def __init__(self):
        pass

    def startLog(self,level):
        if level=='debug':
            return debugLogger()
        elif level=='info':
            return infoLogger()
        else:
            raise Exception("invalid logger:  debug, info ONLY!")

def getLogger(log=None,level='debug'):
    """
    This function is useful as it will instantiate a dummy logger object if log=None
    log may be a logger object or None
    pyFile is the python file that the code was called from (unused)
    """
    if log==None:
        log = dummyLogger()
        logger = log.startLog(level)
    else:
        logger = log
    return logger

if __name__=='__main__':
    # how to use a dummy logger
    logger = dummyLogger()
    log = logger.startLog('debug') # or info
    log.debug('test message')

