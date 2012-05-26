dev = False
if dev:
    import os
    import pysvn
    
    import inspect, os
    print inspect.getfile(inspect.currentframe()) # script filename (usually with path)
    curDir = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
    
    client = pysvn.Client()
    #print dir(client)
    client.update(curDir)
    revision = client.info(curDir).revision.number
    print "revision = ",revision
    changes = client.get_changelist(curDir)
    print "changes = ",changes
    #client.update('./examples/pysvn')
else:
    revision = 'dev790'

__author__  = 'Steven Doyle, Al Danials'
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__copyright__   = 'Copyright 2011-2012, pyNastran; %s' %(__author__)
__license__     = 'LGPLv3'
__releaseDate__ = '2012/5/25'
__releaseDate2__ = 'MAY 25, 2012'
__version__     = '0.5.r%s' %(revision)
__website__     = 'http://code.google.com/p/pynastran/'
