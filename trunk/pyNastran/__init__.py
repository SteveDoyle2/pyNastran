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
    revision = '???'

__author__  = "Steven Doyle"
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF & OP2 File reader/editor/writer/viewer'
__copyright__   = 'Copyright 2011-2012, The pyNastran Project'
__license__     = 'LGPLv3'
__releaseDate__ = '2011/2/??'
__version__     = '0.3.0.r%s' %(revision)
__website__     = 'http://code.google.com/p/pynastran/'
