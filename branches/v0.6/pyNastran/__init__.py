dev = 0
if dev == 1:
    import os
    import pysvn
    import inspect

    print(inspect.getfile(inspect.currentframe())) # script filename (usually with path)
    curDir = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory

    client = pysvn.Client()
    #print(dir(client))
    client.update(curDir)
    revision = client.info(curDir).revision.number
    print("revision = %s" %(revision))
    changes = client.get_changelist(curDir)
    print("changes = %s" %(changes))
    #client.update('./examples/pysvn')
elif dev == 2:
    revision = 'dev'

__author__  = 'Steven Doyle, Al Danials, Marcin Gasiorek'
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__copyright__   = 'Copyright 2011-2013, pyNastran; %s' % __author__
__license__     = 'LGPLv3'
__releaseDate__ = '2014/6/xx'
__releaseDate2__ = 'JUNE xx, 2014'
__version__     = '0.6.2'
__website__     = 'http://code.google.com/p/pynastran/'

isRelease = True  ## turns on skipping of tables that aren't supported