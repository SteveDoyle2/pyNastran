## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
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
__releaseDate__ = '2013/6/12'
__releaseDate2__ = 'JUNE 12, 2013'
__version__     = '0.6.1'
__website__     = 'http://code.google.com/p/pynastran/'

isRelease = True  ## turns on skipping of tables that aren't supported