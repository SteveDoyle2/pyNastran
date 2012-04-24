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
from __future__ import print_function
import sys
from pyNastran.bdf.cards.constraints import constraintObject2

class XrefMesh(object):
    def __init__(self):
        pass

    def crossReference(self,xref=True):
        """
        Links up all the cards to the cards they reference
        """
        if xref:
            self.log.info("Cross Referencing...")
            #for key,e in self.elements.items():
            #    print(e)

            self.crossReference_Nodes()
            self.crossReference_Coordinates()

            self.crossReference_Elements()
            self.crossReference_Properties()
            self.crossReference_Materials()

            self.crossReference_Aero()
            self.crossReference_Constraints()
            #self.crossReference_Loads()
            #self.caseControlDeck.crossReference(self)
        ###
    
    def crossReference_Constraints(self):
        #self.spcObject.crossReference(self)  # enable to output SPCs
        #self.mpcObject.crossReference(self)  # enable to output MPCs
        pass
        
        #self.spcObject2 = constraintObject2()
        for key,spcadd in sorted(self.spcadds.items()):
            self.spcObject2.Add(spcadd)

        for key,spcs in sorted(self.spcs.items()):
            #print spcs
            for spc in spcs:
                self.spcObject2.append(spc)
        ###

        #self.mpcObject2 = constraintObject2()
        for key,mpcadd in sorted(self.mpcadds.items()):
            self.mpcObject2.Add(mpcadd)

        for key,mpcs in sorted(self.mpcs.items()):
            #print spcs
            for mpc in mpcs:
                self.mpcObject2.append(mpc)
        ###
        #self.spcObject.crossReference(self)


    def crossReference_Coordinates(self):
        """
        Links up all the coordinate cards to other coordinate cards and nodes
        """
        for cid,c in self.coords.items(): # CORD2x: links the rid to coordinate systems
            c.crossReference(self)        # CORD1x: links g1,g2,g3 to grid points
        ###
        for cid,c in self.coords.items(): # CORD1x: Since the grid points were already referenced,
            c.resolveCid()                # we can now resolve the coordinate systems.
        ###                               # We couldnt do it in the previous step b/c
                                          # the grid's coordinate system might have been
                                          # unresolved

    def crossReference_Aero(self):
        """
        Links up all the aero cards
        """
        for ID,caero in self.caeros.items():
            caero.crossReference(self)
        ###
        for ID,spline in self.splines.items():
            spline.crossReference(self)
        ###

    def crossReference_Nodes(self):
        """
        Links the nodes to coordinate systems
        """
        gridSet = self.gridSet
        for nid,n in self.nodes.items():
            try:
                n.crossReference(self,gridSet)
            except:
                sys.stderr.write("Couldn't cross reference GRID.  Are all Coordinate Systemes supported?\n%s" %(str(n)))
                raise
        ###
        if self.spoints:
            self.spointi = self.spoints.createSPOINTi()
        ###

    def crossReference_Elements(self):
        """
        Links the elements to nodes, properties (and materials depending on the card)
        """
        for eid,e in self.elements.items():
            try:
                e.crossReference(self)
            except:
                sys.stderr.write("Couldn't cross reference Element.  Are all Properties supported?\n%s" %(str(e)))
                raise
        ###

    def crossReference_Properties(self):
        """
        Links the properties to materials
        """
        for pid,p in self.properties.items():
            #print p
            try:
                p.crossReference(self)
            except:
                sys.stderr.write("Couldn't cross reference Property.  Are all Materials supported?\n%s" %(str(p)))
                raise
        ###

    def crossReference_Materials(self):
        """
        Links the materials to materials (e.g. MAT1, CREEP)
        often this is a pass statement
        """
        for mid,m in self.materials.items(): # MAT1
            try:
                m.crossReference(self)
            except:
                sys.stderr.write("Couldn't cross reference Material\n%s" %(str(m)))
                raise
        ###
        for mid,m in self.materialDeps.items(): # CREEP - depends on MAT1
            try:
                m.crossReference(self)
            except:
                sys.stderr.write("Couldn't cross reference Material\n%s" %(str(m)))
                raise
        ###

    def crossReference_Loads(self):
        """
        Links the loads to nodes, coordinate systems, and other loads
        """
        for lid,sid in self.loads.items():
            self.log.debug("lid=%s sid=%s" %(lid,sid))
            for load in sid:
                try:
                    load.crossReference(self)
                except:
                    sys.stderr.write("Couldn't cross reference Load\n%s" %(str(load)))
                    raise
        ###

