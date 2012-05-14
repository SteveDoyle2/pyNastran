import sys
from struct import pack
from numpy import array,sqrt,dot,abs,degrees

from op2_Objects import scalarObject

try:
    from pylab import xlabel,ylabel,show,grid,legend,plot,title
    import matplotlib.pyplot as plt
except ImportError:
    pass

class TableObject(scalarObject):  # displacement style table
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.gridTypes    = {}
        self.translations = {}
        self.rotations    = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def isImaginary(self):
        return False

    def addF06Data(self,data,transient):
        if transient is None:
            for line in data:
                (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
                self.gridTypes[nodeID]    = gridType
                self.translations[nodeID] = array([t1,t2,t3])
                self.rotations[nodeID]    = array([r1,r2,r3])
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.translations:
            self.updateDt(self.dataCode,dt)

        for line in data:
            (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
            self.gridTypes[nodeID]    = gridType
            self.translations[dt][nodeID] = array([t1,t2,t3])
            self.rotations[dt][nodeID]    = array([r1,r2,r3])
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],dt,self.iSubcase))
            self.dt = dt
            self.addNewTransient()
        ###

    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def addNewTransient(self,dt):
        """initializes the transient variables"""
        self.dt = dt
        self.translations[dt] = {}
        self.rotations[dt]    = {}

    #def addBinary(self,deviceCode,data):
        #(nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        #msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations

        #self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        #self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,dt,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "nodeID=%s gridType=%s v1=%s v2=%s v3=%s" %(nodeID,gridType,v1,v2,v3)
        #print msg
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements,'displacementObject - static failure'
        
        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addSort1(self,dt,out):
        #print "dt=%s out=%s" %(dt,out)
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        if dt not in self.translations:
            self.addNewTransient(dt)
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations[self.dt],'displacementObject - transient failure'

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[dt][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addSort2(self,nodeID,out):
        (dt,gridType,v1,v2,v3,v4,v5,v6) = out
        if dt not in self.translations:
            self.addNewTransient(dt)
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations[self.dt],'displacementObject - transient failure'

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[dt][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    #def writeOp2(self,block3,deviceCode=1):
        #"""
        #creates the binary data for writing the table
        #@warning hasnt been tested...
        #"""
        #msg = block3
        #for nodeID,translation in sorted(self.translations.iteritems()):
            #rotation = self.rotations[nodeID]
            #(dx,dy,dz) = translation
            #(rx,ry,rz) = rotation
            #grid = nodeID*10+deviceCode
            #msg += pack('iffffff',grid,dx,dy,dz,rx,ry,rz)
        ###
        #return msg

    #def writeOp2Transient(self,block3,deviceCode=1):
    #    """
    #    creates the binary data for writing the table
    #    @warning hasnt been tested...
    #    @warning dt slot needs to be fixed...
    #    """
    #    msg = ''
    #    for dt,displacements in sorted(self.displacements.iteritems()):
    #        XXX = 50 ## this isnt correct... @todo update dt
    #        msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
    #        #msg += '%s = %g\n' %(self.dataCode['name'],dt)
    #
    #        for nodeID,displacement in sorted(displacements.iteritems()):
    #            rotation = self.rotations[nodeID]
    #            (dx,dy,dz) = displacement
    #            (rx,ry,rz) = rotation
    #
    #            grid = nodeID*10+deviceCode
    #            msg += pack('iffffff',grid,dx,dy,dz,rx,ry,rz)
    #        ###
    #    ###
    #    return msg

    def writeHeader(self):
        #(mainHeaders,headers) = self.getHeaders()
        mainHeaders = ('nodeID','gridType')
        headers = ('T1','T2','T3','R1','R2','R3')
        msg = '%-10s %8s ' %(mainHeaders)
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        return msg

    def getAsSort1(self):
        return (self.translations,self.rotations)

    def getAsSort2(self):
        """returns translations and rotations in sort2 format"""
        translations2 = {}
        rotations2 = {}
        if self.dt is not None:
            return self.__reprTransient__()

            for dt,translations in sorted(self.translations.iteritems()):
                nodeIDs = translations.keys()
                for nodeID in nodeIDs:
                    translations2[nodeID] = {}
                    rotations2[nodeID] = {}

                for nodeID,translation in sorted(translations.iteritems()):
                    rotation = self.rotations[dt][nodeID]
                    translations2[nodeID][dt] = translation
                    rotations2[nodeID][dt]    = rotation
                ###
        else:
            for nodeID,translation in sorted(self.translations.iteritems()):
                rotation = self.rotations[nodeID]
                translations2[nodeID] = translation
                rotations2[nodeID]    = rotation
            ###
        ###
        return translations2,rotations2
        
    def plot(self,nodeList=None,resultType='Translation',coord=3,markers=None,
                  Title=None,hasLegend=True,Legend=None,xLabel=None,yLabel=None,alphaLegend=0.5):
        """
        @param nodeList a list of the node IDs to plot vs the
               independent variable (default=None; all nodes)
        @param resultType the variable to plot ('Translation','Rotation')
        @param coord the coordinate to plot (for <x,y,z>, x=0,y=1,z=2,Mag=3);
               default=Magnitude
        @param markers  a list of colors/marker shapes for each line
        @param Title title of the plot (default=the object name)
        @param hasLegend should a legend be shown (default=False)
        @param Legend the list of the legend titles (default=No Legend)
        @param xLabel the name of the xAxis (default=the name of
               the independent variable; string)
        @param yLabel the name of the xAxis (default=the name of
               the dependent variable; string)
        @param alphaLegend the transparency of the legend;
               (0.0=solid; 1.0=transparent; default=0.5)
        
        @todo fix alphaLegend; test options more...
        """
        (results,nodeList,markers,Title,xLabel,yLabel) = getPlotData(
                 self,nodeList,resultType,coord,markers,Title,hasLegend,Legend,xLabel,yLabel)

        i = 0
        Labels = []
        #print "nodeList = ",nodeList
        node0 = nodeList[0]
        Xs = results[node0].keys()

        (fig,ax) = plt.subplots(1)
        #leg = plt.legend(loc='best', fancybox=True,alpha=0.5)
        for nodeID in nodeList: # translation/rotation
            result = results[nodeID]
            Ys = []
            for dt,res in sorted(result.iteritems()):
                if coord==3:
                    val = sqrt(res.dot(res))
                else:
                    val = res[coord]
                ###
                Ys.append(val)
            Label = 'Node %s' %(nodeID)
            Labels.append(Label)
            #ax.plot(Xs,Ys,markers[i],label=Label)
            #plot(Xs,Ys,Label)
            plot(Xs,Ys)
            i+=1
        plt.legend(Labels,loc='best',fancybox=True)
        #print dir(leg)
        #leg.get_frame().set_alpha(0.5)
        ax.set_title(Title)
        #ax.set_legend(Labels)
        #if hasLegend and Legend is None:        
        #    plt.legend(Labels)
        #elif hasLegend:
        #    plt.legend(Legend)
        #plt.axes
        #print dir(plt)
        ax.grid(True)
        ax.set_ylabel(yLabel)
        ax.set_xlabel(xLabel)
        #alpha(alphaLegend)
        show()
        
        
class complexTableObject(scalarObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.gridTypes    = {}
        self.translations = {}
        self.rotations    = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###
        
    def isImaginary(self):
        return True

    def addF06Data(self,data,transient):
        if transient is None:
            for line in data:
                (nodeID,gridType,v1,v2,v3,v4,v5,v6) = line
                self.gridTypes[nodeID] = gridType
                self.translations[self.dt][nodeID] = [v1,v2,v3] # dx,dy,dz
                self.rotations[self.dt][nodeID]    = [v4,v5,v6] # rx,ry,rz
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.translations:
            self.updateDt(self.dataCode,dt)

        for line in data:
            (nodeID,gridType,v1,v2,v3,v4,v5,v6) = line
            self.gridTypes[nodeID]    = gridType
            self.translations[self.dt][nodeID] = [v1,v2,v3] # dx,dy,dz
            self.rotations[self.dt][nodeID]    = [v4,v5,v6] # rx,ry,rz
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],dt,self.iSubcase))
            self.addNewTransient()
        ###

    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def addNewTransient(self,dt):
        """initializes the transient variables"""
        self.translations[dt] = {}
        self.rotations[dt]    = {}

    def add(self,dt,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        #msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s" %(dt,nodeID,v1,v2,v3)
        #assert isinstance(nodeID,int),nodeID
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[nodeID] = [v1,v2,v3] # dx,dy,dz
        self.rotations[nodeID]    = [v4,v5,v6] # rx,ry,rz
    ###

    def addSort1(self,dt,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s" %(dt,nodeID,v1,v2,v3)
        print msg
        if dt not in self.translations:
            self.addNewTransient(dt)
        #print msg
        #msg = ''
        #assert isinstance(nodeID,int),nodeID
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[dt][nodeID] = [v1,v2,v3] # dx,dy,dz
        self.rotations[dt][nodeID]    = [v4,v5,v6] # rx,ry,rz

    def addSort2(self,eid,data):
        [dt,gridType,v1,v2,v3,v4,v5,v6] = data

        if dt not in self.translations:
            self.addNewTransient(dt)

        assert isinstance(eid,int),eid
        assert 0<nodeID<1000000000, msg

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[dt][nodeID] = [v1,v2,v3] # dx,dy,dz
        self.rotations[dt][nodeID]    = [v4,v5,v6] # rx,ry,rz
        
    def plot(self,nodeList=None,resultType='Translation',displayType='Real Imag',coord=3,markers=None,
                  Title=None,hasLegend=True,Legend=None,xLabel=None,yLabel=None,alphaLegend=0.5):
        """
        @param nodeList a list of the node IDs to plot vs the
               independent variable (default=None; all nodes)
        @param resultType the variable to plot ('Translation','Rotation')
        @param displayType 'Real Imag' or 'Mag Phase'
        @param coord the coordinate to plot (for <x,y,z>, x=0,y=1,z=2,Mag=3);
               default=Magnitude
        @param markers  a list of colors/marker shapes for each line
        @param Title title of the plot (default=the object name)
        @param hasLegend should a legend be shown (default=False)
        @param Legend the list of the legend titles (default=No Legend)
        @param xLabel the name of the xAxis (default=the name of
               the independent variable; string)
        @param yLabel the name of the xAxis (default=the name of
               the dependent variable; string)
        @param alphaLegend the transparency of the legend;
               (0.0=solid; 1.0=transparent; default=0.5)
        
        @todo fix alphaLegend; test options more...
        """
        displayType = displayType.title()
        (results,nodeList,markers,Title,xLabel,yLabel) = getPlotData(
                 self,nodeList,resultType,coord,markers,Title,hasLegend,Legend,xLabel,yLabel)

        i = 0
        Labels = []
        node0 = nodeList[0]
        Xs = results[node0].keys()

        (fig,ax) = plt.subplots(1)
        leg = plt.legend(loc='best', fancybox=True,alpha=0.5)
        for nodeID in nodeList: # translation/rotation
            result = results[nodeID]
            Ys = []
            for dt,res in sorted(result.iteritems()):
                if coord==3:
                    val = sqrt(res.dot(res))
                else:
                    val = res[coord]
                ###
                mag   = abs(val)
                phase = degrees(val)
                Ys.append(val)
            Label = 'Node %s' %(nodeID)
            Labels.append(Label)
            ax.plot(Xs,Ys,markers[i])
            #plot(Xs,Ys,Label)
            i+=1
        #print dir(leg)
        #leg.get_frame().set_alpha(0.5)
        ax.set_title(Title)
        #ax.set_legend(Labels)
        if hasLegend and Legend is None:        
            plt.legend(Labels)
        elif hasLegend:
            plt.legend(Legend)
        #plt.axes
        #print dir(plt)
        ax.grid(True)
        ax.set_ylabel(yLabel)
        ax.set_xlabel(xLabel)
        #alpha(alphaLegend)
        show()

def getPlotData(obj,nodeList,resultType,coord,markers,Title,hasLegend,Legend,xLabel,yLabel):
    if Title is None:
        label = ''
        subtitle = ''
        if obj.label:
            label = ' - %s' %(obj.label)
        if obj.subtitle:
            subtitle = ' - %s' %(obj.subtitle)
        Title = '%s%s%s - Subcase %s' %(obj.__class__.__name__,label,subtitle,obj.iSubcase)

    resultType  = resultType.title()
    assert coord in [0,1,2,3],'invalid coord...options=[0,1,2,3].  coord=%s' %(coord)
    assert resultType in ['Translation','Rotation'],"invalid resultType...options=['Translation','Rotation']."

    if xLabel is None:
        xLabel = obj.dataCode['name'].title()
    if yLabel is None:
        if   coord==0:  yLabel='X %s' %(resultType)
        elif coord==1:  yLabel='Y %s' %(resultType)
        elif coord==2:  yLabel='Z %s' %(resultType)
        elif coord==3:  yLabel='%s (Magnitude)' %(resultType)
        else: raise RuntimeError('invalid coord...options=[0,1,2,3].  choice=%s' %(coord))

    (translations,rotations) = obj.getAsSort2()

    if resultType=='Translation':
        results = translations
    else:
        results = rotations

    nodeListAll = results.keys()
    if nodeList is None:
        nodeList = nodeListAll

    if hasLegend and Legend is not None:
        assert len(nodeList)==len(Legend),'len(nodeList)=%s len(legend)=%s' %(len(nodeList),len(legend))
    if markers==None:
        markers = ['-']*len(nodeList)
    #print "subcase = ",obj.iSubcase
    return (results,nodeList,markers,Title,xLabel,yLabel)
