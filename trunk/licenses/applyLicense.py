import os
import sys
import copy
import pyNastran

def getFoldersFiles(dirname):
    """
    I'm sure there's an automatic way to create this...
    try using os.walk
    """
    print(dirname)
    filesFolders = os.listdir(dirname)
    
    files = []
    folders = []
    for fileFolder in filesFolders:
        if '.svn' not in fileFolder:
            fullPath = os.path.join(dirname,fileFolder)
            if os.path.isdir(fullPath):
                (subFolders,fullPaths) = getFoldersFiles(fullPath)
                folders.append(fullPath)
                for folder in subFolders:
                    folders.append(folder)
                ###
                for fullPath in fullPaths:
                    #print "A %s" %(fullPath)
                    files.append(fullPath)
                ###
            else:
                if '.pyc' not in fileFolder and '.pyx' not in fileFolder and '.bdf' not in fileFolder and '__init__.py' not in fileFolder:
                    #print "B %s" %(fullPath)
                    files.append(fullPath)
                ###
                else:
                    pass
                    #print "skipping %s" %(fullPath)
                ###
            ###
        ###
    ###
    return (folders,files)

def printTree(files):
    """
    I'm sure there's an automatic way to print this...
    """
    levelsPrinted = {}
    levelOrder = []
    for fname in files:
        fullDir = [os.path.basename(fname)]
        while os.path.dirname(fname):
            fullDir.append(os.path.basename(os.path.dirname(fname)))
            fname = os.path.dirname(fname)
            print "fname = ",fname
            
        fullDir.reverse()
        
        nLevels = len(fullDir)
        for i,level in enumerate(fullDir):
            if level not in levelsPrinted:
                levelsPrinted[level] = i
                levelOrder.append(level)
            
        #print levelOrder
        #print fullDir
        #sys.exit()
    ###
    for level in levelOrder:
        nLevels = levelsPrinted[level]
        if '.py' in level:
            print "    "*nLevels+level
        else:
            print "    "*nLevels+level+'/'
        ###
    ###
    #sys.exit()

def cleanHeader(lines):
    """
    removes the copyright lines from the file
    used when the copyright statement is updated.
    """
    markerWord = 'changes to it will be lost.' # this is the last line with header text
    lines2 = copy.deepcopy(lines)
    for i,line in enumerate(lines[0:50]):
        line2 = lines[i]
        print "line2[%s] = |%r|" %(i,line2)
        if markerWord in line2:
            lines2.pop(0) # removes last ## line
            lines2.pop(0) # removes last ## line
            print "breaking on line %s" %(i)
            break
        elif line2 and line2[0:3]=='## ':
            lines2.pop(0)
        else:
            break
        ###
    ###
    #print "i = ",i
    return lines2

def updateCopyright():
    # must be in ../licenses
    dirname = os.path.relpath(os.path.dirname(pyNastran.__file__))
    print("dirname = %s" %(dirname))
    (folders,files) = getFoldersFiles(dirname)
    printTree(files)

    header = open(os.path.join('licenses','header.txt'),'r')
    headerLines = header.readlines()
    header.close()

    for fname in files:
        if '.py' in fname:
            print "adding fname=%s" %(fname)
            f = open(fname,'r')
            lines = f.readlines()
            f.close()

            for line in lines[0:20]:
                print line

            f.close()
            lines2 = cleanHeader(lines)

            #break
            lines3 = headerLines + lines2

            fout = open(fname,'wb')
            for line in lines3:
                fout.write(line)
            ###
            #if 'bdf.py' in fname:
            #    break
        ###
    ###    
    #print "---FOLDERS---"
    #print "folders = ",'\n'.join(folders)+'\n'

    print "---FILES---"
    print "files = ",'\n'.join(files)
###

if __name__=='__main__':
    updateCopyright()
