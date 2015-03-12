# Introduction #

The developer docs are distributed with pyNastran.  They're the way to find out what method you'll need to call assuming you know what you want to do.

The documentation is a bit confusing the first time you use them, but once you get the hang of them, they're pretty easy.  The documentation is an extensive auto-generated html document that lists ALL the methods and ALL the variables that you can access, as well as all the source code.  If you're confused, open the QRG (for the BDF) or the DMAP guide (for the OP2) to help you.

Most of the BDF/OP2/F06 methods that you'll want to use are documented.  Methods that are not documented, you probably want to stay away from as they're probably internal methods that just haven't been documented yet.

For the BDF cards, those are largely undocumented.  There are tons of variables and methods, so it's a bit of work.  Having the QRG handy will help.  As always, if you find a bug, make a ticket.

## Building the docs ##
If you downloaded the software, skip this step.

If you checked it out, you'll need Doxygen.  If you're not on Windows, with Python 2.7, open **Doxyfile**, search for **Python27** and fix the path to python.  Then call **doxygen Doxyfile**.  This should make a docs folder assuming you had no errors.

# Using the docs/html to make use of the BDF class #
There are tons of methods or object attributes that you can access.  The docs are laid out for a developer (someone using pyNastran in their code).

Open the **docs/html** folder:
click on **any html file**.  Your browser will open.

## Investigating BDF methods ##
Let's look at the main bdf files.

Click the **Directories** tab.  Click **pyNastran/bdf**.  You should see things like **bdf.py**, **fieldWriter.py**, **caseControl.py**, etc.  Files in this folder have functions that will probably be useful to you.  Files in **pyNastran/bdf/bdfInterface** deal with the reader and creating objects.  Generally you don't need to look at those.

Click on **bdf.py**.  It will give you a list of classes.  Click **pyNastran.bdf.bdf.BDF** and it will bring you to the **Classes** tab.

You'll see the list of public member functions with a brief docstring.  Click on **readBDF** or just scroll down to the method for a more detailed description.  You'll see:
```
def pyNastran.bdf.bdf.BDF.readBDF( self,
		  	           infilename,
		  	           includeDir = None,
		  	           xref = True)

main read method for the bdf

Parameters:
    infilename	the input bdf
    includeDir	the relative path to any include files (default=None if no include files)
    xref	should the bdf be cross referenced (default=True)

Definition at line 442 of file bdf.py.

```

Scroll back up to the top and find the updateSolution method.

```
def pyNastran.bdf.bdf.BDF.updateSolution(self,
		  	                 sol,
		  	                 method = None) 		

updates the overall solution type (e.g. 101,200,600)

Parameters:
    self	the object pointer
    sol  	the solution type (101,103, etc)
    method	the solution method (only for SOL=600), default=None

Definition at line 540 of file bdf.py.
```

## Investigating BDF Card methods ##
Go back to Directories.  Click **pyNastran/bdf/cards**.  The files in this sub-folder and any levels lower are the interfaces to all the BDF objects.  You should be referencing this directory and it's sub-directories frequently.

Go back to Directories.  Click **pyNastran/bdf/cards/plates** to see things like **elementsShell.py** and **propertiesShell.py**.  Click on "propertiesShell.py" to see the ShellProperty, PSHELL, PCOMP, PCOMPG classes.  Click on **ShellProperty** and look at the class diagram.  The ShellProperty class points to all the other shell classes and is useful when coding.  Click on the **PCOMP** box in the class diagram.

Again you'll see the Public Member Functions.  Find and click on **MassPerArea**.  You'll see:
```
def pyNastran.bdf.cards.plates.propertiesShell.PCOMP.MassPerArea(self,
      		                                                 iPly = 'all') 		

mass = rho*A*t but area comes from the element mass/A = rho*t for the various layers
the final mass calculation will be done later

Parameters:
    self	the object pointer
    iPly	the string 'all' (default) or the mass per area of the ith ply

Definition at line 449 of file propertiesShell.py.
```
Click on "line 449" and you'll see the following.  You can verify the method or learn how the function works if there is a missing docstring.
```
00425     def MassPerArea(self,iPly='all'):
00426         """
00427         mass = rho*A*t
00428         but area comes from the element
00429         mass/A = rho*t for the various layers
00430         the final mass calculation will be done later
00431         @param self the object pointer
00432         @param iPly the string 'all' (default) or the mass per area of the ith ply
00433         """
00434         if iPly=='all': # get all layers
00435             massPerArea = 0.
00436             for (iply,ply) in enumerate(self.plies):
00437                 massPerArea += self.MassPerArea(iply)
00438             ###
00439             if self.isSymmetrical():
00440                 if self.hasCoreLayer():
00441                     massPerArea -= self.MassPerArea(0)/2.  # cut out the thickness of half the core layer
00442                     
00443                 ###
00444                 return massPerArea*2.
00445             ###
00446             return massPerArea
00447         ###
00448         else:
00449             rho = self.Rho(iPly)
00450             t = self.plies[iPly][1]
00451             return rho*t+self.nsm
00452         ###
```
In this version of the docs, you'll see the docstring is wrong and that non-structural mass is considered in the massPerArea method.  The docstring also made no mention of symmetry or core layers, but the function handles that.