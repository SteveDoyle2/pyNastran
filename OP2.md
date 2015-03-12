# Introduction #

This is how the version 0.4 OP2 reading software works.  Make sure you install the software first!  See InstallationGuide for instructions.

# Features #
  * multiple subcases
  * ability to set list of subcases to extract using setSubcases
  * extract specific times from a transient solution using setTransientTimes
  * lots of transient/static results
    * displacements (real/complex)
    * temperature
    * eigenvectors (real/complex)
    * SPC/MPC Forces (real/complex)
    * linear stress/strain
      * spring (CELASi)
      * rod (CROD, CONROD)
      * beam (CBEAM)
      * plate (CQUAD4 bilinear, CQUAD4 centroidal, CTRIA3)
        * maxShear/vonMises, fiberDistance/plateCurvature support
      * composite plate (CQUAD4 centroidal, CTRIA3)
        * maxShear/vonMises support
      * solid elements (CTETRA, CPENTA, CHEXA)
    * nonlinear stress/strain
      * CQUAD4, CTRIA3, CROD, CONROD, CTUBE, CTETRA
    * element strain energy support (limited)
  * automated printing of tables

## Limitations ##
  * no SORT2 support (only SORT1)
  * **skippedCards.out** (run **test\_op2 op2Name.op2**) shows the results that were skipped
  * limited thermal support (only temperature)
  * no heat flux
  * no applied forces
  * limited nonlinear/hyperelastic element support (e.g. SOL 106, SOL 129)
  * VERY limited frequency support
  * see 

&lt;A href="https://pynastran.googlecode.com/svn/tags/0.4.0/op2\_reading.txt"&gt;

 v0.4 OP2 Support

&lt;/A&gt;



## Results Table ##
The following tables are supported, unless otherwise stated.
| **Result Type**      | **BDF Line**        | **F06/OP2 Table** |
|:---------------------|:--------------------|:------------------|
| Displacement       | DISPLACEMENT=ALL  | OUG  |
| Velocity           | VELOCITY=ALL      | OUG  |
| Acceleration       | ACCELERATION=ALL  | OUG  |
| Eigenvector        | DISPLACEMENT/VECTOR=ALL  | OUG  |
| Temperature        | THERMAL=ALL       | OUG  |
| Eigenvalue Order   | ESE=ALL           | LAMA |
| SPC Forces         | SPCFORCES=ALL     | OQG  |
| MPC Forces         | MPCFORCES=ALL     | OQG  |
| Grid Point Forces  | GPFORCE=ALL       | OGF  |
| Applied Forces     | ???               | OGP  |
| Load Vector        | ???               | OGP  |
| Stress             | STRESS=ALL        | OES  |
| Strain             | STRAIN=ALL        | OES  |
| Element Strain Energy | ESE=ALL ???    | OEE  |
| Forces             | FORCE/ELFORCE=ALL | OEF (not supported)  |
| Heat Flux          | FLUX=ALL          | OEF (not supported)  |


# Making an OP2 #
If you cant make an op2, check your bdf:
```
CEND
$ case control deck...
BEGIN BULK
$PARAM,POST,0      $ makes an xdb for Patran
PARAM,POST,-1      $ the line that makes a good op2
PARAM,POST,-2      $ makes an op2 that WILL fail
$NODES
GRID,1,0,1.,2.,3.
$ rest of bulk data deck
```


# Getting Started #
First install pyNastran v0.4 and on the command line type:
```
>> test_op2 fem.op2 > test.out
```

The "hardest" test would be the following:
```
>> test_op2 -fgw fem.op2 > test.out
```
This will read all sections of the op2, print out a copy of the geometry to fem.bdf.out, and print an F06 of the majority of the OP2 results to fem.f06.out

This following will run the model and print out all the results to the screen.  Additional options can be seen by typing:
```
>> test_op2 --help

usage: test_op2-script.py [-h] [-q] [-g] [-w] [-f] [-v] op2FileName

Tests to see if an OP2 will work with pyNastran.

positional arguments:
  op2FileName     path to OP2 file

optional arguments:
  -h, --help      show this help message and exit
  -q, --quiet     Prints debug messages (default=True)
  -g, --geometry  Reads the OP2 for geometry, which can be written out
  -w, --writeBDF  Writes the bdf to fem.bdf.out
  -f, --writeF06  Writes the f06 to fem.f06.out
  -v, --version   show program version number and exits
```

# Using the OP2 Reader #

```

from pyNastran.op2.op2 import OP2

op2Obj = OP2('fem.op2',makeGeom=False,debug=False,log=None)
# the OP2 filename is loaded
# makeGeom=True is required in order to write out the BDF
# debug=True prints error messages
# log is None (dummy logger) or a python logging module object
# see http://docs.python.org/library/logging.html

# limits reading to subcase 100 & 204; be careful with SOL 200, it may ONLY has 1 subcase, despite analyzing 2 cases
op2.setSubcases([100,204])

# limits reading of transient results (good for large models)
# the nearest times will be found, if no time is close, the time will be skipped
times = {100:[0.1,0.2,0.3,1.0], 204:[] }  # 1:[]=all (just like setSubcases)
op2.setTransientTimes(times)  # added in v0.4

# read the OP2
op2Obj.readOP2()

# writes the bdf if makeGeom is True, otherwise the file is empty
op2Obj.writeBDFAsPatran('fem.bdf.out')

# print the results as an F06 (not supported across all results)
print op2Obj.writeF06(f06Name) # added in v0.4

# lets print the results
print op2Obj.printResults() # prints ALL results from ALL tables

# It prints the following formatted results
---ALL RESULTS---
iSubcase = 1
---DISPLACEMENTS---
nodeID   gridType         Dx         Dy         Dz         Rx         Ry         Rz 
1        G         1.523e-02 -1.641e-06 -2.276e-03          0          0          0 
2        G         1.530e-02          0 -2.300e-03          0          0          0 


iSubcase = 200
---DISPLACEMENTS---
nodeID   gridType         Dx         Dy         Dz         Rx         Ry         Rz 
1        G         1.523e-02 -1.641e-06 -2.276e-03          0          0          0 
2        G         1.530e-02          0 -2.300e-03          0          0          0 

iSubcase = 1
---CELASx STRESSES---
EID     eType     stress 
10001  CELAS1          0 
10002  CELAS2          0 

iSubcase = 200
---CELASx STRESSES---
EID     eType     stress 
10001  CELAS1          0 
10002  CELAS2          0 

iSubcase = 1
---CELASx STRAINS---
EID     eType     strains
10001  CELAS1          0 
10002  CELAS2          0 

iSubcase = 200
---CELASx STRAINS---
EID     eType     strains
10001  CELAS1          0 
10002  CELAS2          0 



```
# Static Table #
All static tables follow this format, so if you can access results from this example, you can access results from any static table.

```
iSubcase = 1
---ISOTROPIC PLATE STRESS---
EID     eType   nodeID  iLayer  curvature        oxx        oyy        txy     majorP     minorP   maxShear 
110001 CQUAD4        C       1      -0.05     -17683      -7216      -9280      -1795     -23103      22260 
110001 CQUAD4        C       2       0.05     -17762      -7252      -9275      -1846     -23168      22302 
110002 CQUAD4        C       1      -0.05      -4872        750       -904        891      -5014       5514 
110002 CQUAD4        C       2       0.05      -4852        740       -907        883      -4996       5492 
110003 CQUAD4        C       1      -0.05      -1870       -362       -284       -310      -1922       1787 

# "C" for nodeID means Centroid.
# if nodeID was a number, that would refer to the actual node ID
# iLayer=1,2 for an isotropic material means Z1/Z2 (Top Half/Bottom Half).


# lets print all the plate stress only (static loading only)
# and then print oxx and ovm (x-stress and von mises stress)
# you can use this instead to process the data (e.g. calculate margins)
# generally all tables follow this pattern
for iSubcase,stressCase in op2Obj.plateStress.items():
    print stressCase
    for eid,oxxNodes in sorted(stressCase.oxx.items()):
        eType = stressCase.eType[eid]
        for nid in sorted(oxxNodes):
            for iLayer in range(len(stressCase.oxx[eid][nid])):
                fd    = stressCase.fiberCurvature[eid][nid][iLayer] # fiber/curvature
                oxx   = stressCase.oxx[eid][nid][iLayer]
                oyy   = stressCase.oyy[eid][nid][iLayer]
                txy   = stressCase.txy[eid][nid][iLayer]
                angle = stressCase.angle[eid][nid][iLayer]
                major = stressCase.majorP[eid][nid][iLayer]
                minor = stressCase.minorP[eid][nid][iLayer]
                ovm   = stressCase.ovmShear[eid][nid][iLayer] # VonMises/MaxShear


# To get ONLY plate stress for subcase 1:
stressSubcase1 = op2Obj.plateStress[1]

# to check if the results contains Von Mises Stress or max Shear Stress
# look at the printout or in code...
# NOTE:  this will be consistent across the entire subcase
isVM = stressSubcase1.isVonMises() # returns True/False
isMaxShear = stressSubcase1.isMaxShear() # returns True/False

# to check if the results contains Fiber / Curvature (composites only)
# look at the printout or in code...
isFiber = stressSubcase1.isFiberDistance() # returns True/False
isCurvature = stressSubcase1.isCurvature() # returns True/False

# we know by looking at the printout or by knowing what was run, but we
# can also query the object to find out if is transient stress or not; 
# this works for all result objects
isTransient = stressSubcase1.isTransient # returns True/False

```
# Transient Table #
```

# If stress looked like this instead...
---ISOTROPIC PLATE STRESS---
EID     eType   nodeID  iLayer  curvature        oxx        oyy        txy     majorP     minorP   maxShear 
dt = 0.0
11     CQUAD4        C       1     -0.025          0          0          0          0          0          0 
11     CQUAD4        C       2      0.025          0          0          0          0          0          0 
12     CQUAD4        C       1     -0.025          0          0          0          0          0          0 
12     CQUAD4        C       2      0.025          0          0          0          0          0          0 
dt = 0.0249999985099
11     CQUAD4        C       1     -0.025       7140      22854       -117      22855       7139      20252 
11     CQUAD4        C       2      0.025      -7140     -22854        117      -7139     -22855      20252 
12     CQUAD4        C       1     -0.025       7140      22854        117      22855       7139      20252 
12     CQUAD4        C       2      0.025      -7140     -22854       -117      -7139     -22855      20252 

# I'll use the transient version instead (print still works the same)
# only accessing the data changes
for iSubcase,stressCase in op2Obj.plateStress.items():
    print stressCase
    for dt,oxxs in sorted(stressCase.oxx.items()):
        msg += '%s = %s\n' %(stressCase.dataCode['name'],dt)
        for eid,oxxNodes in sorted(oxxs.items()):
            eType = stressCase.eType[eid]
            for nid in sorted(oxxNodes):
                for iLayer in range(len(self.oxx[dt][eid][nid])):
                    fd    = stressCase.fiberCurvature[dt][eid][nid][iLayer]
                    oxx   = stressCase.oxx[   dt][eid][nid][iLayer]
                    oyy   = stressCase.oyy[   dt][eid][nid][iLayer]
                    txy   = stressCase.txy[   dt][eid][nid][iLayer]
                    angle = stressCase.angle[ dt][eid][nid][iLayer]
                    major = stressCase.majorP[dt][eid][nid][iLayer]
                    minor = stressCase.minorP[dt][eid][nid][iLayer]
                    ovm   = stressCase.ovmShear[dt][eid][nid][iLayer]

# NOTE:  dataCode contains important information about the case that
# was analyzed (e.g. loadSet, method, the name of the "transient"
# variable, etc.)
# valid transient variable names include:  None, freq, dt, loadStep



```
# More Result Tables #
```

iSubcase = 1
---COMPOSITE PLATE STRESS---
EID       eType   iLayer        o11        o22        t12        t1z        t2z  oVonMises 
1001    QUAD4LC        1       1729       1273        667          0          0        705 
1211    TRIA3LC       20      -5397      -1694      -1453          0          0       2353 

## for a QUAD4LC, the stress only exists at the center of the element

iSubcase = 1
---COMPOSITE PLATE STAIN---
EID       eType   iLayer        e11        e22        e12        e1z        e2z  eVonMises 
1001    QUAD4LC        1    0.00011   9.19e-05   6.02e-05          0          0   6.28e-05 

iSubcase = 1
---COMPOSITE PLATE STAIN---
EID       eType   iLayer        e11        e22        e12        e1z        e2z  eMaxShear 
1001    QUAD4LC        1    0.00011   9.19e-05   6.02e-05          0          0   6.28e-05 

# note that in the top first composite strain example, eVonMises is
# used, while in the second eMaxShear is used.  isVonMises() will 
# evaluate to False

iSubcase = 1
---BAR STRAIN---
EID       eType         e1         e2         e3         e4      Axial       eMax       eMin 
22052      CBAR   -0.00062    0.00017   0.000534  -0.000583  -0.000525   8.99e-06   -0.00114 
                  -0.00115   0.000809   0.000555   -0.00118              0.000284    -0.0017 
iSubcase = 1
---BAR STRESS---
EID     eType     s1     s2     s3     s4  Axial   sMax   sMin 
22052    CBAR  -6383   1754   5502  -6008  -5409     92 -11793 
              -11870   8336   5718 -12132          2927 -17542 

iSubcase = 1
---APPLIED LOADS---
nodeID   eID        source         Fx         Fy         Fz         Mx         My         Mz 
----------------------------------------------------------------------------------------------------
1        1        HEXA         28.956      5.041   2576.704          0          0          0 
1        6        QUAD4        22.757     -3.334  -1276.954      0.007     -0.074      0.624 
1        10       TRIA3       -37.645     -1.545    -14.417      0.685      0.038     -1.486 
1        11       TRIA3       -14.068     -0.163  -1285.333     -0.692      0.036      0.862 
1        =        *TOTALS*          0          0          0          0          0          0 
----------------------------------------------------------------------------------------------------

iSubcase = 1
---EIGENVECTORS---
-eigenvalues-
0  1.000000
1  3.000000

eigenvalueReal = 1.000000
nodeID   GridType         Tx         Ty         Tz         Rx         Ry         Rz 
100000          G          0          0          0          0          0          0 
100001          L          0          0          0          0          0          0 
101             S  -0.000399          0          0          0          0          0 

eigenvalueReal = 3.000000
nodeID   GridType         Tx         Ty         Tz         Rx         Ry         Rz 
100000          G          0          0          0          0          0          0 
100001          L          0          0          0          0          0          0 
101             S  -0.000399          0          0          0          0          0 

iSubcase = 1
---TRANSIENT TEMPERATURE---
NodeID     GridType Temperature
dt = 0
         1        S        110
         2        S        110
         3        S        110
dt = 110
         1        S    112.945
         2        S    110.147
         3        S    109.838


iSubcase = 1
---BEAM STRESSES---
EID     eType    NID    xxb      sMax       sMin MS_tension MS_compression 
1       CBEAM      1 0.00    -249940    -249940          0          0 
1       CBEAM      3 1.00     499251     499251          0          0 


iSubcase = 2
---TRANSIENT DISPLACEMENTS---
NodeID   GridType         Dx         Dy         Dz         Rx         Ry         Rz 
lftsfq = 1.0
1        G                 0          0          0          0          0          0 
2        G                 0          0          0          0          0          0 
3        G                 0          0          0          0          0          0 

iSubcase = 1
---TRANSIENT DISPLACEMENTS---
NodeID   GridType         Dx         Dy         Dz         Rx         Ry         Rz 
dt = 0
1        G                 0          0          0          0          0          0 
2        G                 0          0          0          0          0          0 
3        G                 0          0          0          0          0          0 

# NOTE the difference in the variable names (dt vs lftsfq)

```


# Table Data #

Static Tables
  * Temperature
  * Displacement
  * Strain Energy Density
  * Applied Loads Summary
  * Stress/Strain (rod-type, bar-type, beam-type, compositeShell-type, shell-type, solid-type)


Transient Tables
  * Temperature
  * Displacement
  * Eigenvectors
  * Strain Energy Density
  * Applied Loads Summary
  * Stress/Strain (rod-type, bar-type, beam-type, compositeShell-type, shell-type, solid-type)

The tables are organized into cascading dictionaries with the final value being a float unless otherwise noted.

### Notes ###
  * In order to get the static table, drop the **[dt](dt.md)** term from the transient table.
  * The actual name of dt can be found by extracting the object.dataCode['name'] string.
  * eType means element Type
  * eid means element ID
  * nid means node ID

Possible values of gridType include:
  * G (GRID)
  * S (SPOINT)
  * L (Rigid Element Point - not 100% sure)

The list of data members in the various results:
```

 * Temperature
    * gridType[nid]
    * T[dt][nid]

 * Displacement / Velocity / Acceleration / SPC Forces / MPC Forces (Real)
   * gridType[nid]
   * translations[dt][nid]
     * numpy array of length 3
   * rotations[dt][nid]
     * numpy array of length 3

 * Displacement / SPC Forces / MPC Forces (Complex)
   * gridType[nid]
   * translations[dt][nid]
     * numpy array of length 6
   * rotations[dt][nid]
     * numpy array of length 6

 * Real Eigenvalue Order
   * extractionOrder[modeNum]
   * eigenvalues[modeNum]
   * radians[modeNum]
   * cycles[modeNum]
   * generalizedMass[modeNum]
   * generalizedStiffness[modeNum]

 * Eigenvectors
   * gridType[nid]
   * displacements[dt][nid]
     * numpy array of length 3
   * rotations[dt][nid]
     * numpy array of length 3

 * Complex Eigenvectors
   * gridType[nid]
   * translations[dt][nid]
     * numpy array of length 6
   * rotations[dt][nid]
     * numpy array of length 6

 * Strain Energy Density
    * energy[dt][eid]
    * percent[dt][eid]
    * density[dt][eid]

 * Applied Loads Summary
    * eids[dt][nid][iLoad]
      * iLoad represents a list
    * source[dt][nid][iLoad]
    * forces[dt][nid][iLoad]
      * numpy array of length 3
    * moments[dt][nid][iLoad]
      * numpy array of length 3

 * Stress/Strain
   * rod-type
     * eType[eid]
     * axial[dt][eid]
     * torsion[dt][eid]
     * MS_axial[dt][eid]
     * MS_torsion[dt][eid]

   * bar-type (replace "s" for "e" for strain)
     * eType[eid]
     * s1[dt][eid][iPoint]
       * NOTE:  iPoint is a length=2 list of points A and B of the bar
     * s2[dt][eid][iPoint]
     * s3[dt][eid][iPoint]
     * s4[dt][eid][iPoint]
     * axial[dt][eid]
     * smax[dt][eid][iPoint]
     * smin[dt][eid][iPoint]

   * beam-type (replace "s" with "e" when dealing with strain)
     * eType[eid]
     * xxb[dt][eid][iPoint]
       * NOTE:  iPoint is a variable length list of the number of points along the beam
     * grids[dt][eid][iPoint]
     * smax[dt][eid][iPoint]
     * smin[dt][eid][iPoint]
     * MS_tension[dt][eid][iPoint]
     * MS_compression[dt][eid][iPoint]

   * compositePlate-type (replace "e" with "o" and "t" for strain)
      * eType[eid]
      * o11[dt][eid][iLayer]
        * NOTE: iLayer is a variable length list depending on the number of plies
      * o22[dt][eid][iLayer]
      * t12[dt][eid][iLayer]
      * t1z[dt][eid][iLayer]
      * t2z[dt][eid][iLayer]
      * angle[dt][eid][iLayer]
      * majorP[dt][eid][iLayer]
      * minorP[dt][eid][iLayer]
      * ovmShear[dt][eid][iLayer]

   * plate-type
     * NOTE:  same as compositePlate-type, but 1 added parameter
     * fiberCurvature[dt][eid][iLayer]
      * iLayer refers to Z_1 or Z_2

   * solid-type (same rule for strain as compositePlate-type)
      * eType[eid]
      * oxx[dt][eid][nid]
      * oyy[dt][eid][nid]
      * ozz[dt][eid][nid]
      * txy[dt][eid][nid]
      * tyz[dt][eid][nid]
      * txz[dt][eid][nid]
      * ovmShear[dt][eid][nid]

```