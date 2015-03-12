## Overview ##
The pyNastran BDF interface is intended for use in an optimization process.
Users can take advantage of this to do more complicated things.

This functionality requires version 0.4.**## Mass Properties ##
Mass is usually the minimization goal in most optimization problems.
pyNastran handles most elements/properties.
```
bdf = BDF()
bdf.readBDF(bdfName)
(mass,cg,I) = bdf.MassProperties()
(Ixx,Iyy,Izz,Ixy,Ixz,Iyz) = I
(xCG,yCG,zCG) = cg

#or
mass = bdf.Mass()
```**

## API Parametrization Method ##
pyNastran can be used in an optimization process.
```
bdf = BDF()
bdf.readBDF(bdfName)

# double the thickness on pshell pid=10
p10 = bdf.Property(10)
p10.t *= 2

# change some material properties (Young's Modulus, Poisson's ratio)
m10 = bdf.Material(100)
m10.e = 50E6
m10.nu = 0.29

# triple the thickness of the 0th layer (n-1) of the PCOMP
p11 = bdf.Property(11)

n=0
p.layers[n][1] *=3  # [mid,t,sout,ge]  1->t

# write out the model
bdf.writeBDF(bdfNameOut)
```
See the [Developer Docs](http://code.google.com/p/pynastran/wiki/DeveloperDocs) for more information on how to interface with a BDF.

## User Parametrization Method ##
The user may also parametrize the file themselves and provide a code
with a file that will be parsed by the code that interfaces to pyNastran.
The developer only needs to feed pyNastran a dictionary of the key-value pairs to update the BDF.

Cards that are not supported in pyNastran will not be deparametrized.
Luckily many cards are supported, so this shouldn't be much of an issue.
```
---original line---
FORCE    1       9       0      1000.   1.       0.      0.

--parametrized line using small field format---
FORCE    1       9       0      %oneK   1.       0.      0.

--parametrized line using CSV format---
FORCE,1,9,0,%oneK,1.,0.,0.

```
Large field format and tabs are supported as well.  Variable names have a **maximum of 7-characters** excluding the % sign.
Now we call the code to update the model.  Note that variable names are NOT case sensitive.
```

dictOfVars = {'oneK': 2000.}
bdf = BDF()
bdf.setDynamicSyntax(dictOfVars)
bdf.readBDF(bdfName)
bdf.writeBDF(bdfName2)

---new---
FORCE    1       9       0      2000.   1.       0.      0.
```


## Choosing Which Method To Use ##

The **API Method** is better if you're doing a complicated BDF manipulation or you are hooking into most **existing** systems.

1.  The API Method is probably better for a GUI.  Keeping a valid BDF is always a good thing.

2. It's bad for user customization.

The **User Parametrization Method** is preferred for a simple scripting process.

1. It's more likely to be compatible with a large number of cards that the developer never intended to support.

2. It's follows the syntax of [OpenMDAO's](http://openmdao.org/) [NastranWrapper](https://github.com/OpenMDAO-Plugins/nastranwrapper).

3. It's bad for changing more complicated parameters, such as the number of layers in a PCOMP or adding materials in order to do a Progressive Failure analysis.

4. It's easy to go outside of the 8-character field.  The program will complain if you do.

5.  It only works for the Bulk Data Deck (not the Executive Control or Case Control decks)

6. Not affected by API changes.