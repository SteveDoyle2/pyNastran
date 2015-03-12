# Introduction #

fieldWriter.py is a convenient module that allows you to write 8-character width Nastran BDF cards without having to worry about formatting.  You'll get the highest precision possible without having to worry about if the value should be written as a float or in scientific notation.  It will also write values in as little space as possible, so 1.01000 will be written as 1.01.  All BDF cards in pyNastran use the printCard method.

# Details #

Lets say you want to write a CQUAD4 card
```
# helps to avoid errors
card = "CQUAD4,%i,%i,%i,%i,%i\n" %(eid,n1,n2,n3,n4)
print card

or
# easier to read
card = "CQUAD4  %8i%8i%8i%8i%8i\n" %(eid,n1,n2,n3,n4)
print card
```
That was easy because a CQUAD4 has integer fields.  The second method is preferred because the first results in messy BDFs.


Lets say you want to write a GRID card
```
card = "GRID    %8i%8i%8f%8f%8f\n" %(nid,cp,x,y,z)
print card
```
This will work in 90% of cases, but has problems when values are large or small.  You can overbound the 8-character field, which will lead to crashes when running Nastran.

Try it this way:
```
from pyNastran.bdf.fieldWriter import printCard
fields = ['GRID',nid,cp,x,y,z]
card = printCard(fields)
```

What about blank fields:
```
fields = ['GRID',nid,'',x,y,z]
card = printCard(fields)

or

fields = ['GRID',nid,None,x,y,z]
card = printCard(fields)
```

Writing the fields gets even harder if you try to write a multi-line card.  You have to worry about the leading spaces on the 2nd line.  In pyNastran, the printCard method does it for you:
```
# error prone way
card  = 'CHEXA   %8i%8i%8i%8i%8i%8i%8i%8i\n' %(eid,pid,n1,n2,n3,n4,n5,n6)
card += '        %8i%8i\n'                   %(n7,n8)

# note that return lines are added automatically 
fields = ['CHEXA',eid,pid,n1,n2,n3,n4,n5,n6,n7,n8]
card = printCard(fields)
```

## Future Enhancements ##
If you want to write a double precision card (16-character field width):
```
card = printCard(fields,size=16)
```
The default is the same.
```
card = printCard(fields,size=8)
card = printCard(fields)
```