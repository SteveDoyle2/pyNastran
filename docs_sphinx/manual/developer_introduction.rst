==================================
Developer: Getting Started
==================================


Code Formatting - PEP-8
--------------------------


:mod:`pyNastran` is attempting to implement Python's `PEP-8
<http://www.python.org/dev/peps/pep-0008/>`_ formatting standard after not
following it for multiple versions.  It's a work in progress to upgrade the
code, so try to follow it. The biggest challenge has been to come up with a
clear way to reference methods and very similar variable names.
For example, the element.Area() method uses the element.area variable in the
case of the CQUAD4, but references element.areas in the case of a
multi-segmented CBAR.  It's recommended that developers use `Spyder
<http://pypi.python.org/pypi/spyder/>`_ with `pylint
<http://pypi.python.org/pypi/pylint/>`_ integration to help with identifying
`PEP-8 <http://www.python.org/dev/peps/pep-0008/>`_ deficiencies.


Most important points of Python PEP-8 (http://www.python.org/dev/peps/pep-0008/):
     
#. Use spaces (not tabs) when indenting code
#. Try to limit code to **80**-characters width
#. Files are named **as_follows.py**
#. Classes are named **AsFollows**
#. Functions are are named **as_follows()** and private functions are named 
   **_as_follows()**
#. Variable names**_have_embedded_underscores** or without any (e.g. self.elements).  
   A private variable would look **_similar_to_a_function**
#. Don't use **import \***
#. **Break the rules if it's clearer**
#. Document functions/classes/variables


Additional Guidelines
----------------------

#. :mod:`pyNastran` liberally uses classes (it has too many!), so try not to
   write more classes than necessary. Functions are often better to avoid having
   extra functions in objects the    user has immediate access to (imagine having
   100 functions to sift through)    in an IDE.  Watch `this video 
   <http://www.youtube.com/watch?v=o9pEzgHorH0>`_    for a better explanation.
   Again, it's a work in progress.

#. When writing classes, group similar classes in the same file to reduce having
   to search for a class across multiple files (very common in Java).

#. When deprecating functions, throw a warning using the Python warnings module.
   It's easier to make sure the rest of
   the code isn't using it this way.  This is not nearly as important if it's
   not a function the user is likely to interact with.

Documentation
--------------


:mod:`pyNastran` uses `Doxygen <http://www.doxygen.org/>`_ for docstrings. 
Doxygen is nice because it allows LaTeX to be embedded in the documentation and
HTML docs are easily generated.  Sphinx is another documentation format, It is
known that Sphinx also allows LaTeX and that it is much more commonly used.

Some guidelines for Doxygen:

:@param:  in\_variable\_name a description goes here
:@retval: out\_variable\_name another description
:@see:    look at AsFollows class or as\_follows function
:@note:   a message
:@warning: a warning
:@todo:   something that's not done
