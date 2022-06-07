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
<http://pypi.python.org/pypi/spyder/>`_ or `WingIDE <https://wingware.com/>`_
with `pylint <http://pypi.python.org/pypi/pylint/>`_ integration to help with
identifying `PEP-8 <http://www.python.org/dev/peps/pep-0008/>`_ deficiencies.


Most important points of Python PEP-8 (http://www.python.org/dev/peps/pep-0008/):

#. Use spaces (not tabs) when indenting code
#. Try to limit code to **80**-characters width
#. Files are named **as_follows.py**
#. Classes are named **AsFollows**
#. Functions are are named **as_follows()** and private functions are named
   **_as_follows()**
#. Variable names**_have_embedded_underscores** or without any (e.g. self.elements).
   A private variable would look **_similar_to_a_function**
#. Don't use **import \*** except for tests
#. **Break the rules if it's clearer**
#. Document functions/classes/variables


Additional Guidelines
----------------------

#. :mod:`pyNastran` liberally uses classes (it has too many!), so try not to
   write more classes than necessary. Functions are often better to avoid having
   extra functions in objects the user has immediate access to (imagine having
   100 functions to sift through) in an IDE.  Watch `this video
   <http://www.youtube.com/watch?v=o9pEzgHorH0>`_ for a better explanation.
   Again, it's a work in progress.

#. When writing classes, group similar classes in the same file to reduce having
   to search for a class across multiple files.

#. When deprecating functions, throw a warning using the Python warnings module.
   It's easier to make sure the rest of the code isn't using it this way.  This
   is not nearly as important if it's not a function the user is likely to interact with.

Documentation
--------------

:mod:`pyNastran` uses `Sphinx <http://www.sphinx-doc.org/en/master/>`_ for docstrings with
Numpy's `style guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
for documenting variables.  This style is nice because:

 1. Sphinx is standard for documenting Python code.  As such, it works with
    `ReadTheDocs <http://pynastran-git.readthedocs.org/en/main/>`_, which
    auto-generates documentation.
 2. Numpy's style makes it easy to read in the source code
 3. Data types are clearly defined
 4. It's supported in IDE's like `WingIDE <https://wingware.com/>`_.
 5. It supports embedded LaTeX

An example is shown below.

.. code-block:: python

  def read_op2(op2_filename=None, combine=True,
               log=None, debug=True, debug_file=None, build_dataframe=False,
               skip_undefined_matrices=True, mode='msc'):
      """
      Creates the OP2 object without calling the OP2 class.

      Parameters
      ----------
      op2_filename : str (default=None -> popup)
          the op2_filename
      combine : bool; default=True
          True : objects are isubcase based
          False : objects are (isubcase, subtitle) based;
                  will be used for superelements regardless of the option
      build_dataframe : bool; default=False
          builds a pandas DataFrame for op2 objects
      skip_undefined_matrices : bool; default=False
           True : prevents matrix reading crashes
      debug : bool; default=False
          enables the debug log and sets the debug in the logger
      log : Log()
          a logging object to write debug messages to
       (.. seealso:: import logging)
      debug_file : str; default=None (No debug)
          sets the filename that will be written to

      Returns
      -------
      model : OP2()
          an OP2 object

      .. todo :: creates the OP2 object without all the read methods

      .. note :: this method will change in order to return an object that
                 does not have so many methods
      """
      # code goes here
