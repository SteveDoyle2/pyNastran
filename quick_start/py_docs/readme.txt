The python files in py_docs are valid python code.  However, because it's annoying to write documentation and worry about formatting, a new style is defined.  I'm sure this has been done, but switching shouldn't be hard if it's found.

Rules:

1.  All imports before the first function are not written
2.  Functions are used for flagging sections
3.  The section title is in the first line of the docstring along with the level of the title
4.  Unindented comments are written with ('# ')
5.  Code is just written as normal
6.  #### indicates to write exactly as is written.  It allows for fancy things like URLs, tables, and blank lines to be embedded.
7.  The 'def main():' signals the end of the parsing
8.  All code must have a 'def main():'
9.  Classes are not supported
10. The order in the file is the order that the rst will be written.  Put it in order.
11. The code in main() is executed when py_to_rst.py is run.


A sample function:

def func():
    """
    Section 1 - Level 2
    """
    # this marks unindented comments
    
    this must be valid python code and will be written as:
    >>> this must be valid python code and will be written as:
   
    ##
    ## this is the output of a code line
    
    #### this is a give up, just write the code exactly as specified

#=================================

def func_b():
    """
    Func Title - Level 2
    """
    
    # this is a comment
    import os
    fname = os.path.join('test', 'spike.out')
    print(fname)
    ####
    ## test/spike.out

    # closing comment

def final():
    """"
    Final Block - Level 3
    """"
    # Wasn't that easy?

def main():
    func_b()
    final()

if __name__ == '__main__':
    main()

#=================================

func_b becomes:

==========
Func Title
==========

this is a comment

..code-block:: python

    >>> import os
    >>> fname = os.path.join('test', 'spike.out')
    >>> print(fname)

    test/spike.out

closing comment

-----------
Final Block
-----------

Wasn't that easy?