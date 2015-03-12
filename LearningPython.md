# Introduction #

Python is a one of the easiest languages there is to program in.  It's easier than C++, Matlab, Perl, Fortran (all languages I've used).  It's also free.  It comes with a lot of built-in functions, which means you won't ever have to code those very common functions (path manipulation, file copying, etc.).  The book [Learn Python the Hard Way](http://learnpythonthehardway.org/book/) (free for the HTML version) is a good place to start.  [Dive Into Python](http://www.diveintopython.net/) is another great HTML book.

## Basic Math (Python 2.x) ##
Python distinguishes between integers and floats.
When you first learned division:
5/3 = 1 remainder 2

Eventually, you started doing:
5/3 = 1.666

Python treats:
5/3 = 5//3 = 1 (the truncated result)

and:
5%3 = 2 (the remainder)

5/3.0 = 5.0/3 = 1.666 (an input is a float, so the outputs is a float)

It's strongly recommended that all integer division in Python be done as:
5//3 and not 5/3

because:
5/3 = 1.666
in Python 3.

## Recommendations ##

[Anaconda](https://store.continuum.io/cshop/anaconda/) - Python distribution

[Wing IDE](http://www.wingware.com/) - development environment with an amazing debugger

[![](http://wingware.com/images/wingware-logo-180x58.png)](http://www.wingware.com/)

[Numpy](http://www.numpy.org/), [Scipy](http://scipy.org/) - numerical libraries

[Matplotlib](http://matplotlib.org/) - plotting library