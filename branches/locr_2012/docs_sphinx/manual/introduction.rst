========================
Brief Project Overview
========================



Since the 1960's `NASTRAN <http://en.wikipedia.org/wiki/Nastran>`_ (NASA
Structural ANalysis) has been used to solve
structural/thermal/aerodynamic/dynamics/etc. problems.  The file formats were
originally developed by MSC for a product now called MSC Nastran.  There have
been many spinoff version of NASTRAN that have been created based on the 2001
source code release of MSC Nastran in 2002 after settlement with the FTC
(Federal Trade Commisson).  There is now NX Nastran and NEi Nastran, which are
developed independently.

:mod:`pyNastran` is at it's core an API (Application Programming Interface) to
the legacy formats used by Nastran.  These files include the BDF, F06, OP2, OP4,
and PCH files.  Other code has been added to :mod:`pyNastran` in order to drive
development of the software in order to be able to solve other engineering
problems.  For example, `Code\_Aster <http://www.code-aster.org>`_, an
open-source finite element code developed by the EDF (Electricity of France),
has a Nastran to `Code\_Aster <http://www.code-aster.org>`_ converter that is in
development.  The development has helped to define the API for the loads in
order to be able to extract them in a way that makes sense.  However, this is
not the focus of the software.



Target Audience
-----------------


:mod:`pyNastran` target audience are users of Nastran and therefore are expected
to be familiar with the software.  This has greatly reduced the necessity of
documenting every variable exhaustively as users can easily reference existing
Nastran documentation. The BDF file has roughly 700 cards availble to a user
with 200 being currently supported by pyNastran.  The majority of the cards,
defined as separate Python classes, are not documented.  However, the Quick
Reference Guide (QRG) defines each input to the card.  A user with the QRG
should have little effort in understanding what the various objects do. 
However, for convenience, it's still good to document variables.

:mod:`pyNastran` target audience largely uses `MATLAB
<http://www.mathworks.com/products/matlab/>`_, a matrix based programming
language, and typically has little experience with general programming.  There
are also users that know Python, but have never used a class or a dictionary,
which makes an API seems complicated.  Additionally, while not difficult to
build, the `Cython <http://www.cython.org/>`_ OP4 reader, which requires a C
compiler, has been a source of frustration for many users.  That said, most
operations should be relatively easy to accomplish.

