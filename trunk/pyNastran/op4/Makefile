# CFLAGS='-isysroot /Developer/SDKs/MacOSX10.6.sdk -arch i386'
# LDFLAGS='-arch i386'
CFLAGS= -g
LDFLAGS= -g

all: op4.pyx libop4.c op4.h
	LDFLAGS=$(LDFLAGS) CFLAGS=$(CFLAGS) python setup.py build_ext --inplace
clean:
	rm -rf build *.so op4.c libop4.o

test: ./test_loadop4.py
	./test_loadop4.py


# if that doesn't work try one of these...
# Windows
# -------
#    1.  python setup.py build_ext --inplace
#    2.  Install MinGW from http://www.mingw.org/
#        python setup.py build_ext --inplace --compiler=mingw32
#
# Linux
# -------
#    Google :)
#
# If you have other suggestions, please send an email.

