# CFLAGS='-isysroot /Developer/SDKs/MacOSX10.6.sdk -arch i386'
# LDFLAGS='-arch i386'
CFLAGS=
LDFLAGS=

all:
	CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) python setup.py build_ext --inplace

clean:
	rm -rf build *.so loadop4.c
