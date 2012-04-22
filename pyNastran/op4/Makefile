# CFLAGS='-isysroot /Developer/SDKs/MacOSX10.6.sdk -arch i386'
# LDFLAGS='-arch i386'
# CFLAGS=
# LDFLAGS=

all:
	python setup.py build_ext --inplace
	# CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) python setup.py build_ext --inplace

clean:
	rm -rf build *.so op4.c libop4.o

test: ./test_loadop4.py
	./test_loadop4.py
