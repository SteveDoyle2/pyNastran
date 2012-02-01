# we cimport numpy.pxd.
cimport numpy as cnp

# we cimport loadop4.pxd.
cimport loadop4

# so that we can call free.
# cimport stdlib
cimport libc.stdlib

DEF ND = 2

# must call this at the toplevel to use
# PyArray_SimpleNewFromData inside read() 
cnp.import_array()

def read(filename):

    # the number of dimensions.
    cdef int nd = ND
    # the dims array.
    cdef int dims[ND]
    dims[0] = loadop4.DIM0; dims[1] = loadop4.DIM1

    # call load_K_matrix()
    cdef float *data = NULL
    data = loadop4.load_K_matrix(<char*>filename)
#   cdef float *data = loadop4.load_K_matrix(<char*>filename)

    # ensure the call succeeded.
    if data == NULL:
        raise MemoryError()

    # create a numpy array from the data pointer.
    arr = cnp.PyArray_SimpleNewFromData(nd, dims, cnp.NPY_FLOAT32, data) 

    # the 'arr' array does not call free() on the data when it is garbage
    # collected. We have to create a _dealloc_object() below and set it to the
    # base of the array.  The _dealloc_object is responsible for freeing the
    # data, which is done in the __dealloc__ call.
    do = _dealloc_object()
    do.ptr = <void*>data
    cnp.set_array_base(arr, do)
    return arr

cdef class _dealloc_object:

    # a pointer to the data
    cdef void *ptr

    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        print "deallocating"
        if self.ptr:
            libc.stdlib.free(self.ptr)
