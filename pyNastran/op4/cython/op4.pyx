"""
Based on cython_wrapper.pyx from https://gist.github.com/1249305 
by Gael Varoquaux
"""

# prototypes of the C function we are interested in calling
cdef extern from "libop4.c":
    float *op4_load(int size)
    int    Op4_scan(char *filename)
    int    op4_scan(char *filename  , #/* in                          */
                    int  *n_mat     , #/* out number of matrices      */
                    char  name[][9] , #/* out matrix names            */
                    int  *storage   , #/* out 0=dn; 1=sp1; 2=sp2      */
                    int  *nRow      , #/* out number of rows          */
                    int  *nCol      , #/* out number of columns       */
                    int  *nStr      , #/* out number of strings       */
                    int  *nNnz      , #/* out number of nonzero terms */
                    int  *Type      , #/* out 1=RS; 2=RD; 3=CS; 4=CD  */
                    int  *form      , #/* out matrix form 6=symm, etc */
                    int  *digits    , #/* out size of mantissa        */
                    long *offset      #/* out byte offset to matrix   */
                    )

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
import os

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# We need to build an array-wrapper class to deallocate our array when
# the Python object is deleted.

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array

        This cannot be done in the constructor as it must recieve C-level
        arguments.

        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data            

        """
        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_INT, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)

def File(char *filename, 
         char *mode    ):
    cdef int  Max_Matrices_Per_OP4 = 100
    cdef int  n_mat_in_file
#   cdef char name[Max_Matrices_Per_OP4][9]   # not allowed
    cdef char name[100][9]
    cdef int  storage[100]
    cdef int  nRow[100]   
    cdef int  nCol[100]   
    cdef int  nStr[100]   
    cdef int  nNnz[100]   
    cdef int  Type[100]   
    cdef int  form[100]   
    cdef int  digits[100]
    cdef long offset[100]
#   print('op4.File got [%s], [%c]' % (filename, mode))
    if not os.path.exists(filename):
        print('op4.File: no such file %s' % (filename))
        raise IOError
    elif not os.path.isfile(filename):
        print('op4.File: not a file %s' % (filename))
        raise IOError
    fh = { 'File' : filename }
    # Scan the file for number of matrices and headers for each matrix.
    try:
        rv = op4_scan(filename , # in                            
                 &n_mat_in_file, # out number of matrices        
                  name         , # out matrix names              
                  storage      , # out 0=dn; 1=sp1; 2=sp2        
                  nRow         , # out number of rows            
                  nCol         , # out number of columns         
                  nStr         , # out number of strings         
                  nNnz         , # out number of nonzero terms   
                  Type         , # out 1=RS; 2=RD; 3=CS; 4=CD    
                  form         , # out matrix form 6=symm, etc   
                  digits       , # out size of mantissa          
                  offset)        # out byte offset to matrix     

        if n_mat_in_file > Max_Matrices_Per_OP4:
            print('op4.File: too many matrices in %s' % (filename))
            raise IOError

        fh['nMat']   = n_mat_in_file
        fh['name']   = []
        fh['nRow']   = []
        fh['nCol']   = []
        fh['nStr']   = []
        fh['nNnz']   = []
        fh['type']   = []
        fh['form']   = []
        fh['digits'] = []
        fh['offset'] = []
        for i in range(n_mat_in_file):
            fh['name'].append( name[i] )
            fh['nRow'].append( nRow[i] )
            fh['nCol'].append( nCol[i] )
            fh['nStr'].append( nStr[i] )
            fh['nNnz'].append( nNnz[i] )
            fh['type'].append( Type[i] )
            fh['form'].append( form[i] )
            fh['digits'].append( digits[i] )
            fh['offset'].append( digits[i] )
    except:
        print('op4.File: failed to scan %s' % (filename))
        raise IOError

    return fh

def load(int size):
    """ Python binding of the 'compute' function in 'libop4.c' that does
        not copy the data allocated in C.
    """
    cdef float *array
    cdef np.ndarray ndarray
    # Call the C function
    array = op4_load(size)

    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(size, <void*> array) 
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)

    ndarray = np.reshape(ndarray, (2,5)) 

    return ndarray
