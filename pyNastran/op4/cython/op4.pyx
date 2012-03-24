"""
Based on cython_wrapper.pyx from https://gist.github.com/1249305 
by Gael Varoquaux
"""

# prototypes of the C function we are interested in calling

cdef extern from "stdio.h":
    cdef int SEEK_SET = 0
    ctypedef void* FILE
    int   fseek (FILE *stream, long int off, int whence)
    char *fgets (char *s, int n, FILE *stream)

cdef extern from "libop4.c":
    ctypedef struct str_t:
        int     len         #/* Number of terms in the string (a complex       */
                            #/* number counts as a single term).               */
        int     start_row   #/* Zero based, so first row has start_row = 0.    */
        int     N_idx       #/* Index into N[] to first numeric term for this  */
                            #/* string.                                        */
    float  *op4_load_S(FILE   *fp         ,
                       int     nRow       ,  #/* in  # rows    in matrix        */
                       int     nCol       ,  #/* in  # columns in matrix        */
                       char   *fmt_str    ,  #/* in  eg "%23le%23le%23le"       */
                       int     col_width  ,  #/* in  # characters in format str */
                       int     storage    ,  #/* in  0=dense  1,2=sparse  3=ccr */
                       int     complx     ,  #/* in  0=real   1=complex         */
                       int    *n_str      ,  #/* out # strings   (s_o) = 1,2    */
                       str_t  *str_data   ,  #/* out string data (s_o) = 1,2    */
                       int    *N_index    ,  #/* in/out          (s_o) = 1,2    */
                       )
    double *op4_load_D(FILE   *fp         ,
                       int     nRow       ,  #/* in  # rows    in matrix        */
                       int     nCol       ,  #/* in  # columns in matrix        */
                       char   *fmt_str    ,  #/* in  eg "%23le%23le%23le"       */
                       int     col_width  ,  #/* in  # characters in format str */
                       int     storage    ,  #/* in  0=dense  1,2=sparse  3=ccr */
                       int     complx     ,  #/* in  0=real   1=complex         */
                       int    *n_str      ,  #/* out # strings   (s_o) = 1,2    */
                       str_t  *str_data   ,  #/* out string data (s_o) = 1,2    */
                       int    *N_index    ,  #/* in/out          (s_o) = 1,2    */
                       )
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
    int  op4_filetype(char *filename)
    FILE* op4_open_r(char *filename, long offset)
    int op4_read_col_t(FILE   *fp         ,
                       int     c_in       ,  #/* in  requested column to read   */
                       int     nRow       ,  #/* in  # rows    in matrix        */
                       int     nCol       ,  #/* in  # columns in matrix        */
                       char   *fmt_str    ,  #/* in  eg "%23le%23le%23le"       */
                       int     col_width  ,  #/* in  # characters in format str */
                       int     storage    ,  #/* in  0=dense  1,2=sparse  3=ccr */
                       int     complx     ,  #/* in  0=real   1=complex         */
                       int    *n_str      ,  #/* out # strings   (s_o) = 1,2    */
                       str_t  *str_data   ,  #/* out string data (s_o) = 1,2    */
                       int    *N_index    ,  #/* in/out          (s_o) = 1,2    */
                       double *N             #/* out numeric data               */
                      )
    int  op4_read_col_b(FILE   *fp         ,
                        int     endian     ,  #/* in  0=native   1=flipped    */
                        int     c_in       ,  #/* in  requested column to read   */
                        int     nRow       ,  #/* in  # rows    in matrix        */
                        int     nCol       ,  #/* in  # columns in matrix        */
                        int     nType      ,  #/* in  1=RS 2=RD 3=CS 4=CD        */
                        int     storage    ,  #/* in  0=dn; 1=sp1; 2=sp2         */
                        int    *n_str      ,  #/* in/out idx str_data[] (s_o)=1,2*/
                        str_t  *S          ,  #/* out string data   (s_o)=1,2    */
                        int    *N_index    ,  #/* in/out idx N[]    (s_o)=1,2    */
                        double *N             #/* out numeric data               */
                       )

cdef int OP4_TXT_LINE_SIZE = 82

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

# types at http://cython.org/release/Cython-0.13/Cython/Includes/numpy.pxd
cdef class ArrayWrapper_RS:                                # {{{1
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
                                               np.NPY_FLOAT, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
# 1}}}
cdef class ArrayWrapper_RD:                                # {{{1
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
                                               np.NPY_DOUBLE, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
# 1}}}
cdef class ArrayWrapper_CS:                                # {{{1
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
                                               np.NPY_CFLOAT, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
# 1}}}
cdef class ArrayWrapper_CD:                                # {{{1
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
                                               np.NPY_CDOUBLE, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
# 1}}}
def File(char *filename,                                   # {{{1
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

        fh['nMat']    = n_mat_in_file
        fh['name']    = []
        fh['storage'] = []
        fh['nRow']    = []
        fh['nCol']    = []
        fh['nStr']    = []
        fh['nNnz']    = []
        fh['type']    = []
        fh['form']    = []
        fh['digits']  = []
        fh['offset']  = []
        for i in range(n_mat_in_file):
            fh['name'].append(    name[i]    )
            fh['storage'].append( storage[i] )
            fh['nRow'].append(    nRow[i]    )
            fh['nCol'].append(    nCol[i]    )
            fh['nStr'].append(    nStr[i]    )
            fh['nNnz'].append(    nNnz[i]    )
            fh['type'].append(    Type[i]    )
            fh['form'].append(    form[i]    )
            fh['digits'].append(  digits[i]  )
            fh['offset'].append(  offset[i]  )
    except:
        print('op4.File: failed to scan %s' % (filename))
        raise IOError

    return fh
# 1}}}
def Load(op4_handle,    # in, as created by File()           {{{1
         int n_mat=1 ,  # in, number of matrices to return
         int n_skip=0): # in, number of matrices to skip before loading
    cdef FILE   *fp
    cdef char    line[83]        # not liking line[OP4_TXT_LINE_SIZE+1]
    cdef int     size
    cdef int    *unused
    cdef str_t  *unused_s
    cdef float  *array_RS
    cdef float  *array_CS
    cdef double *array_RD
    cdef double *array_CD
    cdef np.ndarray ndarray
    """ 
    Skip over the first n_skip matrices then load n_mat matrices from the 
    previously opened op4 file.
    """
    if n_mat < 1 or \
       n_mat > len(op4_handle['digits']):
        print('op4.Load: 1 <= n_mat <= %d' % (len(op4_handle['digits'])))
        raise IOError

    if n_skip < 0 or \
       n_skip > len(op4_handle['digits'])-1:
        print('op4.Load: 0 <= n_skip <= %d' % (len(op4_handle['digits']) - 1))
        raise IOError

    print('Will return %d matrices after skipping %d' % (n_mat, n_skip))
    filetype = op4_filetype(op4_handle['File'])
    print('File type of %s is %d' % (op4_handle['File'], filetype))
    for i in range(n_mat):
        offset = i + n_skip
        print('Fetching %d. %-8s  %5d x %5d' % (
            offset + 1, 
            op4_handle['name'][offset],
            op4_handle['nRow'][offset],
            op4_handle['nCol'][offset], ))

    # print(op4_handle)

    if os.path.basename(op4_handle['File']) != 'mat_t_dn.op4': return

    # print('Attempting file open of %s' % (op4_handle['File']))
    fp = op4_open_r(op4_handle['File'], 0)
    if fp == NULL:
        print('unable to open %s' % (op4_handle['File']))
        raise IOError

    print(' open successful')
    for i in range(n_skip, n_skip+n_mat):
        complx = False
        if op4_handle['type'][i] > 2: complx = True
        # print('%s complx=%d' % (op4_handle['name'][i], complx))
        if op4_handle['storage'][i]:
            print('%s is sparse, skipping for now' % (op4_handle['name'][i]))
        else:

            # now read the file contents
            fseek(fp, op4_handle['offset'][i], SEEK_SET)
            if filetype == 1: # text
                fgets(line, OP4_TXT_LINE_SIZE, fp)   # skip the header line
                # create the scanf format string
                col_width = op4_handle['digits'][i] + 7
                fmt_one   = '%%%dle' % col_width
                nNum_cols = 80 // col_width
                fmt_str   = fmt_one*nNum_cols
                # print('fmt_str=[%s]' % (fmt_str))

            if filetype == 1: # text
                if op4_handle['storage'][i]:     # sparse
                    pass
                else:                            # dense
                    size = op4_handle['nRow'][i] * op4_handle['nCol'][i]
                    if   op4_handle['type'][i] == 1: 
                        array_RS = op4_load_S(fp, 
                                              op4_handle['nRow'][i], 
                                              op4_handle['nCol'][i], 
                                              fmt_str, 
                                              col_width,
                                              op4_handle['storage'][i], # in 0=dn  1,2=sp1,2  3=ccr  
                                              complx    , # in  0=real   1=complex     
                                              unused    , # in/out index m.S (if sp1,2)
                                              unused_s  , # out string data (if sp1,2) 
                                              unused    ) # in/out index m.N (sp 1,2)  

                        array_wrapper_RS = ArrayWrapper_RS()
                        array_wrapper_RS.set_data(size, <void*> array_RS) 
                        ndarray = np.array(array_wrapper_RS, copy=False)
                        ndarray.base = <PyObject*> array_wrapper_RS
                        Py_INCREF(array_wrapper_RS)

                    elif op4_handle['type'][i] == 2: 
                        array_RD = op4_load_D(fp, 
                                              op4_handle['nRow'][i], 
                                              op4_handle['nCol'][i], 
                                              fmt_str, 
                                              col_width,
                                              op4_handle['storage'][i], # in 0=dn  1,2=sp1,2  3=ccr  
                                              complx    , # in  0=real   1=complex     
                                              unused    , # in/out index m.S (if sp1,2)
                                              unused_s  , # out string data (if sp1,2) 
                                              unused    ) # in/out index m.N (sp 1,2)  
                        array_wrapper_RD = ArrayWrapper_RD()
                        array_wrapper_RD.set_data(size, <void*> array_RD) 
                        ndarray = np.array(array_wrapper_RD, copy=False)
                        ndarray.base = <PyObject*> array_wrapper_RD
                        Py_INCREF(array_wrapper_RD)

                    elif op4_handle['type'][i] == 3: 
                        array_CS = op4_load_S(fp, 
                                              op4_handle['nRow'][i], 
                                              op4_handle['nCol'][i], 
                                              fmt_str, 
                                              col_width,
                                              op4_handle['storage'][i], # in 0=dn  1,2=sp1,2  3=ccr  
                                              complx    , # in  0=real   1=complex     
                                              unused    , # in/out index m.S (if sp1,2)
                                              unused_s  , # out string data (if sp1,2) 
                                              unused    ) # in/out index m.N (sp 1,2)  
                        array_wrapper_CS = ArrayWrapper_CS()
                        array_wrapper_CS.set_data(size, <void*> array_CS) 
                        ndarray = np.array(array_wrapper_CS, copy=False)
                        ndarray.base = <PyObject*> array_wrapper_CS
                        Py_INCREF(array_wrapper_CS)

                    elif op4_handle['type'][i] == 4: 
                        array_CD = op4_load_D(fp, 
                                              op4_handle['nRow'][i], 
                                              op4_handle['nCol'][i], 
                                              fmt_str, 
                                              col_width,
                                              op4_handle['storage'][i], # in 0=dn  1,2=sp1,2  3=ccr  
                                              complx    , # in  0=real   1=complex     
                                              unused    , # in/out index m.S (if sp1,2)
                                              unused_s  , # out string data (if sp1,2) 
                                              unused    ) # in/out index m.N (sp 1,2)  
                        array_wrapper_CD = ArrayWrapper_CD()
                        array_wrapper_CD.set_data(size, <void*> array_CD) 
                        ndarray = np.array(array_wrapper_CD, copy=False)
                        ndarray.base = <PyObject*> array_wrapper_CD
                        Py_INCREF(array_wrapper_CD)

                    else:
                        print('op4.Load type failure  %s: %d' % (
                            op4_handle['name'][i], op4_handle['type'][i]))
                        raise IOError

                    ndarray = np.reshape(ndarray, (op4_handle['nCol'][i], op4_handle['nRow'][i])).T

                    return ndarray

            else:             # binary
                if op4_handle['storage'][i]:     # sparse
                    pass
                else:                            # dense
                    pass
# 1}}}
