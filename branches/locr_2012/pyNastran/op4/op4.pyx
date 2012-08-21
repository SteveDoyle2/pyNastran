# {{{1 GNU Lesser General Public License
# 
# Copyright (C) 1999-2012  <al.danial@gmail.com>
# 
# This program is free software; you can redistribute it and/or   
# modify it under the terms of the GNU Lesser General Public      
# License as published by the Free Software Foundation;           
# version 3 of the License.                                       
#                                                                 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of  
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    
# GNU Lesser General Public License for more details.             
#                                                                 
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the           
# Free Software Foundation, Inc.,                                 
# 675 Mass Ave, Cambridge, MA 02139, USA.                         
# 1}}}

""" 
Load matrices from NASTRAN Output 4 (.op4) files.
Example:

  >>> fh = op4.File('sol103.op4', 'r')

    Creates a read object associated with the file sol103.op4.

  >>> fh.nmat

    Returns the number of matrices contained in the file.

  >>> fh.print_headers()

    Prints detailed information on the matrices stored in the file.

  >>> K, M = fh.Load(skip=1, nmat=2)

    Loads the second and third matrices from the file sol103.op4 into K and M.

Limitations (these will be resolved in future releases):
    1. only handles dense matrices in text and binary op4 files
    2. unable to byte-swap files created on a different-endian machine
    3. unable to save Python arrays to new op4 files
"""

# Early version of op4.pyx was based on cython_wrapper.pyx from
# https://gist.github.com/1249305 by Gael Varoquaux

cdef extern from "stdio.h":
    cdef int SEEK_SET = 0
    ctypedef void* FILE
    int   fseek (FILE *stream, long int off, int whence)
    char *fgets (char *s, int n, FILE *stream)
    cdef int fclose(FILE *fp)
    FILE *fopen(char *path, char *mode)

cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)

cdef extern from "string.h":
    void* memset(void *dest, int c, size_t n)
#   void* memcpy(void *dest, void *src, size_t n)

# Import the Python-level symbols of numpy
import numpy as np
from   scipy import sparse

# Import the C-level symbols of numpy
cimport numpy as np

cdef extern from "libop4.c":
    ctypedef np.int_t      itype_t   # can't get this to do anything useful
    ctypedef struct str_t:
        int     len         # Number of terms in the string (a complex     
                            # number counts as a single term).             
        int     start_row   # Zero based, so first row has start_row = 0.  
        int     N_idx       # Index into N[] to first numeric term for this
                            # string.                                      
    ctypedef np.float64_t  dtype_t
    float  *op4_load_S(FILE   *fp         ,
                       int     filetype   ,  # in  1=text, other is binary  
                       int     nRow       ,  # in  # rows    in matrix       
                       int     nCol       ,  # in  # columns in matrix       
                       char   *fmt_str    ,  # in  eg "%23le%23le%23le"      
                       int     col_width  ,  # in  # characters in format str
                       int     storage    ,  # in  0=dense  1,2=sparse  3=ccr
                       int     complx     ,  # in  0=real   1=complex        
                       int     n_Nnz      ,  # in  number of nonzero terms   
                       double *I_coo      ,  # out sparse row ind            
                       double *J_coo      ,  # out sparse col ind            
                       )
    double *op4_load_D(FILE   *fp         ,
                       int     filetype   ,  # in  1=text, other is binary   
                       int     nRow       ,  # in  # rows    in matrix       
                       int     nCol       ,  # in  # columns in matrix       
                       char   *fmt_str    ,  # in  eg "%23le%23le%23le"      
                       int     col_width  ,  # in  # characters in format str
                       int     storage    ,  # in  0=dense  1,2=sparse  3=ccr
                       int     complx     ,  # in  0=real   1=complex        
                       int     n_Nnz      ,  # in  number of nonzero terms   
                       double *I_coo      ,  # out sparse row ind            
                       double *J_coo      ,  # out sparse col ind            
                       )
    int    op4_scan(char *filename  , # in                         
                    int  *nmat      , # out number of matrices     
                    char  name[][9] , # out matrix names           
                    int  *storage   , # out 0=dn; 1=sp1; 2=sp2     
                    int  *nRow      , # out number of rows         
                    int  *nCol      , # out number of columns      
                    int  *nStr      , # out number of strings      
                    int  *nNnz      , # out number of nonzero terms
                    int  *Type      , # out 1=RS; 2=RD; 3=CS; 4=CD 
                    int  *form      , # out matrix form 6=symm, etc
                    int  *digits    , # out size of mantissa       
                    long *offset      # out byte offset to matrix  
                    )
    int  op4_filetype(char *filename)
    FILE* op4_open_r(char *filename, long offset)
    int op4_read_col_t(FILE   *fp         ,
                       int     c_in       ,  # in  requested column to read  
                       int     nRow       ,  # in  # rows    in matrix       
                       int     nCol       ,  # in  # columns in matrix       
                       char   *fmt_str    ,  # in  eg "%23le%23le%23le"      
                       int     col_width  ,  # in  # characters in format str
                       int     storage    ,  # in  0=dense  1,2=sparse  3=ccr
                       int     complx     ,  # in  0=real   1=complex        
                       int    *n_str      ,  # out # strings   (s_o) = 1,2   
                       str_t  *str_data   ,  # out string data (s_o) = 1,2   
                       int    *N_index    ,  # in/out          (s_o) = 1,2   
                       double *N             # out numeric data              
                      )
    int  op4_read_col_b(FILE   *fp         ,
                        int     endian     , # in  0=native   1=flipped       
                        int     c_in       , # in  requested column to read   
                        int     nRow       , # in  # rows    in matrix        
                        int     nCol       , # in  # columns in matrix        
                        int     nType      , # in  1=RS 2=RD 3=CS 4=CD        
                        int     storage    , # in  0=dn; 1=sp1; 2=sp2         
                        int    *n_str      , # in/out idx str_data[] (s_o)=1,2
                        str_t  *S          , # out string data   (s_o)=1,2    
                        int    *N_index    , # in/out idx N[]    (s_o)=1,2    
                        double *N            # out numeric data               
                       )
    int  op4_wrt_header(FILE   *fp         ,
                        int     endian     , # in  0=native   1=flipped    
                        char   *name       , # in  matrix name             
                        int     nRow       , # in  # rows    in matrix     
                        int     nCol       , # in  # columns in matrix     
                        int     nType      , # in  1=RS 2=RD 3=CS 4=CD     
                        int     Form       , # in  1=rect; 2=square        
                        int     sparse     , # in  1=sp2; 0=dn             
                        int     digits       # in -1=flipped endian binary 
                                             #     0=native  endian binary 
                                             #    >2=text; number of DIGITS
                       )
    int  op4_wrt_col_dn(FILE   *fp         ,
                        int     column     , # first column is 0           
                        int     nRows      , # number of rows              
                        int     nCols      , # number of columns           
                        double *A          , # array of terms to write     
                        int     complx     , # 1=complex  0=real           
                        int     digits       # in -1=flipped endian binary 
                                             #     0=native  endian binary 
                                             #    >2=text; number of DIGITS
                        )
    int  op4_wrt_trailer(FILE *fp          ,
                         int   column      , # first column is 0
                         int   digits        # -1=flipped; 0=native; >0=digits
                        )

cdef int OP4_TXT_LINE_SIZE = 82
DTYPE = np.float64
FTYPE = np.float32
ctypedef np.float64_t DTYPE_t

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
import os

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
class OP4:                                                 # {{{1
    type_map = { 1 : 'RS' ,   # real    single precision
                 2 : 'RD' ,   # real    double precision
                 3 : 'CS' ,   # complex single precision
                 4 : 'CD' ,}  # complex double precision
    def __init__(self          ,                # {{{2
                 char *filename,
                 mode          ):
        cdef int  Max_Matrices_Per_OP4 = 100
        cdef int  nmat_in_file
#       cdef char name[Max_Matrices_Per_OP4][9]   # not allowed
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
        cdef np.ndarray ndfh
#       cdef FILE *fh[1]
        cdef double fh[1]

#       print('op4.File got [%s], [%c]' % (filename, mode))
        if mode == 'r' or mode == 'R':
            if not os.path.exists(filename):
                print('op4.File: no such file %s' % (filename))
                raise IOError
            elif not os.path.isfile(filename):
                print('op4.File: not a file %s' % (filename))
                raise IOError
        self.filename = os.path.abspath(filename)
        if mode in ['w', 'W', 'a', 'A']:
#           self.fh = fopen(filename, "wb")
#           # self.fh_txt = fopen(filename, "w")
#           fh[0] = <void *> fopen(filename, "wb")
#           fh_wrapper_RD = ArrayWrapper_RD()
#           fh_wrapper_RD.set_data(1, <void*> fh)   # crash
#           ndfh = np.array(fh_wrapper_RD, copy=False)
#           ndfh.base = <PyObject*> fh_wrapper_RD
#           Py_INCREF(fh_wrapper_RD)
#           self.fh = fh_wrapper_RD
            return

        # Scan the file for number of matrices and headers for each matrix.
        try:
            rv = op4_scan(filename , # in                            
                     &nmat_in_file , # out number of matrices        
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

            if nmat_in_file > Max_Matrices_Per_OP4:
                print('op4.File: too many matrices in %s; max allowed is %d' % (
                    filename, Max_Matrices_Per_OP4))
                raise IOError

            self.nmat      = nmat_in_file
            self.name      = []
            self.storage   = []
            self.nrow      = []
            self.ncol      = []
            self.nstr      = []
            self.nnz       = []
            self.type      = []
            self.form      = []
            self.digits    = []
            self.offset    = []
            for i in range(nmat_in_file):
                self.name.append(    name[i]    )
                self.storage.append( storage[i] )
                self.nrow.append(    nRow[i]    )
                self.ncol.append(    nCol[i]    )
                self.nstr.append(    nStr[i]    )
                self.nnz.append(     nNnz[i]    )
                self.type.append(    Type[i]    )
                self.form.append(    form[i]    )
                self.digits.append(  digits[i]  )
                self.offset.append(  offset[i]  )
        except:
            print('op4.File: failed to scan %s' % (filename))
            raise IOError
    # 2}}}
    def print_instance(self, i):                # {{{2 
        print('%2d. %-8s %5d %5d %8d %8d %2s %2d %2d %9d' % (
               i+1          ,
               self.name[i],
               self.nrow[i],
               self.ncol[i],
               self.nstr[i],
               self.nnz[i],
               self.type_map[self.type[i]],
               self.form[i],
               self.digits[i],
               self.offset[i],))
    # 2}}}
    def print_header(self, index=None):         # {{{2 

        print('    %-8s %5s %5s %8s %8s %2s %2s %2s %9s' % (
              'Name', 'nRow', 'nCol', 'nStr', 'nNnz', 'T', 'Fr', 'Dg', 'Offset'))
        if index is None:
            for i in range(self.nmat):
                self.print_instance(i)
        else:
            self.print_instance(index)
    # 2}}}
    def Load(self ,                              # {{{2
             int nmat=1 ,  # in, number of matrices to return
             int skip=0):  # in, number of matrices to skip before loading
        cdef FILE   *fp
        cdef char    line[83]        # not liking line[OP4_TXT_LINE_SIZE+1]
        cdef int     size
        cdef int     filetype
        cdef double *unused
        cdef int     n_str
        cdef str_t  *unused_s
        cdef str_t  *str_data
        cdef int     n_nnz = self.nnz[skip]
        cdef int    *N_index
        cdef float  *array_RS
        cdef float  *array_CS
        cdef double *array_RD
        cdef double *array_CD
        cdef double *I_coo_C
        cdef double *J_coo_C
        cdef np.ndarray ndarray
        cdef np.ndarray ndarray_I
        cdef np.ndarray ndarray_J
        All_Matrices = []

        if nmat < 1 or \
           nmat > len(self.digits):
            print('op4.Load: 1 <= nmat <= %d' % (len(self.digits)))
            raise IOError

        if skip < 0 or \
           skip > len(self.digits)-1:
            print('op4.Load: 0 <= skip <= %d' % (len(self.digits) - 1))
            raise IOError

        filetype = op4_filetype(self.filename)
        if False: # set to True for debug output
            print('Will return %d matrices after skipping %d' % (nmat, skip))
            print('File type of %s is %d' % (self.filename, filetype))
            for i in range(nmat):
                offset = i + skip
                print_header(self, index=offset)
                print('Fetching %d. %-8s  %5d x %5d from %s' % (
                        offset + 1, 
                        self.name[offset],
                        self.nrow[offset],
                        self.ncol[offset], 
                        self.filename))

#       print('Attempting file open of %s' % (self.filename))
        fp = op4_open_r(self.filename, 0)
        if fp == NULL:
            print('unable to open %s' % (self.filename))
            raise IOError

        # print(' open successful')
        for i in range(skip, skip+nmat):
            complx = False
            if self.type[i] > 2: complx = True
            # print('%s complx=%d' % (self.name[i], complx))

            # now read the file contents
            fseek(fp, self.offset[i], SEEK_SET)
            if filetype == 1: # text
                fgets(line, OP4_TXT_LINE_SIZE, fp)   # skip the header line
                # create the scanf format string
                col_width = self.digits[i] + 7
                fmt_one   = '%%%dle' % col_width
                nNum_cols = 80 // col_width
                fmt_str   = fmt_one*nNum_cols
                # print('fmt_str=[%s]' % (fmt_str))
            else: # binary
                col_width = 0
                fmt_str   = ''

            if self.storage[i]:          # sparse

                I_coo_C = <DTYPE_t*>malloc(sizeof(double) * self.nnz[i])
                memset(I_coo_C, 0,         sizeof(double) * self.nnz[i])
                J_coo_C = <DTYPE_t*>malloc(sizeof(double) * self.nnz[i])
                memset(J_coo_C, 0,         sizeof(double) * self.nnz[i])

                if   self.type[i] == 1:  # sparse real single precision
                    array_RS = op4_load_S(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   , 
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          I_coo_C   , # out sparse row ind
                                          J_coo_C   ) # out sparse col ind

                    array_wrapper_RS = ArrayWrapper_RS()
                    array_wrapper_RS.set_data(self.nnz[i], <void*> array_RS) 
                    ndarray = np.array(array_wrapper_RS, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_RS
                    Py_INCREF(array_wrapper_RS)

                elif self.type[i] == 2:   # sparse real double precision
                    array_RD = op4_load_D(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   , 
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          I_coo_C   , # out sparse row ind
                                          J_coo_C   ) # out sparse col ind

                    array_wrapper_RD = ArrayWrapper_RD()
                    array_wrapper_RD.set_data(self.nnz[i], <void*> array_RD) 
                    ndarray = np.array(array_wrapper_RD, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_RD
                    Py_INCREF(array_wrapper_RD)

                elif self.type[i] == 3:   # dense complex single precision
                    array_CS = op4_load_S(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   , 
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          I_coo_C   , # out sparse row ind
                                          J_coo_C   ) # out sparse col ind

                    array_wrapper_CS = ArrayWrapper_CS()
                    array_wrapper_CS.set_data(self.nnz[i], <void*> array_CS) 
                    ndarray = np.array(array_wrapper_CS, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_CS
                    Py_INCREF(array_wrapper_CS)

                elif self.type[i] == 4:   # dense complex double precision
                    array_CD = op4_load_D(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   , 
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          I_coo_C   , # out sparse row ind
                                          J_coo_C   ) # out sparse col ind

                    array_wrapper_CD = ArrayWrapper_CD()
                    array_wrapper_CD.set_data(self.nnz[i], <void*> array_CD) 
                    ndarray = np.array(array_wrapper_CD, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_CD
                    Py_INCREF(array_wrapper_CD)
                else:
                    print('op4.Load sparse type failure  %s: %d' % (
                        self.name[i], self.type[i]))
                    raise IOError
###############     print("matrix %d is sparse and not single prec, skipping" % i)

                array_wrapper_I  = ArrayWrapper_RD()
                array_wrapper_I.set_data(self.nnz[i], <void*> I_coo_C) 
                ndarray_I = np.array(array_wrapper_I , copy=False)
                ndarray_I.base = <PyObject*> array_wrapper_I
                Py_INCREF(array_wrapper_I)

                array_wrapper_J  = ArrayWrapper_RD()
                array_wrapper_J.set_data(self.nnz[i], <void*> J_coo_C) 
                ndarray_J = np.array(array_wrapper_J , copy=False)
                ndarray_J.base = <PyObject*> array_wrapper_J
                Py_INCREF(array_wrapper_J)

                All_Matrices.append( sparse.coo_matrix(
                    (ndarray, (ndarray_I ,ndarray_J)), 
                     shape=(self.nrow[i], self.ncol[i])) )

            else:                        # dense
                size = self.nrow[i] * self.ncol[i]
                if   self.type[i] == 1:  # dense real single precision
                    array_RS = op4_load_S(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   , 
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          unused    , # out sparse row ind
                                          unused    ) # out sparse col ind

                    array_wrapper_RS = ArrayWrapper_RS()
                    array_wrapper_RS.set_data(size, <void*> array_RS) 
                    ndarray = np.array(array_wrapper_RS, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_RS
                    Py_INCREF(array_wrapper_RS)

                elif self.type[i] == 2:   # dense real double precision
                    array_RD = op4_load_D(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   ,
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          unused    , # out sparse row ind
                                          unused    ) # out sparse col ind
                    array_wrapper_RD = ArrayWrapper_RD()
                    array_wrapper_RD.set_data(size, <void*> array_RD) 
                    ndarray = np.array(array_wrapper_RD, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_RD
                    Py_INCREF(array_wrapper_RD)

                elif self.type[i] == 3:   # dense complex single precision
                    array_CS = op4_load_S(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   ,
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          unused    , # out sparse row ind
                                          unused    ) # out sparse col ind
                    array_wrapper_CS = ArrayWrapper_CS()
                    array_wrapper_CS.set_data(size, <void*> array_CS) 
                    ndarray = np.array(array_wrapper_CS, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_CS
                    Py_INCREF(array_wrapper_CS)

                elif self.type[i] == 4:   # dense complex double precision
                    array_CD = op4_load_D(fp        ,
                                          filetype  , # 1 = text
                                          self.nrow[i], 
                                          self.ncol[i], 
                                          fmt_str   ,
                                          col_width ,
                                          self.storage[i], # in 0=dn  1,2=sp1,2  3=ccr  
                                          complx    , # in  0=real   1=complex     
                                          self.nnz[i],# in number of nonzero terms
                                          unused    , # out sparse row ind
                                          unused    ) # out sparse col ind
                    array_wrapper_CD = ArrayWrapper_CD()
                    array_wrapper_CD.set_data(size, <void*> array_CD) 
                    ndarray = np.array(array_wrapper_CD, copy=False)
                    ndarray.base = <PyObject*> array_wrapper_CD
                    Py_INCREF(array_wrapper_CD)

                else:
                    print('op4.Load dense type failure  %s: %d' % (
                        self.name[i], self.type[i]))
                    raise IOError

                ndarray = np.reshape(ndarray, (self.ncol[i], self.nrow[i])).T

                All_Matrices.append(ndarray)

        if len(All_Matrices) == 1:
            return All_Matrices[0]
        else:
            return All_Matrices

    # 2}}}
# 1}}}
def Save(                                                  # {{{1
         filename ,
         *args    ,  # unnamed matrices
         **kwargs ): # named matrices and other options like digits
    """

    cop4.Save('kd.op4', Kxx=K, D)
        Write matrices K and D to x.op4 with names 
        "Kxx" and "M0000001" using native binary format.
     
    cop4.Save('kd.op4', Kxx=K, Dxx=D, digits=5)
        Write matrices K and D to x.op4 with names 
        "Kxx" and "Dxx" as a text file with 5 mantissa digits.
     
    """

    cdef FILE* fp
    cdef int endian = 0
    cdef int sparse = 0
    cdef double *one_column

    if 'digits' in kwargs: 
        digits = kwargs['digits']
        if digits and digits <  2: digits =  2
        if digits and digits > 26: digits = 26
        del kwargs['digits']
    else:
        digits = 0

#   print('write to %s with digits=%d' % (     filename, digits))

    # assign a name to each unnamed matrix: M0000001 through M9999999
    for mat in args:
        i = 1
        while True:
            name = 'M%07d' % i
            if name not in kwargs: break
            i += 1
        kwargs[name] = mat

    if digits: fp = fopen(filename, 'w' )  # text
    else:      fp = fopen(filename, 'wb')  # binary

    for name in kwargs:  # loop over the matrices

        try:
#           print('a %s' % name)
            if sparse.issparse(kwargs[name]): 
#               print('b')
                sparse = 1
            else:                             
#               print('c')
                sparse = 0
        except:
#           print('d')
            sparse = 0

        if (not sparse) and (type(kwargs[name]) is not np.ndarray):
            print('OP4.Save skipping %s, wrong type (%s)' % (
                name, type(kwargs[name])))
            continue
        if   kwargs[name].ndim == 1:
            nR, nC = kwargs[name].shape , 1
        elif kwargs[name].ndim == 2:
            nR, nC = kwargs[name].shape
        else:
            print('OP4.Save skipping %s, not 1D or 2D' % (name))
            continue

        if   kwargs[name].dtype == 'float32':
            op4_type   = 1
            op4_complx = 0
            npt        = 1
        elif kwargs[name].dtype == 'float64':
            op4_type   = 2
            op4_complx = 0
            npt        = 1
        elif kwargs[name].dtype == 'complex64':
            op4_type   = 3
            op4_complx = 1
            npt        = 2
        elif kwargs[name].dtype == 'complex128':
            op4_type   = 4
            op4_complx = 1
            npt        = 2
        else: # kwargs[name].dtype is 'complex128':
            print('   type %s not supported' % (kwargs[name].dtype))
            raise IOError

        if nR == nC: op4_form = 1   # square
        else:        op4_form = 2   # rectangular

        # text with few digits; demote to single precision
        if (op4_type in [2,4]) and (2 <= digits <= 9): op4_type -= 1

        op4_wrt_header(fp, endian, name, nR, nC, op4_type, 
                       op4_form, sparse, digits)

        one_column = <DTYPE_t*>malloc(sizeof(double) * nR * npt)

        # write each column
        for c in xrange(nC):
            if sparse:
                print('sparse column %d' % c)
                print kwargs[name].get_col(c)
            else:  # dense
                # probably a better way to make the column copies below
                if   op4_type in [1,2]:
                    for r in xrange(nR): one_column[r] = kwargs[name][r,c]
                else:
                    for r in xrange(nR): one_column[2*r], one_column[2*r+1] = \
                            kwargs[name][r,c].real, kwargs[name][r,c].imag

                op4_wrt_col_dn(fp, c, nR, nC, one_column, op4_complx, digits)

        free(one_column)

    fclose(fp)

    return
# 1}}}
