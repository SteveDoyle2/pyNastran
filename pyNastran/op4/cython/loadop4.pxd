
cdef extern from "loadop4.h":
    float *load_K_matrix(char *filename)
    int DIM0
    int DIM1
