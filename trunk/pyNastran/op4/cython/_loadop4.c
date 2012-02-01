#include "loadop4.h"
#include <stdlib.h>

float *
load_K_matrix(char *filename) {

    /* we don't actually use the filename here... */

    float *k_matrix = NULL;
    int i;

    k_matrix = (float*)malloc(DIM0*DIM1*sizeof(float));
    if(!k_matrix)
        goto fail;

    /* simple initialization */
    for(i=0; i<DIM0*DIM1; i++)
        k_matrix[i] = (float)(i*i);

fail:
    /* if allocation failed we return the NULL value. */
    return k_matrix;
}
