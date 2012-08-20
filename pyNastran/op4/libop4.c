/* {{{1 GNU Lesser General Public License

Copyright (C) 1999-2012  <al.danial@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation;
version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the
Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
1}}} */
#ifndef __STRICT_ANSI__
   #define __STRICT_ANSI__
#endif
#define _XOPEN_SOURCE 500 /* snprintf */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>      /* fabs */
#include <ctype.h>     /* isdigit, isalpha */
#include "op4.h"


#define  MAX(a,b)              ((a) >  (b)) ? (a)  : (b)
void get_file_suffix(char *file, char *ext);
/* constants {{{1
   would like to have the constant arrays below defined in op4.h
   but they cause 'multiple definition' errors
 */

const char op4_store_str[3][3] = { {'d',  'n',  0},
                                   {'s',  '1',  0},
                                   {'s',  '2',  0}
                           };
const int op4_words_per_term[5] = {-1, 1, 2, 2, 4};
const int Nm_per_term[2]  = {1,  /* numbers per term[real data type]    = 1  */
                             2}; /* numbers per term[complex data type] = 2  */
const char op4_type_str[4][3] = { {'R',  'S',  0},
                                  {'R',  'D',  0},
                                  {'C',  'S',  0},
                                  {'C',  'D',  0}
                           };
 /* 1}}} */
float  *op4_load_S(FILE   *fp         ,  /* {{{1 */
                   int     filetype   ,  /* in  1=text, other is binary    */
                   int     nRow       ,  /* in  # rows    in matrix        */
                   int     nCol       ,  /* in  # columns in matrix        */
                   char   *fmt_str    ,  /* in  eg "%23le%23le%23le"       */
                   int     col_width  ,  /* in  # characters in format str */
                   int     storage    ,  /* in  0=dense  1,2=sparse  3=ccr */
                   int     complx     ,  /* in  0=real   1=complex         */
                   int     n_Nnz      ,  /* in  number of nonzero terms    */
                   double *I_coo      ,  /* out: sparse row ind            */
                   double *J_coo      )  /* out: sparse col ind            */
{
    int     DEBUG = 0;
    float  *array;
    double *column;
    int     r, c, s, i, j, size, NPT, nnz, nType, s_ptr, n_ptr;
    str_t  *str_data;

    int     endian = 0;      /* FIX THIS */

    if (DEBUG) { printf("74 -> op4_load_S\n"); fflush(stdout); }
    i     = 0;
    nnz   = 0;
    size  = nRow*nCol;  /* only used for dense */
    NPT   = 1;
    nType = 1;
    if (complx) {
        NPT   *= 2;
        nType += 2;
    }
    /* op4_read_col_text() only populates arrays of type double.
     * To minimize memory use, only allocate memory for a single
     * column of type double and copy it to the float array each
     * pass through the columns.
     */
    if (!storage) { /* dense; only need one column at a time */
        array    = malloc(sizeof(float )*size*NPT);
        column   = malloc(sizeof(double)*nRow*NPT);
    } else {        /* sparse; must malloc all numeric data */
        array    = malloc(sizeof(float )*n_Nnz*NPT);
        column   = malloc(sizeof(double)*n_Nnz*NPT);
        str_data = malloc(sizeof(str_t)*(nRow+1)/2);  /* max # strings in a row is nR/2 + 1*/
    }
    if (DEBUG) { printf(" 97 -> op4_load_S after malloc\n"); fflush(stdout); }
    for (c = 0; c < nCol; c++) {
        if (DEBUG) { printf(" 99 - op4_load_S column %d\n", c); fflush(stdout); }
        if (storage) { /* sparse */
            s_ptr                 = 0;
            n_ptr                 = 0;
            str_data[0].start_row = 0;
            str_data[0].len       = 0;
            str_data[0].N_idx     = 0;
        }
        if (filetype == 1) {
            nnz = op4_read_col_text(fp, c+1, nRow, nCol, fmt_str, col_width,
                                   storage   ,  /* in 0=dn  1,2=sp1,2  3=ccr  */
                                   complx    ,  /* in  0=real   1=complex     */
                                  &s_ptr     ,  /* in/out index m.S (if sp1,2)*/
                                   str_data  ,  /* out string data (if sp1,2) */
                                  &n_ptr     ,  /* in/out index m.N (sp 1,2)  */
                                   column    ); /* out numeric data */
        } else {
            nnz = op4_read_col_binary(  fp        ,
                                   endian    ,  /* in  0=native   1=flipped    */
                                   c+1       ,  /* in  requested column to read   */
                                   nRow      ,  /* in  # rows    in matrix        */
                                   nCol      ,  /* in  # columns in matrix        */
                                   nType     ,  /* in  1=RS 2=RD 3=CS 4=CD        */
                                   storage   ,  /* in  0=dn; 1=sp1; 2=sp2         */
                                  &s_ptr     ,  /* in/out idx str_data[] (s_o)=1,2*/
                                   str_data  ,  /* out string data   (s_o)=1,2    */
                                  &n_ptr     ,  /* in/out idx N[]    (s_o)=1,2    */
                                   column    ); /* out numeric data               */
        }
        if (DEBUG) { printf("128 - op4_load_S column %d populating array\n", c); fflush(stdout); }
        if (storage) { /* sparse */
             for (s = 0; s < s_ptr; s++) {
                if (DEBUG) { printf("131 - op4_load_S    str_data[%d].len=%d  .start_row=%d  .N_idx=%d\n",
                        s, str_data[s].len, str_data[s].start_row, str_data[s].N_idx); fflush(stdout); }
                for (j = 0; j < str_data[s].len; j++) {
                     I_coo[i] = str_data[s].start_row + j;
                     J_coo[i] = c;
                     if (DEBUG) { printf("136 - op4_load_S    I,J[%d]=%d,%d\n",
                            i, (int) I_coo[i], (int) J_coo[i]); fflush(stdout); }
                     array[i*NPT] = column[str_data[s].N_idx + j*NPT];
                     if (NPT > 1) {
                        array[i*NPT + 1] = column[str_data[s].N_idx + j*NPT + 1];
                     }
                     i++;
                     // printf("  row %3d  Real=%e\n", S[s].start_row + j, N[S[s].N_idx + j*NPT]);
                }
             }
        } else {       /* dense  */
            for (r = 0; r < nRow*NPT; r++) {
                array[c*nRow*NPT + r] = column[r];
            }
        }
    }

    if (storage) { /* sparse */
        free(str_data);
    }
    free(column);
    if (DEBUG) { printf("157 <- op4_load_S\n"); fflush(stdout); }
    return array;
} /* 1}}} */
double *op4_load_D(FILE   *fp         ,  /* {{{1 */
                   int     filetype   ,  /* in  1=text, other is binary    */
                   int     nRow       ,  /* in  # rows    in matrix        */
                   int     nCol       ,  /* in  # columns in matrix        */
                   char   *fmt_str    ,  /* in  eg "%23le%23le%23le"       */
                   int     col_width  ,  /* in  # characters in format str */
                   int     storage    ,  /* in  0=dense  1,2=sparse  3=ccr */
                   int     complx     ,  /* in  0=real   1=complex         */
                   int     n_Nnz      ,  /* in  number of nonzero terms    */
                   double *I_coo      ,  /* out: sparse row ind            */
                   double *J_coo      )  /* out: sparse col ind            */
{
    int     DEBUG = 0;
    double *array, *column;
    int     r, c, s, i, j, size, NPT, nnz, nType, s_ptr, n_ptr;
    str_t  *str_data;

    int     endian = 0;      /* FIX THIS */

    if (DEBUG) { printf("179 -> op4_load_D:  %d x %d\n", nRow, nCol); fflush(stdout); }
    i     = 0;
    nnz   = 0;
    size  = nRow*nCol;  /* only used for dense */
    NPT   = 1;
    nType = 2;
    if (complx) {
        NPT   *= 2;
        nType += 2;
    }
    if (!storage) { /* dense */
        array    = malloc(sizeof(double)*size*NPT);
        column   = malloc(sizeof(double)*nRow*NPT); /* only need one column at a time */
    } else {        /* sparse; must malloc all numeric data */
        array    = malloc(sizeof(double)*n_Nnz*NPT);
        column   = malloc(sizeof(double)*n_Nnz*NPT);
        str_data = malloc(sizeof(str_t)*(nRow+1)/2);  /* max # strings in a row is nR/2 + 1*/
        if (DEBUG) { printf("op4_load_D malloc'ed array, column, str_data\n"); fflush(stdout); }
    }
    for (c = 0; c < nCol; c++) {
        if (DEBUG) { printf("op4_load_D column %d\n", c); fflush(stdout); }
        if (storage) { /* sparse */
            s_ptr                 = 0;
            n_ptr                 = 0;
            str_data[0].start_row = 0;
            str_data[0].len       = 0;
            str_data[0].N_idx     = 0;
        }
        if (filetype == 1) {
            nnz = op4_read_col_text(fp, c+1, nRow, nCol, fmt_str, col_width,
                                   storage   , /* in 0=dn  1,2=sp1,2  3=ccr  */
                                   complx    , /* in  0=real   1=complex     */
                                  &s_ptr     , /* in/out index m.S (if sp1,2)*/
                                   str_data  , /* out string data (if sp1,2) */
                                  &n_ptr     , /* in/out index m.N (sp 1,2)  */
//                                &array[c*nRow*NPT]
                                   column
                                  ); /* out numeric data */
        } else {
            nnz = op4_read_col_binary(fp        ,
                                   endian    ,  /* in  0=native   1=flipped    */
                                   c+1       ,  /* in  requested column to read   */
                                   nRow      ,  /* in  # rows    in matrix        */
                                   nCol      ,  /* in  # columns in matrix        */
                                   nType     ,  /* in  1=RS 2=RD 3=CS 4=CD        */
                                   storage   ,  /* in  0=dn; 1=sp1; 2=sp2         */
                                  &s_ptr     ,  /* in/out idx str_data[] (s_o)=1,2*/
                                   str_data  ,  /* out string data   (s_o)=1,2    */
                                  &n_ptr     ,  /* in/out idx N[]    (s_o)=1,2    */
//                                &array[c*nRow*NPT]     /* out numeric data               */
                                   column
                                  );
        }
        if (DEBUG) { printf("op4_load_D column %d populating array\n", c); fflush(stdout); }
        if (storage) { /* sparse */
             for (s = 0; s < s_ptr; s++) {
                if (DEBUG) { printf("op4_load_D    str_data[%d].len=%d  .start_row=%d  .N_idx=%d\n",
                        s, str_data[s].len, str_data[s].start_row, str_data[s].N_idx); fflush(stdout); }
                for (j = 0; j < str_data[s].len; j++) {
                     I_coo[i] = str_data[s].start_row + j;
                     J_coo[i] = c;
                     if (DEBUG) { printf("op4_load_D    I,J[%d]=%d,%d\n",
                            i, (int) I_coo[i], (int) J_coo[i]); fflush(stdout); }
                     array[i*NPT] = column[str_data[s].N_idx + j*NPT];
                     if (NPT > 1) {
                        array[i*NPT + 1] = column[str_data[s].N_idx + j*NPT + 1];
                     }
                     i++;
                }
             }
        } else {       /* dense  */
            for (r = 0; r < nRow*NPT; r++) {
                if (DEBUG) { printf("   op4_load_D array[%d] = column[%d] = %e\n", c*nRow*NPT + r, r, column[r]); fflush(stdout); }
                array[c*nRow*NPT + r] = column[r];
            }
        }
    }
    if (storage) { /* sparse */
        free(str_data);
    }
    free(column);
    if (DEBUG) { printf("260 <- op4_load_D\n"); fflush(stdout); }
    return array;
} /* 1}}} */
int  op4_filetype(const char *filename)    /* {{{1 */
{
    /* Returns
           1 if a text file
           2 if a binary file native  endian
           3 if a binary file flipped endian
           0 if none of the above or file read error
     */
	int DEBUG = 0;
    FILE *fp;
    int   word, type, nRead;
    char *w;

    w  = (char *) &word;
    fp = fopen(filename, "rb"); /* open in binary mode */
    if (fp == NULL) {
        return 0;
    }
    nRead = fread(&word, BYTES_PER_WORD, 1, fp);
if (DEBUG)
printf("283 - op4_filetype first word is %d  chars:[%c][%c][%c][%c]\n",
word, w[0], w[1], w[2], w[3]);
    if (nRead == 1) { /* in units of words */
        if (word >= 538976288) {
            /* If this is a text file the smallest possible value for
               the integer comprising the first four bytes is
               538976288 which corresponds to four blank spaces.
             */
            if ((ISBLANK(      w[0]) && ISBLANK(      w[1]) &&   /* \s\s\s\s */
                 ISBLANK(      w[2]) && ISBLANK(      w[3])) ||
                (ISBLANK(      w[0]) && ISBLANK(      w[1]) &&   /* \s\s\s\d */
                 ISBLANK(      w[2]) && isdigit((int) w[3])) ||
                (ISBLANK(      w[0]) && ISBLANK(      w[1]) &&   /* \s\s\d\d */
                 isdigit((int) w[2]) && isdigit((int) w[3])) ||
                (ISBLANK(      w[0]) && isdigit((int) w[1]) &&   /* \s\d\d\d */
                 isdigit((int) w[2]) && isdigit((int) w[3])) ||
                (isdigit((int) w[0]) && isdigit((int) w[1]) &&   /* \d\d\d\d */
                 isdigit((int) w[2]) && isdigit((int) w[3]))) {
                type = 1;
            } else {
                type = 0;
            }

        } else if (word == 24) {        /* Fortran binary record header,     */
            type = 2;                   /* native endian, for 6 word block   */

        } else if (word == 402653184) { /* Fortran binary record header,     */
            type = 3;                   /* opposite endian, for 6 word block */

        } else {    /* the first word indicates the file is not an .op4 */
            type = 0;
        }
    } else {        /* unable to read four bytes from the file */
        type = 0;
    }
    fclose(fp);
    return type;
} /* 1}}} */
int  op4_scan(  const char *filename  ,  /* in  {{{1                    */
                int        *n_mat     ,  /* out number of matrices      */
                char        name[][9] ,  /* out matrix names            */
                int        *storage   ,  /* out 0=dn; 1=sp1; 2=sp2      */
                int        *nRow      ,  /* out number of rows          */
                int        *nCol      ,  /* out number of columns       */
                int        *nStr      ,  /* out number of strings       */
                int        *nNnz      ,  /* out number of nonzero terms */
                int        *type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD  */
                int        *form      ,  /* out matrix form 6=symm, etc */
                int        *digits    ,  /* out size of mantissa        */
                long       *offset       /* out byte offset to matrix   */
                )
{
    /* Count the number and determine storage type of matrices
       in the text op4 file by scanning for matrix, column & string headers.
       Driver to text and binary versions.
     */

    int i, result /* = 0  causes segfault */, filetype, endian;

    filetype = op4_filetype(filename); /* 1 if a text file
                                          2 if a binary file native  endian
                                          3 if a binary file flipped endian
                                          0 if none of the above or file error
                                        */
    switch ( filetype ) {
        case 1  :
            result =
            op4_scan_text(filename  ,  /* in                          */
                       n_mat     ,  /* out number of matrices      */
                       name      ,  /* out matrix names            */
                       storage   ,  /* out 0=dn; 1=sp1; 2=sp2      */
                       nRow      ,  /* out number of rows          */
                       nCol      ,  /* out number of columns       */
                       nStr      ,  /* out number of strings       */
                       nNnz      ,  /* out number of nonzero terms */
                       type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD  */
                       form      ,  /* out matrix form 6=symm, etc */
                       digits    ,  /* out size of mantissa        */
                       offset);     /* out byte offset to matrix   */
            break;
        case 2:
        case 3:
            endian = filetype - 2;
            result =
            op4_scan_binary(filename  ,  /* in                          */
                       endian    ,  /* in  0=native   1=flipped    */
                       n_mat     ,  /* out number of matrices      */
                       name      ,  /* out matrix names            */
                       storage   ,  /* out 0=dn; 1=sp1; 2=sp2      */
                       nRow      ,  /* out number of rows          */
                       nCol      ,  /* out number of columns       */
                       nStr      ,  /* out number of strings       */
                       nNnz      ,  /* out number of nonzero terms */
                       type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD  */
                       form      ,  /* out matrix form 6=symm, etc */
                       offset);     /* out byte offset to matrix   */
            if (endian)
                for (i = 0; i < *n_mat; i++)
                    digits[i] = -1;
            else
                for (i = 0; i < *n_mat; i++)
                    digits[i] =  0;
            break;
        case 0  :
            result = 0;
            break;
          default :
            result = 0;
            break;
    }
    return result;
} /* 1}}} */
int  op4_scan_text(const char *filename  ,  /* in  {{{1                    */
                int        *n_mat     ,  /* out number of matrices      */
                char        name[][9] ,  /* out matrix names            */
                int        *storage   ,  /* out 0=dn; 1=sp1; 2=sp2      */
                int        *nRow      ,  /* out number of rows          */
                int        *nCol      ,  /* out number of columns       */
                int        *nStr      ,  /* out number of strings       */
                int        *nNnz      ,  /* out number of nonzero terms */
                int        *type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD  */
                int        *form      ,  /* out matrix form 6=symm, etc */
                int        *digits    ,  /* out size of mantissa        */
                long       *offset       /* out byte offset to matrix   */
               )
{
    /* Count the number and determine storage type of matrices
       in the text op4 file by scanning for matrix, column & string headers.
       Text version.
       Keep track of byte location, matrix name and matrix dims.
       Here are some typical headers:
      10      10       6       2EYE     1P,3E23.16
      10     -10       6       2EYE     1P,3E23.16
123456789x123456789x123456789x123456789x123456789x123456789x
                                 |-> optional from here on

       The first is for dense and sparse form 1, and the
       second is for sparse form 2.
       The headers are distinguished from numeric data in that
       there is no decimal point in the third column--a quality
       all text based numeric data has, as in
-3.937729118E+02-4.689545312E+07 3.766712501E+07-6.581960491E+06 0.000000000E+00
     */
	int DEBUG = 0;
    FILE *fp;
    char  line[OP4_TXT_LINE_SIZE + 1];
    int   n_text_cols, text_col_width, mixed, nwords, start_row, nRead;
    long  location;
    char  junk1[12], junk2[12], junk3[12];

if (DEBUG)
printf("434 - op4_scan_text file=[%s]\n", filename);
    *n_mat = 0;
    fp    = fopen(filename, "r");
    if (fp == NULL) {
        return 0;
    }
if (DEBUG)
printf("441 - op4_scan_text opened ok\n");
    while (!feof(fp)) {
        location = ftell(fp);
        fgets(line, OP4_TXT_LINE_SIZE, fp);
        if (strlen(line) >= 8) {
            /* actuall the op4 file is illegal if line length < 8 */
            switch (op4_line_type(line)) {
                case OP4_TEXT_NUMERIC   : /* {{{2 */
                    break; /* 2}}} */
                case OP4_TEXT_HEADER    : /* {{{2 */
                    if (*n_mat && !storage[*n_mat-1]) {
                        /* previous matrix stored in dense scheme */
                        nNnz[*n_mat-1] = nRow[*n_mat-1]*nCol[*n_mat-1];
                    }
                    storage[*n_mat] = 0;  /* dense unless proven otherwise */
                    nStr[*n_mat]    = 0;
                    nNnz[*n_mat]    = 0;
                    if (line[41] == 'P') {
                        nRead = sscanf(line,
                                      "%8d%8d%8d%8d%8s%[ 1P,]%1d%1[E]%d%1[.]%d",
                                      &nCol[*n_mat],
                                      &nRow[*n_mat],
                                      &form[*n_mat],
                                      &type[*n_mat],
                                       name[*n_mat],
                                       junk1,           /*   1P,  */
                                      &n_text_cols,
                                       junk2,           /*   E    */
                                      &text_col_width,
                                       junk3,           /*   .    */
                                      &digits[*n_mat]);
                        if (nRead != 11) {
                            fclose(fp);
                            return 0;
                        }
                    } else {
                        n_text_cols     =  5;
                        text_col_width  = 15;
                        digits[*n_mat]  =  9;
                        nRead = sscanf(line, "%8d%8d%8d%8d%8s",
                                      &nCol[*n_mat],
                                      &nRow[*n_mat],
                                      &form[*n_mat],
                                      &type[*n_mat],
                                       name[*n_mat]);
                        if (nRead != 5) {
                            fclose(fp);
                            return 0;
                        }
                    }
                    nRow[*n_mat]   = fabs(nRow[*n_mat]);
                    offset[*n_mat] = location;
                    ++(*n_mat);
                    break; /* 2}}} */
                case OP4_TEXT_NEW_COLUMN: /* {{{2 */
if (DEBUG)
printf("497 - new col: %s", line);
                    break; /* 2}}} */
                case OP4_TEXT_STRING_2  : /* {{{2 */
                    sscanf(line, "%8d%8d", &nwords, &start_row);
                    --nwords;
                    nNnz[   *n_mat-1] += nwords/op4_words_per_term[type[*n_mat-1]];
if (DEBUG)
printf("504 - s2:  nwords=%d  nNnz=%d\n", nwords, nNnz[   *n_mat-1]);
                    storage[*n_mat-1] = 2;
                    ++nStr[ *n_mat-1];
                    break; /* 2}}} */
                case OP4_TEXT_STRING_1  : /* {{{2 */
if (DEBUG)
printf("510 - s1:  %s", line);
                    sscanf(line, "%8d", &mixed);
                    nwords    = (mixed/65536) - 1;
                    start_row = mixed - 65536 * (nwords + 1);
                    nNnz[   *n_mat-1] += nwords/op4_words_per_term[type[*n_mat-1]];
                    storage[*n_mat-1] = 1;
                    ++nStr[*n_mat-1];
                    break; /* 2}}} */
                case OP4_TEXT_ERROR     : /* {{{2 */
                    fclose(fp);
                    return 0;
                    break; /* 2}}} */
            }
        }
    }
    offset[*n_mat] = ftell(fp);  /* last byte in the file */
    fclose(fp);
    if (!storage[*n_mat-1]) {
        /* last matrix stored in dense scheme */
        nNnz[*n_mat-1] = nRow[*n_mat-1]*nCol[*n_mat-1];
    }
    return 1;
} /* 1}}} */
int  op4_scan_binary(const char *filename  ,  /* in  {{{1                    */
                int         endian    ,  /* in  0=native   1=flipped    */
                int        *n_mat     ,  /* out number of matrices      */
                char        name[][9] ,  /* out matrix names            */
                int        *storage   ,  /* out 0=dn; 1=sp1; 2=sp2      */
                int        *nRow      ,  /* out number of rows          */
                int        *nCol      ,  /* out number of columns       */
                int        *nStr      ,  /* out number of strings       */
                int        *nNnz      ,  /* out number of nonzero terms */
                int        *type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD  */
                int        *form      ,  /* out matrix form 6=symm, etc */
                long       *offset       /* out byte offset to matrix   */
               )
{
    /* Count the number and determine storage type of matrices
       in the text op4 file by scanning for matrix, column & string headers.
       Binary version.
     */
	int DEBUG = 0;
    FILE *fp;
    int   record_length, is_header, result /* = 0  causes segfault */;
    char  my_name[9], next_byte;
    long  location;

if (DEBUG)
printf("558 - op4_scan_binary file=[%s]\n", filename);
    *n_mat = 0;
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        return 0;
    }
if (DEBUG)
printf("565 - op4_scan_binary opened ok\n");
    while (!feof(fp)) {
        location  = ftell(fp);
if (DEBUG) printf("568 - op4_scan_binary location A =%ld\n", location);
        /* Determining EOF is tricky because positioning the file pointer
           via Fortran record lengths ends up putting the FP at the very
           last byte in the file... without triggering EOF.  Have to
           explicitly look ahead one byte to get feof() to register .
         */
        next_byte = fgetc(fp);
        if (feof(fp)) /* ok I'm really done with this file */
            break;
        else          /* oops, am not at EOF, put the byte back */
            ungetc(next_byte, fp);

        op4_is_mat_header_binary(fp              ,
                            endian          ,
                           &record_length   ,
                           &is_header       ,
                            my_name         ,
                           &storage[*n_mat] ,
                           &nRow[*n_mat]    ,
                           &nCol[*n_mat]    ,
                           &nStr[*n_mat]    ,
                           &nNnz[*n_mat]    ,
                           &type[*n_mat]    ,
                           &form[*n_mat]    ,
                           &offset[*n_mat]);
if (DEBUG) printf("593 - op4_scan_binary got RL=%d\n", record_length);
        fseek(fp, record_length + 2*BYTES_PER_WORD, SEEK_CUR); /* skip header */
if (DEBUG) printf("595 - op4_scan_binary loc=%ld is_header=%d\n", ftell(fp), is_header);
        if (is_header) {
            strncpy(name[*n_mat], my_name, 9);
            if (storage[*n_mat]) { /* sparse 1 or 2; need to count strings */
                result = op4_count_str_binary(fp, endian, type[*n_mat],
                                         storage[*n_mat],
                                         nRow[*n_mat], nCol[*n_mat],
                                        &nStr[*n_mat], &nNnz[*n_mat]);
                if (!result) {
                    fclose(fp);
                    *n_mat = 0;
                    return result;
                }
            }
            ++(*n_mat);
        }
if (DEBUG) {
location  = ftell(fp);
printf("613 - op4_scan_binary location B =%ld nStr=%d\n", location, nStr[*n_mat-1]);
}
    }

    fseek(fp, 0, SEEK_END);      /* Set the file pointer to the end of file. */
    offset[*n_mat] = ftell(fp);  /* If this isn't done ftell() returns bogus.*/
if (DEBUG)
printf("620 - op4_scan_binary end byte position=%ld\n", offset[*n_mat]);
    fclose(fp);
    if (!storage[*n_mat-1]) {
        /* last matrix stored in dense scheme */
        nNnz[*n_mat-1] = nRow[*n_mat-1]*nCol[*n_mat-1];
    }

    return 1;
} /* 1}}} */
void op4_is_mat_header_binary(FILE *fp,        /* {{{1 */
                         int   endian    ,     /* in  0=native   1=flipped  */
                         int  *record_length,  /* out */
                         int  *is_header ,  /* out 1=is matrix hdr; 0=isn't */
                         char  my_name[] ,  /* out matrix name              */
                         int  *storage   ,  /* out 0=dn; 1=sp1; 2=sp2       */
                         int  *nRow      ,  /* out number of rows           */
                         int  *nCol      ,  /* out number of columns        */
                         int  *nStr      ,  /* out number of strings        */
                         int  *nNnz      ,  /* out number of nonzero terms  */
                         int  *Type      ,  /* out 1=RS; 2=RD; 3=CS; 4=CD   */
                         int  *Form      ,  /* out matrix form 6=symm, etc  */
                         long *offset)      /* out byte offset to matrix    */
/*
   If this Fortran record contains a matrix header, reads the appropriate
   entries.  The file pointer is rewound to where it was when
   this routine was called.
 */
{
	int DEBUG = 0;
    int   is_ascii, nCol_, nRow_, Form_, Type_, rec_len,
          column_id, start_row, n_words, header[7];

    long  location;
    char *name_ptr;

    *is_header = 0;

    location   = ftell(fp);
    fread(record_length, BYTES_PER_WORD, 1, fp);
	if (DEBUG)
		printf("660 - op4_is_mat_header_binary BYTES_PER_WORD=%d endian=%s RL=%d location=%ld\n", BYTES_PER_WORD,endian,*record_length, location);
    if (endian)
        *record_length = flip_bytes_int(*record_length);
	if (DEBUG)
		printf("664 - op4_is_mat_header_binary RL=%d location=%ld\n", *record_length, location);

    if (*record_length == 24) {
        /* candidate for a matrix header */
        fread(header, BYTES_PER_WORD, 7, fp); /* + trailing rec marker */
        /* words header[4],header[5] make up the ASCII name--is it valid? */
        name_ptr = (char *) &header[4];
        is_ascii = 1;
        if (op4_valid_name(name_ptr)) {
            strncpy(my_name, name_ptr,  9);
        } else {
            strncpy(my_name, "unnamed", 9);
        }

        nCol_  = header[0];
        nRow_  = header[1];
        Form_  = header[2];
        Type_  = header[3];
        if (endian) {
            nCol_ = flip_bytes_int(nCol_);
            nRow_ = flip_bytes_int(nRow_);
            Form_ = flip_bytes_int(Form_);
            Type_ = flip_bytes_int(Type_);
        }
		if (DEBUG) {
			printf("689 - op4_scan_binary name=[%s] ", my_name);
			printf("690 is_ascii=%d nCols=%d nRows=%d Form=%d Type=%d\n",
			is_ascii,nCol_,nRow_,Form_,Type_);
		}
        if (is_ascii && nCol_ > 0 && nRow_ && Form_ < 50 && Type_ < 10) {
            /* we have a winner; this is a matrix header */
            *is_header = 1;
            *nRow      = fabs(nRow_);
            *nCol      = nCol_;
            *offset    = location;
            *Type      = Type_;
            *Form      = Form_;
            *nStr      = 0;
            *nNnz      = 0;

            fread(&rec_len  , BYTES_PER_WORD, 1, fp);
            fread(&column_id, BYTES_PER_WORD, 1, fp);
            fread(&start_row, BYTES_PER_WORD, 1, fp);
            fread(&n_words  , BYTES_PER_WORD, 1, fp);
            if (endian) {
                rec_len   = flip_bytes_int(rec_len);
                column_id = flip_bytes_int(column_id);
                start_row = flip_bytes_int(start_row);
                n_words   = flip_bytes_int(n_words);
            }
if (DEBUG)
printf("715 - op4_scan_binary rl=%2d c=%2d r=%2d nw=%2d\n",
rec_len, column_id, start_row, n_words);
            if (nRow_ < 0) {
                *storage = 2;   /* sparse 2 */
            } else if (!start_row ||
                       ((start_row == 1) && (n_words == 1))) { /* null */
                *storage = 1;   /* sparse 1 */
            } else {                   /* dense    */
                *storage = 0;
            }
        }
    }
    fseek(fp, location, SEEK_SET);
} /* 1}}} */
int  op4_count_str_binary(FILE *fp     ,  /* in  {{{1 */
                     int   endian ,  /* in  */
                     int   nType  ,  /* in  */
                     int   storage,  /* in  */
                     int   nRows  ,  /* in  */
                     int   nCols  ,  /* in  */
                     int  *nStr   ,  /* out */
                     int  *nNnz)     /* out */
{
	int DEBUG = 0;
    int col, row, n_words, length_in_words, record_length, encoded,
        word_count, end_of_data;

#define ErrStrSize 200
    char error_msg[ErrStrSize - 1];
    *nStr = 0;
    col   = 1;

    while (col <= nCols) {

		if (DEBUG)
		printf("750 - op4_count_str_binary top of loop loc=%ld\n", ftell(fp));
        fread(&record_length, BYTES_PER_WORD, 1, fp);
        fread(&col,           BYTES_PER_WORD, 1, fp);
        fread(&end_of_data,   BYTES_PER_WORD, 1, fp);
        fread(&n_words,       BYTES_PER_WORD, 1, fp);
        if (endian) {
            record_length = flip_bytes_int(record_length);
            col           = flip_bytes_int(col);
            end_of_data   = flip_bytes_int(end_of_data);
            n_words       = flip_bytes_int(n_words);
        }
		if (DEBUG) {
		printf("762 - op4_count_str_binary record_length  =%d\n", record_length);
		printf("763 - op4_count_str_binary col            =%d\n", col          );
		printf("764 - op4_count_str_binary end_of_data    =%d\n", end_of_data  );
		printf("765 - op4_count_str_binary n_words        =%d\n", n_words      );
		}
        /* make sure values aren't bogus {{{2 */
        if ((end_of_data != 0) && (end_of_data != 1)) {
            *nNnz = 0;
            snprintf(error_msg, ErrStrSize,
                    "end_of_data=%d out of range, byte offset %ld\n",
                     end_of_data, ftell(fp));

            /* propagate(error_msg); */
            return 0;
        }
        if ((col < 1)          || (col > (nCols + 2))) {
            *nNnz = 0;
            snprintf(error_msg, ErrStrSize,
                    "col=%d out of range, byte offset %ld\n",
                     col, ftell(fp));

            /* propagate(error_msg); */
            return 0;
        }
        if ((record_length < 4*BYTES_PER_WORD) ||
            (record_length > BYTES_PER_WORD *
                                (5 + op4_words_per_term[nType]*nRows))) {
                              /* 5 = 3 column header ints +
                                     2 string header ints */
            *nNnz = 0;
            snprintf(error_msg, ErrStrSize,
                    "RL=%d out of range (should be < %d), byte offset %ld\n",
                     record_length,
                     BYTES_PER_WORD*(5 + op4_words_per_term[nType]*nRows),
                     ftell(fp));

            /* propagate(error_msg); */
            return 0;
        }
        /* 2}}} */
		if (DEBUG)
			printf("803 - op4_count_str_binary A RL=%d col=%d enddata=%d n_words=%d  loc=%ld\n",
			record_length, col, end_of_data, n_words, ftell(fp));

        if (end_of_data) {
            fseek(fp, record_length - 2*BYTES_PER_WORD, SEEK_CUR);
            break;
        }
        word_count = 0;
        while (word_count < n_words) {
            if (storage == 1) {
                fread(&encoded, BYTES_PER_WORD, 1, fp);
                word_count += 1;
                if (endian) {
                    encoded = flip_bytes_int(encoded);
                }
                length_in_words = (encoded/65536) -1;
				if (DEBUG)
					printf("820 - op4_count_str_binary B sp1 length_words=%d encoded=%d loc=%ld\n",
					length_in_words, encoded, ftell(fp));
            } else {  /* storage == 2 */
                fread(&length_in_words, BYTES_PER_WORD, 1, fp);
                fread(&row,             BYTES_PER_WORD, 1, fp); /* unused */
                word_count += 2;
                if (endian) {
                    length_in_words = flip_bytes_int(length_in_words);
                }
                --length_in_words;
				if (DEBUG)
					printf("831 - op4_count_str_binary C sp2 length_words=%d row=%d  loc=%ld\n",
					length_in_words, row, ftell(fp));
            }
            if (length_in_words < 0) {
                snprintf(error_msg, ErrStrSize,
                        "length_in_words=%d out of range, byte offset %ld\n",
                         length_in_words, ftell(fp));

                /* propagate(error_msg); */

            }
            *nNnz += length_in_words;
            ++(*nStr);

            /* skip numeric data */
            fseek(fp, length_in_words*BYTES_PER_WORD, SEEK_CUR);
            word_count += length_in_words;
			if (DEBUG)
				printf("849 - op4_count_str_binary D word_count=%d  loc=%ld\n", word_count, ftell(fp));
        }
        fread(&record_length, BYTES_PER_WORD, 1, fp);
        if (endian) {
            record_length = flip_bytes_int(record_length);
        }
		if (DEBUG)
			printf("856 - op4_count_str_binary E end RL=%d  loc=%ld\n", record_length, ftell(fp));
    }
    *nNnz /= op4_words_per_term[nType];
    return 1;
}
/* 1}}} */
int  op4_line_type(const char *line) {  /* in {{{1 */
    int length, type;

	int DEBUG = 0;
    length = strlen(line);
    if (length > 2 && line[2] == '.') {  /* most common case:  numeric data */
        type = OP4_TEXT_NUMERIC;
    } else if (length >= 34) {  /* header string must be at least 34 chars  */
        type = OP4_TEXT_HEADER;
    } else if (length >= 24 && line[23] >= '0' && line[23] <= '9') {
        type = OP4_TEXT_NEW_COLUMN;
    } else if (length >= 16 && line[15] >= '0' && line[15] <= '9') {
        type = OP4_TEXT_STRING_2;
    } else if (length >=  8 && line[ 7] >= '0' && line[ 7] <= '9') {
        type = OP4_TEXT_STRING_1;
    } else {
if (DEBUG)
printf("879 - op4_line_type Error on [%s]\n", line);
        type = OP4_TEXT_ERROR;
    }
if (DEBUG)
printf("883 - T=%d:[%s]\n", type, line);
    return type;
} /* 1}}} */
int  op4_read_col_text(FILE   *fp         ,  /* {{{1 */
                    int     c_in       ,  /* in  requested column to read   */
                    int     nRow       ,  /* in  # rows    in matrix        */
                    int     nCol       ,  /* in  # columns in matrix        */
                    char   *fmt_str    ,  /* in  eg "%23le%23le%23le"       */
                    int     col_width  ,  /* in  # characters in format str */
                    int     storage    ,  /* in  0=dense  1,2=sparse        */
                    int     complx     ,  /* in  0=real   1=complex         */
                    int    *n_str      ,  /* in/out idx str_data[] (s_o)=1,2*/
                    str_t  *S          ,  /* out string data   (s_o)=1,2    */
                    int    *N_index    ,  /* in/out idx N[]    (s_o)=1,2    */
                    double *N             /* out numeric data               */
                   )
{
	int DEBUG = 0;
    double      x[OP4_MAX_TEXT_COLS];
    int         nRead, j, nNum_cols, still_reading, location, col, row,
                nNum, got_this_col, n_nnz, mixed, nwords, start_row,
                line_type, NPT,
                n_read_this_str = 0;
    char  line[OP4_TXT_LINE_SIZE + 1];

    still_reading = 1;
    got_this_col  = 0;
    n_nnz         = 0;  /* count of numeric values loaded */
    NPT           = Nm_per_term[complx]; /* numbers per term; real=1 complex=2*/
    if (!storage) {
        /* zero out dense column */
        for (j = 0; j < NPT*nRow; j++)
            N[j] = 0.0;
    }
    while (still_reading) {
        location = ftell(fp);
        fgets(line, OP4_TXT_LINE_SIZE, fp);
        line_type = op4_line_type(line);
if (DEBUG) {
printf("922 - op4_read_col_text line_type=%d\n", line_type); fflush(stdout); }
        if (strlen(line) >= 8) {
            switch (line_type) {
                case OP4_TEXT_NUMERIC   : /* {{{2 */
                    nNum_cols = (int) (strlen(line) / col_width);
                    nRead = sscanf(line, fmt_str, &x[0], &x[1], &x[2],
                                                  &x[3], &x[4], &x[5],
                                                  &x[6], &x[7], &x[8]);
                    if (storage) { /* sparse  */
if (DEBUG) {
printf("932 - op4_read_col_text inserting %d terms at row=%d\n",
nRead, S[*n_str-1].start_row); fflush(stdout); }
                        n_read_this_str += nRead;
                        S[*n_str-1].len = n_read_this_str/NPT;
                        for (j = 0; j < nRead; j++) {
                            N[(*N_index)++] = x[j]; }
                    } else {       /* dense  */
if (DEBUG) { printf("939 - op4_read_col_text inserting %d terms at row=%d\n", nRead, row); fflush(stdout); }
                        for (j = 0; j < nRead; j++) {
if (DEBUG) {
printf("942 - N[row=%d + *N_index=%d]=%e\n", row, (*N_index), x[j]); fflush(stdout); }
                            N[row + j] = x[j];
                        }
                        row += nRead;
                    }
                    n_nnz += nRead;
if (DEBUG) {
printf("949 - op4_read_col_text nR=%d n_nnz now=%d\n", nRead, n_nnz); fflush(stdout); }
                    break; /* 2}}} */
                case OP4_TEXT_HEADER    : /* {{{2 */


                /* propagate("op4_read_col_text new matrix header while reading column"); */

                    break; /* 2}}} */
                case OP4_TEXT_NEW_COLUMN: /* {{{2 */
if (DEBUG)
printf("959 - op4_read_col_text start new column (got_this_col=%d) loc=%ld\n",
got_this_col, ftell(fp));
                    if (got_this_col) {
                        /* have read next column's header; back up & exit */
                        fseek(fp, location, SEEK_SET);
                        still_reading = 0;
                    } else {
                        nRead = sscanf(line, "%8d%8d%8d", &col, &row, &nNum);
if (DEBUG)
printf("968 - requested col=%d  this col=%d  this row=%d\n", c_in, col, row);
                        if (c_in != col) {
                            /* this is not the column of interest            */
                            fseek(fp, location, SEEK_SET);
                            return n_nnz;
                        }
                        --col; --row; /* internal indexing is 0 based */
                        if (!storage) /* dense  */
                            row *= NPT;
if (DEBUG)
printf("\n978 - new column %2d row %2d complx=%d\n", col, row, complx);
                        got_this_col = 1;
                        if ((col > nCol) || (nRead != 3)) {
                            /* have reached the spurious end of matrix entry */
                            /* or have encountered bad data                  */
                            return n_nnz;
                        }
                    }
                    break; /* 2}}} */
                case OP4_TEXT_STRING_1  :
                case OP4_TEXT_STRING_2  : /* {{{2 */
                    n_read_this_str = 0;
                    if (line_type == OP4_TEXT_STRING_1) {
                        sscanf(line, "%8d", &mixed);
                        nwords    = (mixed/65536) - 1;
                        start_row = mixed - 65536 * (nwords + 1);
                    } else { /* type 2 */
                        sscanf(line, "%8d%8d", &nwords, &start_row);
                        --nwords;
                    }
                    S[*n_str].start_row = start_row - 1;
                    S[*n_str].len       = 0;
                    S[*n_str].N_idx     = *N_index;
                    ++(*n_str);
                    break; /* 2}}} */
                case OP4_TEXT_ERROR     : /* {{{2 */


                    /* propagate("op4_read_col_text error reading a column "); */

                    break; /* 2}}} */
            }
        }
    }
    return n_nnz;
} /* 1}}} */
int  op4_read_col_binary(FILE   *fp         ,  /* {{{1 */
                    int     endian     ,  /* in  0=native   1=flipped    */
                    int     c_in       ,  /* in  requested column to read   */
                    int     nRow       ,  /* in  # rows    in matrix        */
                    int     nCol       ,  /* in  # columns in matrix        */
                    int     nType      ,  /* in  1=RS 2=RD 3=CS 4=CD        */
                    int     storage    ,  /* in  0=dn; 1=sp1; 2=sp2         */
                    int    *n_str      ,  /* in/out idx str_data[] (s_o)=1,2*/
                    str_t  *S          ,  /* out string data   (s_o)=1,2    */
                    int    *N_index    ,  /* in/out idx N[]    (s_o)=1,2    */
                    double *N             /* out numeric data               */
                   )
{
	int DEBUG = 0;
    int    complx, NPT, WPN, BPN, i, record_length, col, start_row,
           n_numbers, got_this_col, n_nnz, is_header, unused, mixed,
           n_w_this_col, n_w_this_str;
    long   l_unused;
    char   unused_name[9];
    float  xs;

	if (DEBUG)
		printf("\n1036 - op4_read_col_binary 1 location=%ld  endian=%d want col %d\n",
		ftell(fp), endian, c_in);

    complx        = nType > 2 ? 1 : 0;
    NPT           = Nm_per_term[complx]; /* numbers per term; real=1 complex=2*/
    WPN           = nType % 2 ? 1 : 2;   /* words per number; float=1 double=2*/
    BPN           = BYTES_PER_WORD * WPN;/* bytes per number */
    got_this_col  = 0;
    n_nnz         = 0;  /* count of numeric values loaded */
    op4_is_mat_header_binary(fp              , /* in  */
                        endian          , /* in  */
                       &record_length   , /* out */
                       &is_header       , /* out */
                        unused_name     , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &unused          , /* out */
                       &l_unused);        /* out */
	if (DEBUG)
		printf("1059 - op4_read_col_binary 2 complex=%d NPT=%d WPN=%d BPN=%d RL=%d location=%ld\n", complx,NPT,WPN,BPN,record_length, ftell(fp));
    if (is_header) /* yes?  skip it */
        fseek(fp, record_length+2*BYTES_PER_WORD, SEEK_CUR);
	if (DEBUG)
		printf("1063 - op4_read_col_binary 3 RL=%d location=%ld\n", record_length, ftell(fp));

    if (!storage) {
        /* zero out dense column */
        for (i = 0; i < NPT*nRow; i++)
            N[i] = 0.0;
    }

    fread(&record_length, BYTES_PER_WORD, 1, fp);
    fread(&col,           BYTES_PER_WORD, 1, fp);
    fread(&start_row,     BYTES_PER_WORD, 1, fp);
    fread(&n_w_this_col,  BYTES_PER_WORD, 1, fp);
    if (endian) {
        record_length = flip_bytes_int(record_length);
        col           = flip_bytes_int(col);
        start_row     = flip_bytes_int(start_row);
        n_w_this_col  = flip_bytes_int(n_w_this_col);
    }
	if (DEBUG)
			printf("1082 - op4_read_col_binary: requested col=%d  RL=%d got col=%d row=%d nW=%d\n",
			c_in, record_length, col, start_row, n_w_this_col);
    n_numbers = n_w_this_col/WPN;
    if        (c_in > col) { /* requested column ahead of me; skip record */
        fseek(fp, record_length-2*BYTES_PER_WORD, SEEK_CUR);
		if (DEBUG) {
			printf("1088 - op4_read_col_binary: column mismatch, advance to location %ld & exit\n",
			ftell(fp));
		}
        return 0;
    } else if (c_in < col) { /* requested column behind me; back up */
        fseek(fp,              -4*BYTES_PER_WORD, SEEK_CUR);
		if (DEBUG) {
			printf("1095 - op4_read_col_binary: column mismatch, back up to location %ld & exit\n",ftell(fp));
		}
        return 0;
    }

if (DEBUG)
printf("1102 - op4_read_col_binary: going to read col=%d\n", col);
    if (!storage) {             /* dense; read entire column */
        if (WPN == 2) {
            /* double precision:  can do bulk read */
if (DEBUG) printf("1106 - op4_read_col_binary: double prec\n");
            fread(&N[NPT*(start_row-1)], BPN, n_numbers, fp);
            if (endian) { /* opposite endian, need to flip terms */
                for (i = 0; i < n_numbers; i++)
                    N[NPT*(start_row-1) + i] =
                        flip_bytes_double(N[NPT*(start_row-1) + i]);
            }
        } else if (!endian) {
            /* single precision native endian */
if (DEBUG) printf("1115 - op4_read_col_binary: single prec native endian\n");
            for (i = 0; i < n_numbers; i++) {
                fread(&xs, BPN, 1, fp);
                N[NPT*(start_row-1) + i] = xs;
            }
        } else {
            /* single precision opposite endian */
if (DEBUG) printf("1122 - op4_read_col_binary: single prec opposite endian\n");
            for (i = 0; i < n_numbers; i++) {
                fread(&xs, BPN, 1, fp);
                N[NPT*(start_row-1) + i] = flip_bytes_float(xs);
            }
        }
if (DEBUG) {
for (i = 0; i < n_numbers; i++)
printf("1130 - op4_read_col_binary: N[%d]=%le\n",
NPT*(start_row-1) + i,N[NPT*(start_row-1) + i]);
}
        fseek(fp, BYTES_PER_WORD, SEEK_CUR);  /* skip trailing RL */
if (DEBUG)
printf("1135 - op4_read_col_binary: finished numeric data location=%ld\n", ftell(fp));

    } else {  /* sparse; loop over all strings */
        record_length -= 2*BYTES_PER_WORD; /* don't count start, end rec mark */
        while (record_length >= 8) {
            if (storage == 1) {  /* sparse 1 */
                fread(&mixed, BYTES_PER_WORD, 1, fp);
                if (endian) {
                    mixed = flip_bytes_int(mixed);
                }
                n_w_this_str  = (mixed/65536) - 1;
                start_row     = mixed - 65536 * (n_w_this_str + 1);
                record_length -= 1*BYTES_PER_WORD;
            } else {             /* sparse 2 */
                fread(&n_w_this_str, BYTES_PER_WORD, 1, fp);
                fread(&start_row,    BYTES_PER_WORD, 1, fp);
                record_length -= 2*BYTES_PER_WORD;
                if (endian) {
                    n_w_this_str = flip_bytes_int(n_w_this_str);
                    start_row    = flip_bytes_int(start_row);
                }
                --n_w_this_str;
            }
            n_numbers = n_w_this_str/WPN;
if (DEBUG) printf("\n1159 - op4_read_col_binary top of str: start_row=%d n_w_this_str=%d n_numbers=%d N_index=%d\n",
start_row, n_w_this_str, n_numbers, *N_index);

            S[*n_str].start_row = start_row - 1;
            S[*n_str].len       = n_numbers/(complx + 1);
            S[*n_str].N_idx     = *N_index;
            ++(*n_str);
            if (WPN == 2) {
                /* double precision:  can do bulk read */
if (DEBUG) printf("1168 - op4_read_col_binary sp: double prec\n");
                fread(&N[*N_index], BPN, n_numbers, fp);
                if (endian) { /* opposite endian, need to flip terms */
                    for (i = 0; i < n_numbers; i++)
                        N[*N_index + i] = flip_bytes_double(N[*N_index + i]);
                }
            } else if (!endian) {
                /* single precision native endian */
if (DEBUG) printf("1176 - op4_read_col_binary sp: single prec native endian\n");
                for (i = 0; i < n_numbers; i++) {
                    fread(&xs, BPN, 1, fp);
                    N[*N_index + i] = xs;
                }
            } else {
                /* single precision opposite endian */
if (DEBUG) printf("1183 - op4_read_col_binary sp: single prec opposite endian\n");
                for (i = 0; i < n_numbers; i++) {
                    fread(&xs, BPN, 1, fp);
                    N[*N_index + i] = flip_bytes_float(xs);
                }
            }
            record_length -= BYTES_PER_WORD*n_w_this_str;
if (DEBUG) {
printf("1191 - op4_read_col_binary sp end string  RL=%d loc=%ld\n",
record_length, ftell(fp));
for (i = 0; i < n_numbers; i++) {
printf("1194 - op4_read_col_binary N[%3d] = %e\n", *N_index+i, N[*N_index + i]);
}
}
            *N_index += n_numbers;
            n_nnz    += n_numbers;
        }
        fread(&record_length, BYTES_PER_WORD, 1, fp);
    }
    return n_nnz;
} /* 1}}} */
FILE* op4_open_r(const char *filename, long offset) {  /* in {{{1 */
    FILE *fp;
    /* Open an op4 file for reading, optionally positioning the file
     * handle to a location past the start of the file.
     */

    fp = fopen(filename, "rb");
    if (fp == NULL)
        return NULL;
    if (offset)
        fseek(fp, offset, SEEK_SET);
    return fp;
} /* 1}}} */

int  op4_wrt_header(FILE   *fp         ,  /* {{{1 */
                    int     endian     ,  /* in  0=native   1=flipped    */
                    char   *name       ,  /* in  matrix name                */
                    int     nRow       ,  /* in  # rows    in matrix        */
                    int     nCol       ,  /* in  # columns in matrix        */
                    int     nType      ,  /* in  1=RS 2=RD 3=CS 4=CD        */
                    int     Form       ,  /* in  1=rect; 2=square           */
                    int     sparse     ,  /* in  1=sp2; 0=dn                */
                    int     digits        /* in -1=flipped endian binary    */
                                          /*     0=native  endian binary    */
                                          /*    >2=text; number of DIGITS   */
                   )
{
/*
int DEBUG = 0;
*/
    int width, text_cols, RL = 24, 
        F_RL, F_nCol, F_nRow, F_Form, F_nType;


    if (sparse)
        nRow = -nRow;

    if (2 <= digits && digits <= 9) { /* single precision text */
        if (nType <= 2) nType = 1;
        else            nType = 3;
    } else {                          /* all other cases:  force double prec */
        if (nType <= 2) nType = 2;
        else            nType = 4;
    }

    if (digits > 1) {          /* text */
        width     = digits + 7;
        text_cols = (int) 80/width;
        fprintf(fp, "%8d%8d%8d%8d%-8s1P,%dE%d.%d\n",
                    nCol, nRow, Form, nType, name, text_cols, width, digits);

    } else if (digits == -1) { /* binary, opposite endian */
        F_RL    = flip_bytes_int(RL   );
        F_nCol  = flip_bytes_int(nCol );
        F_nRow  = flip_bytes_int(nRow );
        F_Form  = flip_bytes_int(Form );
        F_nType = flip_bytes_int(nType);

        fwrite(&F_RL,      BYTES_PER_WORD, 1, fp);
        fwrite(&F_nCol,    BYTES_PER_WORD, 1, fp);
        fwrite(&F_nRow,    BYTES_PER_WORD, 1, fp);
        fwrite(&F_Form,    BYTES_PER_WORD, 1, fp);
        fwrite(&F_nType,   BYTES_PER_WORD, 1, fp);
        fwrite( name,    2*BYTES_PER_WORD, 1, fp);
        fwrite(&F_RL,      BYTES_PER_WORD, 1, fp);

    } else {                   /* binary, native endian   */
        fwrite(&RL,        BYTES_PER_WORD, 1, fp);
        fwrite(&nCol,      BYTES_PER_WORD, 1, fp);
        fwrite(&nRow,      BYTES_PER_WORD, 1, fp);
        fwrite(&Form,      BYTES_PER_WORD, 1, fp);
        fwrite(&nType,     BYTES_PER_WORD, 1, fp);
        fwrite( name,    2*BYTES_PER_WORD, 1, fp);
        fwrite(&RL,        BYTES_PER_WORD, 1, fp);
    }
    return 1;
} /* 1}}} */
int  op4_wrt_col_dn(FILE   *fp    , /* {{{1 */
                    int     column, /* first column is 0       */
                    int     nRows , /* number of rows          */
                    int     nCols , /* number of columns       */
                    double *A     , /* array of terms to write */
                    int     complx, /* 1=complex  0=real       */
                    int     digits  /* -1=flipped endian; 0=native; >0=digits */
                    )
{
int DEBUG = 0;
    int    i, r, first_nonzero_row, last_nonzero_row, nNums, width, text_cols,
           RL, A_start, F_RL, F_col, F_fnz, F_nNums;
    double F_x;
    char   format[10];

    ++column; /* go from 0- to 1-based indexing */
    first_nonzero_row = nRows;
    last_nonzero_row  =   -1;
    /* find the first and last rows in this column with nonzero's */
    for (r = 0;       r < nRows; r++) {
        if (complx && ((A[2*r] != 0.0) || (A[2*r+1] != 0.0))) {
            first_nonzero_row = r;
            break;
        } else if (A[r] != 0.0) {
            first_nonzero_row = r;
            break;
        }
    }
    for (r = nRows-1; r >= 0;    r--) {
        if (complx && ((A[2*r] != 0.0) || (A[2*r+1] != 0.0))) {
            last_nonzero_row = r;
            break;
        } else if (A[r] != 0.0) {
            last_nonzero_row = r;
            break;
        }
    }
    if (last_nonzero_row < first_nonzero_row) {
        /* column is all zeros; do nothing unless this is the last column */
        if (column == nCols)
            op4_wrt_trailer(fp, column, digits);
        return 1;
    } else {
        nNums = Nm_per_term[complx]*(last_nonzero_row - first_nonzero_row + 1);
    }
if (DEBUG)
printf("op4_wrt_col_dn nNums=%d 1st nonzero_row=%d last_row=%d Nm_per_term[complx]=%d\n", nNums, first_nonzero_row, last_nonzero_row, Nm_per_term[complx]);

    /* write terms of the column from the first to the last nonzero */
    A_start = first_nonzero_row*Nm_per_term[complx];
    if (digits > 1) {          /* text */
        fprintf(fp, "%8d%8d%8d\n", column, first_nonzero_row+1, nNums);
        sprintf(format, "%%%d.%dE", digits+7, digits);
        width     = digits + 7;
        text_cols = (int) 80/width;
        for (i = 0; i < nNums; i++) {
            fprintf(fp, format, A[i + A_start]);
            if (!((i+1) % text_cols))
                fprintf(fp, "\n");
        }
        if (nNums % text_cols)
            fprintf(fp, "\n");

    } else if (digits == -1) { /* binary, opposite endian */
        RL = BYTES_PER_WORD*(3 + 2*nNums);
        first_nonzero_row++;  /* convert to 1-based indexing */
        nNums *= 2;           /* convert to # of words */
        F_RL    = flip_bytes_int(RL);
        F_col   = flip_bytes_int(column);
        F_fnz   = flip_bytes_int(first_nonzero_row);
        F_nNums = flip_bytes_int(nNums);

        fwrite(&F_RL   ,           BYTES_PER_WORD,     1, fp);
        fwrite(&F_col  ,           BYTES_PER_WORD,     1, fp);
        fwrite(&F_fnz  ,           BYTES_PER_WORD,     1, fp);
        fwrite(&F_nNums,           BYTES_PER_WORD,     1, fp); /* # words     */
        for (i = 0; i < nNums/2; i++) {
            F_x = flip_bytes_double(A[i + A_start]);
            fwrite(&F_x,  2*BYTES_PER_WORD, 1, fp);
        }
        fwrite(&F_RL,              BYTES_PER_WORD,     1, fp);

    } else {                   /* binary, native endian   */
        RL = BYTES_PER_WORD*(3 + 2*nNums);
        first_nonzero_row++;  /* convert to 1-based indexing */
        nNums *= 2;           /* convert to # of words */
if (DEBUG)
printf("op4_wrt_col_dn loc 1 = %ld\n", ftell(fp));
        fwrite(&RL,                BYTES_PER_WORD,     1, fp);
        fwrite(&column,            BYTES_PER_WORD,     1, fp);
        fwrite(&first_nonzero_row, BYTES_PER_WORD,     1, fp);
        fwrite(&nNums,             BYTES_PER_WORD,     1, fp); /* # words     */
if (DEBUG)
printf("op4_wrt_col_dn nNums=%d words starting at A[%d]\n", nNums, A_start);
if (DEBUG)
printf("op4_wrt_col_dn loc 2 = %ld\n", ftell(fp));
        fwrite(&A[A_start],        BYTES_PER_WORD, nNums, fp);
if (DEBUG)
printf("op4_wrt_col_dn loc 3 = %ld\n", ftell(fp));
        fwrite(&RL,                BYTES_PER_WORD,     1, fp);
if (DEBUG)
printf("op4_wrt_col_dn loc 4 = %ld\n", ftell(fp));
    }

    if (column == nCols)
        op4_wrt_trailer(fp, column, digits);
    return 1;
} /* 1}}} */
int  op4_wrt_col_sp(FILE         *fp    , /* {{{1 */
                    int           column, /* first column is 0    */
                    int           A_col , /* column index to A[]; this differs
                                             from 'column' if A[] contains only
                                             part of the entire matrix */
                    int           nCols , /* number of columns    */
                    SparseMatrix  A     , /* entire sparse matrix */
                    int           complx, /* 1=complex  0=real    */
                    int           digits  /* -1=flipped; 0=native; >0=digits */
                    )
{
int DEBUG = 0;
    int    i, n, s, nStr, nNums, nTerms, s_ptr, n_ptr, width, text_cols, 
           RL, n_words_str, n_words_col, type, start_row,
           F_RL, F_col, F_zero, F_nnw, F_n_words_str, F_start_row,
           zero = 0;
    char   format[10];
    double x, F_x;

    nStr  = A.S_start[A_col + 1] - A.S_start[A_col];
    nNums = A.N_start[A_col + 1] - A.N_start[A_col]; /* in this column */
    s_ptr = A.S_start[A_col];
    n_ptr = A.N_start[A_col];

    ++column; /* go from 0- to 1-based indexing */
    if (!nStr) { /* empty column */
        if (column == nCols)
            op4_wrt_trailer(fp, column, digits);
        return 1;
    }

    n     = 0;
    type  = 2;     /* double precision real */
    if (complx) {
        type = 4;  /* double precision complex */
        nTerms = nNums/2;
    } else {
        nTerms = nNums;
    }
    if ((digits > 2) && (digits < 10)) /* demote to single precision */
        --type;

if (DEBUG)
printf("op4_wrt_col_sp a col=%d nStr=%d nNums=%d\n", column, nStr, nNums);
    n_words_col = op4_words_per_term[type]*nTerms + 2*nStr;
    if (digits > 1) {          /* text */
        fprintf(fp, "%8d%8d%8d\n", column, 0, n_words_col);
        sprintf(format, "%%%d.%dE", digits+7, digits);
        width     = digits + 7;
        text_cols = (int) 80/width;
        for (s = 0; s < nStr; s++) {
            nNums = Nm_per_term[complx]*A.S[s_ptr + s].len; /* in this string */
            n_words_str = op4_words_per_term[type]*A.S[s_ptr + s].len + 1;
            fprintf(fp, "%8d%8d\n", n_words_str, A.S[s_ptr + s].start_row + 1);
if (DEBUG)
printf("op4_wrt_col_sp b s=%d\n", s);
            for (i = 0; i < nNums; i++) {
                fprintf(fp, format, A.N[n_ptr + n++]);
                if (!((i+1) % text_cols))
                    fprintf(fp, "\n");
            }
            if (nNums % text_cols)
                fprintf(fp, "\n");
        }

    } else if (digits == -1) { /* binary, opposite endian */

if (DEBUG)
printf("flipping for column %d\n", column);
        RL = BYTES_PER_WORD*(3 + n_words_col);

        F_RL    = flip_bytes_int(RL);
        F_col   = flip_bytes_int(column);
        F_zero  = flip_bytes_int(zero);   /* not really necessary */
        F_nnw   = flip_bytes_int(n_words_col);

        fwrite(&F_RL  , BYTES_PER_WORD, 1, fp);
        fwrite(&F_col , BYTES_PER_WORD, 1, fp);
        fwrite(&F_zero, BYTES_PER_WORD, 1, fp);
        fwrite(&F_nnw , BYTES_PER_WORD, 1, fp);
if (DEBUG)
printf("column header RL=%d col=%d nnw=%d\n", RL, column, n_words_col);
        for (s = 0; s < nStr; s++) {
            n_words_str = op4_words_per_term[type]*A.S[s_ptr + s].len + 1;
            start_row   = A.S[s_ptr + s].start_row + 1;
if (DEBUG)
printf("string s=%d nwords=%d start_row=%d loc=%ld  Fnw=%d  Fsr=%d\n", 
s, n_words_str, start_row, ftell(fp), F_n_words_str, F_start_row);

            F_n_words_str = flip_bytes_int(n_words_str);
            F_start_row   = flip_bytes_int(start_row);

            fwrite(&F_n_words_str,   BYTES_PER_WORD, 1, fp);
            fwrite(&F_start_row, BYTES_PER_WORD, 1, fp);

            nNums = Nm_per_term[complx]*A.S[s_ptr + s].len; /* in this string */
            for (i = 0; i < nNums; i++) {
if (DEBUG)
printf("N[%d]=%le\n", n_ptr + n, A.N[n_ptr + n]);
                x   = A.N[n_ptr + n++];
                F_x = flip_bytes_double(x);
                fwrite(&F_x, 2*BYTES_PER_WORD, 1, fp);
            }
        }
        fwrite(&F_RL  , BYTES_PER_WORD, 1, fp);

    } else {                   /* binary, native endian   */

        RL = BYTES_PER_WORD*(3 + n_words_col);

        fwrite(&RL         ,  BYTES_PER_WORD,     1, fp);
        fwrite(&column     ,  BYTES_PER_WORD,     1, fp);
        fwrite(&zero       ,  BYTES_PER_WORD,     1, fp);
        fwrite(&n_words_col,  BYTES_PER_WORD,     1, fp);
        for (s = 0; s < nStr; s++) {
            nNums = Nm_per_term[complx]*A.S[s_ptr + s].len; /* in this string */
            n_words_str = op4_words_per_term[type]*A.S[s_ptr + s].len + 1;
            start_row = A.S[s_ptr + s].start_row + 1;
            fwrite(&n_words_str   ,   BYTES_PER_WORD,     1, fp);
            fwrite(&start_row     ,   BYTES_PER_WORD,     1, fp);
            fwrite(&A.N[n_ptr + n], 2*BYTES_PER_WORD, nNums, fp);
            n += nNums;
        }
        fwrite(&RL         ,  BYTES_PER_WORD,     1, fp);
    }

    if (column == nCols)
        op4_wrt_trailer(fp, column, digits);
    return 1;
} /* 1}}} */
int  op4_wrt_trailer(FILE *fp    , /* {{{1 */
                     int   column, /* first column is 0    */
                     int   digits  /* -1=flipped; 0=native; >0=digits */
                    )
{
int DEBUG = 0;
    int    RL  = 20,  /* 16 if single precision, 20 for double */ 
           one =  1, F_RL, F_col, F_one;
    double x   = 1.0, F_x;
    char  format[10];

    ++column; 
    if (digits > 1) {          /* text */
        sprintf(format, "%%%d.%dE", digits+7, digits);
        fprintf(fp, "%8d%8d%8d\n", column, 1, 1);
        fprintf(fp, format, x);
        fprintf(fp, "\n");
    } else if (digits == -1) { /* binary, opposite endian */
        F_RL  = flip_bytes_int(RL);
        F_one = flip_bytes_int(one);
        F_col = flip_bytes_int(column);
        F_x   = flip_bytes_double(x);
if (DEBUG)
printf("op4_wrt_trailer loc=%ld\n", ftell(fp));
        fwrite(&F_RL,              BYTES_PER_WORD,     1, fp);
        fwrite(&F_col,             BYTES_PER_WORD,     1, fp);
        fwrite(&F_one,             BYTES_PER_WORD,     1, fp);
        fwrite(&F_one,             BYTES_PER_WORD,     1, fp);
        fwrite(&F_x,             2*BYTES_PER_WORD,     1, fp);
        fwrite(&F_RL,              BYTES_PER_WORD,     1, fp);
    } else {                   /* binary, native endian   */
        fwrite(&RL,                BYTES_PER_WORD,     1, fp);
        fwrite(&column,            BYTES_PER_WORD,     1, fp);
        fwrite(&one,               BYTES_PER_WORD,     1, fp);
        fwrite(&one,               BYTES_PER_WORD,     1, fp);
        fwrite(&x,               2*BYTES_PER_WORD,     1, fp);
        fwrite(&RL,                BYTES_PER_WORD,     1, fp);
    }
    return 1;
} /* 1}}} */
void strings_in_list( /* {{{1 */
                     int  n_terms   ,     /* in  length of list[]            */
                     int *list      ,     /* in                              */
                     int *n_str     ,     /* out number of strings in list[] */
                     int *start_ind ,     /* out index of 1st string terms   */
                     int *str_len   ,     /* out length of each string       */
                     int *max_length)     /* out longest string in list[]    */
     /*  list[] contains non-negative integers in ascending order, eg:
      *     n_terms = 7
      *     list[0..6]  =  2, 5, 6, 7, 10, 12, 13
      *  This routine counts the number of blocks of consecutive terms
      *  and also returns the length of the longest string.  In the example
      *  the strings are   (2)  (5,6,7)  (10)  (12,13)
      *  so:
      *      n_str                   = 4 
      *      start_ind[0..(n_str-1)] =  0, 1, 4, 5
      *      str_len[  0..(n_str-1)] =  1, 3, 1, 2
      *      max_length              = 3
      */
{
int DEBUG = 0;
    int i, length_current;

    *n_str         = 0;
    *max_length    = 0;
    length_current = 0;
    if (n_terms > 0) {
        length_current    = 1;
        start_ind[*n_str] = 0;
        str_len[  *n_str] = length_current;
        ++(*n_str);
    } else {
        return;
    }
    for (i = 1; i < n_terms; i++) {
        if (list[i] != (list[i-1] + 1)) {
            length_current    = 1;
            start_ind[*n_str] = i;
            str_len[  *n_str] = length_current;
if (DEBUG) {
printf("list[%d]=%d  list[%d]=%d  nstr=%d   len=%d\n", 
i, list[i], i-1, list[i-1], *n_str, str_len[*n_str]);
}
            ++(*n_str);
        } else {
            ++length_current;
            str_len[*n_str-1] = length_current;
if (DEBUG) {
printf("           nstr=%d   len=%d\n", *n_str-1, str_len[*n_str-1]);
}
        }
        *max_length = MAX(*max_length, length_current);
    }
} /* 1}}} */

int    flip_bytes_int(int x) {  /* {{{1 */
    int   y;
    char *c_x, *c_y;

    c_x = (char *) &x;
    c_y = (char *) &y;

    c_y[0] = c_x[3];
    c_y[1] = c_x[2];
    c_y[2] = c_x[1];
    c_y[3] = c_x[0];

    return y;
} /* 1}}} */
float  flip_bytes_float(float x) {  /* {{{1 */
    /* Looks identical to flip_bytes_int() and it practice it might
       be possible to use this function and flip_bytes_int()
       interchangably.  The combination of gcc and x86 hardware may
       cause problems though since floats might be assigned extra
       bytes internally.
     */
    float   y;
    char    *c_x, *c_y;

    c_x = (char *) &x;
    c_y = (char *) &y;

    c_y[0] = c_x[3];
    c_y[1] = c_x[2];
    c_y[2] = c_x[1];
    c_y[3] = c_x[0];

    return y;
} /* 1}}} */
double flip_bytes_double(double x) {  /* {{{1 */
    double   y;
    char    *c_x, *c_y;

    c_x = (char *) &x;
    c_y = (char *) &y;

    c_y[0] = c_x[7];
    c_y[1] = c_x[6];
    c_y[2] = c_x[5];
    c_y[3] = c_x[4];
    c_y[4] = c_x[3];
    c_y[5] = c_x[2];
    c_y[6] = c_x[1];
    c_y[7] = c_x[0];

    return y;
} /* 1}}} */
int  op4_valid_name(char   *name)         /* {{{1 */
{
     /* Returns:
      *    1 if the name passed in is a valid Nastran op4 matrix name
      *   -1 if the name is null
      *    0 if the name is not valid
      */
	int DEBUG = 0;

    int i, len, fixed_len;
	if (DEBUG) printf("------------------------------------------------------\n");
	if (DEBUG) printf("1281 - op4_valid_name with [%s]\n", name);

    len = strlen(name);

    /* strip trailing blanks */
    fixed_len = len;
    for (i = len-1; i >= 0; i--) {
        if (ISBLANK(name[i]) || (name[i] == 24)) {
            /* binary op4 files use ASCII ord 24 for spaces */
            name[i] = 0;
            --fixed_len;
        } else {
            break;
        }
    }
    len = fixed_len;

    if (!len) {                            /* null string? */
        return -1;
    } else if (len > 8) {                  /* too long?    */
        return  0;
    } else if (!isalpha((int) name[0])) {  /* first char must be [a-zA-Z] */
        return  0;
    }

    for (i = 1; i < len; i++) {
        /* this test is too conservative; Nastran allows names like 'a+' */
        if (!(ISUNERDSCORE( name[i]) ||
              isalpha((int) name[i]) ||
              isdigit((int) name[i]))) {
            return 0;
        }
    }

    return 1;
} /* 1}}} */
