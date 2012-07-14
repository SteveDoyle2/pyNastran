/* {{{1 GNU Lesser General Public License

Copyright (C) 2001-2012  <al.danial@gmail.com>

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

#define OP4_TXT_LINE_SIZE    82  /* have to include the \n */
#define Max_Matrices_Per_OP4 40
#define OP4_FMT_STR_SIZE     10  /* max chars in, eg, 1PE23.16 */
#define OP4_MAX_TEXT_COLS    10  /* at most 10 text columns */
#define OP4_FMT_LINE_SIZE    OP4_MAX_TEXT_COLS*OP4_FMT_STR_SIZE
#define BYTES_PER_WORD        4
#define ISBLANK(a)      ((a) == ' ') ? 1 : 0
#define ISUNERDSCORE(a) ((a) == '_') ? 1 : 0

typedef struct {
    int     len       ; /* Number of terms in the string (a complex       */
                        /* number counts as a single term).               */
    int     start_row ; /* Zero based, so first row has start_row = 0.    */
    int     N_idx     ; /* Index into N[] to first numeric term for this  */
                        /* string.                                        */
} str_t;

typedef struct {
    int     *magic    ;
    int     *H        ;         /* size = SPARSE_HDR_SIZE  */
    int     *S_start  ;         /* size = # columns + 1    */
    int     *N_start  ;         /* size = # columns + 1    */
    str_t   *S        ;         /* size = n_strings_total  */
    double  *N        ;         /* size = n_nonzeros_total */

    char    *data     ;         /* main data block         */

} SparseMatrix;

int op4_filetype(const char *filename);
int Op4_scan(  const char *filename  ); /* TESTING ONLY                */
int op4_scan(  const char *filename  ,  /* in                          */
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
              );
int op4_scan_text(const char *filename  ,  /* in                          */
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
              );
int op4_scan_binary(const char *filename  ,  /* in                          */
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
              );
void op4_is_mat_header_binary(FILE *fp        ,
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
                         long *offset);     /* out byte offset to matrix    */
int op4_count_str_binary(FILE *fp     ,  /* in  */
                    int   endian ,  /* in  */
                    int   nType  ,  /* in  */
                    int   storage,  /* in  */
                    int   nRows  ,  /* in  */
                    int   nCols  ,  /* in  */
                    int  *nStr   ,  /* out */
                    int  *nNnz);    /* out */
int op4_line_type(const char *line);
FILE* op4_open_r(const char *filename, long offset);
int op4_read_col_text(FILE   *fp         ,
                   int     c_in       ,  /* in  requested column to read   */
                   int     nRow       ,  /* in  # rows    in matrix        */
                   int     nCol       ,  /* in  # columns in matrix        */
                   char   *fmt_str    ,  /* in  eg "%23le%23le%23le"       */
                   int     col_width  ,  /* in  # characters in format str */
                   int     storage    ,  /* in  0=dense  1,2=sparse  3=ccr */
                   int     complx     ,  /* in  0=real   1=complex         */
                   int    *n_str      ,  /* out # strings   (s_o) = 1,2    */
                   str_t  *str_data   ,  /* out string data (s_o) = 1,2    */
                   int    *N_index    ,  /* in/out          (s_o) = 1,2    */
                   double *N             /* out numeric data               */
                  );
int  op4_read_col_binary(FILE   *fp         ,
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
                   );
int    flip_bytes_int(   int    x);
float  flip_bytes_float( float  x);
double flip_bytes_double(double x);
int  op4_wrt_header(FILE   *fp         ,
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
                   );
int  op4_wrt_col_dn(FILE   *fp    ,
                    int     column, /* first column is 0       */
                    int     nRows , /* number of rows          */
                    int     nCols , /* number of columns       */
                    double *A     , /* array of terms to write */
                    int     complx, /* 1=complex  0=real       */
                    int     digits  /* -1=flipped endian; 0=native; >0=digits */
                    );
int  op4_wrt_col_sp(FILE         *fp    ,
                    int           column, /* first column is 0    */
                    int           A_col , /* column index to A[]; this differs
                                             from 'column' if A[] contains only
                                             part of the entire matrix */
                    int           nCols , /* number of columns    */
                    SparseMatrix  A     , /* entire sparse matrix */
                    int           complx, /* 1=complex  0=real    */
                    int           digits  /* -1=flipped; 0=native; >0=digits */
                    );
int  op4_wrt_trailer(FILE *fp    ,
                     int   column, /* first column is 0    */
                     int   digits  /* -1=flipped; 0=native; >0=digits */
                    );

enum { OP4_TEXT_NUMERIC = 0,
       OP4_TEXT_HEADER     ,
       OP4_TEXT_NEW_COLUMN ,
       OP4_TEXT_STRING_2   ,
       OP4_TEXT_STRING_1   ,
       OP4_TEXT_ERROR
     };

/* The items of this enum comprise the sparse matrix header. */
enum { ROWS              = 0,

       COLS                 ,

       n_STR                ,

       n_NONZ               ,   /* number of non-zero *terms*; a complex *
                                 * value is a single term                */
       COMPLX               ,   /* 0 => real          1 => complex       */

       DATA_SIZE                /* in bytes; total length of .data[]        */
                                /* This should always stay as the last term */
                                /* of the header since the constant         */
                                /* SPARSE_HDR_SIZE is defined from it.      */
};
