#ifndef __UTIL_H
#define __UTIL_H

#include <stdio.h>

#if SUPERLU > 0
/* Lib for SuperLU */
#include "slu_ddefs.h"
#endif

#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))
#define MAX3(a,b,c)		MAX(MAX(a,b),c)
#define MIN3(a,b,c)		MIN(MIN(a,b),c)
#define MID3(a,b,c)		((MIN(a,b)<(c))?(MIN(MAX(a,b),c)):(MAX(MIN(a,b),c)))
#define MAX4(a,b,c,d)	MAX(MAX(MAX(a,b),c),d)
#define MIN4(a,b,c,d)	MIN(MIN(MIN(a,b),c),d)
#define DELTA			1.0e-6
#define LARGENUM		1.0e100
#define NULLFILE		"(null)"


#define TRUE			1
#define	FALSE			0

#define RAND_SEED		1500450271

#define STR_SIZE		512
#define LINE_SIZE		65536
#define MAX_ENTRIES		512

int eq(double x, double y);
int le(double x, double y);
int ge(double x, double y);
void fatal (char *s);
void warning (char *s);
void swap_ival (int *a, int *b);
void swap_dval (double *a, double *b);
int tolerant_ceil(double val);
int tolerant_floor(double val);

/* vector routines	*/
double 	*dvector(int n);
void free_dvector(double *v);
void dump_dvector(double *v, int n);
void copy_dvector (double *dst, double *src, int n);
void zero_dvector (double *v, int n);

int *ivector(int n);
void free_ivector(int *v);
void dump_ivector(int *v, int n);
void copy_ivector (int *dst, int *src, int n);
void zero_ivector (int *v, int n);
/* sum of the elements	*/
double sum_dvector (double *v, int n);

/* matrix routines - Thanks to Greg Link
 * from Penn State University for the
 * memory allocators/deallocators
 */
double **dmatrix(int nr, int nc);
void free_dmatrix(double **m);
void dump_dmatrix(double **m, int nr, int nc);
void copy_dmatrix(double **dst, double **src, int nr, int nc);
void zero_dmatrix(double **m, int nr, int nc);
void resize_dmatrix(double **m, int nr, int nc);
/* mirror the lower triangle to make 'm' fully symmetric	*/
void mirror_dmatrix(double **m, int n);

int **imatrix(int nr, int nc);
void free_imatrix(int **m);
void dump_imatrix(int **m, int nr, int nc);
void copy_imatrix(int **dst, int **src, int nr, int nc);
void resize_imatrix(int **m, int nr, int nc);

/* routines for 3-d matrix with tail	*/
/* allocate 3-d matrix with 'nr' rows, 'nc' cols,
 * 'nl' layers	and a tail of 'xtra' elements
 */
double ***dcuboid_tail(int nr, int nc, int nl, int xtra);
/* destructor	*/
void free_dcuboid(double ***m);

/* initialize random number generator	*/
void init_rand(void);
/* random number within the range [0, max-1]	*/
int rand_upto(int max);
/* random number in the range [0, 1)	*/
double rand_fraction(void);

/* a table of name value pairs	*/
typedef struct str_pair_st
{
	char name[STR_SIZE];
	char value[STR_SIZE];
}str_pair;
/*
 * reads tab-separated name-value pairs from file into
 * a table of size max_entries and returns the number
 * of entries read successfully
 */
int read_str_pairs(str_pair *table, int max_entries, char *file);
/* same as above but from command line instead of a file	*/
int parse_cmdline(str_pair *table, int max_entries, int argc, char **argv);
/* append the table onto a file	*/
void dump_str_pairs(str_pair *table, int size, char *file, char *prefix);
/* table lookup	for a name */
int get_str_index(str_pair *table, int size, char *str);
/*
 * remove duplicate names in the table - the entries later
 * in the table are discarded. returns the new size of the
 * table
 */
int str_pairs_remove_duplicates(str_pair *table, int size);
/*
 * binary search a sorted double array 'arr' of size 'n'. if found,
 * the 'loc' pointer has the address of 'ele' and the return
 * value is TRUE. otherwise, the return value is FALSE and 'loc'
 * points to the 'should have been' location
 */
int bsearch_double(double *arr, int n, double ele, double **loc);
/*
 * binary search and insert an element into a partially sorted
 * double array if not already present. returns FALSE if present
 */
int bsearch_insert_double(double *arr, int n, double ele);

/* search if an array contains a value
 * return the index if so and -1 otherwise
 */
int contains(int *array, int size, int value);

/*
 * population count of an 8-bit integer - using pointers from
 * http://aggregate.org/MAGIC/
 */
unsigned int ones8(register unsigned char n);
/*
 * find the number of non-empty, non-comment lines
 * in a file open for reading
 */
int count_significant_lines(FILE *fp);

/* For Matrix format conversion */
int coo2csc(int size, int nnz,
            int *cooX, int *cooY, double *cooV,
            int *cscRowInd, int *cscColPtr, double *cscV);
int c2c_cmp( const void *a , const void *b);

void gaussj(double **a, int n, double *b);

#if SUPERLU > 0
/* Some Basic Matrix Data Structures and Operations */

/*
 * Structure to store a diagonal matrix
 */
typedef struct diagonal_matrix_t_st
{
  int n;
  double *vals;
} diagonal_matrix_t;

/*
 * computes A = c*diag + A
 */
int diagonal_add_SparseMatrix(double c, diagonal_matrix_t *diag, SuperMatrix *A);

/*
 * computes vector = c*diag*vector
 */
int diagonal_mul_vector(double c, diagonal_matrix_t *diag, double **vector);

/*
 * computes vector2 = c1*vector[1] + c2*vector[2]
 */
int vector_add_vector(int n, double c1, double *vector1, double c2, double *vector2);

/*
 * computes vector = A*vector
 */
int SparseMatrix_mul_vector(SuperMatrix *A, double *vector);

void cooTocsv(char *filename, int size, int nnz, int *cooX, int *cooY, double *cooV);
void diagTocsv(char *filename, diagonal_matrix_t *diag);
void vectorTocsv(char *filename, int size, double *vector);
#endif

#endif
