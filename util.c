#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include <strings.h>
#endif
#include <assert.h>

#include "util.h"

#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;}

int eq(double x, double y)
{
	return (fabs(x-y) <  DELTA);
}

int le(double x, double y)
{
	return ((x < y) || eq(x,y));
}

int ge(double x, double y)
{
	return ((x > y) || eq(x,y));
}

void fatal(char *s)
{
	fprintf(stderr, "error: %s", s);
	exit(1);
}

void warning(char *s)
{
	fprintf(stderr, "warning: %s", s);
}

void swap_ival (int *a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

void swap_dval (double *a, double *b)
{
	double t = *a;
	*a = *b;
	*b = t;
}

int tolerant_ceil(double val)
{
	double nearest = floor(val+0.5);
	/* numbers close to integers	*/
	if (eq(val, nearest))
		return ((int) nearest);
	/* all others	*/
	else
		return ((int) ceil(val));
}

int tolerant_floor(double val)
{
	double nearest = floor(val+0.5);
	/* numbers close to integers	*/
	if (eq(val, nearest))
		return ((int) nearest);
	/* all others	*/
	else
		return ((int) floor(val));
}

double *dvector(int n)
{
	double *v;

	v=(double *)calloc(n, sizeof(double));
	if (!v) fatal("allocation failure in dvector()\n");

	return v;
}

void free_dvector(double *v)
{
	free(v);
}

void dump_dvector (double *v, int n)
{
	int i;
	for (i=0; i < n; i++)
		fprintf(stdout, "%.5f\t", v[i]);
	fprintf(stdout, "\n");
}

void copy_dvector (double *dst, double *src, int n)
{
	memmove(dst, src, sizeof(double) * n);
}

void zero_dvector (double *v, int n)
{
	memset(v, 0, sizeof(double) * n);
}

/* sum of the elements	*/
double sum_dvector (double *v, int n)
{
	double sum = 0;
	int i;
	for(i=0; i < n; i++)
		sum += v[i];
	return sum;
}

int *ivector(int n)
{
	int *v;

	v = (int *)calloc(n, sizeof(int));
	if (!v) fatal("allocation failure in ivector()\n");

	return v;
}

void free_ivector(int *v)
{
	free(v);
}

void dump_ivector (int *v, int n)
{
	int i;
	for (i=0; i < n; i++)
		fprintf(stdout, "%d\t", v[i]);
	fprintf(stdout, "\n\n");
}

void copy_ivector (int *dst, int *src, int n)
{
	memmove(dst, src, sizeof(int) * n);
}

void zero_ivector (int *v, int n)
{
	memset(v, 0, sizeof(int) * n);
}

/*
 * Thanks to Greg Link from Penn State University
 * for these memory allocators/deallocators
 */
double **dmatrix(int nr, int nc)
{
	int i;
	double **m;

	m = (double **) calloc (nr, sizeof(double *));
	assert(m != NULL);
	m[0] = (double *) calloc (nr * nc, sizeof(double));
	assert(m[0] != NULL);

	for (i = 1; i < nr; i++)
    	m[i] =  m[0] + nc * i;

	return m;
}

void free_dmatrix(double **m)
{
	free(m[0]);
	free(m);
}

int **imatrix(int nr, int nc)
{
	int i;
	int **m;

	m = (int **) calloc (nr, sizeof(int *));
	assert(m != NULL);
	m[0] = (int *) calloc (nr * nc, sizeof(int));
	assert(m[0] != NULL);

	for (i = 1; i < nr; i++)
		m[i] = m[0] + nc * i;

	return m;
}

void free_imatrix(int **m)
{
	free(m[0]);
	free(m);
}

void dump_dmatrix (double **m, int nr, int nc)
{
	int i;
	for (i=0; i < nr; i++)
		dump_dvector(m[i], nc);
	fprintf(stdout, "\n");
}

void copy_dmatrix (double **dst, double **src, int nr, int nc)
{
	memmove(dst[0], src[0], sizeof(double) * nr * nc);
}

void zero_dmatrix(double **m, int nr, int nc)
{
	memset(m[0], 0, sizeof(double) * nr * nc);
}

void resize_dmatrix(double **m, int nr, int nc)
{
	int i;
	for (i = 1; i < nr; i++)
		m[i] = m[0] + nc * i;
}

/* allocate 3-d matrix with 'nr' rows, 'nc' cols,
 * 'nl' layers	and a tail of 'xtra' elements
 */
double ***dcuboid_tail(int nr, int nc, int nl, int xtra)
{
	int i, j;
	double ***m;

	/* 1-d array of pointers to the rows of the 2-d array below	*/
	m = (double ***) calloc (nl, sizeof(double **));
	assert(m != NULL);
	/* 2-d array of pointers denoting (layer, row)	*/
	m[0] = (double **) calloc (nl * nr, sizeof(double *));
	assert(m[0] != NULL);
	/* the actual 3-d data array	*/
	m[0][0] = (double *) calloc (nl * nr * nc + xtra, sizeof(double));
	assert(m[0][0] != NULL);

	/* remaining pointers of the 1-d pointer array	*/
	for (i = 1; i < nl; i++)
    	m[i] =  m[0] + nr * i;

	/* remaining pointers of the 2-d pointer array	*/
	for (i = 0; i < nl; i++)
		for (j = 0; j < nr; j++)
			/* to reach the jth row in the ith layer,
			 * one has to cross i layers i.e., i*(nr*nc)
			 * values first and then j rows i.e., j*nc
			 * values next
			 */
    		m[i][j] =  m[0][0] + (nr * nc) * i + nc * j;

	return m;
}

void free_dcuboid(double ***m)
{
	free(m[0][0]);
	free(m[0]);
	free(m);
}

/* mirror the lower triangle to make 'm' fully symmetric	*/
void mirror_dmatrix(double **m, int n)
{
	int i, j;
	for(i=0; i < n; i++)
		for(j=0; j < i; j++)
			m[j][i] = m[i][j];
}

void dump_imatrix (int **m, int nr, int nc)
{
	int i;
	for (i=0; i < nr; i++)
		dump_ivector(m[i], nc);
	fprintf(stdout, "\n");
}

void copy_imatrix (int **dst, int **src, int nr, int nc)
{
	memmove(dst[0], src[0], sizeof(int) * nr * nc);
}

void resize_imatrix(int **m, int nr, int nc)
{
	int i;
	for (i = 1; i < nr; i++)
		m[i] = m[0] + nc * i;
}

/* initialize random number generator	*/
void init_rand(void)
{
	srand(RAND_SEED);
}

/* random number within the range [0, max-1]	*/
int rand_upto(int max)
{
	return (int) (max * (double) rand() / (RAND_MAX+1.0));
}

/* random number in the range [0, 1)	*/
double rand_fraction(void)
{
	return ((double) rand() / (RAND_MAX+1.0));
}

/*
 * reads tab-separated name-value pairs from file into
 * a table of size max_entries and returns the number
 * of entries read successfully
 */
int read_str_pairs(str_pair *table, int max_entries, char *file)
{
	int i=0;
	char str[LINE_SIZE], copy[LINE_SIZE];
	char name[STR_SIZE];
	char *ptr;
	FILE *fp = fopen (file, "r");
	if (!fp) {
		sprintf (str,"error: %s could not be opened for reading\n", file);
		fatal(str);
	}
	while(i < max_entries) {
		fgets(str, LINE_SIZE, fp);
		if (feof(fp))
			break;
		strcpy(copy, str);

		/* ignore comments and empty lines  */
		ptr = strtok(str, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		if ((sscanf(copy, "%s%s", name, table[i].value) != 2) || (name[0] != '-'))
			fatal("invalid file format\n");
		/* ignore the leading "-"	*/
		strcpy(table[i].name, &name[1]);
		i++;
	}
	fclose(fp);
	return i;
}

/*
 * same as above but from command line instead of a file. the command
 * line is of the form <prog_name> <name-value pairs> where
 * <name-value pairs> is of the form -<variable> <value>
 */
int parse_cmdline(str_pair *table, int max_entries, int argc, char **argv)
{
	int i, count;
	for (i=1, count=0; i < argc && count < max_entries; i++) {
		if (i % 2) {	/* variable name	*/
			if (argv[i][0] != '-')
				fatal("invalid command line. check usage\n");
			/* ignore the leading "-"	*/
			strncpy(table[count].name, &argv[i][1], STR_SIZE-1);
			table[count].name[STR_SIZE-1] = '\0';
		} else {		/* value	*/
			strncpy(table[count].value, argv[i], STR_SIZE-1);
			table[count].value[STR_SIZE-1] = '\0';
			count++;
		}
	}
	return count;
}

/* append the table onto a file	*/
void dump_str_pairs(str_pair *table, int size, char *file, char *prefix)
{
	int i;
	char str[STR_SIZE];
	FILE *fp = fopen (file, "w");
	if (!fp) {
		sprintf (str,"error: %s could not be opened for writing\n", file);
		fatal(str);
	}
	for(i=0; i < size; i++)
		fprintf(fp, "%s%s\t%s\n", prefix, table[i].name, table[i].value);
	fclose(fp);
}

/* table lookup	for a name */
int get_str_index(str_pair *table, int size, char *str)
{
	int i;

	if (!table)
		fatal("null pointer in get_str_index\n");

	for (i = 0; i < size; i++)
		if (!strcasecmp(str, table[i].name))
			return i;
	return -1;
}

/* delete entry at 'at'	*/
void delete_entry(str_pair *table, int size, int at)
{
	int i;
	/*
	 * overwrite this entry using the next and
	 * shift all later entries once
	 */
	for (i=at+1; i < size; i++) {
		strcpy(table[i-1].name, table[i].name);
		strcpy(table[i-1].value, table[i].value);
	}
}

/*
 * remove duplicate names in the table - the entries later
 * in the table are discarded. returns the new size of the
 * table
 */
int str_pairs_remove_duplicates(str_pair *table, int size)
{
	int i, j;

	for(i=0; i < size-1; i++)
		for(j=i+1; j < size; j++)
			if (!strcasecmp(table[i].name, table[j].name)) {
				delete_entry(table, size, j);
				size--;
				j--;
			}
	return size;
}

/* debug	*/
void print_str_pairs(str_pair *table, int size)
{
	int i;
	fprintf(stdout, "printing string table\n");
	for (i=0; i < size; i++)
		fprintf(stdout, "%s\t%s\n", table[i].name, table[i].value);
}

/*
 * binary search a sorted double array 'arr' of size 'n'. if found,
 * the 'loc' pointer has the address of 'ele' and the return
 * value is TRUE. otherwise, the return value is FALSE and 'loc'
 * points to the 'should have been' location
 */
int bsearch_double(double *arr, int n, double ele, double **loc)
{
	if(n < 0)
		fatal("wrong index in binary search\n");

	if(n == 0) {
		*loc = arr;
		return FALSE;
	}

	if(eq(arr[n/2], ele)) {
		*loc = &arr[n/2];
		return TRUE;
	} else if (arr[n/2] < ele) {
		return bsearch_double(&arr[n/2+1], (n-1)/2, ele, loc);
	} else
		return bsearch_double(arr, n/2, ele, loc);

}

/*
 * binary search and insert an element into a partially sorted
 * double array if not already present. returns FALSE if present
 */
int bsearch_insert_double(double *arr, int n, double ele)
{
	double *loc;
	int i;

	/* element found - nothing more left to do	*/
	if (bsearch_double(arr, n, ele, &loc))
		return FALSE;
	else {
		for(i=n-1; i >= (loc-arr); i--)
			arr[i+1] = arr[i];
		arr[loc-arr] = ele;
	}
	return TRUE;
}

/* search if an array contains a value
 * return the index if so and -1 otherwise
 */
int contains(int *array, int size, int value)
{
  int i;
  for(i = 0; i < size; i++)
  {
    if(array[i] == value)
      return i;
  }

  return -1;
}

/*
 * population count of an 8-bit integer - using pointers from
 * http://aggregate.org/MAGIC/
 */
unsigned int ones8(register unsigned char n)
{
	/* group the bits in two and compute the no. of 1's within a group
	 * this works because 00->00, 01->01, 10->01, 11->10 or
	 * n = n - (n >> 1). the 0x55 masking prevents bits flowing across
	 * group boundary
	 */
	n -= ((n >> 1) & 0x55);
	/* add the 2-bit sums into nibbles */
	n = ((n & 0x33) + ((n >> 2) & 0x33));
	/* add both the nibbles */
	n = ((n + (n >> 4)) & 0x0F);
	return n;
}

/*
 * find the number of non-empty, non-comment lines
 * in a file open for reading
 */
int count_significant_lines(FILE *fp)
{
    char str[LINE_SIZE], *ptr;
    int count = 0;

	fseek(fp, 0, SEEK_SET);
	while(!feof(fp)) {
		fgets(str, LINE_SIZE, fp);
		if (feof(fp))
			break;

		/* ignore comments and empty lines	*/
		ptr = strtok(str, " \r\t\n");
		if (!ptr || ptr[0] == '#')
			continue;

		count++;
	}
	return count;
}

/* Ke's code: Coo2CSC */
struct coo_elem
{
  int x;
  int y;
  double val;
};

int c2c_cmp( const void *a , const void *b )
{
  struct coo_elem *c = (struct coo_elem *)a;
  struct coo_elem *d = (struct coo_elem *)b;
  if(c->y != d->y) return c->y - d->y;
  else return c->x - d->x;
}

int coo2csc(int size, int nnz,
            int *cooX, int *cooY, double *cooV,  // input COO array
            int *cscRowInd, int *cscColPtr, double *cscV) //output CSC array
{
  int i, j;
  int prev_x, prev_y;
  // Init struct array
  struct coo_elem *cooArray;
  cooArray = (struct coo_elem *) calloc (nnz, sizeof(struct coo_elem));

  // Copy in
  for (i =0; i <nnz; i++) {
      cooArray[i].x = cooX[i];
      cooArray[i].y = cooY[i];
      cooArray[i].val = cooV[i];
  }

  // Sort in col major
  qsort(cooArray, nnz, sizeof(cooArray[0]), c2c_cmp);

  // Copy out, check duplicate
  j = -1;
  prev_x = -1;
  prev_y = -1;
  for (i =0; i <nnz; i++) {
      cscRowInd[i]=cooArray[i].x;
      cscV[i]=cooArray[i].val;
      while(j<cooArray[i].y){
          j++;
          cscColPtr[j]=i;
      }
      if((cooArray[i].x == prev_x) &&
         (cooArray[i].y == prev_y))
        printf("Warning: Duplicate elements in Matrix!\n");

      prev_x = cooArray[i].x;
      prev_y = cooArray[i].y;
  }
  cscColPtr[j+1]=i;

  free(cooArray);

  return 1;
}

/*
 * Gauss-Jordan elimination with full pivoting
 * Taken from Numerical Recipes in C
 * gaussj is written assuming that b is a vector, but it is trivial to instead
 * make b an nxm matrix in order to solve Ax=b for m different values of b
 */
void gaussj(double **a, int n, double *b) {
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;

  indxc = calloc(n, sizeof(int));
  indxr = calloc(n, sizeof(int));
  ipiv  = calloc(n, sizeof(int));

  // Main Loop
  for(i = 0; i < n; i++) {
    big = 0.0;

    // Search for a pivot element (choose the largest)
    for(j = 0; j < n; j++) {
      if (ipiv[j] != 1)
        for(k = 0; k < n; k++) {
          if(ipiv[k] == 0) {
            if(fabs(a[j][k]) >= big) {
              big = fabs(a[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
    }
    ++(ipiv[icol]);

    /*
     * We've found the pivot, so we interchange rows if necessary to put
     * the pivot on the diagonal
     * indxc[i] = the column of the ith pivot element
     * indxr[i] = the row in which that pivot element was originally located
     * if indxr[i] != indxc[i], there is an implied column intercahnge
     */
    if(irow != icol) {
      for(l = 0; l < n; l++)
        SWAP(a[irow][l], a[icol][l])
      SWAP(b[irow], b[icol])
    }

    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.0)
      fatal("gaussj: Singular Pressure Matrix\n");

    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;

    // Divide the pivot row by the pivot element
    for(l = 0; l < n; l++)
      a[icol][l] *= pivinv;
    b[icol] *= pivinv;

    // Reduce all rows except the row containing the pivot
    for(ll = 0; ll < n; ll++) {
      if(ll != icol) {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for(l = 0; l < n; l++)
          a[ll][l] -= a[icol][l] * dum;
        b[ll] -= b[icol] * dum;
      }
    }
  } // end of Main Loop

  // Put columns back in original order
  for(l = n-1; l >= 0; l--) {
    if(indxr[l] != indxc[l]) {
      for(k = 1; k < n; k++) {
        SWAP(a[k][indxr[l]], a[k][indxc[l]]);
      }
    }
  }

  free(ipiv);
  free(indxr);
  free(indxc);
}

#if SUPERLU > 0
/*
 * computes A = c*diag + A
 * NOTE: Assumes that A contains only nonzero elements on its diagonal
 */
int diagonal_add_SparseMatrix(double c, diagonal_matrix_t *diag, SuperMatrix *A) {
  NCformat *Astore;
  int i, j;
  int n;
  double *a;
  int *asub, *xa;
  int flag = 1;

  Astore = A->Store;
  n = A->ncol;
  a = Astore->nzval;
  asub = Astore->rowind;
  xa = Astore->colptr;

  double *diag_vals = diag->vals;

  for(i=0; i<n; i++){
      flag = 1;
      j = xa[i];
      while(j<xa[i+1]){
          if(asub[j] == i){
              a[j] += c*diag_vals[i];
              flag = 0;
              fprintf(stderr, "A[%d] = %e\n", j, a[j]);
          }
          j++;
      }
      if(flag)
        fatal("Matrix missing diagonal element! Cannot support that yet\n");
  }
  return 1;
}

/*
 * computes vector = c*diag*vector
 */
int diagonal_mul_vector(double c, diagonal_matrix_t *diag, double **vector) {
  int i, n;
  double *diag_vals;

  n = diag->n;
  diag_vals = diag->vals;

  for(i = 0; i < n; i++) {
    (*vector)[i] *= c*diag_vals[i];
  }

  return 1;
}

/*
 * computes vector2 += vector1
 */
int vector_add_vector(int n, double c1, double *vector1, double c2, double *vector2) {
  int i;

  for(i = 0; i < n; i++) {
    //fprintf(stderr, "vector2[%d] = vector1[%d] + vector[%d] = %e + %e = ", i, i, i, vector1[i], vector2[i]);
    vector2[i] = c1*vector1[i] + c2*vector2[i];
    //fprintf(stderr, "%e\n", vector2[i]);
  }

  return 1;
}

/*
 * computes A = A*vector
 */
int SparseMatrix_mul_vector(SuperMatrix *A, double *vector) {
  NCformat *Astore;
  double *result;
  int i, j, row_index;
  int m, n;
  double *a;
  int *asub, *xa;

  m = A->nrow;
  n = A->ncol;
  Astore = A->Store;
  a = Astore->nzval;
  asub = Astore->rowind;
  xa = Astore->colptr;

  if ( !(result = (double *) calloc(m, sizeof(double))) )
    fatal("Malloc fails for result[].\n");

  for(i=0; i<m; i++)
    result[i] = 0;

  for(i=0; i<n; i++){
      j = xa[i];
      while(j<xa[i+1]){
          row_index = asub[j];
          result[row_index] += a[j] * vector[i];
          j++;
      }
  }

  copy_dvector(vector, result, m);

  free(result);

  return 1;
}

void cooTocsv(char *filename, int size, int nnz, int *cooX, int *cooY, double *cooV) {
  int i, j, k;
  double **matrix;

  matrix = calloc(size, sizeof(double *));

  for(i = 0; i < size; i++)
    matrix[i] = calloc(size, sizeof(double));

  for(i = 0; i < nnz; i++) {
    matrix[cooX[i]][cooY[i]] = cooV[i];
  }

  FILE *fp = fopen(filename, "w");

  fprintf(fp, ",");
  for(i = 0; i < size-1; i++)
    fprintf(fp, "%d, ", i);

  fprintf(fp, "%d\n", size-1);

  for(i = 0; i < size; i++) {
    fprintf(fp, "%d, ", i);
    for(j = 0; j < size-1; j++) {
      fprintf(fp, "%e, ", matrix[i][j]);
    }
    fprintf(fp, "%e\n", matrix[i][size-1]);
  }

  for(i = 0; i < size; i++)
    free(matrix[i]);

  free(matrix);

  fclose(fp);
}

void diagTocsv(char *filename, diagonal_matrix_t *diag) {
  int n = diag->n;
  double *vals = diag->vals;
  int i, j;

  FILE *fp = fopen(filename, "w");

  fprintf(fp, ",");
  for(i = 0; i < n-1; i++)
    fprintf(fp, "%d, ", i);

  fprintf(fp, "%d\n", n-1);

  for(i = 0; i < n; i++) {
    fprintf(fp, "%d, ", i);
    for(j = 0; j < n-1; j++) {
      if(i != j)
        fprintf(fp, "0, ");
      else
        fprintf(fp, "%e, ", vals[i]);
    }

    if(i == n-1)
      fprintf(fp, "%e\n", vals[n-1]);
    else
      fprintf(fp, "0\n");
  }

  fclose(fp);
}

void vectorTocsv(char *filename, int size, double *vector) {
  FILE *fp = fopen(filename, "w");
  int i;

  fprintf(fp, ",0\n");
  for(i = 0; i < size; i++)
    fprintf(fp, "%d, %e\n", i, vector[i]);

  fclose(fp);
}

#endif
