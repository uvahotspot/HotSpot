/*
 * Thanks to Greg Link from Penn State University
 * for his math acceleration engine. Where available,
 * the modified version of the engine found here, uses
 * the fast, vendor-provided linear algebra routines
 * from the BLAS and LAPACK packages in  lieu of
 * the vanilla C code present in the matrix functions
 * of the previous versions of HotSpot.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "temperature.h"
#include "flp.h"
#include "util.h"

/* thermal resistance calculation	*/
double getr(double conductivity, double thickness, double area)
{
	return thickness / (conductivity * area);
}

/* thermal capacitance calculation	*/
double getcap(double sp_heat, double thickness, double area)
{
	/* include lumped vs. distributed correction	*/
	return C_FACTOR * sp_heat * thickness * area;
  //return sp_heat * thickness * area;
}

/*
 * LUP decomposition from the pseudocode given in the CLR
 * 'Introduction to Algorithms' textbook. The matrix 'a' is
 * transformed into an in-place lower/upper triangular matrix
 * and the vector'p' carries the permutation vector such that
 * Pa = lu, where 'P' is the matrix form of 'p'. The 'spd' flag
 * indicates that 'a' is symmetric and positive definite
 */

void lupdcmp(double**a, int n, int *p, int spd)
{
	#if(MATHACCEL == MA_INTEL)
	int info = 0;
	if (!spd)
		dgetrf(&n, &n, a[0], &n, p, &info);
	else
		dpotrf("U", &n, a[0], &n, &info);
	assert(info == 0);
	#elif(MATHACCEL == MA_AMD)
	int info = 0;
	if (!spd)
		dgetrf_(&n, &n, a[0], &n, p, &info);
	else
		dpotrf_("U", &n, a[0], &n, &info, 1);
	assert(info == 0);
	#elif(MATHACCEL == MA_APPLE)
	int info = 0;
	if (!spd)
		dgetrf_((__CLPK_integer *)&n, (__CLPK_integer *)&n, a[0],
				(__CLPK_integer *)&n, (__CLPK_integer *)p,
				(__CLPK_integer *)&info);
	else
		dpotrf_("U", (__CLPK_integer *)&n, a[0], (__CLPK_integer *)&n,
				(__CLPK_integer *)&info);
	assert(info == 0);
	#elif(MATHACCEL == MA_SUN)
	int info = 0;
	if (!spd)
		dgetrf_(&n, &n, a[0], &n, p, &info);
	else
		dpotrf_("U", &n, a[0], &n, &info);
	assert(info == 0);
	#else
	int i, j, k, pivot=0;
	double max = 0;

	/* start with identity permutation	*/
	for (i=0; i < n; i++)
		p[i] = i;

	for (k=0; k < n-1; k++)	 {
		max = 0;
		for (i = k; i < n; i++)	{
			if (fabs(a[i][k]) > max) {
				max = fabs(a[i][k]);
				pivot = i;
			}
		}
		if (eq (max, 0))
			fatal ("singular matrix in lupdcmp\n");

		/* bring pivot element to position	*/
		swap_ival (&p[k], &p[pivot]);
		for (i=0; i < n; i++)
			swap_dval (&a[k][i], &a[pivot][i]);

		for (i=k+1; i < n; i++) {
			a[i][k] /= a[k][k];
			for (j=k+1; j < n; j++)
				a[i][j] -= a[i][k] * a[k][j];
		}
	}
	#endif
}

/*
 * the matrix a is an in-place lower/upper triangular matrix
 * the following macros split them into their constituents
 */

#define LOWER(a, i, j)		((i > j) ? a[i][j] : 0)
#define UPPER(a, i, j)		((i <= j) ? a[i][j] : 0)

/*
 * LU forward and backward substitution from the pseudocode given
 * in the CLR 'Introduction to Algorithms' textbook. It solves ax = b
 * where, 'a' is an in-place lower/upper triangular matrix. The vector
 * 'x' carries the solution vector. 'p' is the permutation vector. The
 * 'spd' flag indicates that 'a' is symmetric and positive definite
 */

void lusolve(double **a, int n, int *p, double *b, double *x, int spd)
{
	#if(MATHACCEL == MA_INTEL)
	int one = 1, info = 0;
	cblas_dcopy(n, b, 1, x, 1);
	if (!spd)
		dgetrs("T", &n, &one, a[0], &n, p, x, &n, &info);
	else
		dpotrs("U", &n, &one, a[0], &n, x, &n, &info);
	assert(info == 0);
	#elif(MATHACCEL == MA_AMD)
	int one = 1, info = 0;
	dcopy(n, b, 1, x, 1);
	if (!spd)
		dgetrs_("T", &n, &one, a[0], &n, p, x, &n, &info, 1);
	else
		dpotrs_("U", &n, &one, a[0], &n, x, &n, &info, 1);
	assert(info == 0);
	#elif(MATHACCEL == MA_APPLE)
	int one = 1, info = 0;
	cblas_dcopy(n, b, 1, x, 1);
	if (!spd)
		dgetrs_("T", (__CLPK_integer *)&n, (__CLPK_integer *)&one, a[0],
				(__CLPK_integer *)&n, (__CLPK_integer *)p, x,
				(__CLPK_integer *)&n, (__CLPK_integer *)&info);
	else
		dpotrs_("U", (__CLPK_integer *)&n, (__CLPK_integer *)&one, a[0],
				(__CLPK_integer *)&n, x, (__CLPK_integer *)&n,
				(__CLPK_integer *)&info);
	assert(info == 0);
	#elif(MATHACCEL == MA_SUN)
	int one = 1, info = 0;
	dcopy(n, b, 1, x, 1);
	if (!spd)
		dgetrs_("T", &n, &one, a[0], &n, p, x, &n, &info);
	else
		dpotrs_("U", &n, &one, a[0], &n, x, &n, &info);
	assert(info == 0);
	#else
	int i, j;
	double *y = dvector (n);
	double sum;

	/* forward substitution	- solves ly = pb	*/
	for (i=0; i < n; i++) {
		for (j=0, sum=0; j < i; j++)
			sum += y[j] * LOWER(a, i, j);
		y[i] = b[p[i]] - sum;
	}

	/* backward substitution - solves ux = y	*/
	for (i=n-1; i >= 0; i--) {
		for (j=i+1, sum=0; j < n; j++)
			sum += x[j] * UPPER(a, i, j);
		x[i] = (y[i] - sum) / UPPER(a, i, i);
	}

	free_dvector(y);
	#endif
}

#if SUPERLU > 0
// Builds the matrix A = (1/h)C + G for the Backward Euler Method
int build_A_matrix(SuperMatrix *G, diagonal_matrix_t *C, double h, SuperMatrix *A)
{
  NCformat *Astore;
  int i, j;
  int n = C->n;
  double *a;
  int *asub, *xa;
  int flag = 1;

  NCformat *Gstore = G->Store;
  int nrow = G->nrow;
  int ncol = G->ncol;
  int nnz = Gstore->nnz;

  // Copy G into A
  double *nzval;
  int *rowind;
  int *colptr;
  if ( !(nzval = doubleMalloc(nnz)) ) fatal("Malloc fails for a[].\n");
  if ( !(rowind = intMalloc(nnz)) ) fatal("Malloc fails for asub[].\n");
  if ( !(colptr = intMalloc(n+1)) ) fatal("Malloc fails for xa[].\n");


  dCreate_CompCol_Matrix(A, nrow, ncol, nnz, nzval, rowind,
                         colptr, SLU_NC, SLU_D, SLU_GE);


  dCopy_CompCol_Matrix(G, A);

  Astore = A->Store;
  n = A->ncol;
  a = Astore->nzval;
  asub = Astore->rowind;
  xa = Astore->colptr;

  double *C_vals = C->vals;

  for(i=0; i<n; i++){
      flag = 1;
      j = xa[i];
      while(j<xa[i+1]){
          if(asub[j] == i){
              //fprintf(stderr, "i = %d, j = %d, a[%d]  = %e + (1.0 / %e) * %e = ", i, j, j, a[j], h, C_vals[i]);
              a[j] += (1.0 / h)*C_vals[i];
              //fprintf(stderr, "%e\n", a[j]);
              flag = 0;
          }
          j++;
      }
      if(flag)
        return 0;
  }
  return 1;
}

// Builds the matrix B = (1/h)CT + P for the Backward Euler method
int build_B_matrix(diagonal_matrix_t *C, double *T, double *P, double h, SuperMatrix *B)
{
  int i, n = C->n;
  double *B_temp = dvector(n);
  double *C_vals = C->vals;

  for(i = 0; i < n; i++) {
    B_temp[i] = (1.0 / h) * C_vals[i] * T[i] + P[i];
  }

  dCreate_Dense_Matrix(B, n, 1, B_temp, n, SLU_DN, SLU_D, SLU_GE);

  return 1;
}

/*
 * Backward Euler solver with adaptive step sizing.
 * It takes as input the conductance matrix G, capacitances C,
 * temperatures T and power P at time t, and solves the ODE GT + CdT = P
 * between t and t+h. It returns the correct step size to be used next time.
 */
#define BACKWARD_EULER_SAFETY    0.95
#define BACKWARD_EULER_MAXUP     5.0
#define BACKWARD_EULER_MAXDOWN   10.0
#define BACKWARD_EULER_PRECISION 0.01
double backward_euler(SuperMatrix *G, diagonal_matrix_t *C, double *T, double *P, double *h, double *Tout)
{
  ////////////////// Set Options and Statistics //////////////////////////
  SuperMatrix L, U;
  int *perm_r; // row permutations from partial pivoting
  int *perm_c; // column permutation vector
  int info, n = C->n;
  double max, new_h = (*h);
  superlu_options_t options;
  SuperLUStat_t stat;
  DNformat     *Bstore;
  double       *dp;

  // TODO: Implement re-use of perm_c
  if ( !(perm_r = intMalloc(n)) ) fatal("Malloc fails for perm_r[].\n");
  if ( !(perm_c = intMalloc(n)) ) fatal("Malloc fails for perm_c[].\n");

  set_default_options(&options);

  // Initialize the statistics variables.
  StatInit(&stat);

  SuperMatrix *A, *B;
  int i;
  double *Ttemp = dvector(n);
  double *T1    = dvector(n);
  double *T2    = dvector(n);


  A = calloc(1, sizeof(SuperMatrix));
  B = calloc(1, sizeof(SuperMatrix));

  int flag = 1;

  if(build_A_matrix(G, C, *h, A) != 1)
    fatal("error\n");

  if(build_B_matrix(C, T, P, *h, B) != 1)
    fatal("error\n");

  dgssv(&options, A, perm_c, perm_r, &L, &U, B, &stat, &info);

  // Copy results into Ttemp
  Bstore = (DNformat *) B->Store;
  dp = (double *) Bstore->nzval;
  for(i = 0; i < n; i++) {
    Ttemp[i] = dp[i];
    //fprintf(stderr, "  Ttemp[%d] = %e\n", i, dp[i]);
  }

  // TEST
  /*
  do {
    (*h) = new_h;

    ////////////////////// Evaluate once with step size h ////////////////////
    // Build the matrix A = (1/h)C + G
    if(build_A_matrix(G, C, *h, A) != 1)
      fatal("In backward_euler(): Failed to build A matrix\n");

    // Build the matrix B = (1/h)CT + P
    if(build_B_matrix(C, T, P, *h, B) != 1)
      fatal("In backward_euler(): Failed to build B matrix\n");

    // Solve the linear system Ax = B
    dgssv(&options, A, perm_c, perm_r, &L, &U, B, &stat, &info);

    // Copy results into Ttemp
    Bstore = (DNformat *) B->Store;
    dp = (double *) Bstore->nzval;
    for(i = 0; i < n; i++) {
      Ttemp[i] = dp[i];
      //fprintf(stderr, "  Ttemp[%d] = %e\n", i, dp[i]);
    }

    ////////////////////// Repeat as two steps of size h/2 ////////////////////
    // Build the matrix A = (2/h)C + G
    if(build_A_matrix(G, C, (*h)/2, A) != 1)
      fatal("In backward_euler(): Failed to build A matrix\n");

    // Build the matrix B = (2/h)CT + P
    if(build_B_matrix(C, T, P, (*h)/2, B) != 1)
      fatal("In backward_euler(): Failed to build B matrix\n");

    // Solve the linear system Ax = B
    dgssv(&options, A, perm_c, perm_r, &L, &U, B, &stat, &info);

    // Copy results into T1
    Bstore = (DNformat *) B->Store;
    dp = (double *) Bstore->nzval;
    for(i = 0; i < n; i++) {
      T1[i] = dp[i];
      //fprintf(stderr, "  T1[%d] = %e\n", i, T1[i]);
    }

    // Rebuild B = (2/h)CT + P using T1
    if(build_B_matrix(C, T1, P, (*h)/2, B) != 1)
      fatal("In backward_euler(): Failed to build B matrix\n");

    // Solve the linear system Ax = B
    dgssv(&options, A, perm_c, perm_r, &L, &U, B, &stat, &info);

    // Copy results into T2
    Bstore = (DNformat *) B->Store;
    dp = (double *) Bstore->nzval;
    for(i = 0; i < n; i++) {
      T2[i] = dp[i];
      //fprintf(stderr, "  T2[%d] = %e\n", i, T2[i]);
    }

    // Find maximum difference between Ttemp and T2
    for(i = 0; i < n; i++) {
      T1[i] = fabs(Ttemp[i] - T2[i]);
      //fprintf(stderr, "    delta[%d] = |%e - %e| = %e\n", i, Ttemp[i], T2[i], T1[i]);
    }
    max = T1[0];
    for(i = 1; i < n; i++) {
      if(max < T1[i]) {
        max = T1[i];
      }
    }
    //fprintf(stderr, "    max delta = %e\n", max);

    // Accuracy OK. Increase step size
    if(max <= BACKWARD_EULER_PRECISION) {
      new_h = BACKWARD_EULER_SAFETY * (*h) * pow(fabs(BACKWARD_EULER_PRECISION/max), 0.2);
      if (new_h > BACKWARD_EULER_MAXUP * (*h))
        new_h = BACKWARD_EULER_MAXUP * (*h);
    // Inaccuracy error. Decrease step size	and compute again
      //fprintf(stderr, "  max delta OK; increasing step size to h = %e\n", new_h);
      flag = 0;
    }
    else {
      new_h = BACKWARD_EULER_SAFETY * (*h) * pow(fabs(BACKWARD_EULER_PRECISION/max), 0.25);
        if (new_h < (*h) / BACKWARD_EULER_MAXDOWN)
          new_h = (*h) / BACKWARD_EULER_MAXDOWN;

      //fprintf(stderr, "  max delta NOT OK; decreasing step size to h = %e\n", new_h);
    }

  } while(new_h < (*h));
*/
  // Copy temperatures over to output
  for(i = 0; i < n; i++) {
    Tout[i] = Ttemp[i];
    //fprintf(stderr, "  Returning T[%d] = %e\n", i, Tout[i]);
  }

  free_dvector(Ttemp);
  free_dvector(T1);
  free_dvector(T2);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(A);
  Destroy_SuperMatrix_Store(B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);

  return new_h;
}
#endif

/* core of the 4th order Runge-Kutta method, where the Euler step
 * (y(n+1) = y(n) + h * k1 where k1 = dydx(n)) is provided as an input.
 * to evaluate dydx at different points, a call back function f (slope
 * function) is also passed as a parameter. Given values for y, and k1,
 * this function advances the solution over an interval h, and returns
 * the solution in yout. For details, see the discussion in "Numerical
 * Recipes in C", Chapter 16, from
 * http://www.nrbook.com/a/bookcpdf/c16-1.pdf
 */
void rk4_core(void *model, double *y, double *k1, void *p, int n, double h, double *yout, slope_fn_ptr f)
{
	int i;
	double *t, *k2, *k3, *k4;
	k2 = dvector(n);
	k3 = dvector(n);
	k4 = dvector(n);
	t = dvector(n);

	/* k2 is the slope at the trial midpoint (t) found using
	 * slope k1 (which is at the starting point).
	 */
	/* t = y + h/2 * k1 (t = y; t += h/2 * k1) */
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dcopy(n, y, 1, t, 1);
	cblas_daxpy(n, h/2.0, k1, 1, t, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	dcopy(n, y, 1, t, 1);
	daxpy(n, h/2.0, k1, 1, t, 1);
	#else
	for(i=0; i < n; i++)
		t[i] = y[i] + h/2.0 * k1[i];
	#endif
	/* k2 = slope at t */
	(*f)(model, t, p, k2);

	/* k3 is the slope at the trial midpoint (t) found using
	 * slope k2 found above.
	 */
	/* t =  y + h/2 * k2 (t = y; t += h/2 * k2) */
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dcopy(n, y, 1, t, 1);
	cblas_daxpy(n, h/2.0, k2, 1, t, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	dcopy(n, y, 1, t, 1);
	daxpy(n, h/2.0, k2, 1, t, 1);
	#else
	for(i=0; i < n; i++)
		t[i] = y[i] + h/2.0 * k2[i];
	#endif
	/* k3 = slope at t */
	(*f)(model, t, p, k3);

	/* k4 is the slope at trial endpoint (t) found using
	 * slope k3 found above.
	 */
	/* t =  y + h * k3 (t = y; t += h * k3) */
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dcopy(n, y, 1, t, 1);
	cblas_daxpy(n, h, k3, 1, t, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	dcopy(n, y, 1, t, 1);
	daxpy(n, h, k3, 1, t, 1);
	#else
	for(i=0; i < n; i++)
		t[i] = y[i] + h * k3[i];
	#endif
	/* k4 = slope at t */
	(*f)(model, t, p, k4);

	/* yout = y + h*(k1/6 + k2/3 + k3/3 + k4/6)	*/
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	/* yout = y	*/
	cblas_dcopy(n, y, 1, yout, 1);
	/* yout += h*k1/6	*/
	cblas_daxpy(n, h/6.0, k1, 1, yout, 1);
	/* yout += h*k2/3	*/
	cblas_daxpy(n, h/3.0, k2, 1, yout, 1);
	/* yout += h*k3/3	*/
	cblas_daxpy(n, h/3.0, k3, 1, yout, 1);
	/* yout += h*k4/6	*/
	cblas_daxpy(n, h/6.0, k4, 1, yout, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	dcopy(n, y, 1, yout, 1);
	/* yout += h*k1/6	*/
	daxpy(n, h/6.0, k1, 1, yout, 1);
	/* yout += h*k2/3	*/
	daxpy(n, h/3.0, k2, 1, yout, 1);
	/* yout += h*k3/3	*/
	daxpy(n, h/3.0, k3, 1, yout, 1);
	/* yout += h*k4/6	*/
	daxpy(n, h/6.0, k4, 1, yout, 1);
	#else
	for (i =0; i < n; i++)
		yout[i] = y[i] + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
	#endif

	free_dvector(k2);
	free_dvector(k3);
	free_dvector(k4);
	free_dvector(t);
}

/*
 * 4th order Runge Kutta solver	with adaptive step sizing.
 * It integrates and solves the ODE dy + cy = p between
 * t and t+h. It returns the correct step size to be used
 * next time. slope function f is the call back used to
 * evaluate the derivative at each point
 */
#define RK4_SAFETY		0.95
#define RK4_MAXUP		5.0
#define RK4_MAXDOWN		10.0
#define RK4_PRECISION	0.01
double rk4(void *model, double *y, void *p, int n, double *h, double *yout, slope_fn_ptr f)
{
	int i;
	double *k1, *t1, *t2, *ytemp, max, new_h = (*h);

	k1 = dvector(n);
	t1 = dvector(n);
	t2 = dvector(n);
	ytemp = dvector(n);

	/* evaluate the slope k1 at the beginning */
	(*f)(model, y, p, k1);

	/* try until accuracy is achieved	*/
  do {
		(*h) = new_h;

		/* try RK4 once with normal step size	*/
		rk4_core(model, y, k1, p, n, (*h), ytemp, f);

		/* repeat it with two half-steps	*/
		rk4_core(model, y, k1, p, n, (*h)/2.0, t1, f);

		/* y after 1st half-step is in t1. re-evaluate k1 for this	*/
		(*f)(model, t1, p, k1);

		/* get output of the second half-step in t2	*/
		rk4_core(model, t1, k1, p, n, (*h)/2.0, t2, f);

		/* find the max diff between these two results:
		 * use t1 to store the diff
		 */
		#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
		/* t1 = ytemp	*/
		cblas_dcopy(n, ytemp, 1, t1, 1);
		/* t1 = -1 * t2 + t1 = ytemp - t2	*/
		cblas_daxpy(n, -1.0, t2, 1, t1, 1);
		/* max = |t1[max_abs_index]|	*/
		max = fabs(t1[cblas_idamax(n, t1, 1)]);
		#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
		/* t1 = ytemp	*/
		dcopy(n, ytemp, 1, t1, 1);
		/* t1 = -1 * t2 + t1 = ytemp - t2	*/
		daxpy(n, -1.0, t2, 1, t1, 1);
		/*
		 * max = |t1[max_abs_index]| note: FORTRAN BLAS
		 * indices start from 1 as opposed to CBLAS where
		 * indices start from 0
		 */
		max = fabs(t1[idamax(n, t1, 1)-1]);
		#else
		for(i=0; i < n; i++) {
			t1[i] = fabs(ytemp[i] - t2[i]);
    }
		max = t1[0];
		for(i=1; i < n; i++)
			if (max < t1[i])
				max = t1[i];
		#endif

		/*
		 * compute the correct step size: see equation
		 * 16.2.10  in chapter 16 of "Numerical Recipes
		 * in C"
		 */
		/* accuracy OK. increase step size	*/
		if (max <= RK4_PRECISION) {
			new_h = RK4_SAFETY * (*h) * pow(fabs(RK4_PRECISION/max), 0.2);
			if (new_h > RK4_MAXUP * (*h))
				new_h = RK4_MAXUP * (*h);
		/* inaccuracy error. decrease step size	and compute again */
		} else {
			new_h = RK4_SAFETY * (*h) * pow(fabs(RK4_PRECISION/max), 0.25);
			if (new_h < (*h) / RK4_MAXDOWN)
				new_h = (*h) / RK4_MAXDOWN;
		}

	} while (new_h < (*h));

	/* commit ytemp to yout	*/
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dcopy(n, ytemp, 1, yout, 1);
	#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
	dcopy(n, ytemp, 1, yout, 1);
	#else
	copy_dvector(yout, ytemp, n);
	#endif

	/* clean up */
	free_dvector(k1);
	free_dvector(t1);
	free_dvector(t2);
	free_dvector(ytemp);

	/* return the step-size	*/
	return new_h;
}

/* matmult: C = AB, A, B are n x n square matrices	*/
void matmult(double **c, double **a, double **b, int n)
{
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				n, n, n, 1.0, a[0], n, b[0], n, 0.0, c[0], n);
	#elif (MATHACCEL == MA_AMD  || MATHACCEL == MA_SUN)
	/* B^T * A^T = (A * B)^T	*/
	dgemm('N', 'N', n, n, n, 1.0, b[0], n, a[0], n, 0.0, c[0], n);
	#else
	int i, j, k;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	#endif
}

/* same as above but 'a' is a diagonal matrix stored as a 1-d array	*/
void diagmatmult(double **c, double *a, double **b, int n)
{
	int i;
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	zero_dmatrix(c, n, n);
	for(i=0; i < n; i++)
		cblas_daxpy(n, a[i], b[i], 1, c[i], 1);
	#elif (MATHACCEL == MA_AMD  || MATHACCEL == MA_SUN)
	zero_dmatrix(c, n, n);
	for(i=0; i < n; i++)
		daxpy(n, a[i], b[i], 1, c[i], 1);
	#else
	int j;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			c[i][j] = a[i] * b[i][j];
	#endif
}

/* mult of an n x n matrix and an n x 1 column vector	*/
void matvectmult(double *vout, double **m, double *vin, int n)
{
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, m[0],
				n, vin, 1, 0.0, vout, 1);
	#elif (MATHACCEL == MA_AMD  || MATHACCEL == MA_SUN)
	dgemv('T', n, n, 1.0, m[0], n, vin, 1, 0.0, vout, 1);
	#else
	int i, j;

	for (i = 0; i < n; i++) {
		vout[i] = 0;
		for (j = 0; j < n; j++)
			vout[i] += m[i][j] * vin[j];
	}
	#endif
}

/* same as above but 'm' is a diagonal matrix stored as a 1-d array	*/
void diagmatvectmult(double *vout, double *m, double *vin, int n)
{
	#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
	cblas_dsbmv(CblasRowMajor, CblasUpper, n, 0, 1.0, m, 1, vin,
				1, 0.0, vout, 1);
	#elif (MATHACCEL == MA_AMD  || MATHACCEL == MA_SUN)
	dsbmv('U', n, 0, 1.0, m, 1, vin, 1, 0.0, vout, 1);
	#else
	int i;

	for (i = 0; i < n; i++)
		vout[i] = m[i] * vin[i];
	#endif
}

/*
 * inv = m^-1, inv, m are n by n matrices.
 * the spd flag indicates that m is symmetric
 * and positive definite
 */
#define BLOCK_SIZE		256
void matinv(double **inv, double **m, int n, int spd)
{
	int *p, lwork;
	double *work;
	int i, j;
	#if (MATHACCEL != MA_NONE)
	int info;
	#endif
	double *col;

	p = ivector(n);
	lwork = n * BLOCK_SIZE;
	work = dvector(lwork);

	#if (MATHACCEL == MA_INTEL)
	info = 0;
	cblas_dcopy(n*n, m[0], 1, inv[0], 1);
	if (!spd) {
		dgetrf(&n, &n, inv[0], &n, p, &info);
		assert(info == 0);
		dgetri(&n, inv[0], &n, p, work, &lwork, &info);
		assert(info == 0);
	} else {
		dpotrf("U", &n, inv[0], &n, &info);
		assert(info == 0);
		dpotri("U", &n, inv[0], &n, &info);
		assert(info == 0);
		mirror_dmatrix(inv, n);
	}
	#elif(MATHACCEL == MA_AMD)
	info = 0;
	dcopy(n*n, m[0], 1, inv[0], 1);
	if (!spd) {
		dgetrf_(&n, &n, inv[0], &n, p, &info);
		assert(info == 0);
		dgetri_(&n, inv[0], &n, p, work, &lwork, &info);
		assert(info == 0);
	} else {
		dpotrf_("U", &n, inv[0], &n, &info, 1);
		assert(info == 0);
		dpotri_("U", &n, inv[0], &n, &info, 1);
		assert(info == 0);
		mirror_dmatrix(inv, n);
	}
	#elif (MATHACCEL == MA_APPLE)
	info = 0;
	cblas_dcopy(n*n, m[0], 1, inv[0], 1);
	if (!spd) {
		dgetrf_((__CLPK_integer *)&n, (__CLPK_integer *)&n,
				inv[0], (__CLPK_integer *)&n, (__CLPK_integer *)p,
				(__CLPK_integer *)&info);
		assert(info == 0);
		dgetri_((__CLPK_integer *)&n, inv[0], (__CLPK_integer *)&n,
				(__CLPK_integer *)p, work, (__CLPK_integer *)&lwork,
				(__CLPK_integer *)&info);
		assert(info == 0);
	} else {
		dpotrf_("U", (__CLPK_integer *)&n, inv[0], (__CLPK_integer *)&n,
				(__CLPK_integer *)&info);
		assert(info == 0);
		dpotri_("U", (__CLPK_integer *)&n, inv[0], (__CLPK_integer *)&n,
				(__CLPK_integer *)&info);
		assert(info == 0);
		mirror_dmatrix(inv, n);
	}
	#elif(MATHACCEL == MA_SUN)
	info = 0;
	dcopy(n*n, m[0], 1, inv[0], 1);
	if (!spd) {
		dgetrf_(&n, &n, inv[0], &n, p, &info);
		assert(info == 0);
		dgetri_(&n, inv[0], &n, p, work, &lwork, &info);
		assert(info == 0);
	} else {
		dpotrf_("U", &n, inv[0], &n, &info);
		assert(info == 0);
		dpotri_("U", &n, inv[0], &n, &info);
		assert(info == 0);
		mirror_dmatrix(inv, n);
	}
	#else
	col = dvector(n);

	lupdcmp(m, n, p, spd);

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) col[i]=0.0;
		col[j] = 1.0;
		lusolve(m, n, p, col, work, spd);
		for (i = 0; i < n; i++) inv[i][j]=work[i];
	}

	free_dvector(col);
	#endif

	free_ivector(p);
	free_dvector(work);
}

/* dst = src1 + scale * src2	*/
void scaleadd_dvector (double *dst, double *src1, double *src2, int n, double scale)
{
	#if (MATHACCEL == MA_NONE)
	int i;
	for(i=0; i < n; i++)
		dst[i] = src1[i] + scale * src2[i];
	#else
	/* dst == src2. so, dst *= scale, dst += src1	*/
	if (dst == src2 && dst != src1) {
		#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
		cblas_dscal(n, scale, dst, 1);
		cblas_daxpy(n, 1.0, src1, 1, dst, 1);
		#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
		dscal(n, scale, dst, 1);
		daxpy(n, 1.0, src1, 1, dst, 1);
		#else
		fatal("unknown math acceleration\n");
		#endif
	/* dst = src1; dst += scale * src2	*/
	} else {
		#if (MATHACCEL == MA_INTEL || MATHACCEL == MA_APPLE)
		cblas_dcopy(n, src1, 1, dst, 1);
		cblas_daxpy(n, scale, src2, 1, dst, 1);
		#elif (MATHACCEL == MA_AMD || MATHACCEL == MA_SUN)
		dcopy(n, src1, 1, dst, 1);
		daxpy(n, scale, src2, 1, dst, 1);
		#else
		fatal("unknown math acceleration\n");
		#endif
	}
	#endif
}
