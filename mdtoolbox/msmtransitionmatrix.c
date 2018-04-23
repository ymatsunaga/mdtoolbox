#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define x(m,n) (X[m + n*nstate])
#define x_new(m,n) (X_NEW[m + n*nstate])
#define c(m,n) (C[m + n*nstate])
#define x_rs(m) (X_RS[m])
#define c_rs(m) (C_RS[m])

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* inputs */
  double  *C;
  double  *tol_pointer;
  double  tol;

  /* outputs */
  double  *T;
  double  *pi_i;

  /* working variables */
  int     nstate;
  int     iter;
  int     i, j, k, l;
  int     m, n;
  double  *X;
  double  *X_NEW;
  double  *X_RS;
  double  *C_RS;
  double  logl, oldlogl;
  double  a, b, c;
  double  v;
  double  tmp, denom, pi_sum;
  int     is_error;

  /* get inputs */
  if (mxIsSparse(prhs[0])) {
    mexErrMsgIdAndTxt("MDToolbox:msmtransitionmatrix", "Mex Error: sparse matrix is nor supported");
  }
  
  C = mxGetPr(prhs[0]);

  if ( nrhs < 2 ) {
    tol = 0.0001;
  } else {
    tol_pointer = mxGetPr(prhs[1]);
    tol = tol_pointer[0];
  }

  #ifdef DEBUG
    mexPrintf("MEX: C[0, 0] = %f\n", C[0]);
    mexPrintf("MEX: TOLERANCE = %f\n", tol);
  #endif

  /* setup: working variables */
  nstate = mxGetM(prhs[0]);
  X = (double*) malloc(nstate*nstate*sizeof(double));
  X_NEW = (double*) malloc(nstate*nstate*sizeof(double));
  X_RS = (double*) malloc(nstate*sizeof(double));
  C_RS = (double*) malloc(nstate*sizeof(double));

  #ifdef DEBUG
    mexPrintf("MEX: step 1\n");
  #endif

  /* allocate output variables */
  plhs[0] = mxCreateDoubleMatrix(nstate, nstate, mxREAL);
  T = mxGetPr(plhs[0]);
  for (m = 0; m < nstate; m++) {
    for (n = 0; n < nstate; n++) {
      T[m + n*nstate] = 0.0;
    }
  }

  #ifdef DEBUG
    mexPrintf("MEX: step 2\n");
  #endif

  plhs[1] = mxCreateDoubleMatrix(1, nstate, mxREAL);
  pi_i = mxGetPr(plhs[1]);
  for (m = 0; m < nstate; m++) {
    pi_i[m] = 0.0;
  }

  #ifdef DEBUG
    mexPrintf("MEX: step 3\n");
  #endif

  /* initialize X */
  for (i = 0; i < nstate; i++)
    for (j = 0; j < nstate; j++)
      x(i,j) = c(i,j) + c(j,i);

  #ifdef DEBUG
    mexPrintf("MEX: step 4\n");
  #endif

  /* initialize x_rs and c_rs */
  for (i = 0; i < nstate; i++) {
    x_rs(i) = 0;
    c_rs(i) = 0;
    for (j = 0; j < nstate; j++) {
      x_rs(i) += x(i,j);
      c_rs(i) += c(i,j);
    }

    if (x_rs(i) <= 0 || c_rs(i) <= 0) {
      mexErrMsgIdAndTxt("MDToolbox:msmtransitionmatrix", "Mex Error: domain error. we can't have rows with sum=0");
    }
  }

  #ifdef DEBUG
    mexPrintf("MEX: step 5\n");
  #endif

  oldlogl = 2*tol;
  logl = 0;
  is_error = 0;

  for (iter=0; fabs(oldlogl - logl) >= tol; iter++) {
    oldlogl = logl;
    logl = 0;

    #ifdef DEBUG
      mexPrintf("MEX: step 6\n");
    #endif

    /* fixed-point method */
    #pragma omp parallel for \
    default(none) \
    private(i, j, a, b, c, v, denom)                      \
    shared(nstate, X, X_NEW, C, X_RS, C_RS, is_error)
    for (i = 0; i < nstate; i++) {
      for (j = 0; j < nstate; j++) {
        denom = (c_rs(i)/x_rs(i)) + (c_rs(j)/x_rs(j));
        if (denom <= 0) is_error = 1;
        x_new(i, j) = (c(i, j) + c(j, i)) / denom;
      }
    }
    if (is_error)
      mexErrMsgIdAndTxt("MDToolbox:msmtransitionmatrix", "Mex Error: denominator becomes non-positive");

    #ifdef DEBUG
      mexPrintf("MEX: step 8\n");
    #endif

    /* update */
    #pragma omp parallel for \
    default(none) \
    reduction(+:logl) \
    private(i, j) \
    shared(nstate, X, X_NEW, C, X_RS, C_RS, is_error)
    for (i = 0; i < nstate; i++) {
      x_rs(i) = 0;
      for (j = 0; j < nstate; j++) {
        x_rs(i) += x_new(i,j);
      }
      if (x_rs(i) <= 0) {
        is_error = 1;
      }
      for (j = 0; j < nstate; j++) {
        x(i,j) = x_new(i,j);
      }
      for (j = 0; j < nstate; j++) {
        if (x(i,j) > 0)
          logl += c(i,j) * log(x(i,j) / x_rs(i));
      }
    }
    if (is_error)
      mexErrMsgIdAndTxt("MDToolbox:msmtransitionmatrix", "Mex Error: x_rs(i) <= 0. No transition counts");

    #ifdef DEBUG
      mexPrintf("MEX: step 9\n");
    #endif

    if ((iter+1)%10 == 0) 
      mexPrintf("MEX: iteration %d  LogLikelihood = %f  delta = %8.5e\n", iter+1,  logl, fabs(oldlogl - logl));
    if (logl != logl)
      mexErrMsgIdAndTxt("MDToolbox:msmtransitionmatrix", "Mex Error: LogLikelihood is something wrong");
  }

  pi_sum = 0;
  for (i = 0; i < nstate; i++) {
    pi_sum += x_rs(i);
    for (j = 0; j < nstate; j++) {
      T[i + j*nstate] = x(i,j) / x_rs(i);
    }
  }
  for (i = 0; i < nstate; i++)
    pi_i[i] = x_rs(i) / pi_sum;

  if (X != NULL) {
    free(X);
  }

  if (X_NEW != NULL) {
    free(X_NEW);
  }

  if (X_RS != NULL) {
    free(X_RS);
  }

  if (C_RS != NULL) {
    free(C_RS);
  }

  /* exit(EXIT_SUCCESS); */
}

