#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* inputs */
  double  *N_k_pr;
  int     *N_k;
  double  *f_k;
  double  *u_kln;
  double  *u_kn;
  double  *K_pr;
  int      K;
  double  *N_max_pr;
  int      N_max;

  /* outputs */
  double  *log_wi_jn;

  /* working variables */
  int     i, j;
  int     k, l;
  int     n;
  int     mrows;
  int     ncols;

  double  *FlogN;
  double  *log_term;
  double  term_sum;
  double  max_log_term;
  double  log_sum;
  double  u_l, u_k;

  /* check inputs and outputs */
  if ( nrhs < 6 ) {
    mexErrMsgTxt("MEX: Not enough input arguments.");
  }

  /* get inputs */
  N_k_pr   = mxGetPr(prhs[0]);
  f_k      = mxGetPr(prhs[1]);
  u_kln    = mxGetPr(prhs[2]);
  u_kn     = mxGetPr(prhs[3]);
  K_pr     = mxGetPr(prhs[4]);
  N_max_pr = mxGetPr(prhs[5]);

  K     = (int) (K_pr[0] + 0.5);
  N_max = (int) (N_max_pr[0] + 0.5);

  N_k = (int *) malloc(K*sizeof(int));
  for (k = 0; k < K; k++) {
    N_k[k] = (int) (N_k_pr[k] + 0.5);
  }

  #ifdef DEBUG
    mexPrintf("MEX: K = %d\n", K);
    mexPrintf("MEX: N_max = %d\n", N_max);
  #endif

  /* setup: working variables */
  FlogN = (double *) malloc(K*sizeof(double));
  log_term = (double *) malloc(K*sizeof(double));

  for (k = 0; k < K; k++) {
    FlogN[k] = log((double)N_k[k])+f_k[k];
  }

  /* allocate output variables */
  plhs[0] = mxCreateDoubleMatrix(K, N_max, mxREAL);
  log_wi_jn = mxGetPr(plhs[0]);

  /* calculation */
  for (k = 0; k < K; k++) {

    for (n = 0; n < N_k[k]; n++) {
      max_log_term = -1e100;
      u_k = u_kn[k + n*K];
      
      for (l = 0; l < K; l++) {
        u_l = u_kln[k + l*(K) + n*(K*K)];
        log_term[l] = FlogN[l] - (u_l - u_k);
        if (log_term[l] > max_log_term) {max_log_term = log_term[l];}
      }

      term_sum = 0.0;
      for (l = 0; l < K; l++) {
        term_sum += exp(log_term[l]-max_log_term);
      }
      log_sum = log(term_sum) + max_log_term;
      log_wi_jn[k + n*K] = -log_sum;
    }
  }

  if (N_k != NULL) {
    free(N_k);
  }

  if (FlogN != NULL) {
    free(FlogN);
  }

  if (log_term != NULL) {
    free(log_term);
  }

  /* exit(EXIT_SUCCESS); */
}

