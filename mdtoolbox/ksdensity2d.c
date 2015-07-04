#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* inputs */
  double  *data;
  double  *grid_x;
  double  *grid_y;
  double  *bandwidth;
  double  *box;
  double  *weight;

  /* outputs */
  double  *f;

  /* working variables */
  int     nstep;
  int     nx, ny;
  double  dx, dy;
  double  dgrid_x, dgrid_y;
  int     dims[2];
  int     i, j, istep;
  int     ix, iy;
  int     ix_min, ix_max;
  int     iy_min, iy_max;
  double  rx, ry;
  double  *gaussx, *gaussy;
  double  *f_private;
  int     is_box;
  int     alloc_bandwidth;
  int     alloc_weight;

  int     *ix_array;
  int     *iy_array;
  int     ix_count;
  int     iy_count;

  /* check inputs and outputs */
  if (nrhs < 4) {
    mexErrMsgTxt("MEX: Not enough input arguments. MEX version requires both grids and bandwidth.");
  }

  if (nlhs > 1) {
    mexErrMsgTxt("MEX: Too many output arguments.");
  }

  /* get inputs */
  data      = mxGetPr(prhs[0]);
  grid_x    = mxGetPr(prhs[1]);
  grid_y    = mxGetPr(prhs[2]);

  nstep = mxGetM(prhs[0]);
  nx    = mxGetNumberOfElements(prhs[1]);
  ny    = mxGetNumberOfElements(prhs[2]);
  dims[0] = nx;
  dims[1] = ny;

  #ifdef DEBUG
    mexPrintf("MEX: nstep = %d\n", nstep);
    mexPrintf("MEX: nx    = %d\n", nx);
    mexPrintf("MEX: ny    = %d\n", ny);
  #endif

  /* setup: bandwidth */
  dgrid_x = grid_x[1] - grid_x[0];
  dgrid_y = grid_y[1] - grid_y[0];

  /* setup: bandwidth */
  bandwidth = NULL;
  if (nrhs > 3) {
    if (mxGetNumberOfElements(prhs[3]) != 0) {
      alloc_bandwidth = 0; /* not allocated */
      bandwidth = mxGetPr(prhs[3]);
      rx = bandwidth[0]*5.0;
      ry = bandwidth[1]*5.0;
    }
  }
  if (bandwidth == NULL) {
    /* ................TODO................ */
    /* alloc_bandwidth = 1; */
    /* bandwidth = (double *) malloc(3*sizeof(double)); */
    /* mexErrMsgTxt("MEX: In MEX version, please specify bandwidth explicitly"); */
  }

  #ifdef DEBUG
    mexPrintf("MEX: bandwidth in x-axis = %f\n", bandwidth[0]);
    mexPrintf("MEX: bandwidth in y-axis = %f\n", bandwidth[1]);
  #endif

  /* setup: box */
  box = NULL;
  is_box = 0;
  if (nrhs > 4) {
    if (mxGetNumberOfElements(prhs[4]) != 0) {
      is_box = 1; /* box is turned on */
      box = mxGetPr(prhs[4]);
    }
  }

  /* setup: weight */
  weight = NULL;
  if (nrhs > 5) {
    if (mxGetNumberOfElements(prhs[5]) != 0) {
      alloc_weight = 0; /* not allocated */
      weight = mxGetPr(prhs[5]);
    }
  }
  if (weight == NULL) {
    alloc_weight = 1; /* allocated */
    weight = (double *) malloc(nstep*sizeof(double));
    for (istep = 0; istep < nstep; istep++) {
      weight[istep] = 1.0/(double)nstep;
    }
  }

  /* setup: f */
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  f = mxGetPr(plhs[0]);

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      f[iy*nx + ix] = 0.0;
    }
  }

  /* calculation */
  if (is_box) {

    #pragma omp parallel \
    default(none) \
      private(istep, ix, ix_min, ix_max, iy, iy_min, iy_max, dx, dy, gaussx, gaussy, f_private, ix_count, iy_count, ix_array, iy_array) \
      shared(nstep, grid_x, grid_y, dgrid_x, dgrid_y, rx, ry, data, bandwidth, weight, f, nx, ny, alloc_bandwidth, alloc_weight, box)
    {
      f_private = (double *) malloc(nx*ny*sizeof(double));
      gaussx    = (double *) malloc(nx*sizeof(double));
      gaussy    = (double *) malloc(ny*sizeof(double));
      ix_array  = (int *) malloc(nx*sizeof(int));
      iy_array  = (int *) malloc(ny*sizeof(int));

      for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
          f_private[iy*nx + ix] = 0.0;
        }
      }

      for (ix = 0; ix < nx; ix++) {
        gaussx[ix] = 0.0;
      }
      for (iy = 0; iy < ny; iy++) {
        gaussy[iy] = 0.0;
      }

      #pragma omp for
      for (istep = 0; istep < nstep; istep++) {

        ix_count = 0;
        for (ix = 0; ix < nx; ix++) {
          dx = data[istep + nstep*0] - grid_x[ix];
          dx = dx - round(dx/box[0])*box[0];
          if (fabs(dx) < rx) {
            dx = dx/bandwidth[0];
            gaussx[ix_count] = exp(-0.5*dx*dx)/(sqrt(2*M_PI)*bandwidth[0]);
            ix_array[ix_count] = ix;
            ix_count++;
          }
        }

        iy_count = 0;
        for (iy = 0; iy < ny; iy++) {
          dy = data[istep + nstep*1] - grid_y[iy];
          dy = dy - round(dy/box[1])*box[1];
          if (fabs(dy) < ry) {
            dy = dy/bandwidth[1];
            gaussy[iy_count] = exp(-0.5*dy*dy)/(sqrt(2*M_PI)*bandwidth[1]);
            iy_array[iy_count] = iy;
            iy_count++;
          }
        }

        for (ix = 0; ix < ix_count; ix++) {
          for (iy = 0; iy < iy_count; iy++) {
            f_private[iy_array[iy]*nx + ix_array[ix]] += weight[istep]*gaussx[ix]*gaussy[iy];
          }
        }

        for (ix = 0; ix < nx; ix++) {
          gaussx[ix] = 0.0;
        }
        for (iy = 0; iy < ny; iy++) {
          gaussy[iy] = 0.0;
        }

      } /* pragma omp for */

      #pragma omp critical
      for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
          f[iy*nx + ix] += f_private[iy*nx + ix];
        }
      }

      if (f_private != NULL) {
        free(f_private);
      }
      if (gaussx != NULL) {
        free(gaussx);
      }
      if (gaussy != NULL) {
        free(gaussy);
      }
      if (ix_array != NULL) {
        free(ix_array);
      }
      if (iy_array != NULL) {
        free(iy_array);
      }

    } /* pragma omp parallel */

  } else {
    
    #pragma omp parallel \
    default(none) \
      private(istep, ix, ix_min, ix_max, iy, iy_min, iy_max, dx, dy, gaussx, gaussy, f_private) \
      shared(nstep, grid_x, grid_y, dgrid_x, dgrid_y, rx, ry, data, bandwidth, weight, f, nx, ny, alloc_bandwidth, alloc_weight)
    {
      f_private = (double *) malloc(nx*ny*sizeof(double));
      gaussx    = (double *) malloc(nx*sizeof(double));
      gaussy    = (double *) malloc(ny*sizeof(double));

      for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
          f_private[iy*nx + ix] = 0.0;
        }
      }

      for (ix = 0; ix < nx; ix++) {
        gaussx[ix] = 0.0;
      }
      for (iy = 0; iy < ny; iy++) {
        gaussy[iy] = 0.0;
      }

      #pragma omp for
      for (istep = 0; istep < nstep; istep++) {

        dx = data[istep + nstep*0] - grid_x[0];
        ix_min = (int) ((dx - rx)/dgrid_x);
        ix_min = ix_min > 0 ? ix_min : 0;
        ix_max = ((int) ((dx + rx)/dgrid_x)) + 1;
        ix_max = ix_max < nx ? ix_max : nx;
        for (ix = ix_min; ix < ix_max; ix++) {
          dx = (grid_x[ix] - data[istep + nstep*0])/bandwidth[0];
          gaussx[ix] = exp(-0.5*dx*dx)/(sqrt(2*M_PI)*bandwidth[0]);
        }

        dy = data[istep + nstep*1] - grid_y[0];
        iy_min = (int) ((dy - ry)/dgrid_y);
        iy_min = iy_min > 0 ? iy_min : 0;
        iy_max = ((int) ((dy + ry)/dgrid_y)) + 1;
        iy_max = iy_max < ny ? iy_max : ny;
        for (iy = iy_min; iy < iy_max; iy++) {
          dy = (grid_y[iy] - data[istep + nstep*1])/bandwidth[1];
          gaussy[iy] = exp(-0.5*dy*dy)/(sqrt(2*M_PI)*bandwidth[1]);
        }

        for (ix = ix_min; ix < ix_max; ix++) {
          for (iy = iy_min; iy < iy_max; iy++) {
            f_private[iy*nx + ix] += weight[istep]*gaussx[ix]*gaussy[iy];
          }
        }

        for (ix = ix_min; ix < ix_max; ix++) {
          gaussx[ix] = 0.0;
        }
        for (iy = iy_min; iy < iy_max; iy++) {
          gaussy[iy] = 0.0;
        }

      } /* pragma omp for */

      #pragma omp critical
      for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
          f[iy*nx + ix] += f_private[iy*nx + ix];
        }
      }

      if (f_private != NULL) {
        free(f_private);
      }
      if (gaussx != NULL) {
        free(gaussx);
      }
      if (gaussy != NULL) {
        free(gaussy);
      }

    } /* pragma omp parallel */

  } /* is_box */

  if (alloc_bandwidth) {
    free(bandwidth);
  }
  if (alloc_weight) {
    free(weight);
  }

  /* exit(EXIT_SUCCESS); */
}

