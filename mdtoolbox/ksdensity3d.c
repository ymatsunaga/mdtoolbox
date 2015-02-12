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
  double  *grid_z;
  double  *bandwidth;
  double  *weight;

  /* outputs */
  double  *f;

  /* working variables */
  int     nstep;
  int     nx, ny, nz;
  double  dx, dy, dz;
  int     dims[3];
  int     i, j, istep;
  int     ix, iy, iz;
  int     ix_min, ix_max;
  int     iy_min, iy_max;
  int     iz_min, iz_max;
  double  *gaussx, *gaussy, *gaussz;
  double  *f_private;
  int     alloc_bandwidth;
  int     alloc_weight;

  /* check inputs and outputs */
  if (nrhs < 5) {
    mexErrMsgTxt("MEX: Not enough input arguments. MEX version requires both grids and bandwidth.");
  }

  if (nlhs > 1) {
    mexErrMsgTxt("MEX: Too many output arguments.");
  }

  /* get inputs */
  data      = mxGetPr(prhs[0]);
  grid_x    = mxGetPr(prhs[1]);
  grid_y    = mxGetPr(prhs[2]);
  grid_z    = mxGetPr(prhs[3]);
  bandwidth = mxGetPr(prhs[4]);

  nstep = mxGetM(prhs[0]);
  nx    = mxGetNumberOfElements(prhs[1]);
  ny    = mxGetNumberOfElements(prhs[2]);
  nz    = mxGetNumberOfElements(prhs[3]);
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = nz;

  #ifdef DEBUG
    mexPrintf("MEX: nstep = %d\n", nstep);
    mexPrintf("MEX: nx    = %d\n", nx);
    mexPrintf("MEX: ny    = %d\n", ny);
    mexPrintf("MEX: nz    = %d\n", nz);
  #endif

  /* setup: bandwidth */
  bandwidth = NULL;
  if (nrhs > 4) {
    if (mxGetNumberOfElements(prhs[4]) != 0) {
      alloc_bandwidth = 0; /* not allocated */
      bandwidth = mxGetPr(prhs[4]);
    }
  }
  if (bandwidth == NULL) {
    /* ................TODO................ */
    /* alloc_bandwidth = 1; */
    /* bandwidth = (double *) malloc(3*sizeof(double)); */
    /* mexErrMsgTxt("MEX: In MEX version, please specify bandwidth explicitly"); */
  }
  mexPrintf("MEX: bandwidth in x-axis: %f\n", bandwidth[0]);
  mexPrintf("MEX: bandwidth in y-axis: %f\n", bandwidth[1]);
  mexPrintf("MEX: bandwidth in z-axis: %f\n", bandwidth[2]);

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
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  f = mxGetPr(plhs[0]);

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      for (iz = 0; iz < nz; iz++) {
        f[iz*ny*nx + iy*nx + ix] = 0.0;
      }
    }
  }

  /* calculation */
  #pragma omp parallel \
  default(none) \
    private(istep, ix, ix_min, ix_max, iy, iy_min, iy_max, iz, iz_min, iz_max, dx, dy, dz, gaussx, gaussy, gaussz, f_private) \
    shared(nstep, grid_x, grid_y, grid_z, data, bandwidth, weight, f)
  {
    f_private = (double *) malloc(nx*ny*nz*sizeof(double));
    gaussx    = (double *) malloc(nx*sizeof(double));
    gaussy    = (double *) malloc(ny*sizeof(double));
    gaussz    = (double *) malloc(nz*sizeof(double));

    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {
          f_private[iz*ny*nx + iy*nx + ix] = 0.0;
        }
      }
    }

    for (ix = 0; ix < nx; ix++) {
      gaussx[ix] = 0.0;
    }
    for (iy = 0; iy < ny; iy++) {
      gaussy[iy] = 0.0;
    }
    for (iz = 0; iz < nz; iz++) {
      gaussz[iz] = 0.0;
    }

    #pragma omp for
    for (istep = 0; istep < nstep; istep++) {

      for (ix = 0; ix < nx; ix++) {
        dx = (grid_x[ix] - data[istep + nstep*0])/bandwidth[0];
        gaussx[ix] = exp(-0.5*dx*dx)/(sqrt(2*M_PI)*bandwidth[0]);
      }
      for (iy = 0; iy < ny; iy++) {
        dy = (grid_y[iy] - data[istep + nstep*1])/bandwidth[1];
        gaussy[iy] = exp(-0.5*dy*dy)/(sqrt(2*M_PI)*bandwidth[1]);
      }
      for (iz = 0; iz < nz; iz++) {
        dz = (grid_z[iz] - data[istep + nstep*2])/bandwidth[2];
        gaussz[iz] = exp(-0.5*dz*dz)/(sqrt(2*M_PI)*bandwidth[2]);
      }

      for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
          for (iz = 0; iz < nz; iz++) {
            f_private[iz*ny*nx + iy*nx + ix] += weight[istep]*gaussx[ix]*gaussy[iy]*gaussz[iz];
          }
        }
      }

      for (ix = 0; ix < nx; ix++) {
        gaussx[ix] = 0.0;
      }
      for (iy = 0; iy < ny; iy++) {
        gaussy[iy] = 0.0;
      }
      for (iz = 0; iz < nz; iz++) {
        gaussz[iz] = 0.0;
      }

    } /* pragma omp for */

    #pragma omp critical
    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {
          f[iz*ny*nx + iy*nx + ix] += f_private[iz*ny*nx + iy*nx + ix];
        }
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
    if (gaussz != NULL) {
      free(gaussz);
    }
    if (alloc_bandwidth) {
      free(bandwidth);
    }
    if (alloc_weight) {
      free(weight);
    }

  } /* omp pragma parallel */

  /* exit(EXIT_SUCCESS); */
}

