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
  double  *box;
  double  *weight;

  /* outputs */
  double  *f;

  /* working variables */
  int     nstep;
  int     nx, ny, nz;
  double  dx, dy, dz;
  double  dgrid_x, dgrid_y, dgrid_z;
  int     dims[3];
  int     i, j, istep;
  int     ix, iy, iz;
  int     ix_min, ix_max;
  int     iy_min, iy_max;
  int     iz_min, iz_max;
  double  rx, ry, rz;
  double  *gaussx, *gaussy, *gaussz;
  double  *f_private;
  int     is_box;
  int     alloc_bandwidth;
  int     alloc_weight;

  int     *ix_array;
  int     *iy_array;
  int     *iz_array;
  int     ix_count;
  int     iy_count;
  int     iz_count;

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
  dgrid_x = grid_x[1] - grid_x[0];
  dgrid_y = grid_y[1] - grid_y[0];
  dgrid_z = grid_z[1] - grid_z[0];

  /* setup: bandwidth */
  bandwidth = NULL;
  if (nrhs > 4) {
    if (mxGetNumberOfElements(prhs[4]) != 0) {
      alloc_bandwidth = 0; /* not allocated */
      bandwidth = mxGetPr(prhs[4]);
      rx = bandwidth[0]*5.0;
      ry = bandwidth[1]*5.0;
      rz = bandwidth[2]*5.0;
    }
  }
  if (bandwidth == NULL) {
    /* ................TODO................ */
    /* alloc_bandwidth = 1; */
    /* bandwidth = (double *) malloc(3*sizeof(double)); */
    /* mexErrMsgTxt("MEX: In MEX version, please specify bandwidth explicitly"); */
  }

  #ifdef DEBUG
    mexPrintf("MEX: bandwidth in x-axis: %f\n", bandwidth[0]);
    mexPrintf("MEX: bandwidth in y-axis: %f\n", bandwidth[1]);
    mexPrintf("MEX: bandwidth in z-axis: %f\n", bandwidth[2]);
  #endif

  /* setup: box */
  box = NULL;
  is_box = 0;
  if (nrhs > 5) {
    if (mxGetNumberOfElements(prhs[5]) != 0) {
      is_box = 1; /* box is turned on */
      box = mxGetPr(prhs[5]);
    }
  }

  /* setup: weight */
  weight = NULL;
  if (nrhs > 6) {
    if (mxGetNumberOfElements(prhs[6]) != 0) {
      alloc_weight = 0; /* not allocated */
      weight = mxGetPr(prhs[6]);
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
  if (is_box) {

    #pragma omp parallel \
    default(none) \
      private(istep, ix, ix_min, ix_max, iy, iy_min, iy_max, iz, iz_min, iz_max, dx, dy, dz, gaussx, gaussy, gaussz, f_private, ix_count, iy_count, iz_count, ix_array, iy_array, iz_array) \
      shared(nstep, grid_x, grid_y, grid_z, dgrid_x, dgrid_y, dgrid_z, rx, ry, rz, data, bandwidth, weight, f, nx, ny, nz, alloc_bandwidth, alloc_weight, box)
    {
      f_private = (double *) malloc(nx*ny*nz*sizeof(double));
      gaussx    = (double *) malloc(nx*sizeof(double));
      gaussy    = (double *) malloc(ny*sizeof(double));
      gaussz    = (double *) malloc(nz*sizeof(double));
      ix_array  = (int *) malloc(nx*sizeof(int));
      iy_array  = (int *) malloc(ny*sizeof(int));
      iz_array  = (int *) malloc(nz*sizeof(int));

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

        ix_count = 0;
        for (ix = 0; ix < nx; ix++) {
          dx = data[istep + nstep*0] - grid_x[ix];
          dx = dx - round(dx/box[0])*box[0];
          if (abs(dx) < rx) {
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
          if (abs(dy) < ry) {
            dy = dy/bandwidth[1];
            gaussy[iy_count] = exp(-0.5*dy*dy)/(sqrt(2*M_PI)*bandwidth[1]);
            iy_array[iy_count] = iy;
            iy_count++;
          }
        }

        iz_count = 0;
        for (iz = 0; iz < nz; iz++) {
          dz = data[istep + nstep*2] - grid_z[iz];
          dz = dz - round(dz/box[2])*box[2];
          if (abs(dz) < rz) {
            dz = dz/bandwidth[2];
            gaussz[iz_count] = exp(-0.5*dz*dz)/(sqrt(2*M_PI)*bandwidth[2]);
            iz_array[iz_count] = iz;
            iz_count++;
          }
        }

        for (ix = 0; ix < ix_count; ix++) {
          for (iy = 0; iy < iy_count; iy++) {
            for (iz = 0; iz < iz_count; iz++) {
              f_private[iz_array[iz]*ny*nx + iy_array[iy]*nx + ix_array[ix]] += weight[istep]*gaussx[ix]*gaussy[iy]*gaussz[iz];
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
      if (ix_array != NULL) {
        free(ix_array);
      }
      if (iy_array != NULL) {
        free(iy_array);
      }
      if (iz_array != NULL) {
        free(iz_array);
      }

    } /* pragma omp parallel */

  } else {

    #pragma omp parallel \
    default(none) \
      private(istep, ix, ix_min, ix_max, iy, iy_min, iy_max, iz, iz_min, iz_max, dx, dy, dz, gaussx, gaussy, gaussz, f_private) \
      shared(nstep, grid_x, grid_y, grid_z, dgrid_x, dgrid_y, dgrid_z, rx, ry, rz, data, bandwidth, weight, f, nx, ny, nz, alloc_bandwidth, alloc_weight)
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

        dx = data[istep + nstep*0] - grid_x[0];
        ix_min = (int) ((dx - rx)/dgrid_x);
        ix_min = ix_min > 0 ? ix_min : 0;
        ix_max = ((int) ((dx + rx)/dgrid_x)) + 1;
        ix_max = ix_max < nx ? ix_max : nx;
        /* mexPrintf("MEX: ix_min = %d\n", ix_min); */
        /* mexPrintf("MEX: ix_max = %d\n", ix_max); */
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

        dz = data[istep + nstep*2] - grid_z[0];
        iz_min = (int) ((dz - rz)/dgrid_z);
        iz_min = iz_min > 0 ? iz_min : 0;
        iz_max = ((int) ((dz + rz)/dgrid_z)) + 1;
        iz_max = iz_max < nz ? iz_max : nz;
        for (iz = iz_min; iz < iz_max; iz++) {
          dz = (grid_z[iz] - data[istep + nstep*2])/bandwidth[2];
          gaussz[iz] = exp(-0.5*dz*dz)/(sqrt(2*M_PI)*bandwidth[2]);
        }

        for (ix = ix_min; ix < ix_max; ix++) {
          for (iy = iy_min; iy < iy_max; iy++) {
            for (iz = iz_min; iz < iz_max; iz++) {
              f_private[iz*ny*nx + iy*nx + ix] += weight[istep]*gaussx[ix]*gaussy[iy]*gaussz[iz];
            }
          }
        }

        for (ix = ix_min; ix < ix_max; ix++) {
          gaussx[ix] = 0.0;
        }
        for (iy = iy_min; iy < iy_max; iy++) {
          gaussy[iy] = 0.0;
        }
        for (iz = iz_min; iz < iz_max; iz++) {
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

