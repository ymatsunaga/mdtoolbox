/*******************************************************************************
 *  -/_|:|_|_\- 
 *
 *  File:           qcprot.c
 *  Version:        1.4
 *
 *  Function:       Rapid calculation of the least-squares rotation using a 
 *                  quaternion-based characteristic polynomial and 
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *                  
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 * 
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for macromolecular 
 *      superpositions."
 *      in press, Journal of Computational Chemistry 
 *
 *
 *  Copyright (c) 2009-2012 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted 
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of 
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list 
 *    of conditions and the following disclaimer in the documentation and/or other materials 
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to 
 *    endorse or promote products derived from this software without specific prior written 
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *                    If trying all rows of the adjoint still gives too small
 *                    qsqr, then just return identity matrix. (DLT)
 *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
 *                    invalid mem access
 *    2011/02/21      Made CenterCoords use weights
 *    2011/05/02      Finally changed CenterCoords declaration in qcprot.h
 *                    Also changed some functions to static
 *    2011/07/08      put in fabs() to fix taking sqrt of small neg numbers, fp error
 *    2012/07/26      minor changes to comments and main.c, more info (v.1.4)
 *  
 ******************************************************************************/

/* #include "qcprot.h" */
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


static double
InnerProduct(double *A, double **coords1, double **coords2, const int len, const double *weight)
{
    double          x1, x2, y1, y2, z1, z2;
    int             i;
    const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
    const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
    double          G1 = 0.0, G2 = 0.0;

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

    if (weight != NULL)
    {
        for (i = 0; i < len; ++i)
        {
            x1 = weight[i] * fx1[i];
            y1 = weight[i] * fy1[i];
            z1 = weight[i] * fz1[i];

            G1 += x1 * fx1[i] + y1 * fy1[i] + z1 * fz1[i];

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);   
        }
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            x1 = fx1[i];
            y1 = fy1[i];
            z1 = fz1[i];

            G1 += x1 * x1 + y1 * y1 + z1 * z1;

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);  
        }
    }

    return (G1 + G2) * 0.5;
}


int
FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double wsum, double minScore)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double C[4];
    int i;
    double mxEigenV; 
    double oldg = 0.0;
    double b, a, delta, rms, qsqr;
    double q1, q2, q3, q4, normq;
    double a11, a12, a13, a14, a21, a22, a23, a24;
    double a31, a32, a33, a34, a41, a42, a43, a44;
    double a2, x2, y2, z2; 
    double xy, az, zx, ay, yz, ax; 
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
    double evecprec = 1e-6;
    double evalprec = 1e-11;

    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    mxEigenV = E0;
    for (i = 0; i < 50; ++i)
    {
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        /* printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV); */
        if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
            break;
    }

    if (i == 50) 
       fprintf(stderr,"\nMore than %d iterations needed!\n", i);

    /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
    /* rms = sqrt(fabs(2.0 * (E0 - mxEigenV)/len)); */
    rms = sqrt(fabs(2.0 * (E0 - mxEigenV)/wsum));
    (*rmsd) = rms;
    /* printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); */

    if (minScore > 0) 
        if (rms < minScore)
            return (-1);

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

/* The following code tries to calculate another column in the adjoint matrix when the norm of the 
   current column is too small.
   Usually this commented block will never be activated.  To be absolutely safe this should be
   uncommented, but it is most likely unnecessary.  
*/
    if (qsqr < evecprec)
    {
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

        if (qsqr < evecprec)
        {
            double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
            double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
            double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr < evecprec)
                {
                    /* if qsqr is still too small, return the identity matrix. */
                    rot[0] = rot[4] = rot[8] = 1.0;
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

                    return(0);
                }
            }
        }
    }

    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;

    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;

    rot[0] = a2 + x2 - y2 - z2;
    rot[1] = 2 * (xy + az);
    rot[2] = 2 * (zx - ay);
    rot[3] = 2 * (xy - az);
    rot[4] = a2 - x2 + y2 - z2;
    rot[5] = 2 * (yz + ax);
    rot[6] = 2 * (zx + ay);
    rot[7] = 2 * (yz - ax);
    rot[8] = a2 - x2 - y2 + z2;

    return (1);
}


static void
CenterCoords(double **coords, double *xsum, double *ysum, double *zsum, const int len, const double *weight)
{
    int            i;
    double         wsum;
    double         *x = coords[0], *y = coords[1], *z = coords[2];

    *xsum = *ysum = *zsum = 0.0;

    if (weight != NULL)
    {
        wsum = 0.0;
        for (i = 0; i < len; ++i)
        {
            *xsum += weight[i] * x[i];
            *ysum += weight[i] * y[i];
            *zsum += weight[i] * z[i];
            wsum += weight[i];
        }
        *xsum /= wsum;
        *ysum /= wsum;
        *zsum /= wsum;
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            *xsum += x[i];
            *ysum += y[i];
            *zsum += z[i];
        }
        *xsum /= len;
        *ysum /= len;
        *zsum /= len;
    }

    for (i = 0; i < len; ++i)
    {
        x[i] -= *xsum;
        y[i] -= *ysum;
        z[i] -= *zsum;
    }
}


/* Superposition coords2 onto coords1 -- in other words, coords2 is rotated, coords1 is held fixed */
double
CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight)
{
    int i;
    double wsum, xsum, ysum, zsum;
    double A[9], rmsd;

    /* center the structures -- if precentered you can omit this step */
    /* CenterCoords(coords1, &xsum, &ysum, &zsum, len, weight); */
    /* CenterCoords(coords2, &xsum, &ysum, &zsum, len, weight); */

    /* calculate the (weighted) inner product of two structures */
    double E0 = InnerProduct(A, coords1, coords2, len, weight);

    /* calculate the RMSD & rotational matrix */
    if (weight != NULL) {
      wsum = 0.0;
      for (i =0; i < len; i++) {
        wsum += weight[i];
      }
    } else {
      wsum = (double) len;
    }
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, len, wsum, -1);

    return rmsd;
}


double **MatInit(const int rows, const int cols)
{
    int             i;
    double        **matrix = NULL;
    double         *matspace = NULL;

    matspace = (double *) calloc((rows * cols), sizeof(double));
    if (matspace == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate matrix space in MatInit(): (%d x %d)\n", rows, cols);
        exit(EXIT_FAILURE);
    }

    /* allocate room for the pointers to the rows */
    matrix = (double **) malloc(rows * sizeof(double *));
    if (matrix == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate room for row pointers in MatInit(): (%d)\n", rows);
        exit(EXIT_FAILURE);
    }

    /*  now 'point' the pointers */
    for (i = 0; i < rows; i++)
        matrix[i] = matspace + (i * cols);

    return(matrix);
}


void MatDestroy(double ***matrix_ptr)
{
    double **matrix = *matrix_ptr;

    if (matrix != NULL)
    {
        if (matrix[0] != NULL)
        {
            free(matrix[0]);
            matrix[0] = NULL;
        }

        free(matrix);
        *matrix_ptr = NULL;
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* inputs */
  double  *ref;
  double  *trj;
  bool    *bindex;
  double  *dindex;
  double  *weight;
  double  *vel;

  /* outputs */
  double  *rmsd;
  double  *trj_fitted;
  double  *vel_fitted;

  /* working variables */
  int     i, j, istep;
  int     mrows;
  int     ncols;

  double  **frag_a, **frag_b;
  double  **frag_a_sub, **frag_b_sub;
  int     *index;
  double  *weight_sub;
  double  **frag_v;
  double  rotmat[9];

  mxArray *mxIsLogical;
  bool    *isLogical;

  int     natom, natom3;
  int     nstep;
  int     len;
  double  xsum_a, ysum_a, zsum_a;
  double  xsum_b, ysum_b, zsum_b;
  double  x, y, z;

  /* check inputs and outputs */
  if ( nrhs < 2 ) {
    mexErrMsgTxt("MEX: Not enough input arguments.");
  }

  if ( nlhs > 3 ) {
    mexErrMsgTxt("MEX: Too many output arguments.");
  }

  /* get inputs: array sizes */
  ncols = mxGetN(prhs[0]);
  natom3 = ncols;
  natom = ncols/3;
  mrows = mxGetM(prhs[1]);
  nstep = mrows;

  #ifdef DEBUG
    mexPrintf("MEX: *** stage 1 ***\n");
    mexPrintf("MEX: natom3 = %d\n", natom3);
    mexPrintf("MEX: natom = %d\n", natom);
    mexPrintf("MEX: nstep = %d\n", nstep);
    mexPrintf("\n", len);
  #endif

  /* setup: two coordinates */
  frag_a = NULL;
  frag_b = NULL;
  ref = mxGetPr(prhs[0]);
  frag_a = MatInit(3, natom);
  for ( i = 0; i < natom; i++) {
    for ( j = 0; j < 3; j++) {
      frag_a[j][i] = ref[3*i+j];
    }
  }
  trj = mxGetPr(prhs[1]);
  frag_b = MatInit(3, natom);

  /* setup: index */
  len = natom;
  index = NULL;
  if ( nrhs > 2 ) {
    if ( mxGetNumberOfElements(prhs[2]) != 0 ) {
      mexCallMATLAB(1, &mxIsLogical, 1, &prhs[2], "islogical");
      isLogical = (bool *) mxGetData(mxIsLogical);

      if( isLogical[0] ){
        /* logical index */
        bindex = (bool*) mxGetData(prhs[2]);
        len = 0;
        for (i = 0; i < mxGetNumberOfElements(prhs[2]); i++) {
          if ( bindex[i] ) {
            len++;
          }
        }
        index = (int *) malloc(len*sizeof(int));
        j = 0;
        for (i = 0; i < mxGetNumberOfElements(prhs[2]); i++) {
          if ( bindex[i] ) {
            index[j] = i;
            j++;
          }
        }
      } else {
        /* integer index */
        dindex = (double *) mxGetData(prhs[2]);
        len = mxGetNumberOfElements(prhs[2]);
        index = (int *) malloc(len*sizeof(int));
        for (i = 0; i < len; i++) {
          index[i] = (int) (dindex[i] - 0.5);
        }
      }

    }
  }

  #ifdef DEBUG
    mexPrintf("MEX: *** stage 2 ***\n");
    mexPrintf("MEX: natom = %d\n", natom);
    mexPrintf("MEX: len = %d\n", len);
    mexPrintf("\n", len);
  #endif

  /* setup: mass */
  weight = NULL;
  if ( nrhs > 3 ) {
    if ( mxGetNumberOfElements(prhs[3]) != 0 ) {
      weight = mxGetPr(prhs[3]);
    }
  }

  /* setup: velocity */
  frag_v = NULL;
  if ( nrhs > 4 ) {
    if ( mxGetNumberOfElements(prhs[4]) != 0 ) {
      vel = mxGetPr(prhs[4]);
      frag_v = MatInit(3, natom);
    }
  }

  /* setup: allocate subsets if index is given */
  frag_a_sub = NULL;
  frag_b_sub = NULL;
  weight_sub = NULL;
  if ( index != NULL ) {
    frag_a_sub = MatInit(3, len);
    for ( i = 0; i < len; i++) {
      for ( j = 0; j < 3; j++) {
        frag_a_sub[j][i] = frag_a[j][index[i]];
      }
    }
    frag_b_sub = MatInit(3, len);
    if ( weight != NULL ) {
      weight_sub = (double *) malloc(len*sizeof(double));
      for (i = 0; i < len; i++) {
        weight_sub[i] = weight[index[i]];
      }
    }
  }

  /* allocate output variables */
  plhs[0] = mxCreateDoubleMatrix(nstep, 1, mxREAL);
  rmsd = mxGetPr(plhs[0]);

  if ( nlhs > 1 ) {
    plhs[1] = mxCreateDoubleMatrix(nstep, natom3, mxREAL);
    trj_fitted = mxGetPr(plhs[1]);
  }

  if ( nlhs > 2 ) {
    plhs[2] = mxCreateDoubleMatrix(nstep, natom3, mxREAL);
    vel_fitted = mxGetPr(plhs[2]);
  }

  /* center coordinates */
  if( index == NULL ) {
    CenterCoords(frag_a, &xsum_a, &ysum_a, &zsum_a, natom, weight);
  } else {
    CenterCoords(frag_a_sub, &xsum_a, &ysum_a, &zsum_a, len, weight_sub);
  }

  /* calculate RMSD and Rotaion matrix */
  if( index == NULL ) {
    for (istep = 0; istep < nstep; istep++) {

      for (i = 0; i < natom; i++) {
        for (j = 0; j < 3; j++) {
          frag_b[j][i] = trj[istep + 3*nstep*i + nstep*j];
        }
      }

      CenterCoords(frag_b, &xsum_b, &ysum_b, &zsum_b, natom, weight);
      rmsd[istep] = CalcRMSDRotationalMatrix((double **) frag_a, (double **) frag_b, natom, rotmat, weight);

      if ( nlhs > 1 ) {
        /* apply rotation matrix */
        for (i = 0; i < natom; ++i) {
          x = rotmat[0]*frag_b[0][i] + rotmat[1]*frag_b[1][i] + rotmat[2]*frag_b[2][i];
          y = rotmat[3]*frag_b[0][i] + rotmat[4]*frag_b[1][i] + rotmat[5]*frag_b[2][i];
          z = rotmat[6]*frag_b[0][i] + rotmat[7]*frag_b[1][i] + rotmat[8]*frag_b[2][i];
        
          frag_b[0][i] = x + xsum_a;
          frag_b[1][i] = y + ysum_a;
          frag_b[2][i] = z + zsum_a;
        }

        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            trj_fitted[istep + 3*nstep*i + nstep*j] = frag_b[j][i];
          }
        }
      }

      if ( nlhs > 2 ) {
        /* apply rotation matrix to velocity*/
        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            frag_v[j][i] = vel[istep + 3*nstep*i + nstep*j];
          }
        }

        for (i = 0; i < natom; ++i) {
          x = rotmat[0]*frag_v[0][i] + rotmat[1]*frag_v[1][i] + rotmat[2]*frag_v[2][i];
          y = rotmat[3]*frag_v[0][i] + rotmat[4]*frag_v[1][i] + rotmat[5]*frag_v[2][i];
          z = rotmat[6]*frag_v[0][i] + rotmat[7]*frag_v[1][i] + rotmat[8]*frag_v[2][i];
          frag_v[0][i] = x;
          frag_v[1][i] = y;
          frag_v[2][i] = z;
        }

        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            vel_fitted[istep + 3*nstep*i + nstep*j] = frag_v[j][i];
          }
        }
      }

    }
  } else {
    for (istep = 0; istep < nstep; istep++) {

      for (i = 0; i < natom; i++) {
        for (j = 0; j < 3; j++) {
          frag_b[j][i] = trj[istep + 3*nstep*i + nstep*j];
        }
      }

      for (i = 0; i < len; i++) {
        for (j = 0; j < 3; j++) {
          frag_b_sub[j][i] = frag_b[j][index[i]];
        }
      }

      CenterCoords(frag_b_sub, &xsum_b, &ysum_b, &zsum_b, len, weight_sub);
      rmsd[istep] = CalcRMSDRotationalMatrix((double **) frag_a_sub, (double **) frag_b_sub, len, rotmat, weight_sub);

      if ( nlhs > 1 ) {
        /* apply rotation matrix */
        for (i = 0; i < natom; ++i) {
          frag_b[0][i] = frag_b[0][i] - xsum_b;
          frag_b[1][i] = frag_b[1][i] - ysum_b;
          frag_b[2][i] = frag_b[2][i] - zsum_b;

          x = rotmat[0]*frag_b[0][i] + rotmat[1]*frag_b[1][i] + rotmat[2]*frag_b[2][i];
          y = rotmat[3]*frag_b[0][i] + rotmat[4]*frag_b[1][i] + rotmat[5]*frag_b[2][i];
          z = rotmat[6]*frag_b[0][i] + rotmat[7]*frag_b[1][i] + rotmat[8]*frag_b[2][i];
        
          frag_b[0][i] = x + xsum_a;
          frag_b[1][i] = y + ysum_a;
          frag_b[2][i] = z + zsum_a;
        }

        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            trj_fitted[istep + 3*nstep*i + nstep*j] = frag_b[j][i];
          }
        }
      }

      if ( nlhs > 2 ) {
        /* apply rotation matrix to velocity*/
        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            frag_v[j][i] = vel[istep + 3*nstep*i + nstep*j];
          }
        }

        for (i = 0; i < natom; ++i) {
          x = rotmat[0]*frag_v[0][i] + rotmat[1]*frag_v[1][i] + rotmat[2]*frag_v[2][i];
          y = rotmat[3]*frag_v[0][i] + rotmat[4]*frag_v[1][i] + rotmat[5]*frag_v[2][i];
          z = rotmat[6]*frag_v[0][i] + rotmat[7]*frag_v[1][i] + rotmat[8]*frag_v[2][i];
          frag_v[0][i] = x;
          frag_v[1][i] = y;
          frag_v[2][i] = z;
        }

        for (i = 0; i < natom; i++) {
          for (j = 0; j < 3; j++) {
            vel_fitted[istep + 3*nstep*i + nstep*j] = frag_v[j][i];
          }
        }
      }

    }
  }

  MatDestroy(&frag_a);
  MatDestroy(&frag_b);
  MatDestroy(&frag_a_sub);
  MatDestroy(&frag_b_sub);
  MatDestroy(&frag_v);

  if (index != NULL) {
    free(index);
  }

  if (weight_sub != NULL) {
    free(weight_sub);
  }

  /* exit(EXIT_SUCCESS); */
}

