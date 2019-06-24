/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <stdio.h>
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "SMsmooth.h"

#undef __FUNC__
#define __FUNC__ "SMformGrammian"
int SMformGrammian(int num_active, double **vec, double **G, int dimension)
{
   int ierr;
   int i, j;

   if (num_active > OPTMS_MAX_G_NUM) {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Exceeded maximum allowed active values");
   }
   /* form the grammian with the dot products of the gradients */
   for (i=0; i<num_active; i++) {
      for (j=i; j<num_active; j++) {
         G[i][j] = 0.;
	 OPTMS_DOT(G[i][j],vec[i],vec[j],dimension);
	 G[j][i] = G[i][j];
      }
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMformPDGrammian"
int SMformPDGrammian(SMoptimal *opt_info)
{
    int ierr;
    int    i,j,k,g_ind_1;
    int    dimension;
    int    num_active, *active_ind, *PDG_ind;
    int    singular;
    double **G, **PDG;

    OPTMS_CHECK_NULL(opt_info);

    dimension = opt_info->dimension;

    num_active = opt_info->active->num_active;
    active_ind = opt_info->active->active_ind;
    PDG_ind = opt_info->PDG_ind;
    G = opt_info->G;
    PDG = opt_info->PDG;

    /* this assumes the grammian has been formed */
    for (i=0;i<num_active;i++) {
      for (j=0;j<num_active;j++) {
        if (G[i][j]==-1) OPTMS_SETERR(OPTMS_INIT_ERR,0,"Grammian not computed properly");
      }
    }

    /* use the first gradient in the active set */
    g_ind_1 = 0;
    PDG[0][0] = G[0][0];
    PDG_ind[0] = active_ind[0];

    /* test the rest and add them as appropriate */
    k = 1; i = 1;
    while( (k<dimension) && (i < num_active) ) {
        PDG[0][k] = PDG[k][0] = G[0][i];
        PDG[k][k] = G[i][i];
        if ( k == 2) { /* add the dot product of g1 and g2 */
           PDG[1][k] = PDG[k][1] = G[g_ind_1][i];
        }
        ierr = SMsingularTest(k+1,PDG,&singular); OPTMS_CHKERR(ierr);
        if (!singular) {
           PDG_ind[k] = active_ind[i];
           if (k==1) g_ind_1 = i;
           k++;
        }
        i++;
    }
    opt_info->num_LI = k;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMformVerticalMatrix"
int SMformVerticalMatrix(int num, double **PDG, double ***N)
{
    int ierr;
    int i,j;

    OPTMS_MALLOC((*N),(double **),sizeof(double *)*num,1);
    for (i=0;i<num;i++) OPTMS_MALLOC((*N)[i],(double *),sizeof(double)*num,1);

    for (i=0;i<num;i++) {
      (*N)[i][i] = PDG[i][i] + 1;
      for (j=i+1;j<num;j++) {
         (*N)[i][j] = (*N)[j][i] =  PDG[i][j] + 1;
      }
    }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMformReducedMatrix"
int SMformReducedMatrix(int num_active,double **G, double ***P)
{
    int ierr;
    int i,j;

    OPTMS_MALLOC((*P),(double **),sizeof(double *)*(num_active-1),1);
    for (i=0; i<num_active-1; i++)
        OPTMS_MALLOC((*P)[i],(double *),sizeof(double)*(num_active-1),1);

    for (i=0;i<num_active-1;i++) {
        (*P)[i][i] = G[0][0] - 2*G[0][i+1] + G[i+1][i+1];
        for (j=i+1;j<num_active-1;j++) {
            (*P)[i][j] = G[0][0] - G[0][j+1] - G[i+1][0] + G[i+1][j+1];
            (*P)[j][i] = (*P)[i][j];
        }
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMsingularTest"
int SMsingularTest(int n, double **A, int *singular)
{
    int ierr;
    double cond;

    if ((n>3) || (n<1)) {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Singular test works only for n=1 to n=3");
    }

    (*singular)=OPTMS_TRUE;
    switch(n) {
    case 1:
        if (A[0][0] > 0) (*singular) = OPTMS_FALSE;
        break;
    case 2:
        if (fabs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) > OPTMS_MACHINE_EPS)
            (*singular) = OPTMS_FALSE;
        break;
    case 3:
        /* calculate the condition number */
         ierr = SMcondition3x3(A, &cond); OPTMS_CHKERR(ierr);
 	 if (cond < 1E14) (*singular)=OPTMS_FALSE;
         break;
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMsolve2x2"
int SMsolve2x2(double a11, double a12, double a21, double a22,
		   double b1, double b2, double **x)
{
    int ierr;
    double factor;

    /* if the system is not singular, solve it */
    if (fabs(a11*a22 - a21*a12) > OPTMS_MACHINE_EPS) {
	OPTMS_MALLOC((*x),(double *),sizeof(double)*2,1);
	if (fabs(a11) > OPTMS_MACHINE_EPS) {
	    factor = (a21/a11);
	    (*x)[1] = (b2 - factor*b1)/(a22 - factor*a12);
	    (*x)[0] = (b1 - a12*(*x)[1])/a11;
	} else if (fabs(a21) > OPTMS_MACHINE_EPS) {
	    factor = (a11/a21);
	    (*x)[1] = (b1 - factor*b2)/(a12 - factor*a22);
	    (*x)[0] = (b2 - a22*(*x)[1])/a21;
	}
    } else {
	(*x) = NULL;
    }
    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMsolveSPD3x3"
int SMsolveSPD3x3(double **A, double *B, double **x)
{
    int ierr;
    double L11, L21, L31, L22, L32, L33;
    double Y1, Y2, Y3;

    /* solves a 3x3 symmetric positive definite system
       using cholesky factorization */

    OPTMS_MALLOC((*x),(double *),sizeof(double)*3,1);

    /* Find the cholesky factor L */
    L11 = sqrt(A[0][0]);
    L21 = A[1][0] / L11;
    L22 = sqrt(A[1][1] - L21*L21);
    L31 = A[2][0] / L11;
    L32 = (A[2][1] - L31*L21) / L22;
    L33 = sqrt(A[2][2] - L31*L31 - L32*L32);

    if (fabs(L11) < OPTMS_MACHINE_EPS) OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    if (fabs(L22) < OPTMS_MACHINE_EPS) OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    if (fabs(L33) < OPTMS_MACHINE_EPS) OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);

    /* Find the solution to the lower triangular system */
    Y1 = B[0] / L11;
    Y2 = (B[1] - L21*Y1) / L22;
    Y3 = (B[2] - L31*Y1 - L32*Y2) / L33;

    /* Find the solution of the upper triangular system */
    (*x)[2] = Y3 / L33;
    (*x)[1] = (Y2 - L32*(*x)[2]) / L22;
    (*x)[0] = (Y1 - L21*(*x)[1] - L31*(*x)[2]) / L11;

    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMtransposeMatrix"
int  SMtransposeMatrix(double *mat, int rows, int cols, double *mat_T)
{
   int ierr;
   int i,j,k;

   OPTMS_CHECK_NULL(mat);
   OPTMS_CHECK_NULL(mat_T);

   k = 0;
   for (j=0;j<cols;j++) {
     for (i=0;i<rows;i++)   mat_T[k++] = mat[i*rows+j];
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMtransposeMatrix2"
int  SMtransposeMatrix2(double **mat, int rows, int cols, double **mat_T)
{
   int ierr;
   int i,j;

   OPTMS_CHECK_NULL(mat);
   OPTMS_CHECK_NULL(mat_T);

   for (j=0;j<cols;j++) {
     for (i=0;i<rows;i++)  mat_T[j][i] = mat[i][j];
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcondition3x3"
int SMcondition3x3(double **A, double *cond)
{
   int ierr;
   double a11, a12, a13;
   double a21, a22, a23;
   double a31, a32, a33;
   double s1, s2, s3;
   double denom;
   double temp;
   int zero_denom = OPTMS_TRUE;

   a11 = A[0][0]; a12=A[0][1]; a13=A[0][2];
   a21 = A[1][0]; a22=A[1][1]; a23=A[1][2];
   a31 = A[2][0]; a32=A[2][1]; a33=A[2][2];

   denom = -a11*a22*a33+a11*a23*a32+a21*a12*a33-a21*a13*a32-
            a31*a12*a23+a31*a13*a22;

   if ( (fabs(a11) > OPTMS_MACHINE_EPS) &&
        (fabs(denom/a11) > OPTMS_MACHINE_EPS)) {
         zero_denom = OPTMS_FALSE;
   }
   if ( (fabs(a22) > OPTMS_MACHINE_EPS) &&
        (fabs(denom/a22) > OPTMS_MACHINE_EPS)) {
         zero_denom = OPTMS_FALSE;
   }
   if ( (fabs(a33) > OPTMS_MACHINE_EPS) &&
        (fabs(denom/a33) > OPTMS_MACHINE_EPS)) {
         zero_denom = OPTMS_FALSE;
   }

   if (zero_denom) {
     (*cond) = OPTMS_BIG_POS_NMBR;
   } else {
     s1 = sqrt(a11*a11 + a12*a12 + a13*a13 +
               a21*a21 + a22*a22 + a23*a23 +
               a31*a31 + a32*a32 + a33*a33);

     temp = (-a22*a33+a23*a32)/denom;
     s3 = temp*temp;
     temp =(a12*a33-a13*a32)/denom;
     s3 += temp*temp;
     temp = (a12*a23-a13*a22)/denom;
     s3 += temp*temp;
     temp = (a21*a33-a23*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a33-a13*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a23-a13*a21)/denom;
     s3 += temp*temp;
     temp = (a21*a32-a22*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a32+a12*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a22+a12*a21)/denom;
     s3 += temp*temp;

     s2 = sqrt(s3);
     (*cond) = s1*s2;
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMdeterminant2x2"
int SMdeterminant2x2(double a1[2], double a2[2], double *det)
{
  int ierr;
  *det = a1[0] * a2[1] - a1[1] * a2[0];
  return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfrobenius_norm_squared2x2"
int SMfrobenius_norm_squared2x2(double a1[2], double a2[2], double *norma)
{
  int ierr;

  *norma  = a1[0]*a1[0] + a1[1]*a1[1];
  *norma += a2[0]*a2[0] + a2[1]*a2[1];

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMadjoint2x2"
int SMadjoint2x2(double a1[2], double a2[2],
                 double b1[2], double b2[2])
{
  int ierr;
  b1[0] = a2[1];
  b1[1] = -a1[1];
  b2[0] = -a2[0];
  b2[1] = a1[0];
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmultiply2x2"
int SMmultiply2x2(double a1[2], double a2[2],
                   double b1[2], double b2[2],
                   double r1[2], double r2[2])
{
  int ierr;
  r1[0] = a1[0] * b1[0] + a2[0] * b1[1];
  r1[1] = a1[1] * b1[0] + a2[1] * b1[1];

  r2[0] = a1[0] * b2[0] + a2[0] * b2[1];
  r2[1] = a1[1] * b2[0] + a2[1] * b2[1];
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMtranspose2x2"
int SMtranspose2x2(double a1[2], double a2[2],
		    double b1[2], double b2[2])
{
  int ierr;
  b1[0] = a1[0];  b1[1] = a2[0];
  b2[0] = a1[1];  b2[1] = a2[1];
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMdeterminant3x3"
int SMdeterminant3x3(double a1[3], double a2[3], double a3[3], double *determinant)
{
  int ierr;

  *determinant  = a1[0] * (a2[1]*a3[2]-a3[1]*a2[2]);
  *determinant -= a2[0] * (a1[1]*a3[2]-a3[1]*a1[2]);
  *determinant += a3[0] * (a1[1]*a2[2]-a2[1]*a1[2]);

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmultiply3x3"
int SMmultiply3x3(double a1[3], double a2[3], double a3[3],
                   double b1[3], double b2[3], double b3[3],
                   double r1[3], double r2[3], double r3[3])
{
  int ierr;
  r1[0] = a1[0] * b1[0] + a2[0] * b1[1] + a3[0] * b1[2];
  r1[1] = a1[1] * b1[0] + a2[1] * b1[1] + a3[1] * b1[2];
  r1[2] = a1[2] * b1[0] + a2[2] * b1[1] + a3[2] * b1[2];

  r2[0] = a1[0] * b2[0] + a2[0] * b2[1] + a3[0] * b2[2];
  r2[1] = a1[1] * b2[0] + a2[1] * b2[1] + a3[1] * b2[2];
  r2[2] = a1[2] * b2[0] + a2[2] * b2[1] + a3[2] * b2[2];

  r3[0] = a1[0] * b3[0] + a2[0] * b3[1] + a3[0] * b3[2];
  r3[1] = a1[1] * b3[0] + a2[1] * b3[1] + a3[1] * b3[2];
  r3[2] = a1[2] * b3[0] + a2[2] * b3[1] + a3[2] * b3[2];
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMtranspose3x3"
int SMtranspose3x3(double a1[3], double a2[3], double a3[3],
		    double b1[3], double b2[3], double b3[3])
{
  int ierr;
  b1[0] = a1[0];  b1[1] = a2[0];  b1[2] = a3[0];
  b2[0] = a1[1];  b2[1] = a2[1];  b2[2] = a3[1];
  b3[0] = a1[2];  b3[1] = a2[2];  b3[2] = a3[2];
  return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMadjoint3x3"
int SMadjoint3x3(double a1[3], double a2[3], double a3[3],
                  double b1[3], double b2[3], double b3[3])
{
  int ierr;
  b1[0] = a2[1]*a3[2]-a3[1]*a2[2];
  b1[1] = a3[1]*a1[2]-a1[1]*a3[2];
  b1[2] = a1[1]*a2[2]-a2[1]*a1[2];

  b2[0] = a3[0]*a2[2]-a2[0]*a3[2];
  b2[1] = a1[0]*a3[2]-a3[0]*a1[2];
  b2[2] = a2[0]*a1[2]-a1[0]*a2[2];

  b3[0] = a2[0]*a3[1]-a3[0]*a2[1];
  b3[1] = a3[0]*a1[1]-a1[0]*a3[1];
  b3[2] = a1[0]*a2[1]-a2[0]*a1[1];
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfrobenius_norm_squared3x3"
int SMfrobenius_norm_squared3x3(double a1[3], double a2[3], double a3[3], double *norma)
{
  int ierr;

  *norma  = a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2];
  *norma += a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2];
  *norma += a3[0]*a3[0] + a3[1]*a3[1] + a3[2]*a3[2];

  return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMfrobenius_norm_squared_adj3x3"
int SMfrobenius_norm_squared_adj3x3(double a1[3], double a2[3], double a3[3],
               double *norm_adja)
{
  int ierr;
  double tmp[3];

  tmp[0] =  a1[1]*a2[2] - a1[2]*a2[1];
  tmp[1] = -(a1[0]*a2[1] - a1[2]*a2[0]);
  tmp[2] =  a1[0]*a2[1] - a1[1]*a2[0];
  *norm_adja = tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];

  tmp[0] =  a2[1]*a3[2] - a2[2]*a3[1];
  tmp[1] = -(a2[0]*a3[1] - a2[2]*a3[0]);
  tmp[2] =  a2[0]*a3[1] - a2[1]*a3[0];
  *norm_adja += tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];

  tmp[0] =  a3[1]*a1[2] - a3[2]*a1[1];
  tmp[1] = -(a3[0]*a1[1] - a3[2]*a1[0]);
  tmp[2] =  a3[0]*a1[1] - a3[1]*a1[0];
  *norm_adja += tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];

  return (ierr=0);
}
