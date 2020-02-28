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
#define __FUNC__ "SMcomputeTriCosines"
int SMcomputeTriCosines(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double xlen, ylen;
    double L1, L2, L3;

    xlen = vtx2[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx2[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L1 = xlen*xlen + ylen*ylen;

    xlen = vtx1[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L2 = xlen*xlen + ylen*ylen;

    xlen = vtx1[OPTMS_XDIR] - vtx2[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx2[OPTMS_YDIR];
    L3 = xlen*xlen + ylen*ylen;

    if (L1<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    function[0] = .5*(L2+L3-L1)/(sqrt(L2)*sqrt(L3));
    function[1] = .5*(L1+L3-L2)/(sqrt(L1)*sqrt(L3));
    function[2] = .5*(L1+L2-L3)/(sqrt(L1)*sqrt(L2));
    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorTriCosines"
int SMcomputeInteriorTriCosines(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double xlen, ylen;
    double L1, L2, L3;

    xlen = vtx2[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx2[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L1 = xlen*xlen + ylen*ylen;

    xlen = vtx1[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L2 = xlen*xlen + ylen*ylen;

    xlen = vtx1[OPTMS_XDIR] - vtx2[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx2[OPTMS_YDIR];
    L3 = xlen*xlen + ylen*ylen;

    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    function[0] = .5*(L2+L3-L1)/(sqrt(L2)*sqrt(L3));
    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeNegTriCosines"
int SMcomputeNegTriCosines(double *vtx1, double *vtx2,
                  double *vtx3, double *function, int *num_values)
{
    int ierr;
    int i;
    ierr = SMcomputeTriCosines(vtx1, vtx2, vtx3, function, num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) function[i] = -function[i];
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeTriAngles"
int SMcomputeTriAngles(double *vtx1, double *vtx2, double *vtx3,
                                  double *function, int *num_values)
{
    int ierr;
    int i;

    ierr = SMcomputeTriCosines(vtx1, vtx2, vtx3, function, num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) function[i] = acos(function[i]);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorTriAngles"
int SMcomputeInteriorTriAngles(double *vtx1, double *vtx2, double *vtx3,
                                  double *function, int *num_values)
{
    int ierr;
    int i;

    ierr = SMcomputeInteriorTriCosines(vtx1, vtx2, vtx3, function, num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) function[i] = acos(function[i]);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeNegTriAngles"
int SMcomputeNegTriAngles(double *vtx1, double *vtx2,
                  double *vtx3, double *function, int *num_values)
{
    int ierr;
    int i;
    ierr = SMcomputeTriAngles(vtx1, vtx2, vtx3, function, num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) function[i] = -function[i];
    return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMcomputeTriSines"
int SMcomputeTriSines(double *vtx1, double *vtx2, double *vtx3,
                                 double *function, int *num_values)
{
    int ierr;
    double xlen, ylen;
    double L1, L2, L3;
    double sL1, sL2, sL3;
    double cos0;

    /* Computes the sine of the angles of the triangle, but only returns
        positive values - as it stands this function cannot be used for mesh
        untangling */

    xlen = vtx2[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx2[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L1 = xlen*xlen + ylen*ylen;
    sL1 = sqrt(L1);

    xlen = vtx1[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L2 = xlen*xlen + ylen*ylen;
    sL2 = sqrt(L2);

    xlen = vtx1[OPTMS_XDIR] - vtx2[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx2[OPTMS_YDIR];
    L3 = xlen*xlen + ylen*ylen;
    sL3 = sqrt(L3);

    if (L1<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    cos0 = .5*(L2+L3-L1)/(sL2*sL3);
    function[0] = sqrt(1.0 - cos0*cos0);
    function[1] = (sL2 * function[0]) / sL1;
    function[2] = (sL3 * function[0]) / sL1;
    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorTriSines"
int SMcomputeInteriorTriSines(double *vtx1, double *vtx2, double *vtx3,
                                 double *function, int *num_values)
{
    int ierr;
    double xlen, ylen;
    double L1, L2, L3;
    double sL1, sL2, sL3;
    double cos0;

    /* Computes the sine of the angles of the triangle, but only returns
        positive values - as it stands this function cannot be used for mesh
        untangling */

    xlen = vtx2[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx2[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L1 = xlen*xlen + ylen*ylen;
    sL1 = sqrt(L1);

    xlen = vtx1[OPTMS_XDIR] - vtx3[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    L2 = xlen*xlen + ylen*ylen;
    sL2 = sqrt(L2);

    xlen = vtx1[OPTMS_XDIR] - vtx2[OPTMS_XDIR];
    ylen = vtx1[OPTMS_YDIR] - vtx2[OPTMS_YDIR];
    L3 = xlen*xlen + ylen*ylen;
    sL3 = sqrt(L3);

    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    cos0 = .5*(L2+L3-L1)/(sL2*sL3);
    function[0] = sqrt(1.0 - cos0*cos0);
    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeTriJacobians"
int SMcomputeTriJacobians(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double a, b, c;
    double equal_jac;
    double L1_2;
    double x1, x2, x3, y1, y2, y3;

    /* Computes the square of the difference between the actual jacbian
        and the jacobian of an equilateral triangle based on the edge opposite
        the free vertex.  We need to minimize the maximum difference and
        the function is therefore negated.  The corresponding gradient code
        is in SMcomputeJacobianGradients.  There is one value per triangle */

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1_2 = (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3);
    equal_jac = .5*sqrt(3)*L1_2;

    if (L1_2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    a = x2*y3 - x3*y2;
    b = y2 - y3;
    c = x3 - x2;
    function[0] = a + b*x1 + c*y1;
    function[0] = -((function[0] - equal_jac)/equal_jac)*
                   ((function[0]-equal_jac)/equal_jac);

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeTriArea"
int SMcomputeTriArea(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double a, b, c;
    double x1, x2, x3, y1, y2, y3;

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    a = x2*y3 - x3*y2;
    b = y2 - y3;
    c = x3 - x2;

    function[0] = .5*(a + b*x1 + c*y1);
    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeScaledTriJacobians"
int SMcomputeScaledTriJacobians(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double a, b, c;
    double L1, L2, L3;
    double x1, x2, x3, y1, y2, y3;

    /* This function computes the scaled jacobian for each angle of the triangle
        -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
        Note that this is equivalent to the sine of the angle but is an efficient way to compute
        both positive and negative sines for mesh untangling.  There are 3 function values
        per triangle and we need to maximize the minimum scaled jacobian.  The
        corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
    L2 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
    L3 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

    if (L1<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    a = x2*y3 - x3*y2;    b = y2 - y3;    c = x3 - x2;
    function[0] = a + b*x1 + c*y1;
    function[0] = function[0]/(L2*L3);

    a = x3*y1 - x1*y3;    b = y3 - y1;    c = x1 - x3;
    function[1] = a + b*x2 + c*y2;
    function[1] = function[1]/(L1*L3);

    a = x1*y2 - x2*y1;    b = y1 - y2;    c = x2 - x1;
    function[2] = a + b*x3 + c*y3;
    function[2] = function[2]/(L1*L2);

    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorScaledTriJacobians"
int SMcomputeInteriorScaledTriJacobians(double *vtx1, double *vtx2, double *vtx3,
                                   double *function, int *num_values)
{
    int ierr;
    double a, b, c;
    double L1, L2, L3;
    double x1, x2, x3, y1, y2, y3;

    /* This function computes the scaled jacobian for each angle of the triangle
        -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
        Note that this is equivalent to the sine of the angle but is an efficient way to compute
        both positive and negative sines for mesh untangling.  There are 3 function values
        per triangle and we need to maximize the minimum scaled jacobian.  The
        corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
    L2 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
    L3 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

    if (L2<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    a = x2*y3 - x3*y2;    b = y2 - y3;    c = x3 - x2;
    function[0] = a + b*x1 + c*y1;
    function[0] = function[0]/(L2*L3);

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeAreaLengthRatio"
int SMcomputeAreaLengthRatio(double *vtx1, double *vtx2, double *vtx3,
                         double *function, int *num_values)
{
    int ierr;
    double a, b, c;
    double L1, L2, L3;
    double L1_xdiff, L2_xdiff, L3_xdiff;
    double L1_ydiff, L2_ydiff, L3_ydiff;
    double x1, x2, x3, y1, y2, y3;

    /* comptes the SMcomputeAreaLengthRatio function which
        gives the ratio of the the area of the triangle and the sum of the squares of
        the length of the edges */

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1_xdiff = x2-x3;  L1_ydiff=y2-y3;
    L2_xdiff = x1-x3;  L2_ydiff=y1-y3;
    L3_xdiff = x1-x2;  L3_ydiff=y1-y2;

    L1 = L1_xdiff*L1_xdiff + L1_ydiff*L1_ydiff;
    L2 = L2_xdiff*L2_xdiff + L2_ydiff*L2_ydiff;
    L3 = L3_xdiff*L3_xdiff + L3_ydiff*L3_ydiff;

    if ((L1+L2+L3)<OPTMS_MACHINE_EPS)
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    a = x2*y3 - x3*y2;    b = y2 - y3;    c = x3 - x2;
    function[0] = .5*(a + b*x1 + c*y1);
    function[0] =12./sqrt(3)*function[0]/(L1 + L2 + L3);

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeLengthAreaRatio"
int SMcomputeLengthAreaRatio(double *vtx1, double *vtx2, double *vtx3,
                         double *function, int *num_values)
{
   int ierr;
   int i;
   ierr = SMcomputeAreaLengthRatio(vtx1,vtx2,vtx3,function,num_values);
          OPTMS_CHKERR(ierr);
   for (i=0;i<(*num_values);i++)  function[i]=-1/function[i];
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeAngGradients"
int SMcomputeAngGradients(double *vtx1, double *vtx2, double *vtx3,
                         double **gradient, int *num_values)
{
    int ierr;
     int i;
     int num_vert = 3;
     double *cosines, tmp;

     OPTMS_MALLOC(cosines,(double *),sizeof(double)*num_vert,1);

     ierr = SMcomputeTriCosines(vtx1,vtx2,vtx3,cosines,num_values);
            OPTMS_CHKERR(ierr);
     ierr = SMcomputeCosGradients(vtx1,vtx2,vtx3,gradient,num_values);
            OPTMS_CHKERR(ierr);

     /* compute the gradient of the angle */
     for (i=0;i<*num_values;i++){
         if ((fabs(cosines[i]) - 1)<OPTMS_MACHINE_EPS) {
             OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Linearly aligned vertices");
         }
         tmp = -1. / sqrt( 1 - cosines[i]*cosines[i]);
	 gradient[i][OPTMS_XDIR] *= tmp;
	 gradient[i][OPTMS_YDIR] *= tmp;
    }
    *num_values = num_vert;
    OPTMS_FREE(cosines);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorAngGradients"
int SMcomputeInteriorAngGradients(double *vtx1, double *vtx2, double *vtx3,
                         double **gradient, int *num_values)
{
    int ierr;
     int i;
     int num_vert = 3;
     double *cosines, tmp;

     OPTMS_MALLOC(cosines,(double *),sizeof(double)*num_vert,1);

     ierr = SMcomputeInteriorTriCosines(vtx1,vtx2,vtx3,cosines,num_values);
            OPTMS_CHKERR(ierr);
     ierr = SMcomputeInteriorCosGradients(vtx1,vtx2,vtx3,gradient,num_values);
            OPTMS_CHKERR(ierr);

     /* compute the gradient of the angle */
     for (i=0;i<*num_values;i++){
         if ((fabs(cosines[i]) - 1)<OPTMS_MACHINE_EPS) {
             OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Linearly aligned vertices");
         }
         tmp = -1. / sqrt( 1 - cosines[i]*cosines[i]);
	 gradient[i][OPTMS_XDIR] *= tmp;
	 gradient[i][OPTMS_YDIR] *= tmp;
    }
    *num_values = num_vert;
    OPTMS_FREE(cosines);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeNegAngGradients"
int SMcomputeNegAngGradients(double *vtx1, double *vtx2,
                        double *vtx3, double **gradient, int *num_values)
{
    int ierr;
    int i,j;

    ierr = SMcomputeAngGradients(vtx1,vtx2,vtx3,gradient,num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) {
      for (j=0;j<2;j++) gradient[i][j] = -gradient[i][j];
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeCosGradients"
int SMcomputeCosGradients(double *vtx1, double *vtx2, double *vtx3,
                                     double **gradient, int *num_values)
{
    int ierr;
    int i;
    int num_vert = 3;
    double x_1, x_2, x_3, y_1, y_2, y_3;
    double a, b, c;
    double a2, b2, c2;
    double gx, gy;
    double term1, term2;

    x_1 = vtx1[OPTMS_XDIR];    x_2 = vtx2[OPTMS_XDIR];    x_3 = vtx3[OPTMS_XDIR];
    y_1 = vtx1[OPTMS_YDIR];    y_2 = vtx2[OPTMS_YDIR];    y_3 = vtx3[OPTMS_YDIR];

    a2 = (x_2 - x_1)*(x_2 - x_1) + (y_2 - y_1)*(y_2 - y_1);
    b2 = (x_1 - x_3)*(x_1 - x_3) + (y_1 - y_3)*(y_1 - y_3);
    c2 = (x_3 - x_2)*(x_3 - x_2) + (y_3 - y_2)*(y_3 - y_2);

    a = sqrt(a2);  b = sqrt(b2);  c = sqrt(c2);

    if (a<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (b<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (c<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    /* compute the gradient of the cosine */
    for (i=0;i<num_vert;i++) {
	if (i==0) {
	    term1 = (a2 + c2 - b2) / (2*a2*b) * (1/a);
	    term2 = (b2 + c2 - a2) / (2*a*b2) * (1/b);
	    /* SM_LOG_FLOPS(__SM_GRADIENT__,14); */
	} else if (i==1) {
	    term1 = (a2 + b2 - c2) / (2*a2*c) * (1/a);
            term2 = -1 / (a*c);
	    /* SM_LOG_FLOPS(__SM_GRADIENT__,9); */
        } else {
	    term1 = -1 / (b*c);
	    term2 = (b2 + a2 - c2) / (2*b2*c) * (1/b);
	    /* SM_LOG_FLOPS(__SM_GRADIENT__,9); */
	}

	gx = term1*(x_1 - x_2) + term2*(x_1 - x_3);
	gy = term1*(y_1 - y_2) + term2*(y_1 - y_3);
    	/* SM_LOG_FLOPS(__SM_GRADIENT__,10); */

	gradient[i][OPTMS_XDIR] = gx;
	gradient[i][OPTMS_YDIR] = gy;
    }
    *num_values = num_vert;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorCosGradients"
int SMcomputeInteriorCosGradients(double *vtx1, double *vtx2, double *vtx3,
                                     double **gradient, int *num_values)
{
    int ierr;
    double x_1, x_2, x_3, y_1, y_2, y_3;
    double a, b, c;
    double a2, b2, c2;
    double gx, gy;
    double term1, term2;

    x_1 = vtx1[OPTMS_XDIR];    x_2 = vtx2[OPTMS_XDIR];    x_3 = vtx3[OPTMS_XDIR];
    y_1 = vtx1[OPTMS_YDIR];    y_2 = vtx2[OPTMS_YDIR];    y_3 = vtx3[OPTMS_YDIR];

    a2 = (x_2 - x_1)*(x_2 - x_1) + (y_2 - y_1)*(y_2 - y_1);
    b2 = (x_1 - x_3)*(x_1 - x_3) + (y_1 - y_3)*(y_1 - y_3);
    c2 = (x_3 - x_2)*(x_3 - x_2) + (y_3 - y_2)*(y_3 - y_2);

    a = sqrt(a2);  b = sqrt(b2);  c = sqrt(c2);

    if (a<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (b<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    /* compute the gradient of the cosine */
    term1 = (a2 + c2 - b2) / (2*a2*b) * (1/a);
    term2 = (b2 + c2 - a2) / (2*a*b2) * (1/b);

    gx = term1*(x_1 - x_2) + term2*(x_1 - x_3);
    gy = term1*(y_1 - y_2) + term2*(y_1 - y_3);

    gradient[0][OPTMS_XDIR] = gx;
    gradient[0][OPTMS_YDIR] = gy;
    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeNegCosGradients"
int SMcomputeNegCosGradients(double *vtx1, double *vtx2,
                        double *vtx3, double **gradient, int *num_values)
{
    int ierr;
    int i,j;

    ierr = SMcomputeCosGradients(vtx1,vtx2,vtx3,gradient,num_values);
           OPTMS_CHKERR(ierr);
    for (i=0;i<*num_values;i++) {
      for (j=0;j<2;j++) gradient[i][j] = -gradient[i][j];
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeSineGradients"
int SMcomputeSineGradients(double *vtx1, double *vtx2,
                        double *vtx3, double **gradient, int *num_values)
{
    int ierr;
    int i;
    double sqrtcos0, sqrtcos1, sqrtcos2;
    double *cosine;
    double **cos_gradient;
    int num_vert;

    num_vert = 3;

    OPTMS_MALLOC(cosine,(double *),sizeof(double)*num_vert,1);
    OPTMS_MALLOC(cos_gradient,(double **),sizeof(double *)*num_vert,1);
    for (i=0;i<num_vert;i++)
       OPTMS_MALLOC(cos_gradient[i],(double *),sizeof(double)*2,1);

    ierr = SMcomputeTriCosines(vtx1,vtx2,vtx3,cosine,num_values);
           OPTMS_CHKERR(ierr);
    ierr = SMcomputeCosGradients(vtx1,vtx2,vtx3,cos_gradient,num_values);
           OPTMS_CHKERR(ierr);

    sqrtcos0 = sqrt(1. - cosine[0]*cosine[0]);
    sqrtcos1 = sqrt(1. - cosine[1]*cosine[1]);
    sqrtcos2 = sqrt(1. - cosine[2]*cosine[2]);

    if ( (sqrtcos0<OPTMS_MACHINE_EPS) ||
         (sqrtcos1<OPTMS_MACHINE_EPS) ||
         (sqrtcos2<OPTMS_MACHINE_EPS)) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    }

    gradient[0][OPTMS_XDIR] = (-1./sqrtcos0) * cosine[0] * cos_gradient[0][OPTMS_XDIR];
    gradient[0][OPTMS_YDIR] = (-1./sqrtcos0) * cosine[0] * cos_gradient[0][OPTMS_YDIR];
    gradient[1][OPTMS_XDIR] = (-1./sqrtcos1) * cosine[1] * cos_gradient[1][OPTMS_XDIR];
    gradient[1][OPTMS_YDIR] = (-1./sqrtcos1) * cosine[1] * cos_gradient[1][OPTMS_YDIR];
    gradient[2][OPTMS_XDIR] = (-1./sqrtcos2) * cosine[2] * cos_gradient[2][OPTMS_XDIR];
    gradient[2][OPTMS_YDIR] = (-1./sqrtcos2) * cosine[2] * cos_gradient[2][OPTMS_YDIR];

    for (i=0;i<num_vert;i++) OPTMS_FREE(cos_gradient[i]);
    OPTMS_FREE(cosine);  OPTMS_FREE(cos_gradient);

    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorSineGradients"
int SMcomputeInteriorSineGradients(double *vtx1, double *vtx2,
                        double *vtx3, double **gradient, int *num_values)
{
    int ierr;
    int i;
    double sqrtcos0, sqrtcos1, sqrtcos2;
    double *cosine;
    double **cos_gradient;
    int num_vert;

    num_vert = 1;

    OPTMS_MALLOC(cosine,(double *),sizeof(double)*num_vert,1);
    OPTMS_MALLOC(cos_gradient,(double **),sizeof(double *)*num_vert,1);
    for (i=0;i<num_vert;i++)
       OPTMS_MALLOC(cos_gradient[i],(double *),sizeof(double)*2,1);

    SMcomputeInteriorTriCosines(vtx1,vtx2,vtx3,cosine,num_values);
    SMcomputeInteriorCosGradients(vtx1,vtx2,vtx3,cos_gradient,num_values);

    sqrtcos0 = sqrt(1. - cosine[0]*cosine[0]);
    sqrtcos1 = sqrt(1. - cosine[1]*cosine[1]);
    sqrtcos2 = sqrt(1. - cosine[2]*cosine[2]);

    if (sqrtcos0<OPTMS_MACHINE_EPS)
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    gradient[0][OPTMS_XDIR] = (-1./sqrtcos0) * cosine[0] * cos_gradient[0][OPTMS_XDIR];
    gradient[0][OPTMS_YDIR] = (-1./sqrtcos0) * cosine[0] * cos_gradient[0][OPTMS_YDIR];

    for (i=0;i<num_vert;i++) OPTMS_FREE(cos_gradient[i]);
    OPTMS_FREE(cosine);  OPTMS_FREE(cos_gradient);

    *num_values = 1;
    return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMcomputeJacobianGradients"
int SMcomputeJacobianGradients(double *vtx1, double *vtx2, double *vtx3,
                                   double **gradient, int *num_values)
{
    int ierr;
    double a, b, c;
    double equal_jac;
    double L1;

    L1 = (vtx2[OPTMS_XDIR]-vtx3[OPTMS_XDIR])*(vtx2[OPTMS_XDIR]-vtx3[OPTMS_XDIR]);
    L1 = L1 + (vtx2[OPTMS_YDIR]-vtx3[OPTMS_YDIR])*(vtx2[OPTMS_YDIR]-vtx3[OPTMS_YDIR]);

    if (L1<OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    equal_jac = .5*sqrt(3)*L1;

    a = vtx2[OPTMS_XDIR]*vtx3[OPTMS_YDIR] - vtx3[OPTMS_XDIR]*vtx2[OPTMS_YDIR];
    b = vtx2[OPTMS_YDIR] - vtx3[OPTMS_YDIR];
    c = vtx3[OPTMS_XDIR] - vtx2[OPTMS_XDIR];

    gradient[0][OPTMS_XDIR] = -2*(a+b*vtx1[OPTMS_XDIR]+c*vtx1[OPTMS_YDIR]-equal_jac)*b/equal_jac;
    gradient[0][OPTMS_YDIR] = -2*(a+b*vtx1[OPTMS_XDIR]+c*vtx1[OPTMS_YDIR]-equal_jac)*c/equal_jac;

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeScaledJacobianGradients"
int SMcomputeScaledJacobianGradients(double *vtx1, double *vtx2, double *vtx3,
                                   double **gradient, int *num_values)
{
    int ierr;
    double  f;
    double L1, L2, L3;
    double L1_2, L2_2, L3_2;
    double dL1, dL2, dL3;
    double x1, x2, x3, y1, y2, y3;

    /* This function computes the scaled jacobian for each angle of the
       triangle -- Jac/(L1*L2) where L1 and L2 are the length of the
       edges containing the angle. Note that this is equivalent to the
       sine of the angle but is an efficient way to compute both
       positive and negative sines for mesh untangling.  There are 3
       function values per triangle and we need to maximize the minimum
       scaled jacobian.  The corresponding gradient code is in the function
       SMcomputeScaledJacobianGradients*/

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1_2 = (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3);
    L1 = sqrt(L1_2);
    dL1=L1_2*L1;

    L2_2 = (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3);
    L2 = sqrt(L2_2);
    dL2=L2_2*L2;

    L3_2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
    L3 = sqrt(L3_2);
    dL3=L3_2*L3;


    if (L1 < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L2 < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3 < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    f = x2*y3-x3*y2 + (y2-y3)*x1+(x3-x2)*y1;
    gradient[0][OPTMS_XDIR] = (y2-y3)/(L3*L2) -f*(x1-x2)/(dL3*L2)
                                    -f*(x1-x3)/(L3*dL2);
    gradient[0][OPTMS_YDIR] = (x3-x2)/(L3*L2)-f*(y1-y2)/(dL3*L2)
                                   -f*(y1-y3)/(L3*dL2);

    f = x3*y1-x1*y3+(y3-y1)*x2+(x1-x3)*y2;
    gradient[1][OPTMS_XDIR] = (y2-y3)/(L1*L3)-f*(x1-x2)/(L1*dL3);
    gradient[1][OPTMS_YDIR] = (x3-x2)/(L1*L3)-f*(y1-y2)/(L1*dL3);

    f = x1*y2-x2*y1+(y1-y2)*x3+(x2-x1)*y3;
    gradient[2][OPTMS_XDIR] = (y2-y3)/(L2*L1)-f*(x1-x3)/(dL2*L1);
    gradient[2][OPTMS_YDIR] = (x3-x2)/(L2*L1)-f*(y1-y3)/(dL2*L1);

    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeInteriorScaledJacobianGradients"
int SMcomputeInteriorScaledJacobianGradients(double *vtx1, double *vtx2, double *vtx3,
                                   double **gradient, int *num_values)
{
    int ierr;
    double  f;
    double L1, L2, L3;
    double L1_2, L2_2, L3_2;
    double dL1, dL2, dL3;
    double x1, x2, x3, y1, y2, y3;

    /* This function computes the scaled jacobian for each angle of the
       triangle  Jac/(L1*L2) where L1 and L2 are the length of the
       edges containing the angle.  Note that this is equivalent to
       the sine of the angle but is an efficient way to compute both
       positive and negative sines for mesh untangling.  There are 3
       function values per triangle and we need to maximize the minimum
       scaled jacobian.  The corresponding gradient code is in the function
       SMcomputeScaledJacobianGradients*/

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1_2 = (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3);
    L1 = sqrt(L1_2);
    dL1=L1_2*L1;

    L2_2 = (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3);
    L2 = sqrt(L2_2);
    dL2=L2_2*L2;

    L3_2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
    L3 = sqrt(L3_2);
    dL3=L3_2*L3;

    if (L2 < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    if (L3 < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");

    f = x2*y3-x3*y2 + (y2-y3)*x1+(x3-x2)*y1;
    gradient[0][OPTMS_XDIR] = (y2-y3)/(L3*L2) -f*(x1-x2)/(dL3*L2)
                                    -f*(x1-x3)/(L3*dL2);
    gradient[0][OPTMS_YDIR] = (x3-x2)/(L3*L2)-f*(y1-y2)/(dL3*L2)
                                   -f*(y1-y3)/(L3*dL2);

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeAreaLengthRatioGradients"
int SMcomputeAreaLengthRatioGradients(double *vtx1, double *vtx2, double *vtx3,
                         double **gradient, int *num_values)
{
    int ierr;
    double a, b, c;
    double L1, L2, L3;
    double L1_xdiff, L2_xdiff, L3_xdiff;
    double L1_ydiff, L2_ydiff, L3_ydiff;
    double x1, x2, x3, y1, y2, y3;

    /* comptes the gradient of the SMcomputeAreaLengthRatio function which
       gives the ratio of the the area of the triangle and the sum of
       the squares of the length of the edges */

    x1=vtx1[OPTMS_XDIR];     y1=vtx1[OPTMS_YDIR];
    x2=vtx2[OPTMS_XDIR];     y2=vtx2[OPTMS_YDIR];
    x3=vtx3[OPTMS_XDIR];     y3=vtx3[OPTMS_YDIR];

    L1_xdiff = x2-x3;  L1_ydiff=y2-y3;
    L2_xdiff = x1-x3;  L2_ydiff=y1-y3;
    L3_xdiff = x1-x2;  L3_ydiff=y1-y2;

    L1 = L1_xdiff*L1_xdiff + L1_ydiff*L1_ydiff;
    L2 = L2_xdiff*L2_xdiff + L2_ydiff*L2_ydiff;
    L3 = L3_xdiff*L3_xdiff + L3_ydiff*L3_ydiff;

    a = x2*y3 - x3*y2;    b = y2 - y3;    c = x3 - x2;

    if ((L1+L2+L3)<OPTMS_MACHINE_EPS) {
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    }

    gradient[0][OPTMS_XDIR] = .5*b/(L1+L2+L3) -
                                 (.5*(a+b*x1+c*y1)*(4*x1-2*x3-2*x2))
                                 /((L1+L2+L3)*(L1+L2+L3));
    gradient[0][OPTMS_YDIR] = .5*c/(L1+L2+L3) -
                                 (.5*(a+b*x1+c*y1)*(4*y1-2*y3-2*y2))
                                 /((L1+L2+L3)*(L1+L2+L3));
    gradient[0][OPTMS_XDIR] = 12./sqrt(3)*gradient[0][OPTMS_XDIR];
    gradient[0][OPTMS_YDIR] = 12./sqrt(3)*gradient[0][OPTMS_YDIR];

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMNormJacSquared2D"
int SMNormJacSquared2D(double *vtx1, double *vtx2, double *vtx3,
                     double *function, int *num_values)
{
    int ierr;
    double a1[2], a2[2];

    /* the a matrix contains the Jacobian */
    a1[0] = vtx2[0] - vtx1[0];
    a1[1] = vtx2[1] - vtx1[1];

    a2[0] = vtx3[0] - vtx1[0];
    a2[1] = vtx3[1] - vtx1[1];

    ierr = SMfrobenius_norm_squared2x2(a1,a2,&function[0]); OPTMS_CHKERR(ierr);
    function[0] *= -1.;

    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeNormJacSquaredGradients2D"
int SMcomputeNormJacSquaredGradients2D(double *vtx1, double *vtx2,
             double *vtx3, double **gradient, int *num_values)
{
    int ierr;
    double a1[2], a2[2];

    /* the a matrix contains the Jacobian */
    a1[0] = vtx2[0] - vtx1[0];
    a1[1] = vtx2[1] - vtx1[1];

    a2[0] = vtx3[0] - vtx1[0];
    a2[1] = vtx3[1] - vtx1[1];

    gradient[0][0] = 2*(a1[0] + a2[0]);
    gradient[0][1] = 2*(a1[1] + a2[1]);
    *num_values = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcondition2D"
int SMcondition2D(double *vtx1, double *vtx2, double *vtx3,
                     double *function, int *num_values)
{
    int ierr;
    double a1[2], a2[2];
    double norm_a2;
    double alpha;

    /* the a matrix contains the Jacobian */
    /* note that the jacobian is not unique to the triangle, only
       the det(J) is.. we must therefore compute the condition number
       at each vertex */

    /* vertex 1 */
    a1[0] = vtx2[0] - vtx1[0];
    a1[1] = vtx2[1] - vtx1[1];
    a2[0] = vtx3[0] - vtx1[0];
    a2[1] = vtx3[1] - vtx1[1];

    ierr = SMdeterminant2x2(a1, a2, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared2x2(a1,a2, &norm_a2); OPTMS_CHKERR(ierr);
    function[0] = -.50 * norm_a2 / alpha;

    /* vertex 2 */
    a1[0] = vtx3[0] - vtx2[0];
    a1[1] = vtx3[1] - vtx2[1];
    a2[0] = vtx1[0] - vtx2[0];
    a2[1] = vtx1[1] - vtx2[1];

    ierr = SMfrobenius_norm_squared2x2(a1,a2,&norm_a2); OPTMS_CHKERR(ierr);
    function[1] = -.50 * norm_a2 / alpha;

    /* vertex 2 */
    a1[0] = vtx1[0] - vtx3[0];
    a1[1] = vtx1[1] - vtx3[1];
    a2[0] = vtx2[0] - vtx3[0];
    a2[1] = vtx2[1] - vtx3[1];

    ierr = SMfrobenius_norm_squared2x2(a1,a2,&norm_a2); OPTMS_CHKERR(ierr);
    function[2] = -.50 * norm_a2 / alpha;

    *num_values = 3;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMgradCondition2D"
int SMgradCondition2D(double *vtx1, double *vtx2, double *vtx3,
                               double **gradient, int *num_values)
{
    int ierr;
    double norma;
    double alpha;
    double a1[2], a2[2];
    double dnorma_dx, dnorma_dy;
    double dalpha_dx, dalpha_dy;

    /* the a matrix contains the Jacobian */
    /* vertex 1 */
    a1[0] = vtx2[0] - vtx1[0];
    a1[1] = vtx2[1] - vtx1[1];
    a2[0] = vtx3[0] - vtx1[0];
    a2[1] = vtx3[1] - vtx1[1];

    ierr = SMdeterminant2x2(a1, a2, &alpha);  OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared2x2(a1,a2,&norma); OPTMS_CHKERR(ierr);

    if ((norma < OPTMS_MACHINE_EPS) || (alpha < OPTMS_MACHINE_EPS)){
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Division by zero");
    }

    dnorma_dx = (2*vtx1[0]-vtx2[0]-vtx3[0])*(1/sqrt(norma));
    dnorma_dy = (2*vtx1[1]-vtx2[1]-vtx3[1])*(1/sqrt(norma));
    dalpha_dx = vtx2[1]-vtx3[1];
    dalpha_dy = vtx3[0]-vtx2[0];

    gradient[0][0] = -1*(sqrt(norma)/alpha)*dnorma_dx + .5*norma/(alpha*alpha)*dalpha_dx;
    gradient[0][1] = -1*(sqrt(norma)/alpha)*dnorma_dy + .5*norma/(alpha*alpha)*dalpha_dy;

    /* vertex 2 */
    a1[0] = vtx3[0] - vtx2[0];
    a1[1] = vtx3[1] - vtx2[1];
    a2[0] = vtx1[0] - vtx2[0];
    a2[1] = vtx1[1] - vtx2[1];

    ierr = SMfrobenius_norm_squared2x2(a1,a2,&norma); OPTMS_CHKERR(ierr);

    dnorma_dx = (vtx1[0]-vtx2[0])*(1/sqrt(norma));
    dnorma_dy = (vtx1[1]-vtx2[1])*(1/sqrt(norma));
    dalpha_dx = vtx2[1]-vtx3[1];
    dalpha_dy = vtx3[0]-vtx2[0];

    gradient[1][0] = -1*(sqrt(norma)/alpha)*dnorma_dx + .5*norma/(alpha*alpha)*dalpha_dx;
    gradient[1][1] = -1*(sqrt(norma)/alpha)*dnorma_dy + .5*norma/(alpha*alpha)*dalpha_dy;

    /* vertex 3 */
    a1[0] = vtx1[0] - vtx3[0];
    a1[1] = vtx1[1] - vtx3[1];
    a2[0] = vtx2[0] - vtx3[0];
    a2[1] = vtx2[1] - vtx3[1];

    ierr = SMfrobenius_norm_squared2x2(a1,a2,&norma); OPTMS_CHKERR(ierr);

    dnorma_dx = (vtx1[0]-vtx3[0])/sqrt(norma);
    dnorma_dy = (vtx1[1]-vtx3[1])/sqrt(norma);
    dalpha_dx = vtx2[1]-vtx3[1];
    dalpha_dy = vtx3[0]-vtx2[0];

    gradient[2][0] = -1*(sqrt(norma)/alpha)*dnorma_dx + .5*norma/(alpha*alpha)*dalpha_dx;
    gradient[2][1] = -1*(sqrt(norma)/alpha)*dnorma_dy + .5*norma/(alpha*alpha)*dalpha_dy;

    *num_values = 3;
    return(ierr=0);
}
