/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <assert.h>
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <stdio.h>
#include "SMsmooth.h"
#include "SMdihed_func.h"

#define iFuzzyOne(a) (fabs((a)-1) < 1.e-12)

#undef __FUNC__
#define __FUNC__ "vCross"
void vCross(const double adVecA[3], const double adVecB[3],
		   double adResult[3])
{
  adResult[0] = adVecA[1]*adVecB[2] - adVecA[2]*adVecB[1];
  adResult[1] = adVecA[2]*adVecB[0] - adVecA[0]*adVecB[2];
  adResult[2] = adVecA[0]*adVecB[1] - adVecA[1]*adVecB[0];
}

#undef __FUNC__
#define __FUNC__ "dMagnitude"
double dMagnitude(const double adVec[3])
{ 
  return (sqrt(adVec[0]*adVec[0] +
	       adVec[1]*adVec[1] +
	       adVec[2]*adVec[2]));
}

#undef __FUNC__
#define __FUNC__ "dDot"
double dDot(const double adVecA[3], const double adVecB[3])
{ 
   return (adVecA[0]*adVecB[0] + 
           adVecA[1]*adVecB[1] + 
           adVecA[2]*adVecB[2]); 
}

#undef __FUNC__
#define __FUNC__ "vUnitNormal"
int vUnitNormal(const double adCoord0[3], const double adCoord1[3],
			const double adCoord2[3], double adResult[3])
{
  int ierr;
  double adVecA[3], adVecB[3], dMag;
  int i;
  for (i = 0; i < 3; i++) {
    adVecA[i] = adCoord1[i] - adCoord0[i];
    adVecB[i] = adCoord2[i] - adCoord0[i];
  }

  vCross(adVecA, adVecB, adResult);

  dMag = dMagnitude(adResult);

  if (dMag < OPTMS_MACHINE_EPS) OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Division by zero")

  adResult[0] /= dMag;
  adResult[1] /= dMag;
  adResult[2] /= dMag;

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vSineDihedrals"
int  vSineDihedrals(const double adCoord0[3], 
                    const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
		    double adResult[6], int* const piNResult)
{
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  double adTemp[3];
  int ierr;

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  ierr = vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD); OPTMS_CHKERR(ierr);

  vCross(adNormA, adNormB, adTemp);
  adResult[5] = dMagnitude(adTemp); /* Edge 23 */

  vCross(adNormA, adNormC, adTemp);
  adResult[4] = dMagnitude(adTemp); /* Edge 13 */

  vCross(adNormA, adNormD, adTemp);
  adResult[3] = dMagnitude(adTemp); /* Edge 12 */

  vCross(adNormB, adNormC, adTemp);
  adResult[2] = dMagnitude(adTemp); /* Edge 03 */

  vCross(adNormB, adNormD, adTemp);
  adResult[1] = dMagnitude(adTemp); /* Edge 02 */

  vCross(adNormC, adNormD, adTemp);
  adResult[0] = dMagnitude(adTemp); /* Edge 01 */

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vCosineDihedrals"
int vCosineDihedrals(const double adCoord0[3], 
                      const double adCoord1[3],
		      const double adCoord2[3], const double adCoord3[3],
		      double adResult[6], int* const piNResult)
{
  int ierr;
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  ierr = vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD); OPTMS_CHKERR(ierr);

  adResult[5] = -dDot(adNormA, adNormB); /* Edge 23 */
  adResult[4] = -dDot(adNormA, adNormC); /* Edge 13 */
  adResult[3] = -dDot(adNormA, adNormD); /* Edge 12 */
  adResult[2] = -dDot(adNormB, adNormC); /* Edge 03 */
  adResult[1] = -dDot(adNormB, adNormD); /* Edge 02 */
  adResult[0] = -dDot(adNormC, adNormD); /* Edge 01 */
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vDihedrals"
int vDihedrals(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
  int ierr;
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  double dTemp;

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  ierr = vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);  OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC); OPTMS_CHKERR(ierr);
  ierr = vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD); OPTMS_CHKERR(ierr);

  assert(iFuzzyOne(dMagnitude(adNormA)));
  assert(iFuzzyOne(dMagnitude(adNormB)));
  assert(iFuzzyOne(dMagnitude(adNormC)));
  assert(iFuzzyOne(dMagnitude(adNormD)));

  dTemp = -dDot(adNormA, adNormB);
  adResult[5] = acos(dTemp);        /* Edge 23 */

  dTemp = -dDot(adNormA, adNormC);
  adResult[4] = acos(dTemp);        /* Edge 13 */

  dTemp = -dDot(adNormA, adNormD);
  adResult[3] = acos(dTemp);        /* Edge 12 */

  dTemp = -dDot(adNormB, adNormC);
  adResult[2] = acos(dTemp);        /* Edge 03 */

  dTemp = -dDot(adNormB, adNormD);
  adResult[1] = acos(dTemp);        /* Edge 02 */

  dTemp = -dDot(adNormC, adNormD);
  adResult[0] = acos(dTemp);        /* Edge 01 */
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vNegateDihedrals"
int vNegateDihedrals(const double adCoord0[3], 
                const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
     int ierr;
     int i;
     ierr = vDihedrals(adCoord0, adCoord1, adCoord2, adCoord3,
		adResult, piNResult); OPTMS_CHKERR(ierr);
     for (i=0;i<*piNResult;i++) adResult[i] = -adResult[i];
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vNegateCosineDihedrals"
int vNegateCosineDihedrals(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
     int ierr;
     int i;
     ierr = vCosineDihedrals(adCoord0, adCoord1, adCoord2, adCoord3,
		adResult, piNResult);  OPTMS_CHKERR(ierr);
     for (i=0;i<*piNResult;i++) adResult[i] = -adResult[i];
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vScaledJacobian"
int vScaledJacobian(const double adCoord0[3], 
                    const double adCoord1[3],
  		    const double adCoord2[3], 
                    const double adCoord3[3],
  		    double adResult[6], int* const piNResult)
{
    int ierr;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;
    double L1, L2, L3, L4, L5, L6;
    double L1_squared, L2_squared, L3_squared;
    double L4_squared, L5_squared, L6_squared;
    double jacobian;
    
    x0 = adCoord0[OPTMS_XDIR];  y0 = adCoord0[OPTMS_YDIR];  z0 = adCoord0[OPTMS_ZDIR];
    x1 = adCoord1[OPTMS_XDIR];  y1 = adCoord1[OPTMS_YDIR];  z1 = adCoord1[OPTMS_ZDIR];
    x2 = adCoord2[OPTMS_XDIR];  y2 = adCoord2[OPTMS_YDIR];  z2 = adCoord2[OPTMS_ZDIR];
    x3 = adCoord3[OPTMS_XDIR];  y3 = adCoord3[OPTMS_YDIR];  z3 = adCoord3[OPTMS_ZDIR];

    L1_squared = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    L2_squared = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0);
    L3_squared = (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0);
    L4_squared = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    L5_squared = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1);
    L6_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2);

    L1 = sqrt(L1_squared);            L4 = sqrt(L4_squared);
    L2 = sqrt(L2_squared);            L5 = sqrt(L5_squared);
    L3 = sqrt(L3_squared);            L6 = sqrt(L6_squared);

    if ((L1<OPTMS_MACHINE_EPS) || (L2<OPTMS_MACHINE_EPS) ||
        (L3<OPTMS_MACHINE_EPS) || (L4<OPTMS_MACHINE_EPS) ||
        (L5<OPTMS_MACHINE_EPS) || (L6<OPTMS_MACHINE_EPS)) {
          OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    }

    jacobian = dComputeTetJacobian(adCoord0,adCoord1,adCoord2,adCoord3);

    adResult[0] = jacobian/(L1*L2*L3);
    adResult[1] = jacobian/(L1*L4*L5);
    adResult[2] = jacobian/(L2*L4*L6);
    adResult[3] = jacobian/(L3*L5*L6);

    *piNResult = 4;
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vSMRSVolumeRatio"
int vSMRSVolumeRatio(const double adCoord0[3], 
                     const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
    int ierr;
    double jacobian, volume, srms;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;
    double L1_squared, L2_squared, L3_squared;
    double L4_squared, L5_squared, L6_squared;

    /* computes the ratio of (1/6 sum_i L_i^2) ^ (3/2) / (8.47*Volume)
       this needs to be negated as it ranges from 1 to infinity (equilateral
       triangle is 1, zero volume is infinity) and we therefore need to 
       minimize the maximum ratio */

    x0 = adCoord0[OPTMS_XDIR];  y0 = adCoord0[OPTMS_YDIR];  z0 = adCoord0[OPTMS_ZDIR];
    x1 = adCoord1[OPTMS_XDIR];  y1 = adCoord1[OPTMS_YDIR];  z1 = adCoord1[OPTMS_ZDIR];
    x2 = adCoord2[OPTMS_XDIR];  y2 = adCoord2[OPTMS_YDIR];  z2 = adCoord2[OPTMS_ZDIR];
    x3 = adCoord3[OPTMS_XDIR];  y3 = adCoord3[OPTMS_YDIR];  z3 = adCoord3[OPTMS_ZDIR];

    jacobian = dComputeTetJacobian(adCoord0,adCoord1,adCoord2,adCoord3);
    volume  = jacobian/6.;

/* DPS : changed below to use EPS cubed, since this Jacobian is not scaled. */
    if (fabs(jacobian)<(OPTMS_MACHINE_EPS)*(OPTMS_MACHINE_EPS)*(OPTMS_MACHINE_EPS) ) {

          OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Zero volume tetrahedron");
    }

    L1_squared = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    L2_squared = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0);
    L3_squared = (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0);
    L4_squared = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    L5_squared = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1);
    L6_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2);

    srms = sqrt(1./6.*(L1_squared+L2_squared +L3_squared+
                       L4_squared+L5_squared+L6_squared));

    adResult[0] = - pow(srms, 3) / (8.47967 * volume);
    *piNResult = 1;
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vComputeTetVolume"
int vComputeTetVolume(const double adCoord0[3], 
                       const double adCoord1[3],
		       const double adCoord2[3], const double adCoord3[3],
                       double adResult[6], int* const piNResult)
{
    int ierr;
    double jacobian, volume;

    jacobian = dComputeTetJacobian(adCoord0,adCoord1,adCoord2,adCoord3);
    volume  = jacobian/6.;
    adResult[0] = volume;
    *piNResult = 1;
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "dComputeTetJacobian"
double dComputeTetJacobian(const double adCoord0[3], const double adCoord1[3],
                           const double adCoord2[3], const double adCoord3[3])
{
    double a, b, c ,d;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;

    x0 = adCoord0[OPTMS_XDIR];  y0 = adCoord0[OPTMS_YDIR];  z0 = adCoord0[OPTMS_ZDIR];
    x1 = adCoord1[OPTMS_XDIR];  y1 = adCoord1[OPTMS_YDIR];  z1 = adCoord1[OPTMS_ZDIR];
    x2 = adCoord2[OPTMS_XDIR];  y2 = adCoord2[OPTMS_YDIR];  z2 = adCoord2[OPTMS_ZDIR];
    x3 = adCoord3[OPTMS_XDIR];  y3 = adCoord3[OPTMS_YDIR];  z3 = adCoord3[OPTMS_ZDIR];

    a = x1*(y2*z3-z2*y3) - x2*(y1*z3-y3*z1) + x3*(y1*z2 - z1*y2);
    b = (y2*z3-z2*y3) - (y1*z3-y3*z1) + (y1*z2-y2*z1);
    c = x1*(z3-z2) - x2*(z3-z1)+x3*(z2-z1);
    d = x1*(y2-y3) - x2*(y1-y3)+x3*(y1-y2);

    return (a-b*x0-c*y0-d*z0);
}
/*******************************************************************/
/*                  GRADIENT FUNCTIONS                             */
/*******************************************************************/

#undef __FUNC__
#define __FUNC__ "vGradSineDihedrals"
int vGradSineDihedrals(const double adCoord0[3], 
                       const double adCoord1[3],
                       const double adCoord2[3], const double adCoord3[3],
		       double **adGradient, int* const piNGradient)
{
  int ierr;
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);
  
  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);
  
  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */
  
  /* call the function */
  g_ad_vSineDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);
  
  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    adGradient[i11][0] = DERIV_GRAD(DTResult[i11])[0];
    adGradient[i11][1] = DERIV_GRAD(DTResult[i11])[1];
    adGradient[i11][2] = DERIV_GRAD(DTResult[i11])[2];
  }
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradCosineDihedrals"
int vGradCosineDihedrals(const double adCoord0[3], 
                                          const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double **adGradient, int* const piNGradient)
{
  int ierr;
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);
  
  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);
  
  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */
  
  /* call the function */
  g_ad_vCosineDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);
  
  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    adGradient[i11][0] = DERIV_GRAD(DTResult[i11])[0];
    adGradient[i11][1] = DERIV_GRAD(DTResult[i11])[1];
    adGradient[i11][2] = DERIV_GRAD(DTResult[i11])[2];
  }
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradDihedrals"
int vGradDihedrals(const double adCoord0[3], 
                                          const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double **adGradient, int* const piNGradient)
{
  int ierr;
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);
  
  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);
  
  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */
  
  /* call the function */
  g_ad_vDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);
  
  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    adGradient[i11][0] = DERIV_GRAD(DTResult[i11])[0];
    adGradient[i11][1] = DERIV_GRAD(DTResult[i11])[1];
    adGradient[i11][2] = DERIV_GRAD(DTResult[i11])[2];
  }
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vNegateGradDihedrals"
int vNegateGradDihedrals(const double adCoord0[3], 
                                          const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],  
			double **adGradient, int* const piNGradient)
{
     int ierr;
     int i,j;
     ierr = vGradDihedrals(adCoord0, adCoord1, adCoord2, adCoord3,
			adGradient, piNGradient); OPTMS_CHKERR(ierr);
     for(i=0;i<*piNGradient;i++) {
       for (j=0;j<3;j++) adGradient[i][j] = -adGradient[i][j];
     }
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vNegateGradCosineDihedrals"
int vNegateGradCosineDihedrals(const double adCoord0[3], 
                                          const double adCoord1[3],
			const double adCoord2[3], const double adCoord3[3],
			double **adGradient, int* const piNGradient)
{
     int ierr;
     int i,j;
     ierr = vGradCosineDihedrals(adCoord0, adCoord1, adCoord2, adCoord3,
			adGradient, piNGradient);  OPTMS_CHKERR(ierr);
     for(i=0;i<*piNGradient;i++) {
       for (j=0;j<3;j++) adGradient[i][j] = -adGradient[i][j];
     }
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradScaledJacobian"
int vGradScaledJacobian(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double **adGradient, int* const piNGradient)
{
    int ierr;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;
    double L1, L2, L3, L4, L5, L6;
    double L1_squared, L2_squared, L3_squared;
    double L4_squared, L5_squared, L6_squared;
    double L1_3_2, L2_3_2, L3_3_2;
    double a, b, c, d, J;
    double factor1, factor2, factor3, factor4;

    
    x0 = adCoord0[OPTMS_XDIR];  y0 = adCoord0[OPTMS_YDIR];  z0 = adCoord0[OPTMS_ZDIR];
    x1 = adCoord1[OPTMS_XDIR];  y1 = adCoord1[OPTMS_YDIR];  z1 = adCoord1[OPTMS_ZDIR];
    x2 = adCoord2[OPTMS_XDIR];  y2 = adCoord2[OPTMS_YDIR];  z2 = adCoord2[OPTMS_ZDIR];
    x3 = adCoord3[OPTMS_XDIR];  y3 = adCoord3[OPTMS_YDIR];  z3 = adCoord3[OPTMS_ZDIR];

    L1_squared = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    L2_squared = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0);
    L3_squared = (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0);
    L4_squared = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    L5_squared = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1);
    L6_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2);

    L1 = sqrt(L1_squared);            L4 = sqrt(L4_squared);
    L2 = sqrt(L2_squared);            L5 = sqrt(L5_squared);
    L3 = sqrt(L3_squared);            L6 = sqrt(L6_squared);

    if ((L1<OPTMS_MACHINE_EPS) || (L2<OPTMS_MACHINE_EPS) ||
        (L3<OPTMS_MACHINE_EPS) || (L4<OPTMS_MACHINE_EPS) ||
        (L5<OPTMS_MACHINE_EPS) || (L6<OPTMS_MACHINE_EPS)) {
          OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Colocated vertices");
    }

    L1_3_2 = L1*L1_squared;    
    L2_3_2 = L2*L2_squared;    
    L3_3_2 = L3*L3_squared;    
  
    a = x1*(y2*z3-z2*y3) - x2*(y1*z3-y3*z1) + x3*(y1*z2 - z1*y2);
    b = (y2*z3-z2*y3) - (y1*z3-y3*z1) + (y1*z2-y2*z1);
    c = x1*(z3-z2) - x2*(z3-z1)+x3*(z2-z1);
    d = x1*(y2-y3) - x2*(y1-y3)+x3*(y1-y2);

    J = a-b*x0-c*y0-d*z0;

    factor1 = 1/(L1*L2*L3);
    factor2 = J/(L1_3_2*L2*L3);
    factor3 = J/(L1*L2_3_2*L3);
    factor4 = J/(L1*L2*L3_3_2);

    adGradient[0][0] = -b*factor1 - (x0-x1)*factor2 - 
                                (x0-x2)*factor3 - (x0-x3)*factor4;
    adGradient[0][1] = -c*factor1 - (y0-y1)*factor2 - 
                                (y0-y2)*factor3 - (y0-y3)*factor4;
    adGradient[0][2] = -d*factor1 - (z0-z1)*factor2 - 
                                (z0-z2)*factor3 - (z0-z3)*factor4;

    factor1 = 1/(L1*L4*L5);
    adGradient[1][0] = -b*factor1 - J*(x0-x1)/(L1_3_2*L4*L5);
    adGradient[1][1] = -c*factor1 - J*(y0-y1)/(L1_3_2*L4*L5);
    adGradient[1][2] = -d*factor1 - J*(z0-z1)/(L1_3_2*L4*L5);

    factor1 = 1/(L2*L4*L6);
    adGradient[2][0] = -b*factor1 - J*(x0-x1)/(L2_3_2*L4*L6);
    adGradient[2][1] = -c*factor1 - J*(y0-y1)/(L2_3_2*L4*L6);
    adGradient[2][2] = -d*factor1 - J*(z0-z1)/(L2_3_2*L4*L6);

    factor1 = 1/(L3*L5*L6);
    adGradient[3][0] = -b*factor1 - J*(x0-x1)/(L3_3_2*L5*L6);
    adGradient[3][1] = -c*factor1 - J*(y0-y1)/(L3_3_2*L5*L6);
    adGradient[3][2] = -d*factor1 - J*(z0-z1)/(L3_3_2*L5*L6);
 
    *piNGradient = 4;
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradSMRSVolumeRatio"
int vGradSMRSVolumeRatio(const double adCoord0[3], 
                            const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double **adGradient, int* const piNGradient)
{
    int ierr;
    double jacobian, srms, srms3;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;
    double L1_squared, L2_squared, L3_squared;
    double L4_squared, L5_squared, L6_squared;
    double a, b, c, d;

    x0 = adCoord0[OPTMS_XDIR];  y0 = adCoord0[OPTMS_YDIR];  z0 = adCoord0[OPTMS_ZDIR];
    x1 = adCoord1[OPTMS_XDIR];  y1 = adCoord1[OPTMS_YDIR];  z1 = adCoord1[OPTMS_ZDIR];
    x2 = adCoord2[OPTMS_XDIR];  y2 = adCoord2[OPTMS_YDIR];  z2 = adCoord2[OPTMS_ZDIR];
    x3 = adCoord3[OPTMS_XDIR];  y3 = adCoord3[OPTMS_YDIR];  z3 = adCoord3[OPTMS_ZDIR];
    
    a = x1*(y2*z3-z2*y3) - x2*(y1*z3-y3*z1) + x3*(y1*z2 - z1*y2);
    b = (y2*z3-z2*y3) - (y1*z3-y3*z1) + (y1*z2-y2*z1);
    c = x1*(z3-z2) - x2*(z3-z1)+x3*(z2-z1);
    d = x1*(y2-y3) - x2*(y1-y3)+x3*(y1-y2);

    jacobian = a-b*x0-c*y0-d*z0;
/* DPS : changed below to use EPS cubed, since this Jacobian is not scaled. */
    if (fabs(jacobian)<(OPTMS_MACHINE_EPS)*(OPTMS_MACHINE_EPS)*(OPTMS_MACHINE_EPS) ) {
          OPTMS_SETERR(OPTMS_DIVIDE_BY_ZERO_ERR,0,"Zero volume tetrahedron");
    }

    L1_squared = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    L2_squared = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0);
    L3_squared = (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0);
    L4_squared = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    L5_squared = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1);
    L6_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2);

    srms = sqrt(1./6.*(L1_squared+L2_squared +L3_squared+
                              L4_squared+L5_squared+L6_squared));
    srms3 = pow(srms,3);

    adGradient[0][0] = -( 3./2.*srms*(-2*x1+6*x0-2*x2 -2*x3)/ (8.47967 * jacobian)
                             + 6*srms3*b/(8.47967*jacobian*jacobian));
    adGradient[0][1] = -(3./2.*srms*(-2*y1+6*y0-2*y2 -2*y3)/ (8.47967 * jacobian)
                             + 6*srms3*c/(8.47967*jacobian*jacobian));
    adGradient[0][2] = -(3./2.*srms*(-2*z1+6*z0-2*z2 -2*z3)/ (8.47967 * jacobian)
                             + 6*srms3*d/(8.47967*jacobian*jacobian));
    *piNGradient = 1;
     return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vNormJacSquared"
int vNormJacSquared(const double adCoord0[3], 
                const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
    int ierr;
    double a1[3], a2[3], a3[3];

    /* the a matrix contains the Jacobian */
    a1[0] = adCoord1[0] - adCoord0[0];
    a1[1] = adCoord1[1] - adCoord0[1];
    a1[2] = adCoord1[2] - adCoord0[2];

    a2[0] = adCoord2[0] - adCoord0[0];
    a2[1] = adCoord2[1] - adCoord0[1];
    a2[2] = adCoord2[2] - adCoord0[2];

    a3[0] = adCoord3[0] - adCoord0[0];
    a3[1] = adCoord3[1] - adCoord0[1];
    a3[2] = adCoord3[2] - adCoord0[2];

    ierr = SMfrobenius_norm_squared3x3(a1,a2,a3,&adResult[0]);
           OPTMS_CHKERR(ierr);
    adResult[0] *= -1.;
    
    *piNResult = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradNormJacSquared"
int vGradNormJacSquared(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
    		    double **adGradient, int* const piNGradient)
{
    int ierr;
    double a1[3], a2[3], a3[3];
    /* the a matrix contains the Jacobian */
    a1[0] = adCoord1[0] - adCoord0[0];
    a1[1] = adCoord1[1] - adCoord0[1];
    a1[2] = adCoord1[2] - adCoord0[2];

    a2[0] = adCoord2[0] - adCoord0[0];
    a2[1] = adCoord2[1] - adCoord0[1];
    a2[2] = adCoord2[2] - adCoord0[2];

    a3[0] = adCoord3[0] - adCoord0[0];
    a3[1] = adCoord3[1] - adCoord0[1];
    a3[2] = adCoord3[2] - adCoord0[2];

    adGradient[0][0] = -2*(a1[0] + a2[0] + a3[0]);    
    adGradient[0][1] = -2*(a1[1] + a2[1] + a3[1]);    
    adGradient[0][2] = -2*(a1[2] + a2[2] + a3[2]);    
    *piNGradient = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vCondition"
int vCondition(const double adCoord0[3], 
                const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
    int ierr;
    double norm_a2, norm_adj2;
    double a1[3], a2[3], a3[3];
    double t1[3], t2[3], t3[3];
    double alpha;


    /* create the Jacobian matrix */
    ierr = SMcreateWeightedJacobian(adCoord1,adCoord2,adCoord3,adCoord0,a1,a2,a3);
           OPTMS_CHKERR(ierr);
    /* compute it's determinant and frobenius norm squared */
    ierr = SMdeterminant3x3(a1, a2, a3, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3(a1,a2,a3,&norm_a2); OPTMS_CHKERR(ierr);
    /* compute it's adjoint and its frobenius norm squared */
    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adj2); 
           OPTMS_CHKERR(ierr);

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_SETERR(OPTMS_DIVIDE_BY_ZERO_ERR,0,"Zero volume tetrahedron");
    }

    adResult[0] = -sqrt( norm_a2 * norm_adj2 / (alpha*alpha));
    adResult[0]*=.3333333333333;

    *piNResult = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ " SMcreateJacobian"
int SMcreateJacobian(const double adCoord1[3], 
                    const double adCoord2[3],
		    const double adCoord3[3], const double adCoord0[3],
                    double a1[3], double a2[3], double a3[3])
{
    int ierr;
    /* the a matrix contains the Jacobian */
    a1[0] = adCoord1[0] - adCoord0[0];
    a1[1] = adCoord1[1] - adCoord0[1];
    a1[2] = adCoord1[2] - adCoord0[2];

    a2[0] = adCoord2[0] - adCoord0[0];
    a2[1] = adCoord2[1] - adCoord0[1];
    a2[2] = adCoord2[2] - adCoord0[2];

    a3[0] = adCoord3[0] - adCoord0[0];
    a3[1] = adCoord3[1] - adCoord0[1];
    a3[2] = adCoord3[2] - adCoord0[2];
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcreateWeightedJacobian"
int SMcreateWeightedJacobian(const double adCoord1[3], 
                    const double adCoord2[3],
		    const double adCoord3[3], const double adCoord0[3],
                    double a1[3], double a2[3], double a3[3])
{
    int ierr;
    double inv_root3 = 1./sqrt(3.0);
    double inv_root6 = 1./sqrt(6.0);
    double t1[3];
    double t2[3];
    double t3[3];

    /* the a matrix contains the Jacobian */
    t1[0] = adCoord1[0] - adCoord0[0];
    t1[1] = adCoord1[1] - adCoord0[1];
    t1[2] = adCoord1[2] - adCoord0[2];

    t2[0] = adCoord2[0] - adCoord0[0];
    t2[1] = adCoord2[1] - adCoord0[1];
    t2[2] = adCoord2[2] - adCoord0[2];

    t3[0] = adCoord3[0] - adCoord0[0];
    t3[1] = adCoord3[1] - adCoord0[1];
    t3[2] = adCoord3[2] - adCoord0[2];

    a1[0] = t1[0];
    a1[1] = t1[1];
    a1[2] = t1[2];

    a2[0] = inv_root3*(-t1[0]+2*t2[0]);
    a2[1] = inv_root3*(-t1[1]+2*t2[1]);
    a2[2] = inv_root3*(-t1[2]+2*t2[2]);

    a3[0] = inv_root6*(-t1[0]-t2[0]+3*t3[0]);
    a3[1] = inv_root6*(-t1[1]-t2[1]+3*t3[1]);
    a3[2] = inv_root6*(-t1[2]-t2[2]+3*t3[2]);
    return(ierr=0);

}

#undef __FUNC__
#define __FUNC__ "vGradCondition"
int vGradCondition(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
    		    double **adGradient, int* const piNGradient)
{
    int ierr;
    double norma;
    double alpha;
    double a1[3], a2[3], a3[3];
    double w1[3], w2[3], w3[3];
    double j1[3], j2[3], j3[3];
    double inv_root3 = 1./sqrt(3.0);
    double inv_root6 = 1./sqrt(6.0);
    double denom;
    double dnorma_dx, dnorma_dy, dnorma_dz;
    double dalpha_dx, dalpha_dy, dalpha_dz;
    double dadja_dx, dadja_dy, dadja_dz;
    double sqrt_norma, sqrt_norm_adja, norm_adj2;

    denom = inv_root3*inv_root6;

    ierr = SMcreateJacobian(adCoord1,adCoord2,adCoord3,adCoord0,j1,j2,j3);
           OPTMS_CHKERR(ierr);
    ierr = SMcreateWeightedJacobian(adCoord1,adCoord2,adCoord3,adCoord0,w1,w2,w3);
           OPTMS_CHKERR(ierr);
    ierr = SMdeterminant3x3(w1, w2, w3, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( w1, w2, w3, &norma ); 
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( w1, w2, w3, a1, a2, a3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norm_adj2); 
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adj2);

    /* we got the determinant right */
    dalpha_dx = 6.0*denom*(-j1[1]*j2[2] + j1[1]*j3[2] - j2[1]*j3[2]
                           +j3[1]*j2[2] + j1[2]*j2[1] - j3[1]*j1[2] );
    dalpha_dy = 6.0*denom*( j1[0]*j2[2] - j1[0]*j3[2] - j2[0]*j1[2] +
                            j2[0]*j3[2] - j3[0]*j2[2] + j1[2]*j3[0]);
    dalpha_dz = 6.0*denom*(-j1[0]*j2[1] + j2[0]*j1[1] - j1[1]*j3[0] +
                            j3[0]*j2[1] - j2[0]*j3[1] + j1[0]*j3[1]);

    /* norm is correct*/
    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (-w1[0]-inv_root3*w2[0]-inv_root6*w3[0])/sqrt_norma;
    dnorma_dy = (-w1[1]-inv_root3*w2[1]-inv_root6*w3[1])/sqrt_norma;
    dnorma_dz = (-w1[2]-inv_root3*w2[2]-inv_root6*w3[2])/sqrt_norma;

    /* now let's try the adjoint */
    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (a2[0]*(inv_root3*w3[2]-inv_root6*w2[2]) + 
                a2[1]*(-w3[2]+inv_root6*w1[2]) +
                a2[2]*(-inv_root3*w1[2]+w2[2]) +
                a3[0]*(-inv_root3*w3[1]+inv_root6*w2[1]) +
                a3[1]*(-inv_root6*w1[1]+w3[1]) +
                a3[2]*(-w2[1]+inv_root3*w1[1]))/sqrt_norm_adja;

    dadja_dy = (a1[0]*(-inv_root3*w3[2]+inv_root6*w2[2]) + 
                a1[1]*(w3[2]-inv_root6*w1[2]) +
                a1[2]*(inv_root3*w1[2]-w2[2]) +
                a3[0]*(inv_root3*w3[0]-inv_root6*w2[0]) +
                a3[1]*(inv_root6*w1[0]-w3[0]) +
                a3[2]*(w2[0]-inv_root3*w1[0]))/sqrt_norm_adja;

    dadja_dz = (a1[0]*(-inv_root6*w2[1]+inv_root3*w3[1]) + 
                a1[1]*(inv_root6*w1[1]-w3[1]) +
                a1[2]*(-inv_root3*w1[1]+w2[1]) +
                a2[0]*(-inv_root3*w3[0]+inv_root6*w2[0]) +
                a2[1]*(-inv_root6*w1[0]+w3[0]) +
                a2[2]*(-w2[0]+inv_root3*w1[0]))/sqrt_norm_adja;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    adGradient[0][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[0][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[0][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    adGradient[0][0] *= -.333333333;
    adGradient[0][1] *= -.333333333;
    adGradient[0][2] *= -.333333333;

    *piNGradient = 1;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "vGradConditionOld"
int vGradConditionOld(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
    		    double **adGradient, int* const piNGradient)
{
    int ierr;
    double norma, norm_adja;
    double alpha;
    double a1[3], a2[3], a3[3];
    double t1[3], t2[3], t3[3];
    double dnorma_dx, dnorma_dy, dnorma_dz;
    double dalpha_dx, dalpha_dy, dalpha_dz;
    double dadja_dx, dadja_dy, dadja_dz;
    double sqrt_norma, sqrt_norm_adja;


    /************************ VERTEX 1 ****************************/
    ierr = SMcreateWeightedJacobian(adCoord1,adCoord2,adCoord3,adCoord0,a1,a2,a3);
           OPTMS_CHKERR(ierr);

    ierr = SMdeterminant3x3(a1, a2, a3, &alpha); OPTMS_CHKERR(ierr);

    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norma ); 
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adja );
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adja);

    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (3*adCoord0[0] - adCoord1[0] - adCoord2[0] - adCoord3[0])/
                sqrt_norma;
    dnorma_dy = (3*adCoord0[1] - adCoord1[1] - adCoord2[1] - adCoord3[1])/
                sqrt_norma;
    dnorma_dz = (3*adCoord0[2] - adCoord1[2] - adCoord2[2] - adCoord3[2])/
                sqrt_norma;

    dalpha_dx = (adCoord2[1] - adCoord0[1])*(adCoord1[2] - adCoord3[2]) + 
                (adCoord3[1] - adCoord0[1])*(adCoord2[2] - adCoord1[2]) + 
                (adCoord1[1] - adCoord0[1])*(adCoord3[2] - adCoord2[2]);
    dalpha_dy = (adCoord1[0] - adCoord0[0])*(adCoord2[2] - adCoord3[2]) + 
                (adCoord2[0] - adCoord0[0])*(adCoord3[2] - adCoord1[2]) + 
                (adCoord3[0] - adCoord0[0])*(adCoord1[2] - adCoord2[2]);
    dalpha_dz = (adCoord1[0] - adCoord0[0])*(adCoord3[1] - adCoord2[1]) + 
                (adCoord2[0] - adCoord0[0])*(adCoord1[1] - adCoord3[1]) + 
                (adCoord3[0] - adCoord0[0])*(adCoord2[1] - adCoord1[1]);

    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (t2[0]*(adCoord3[2]-adCoord2[2]) + 
                t2[1]*(adCoord1[2]-adCoord3[2]) + 
                t2[2]*(adCoord2[2]-adCoord1[2]) + 
                t3[0]*(adCoord2[1]-adCoord3[1]) + 
                t3[1]*(adCoord3[1]-adCoord1[1]) + 
                t3[2]*(adCoord1[1]-adCoord2[1]))/sqrt_norm_adja;
    dadja_dy = (t1[0]*(adCoord2[2]-adCoord3[2]) + 
                t1[1]*(adCoord3[2]-adCoord1[2]) + 
                t1[2]*(adCoord1[2]-adCoord2[2]) + 
                t3[0]*(adCoord3[0]-adCoord2[0]) + 
                t3[1]*(adCoord1[0]-adCoord3[0]) + 
                t3[2]*(adCoord2[0]-adCoord1[0]))/sqrt_norm_adja;
    dadja_dz = (t1[0]*(adCoord3[1]-adCoord2[1]) + 
                t1[1]*(adCoord1[1]-adCoord3[1]) + 
                t1[2]*(adCoord2[1]-adCoord1[1]) + 
                t2[0]*(adCoord2[0]-adCoord3[0]) + 
                t2[1]*(adCoord3[0]-adCoord1[0]) + 
                t2[2]*(adCoord1[0]-adCoord2[0]))/sqrt_norm_adja;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    adGradient[0][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[0][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[0][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    adGradient[0][0] *= -.3333333333;
    adGradient[0][1] *= -.3333333333;
    adGradient[0][2] *= -.3333333333;

    /************************ VERTEX 2 ****************************/
    ierr = SMcreateJacobian(adCoord2,adCoord3,adCoord0,adCoord1,a1,a2,a3);
           OPTMS_CHKERR(ierr);
    ierr = SMdeterminant3x3(a1, a2, a3, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norma );
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adja );
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adja);

    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (adCoord0[0] - adCoord1[0])/sqrt_norma;
    dnorma_dy = (adCoord0[1] - adCoord1[1])/sqrt_norma;
    dnorma_dz = (adCoord0[2] - adCoord1[2])/sqrt_norma;

    dalpha_dx = (adCoord2[1] - adCoord1[1])*(adCoord3[2] - adCoord1[2]) -
                (adCoord3[1] - adCoord1[1])*(adCoord2[2] - adCoord1[2]);
    dalpha_dy = (adCoord1[0] - adCoord2[0])*(adCoord3[2] - adCoord1[2]) + 
                (adCoord3[0] - adCoord1[0])*(adCoord2[2] - adCoord1[2]);
    dalpha_dz = (adCoord2[0] - adCoord1[0])*(adCoord3[1] - adCoord1[1]) -
                (adCoord3[0] - adCoord1[0])*(adCoord2[1] - adCoord1[1]);

    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (t2[0]*(adCoord3[2]-adCoord1[2]) + 
                t2[1]*(adCoord1[2]-adCoord2[2]) + 
                t3[0]*(adCoord1[1]-adCoord3[1]) + 
                t3[1]*(adCoord2[1]-adCoord1[1]))/sqrt_norm_adja;
    dadja_dy = (t1[0]*(adCoord1[2]-adCoord3[2]) + 
                t1[1]*(adCoord2[2]-adCoord1[2]) + 
                t3[1]*(adCoord1[0]-adCoord2[0]) + 
                t3[0]*(adCoord3[0]-adCoord1[0]))/sqrt_norm_adja;
    dadja_dz = (t1[0]*(adCoord3[1]-adCoord1[1]) + 
                t1[1]*(adCoord1[1]-adCoord2[1]) + 
                t2[0]*(adCoord1[0]-adCoord3[0]) + 
                t2[1]*(adCoord2[0]-adCoord1[0]))/sqrt_norm_adja;;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    adGradient[1][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[1][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[1][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    
    /************************ VERTEX 2 ****************************/
    ierr = SMcreateJacobian(adCoord2,adCoord3,adCoord0,adCoord1,a1,a2,a3);
           OPTMS_CHKERR(ierr);
    ierr = SMdeterminant3x3(a1, a2, a3, &alpha);  OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norma );
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adja );
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adja);

    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (adCoord0[0] - adCoord1[0])/sqrt_norma;
    dnorma_dy = (adCoord0[1] - adCoord1[1])/sqrt_norma;
    dnorma_dz = (adCoord0[2] - adCoord1[2])/sqrt_norma;

    dalpha_dx = (adCoord2[1] - adCoord1[1])*(adCoord3[2] - adCoord1[2]) -
                (adCoord3[1] - adCoord1[1])*(adCoord2[2] - adCoord1[2]);
    dalpha_dy = (adCoord1[0] - adCoord2[0])*(adCoord3[2] - adCoord1[2]) + 
                (adCoord3[0] - adCoord1[0])*(adCoord2[2] - adCoord1[2]);
    dalpha_dz = (adCoord2[0] - adCoord1[0])*(adCoord3[1] - adCoord1[1]) -
                (adCoord3[0] - adCoord1[0])*(adCoord2[1] - adCoord1[1]);

    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (t2[0]*(adCoord3[2]-adCoord1[2]) + 
                t2[1]*(adCoord1[2]-adCoord2[2]) + 
                t3[0]*(adCoord1[1]-adCoord3[1]) + 
                t3[1]*(adCoord2[1]-adCoord1[1]))/sqrt_norm_adja;
    dadja_dy = (t1[0]*(adCoord1[2]-adCoord3[2]) + 
                t1[1]*(adCoord2[2]-adCoord1[2]) + 
                t3[1]*(adCoord1[0]-adCoord2[0]) + 
                t3[0]*(adCoord3[0]-adCoord1[0]))/sqrt_norm_adja;
    dadja_dz = (t1[0]*(adCoord3[1]-adCoord1[1]) + 
                t1[1]*(adCoord1[1]-adCoord2[1]) + 
                t2[0]*(adCoord1[0]-adCoord3[0]) + 
                t2[1]*(adCoord2[0]-adCoord1[0]))/sqrt_norm_adja;;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    adGradient[1][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[1][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[1][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    /*
    adGradient[1][0] *= -.3333333333;
    adGradient[1][1] *= -.3333333333;
    adGradient[1][2] *= -.3333333333;
    */
    /************************ VERTEX 3 ****************************/
    ierr = SMcreateJacobian(adCoord3,adCoord0,adCoord1,adCoord2,a1,a2,a3);
           OPTMS_CHKERR(ierr);
    ierr = SMdeterminant3x3(a1, a2, a3, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norma ); 
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adja); 
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adja);

    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (adCoord0[0] - adCoord2[0])/sqrt_norma;
    dnorma_dy = (adCoord0[1] - adCoord2[1])/sqrt_norma;
    dnorma_dz = (adCoord0[2] - adCoord2[2])/sqrt_norma;

    dalpha_dx = (adCoord1[1] - adCoord2[1])*(adCoord3[2] - adCoord2[2]) +
                (adCoord2[1] - adCoord3[1])*(adCoord1[2] - adCoord2[2]);
    dalpha_dy = (adCoord3[0] - adCoord2[0])*(adCoord1[2] - adCoord2[2]) + 
                (adCoord2[0] - adCoord1[0])*(adCoord3[2] - adCoord2[2]);
    dalpha_dz = (adCoord2[1] - adCoord1[1])*(adCoord3[0] - adCoord2[0]) +
                (adCoord3[1] - adCoord2[1])*(adCoord1[0] - adCoord2[0]);

    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (t2[0]*(adCoord2[2]-adCoord1[2]) + 
                t2[2]*(adCoord3[2]-adCoord2[2]) + 
                t3[0]*(adCoord1[1]-adCoord2[1]) + 
                t3[2]*(adCoord2[1]-adCoord3[1]))/sqrt_norm_adja;
    dadja_dy = (t1[0]*(adCoord1[2]-adCoord2[2]) + 
                t1[2]*(adCoord2[2]-adCoord3[2]) + 
                t3[0]*(adCoord2[0]-adCoord1[0]) + 
                t3[2]*(adCoord3[0]-adCoord2[0]))/sqrt_norm_adja;
    dadja_dz = (t1[0]*(adCoord2[1]-adCoord1[1]) + 
                t1[2]*(adCoord3[1]-adCoord2[1]) + 
                t2[0]*(adCoord1[0]-adCoord2[0]) + 
                t2[2]*(adCoord2[0]-adCoord3[0]))/sqrt_norm_adja;;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    adGradient[2][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[2][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[2][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    adGradient[2][0] *= -.3333333333;
    adGradient[2][1] *= -.3333333333;
    adGradient[2][2] *= -.3333333333;

    /************************ VERTEX 4 ****************************/
    ierr = SMcreateJacobian(adCoord0,adCoord1,adCoord2,adCoord3,a1,a2,a3);
           OPTMS_CHKERR(ierr);
    ierr = SMdeterminant3x3(a1, a2, a3, &alpha); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( a1, a2, a3, &norma); 
           OPTMS_CHKERR(ierr);
    sqrt_norma = sqrt(norma);

    ierr = SMadjoint3x3( a1, a2, a3, t1, t2, t3 ); OPTMS_CHKERR(ierr);
    ierr = SMfrobenius_norm_squared3x3( t1, t2, t3, &norm_adja); 
           OPTMS_CHKERR(ierr);
    sqrt_norm_adja = sqrt(norm_adja);

    if (sqrt_norma<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dnorma_dx = (adCoord0[0] - adCoord3[0])/sqrt_norma;
    dnorma_dy = (adCoord0[1] - adCoord3[1])/sqrt_norma;
    dnorma_dz = (adCoord0[2] - adCoord3[2])/sqrt_norma;

    dalpha_dx = (adCoord1[1] - adCoord3[1])*(adCoord2[2] - adCoord3[2]) -
                (adCoord2[1] - adCoord3[1])*(adCoord1[2] - adCoord3[2]);
    dalpha_dy = (adCoord3[0] - adCoord1[0])*(adCoord2[2] - adCoord3[2]) + 
                (adCoord2[0] - adCoord3[0])*(adCoord1[2] - adCoord3[2]);
    dalpha_dz = (adCoord1[0] - adCoord3[0])*(adCoord2[1] - adCoord3[1]) -
                (adCoord2[0] - adCoord3[0])*(adCoord1[1] - adCoord3[1]);

    if (sqrt_norm_adja<OPTMS_MACHINE_EPS) {
          OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);
    }
    dadja_dx = (t2[1]*(adCoord2[2]-adCoord3[2]) + 
                t2[2]*(adCoord3[2]-adCoord1[2]) + 
                t3[1]*(adCoord3[1]-adCoord2[1]) + 
                t3[2]*(adCoord1[1]-adCoord3[1]))/sqrt_norm_adja;
    dadja_dy = (t1[1]*(adCoord3[2]-adCoord2[2]) + 
                t1[2]*(adCoord1[2]-adCoord3[2]) + 
                t3[1]*(adCoord2[0]-adCoord3[0]) + 
                t3[2]*(adCoord3[0]-adCoord1[0]))/sqrt_norm_adja;
    dadja_dz = (t1[1]*(adCoord2[1]-adCoord3[1]) + 
                t1[2]*(adCoord3[1]-adCoord1[1]) + 
                t2[1]*(adCoord3[0]-adCoord2[0]) + 
                t2[2]*(adCoord1[0]-adCoord3[0]))/sqrt_norm_adja;;

    if (fabs(alpha)<OPTMS_MACHINE_EPS) {OPTMS_CHKERR(OPTMS_DIVIDE_BY_ZERO_ERR);}

    adGradient[3][0] = sqrt_norm_adja*dnorma_dx/alpha + 
                       sqrt_norma*dadja_dx/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dx/(alpha*alpha);
    adGradient[3][1] = sqrt_norm_adja*dnorma_dy/alpha + 
                       sqrt_norma*dadja_dy/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dy/(alpha*alpha);
    adGradient[3][2] = sqrt_norm_adja*dnorma_dz/alpha + 
                       sqrt_norma*dadja_dz/alpha - 
                       sqrt_norma*sqrt_norm_adja*dalpha_dz/(alpha*alpha);

    /*
    adGradient[3][0] *= -.33333333333;
    adGradient[3][1] *= -.33333333333;
    adGradient[3][2] *= -.33333333333;
    */
    *piNGradient = 4;
    return(ierr=0);

}

