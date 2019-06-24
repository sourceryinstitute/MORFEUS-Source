/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_DIHED_FUNC_H
#define SM_DIHED_FUNC_H 1

void vCross(const double adVecA[3], const double adVecB[3],
		   double adResult[3]);
double dMagnitude(const double adVec[3]);
double dDot(const double adVecA[3], const double adVecB[3]);

int vUnitNormal(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], double adResult[3]);
double dComputeTetJacobian(const double adCoord0[3], const double adCoord1[3],
                           const double adCoord2[3], const double adCoord3[3]);

int vDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vNegateDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vSineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vNegateCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vGradSineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vGradCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vNegateGradCosineDihedrals(const double adCoord0[3], 
               const double adCoord1[3], const double adCoord2[3], 
               const double adCoord3[3], double **adGradient, 
               int* const piNGradient);
int vGradDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vNegateGradDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vScaledJacobian(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vGradScaledJacobian(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vSMRSVolumeRatio(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vGradSMRSVolumeRatio(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vComputeTetVolume(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vNormJacSquared(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int vGradNormJacSquared(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double **adGradient, int* const piNGradient);
int vCondition(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult);
int SMcreateWeightedJacobian(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
               double a1[3], double a2[3], double a3[3]);
int SMcreateJacobian(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
               double a1[3], double a2[3], double a3[3]);
int vGradCondition(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
    	       double **adGradient, int* const piNGradient);
int vGradConditionOld(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
    	       double **adGradient, int* const piNGradient);


int  g_ad_vCross(DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3], DERIV_TYPE adResult[3]);
int  g_ad_dMagnitude(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVec[3]);
int  g_ad_dDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3]);
int  g_ad_dNegDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3]);
int  g_ad_vUnitNormal(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                      DERIV_TYPE adCoord2[3], DERIV_TYPE adResult[3]);
int  g_ad_vSineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                       DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                       DERIV_TYPE adResult[6], int  *piNResult);
int  g_ad_vCosineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                       DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                       DERIV_TYPE adResult[6], int  *piNResult);
int  g_ad_vDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                       DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                       DERIV_TYPE adResult[6], int  *piNResult);

#endif


