/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include "SMderiv.h"
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "SMintrinsic.h"
#include "SMdihed_func.h"

#undef __FUNC__
#define __FUNC__ "g_ad_vCross" 
int  g_ad_vCross(DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3], DERIV_TYPE adResult[3])
{
  int ierr;
  register double *GA0, *GA1, *GA2, *GB0, *GB1, *GB2;
  register double VA0, VA1, VA2, VB0, VB1, VB2;
  register double *GR0, *GR1, *GR2;
  int g_i_;

  GA0 = DERIV_GRAD(adVecA[0]);
  GA1 = DERIV_GRAD(adVecA[1]);
  GA2 = DERIV_GRAD(adVecA[2]);
  GB0 = DERIV_GRAD(adVecB[0]);
  GB1 = DERIV_GRAD(adVecB[1]);
  GB2 = DERIV_GRAD(adVecB[2]);

  VA0 = DERIV_VAL(adVecA[0]);
  VA1 = DERIV_VAL(adVecA[1]);
  VA2 = DERIV_VAL(adVecA[2]);
  VB0 = DERIV_VAL(adVecB[0]);
  VB1 = DERIV_VAL(adVecB[1]);
  VB2 = DERIV_VAL(adVecB[2]);

  GR0 = DERIV_GRAD(adResult[0]);
  GR1 = DERIV_GRAD(adResult[1]);
  GR2 = DERIV_GRAD(adResult[2]);
  for (g_i_ = 2; g_i_ >= 0; g_i_--) {
    GR0[g_i_] = VB2 * GA1[g_i_] + VA1 * GB2[g_i_]
      - VB1 * GA2[g_i_] - VA2 * GB1[g_i_];
    GR1[g_i_] = VB0 * GA2[g_i_] + VA2 * GB0[g_i_]
      - VB2 * GA0[g_i_] - VA0 * GB2[g_i_];
    GR2[g_i_] = VB1 * GA0[g_i_] + VA0 * GB1[g_i_]
      - VB0 * GA1[g_i_] - VA1 * GB0[g_i_];
  }

  DERIV_VAL(adResult[0]) = VA1 * VB2 - VA2 * VB1;
  DERIV_VAL(adResult[1]) = VA2 * VB0 - VA0 * VB2;
  DERIV_VAL(adResult[2]) = VA0 * VB1 - VA1 * VB0;
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_dMagnitude" 
int  g_ad_dMagnitude(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVec[3]) 
{
  int ierr;
  int g_i_;
  DERIV_TYPE g_ad_var_0;
  DERIV_TYPE g_ad_var_1;

  for (g_i_ = 0; g_i_ < PMAX; g_i_++) {
    DERIV_GRAD(g_ad_var_0)[g_i_] =
      2 * DERIV_VAL(adVec[0]) * DERIV_GRAD(adVec[0])[g_i_] +
      2 * DERIV_VAL(adVec[1]) * DERIV_GRAD(adVec[1])[g_i_] +
      2 * DERIV_VAL(adVec[2]) * DERIV_GRAD(adVec[2])[g_i_];
  }

  DERIV_VAL(g_ad_var_0) =
    DERIV_VAL(adVec[0]) * DERIV_VAL(adVec[0]) +
    DERIV_VAL(adVec[1]) * DERIV_VAL(adVec[1]) +
    DERIV_VAL(adVec[2]) * DERIV_VAL(adVec[2]);

  g_ad_sqrt( &g_ad_var_1, g_ad_var_0);

  for (g_i_ = 0; g_i_ < PMAX; g_i_++) {
    DERIV_GRAD( *g_ad_var_)[g_i_] = DERIV_GRAD(g_ad_var_1)[g_i_];
  }

  DERIV_VAL( *g_ad_var_) = DERIV_VAL(g_ad_var_1);
  return (ierr=0);  ;

}

#undef __FUNC__
#define __FUNC__ "g_ad_dDot" 
int  g_ad_dDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3]) 
{
  int ierr;
  double g_ad_adj_3;
  double g_ad_adj_4;
  double g_ad_adj_6;
  double g_ad_adj_7;
  double g_ad_adj_9;
  double g_ad_adj_10;
  int g_i_;
  g_ad_adj_10 = DERIV_VAL(adVecA[2]);
  g_ad_adj_9 = DERIV_VAL(adVecB[2]);
  g_ad_adj_7 = DERIV_VAL(adVecA[1]);
  g_ad_adj_6 = DERIV_VAL(adVecB[1]);
  g_ad_adj_4 = DERIV_VAL(adVecA[0]);
  g_ad_adj_3 = DERIV_VAL(adVecB[0]);

  for (g_i_ = 0; g_i_ < PMAX; g_i_++) {
    DERIV_GRAD( *g_ad_var_)[g_i_] =
      + g_ad_adj_3 * DERIV_GRAD(adVecA[0])[g_i_]
      +  g_ad_adj_4 * DERIV_GRAD(adVecB[0])[g_i_]
      +  g_ad_adj_6 * DERIV_GRAD(adVecA[1])[g_i_]
      +  g_ad_adj_7 * DERIV_GRAD(adVecB[1])[g_i_]
      +  g_ad_adj_9 * DERIV_GRAD(adVecA[2])[g_i_]
      +  g_ad_adj_10 * DERIV_GRAD(adVecB[2])[g_i_];
  }

  DERIV_VAL( *g_ad_var_) =
    DERIV_VAL(adVecA[0]) * DERIV_VAL(adVecB[0]) +
    DERIV_VAL(adVecA[1]) * DERIV_VAL(adVecB[1]) +
    DERIV_VAL(adVecA[2]) * DERIV_VAL(adVecB[2]);

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_dNegDot" 
int  g_ad_dNegDot(DERIV_TYPE  *g_ad_var_, DERIV_TYPE adVecA[3], DERIV_TYPE adVecB[3])
{
  int ierr;
  double g_ad_adj_3;
  double g_ad_adj_4;
  double g_ad_adj_6;
  double g_ad_adj_7;
  double g_ad_adj_9;
  double g_ad_adj_10;
  int g_i_;
  g_ad_adj_10 = DERIV_VAL(adVecA[2]);
  g_ad_adj_9 = DERIV_VAL(adVecB[2]);
  g_ad_adj_7 = DERIV_VAL(adVecA[1]);
  g_ad_adj_6 = DERIV_VAL(adVecB[1]);
  g_ad_adj_4 = DERIV_VAL(adVecA[0]);
  g_ad_adj_3 = DERIV_VAL(adVecB[0]);

  for (g_i_ = 0; g_i_ < PMAX; g_i_++) {
    DERIV_GRAD( *g_ad_var_)[g_i_] = -(
      + g_ad_adj_3 * DERIV_GRAD(adVecA[0])[g_i_]
      +  g_ad_adj_4 * DERIV_GRAD(adVecB[0])[g_i_]
      +  g_ad_adj_6 * DERIV_GRAD(adVecA[1])[g_i_]
      +  g_ad_adj_7 * DERIV_GRAD(adVecB[1])[g_i_]
      +  g_ad_adj_9 * DERIV_GRAD(adVecA[2])[g_i_]
      +  g_ad_adj_10 * DERIV_GRAD(adVecB[2])[g_i_]);
  }

  DERIV_VAL( *g_ad_var_) = -(
    DERIV_VAL(adVecA[0]) * DERIV_VAL(adVecB[0]) +
    DERIV_VAL(adVecA[1]) * DERIV_VAL(adVecB[1]) +
    DERIV_VAL(adVecA[2]) * DERIV_VAL(adVecB[2]));

  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_vUnitNormal" 
int  g_ad_vUnitNormal(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                       DERIV_TYPE adCoord2[3], DERIV_TYPE adResult[3]) 
{
  int ierr;
  DERIV_TYPE adVecA[3], adVecB[3], dMag;
  double g_ad_adj_1, g_ad_adj_20, g_ad_adj_21, g_ad_adj_22;
  int i, g_i_;
  for (i = 0; i < 3; i++) {

    DERIV_GRAD(adVecA[i])[0] =
      DERIV_GRAD(adCoord1[i])[0] - DERIV_GRAD(adCoord0[i])[0];
    DERIV_GRAD(adVecB[i])[0] =
      DERIV_GRAD(adCoord2[i])[0] - DERIV_GRAD(adCoord0[i])[0];

    DERIV_GRAD(adVecA[i])[1] =
      DERIV_GRAD(adCoord1[i])[1] - DERIV_GRAD(adCoord0[i])[1];
    DERIV_GRAD(adVecB[i])[1] =
      DERIV_GRAD(adCoord2[i])[1] - DERIV_GRAD(adCoord0[i])[1];

    DERIV_GRAD(adVecA[i])[2] =
      DERIV_GRAD(adCoord1[i])[2] - DERIV_GRAD(adCoord0[i])[2];
    DERIV_GRAD(adVecB[i])[2] =
      DERIV_GRAD(adCoord2[i])[2] - DERIV_GRAD(adCoord0[i])[2];

    DERIV_VAL(adVecA[i]) = DERIV_VAL(adCoord1[i]) - DERIV_VAL(adCoord0[i]);
    DERIV_VAL(adVecB[i]) = DERIV_VAL(adCoord2[i]) - DERIV_VAL(adCoord0[i]);
     
  }
  g_ad_vCross(adVecA, adVecB, adResult);
  g_ad_dMagnitude( &dMag, adResult);
  g_ad_adj_1 = 1./DERIV_VAL(dMag);
  DERIV_VAL(adResult[0]) *= g_ad_adj_1;
  DERIV_VAL(adResult[1]) *= g_ad_adj_1;
  DERIV_VAL(adResult[2]) *= g_ad_adj_1;
  g_ad_adj_20 = - DERIV_VAL(adResult[0]) * g_ad_adj_1;
  g_ad_adj_21 = - DERIV_VAL(adResult[1]) * g_ad_adj_1;
  g_ad_adj_22 = - DERIV_VAL(adResult[2]) * g_ad_adj_1;

  for (g_i_ = 0; g_i_ < PMAX; g_i_++) {
    DERIV_GRAD(adResult[0])[g_i_] =
      g_ad_adj_1 * DERIV_GRAD(adResult[0])[g_i_] +
      g_ad_adj_20 * DERIV_GRAD(dMag)[g_i_];
    DERIV_GRAD(adResult[1])[g_i_] =
      g_ad_adj_1 * DERIV_GRAD(adResult[1])[g_i_] +
      g_ad_adj_21 * DERIV_GRAD(dMag)[g_i_];
    DERIV_GRAD(adResult[2])[g_i_] =
      g_ad_adj_1 * DERIV_GRAD(adResult[2])[g_i_] +
      g_ad_adj_22 * DERIV_GRAD(dMag)[g_i_];
  }
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_vSineDihedrals" 
int g_ad_vSineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                         DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                         DERIV_TYPE adResult[6], int  *piNResult) 
{
  int ierr;
  DERIV_TYPE adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  DERIV_TYPE adTemp[3];
  *piNResult = 6;
  g_ad_vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  g_ad_vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  g_ad_vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  g_ad_vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);
  g_ad_vCross(adNormA, adNormB, adTemp);
  g_ad_dMagnitude( &adResult[5], adTemp);
  g_ad_vCross(adNormA, adNormC, adTemp);
  g_ad_dMagnitude( &adResult[4], adTemp);
  g_ad_vCross(adNormA, adNormD, adTemp);
  g_ad_dMagnitude( &adResult[3], adTemp);
  g_ad_vCross(adNormB, adNormC, adTemp);
  g_ad_dMagnitude( &adResult[2], adTemp);
  g_ad_vCross(adNormB, adNormD, adTemp);
  g_ad_dMagnitude( &adResult[1], adTemp);
  g_ad_vCross(adNormC, adNormD, adTemp);
  g_ad_dMagnitude( &adResult[0], adTemp);
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_vCosineDihedrals" 
int  g_ad_vCosineDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                                   DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                                   DERIV_TYPE adResult[6], int  *piNResult) 
{
  int ierr;
  DERIV_TYPE adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  *piNResult = 6;
  g_ad_vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  g_ad_vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  g_ad_vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  g_ad_vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);

  g_ad_dNegDot( &adResult[5], adNormA, adNormB);
  g_ad_dNegDot( &adResult[4], adNormA, adNormC);
  g_ad_dNegDot( &adResult[3], adNormA, adNormD);
  g_ad_dNegDot( &adResult[2], adNormB, adNormC);
  g_ad_dNegDot( &adResult[1], adNormB, adNormD);
  g_ad_dNegDot( &adResult[0], adNormC, adNormD);
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "g_ad_vDihedrals" 
int g_ad_vDihedrals(DERIV_TYPE adCoord0[3], DERIV_TYPE adCoord1[3], 
                             DERIV_TYPE adCoord2[3], DERIV_TYPE adCoord3[3], 
                             DERIV_TYPE adResult[6], int  *piNResult) 
{
  int ierr;
  DERIV_TYPE adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  DERIV_TYPE dTemp;
  *piNResult = 6;
  g_ad_vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  g_ad_vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  g_ad_vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  g_ad_vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);
  g_ad_dNegDot( &dTemp, adNormA, adNormB);
  g_ad_acos( &adResult[5], dTemp);
  g_ad_dNegDot( &dTemp, adNormA, adNormC);
  g_ad_acos( &adResult[4], dTemp);
  g_ad_dNegDot( &dTemp, adNormA, adNormD);
  g_ad_acos( &adResult[3], dTemp);
  g_ad_dNegDot( &dTemp, adNormB, adNormC);
  g_ad_acos( &adResult[2], dTemp);
  g_ad_dNegDot( &dTemp, adNormB, adNormD);
  g_ad_acos( &adResult[1], dTemp);
  g_ad_dNegDot( &dTemp, adNormC, adNormD);
  g_ad_acos( &adResult[0], dTemp);
  return(ierr=0);
}


