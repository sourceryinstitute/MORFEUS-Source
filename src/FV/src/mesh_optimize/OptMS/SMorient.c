/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

/* 3D orientation predicate, without exact arithmetic */
#undef __FUNC__
#define __FUNC__ "Orient3D" 
int iOrient3D(const double adLocA[3], const double adLocB[3],
	      const double adLocC[3], const double adLocD[3]) 
     /* Returns 1 if tet ABCD is right-handed, -1 if it's left-handed,
	and 0 if it's essentially a tie. */
{
  double dXa = adLocA[0];
  double dYa = adLocA[1];
  double dZa = adLocA[2];

  double dXb = adLocB[0];
  double dYb = adLocB[1];
  double dZb = adLocB[2];

  double dXc = adLocC[0];
  double dYc = adLocC[1];
  double dZc = adLocC[2];

  double dXd = adLocD[0];
  double dYd = adLocD[1];
  double dZd = adLocD[2];

  double dDX2 = dXb - dXa;
  double dDX3 = dXc - dXa;
  double dDX4 = dXd - dXa;

  double dDY2 = dYb - dYa;
  double dDY3 = dYc - dYa;
  double dDY4 = dYd - dYa;

  double dDZ2 = dZb - dZa;
  double dDZ3 = dZc - dZa;
  double dDZ4 = dZd - dZa;

  /* dDet is proportional to the cell volume */
  double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
    - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;

  /* Compute a length scale based on edge lengths. */
  double dScale = ( sqrt((dXa-dXb)*(dXa-dXb) + (dYa-dYb)*(dYa-dYb) +
			 (dZa-dZb)*(dZa-dZb)) +
		    sqrt((dXa-dXc)*(dXa-dXc) + (dYa-dYc)*(dYa-dYc) +
			 (dZa-dZc)*(dZa-dZc)) +
		    sqrt((dXa-dXd)*(dXa-dXd) + (dYa-dYd)*(dYa-dYd) +
			 (dZa-dZd)*(dZa-dZd)) +
		    sqrt((dXb-dXc)*(dXb-dXc) + (dYb-dYc)*(dYb-dYc) +
			 (dZb-dZc)*(dZb-dZc)) +
		    sqrt((dXb-dXd)*(dXb-dXd) + (dYb-dYd)*(dYb-dYd) +
			 (dZb-dZd)*(dZb-dZd)) +
		    sqrt((dXc-dXd)*(dXc-dXd) + (dYc-dYd)*(dYc-dYd) +
			 (dZc-dZd)*(dZc-dZd)) ) / 6.;

  double dEps = 1.e-13;

  /* Use the length scale to get a better idea if the tet is flat or
     just really small. */
  dDet /= (dScale*dScale*dScale);

  if (dDet > dEps)
    return (1);
  else if (dDet < -dEps) 
    return (-1);
  else
    return (0);
}

