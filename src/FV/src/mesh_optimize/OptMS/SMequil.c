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

#define OPTMS_CCW               1
#define OPTMS_CW                0
#define OPTMS_NO_EQUIL 101
#define OPTMS_CHECK_TOP_DOWN 102
#define OPTMS_CHECK_BOTTOM_UP 103
#define OPTMS_TWO_PT_PLANE 104
#define OPTMS_THREE_PT_PLANE 105
#define OPTMS_CHECK_Y_COORD_DIRECTION 106
#define OPTMS_CHECK_X_COORD_DIRECTION 107
#define OPTMS_CHECK_Z_COORD_DIRECTION 108
#define OPTMS_EQUIL 109
#define OPTMS_HULL_TEST_ERROR 110

#undef __FUNC__
#define __FUNC__ "SMconvexHullTest" 
int SMconvexHullTest(double **vec, int num_vec, int *equil)
{
    int ierr;
    int status, dir_done;
    double pt1[3], pt2[3], pt3[3];
    double normal[3];

    /* tries to determine equilibrium for the 3D case */
    *equil = OPTMS_FALSE;
    status = OPTMS_CHECK_Z_COORD_DIRECTION;
    dir_done = -1;

    if (num_vec <= 2) status = OPTMS_NO_EQUIL;

    while ((status != OPTMS_EQUIL) && (status != OPTMS_NO_EQUIL) && 
           (status != OPTMS_HULL_TEST_ERROR)) {
       if (status == OPTMS_CHECK_Z_COORD_DIRECTION) {
          ierr = SMfindPlanePoints(OPTMS_ZDIR, OPTMS_YDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status);
                 OPTMS_CHKERR(ierr);
          dir_done = 2;
       }
       if (status == OPTMS_CHECK_Y_COORD_DIRECTION) {
          ierr = SMfindPlanePoints(OPTMS_YDIR, OPTMS_XDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status);
                 OPTMS_CHKERR(ierr);
          dir_done = 1;
       }
       if (status == OPTMS_CHECK_X_COORD_DIRECTION) {
          ierr = SMfindPlanePoints(OPTMS_XDIR, OPTMS_ZDIR, 
                          vec, num_vec, pt1, pt2, pt3, &status);
                 OPTMS_CHKERR(ierr);
          dir_done = 0;
       }
       if (status == OPTMS_TWO_PT_PLANE) {
          pt3[0]=0.; pt3[1]=0.; pt3[2]=0.;
       }
       if ((status == OPTMS_TWO_PT_PLANE) || (status == OPTMS_THREE_PT_PLANE)){
           ierr = SMfindPlaneNormal(pt1,pt2,pt3,normal); OPTMS_CHKERR(ierr);
           ierr = SMcheckVectorDots(vec,num_vec,normal,equil); 
                  OPTMS_CHKERR(ierr);
           if (*equil == 1) {
             switch(dir_done){
             case 2:
               *equil = 0; status = OPTMS_CHECK_Y_COORD_DIRECTION;
               break;
             case 1:
               *equil = 0; status = OPTMS_CHECK_X_COORD_DIRECTION;
               break;
             case 0:
               *equil = 1; status = OPTMS_EQUIL;
             }
           } else if (*equil == 0) {
               status = OPTMS_NO_EQUIL;
           } else {
               OPTMS_SETERR(OPTMS_DATA_ERR,0,"equil flag not set to 0 or 1");
           }
       }
    }
    switch (status){
    case OPTMS_NO_EQUIL:
      OPTMS_DEBUG_PRINT(3,"Not an equilibrium point\n");
      *equil = 0; 
      break;
    case OPTMS_EQUIL:
      OPTMS_DEBUG_PRINT(3,"An equilibrium point\n");
      *equil = 1;
      break;
    default:
      OPTMS_DEBUG_ACTION(3,{
        fprintf(stdout,"Failed to determine equil or not; status = %d\n",status);
      });
    }
    return (ierr=0);
}   

#undef __FUNC__
#define __FUNC__ "SMcheckVectorDots" 
int SMcheckVectorDots(double **vec,int num_vec,double *normal, int *equil)
{
    int ierr;
    int i, ind;
    double test_dot, dot;

    *equil = OPTMS_FALSE;
    OPTMS_DOT(test_dot,vec[0],normal,3);
    ind = 1;
    while ((fabs(test_dot) < OPTMS_MACHINE_EPS) && (ind<num_vec)) {
      OPTMS_DOT(test_dot,vec[ind],normal,3);
      ind++;
    }
      
    for (i=ind;i<num_vec;i++) {
       OPTMS_DOT(dot,vec[i],normal,3);
       if ( ((dot>0 && test_dot<0) || (dot<0 && test_dot>0)) &&
            (fabs(dot)>OPTMS_MACHINE_EPS)) {
          *equil = OPTMS_TRUE;
          return(ierr=0);
       }
    }
    return (ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMfindPlaneNormal" 
int SMfindPlaneNormal(double pt1[3], double pt2[3], double pt3[3],
                          double *cross)
{
  int ierr;
  int i;
  double vecA[3], vecB[3];

  for (i=0;i<3;i++) {
    vecA[i] = pt2[i] - pt1[i];
    vecB[i] = pt3[i] - pt1[i];
  }
  cross[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  cross[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  cross[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
  OPTMS_NORMALIZE(cross, 3);
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfindPlanePoints" 
int SMfindPlanePoints(int dir1, int dir2, double **vec, int num_vec,
                      double *pt1, double *pt2, double *pt3, 
                      int *status)
{
    int ierr;
    int i;
    int ind[50], num_min, num_max;
    int rotate=OPTMS_CW;
    int num_rotated=0;
    double pt_1, pt_2;
    double min, inv_slope;
    double min_inv_slope=0.;
    double max; 
    double max_inv_slope=0;
    double inv_origin_slope=0;

    *status = OPTMS_CHECK_BOTTOM_UP;
    /* find the minimum points in dir1 starting at -1 */
    num_min = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; min=1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1]<min) {
	min = vec[i][dir1]; ind[0] = i; num_min = 1;
      } else if (fabs(vec[i][dir1] - min) < OPTMS_MACHINE_EPS) {
	ind[num_min++] = i;
      }
    }
    if (min >= 0) *status = OPTMS_NO_EQUIL;
 
    if (*status != OPTMS_NO_EQUIL) {
      switch(num_min) {
      case 1: /* rotate to find the next point */
	OPTMS_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1]; pt_2 = pt1[dir2];
	if (pt1[dir2] <= 0){rotate=OPTMS_CCW; max_inv_slope=OPTMS_BIG_NEG_NMBR;}
	if (pt1[dir2] > 0){rotate=OPTMS_CW; min_inv_slope=OPTMS_BIG_POS_NMBR;}
	switch(rotate) {
	case OPTMS_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope>max_inv_slope) &&  
		  (fabs(inv_slope - max_inv_slope) > OPTMS_MACHINE_EPS)) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < OPTMS_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case OPTMS_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope<min_inv_slope) && 
		  (fabs(inv_slope - max_inv_slope) > OPTMS_MACHINE_EPS)){
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < OPTMS_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  OPTMS_DEBUG_PRINT(3,"No points in the rotation ... odd\n");
	    *status = OPTMS_HULL_TEST_ERROR;
	  break;
	case 1:
	  OPTMS_DEBUG_PRINT(3,"Found a line in the convex hull\n");
	  OPTMS_COPY_VECTOR(pt2,vec[ind[1]],3); *status = OPTMS_TWO_PT_PLANE;
	  break;
	default:
	  OPTMS_DEBUG_PRINT(3,"Found 2 or more points in the rotation\n");
	    if (fabs(pt_1) > OPTMS_MACHINE_EPS) inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case OPTMS_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=OPTMS_NO_EQUIL;
	    else *status=OPTMS_CHECK_TOP_DOWN;
	    break;
	  case OPTMS_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=OPTMS_NO_EQUIL;
	    else *status=OPTMS_CHECK_TOP_DOWN;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	OPTMS_DEBUG_PRINT(3,"Found two minimum points to define the plane\n");
                OPTMS_COPY_VECTOR(pt1,vec[ind[0]],3);
	OPTMS_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = OPTMS_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	OPTMS_DEBUG_ACTION(3,{fprintf(stdout,"Found 3 or more points in min plane %f\n",min);})
	  if (vec[ind[0]][dir1] >= 0) *status = OPTMS_NO_EQUIL;
	  else *status = OPTMS_CHECK_TOP_DOWN;
    }
    }

    /***************************/
    /*  failed to find any information, checking top/down this coord*/
    /***************************/

    if (*status == OPTMS_CHECK_TOP_DOWN) {
    /* find the maximum points in dir1 starting at 1 */
    num_max = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; max=-1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1] > max) {
	max = vec[i][dir1]; ind[0] = i; num_max = 1;
      } else if (fabs(vec[i][dir1] - max) < OPTMS_MACHINE_EPS) {
	ind[num_max++] = i;
      }
    }
    if (max <= 0) *status = OPTMS_NO_EQUIL;
 
    if (*status != OPTMS_NO_EQUIL) {
      switch(num_max) {
      case 1: /* rotate to find the next point */
	OPTMS_COPY_VECTOR(pt1,vec[ind[0]],3);
	pt_1 = pt1[dir1];  pt_2 = pt1[dir2];
	if (pt1[dir2] < 0){rotate=OPTMS_CW; min_inv_slope=OPTMS_BIG_POS_NMBR;}
	if (pt1[dir2] >= 0){rotate=OPTMS_CCW; max_inv_slope=OPTMS_BIG_NEG_NMBR;}
	switch(rotate) {
	case OPTMS_CCW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope>max_inv_slope) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < OPTMS_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case OPTMS_CW:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope<min_inv_slope) {
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < OPTMS_MACHINE_EPS) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  OPTMS_DEBUG_PRINT(3,"No points in the rotation ... odd\n");
	  *status = OPTMS_HULL_TEST_ERROR;
	  break;
	case 1:
	  OPTMS_DEBUG_PRINT(3,"Found a line in the convex hull\n");
          OPTMS_COPY_VECTOR(pt2,vec[ind[1]],3);
	  *status = OPTMS_TWO_PT_PLANE;
	  break;
	default:
	  OPTMS_DEBUG_PRINT(3,"Found 2 or more points in the rotation\n");
	    /* check to see if rotation got past origin */
	  inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case OPTMS_CCW:
	    if (inv_origin_slope >= max_inv_slope) *status=OPTMS_NO_EQUIL;
	    else if (dir1 == 2) *status=OPTMS_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=OPTMS_CHECK_X_COORD_DIRECTION;
	    else *status=OPTMS_EQUIL;
	    break;
	  case OPTMS_CW:
	    if (inv_origin_slope <= min_inv_slope) *status=OPTMS_NO_EQUIL;
	    else if (dir1 == 2) *status=OPTMS_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) *status=OPTMS_CHECK_X_COORD_DIRECTION;
	    else *status=OPTMS_EQUIL;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	OPTMS_COPY_VECTOR(pt1,vec[ind[0]],3);
	OPTMS_COPY_VECTOR(pt2,vec[ind[1]],3);
	*status = OPTMS_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	OPTMS_DEBUG_ACTION(3,{fprintf(stdout,"Found 3 in max plane %f\n",max);});
	if (vec[ind[0]][dir1] <= 0) *status = OPTMS_NO_EQUIL;
	else if (dir1==2) *status=OPTMS_CHECK_Y_COORD_DIRECTION;
	else if (dir1==1) *status=OPTMS_CHECK_X_COORD_DIRECTION;
	else *status = OPTMS_EQUIL;
      }
    }
  }
  return(ierr=0);
}


