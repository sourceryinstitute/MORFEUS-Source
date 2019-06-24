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
#define __FUNC__ "SMvalidityCheck"
int SMvalidityCheck(SMlocal_mesh *local_mesh, SMparam *smooth_param, int *valid)
{
/*  validity tests : 1. the mesh is still valid based on righthandedness
                     2. the new mesh is indeed better than the original
*/

    int    ierr;
    SMoptimal *opt_info;

    *valid = 1;

    opt_info = local_mesh->opt_info;

    /* check that the mesh is still valid, based on right handedness. 
       Returns a 1 or a 0 */
    ierr = SMvalidMesh(local_mesh, valid); OPTMS_CHKERR(ierr);

    /* check to see that the mesh didn't get worse */
    if (local_mesh->original_value > local_mesh->current_active_value) {
       OPTMS_DEBUG_ACTION(1,{
         printf("The local mesh got worse; initial value %f; final value %f\n",
	       local_mesh->original_value,
 	       local_mesh->current_active_value );
       });
       *valid = 0;
    }

    OPTMS_DEBUG_ACTION(1,{
      if (!(*valid)) {
        int i;
        printf("The mesh is no longer valid\n");
        OPTMS_WRITE_BINARY_ORDERED_PTS(local_mesh);
        for (i=0;i<local_mesh->dimension;i++) 
   	   fprintf(stderr,"free_new[%d] = %f;  ",i,local_mesh->free_vtx[i]);
        fprintf(stderr,"\n");

        OPTMS_PRINT_ORDERED_PTS(local_mesh);
      }
    });
    return(ierr=0);
}    

#undef __FUNC__
#define __FUNC__ "SMvalidMesh"
int SMvalidMesh(SMlocal_mesh *local_mesh, int *valid)
{
    int    i, ierr;
    int    num_tri;
    int    ind1, ind2, ind3;
    double *free_vtx;

    (*valid) = 1;
    num_tri = local_mesh->num_tri;
    free_vtx = local_mesh->free_vtx;

    if (local_mesh->dimension == 2) {
      for (i=0;i<num_tri && (*valid);i++) {
        ind1 = local_mesh->vtx_connectivity[i][0]; 
        ind2 = local_mesh->vtx_connectivity[i][1]; 
        ierr = SMorient2D(free_vtx, 
                           local_mesh->incident_vtx[ind1], 
                           local_mesh->incident_vtx[ind2],valid);
               OPTMS_CHKERR(ierr);
      }
    } else if (local_mesh->dimension == 3) {
      for (i=0;i<num_tri && (*valid);i++) {
         ind1 = local_mesh->vtx_connectivity[i][0]; 
         ind2 = local_mesh->vtx_connectivity[i][1]; 
         ind3 = local_mesh->vtx_connectivity[i][2]; 
         ierr = SMorient3D(free_vtx,local_mesh->incident_vtx[ind1], 
                           local_mesh->incident_vtx[ind2], 
                           local_mesh->incident_vtx[ind3], valid);
                OPTMS_CHKERR(ierr);
	 if (*valid == -1) *valid = 0; /* if it's left handed consider it invalid */
      }
    } else {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension must be 2 or 3");
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMorient2D"
int SMorient2D(double *vtx1, double *vtx2, double *vtx3, int *valid)
{
    int ierr;
    double cross;

    *valid = 1;
    cross = (vtx2[OPTMS_XDIR] - vtx1[OPTMS_XDIR])*(vtx3[OPTMS_YDIR] - vtx2[OPTMS_YDIR])-
            (vtx3[OPTMS_XDIR] - vtx2[OPTMS_XDIR])*(vtx2[OPTMS_YDIR] - vtx1[OPTMS_YDIR]);

    if (cross < 0) *valid = 0;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMorient3D"
int SMorient3D(double *vtx1, double *vtx2, double *vtx3, double *free_vtx, int *valid) 
     /* Returns 1 if tet ABCD is right-handed, -1 if it's left-handed,
	and 0 if it's essentially a tie. */
{
  int ierr;
  double dEps = 1.e-13;

  double dXa = vtx1[0];
  double dYa = vtx1[1];
  double dZa = vtx1[2];

  double dXb = vtx2[0];
  double dYb = vtx2[1];
  double dZb = vtx2[2];

  double dXc = vtx3[0];
  double dYc = vtx3[1];
  double dZc = vtx3[2];

  double dXd = free_vtx[0];
  double dYd = free_vtx[1];
  double dZd = free_vtx[2];

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

  /* Use the length scale to get a better idea if the tet is flat or
     just really small. */
  if (fabs(dScale) < OPTMS_MACHINE_EPS) {
    *valid = 0;
    return(ierr=0);
  } else {
    dDet /= (dScale*dScale*dScale);
  }

  if (dDet > dEps) {
    *valid = 1;
  } else if (dDet < -dEps) {
    *valid = -1;
  } else {
    *valid = 0;
  }
  return (ierr=0);
}







