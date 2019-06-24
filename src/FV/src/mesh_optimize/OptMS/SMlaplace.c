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
#define __FUNC__ "SMlaplaceSmooth"
int SMlaplaceSmooth(SMlocal_mesh *local_mesh, SMparam *smooth_param,
                     SMprocinfo *procinfo)
{
    int    ierr = 0;
    int    i, valid = 0;
    int    num_values;
    double min_value, func_diff;
    double *laplace_function;
    OPTMS_MATLAB_ON({
          FILE *fp;
    });

    SM_LOG_EVENT_BEGIN(__SM_LAP_SMOOTH__);

    OPTMS_CHECK_NULL(local_mesh);
    OPTMS_CHECK_NULL(smooth_param);
    OPTMS_CHECK_NULL(procinfo);

    local_mesh->lap_done = OPTMS_TRUE;

    num_values = local_mesh->num_values;

    ierr = SMinitLap(local_mesh->num_values,local_mesh->lap_info); OPTMS_CHKERR(ierr);

    ierr = SMcentroidSmoothMesh(local_mesh->num_incident_vtx,
        local_mesh->incident_vtx, local_mesh->free_vtx,local_mesh->dimension);
        OPTMS_CHKERR(ierr);

    /* check the validity of the new point */
    ierr = SMvalidMesh(local_mesh, &valid); OPTMS_CHKERR(ierr);
    if (!valid && (local_mesh->validity==OPTMS_VALID_MESH)){

      /* don't use this step if this step makes a previously valid mesh invalid*/
      OPTMS_DEBUG_PRINT(2,"Did Not Accept Laplacian Smoothing\n");
      OPTMS_COPY_VECTOR(local_mesh->free_vtx,local_mesh->original_pt,local_mesh->dimension);
      /* write this out so that we can keep stats */
      local_mesh->lap_info->lap_invalid = OPTMS_TRUE;
      local_mesh->lap_info->lap_accepted = OPTMS_FALSE;
      local_mesh->lap_info->laplacian_value = local_mesh->original_value;
    } else {

      /* compute the new function values */
      ierr = SMcomputeFunction(local_mesh,smooth_param,
                        local_mesh->lap_info->laplacian_function); OPTMS_CHKERR(ierr);
      laplace_function = local_mesh->lap_info->laplacian_function;

      /* find the minimum value */
      min_value = 1E300;
      for (i=0;i<num_values;i++) {
         min_value = OPTMS_MIN(laplace_function[i],min_value);
      }

      /* if this step improves the function, keep it */
      func_diff = min_value - local_mesh->original_value;

      OPTMS_DEBUG_ACTION(2,{
	fprintf(stdout,"The improvement from Laplacian smoothing is %f\n",
	       func_diff);});
      if (func_diff > 0) {

	OPTMS_DEBUG_PRINT(2,"Accepted Laplacian Smoothing \n");
              local_mesh->lap_info->lap_accepted = OPTMS_TRUE;

              /* if the mesh is now has gone from invalid to valid change the status
                  of mesh->validity */
              if (valid) local_mesh->validity=OPTMS_VALID_MESH;

	/* set the iteration counter to 1 */
	smooth_param->iter_count += 1;
	local_mesh->lap_info->lap_improvement = func_diff;
              local_mesh->lap_info->laplacian_value = min_value;
              local_mesh->current_active_value = min_value;

	OPTMS_MATLAB_ON({
	  if (OPTMS_ISROOT(procinfo) && (local_mesh->dimension==2)) {
	    /* plot the vertices of the mesh */
	    OPTMS_DEBUG_PRINT(2,"Plotting the results of centroid smoothing \n");
	    if ((fp = fopen("centroid.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open centroid.m for writing\n");
            }
	    if (smooth_param->new_init_pt_option == OPTMS_CENTROID) {
	      ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	      ierr = SMwriteActiveSet(fp,local_mesh); OPTMS_CHKERR(ierr);
	    }
	    fclose(fp);
	  }
        OPTMS_MATLAB_OFF});

      } else {
	/* don't use this step */
	OPTMS_DEBUG_PRINT(2,"Did Not Accept Laplacian Smoothing\n");
              local_mesh->lap_info->lap_accepted = OPTMS_FALSE;
              local_mesh->lap_info->laplacian_value = local_mesh->original_value;
              OPTMS_COPY_VECTOR(local_mesh->lap_info->laplacian_function,
                           local_mesh->original_function, local_mesh->num_values);
              OPTMS_COPY_VECTOR(local_mesh->free_vtx,local_mesh->original_pt,
                                                  local_mesh->dimension);
      }
    }
    SM_LOG_EVENT_END(__SM_LAP_SMOOTH__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcentroidSmoothMesh"
int SMcentroidSmoothMesh(int num_incident_vtx, double **incident_vtx,
                          double *free_vtx, int dimension)
{
    int ierr;
    int i,j;
    double avg[OPTMS_MAX_DIM];

    OPTMS_CHECK_NULL(incident_vtx);
    OPTMS_CHECK_NULL(free_vtx);

    if (num_incident_vtx==0)
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Number of incident vertex is zero");

    for (j=0;j<dimension;j++) {
       avg[j] = 0.;
       for (i=0;i<num_incident_vtx;i++)  avg[j]+=incident_vtx[i][j];
       free_vtx[j] = avg[j]/num_incident_vtx;
       OPTMS_DEBUG_ACTION(3,{
           fprintf(stdout,"final --  avg[%d] = %f\n",j, free_vtx[j]); });
    }
    return(ierr=0);
}
