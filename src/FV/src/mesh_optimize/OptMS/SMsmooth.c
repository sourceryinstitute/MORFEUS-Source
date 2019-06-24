/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <stdio.h>
#include "SMsmooth.h"

/*@
   SMsmooth - This is the main routine that optimizes a local submesh.  Most of the quality 
      functions available require that the local submesh is initially valid.  
      If the user suspects that the mesh contains invalid elements, a quality 
      assessment should be done, and, if necessary, SMuntangle should be used 
      to try to rectify the problem.

   Input Parameters:
+  num_incident_vtx - the number of incident vertices in the local mesh 
.  num_tri - the number of incident triangles or tetrahedra
.  free_vtx - the coordinates of the free vertex 
           a vector of length equal to the problem dimension in x, y, z order
.  vtx_list - a matrix of the coordinates of the incident vtx;
           matrix dimensions are num_incident_vtx by problem dimension
.  vtx_connectivity - a matrix that gives the connectivity info for the
           incident vertices. matrix dimensions are num_incident_vtx by 
           the problem dimension. Note: this assumes that the connectivity 
           given is for a right handed triangle or tetrahedra the free vertex 
           ordered first
-  ext_smooth_data - data structure for mesh smoothing; created in 
          SMinitSmoothing and cast to the local data structure

   Output Parameter:
.  free_vtx -  contains the new coordinates of the free vertex
 
   Note:
   This function can be called only after SMinitSmoothing has been called to create
    the ext_smooth_data data structure.  Once the mesh has been optimized,
    SMfinalizeSmoothing should be called to release the memory.

.seealso SMinitSmoothing(), SMsetSmoothTechnique(), SMsetSmoothFunction(), 
             SMsetSmoothThreshold(), SMfinializeSmoothing()
@*/
#undef __FUNC__
#define __FUNC__ "SMsmooth"
int SMsmooth(int num_incident_vtx, int num_tri, double *free_vtx, 
              double **vtx_list, int **vtx_connectivity, 
              void *ext_smooth_data)
{
    int ierr;
    int i;
    int dimension;
    int valid;
    double min_value;
    SMsmooth_data *smooth_data;
    SMlocal_mesh *local_mesh;
    SMparam *smooth_param;
    SMstats *smooth_stats;
    SMprocinfo *procinfo;
    SMquality_table *quality_table;
    /*    FILE *fp; */

#ifdef OPTMS_LOCALTEST
    FILE *fp;
#endif

    SM_LOG_EVENT_BEGIN(__SM_TOTAL__);
    SM_LOG_EVENT_BEGIN(__SM_SMOOTH_INIT__);

    OPTMS_DEBUG_PRINT(2,"\nSmoothing a new local submesh\n");

    /* check for null data pointer for Opt-MS */
    OPTMS_CHECK_NULL(ext_smooth_data);

    smooth_data = (SMsmooth_data *) ext_smooth_data;
    local_mesh = smooth_data->local_mesh;
    smooth_param = smooth_data->smooth_param;
    smooth_stats = smooth_data->smooth_stats;
    procinfo = smooth_data->smooth_procinfo;
    quality_table = smooth_data->quality_table;
    dimension = local_mesh->dimension;

    /* initialize the local mesh triangles and the smoothing parameters */
    /* this routine checks for input errors */
    ierr = SMinitLocalMesh(num_incident_vtx,num_tri,free_vtx,vtx_list,
                    vtx_connectivity,local_mesh,smooth_param); OPTMS_CHKERR(ierr);

    SM_LOG_EVENT_END(__SM_SMOOTH_INIT__);

    /* if the technique chosen is Stupid Laplacian, do it and leave */
    SM_LOG_EVENT_BEGIN(__SM_SMOOTH__);
    if (smooth_param->smooth_technique == OPTMS_LAPLACIAN_ONLY) {
        SM_LOG_EVENT_BEGIN(__SM_LAP_SMOOTH__);
        OPTMS_MATLAB_ON({
	if (OPTMS_ISROOT(procinfo) && (dimension==2)) {
	    /* plot out the initial local mesh */
	    OPTMS_DEBUG_PRINT(1,"Plotting the initial local mesh \n");
	    if ((fp = fopen("initlocal.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open initlocal.m for writing\n");
            } 
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
        OPTMS_MATLAB_OFF});

        ierr = SMcentroidSmoothMesh(local_mesh->num_incident_vtx, 
                      local_mesh->incident_vtx, 
                      local_mesh->free_vtx,local_mesh->dimension); 
               OPTMS_CHKERR(ierr);
        OPTMS_COPY_VECTOR(free_vtx,local_mesh->free_vtx,dimension);

        OPTMS_MATLAB_ON({
	if (OPTMS_ISROOT(procinfo) && (dimension==2)) {
	    /* plot out the initial local mesh */
	    OPTMS_DEBUG_PRINT(1,"Plotting the laplaced local mesh \n");
	    if ((fp = fopen("laplace_local.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open local_laplace.m for writing\n");
            } 
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
        OPTMS_MATLAB_OFF});
        SM_LOG_EVENT_END(__SM_LAP_SMOOTH__);
        SM_LOG_EVENT_END(__SM_SMOOTH__);
        SM_LOG_EVENT_END(__SM_TOTAL__);
        SM_LOG_GLOBAL_TIME(__SM_TOTAL__);
        OPTMS_STATS_ON({
             OPTMS_DEBUG_PRINT(2,"Creating stats \n");
             ierr = SMaccumulateStats(local_mesh,smooth_param,smooth_stats);
                    OPTMS_CHKERR(ierr);
        OPTMS_STATS_OFF});
        return(ierr=0);  /* leaving after stupid laplacian smoothing */
     }

    /* check for an invalid mesh; if it's invalid return and ask the user to use SMuntangle */
    ierr = SMvalidityCheck(local_mesh,smooth_param,&valid); OPTMS_CHKERR(ierr);
    if (!valid) {
         ierr = SMwrite_ordered_points(local_mesh); OPTMS_CHKERR(ierr);
         OPTMS_SETERR(OPTMS_INVALID_MESH_ERR,0,"Invalid Mesh: Use SMuntangle to create a valid triangulation");
    }
    local_mesh->validity=OPTMS_VALID_MESH;
      
    /* compute the original function values */
    ierr = SMcomputeFunction(local_mesh,smooth_param,local_mesh->original_function);
           OPTMS_CHKERR(ierr);

    /* find the minimum value */
    min_value = 1E300;
    for (i=0;i<local_mesh->num_values;i++) 
       min_value = OPTMS_MIN(local_mesh->original_function[i],min_value);

    local_mesh->original_value = local_mesh->current_active_value = min_value;
	
    OPTMS_DEBUG_ACTION(3,
        {fprintf(stdout,"The initial minimum value is %f\n",min_value); 
    });

    OPTMS_MATLAB_ON({
	if (OPTMS_ISROOT(procinfo) && (dimension==2)) {
	    /* plot out the initial local mesh */
	    OPTMS_DEBUG_PRINT(1,"Plotting the initial local mesh \n");
	    if ((fp = fopen("initlocal.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open local_laplace.m for writing\n");
            } 
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
    OPTMS_MATLAB_OFF});


    /**********************************************************************************************************
    Methods available in the smoothing code:

    OPTMS_LAPLACIAN_ONLY: take a laplacian step only if it improves the mesh and leave
    OPTMS_OPTIMIZATION_ONLY: for each local submesh do optimization based smoothing.
          This method is not recommended due to it's expense
    OPTMS COMBINED 1: if the minimum value in the local submesh is too small, do 
           optimization based smoothing, otherwise do Laplacian smoothing
    OPTMS COMBINED 2: always do a Laplacian step; if the minimum value of the local
           submesh is too small following this step, perform optimization-based smoothing.
           This is the recommended technique as it gives general improvement to the mesh
           and eliminates extremal angles at roughly twice the cost of Laplacian smoothing
           used alone.
     OPTMS COMBINED 3: if the minimum value of the local submesh exceeds the threshold value,
           do nothing, otherwise use the COMBINED 2 approach.  This is the cheapest of the
           three combined techniques and gives reasonable results (although quality metric
           distribution can clump near the threshold value.
    **********************************************************************************************************/
    switch(smooth_param->smooth_technique) {
    case OPTMS_SMART_LAPLACIAN_ONLY:
        ierr = SMlaplaceSmooth(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        break;
    case OPTMS_OPTIMIZATION_ONLY:
        ierr = SMminmaxOpt(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        break;
    case OPTMS_COMBINED1:
        if (smooth_param->lap_accept_value<=local_mesh->current_active_value){
          /* do Laplacian */
          ierr = SMlaplaceSmooth(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        } else {
          /* do optimization */
          ierr = SMminmaxOpt(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        }
        break;
    case OPTMS_COMBINED2:
        ierr = SMlaplaceSmooth(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
            ierr = SMminmaxOpt(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
            OPTMS_DEBUG_PRINT(2,"Current Active value less than threshold, Doing Optimization\n");
        }
        break;
    case OPTMS_COMBINED3:
        if (smooth_param->lap_accept_value>local_mesh->current_active_value){
          /* do laplacian */
          ierr = SMlaplaceSmooth(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
          /* decide if optimization is necessary */
          if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
            ierr = SMminmaxOpt(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
            OPTMS_DEBUG_PRINT(2,"Current Active value less than threshold, Doing Optimization\n");
          }
        }
        break;
    case OPTMS_COMBINED:
    case OPTMS_FLOATING_THRESHOLD:
        ierr = SMlaplaceSmooth(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
        if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
            ierr = SMminmaxOpt(local_mesh,smooth_param,procinfo); OPTMS_CHKERR(ierr);
            OPTMS_DEBUG_PRINT(2,"Current Active value less than threshold, Doing Optimization\n");
        }
        break;
    default:
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Invalid smoothing technique");
    }

    SM_LOG_EVENT_BEGIN(__SM_SMOOTH_FINAL__);

    /* checking the global minimum value */
    if (smooth_param->global_min_value > local_mesh->current_active_value) {
      smooth_param->global_min_value = local_mesh->current_active_value;
    }

    OPTMS_MATLAB_ON({
	if (OPTMS_ISROOT(procinfo) && (dimension==2)) {
	    /* plot out the local mesh */
	    OPTMS_DEBUG_PRINT(3,"Plotting the local mesh \n");
	    if ((fp = fopen("local.m","w")) == NULL) {
                OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open local.m for writing\n");
            }

	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    if (do_optimization==OPTMS_TRUE) {
	       /* write out the active set if there is one */
               ierr = SMwriteActiveSet(fp,local_mesh);  OPTMS_CHKERR(ierr);
	    }
	    fclose(fp);
	}
    OPTMS_MATLAB_OFF});

    OPTMS_COPY_VECTOR(free_vtx,local_mesh->free_vtx,dimension);
    OPTMS_DEBUG_PRINT(2,"Returning to example program\n");

    SM_LOG_EVENT_END(__SM_SMOOTH_FINAL__);
    SM_LOG_EVENT_END(__SM_SMOOTH__);
    SM_LOG_EVENT_END(__SM_TOTAL__);
    SM_LOG_GLOBAL_TIME(__SM_TOTAL__);

    OPTMS_STATS_ON({
	OPTMS_DEBUG_PRINT(2,"Creating stats \n");
	ierr = SMaccumulateStats(local_mesh,smooth_param,smooth_stats);
        OPTMS_CHKERR(ierr);
    OPTMS_STATS_OFF});

    return(ierr = 0);
}


