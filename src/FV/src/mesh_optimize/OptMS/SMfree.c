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

#undef __FUNC__
#define __FUNC__ "SMfreeLocalMesh" 
int  SMfreeLocalMesh(SMlocal_mesh *local_mesh) 
{
    int     ierr;
    int     i,num_tri,num_pts;

    num_pts = OPTMS_MAX_NUM_VTX;
    num_tri = OPTMS_MAX_NUM_TRI;
    for (i=0;i<num_pts;i++) {
	OPTMS_FREE(local_mesh->incident_vtx[i]);
    }
    for (i=0;i<num_tri;i++) {
	OPTMS_FREE(local_mesh->vtx_connectivity[i]);
    }
    OPTMS_FREE(local_mesh->original_pt);
    OPTMS_FREE(local_mesh->original_function);
    OPTMS_FREE(local_mesh->incident_vtx);
    OPTMS_FREE(local_mesh->vtx_connectivity);
    OPTMS_FREE(local_mesh->lap_info->laplacian_function);
    OPTMS_FREE(local_mesh->lap_info);
    ierr = SMfreeOpt(local_mesh->opt_info);  OPTMS_CHKERR(ierr);
    ierr = SMfreeLP(local_mesh,5,num_tri); OPTMS_CHKERR(ierr);
    OPTMS_FREE(local_mesh); 
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeOpt" 
int SMfreeOpt(SMoptimal *opt_info)
{
    int        ierr;
    int        i;
    int        num_values;

    num_values = OPTMS_MAX_NUM_TRI * OPTMS_DEFAULT_FUNC_PER_TRI;

    for (i=0;i<num_values;i++) {
      OPTMS_FREE(opt_info->gradient[i]);
   }
    for (i=0;i<OPTMS_MAX_G_NUM;i++)  {
      OPTMS_FREE(opt_info->G[i]);
    }
    for (i=0;i<OPTMS_MAX_DIM;i++)  {
      OPTMS_FREE(opt_info->PDG[i]);
    }

    OPTMS_FREE(opt_info->function);
    OPTMS_FREE(opt_info->test_function);
    OPTMS_FREE(opt_info->original_function);
    OPTMS_FREE(opt_info->gradient);
    OPTMS_FREE(opt_info->gs);
    OPTMS_FREE(opt_info->G);
    OPTMS_FREE(opt_info->PDG);
    OPTMS_FREE(opt_info->prev_active_values);

    /* free the active information sets */
    SMfreeActive(opt_info->active);
    SMfreeActive(opt_info->test_active);
    SMfreeActive(opt_info->original_active);
    
    OPTMS_FREE(opt_info);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeLP" 
int SMfreeLP(SMlocal_mesh *local_mesh, int num_active, int num_constraints)
{
    int        ierr;
    int        i;
    SMlp    *lp_info;

    lp_info = local_mesh->lp_info;
    for (i=0;i<num_active;i++)   {
        OPTMS_FREE(lp_info->AAmat[i]);
        OPTMS_FREE(lp_info->Amat_T[i]);
        OPTMS_FREE(lp_info->Amat_T_O[i]);
    }
    for (i=0;i<num_constraints;i++)  {
      OPTMS_FREE(lp_info->Amat[i]);
    }

    OPTMS_FREE(lp_info->ipivot);
    OPTMS_FREE(lp_info->Amat);
    OPTMS_FREE(lp_info->Amat_T);
    OPTMS_FREE(lp_info->Amat_T_O);
    OPTMS_FREE(lp_info->feasible_x);
    OPTMS_FREE(lp_info->free_ind);
    OPTMS_FREE(lp_info->active_ind);
    OPTMS_FREE(lp_info->b);
    OPTMS_FREE(lp_info->c);
    OPTMS_FREE(lp_info->pi);
    OPTMS_FREE(lp_info->s);
    OPTMS_FREE(lp_info->alpha);
    OPTMS_FREE(lp_info->step);
    OPTMS_FREE(lp_info->Bmat);
    OPTMS_FREE(lp_info->Bmat_T);
    OPTMS_FREE(lp_info->AAmat);

    OPTMS_FREE(lp_info);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeActive" 
int SMfreeActive(SMactive *active)
{
    int ierr;
    OPTMS_FREE(active->active_ind);
    OPTMS_FREE(active);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeParam" 
int SMfreeParam(SMparam *smooth_param)
{
    int ierr;
    OPTMS_FREE(smooth_param);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeProcinfo" 
int SMfreeProcinfo(SMprocinfo *procinfo)
{
    int ierr;
    OPTMS_FREE(procinfo);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMfreeQualityTable" 
int SMfreeQualityTable(SMquality_table *quality_table)
{
    int  i, ierr;
    for (i=0;i<quality_table->num_functions;i++) {
      OPTMS_FREE(quality_table->measure[i]);
    }
    /* DPS added */
    OPTMS_FREE(quality_table->measure);
    /* DPS end */

    OPTMS_FREE(quality_table);
    return(ierr=0);
}


