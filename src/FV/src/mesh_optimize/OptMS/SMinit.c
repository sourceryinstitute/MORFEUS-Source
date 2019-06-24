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
#include <string.h>
#include "SMsmooth.h"

#undef __FUNC__
#define __FUNC__ "SMmallocLocalMesh"
int SMmallocLocalMesh(SMlocal_mesh **local_mesh)
{
    int i,num_pts,num_tri,num_values;
    int ierr = 0;

    /* malloc more space than you think you'll ever need, so you only need to do
       it once.  Based on 3D problem */

    OPTMS_MALLOC((*local_mesh),(SMlocal_mesh *),sizeof(SMlocal_mesh),1);

    num_pts = OPTMS_MAX_NUM_VTX;
    num_tri = OPTMS_MAX_NUM_TRI;

    num_values = num_tri * OPTMS_DEFAULT_FUNC_PER_TRI;

    (*local_mesh)->num_values = num_values;

    OPTMS_MALLOC((*local_mesh)->incident_vtx,(double **),sizeof(double *)*num_pts,1);
    OPTMS_MALLOC((*local_mesh)->vtx_connectivity,(int **),sizeof(int *)*num_tri,1);

    for (i=0;i<num_pts;i++) {
        OPTMS_MALLOC((*local_mesh)->incident_vtx[i],(double *),
                     sizeof(double)*OPTMS_MAX_DIM,1);
    }
    for (i=0;i<num_tri;i++) {
        OPTMS_MALLOC((*local_mesh)->vtx_connectivity[i],(int *),
                     sizeof(int)*OPTMS_MAX_DIM,1);
    }

    OPTMS_MALLOC((*local_mesh)->original_pt,(double *),sizeof(double )*OPTMS_MAX_DIM,1);
    OPTMS_MALLOC((*local_mesh)->original_function,(double *),sizeof(double)*num_values,1);

    ierr = SMmallocLap(num_values,&(*local_mesh)->lap_info); OPTMS_CHKERR(ierr);
    ierr = SMmallocOpt(num_values,&(*local_mesh)->opt_info); OPTMS_CHKERR(ierr);
    ierr = SMmallocLP(5,num_tri,&(*local_mesh)->lp_info); OPTMS_CHKERR(ierr);
    ierr = SMinitLap(num_values,(*local_mesh)->lap_info); OPTMS_CHKERR(ierr);
    ierr = SMinitOpt(num_values,(*local_mesh)->opt_info); OPTMS_CHKERR(ierr);

    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmallocLap"
int SMmallocLap(int num_values, SMlap_info **lap_info)
{
    int ierr;
    OPTMS_MALLOC((*lap_info),(SMlap_info *),sizeof(SMlap_info),1);
    OPTMS_MALLOC((*lap_info)->laplacian_function,(double *),sizeof(double)*num_values,1);
    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmallocOpt"
int SMmallocOpt(int num_values, SMoptimal **opt_info)
{
    int i, ierr;

    OPTMS_MALLOC((*opt_info),(SMoptimal *),sizeof(SMoptimal),1);

    OPTMS_MALLOC((*opt_info)->function,(double *),sizeof(double)*num_values,1);
    OPTMS_MALLOC((*opt_info)->test_function,(double *),sizeof(double)*num_values,1);
    OPTMS_MALLOC((*opt_info)->original_function,(double *),sizeof(double)*num_values,1);
    OPTMS_MALLOC((*opt_info)->gradient,(double **),sizeof(double *)*num_values,1);
    OPTMS_MALLOC((*opt_info)->gs,(double *),sizeof(double)*num_values,1);
    OPTMS_MALLOC((*opt_info)->G,(double **),sizeof(double *)*OPTMS_MAX_G_NUM,1);
    OPTMS_MALLOC((*opt_info)->PDG,(double **),sizeof(double *)*OPTMS_MAX_DIM,1);

    for (i=0;i<OPTMS_MAX_DIM;i++) {
       OPTMS_MALLOC((*opt_info)->PDG[i],(double *),sizeof(double)*OPTMS_MAX_DIM,1);
    }
    for (i=0;i<num_values;i++) {
       OPTMS_MALLOC((*opt_info)->gradient[i],(double *),sizeof(double)*OPTMS_MAX_DIM,1);
    }
    for (i=0;i<OPTMS_MAX_G_NUM;i++) {
       OPTMS_MALLOC((*opt_info)->G[i],(double *),sizeof(double)*OPTMS_MAX_G_NUM,1);
    }

/*DPS original    OPTMS_MALLOC((*opt_info)->prev_active_values,(double *),
  sizeof(double)*OPTMS_MAX_OPT_ITER+5,1); */
/* Increased size to avoid stepping on memory in SMminmaxOpt */
    OPTMS_MALLOC((*opt_info)->prev_active_values,(double *),
                 sizeof(double)*OPTMS_MAX_OPT_ITER+8,1);

    /* malloc the active information sets */
    ierr = SMmallocActive(num_values,&((*opt_info)->active)); OPTMS_CHKERR(ierr);
    ierr = SMmallocActive(num_values,&((*opt_info)->test_active)); OPTMS_CHKERR(ierr);
    ierr = SMmallocActive(num_values,&((*opt_info)->original_active)); OPTMS_CHKERR(ierr);

    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmallocLP"
int SMmallocLP(int num_active, int num_constraints, SMlp **lp_info)
{
    int i, ierr;

    OPTMS_MALLOC((*lp_info),(SMlp *),sizeof(SMlp),1);

    OPTMS_MALLOC((*lp_info)->ipivot, (int *),sizeof(int)*num_active,1);
    OPTMS_MALLOC((*lp_info)->Amat,(double **),sizeof(double *)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->Amat_T,(double **),sizeof(double *)*num_active,1);
    OPTMS_MALLOC((*lp_info)->Amat_T_O,(double **),sizeof(double *)*num_active,1);
    OPTMS_MALLOC((*lp_info)->feasible_x,(double *),sizeof(double)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->free_ind,(int *),sizeof(int)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->active_ind,(int *),sizeof(int)*num_active,1);
    OPTMS_MALLOC((*lp_info)->b,(double *),sizeof(double)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->c,(double *),sizeof(double)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->pi,(double *),sizeof(double)*num_active,1);
    OPTMS_MALLOC((*lp_info)->s,(double *),sizeof(double)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->alpha,(double *),sizeof(double)*num_active,1);
    OPTMS_MALLOC((*lp_info)->step,(double *),sizeof(double)*num_constraints,1);
    OPTMS_MALLOC((*lp_info)->Bmat,(double *),sizeof(double)*num_active*num_active,1);
    OPTMS_MALLOC((*lp_info)->Bmat_T,(double *),sizeof(double)*num_active*num_active,1);
    OPTMS_MALLOC((*lp_info)->AAmat,(double **),sizeof(double *)*num_active,1);

    for (i=0;i<num_active;i++) {
        OPTMS_MALLOC((*lp_info)->AAmat[i],(double*),sizeof(double)*num_constraints,1);
        OPTMS_MALLOC((*lp_info)->Amat_T_O[i],(double *),sizeof(double)*num_constraints,1);
        OPTMS_MALLOC((*lp_info)->Amat_T[i],(double *),sizeof(double)*num_constraints,1);
    }
    for (i=0;i<num_constraints;i++)
         OPTMS_MALLOC((*lp_info)->Amat[i],(double *),sizeof(double)*num_active,1);

    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmallocActive"
int SMmallocActive(int num_values, SMactive **active)
{
    int ierr;
    OPTMS_MALLOC((*active),(SMactive *),sizeof(SMactive),1);
    OPTMS_MALLOC((*active)->active_ind,(int *),sizeof(int)*num_values,1);
    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMmallocQualityTable"
int SMmallocQualityTable(SMquality_table **quality_table)
{
    int i, ierr;
    char function_name[10][128];
    double target[10];
    int dimension[10];

    OPTMS_MALLOC((*quality_table),(SMquality_table *),sizeof(SMquality_table),1);

    (*quality_table)->num_functions = 10;
    /* 2D functions */
    strcpy(function_name[0],"Triangle Area");         target[0]=0; dimension[0]=2;
    strcpy(function_name[1],"Angle");                 target[1]=60; dimension[1]=2;
    strcpy(function_name[2],"Deviation from Equil");  target[2]=0; dimension[2]=2;
    strcpy(function_name[3],"Scaled Jacobian");       target[3]=.86; dimension[3]=2;
    strcpy(function_name[4],"Condition Number");      target[4]=1.0; dimension[4]=2;

    /* 3D functions */
    strcpy(function_name[5],"Tet Area");              target[5]=0; dimension[5]=3;
    strcpy(function_name[6],"Dihedral Angle");        target[6]=72; dimension[6]=3;
    strcpy(function_name[7],"Scaled Jacobian");       target[7]=.86; dimension[7]=3;
    strcpy(function_name[8],"Edge Length/Volume");    target[8]=.1; dimension[8]=3;
    strcpy(function_name[9],"Condition Number");      target[9]=1; dimension[9]=3;

    OPTMS_MALLOC((*quality_table)->measure,(SMqualityMeasure **),sizeof(SMqualityMeasure *)*(*quality_table)->num_functions,1);
    for (i=0;i<(*quality_table)->num_functions;i++) {
         OPTMS_MALLOC((*quality_table)->measure[i],(SMqualityMeasure *),sizeof(SMqualityMeasure),1);
         strcpy((*quality_table)->measure[i]->name, function_name[i]);
         (*quality_table)->measure[i]->dimension=dimension[i];
         (*quality_table)->measure[i]->target=target[i];
         (*quality_table)->measure[i]->min_value = OPTMS_BIG_POS_NMBR;
         (*quality_table)->measure[i]->max_value = OPTMS_BIG_NEG_NMBR;
         (*quality_table)->measure[i]->avg_value=0.0;
         (*quality_table)->measure[i]->num_function_values=0;
    }
    (*quality_table)->num_tangled_elements = 0;
    (*quality_table)->initialized=OPTMS_FALSE;
    return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitLocalMesh"
int  SMinitLocalMesh(int num_incident_vtx, int num_tri, double *free_vtx,
                      double **vtx_list, int **vtx_connectivity,
                      SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
    int           i,j;
    int           num_values;
    int           dimension;
    int           ierr;
    double        min[3], max[3];
    double        coord;

    SM_LOG_EVENT_BEGIN(__SM_INIT__);

    /* check for null data */
    OPTMS_CHECK_NULL(free_vtx);
    OPTMS_CHECK_NULL(vtx_list);
    OPTMS_CHECK_NULL(vtx_connectivity);
    OPTMS_CHECK_NULL(local_mesh);
    OPTMS_CHECK_NULL(smooth_param);

    /* check the input argurments for bad data */
    if (num_incident_vtx==0)
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"No incident vertices\n");
    if (num_tri==0)
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"No incident elements\n");

    if (num_incident_vtx>OPTMS_MAX_NUM_TRI) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,
               "Number of incident vertices exceeds OPTMS_MAX_NUM_TRI\n");
    }
    if (num_tri>2*OPTMS_MAX_NUM_TRI) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,
               "Number of incident triangles exceeds 2*OPTMS_MAX_NUM_TRI\n");
    }

    /* check for some simple incorrect input problems */
    dimension = local_mesh->dimension;
    if ((dimension==2) && (num_incident_vtx != num_tri)) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,
               "Number of incident vertices and elements are mismatched \n");
    } else if ((dimension==3) && (num_tri != 2*num_incident_vtx-4)) {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,
               "Number of incident vertices and elements are mismatched \n");
    }

    /* copy data to the local mesh data structure */
    local_mesh->num_incident_vtx = num_incident_vtx;
    local_mesh->num_tri = num_tri;

    num_values = num_tri * smooth_param->function_values_per_tri;
    local_mesh->num_values = num_values;

    OPTMS_COPY_VECTOR(local_mesh->free_vtx,free_vtx,dimension);
    OPTMS_COPY_VECTOR(local_mesh->original_pt,free_vtx,dimension);
    for (i=0;i<num_tri;i++) {
        OPTMS_COPY_VECTOR(local_mesh->vtx_connectivity[i],vtx_connectivity[i],dimension);
    }

    /* find the min and max coordinates in the local mesh */
    for (i=0;i<dimension;i++) { min[i]=1E300; max[i]=-1E300; }

    for (i=0;i<num_incident_vtx;i++) {
      for (j=0;j<dimension;j++) {
	coord = vtx_list[i][j];
	local_mesh->incident_vtx[i][j] = coord;
        if (coord < min[j]) min[j] = coord;
        if (coord > max[j]) max[j] = coord;
      }
    }
    if (local_mesh->free_vtx[OPTMS_XDIR] < min[OPTMS_XDIR]) min[OPTMS_XDIR]=local_mesh->free_vtx[OPTMS_XDIR];
    if (local_mesh->free_vtx[OPTMS_YDIR] < min[OPTMS_YDIR]) min[OPTMS_YDIR]=local_mesh->free_vtx[OPTMS_YDIR];
    if (local_mesh->free_vtx[OPTMS_XDIR] > max[OPTMS_XDIR]) max[OPTMS_XDIR]=local_mesh->free_vtx[OPTMS_XDIR];
    if (local_mesh->free_vtx[OPTMS_YDIR] > max[OPTMS_YDIR]) max[OPTMS_YDIR]=local_mesh->free_vtx[OPTMS_YDIR];
    OPTMS_COPY_VECTOR(local_mesh->min,min,dimension);
    OPTMS_COPY_VECTOR(local_mesh->max,max,dimension);

    /* initialize the optimization function values */
    local_mesh->original_value = 0;
    local_mesh->current_active_value = 0;
    for (i=0;i<local_mesh->num_values;i++) {
       local_mesh->original_function[0] = 0.0;
    }
    local_mesh->opt_info->opt_iter_count = 0;

    local_mesh->lap_done = OPTMS_FALSE;
    local_mesh->opt_done = OPTMS_FALSE;

    OPTMS_DEBUG_ACTION(3,{OPTMS_PRINT_ORDERED_PTS(local_mesh); });
    OPTMS_DEBUG_ACTION(3,{OPTMS_WRITE_BINARY_ORDERED_PTS(local_mesh);});

    SM_LOG_EVENT_END(__SM_INIT__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitLap"
int SMinitLap(int num_values, SMlap_info *lap_info)
{
    int i;
    int ierr;

    SM_LOG_EVENT_BEGIN(__SM_INIT__);

    OPTMS_CHECK_NULL(lap_info);

    lap_info->laplacian_value = 0;
    lap_info->lap_improvement = 0;
    lap_info->lap_invalid = OPTMS_FALSE;
    lap_info->lap_accepted = OPTMS_FALSE;
    for (i=0;i<num_values;i++) {
      lap_info->laplacian_function[i] = 0.;
    }
    SM_LOG_EVENT_END(__SM_INIT__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitLP"
int SMinitLP(SMlocal_mesh *local_mesh)
{
    int ierr;
    int i,j;
    int num_active, num_constraints;
    SMlp *lp_info;

    OPTMS_CHECK_NULL(local_mesh);

    lp_info = local_mesh->lp_info;

    OPTMS_CHECK_NULL(lp_info);

    num_active = 5;
    num_constraints = local_mesh->num_tri;

    for (i=0;i<num_constraints;i++) {
       lp_info->feasible_x[i]=0;
       lp_info->b[i]=0;
       lp_info->c[i]=0;
       lp_info->free_ind[i]=0;
       lp_info->s[i]=0;
       lp_info->step[i]=0;
    }
    for (i=0;i<num_active;i++) {
       lp_info->active_ind[i]=0;
       lp_info->pi[i]=0;
       lp_info->alpha[i]=0;
       lp_info->ipivot[i]=0;
       for (j=0;j<num_constraints;j++) {
           lp_info->AAmat[i][j] = 0;
           lp_info->Amat_T[i][j]=0;
           lp_info->Amat_T_O[i][j]=0;
           lp_info->Amat[j][i]=0;
       }
    }
    for (i=0;i<num_active*num_active;i++) {
       lp_info->Bmat[i]=0;
       lp_info->Bmat_T[i]=0;
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitSmoothParam"
int  SMinitSmoothParam(int technique, int FunctionID,
                       double Threshold, void *ext_smooth_data)
{
    int ierr;

    SMsmooth_data *smooth_data;
    SMparam *smooth_param;

    smooth_data = (SMsmooth_data *)ext_smooth_data;
    OPTMS_CHECK_NULL(smooth_data);
    smooth_param = smooth_data->smooth_param;
    OPTMS_CHECK_NULL(smooth_param);

    smooth_param->iter_count = 0;
    smooth_param->new_init_pt_option = OPTMS_NONE;
    smooth_param->maxit = OPTMS_MAX_OPT_ITER;
    smooth_param->conv_eps = 1E-10;
    smooth_param->active_eps = .00003;
    smooth_param->min_acceptable_imp = 1E-06;
    smooth_param->min_step_size = 1E-6;
    smooth_param->function_id = FunctionID;
    smooth_param->global_min_value = OPTMS_BIG_POS_NMBR;
    ierr = SMsetSmoothTechnique(ext_smooth_data,technique); OPTMS_CHKERR(ierr);
    ierr = SMsetSmoothFunction(ext_smooth_data,FunctionID); OPTMS_CHKERR(ierr);
    ierr = SMsetSmoothThreshold(ext_smooth_data,Threshold); OPTMS_CHKERR(ierr);

    return(ierr = 0);
}

#undef __FUNC__
#define __FUNC__ "SMinitProcinfo"
int SMinitProcinfo(SMprocinfo *procinfo)
{
    int ierr;
#ifdef PARALLEL_LOG
    int nprocs, myid;
#endif

    /* set up procinfo */
#ifdef PARALLEL_LOG
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    procinfo->nprocs  = nprocs;
    procinfo->myid    = myid;
    procinfo->procset = MPI_COMM_WORLD;
#else
    procinfo->nprocs = 1;
    procinfo->myid = 0;
#endif

    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitOpt"
int SMinitOpt(int num_values, SMoptimal *opt_info)
{
    int        ierr;
    int        i, j;

    SM_LOG_EVENT_BEGIN(__SM_INIT__);

    OPTMS_CHECK_NULL(opt_info);
    if (num_values > 900) {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"num_values exceeds 900");
    }

    /* for the purposes of initialization will be set to zero after */
    opt_info->num_values = num_values;
    opt_info->equilibrium_pt = 0;
    opt_info->step_too_small = 0;
    opt_info->step_accepted = 0;
    opt_info->status = 0;
    opt_info->iter_count = 0;
    opt_info->opt_iter_count = 0;
    opt_info->opt_improvement = 0;
    opt_info->maxit = OPTMS_MAX_OPT_ITER;
    opt_info->steepest = 0;
    opt_info->alpha = 0;
    opt_info->max_alpha = 0;

    for (i=0;i<OPTMS_MAX_DIM;i++) {
      opt_info->search[i] = 0;
      opt_info->PDG_ind[i] = -1;
      for (j=0;j<OPTMS_MAX_DIM;j++) opt_info->PDG[i][j] = 0;
    }
    for (i=0;i<num_values;i++) {
       opt_info->function[i] = 0;
       opt_info->test_function[i] = 0;
       opt_info->original_function[i] = 0;
       opt_info->gs[i] = 0;
       for (j=0;j<OPTMS_MAX_DIM;j++) {
           opt_info->gradient[i][j] = 0;
       }
    }
    if (num_values > OPTMS_MAX_G_NUM) {
      for (i=0;i<OPTMS_MAX_G_NUM;i++) {
       for (j=0;j<OPTMS_MAX_G_NUM;j++) opt_info->G[i][j] = -1;
      }
    } else {
      for (i=0;i<num_values;i++) {
       for (j=0;j<num_values;j++) opt_info->G[i][j] = -1;
      }
    }

    for (i=0;i<20;i++) opt_info->prev_active_values[i] = 0;
    SM_LOG_EVENT_END(__SM_INIT__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMinitMaxStepLength"
int SMinitMaxStepLength(SMlocal_mesh *local_mesh)
{
  int ierr;
  int i, j, k;
  double max_diff = 0;
  double diff=0;
  double coord1[3], coord2[3];

  /* check that the input data is correct */
  OPTMS_CHECK_NULL(local_mesh);
  OPTMS_CHECK_NULL(local_mesh->incident_vtx);
  OPTMS_CHECK_NULL(local_mesh->opt_info);
  if (local_mesh->num_incident_vtx==0) OPTMS_SETERR(OPTMS_INIT_ERR,0,"Num incident vtx = 0\n");
  if ((local_mesh->dimension!=2) && (local_mesh->dimension!=3)) {
     OPTMS_SETERR(OPTMS_INIT_ERR,0,"Problem dimension is incorrect\n");
  }

  /* find the maximum distance between two incident vertex locations */
  for (i=0;i<local_mesh->num_incident_vtx;i++) {
    for (j=i;j<local_mesh->num_incident_vtx;j++) {
      diff=0;
      for (k=0;k<local_mesh->dimension;k++) {
        coord1[k] = local_mesh->incident_vtx[i][k];
        coord2[k] = local_mesh->incident_vtx[j][k];
        diff += (coord1[k]-coord2[k])*(coord1[k]-coord2[k]);
      }
      if (max_diff < diff) max_diff=diff;
    }
  }
  max_diff = sqrt(max_diff);
  if (max_diff==0) {
     OPTMS_SETERR(OPTMS_INPUT_ERR,0,
            "Maximum distance between incident vertices = 0\n");
  }
  local_mesh->opt_info->max_alpha = max_diff/100;

  return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMinitStats"
int SMinitStats(SMstats *smooth_stats)
{
    int ierr = 0;
    SM_LOG_EVENT_BEGIN(__SM_INIT_STATS__);

    OPTMS_CHECK_NULL(smooth_stats);
    smooth_stats->stats_initialized = 1;
    smooth_stats->total_cells_smoothed = 0;
    smooth_stats->num_equil = 0;
    smooth_stats->num_started_equil = 0;
    smooth_stats->num_zero_search = 0;
    smooth_stats->num_imp_too_small = 0;
    smooth_stats->num_flat_no_imp = 0;
    smooth_stats->num_step_too_small = 0;
    smooth_stats->num_max_iter_exceeded = 0;
    smooth_stats->num_cells_opted = 0;
    smooth_stats->opt_count = 0.0;
    smooth_stats->opt_iter_count = 0;
    smooth_stats->num_cells_laplaced = 0;
    smooth_stats->num_lap_enough = 0;
    smooth_stats->num_lap_invalid = 0;
    smooth_stats->num_lap_worse = 0;
    smooth_stats->no_improvement = 0;
    smooth_stats->avg_improvement = 0.0;
    smooth_stats->avg_active_val = 0.0;
    smooth_stats->global_minimum_val = OPTMS_BIG_POS_NMBR;

    SM_LOG_EVENT_END(__SM_INIT_STATS__);
    return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMconvertToDegrees"
int SMconvertToDegrees(int function_id, double *value)
{
  int ierr=0;
  switch(function_id) {
  case OPTMS_MAX_MIN_ANGLE:
     *value = *value * 180 / 3.14159; break;
  case OPTMS_MIN_MAX_ANGLE:
     *value = *value * 180 / 3.14159; break;
  case OPTMS_MIN_MAX_COSINE:
     *value = acos(*value) * 180 / 3.14159; break;
  case OPTMS_MAX_MIN_COSINE:
     *value = acos(*value) * 180 / 3.14159; break;
  case OPTMS_MAX_MIN_SINE:
     *value = asin(*value) * 180 / 3.14159; break;
  case OPTMS_MAX_MIN_DIHEDRAL:
     *value = *value * 180 / 3.14159; break;
  case OPTMS_MIN_MAX_DIHEDRAL:
     *value = *value * 180 / 3.14159; break;
  case OPTMS_MIN_MAX_COSINE_DIHEDRAL:
     *value = acos(*value) * 180 / 3.14159; break;
  case OPTMS_MAX_MIN_COSINE_DIHEDRAL:
     *value = acos(*value) * 180 / 3.14159; break;
  case OPTMS_MAX_SINE_DIHEDRAL:
     *value = asin(*value) * 180 / 3.14159; break;
  default:
     *value = asin(*value) * 180 / 3.14159; break;
  }
  OPTMS_DEBUG_ACTION(3,{fprintf(stdout,"Coverting the min value to degrees %f\n",*value);});
  return(ierr);
}
