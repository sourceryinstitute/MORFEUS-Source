/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_DATA_STRUCTS_H
#define SM_DATA_STRUCTS_H 1

typedef struct __SMactive {
    int      num_active;
    int      num_equal;
    int      *active_ind;
    double   true_active_value;    /* the true active value */
} SMactive;

typedef struct __SMoptimal {
    int      dimension;              /* what dimesion problem is it*/
    int      equilibrium_pt;       /* termination cond: equil point found */
    int      step_too_small;       /* termination cond; steps are too small */
    int      step_accepted;        /* accept this optimization step */
    int      status;               /* termination status, starts at 0 */
    int      iter_count;           /* the iteration count */
    int      opt_iter_count;       /* the number of optimization steps*/
    double   opt_improvement;      /* the improvement due to optimization */
    int      maxit;                /* max iteration count */
    int      num_values;
    double   *function;            /* an array of function values */
    double   **gradient;           /* an array of gradient values */
    double   *test_function;
    double   *original_function;
    SMactive *active;
    SMactive *test_active;
    SMactive *original_active;
    int      steepest;
    double   search[OPTMS_MAX_DIM];      /* the search direction */
    double   alpha;                /* the step length */
    double   max_alpha;            /* the maximum step length */
    double   *gs;                  /* projections of gradient on search */
    double   *prev_active_values;  /* a record of true active values */
    double   **G;                  /* the grammian matrix */
    double   **PDG;                /* a grammian of only LI gradients */
    int      PDG_ind[OPTMS_MAX_DIM];     /* the indices of the LI gradients */
    int      num_LI;               /* the number of LI gradients */
} SMoptimal;

typedef struct __SMlp {
    double **Amat;        /* contains the matrix with the constraints */
    double **Amat_T;
    double **Amat_T_O;  /*contains the constraint matrix + a row of ones for the slack variable */
    double *feasible_x;      /* contains the current feasible iterate */
    double *c;                    /* contains the rhs of the constraint equations */
    int *ipivot;             /* an array containing the pivot rows for lapack solve */
    int *active_ind;      /* an array containing the active indices for the LP solve*/
    int *free_ind;         /* an array containing the free indices for the LP solve*/
    double *b;              /* an array containing the rhs (0 0 0 ... 1) */
    double *pi;             /* the array containing x, y, min area */
    double *s;               /* the slack variables */
    double *alpha;        /* array containing the step size before we go infeasible */
    double *step;          /* the pivots */
    double *Bmat;         /* the active submatrix */
    double *Bmat_T;
    double **AAmat;    /* matrix containing the constraints and a row of ones */
 } SMlp;

typedef struct __SMstats {
    int                 stats_initialized;
    int                 total_cells_smoothed;
    int                 num_equil;
    int                 num_started_equil;
    int                 num_zero_search;
    int                 num_imp_too_small;
    int                 num_flat_no_imp;
    int                 num_step_too_small;
    int                 num_max_iter_exceeded;
    int                 num_cells_opted;
    int                 num_cells_laplaced;
    int                 num_lap_enough;
    int                 num_lap_invalid;
    int                 num_lap_worse;
    double           opt_count;
    int                 opt_iter_count;
    int                 no_improvement;
    double           avg_improvement;
    double           avg_active_val;
    double           global_minimum_val;
} SMstats;

typedef struct __SMlap_info {
    double laplacian_value;
    double *laplacian_function;
    double lap_improvement;
    int    lap_invalid;
    int    lap_accepted;
} SMlap_info;

typedef struct __SMlocal_mesh {
    int             dimension;
    int             num_incident_vtx;
    int 	    num_tri;
    double          min[OPTMS_MAX_DIM];
    double          max[OPTMS_MAX_DIM];
    double          free_vtx[OPTMS_MAX_DIM];
    double          **incident_vtx;
    int             **vtx_connectivity;
    double          *original_pt;
    int              num_values;
    double           original_value;
    double          *original_function;
    double           current_active_value;
    int              lap_done;
    int              opt_done;
    int              validity;
    SMoptimal        *opt_info;
    SMlap_info       *lap_info;
    SMlp             *lp_info;
} SMlocal_mesh;

typedef struct __SMparam {
    int         iter_count;
    int         new_init_pt_option;
    int         maxit;
    double      conv_eps;
    double      active_eps;
    double      min_acceptable_imp;
    double      min_step_size;
    int         function_id;
    int         smooth_technique;
    int         function_values_per_tri;
    double      lap_accept_value;
    double      global_min_value;
    SMfunction_ptr2D      ComputeFunctionValues2D;
    SMgradfunc_ptr2D      ComputeGradientValues2D;
    SMfunction_ptr3D      ComputeFunctionValues3D;
    SMgradfunc_ptr3D      ComputeGradientValues3D;
} SMparam;

#ifdef PARALLEL_LOG
typedef struct __SMprocinfo {
   int       nprocs;
   int       myid;
   MPI_Comm  procset;
} SMprocinfo;
#else
typedef struct __SMprocinfo {
   int       nprocs;
   int       myid;
} SMprocinfo;
#endif

typedef struct SMqualityMeasure {
  int    dimension;
  char name[128];
  double target;
  double  min_value;
  double  max_value;
  double avg_value;
  double avg_min_value;
  double avg_max_value;
  int num_function_values;
  int num_elements;
} SMqualityMeasure;

typedef struct SMquality_table {
   int   initialized;
   int   num_functions;
   int num_tangled_elements;
   int mesh_validity;
   SMqualityMeasure **measure;
} SMquality_table;

typedef struct __SMuntangle_param {
   int         untangle_technique;
} SMuntangle_param;

typedef struct __SMsmooth_data {
  int               dimension;
  SMlocal_mesh     *local_mesh;
  SMparam          *smooth_param;
  SMstats          *smooth_stats;
  SMquality_table  *quality_table;
  SMprocinfo       *smooth_procinfo;
  SMuntangle_param *untangle_param;
} SMsmooth_data;

#endif
