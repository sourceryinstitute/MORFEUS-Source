/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_INTERNAL_FUNC_H
#define SM_INTERNAL_FUNC_H 1

/* Function Declarations */

/* INTERNAL FUNCTIONS */
/* routines for mallocing and initializing space */
int   SMmallocLocalMesh(SMlocal_mesh **local_mesh);
int   SMmallocLap(int num_values, SMlap_info **lap_info);
int   SMmallocOpt(int num_values, SMoptimal **opt_info);
int   SMmallocLP(int num_active, int num_tri, SMlp **lp_info);
int   SMmallocActive(int num_values, SMactive **active);
int   SMmallocQualityTable(SMquality_table **quality_table);
int   SMinitLocalMesh(int num_pts, int num_tet,
                      double *free_vtx, double **vtx_list,
                      int **vtx_connectivity, SMlocal_mesh *local_mesh,
                      SMparam *smooth_param);
int   SMinitSmoothParam(int technique, int FunctionID,
                        double AcceptFunction, void *ext_smooth_data);
int   SMinitStats(SMstats *smooth_stats);
int   SMinitProcinfo(SMprocinfo *proc_info);
int   SMinitLap(int num_values, SMlap_info *lap_info);
int   SMinitOpt(int num_values, SMoptimal *opt_info);
int   SMinitMaxStepLength(SMlocal_mesh *local_mesh);
int   SMinitLP(SMlocal_mesh *local_mesh);
int   SMconvertToDegrees(int,double *);

/* function-gradient routines -- 2D */
int   SMcomputeFunction(SMlocal_mesh *local_mesh,
                        SMparam *smooth_param, double *function);
int   SMcomputeGradient(SMlocal_mesh *local_mesh,
                        SMparam *smooth_param, double **gradient);

int   SMcomputeTriCosines(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeCosGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient,int *num_values);

int   SMcomputeInteriorTriCosines(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeInteriorCosGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient,int *num_values);

int   SMcomputeNegTriCosines(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeNegCosGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeTriAngles(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeAngGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient,int *num_values);

int   SMcomputeInteriorTriAngles(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeInteriorAngGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeNegTriAngles(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeNegAngGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient,int *num_values);

int   SMcomputeTriSines(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeSineGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient,int *num_values);

int   SMcomputeInteriorTriSines(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeInteriorSineGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeTriJacobians(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeJacobianGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeScaledTriJacobians(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeScaledJacobianGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeInteriorScaledTriJacobians(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeInteriorScaledJacobianGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeAreaLengthRatio(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeAreaLengthRatioGradients(double *vtx1, double *vtx2,
                         double *vtx3, double **gradient, int *num_values);

int   SMcomputeLengthAreaRatio(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);
int   SMcomputeTriArea(double *vtx1, double *vtx2,
                         double *vtx3, double *function, int *num_values);

int   SMNormJacSquared2D(double *vtx1, double *vtx2, double *vtx3,
		       double *function, int *num_values);
int   SMcomputeNormJacSquaredGradients2D(double *vtx1,
                       double *vtx2, double *vtx3, double **gradient,
                       int *num_values);

int   SMcondition2D(double *vtx1, double *vtx2, double *vtx3,
                    double *function, int *num_values);
int   SMgradCondition2D(double *vtx1, double *vtx2, double *vtx3,
                        double **gradient, int *num_values);

/* function - gradient routines 3D */
#include "SMderiv.h"
#include "SMintrinsic.h"
#include "SMdihed_func.h"

/* Laplacian smoothing routines */
int   SMlaplaceSmooth(SMlocal_mesh *local_mesh,
                      SMparam *smooth_param, SMprocinfo *procinfo);
int   SMcentroidSmoothMesh(int num_incident,double **incident_vtx,
                           double *free_vtx, int dimension);

/* L infinity optimization routines */
int   SMminmaxOpt(SMlocal_mesh *local_mesh, SMparam *smooth_param, SMprocinfo *procinfo);
int   SMfindActiveSet(int num_values, double *function,
                      double active_eps, SMactive *active_info);
int   SMgetActiveDirections(int num_active, double **gradient,
                           int *active_ind, int dimension, double ***dir);
int   SMsearchDirection(SMlocal_mesh *local_mesh);
int   SMsearchEdgesFaces(int num_active, double **G, double **dir,
                         SMoptimal *opt_info);
int   SMprintActiveSet(SMactive *active, double *function);
int   SMstepAcceptance(SMlocal_mesh *local_mesh, SMparam *smooth_param);
int   SMstepToCusp(SMlocal_mesh *local_mesh,SMparam *smooth_param,
                   SMprocinfo *procinfo);
int   SMcomputeVerticalStep(SMlocal_mesh *local_mesh,
                            SMparam *smooth_param);
int   SMformVerticalRHS(int n,double *function, int ind[OPTMS_MAX_DIM],
                         double max_value, double **R);

int   SMgetGradientProjections(SMoptimal *opt_info);
int   SMcomputeAlpha(SMoptimal *opt_info);

int   SMgetMinEstimate(SMoptimal *opt_info, double *est);
int   SMcheckEquilibrium(SMoptimal *opt_info, int *equil);

int   SMconvexHullTest(double **vec, int num_vec, int *equil);
int   SMfindPlaneNormal(double pt1[3], double pt2[3],
                        double pt3[3],double *cross);
int   SMcheckVectorDots(double **vec,int num_vec,double *normal, int *equil);
int   SMfindPlanePoints(int dir1, int dir2, double **vec,
                       int num_vec, double *pt1,
                       double *pt2, double *pt3, int *status);
int   SMcopyActive(SMactive *active1, SMactive *active2);


/* untangling routines */
int   SMuntangle_mesh(SMlocal_mesh *local_mesh, int *degenerate);
int   SMdegenerate(int num_active, int num_constraints, double **A, double *b,
                   int *degenerate);
int   SMcomputeConstraintMatrix(SMlocal_mesh *local_mesh, int num_constriants,
                   double **Amat, double *b);
int   SMphaseOneLP(int num_constraints, int num_active, double **A, double *x,
                   SMlp *lp_info, int *feasible);
int   SMsolveLP(int num_constraints, int num_active, double **A, double *x, double *c, double *pi,
                   SMlp *lp_info, int *solved);
int   SMgetActiveMatrix(double **AAmat,int num_active,int *active_ind, double *Bmat);
int   SMgetActiveRHS(double *c,int num_active,int *active_ind, double *c_active);
int   SMremoveIdenticalVtx(int dimension, int *num_incident_vtx,int *num_tri,
                   double ***vtx_list, int ***vtx_connectivity, int *point_removed);

/* stats routines */
int   SMregisterEvents();
int   SMaccumulateStats(SMlocal_mesh *local_mesh, SMparam *smooth_param,
	                SMstats *smooth_stats);
int   SMprintStats(SMsmooth_data *);
int   SMwriteStatsFile(SMstats *smooth_stats, int smooth_count);

/* matrix routines */
int   SMformGrammian(int num_vecs, double **vecs, double **G, int dimension);
int   SMformPDGrammian(SMoptimal *opt_info);
int   SMformReducedMatrix(int num_active, double **G, double ***P);
int   SMformVerticalMatrix(int num_active, double **PDG, double ***N);
int   SMsolve2x2(double a11, double a12, double a21, double a22,
		 double b1, double b2, double **x);
int   SMsolveSPD3x3(double **A, double *B, double **x);
int   SMsingularTest(int n, double **A, int *singular);
int   SMcondition3x3(double **A, double *cond);
int   SMtransposeMatrix(double *mat, int n, int m, double *mat_T);
int   SMtransposeMatrix2(double **mat, int n, int m, double **mat_T);
int   SMdeterminant2x2(double a1[2], double a2[2], double *det);
int   SMfrobenius_norm_squared2x2(double a1[2], double a2[2], double *norma);
int   SMadjoint2x2(double a1[2], double a2[2], double b1[2], double b2[2]);
int   SMmultiply2x2(double a1[2], double a2[2], double b1[2], double b2[2],
		   double r1[2], double r2[2]);
int   SMtranspose2x2(double a1[2], double a2[2], double b1[2], double b2[2]);
int   SMdeterminant3x3(double a1[3], double a2[3], double a3[3], double *det);
int   SMmultiply3x3(double a1[3], double a2[3], double a3[3],
                   double b1[3], double b2[3], double b3[3],
		   double r1[3], double r2[3], double r3[3]);
int   SMtranspose3x3(double a1[3], double a2[3], double a3[3],
		    double b1[3], double b2[3], double b3[3]);
int   SMadjoint3x3(double a1[3], double a2[3], double a3[3],
                  double b1[3], double b2[3], double b3[3]);
int   SMfrobenius_norm_squared3x3(double a1[3], double a2[3],
                  double a3[3], double *norma );
int   SMfrobenius_norm_squared_adj3x3(double a1[3], double a2[3],
                  double a3[3], double *norm_adja );

/* miscellaneous routines */
int   SMinsertQualityInfo(SMquality_table *quality_table, int measure_id,
                        double *function, int num_values);
int   SMvalidityCheck(SMlocal_mesh *local_mesh, SMparam *smooth_param, int *valid);
int   SMvalidMesh(SMlocal_mesh *local_mesh, int *valid);
int   SMorient2D(double *vtx1, double *vtx2, double *vtx3, int *valid);
int   SMorient3D(double *vtx1, double *vtx2, double *vtx3, double *free_vtx,
                 int *valid);

/* write to matlab file routines */
int   SMwriteLocalMesh(FILE *fp, SMlocal_mesh *local_mesh);
int   SMwriteLocalAxes(FILE *fp, SMlocal_mesh *mesh);
int   SMwriteLocalTriangleList(FILE *fp, SMlocal_mesh *local_mesh);
int   SMwriteActiveSet(FILE *fp,SMlocal_mesh *local_mesh);
int   SMwritePoint(FILE *fp, double x, double y);
int   SMwriteSearch(FILE *fp, SMlocal_mesh *local_mesh);

/* free routines */
int   SMfreeOpt(SMoptimal *opt_info);
int   SMfreeLP(SMlocal_mesh *local_mesh, int num_active, int num_constraints);
int   SMfreeActive(SMactive *active);
int   SMfreeParam(SMparam *smooth_param);
int   SMfreeLocalMesh(SMlocal_mesh *local_mesh);
int   SMfreeProcinfo(SMprocinfo *procinfo);
int   SMfreeQualityTable(SMquality_table *quality_table);

/* error routines */
int   SMwrite_ordered_points(SMlocal_mesh *local_mesh);

/* lapack stuff */
#if defined (rs6000)
#define DGESV dgesv
void dgesv(int *n, int *nrhs, double *A, int *lda, int *IPIV, double *B, int *LDB, int *INFO);
#else
#define DGESV dgesv_
void dgesv_(int *n, int *nrhs, double *A, int *lda, int *IPIV, double *B, int *LDB, int *INFO);
#endif

#endif
