/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_EXTERNAL_FUNC_H
#define SM_EXTERNAL_FUNC_H 1

#include "SMuserDefs.h"
#include "SMerror.h"

/* EXTERNAL FUNCTIONS */
/* functions needed in the user code */

#ifdef __cplusplus
extern "C" {
#endif
  /* Initialization Routines */
int         SMinitSmoothing(int argc, char** argv, int dimension,
                     int technique, int FunctionID, double AcceptFunction,
                     void **smooth_data);
int         SMsetProblemDimension(void *smooth_data, int dimension);
int         SMsetSmoothTechnique(void *smooth_data, int technique);
int         SMinitGlobalMinValue(void *smooth_data);
int         SMsetSmoothThreshold(void *smooth_data, double accept);
int         SMsetSmoothFunction(void *smooth_data, int FunctionID);
int         SMsetUntangleTechnique(void *smooth_data, int technique);

  /* Mesh Improvement Routines */
int         SMsmooth(int num_pts, int num_tet, double *free_vtx,
                     double **vtx_list, int **vtx_connectivity,
                     void *smooth_data);

int         SMuntangle(int num_pts, int num_tet, double *free_vtx,
                       double **vtx_list, int **vtx_connectivity,
                       void *smooth_data);

/* void             SMsetUserQualityFunction2D(void *ext_smooth_data,
                                          int values_per_tri,
                                          SMfunction_ptr2D userQualityFunc,
                                          SMgradfunc_ptr2D userQualityGrad);
*/

  /* Statistics Routines */
int           SMinitSmoothStats(void *smooth_data);
int           SMprintSmoothStats(void *smooth_data);

  /* Quality Routines */
int           SMinitQualityTable(void *smooth_data);
int           SMaccumulateQualityInformation(void *smooth_data, double **vtx);
int           SMprintQualityInformation(void *smooth_data);
int           SMisValidMesh(void *smooth_data, int *valid);
int           SMsetMeshValidity(int mesh_validity, void *ext_smooth_data);

  /* Clean-up routines */
int           SMfinalizeSmoothing(void *smooth_data);

#ifdef __cplusplus
}
#endif

#endif
