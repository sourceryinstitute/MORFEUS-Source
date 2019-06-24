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
#define __FUNC__ "SMcomputeFunction"
int  SMcomputeFunction(SMlocal_mesh *local_mesh, SMparam *smooth_param,
                        double *function)
{
    int ierr = 0;
    int fcn_ind, ind1, ind2, ind3;
    int i,j,num_values,num_tri;
    int dimension;
    double *function1;
    SMfunction_ptr2D ComputeFunctionValues2D;
    SMfunction_ptr3D ComputeFunctionValues3D;

    OPTMS_CHECK_NULL(local_mesh);
    OPTMS_CHECK_NULL(smooth_param);
    OPTMS_CHECK_NULL(function);

    dimension = local_mesh->dimension;

    SM_LOG_EVENT_BEGIN(__SM_FUNCTION__);
    if (dimension==2) {
      ComputeFunctionValues2D = smooth_param->ComputeFunctionValues2D;
      fcn_ind=0; ind1=0; ind2=0;
      OPTMS_MALLOC(function1,(double *),sizeof(double)*smooth_param->function_values_per_tri,1);

      num_tri = local_mesh->num_tri;
      for (i=0;i<num_tri;i++) {
         ind1=local_mesh->vtx_connectivity[i][0];
         ind2=local_mesh->vtx_connectivity[i][1];
         ierr = ComputeFunctionValues2D(local_mesh->free_vtx,
                     local_mesh->incident_vtx[ind1],
                     local_mesh->incident_vtx[ind2],
                     function1,&num_values);  OPTMS_CHKERR(ierr);
          for (j=0;j<num_values;j++)  function[fcn_ind++] = function1[j];
       }
    } else if (dimension ==3) {
      ComputeFunctionValues3D = smooth_param->ComputeFunctionValues3D;
      fcn_ind=0; ind1=0; ind2=0;
      OPTMS_MALLOC(function1,(double *),
                   sizeof(double)*smooth_param->function_values_per_tri,1);

      num_tri = local_mesh->num_tri;
      for (i=0;i<num_tri;i++) {
         ind1 = local_mesh->vtx_connectivity[i][0];
         ind2 = local_mesh->vtx_connectivity[i][1];
         ind3 = local_mesh->vtx_connectivity[i][2];
         ComputeFunctionValues3D(local_mesh->free_vtx,
                   local_mesh->incident_vtx[ind1],
                   local_mesh->incident_vtx[ind2],
                   local_mesh->incident_vtx[ind3],
                   function1,&num_values);
         for (j=0;j<num_values;j++)  function[fcn_ind++] = function1[j];
      }
    } else {
      OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension must be 2 or 3");
    }

    /* free the temporary function memory */
    OPTMS_FREE(function1);
    SM_LOG_EVENT_END(__SM_FUNCTION__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeGradient"
int SMcomputeGradient(SMlocal_mesh *local_mesh, SMparam *smooth_param,
                       double **gradient)
{
    int ierr;
    int i,j,num_values, num_tri;
    int ind1, ind2, ind3, grad_ind;
    int dimension;
    int num_val_per_tri;
    double **gradient1;
    SMgradfunc_ptr2D ComputeGradientValues2D;
    SMgradfunc_ptr3D ComputeGradientValues3D;

    OPTMS_CHECK_NULL(local_mesh);
    OPTMS_CHECK_NULL(smooth_param);
    OPTMS_CHECK_NULL(gradient);

    dimension = local_mesh->dimension;

    SM_LOG_EVENT_BEGIN(__SM_GRADIENT__);
    if (dimension == 2) {
      num_tri = local_mesh->num_tri;
      num_val_per_tri = 3;
      ComputeGradientValues2D = smooth_param->ComputeGradientValues2D;
      OPTMS_MALLOC(gradient1,(double **),sizeof(double *)*num_val_per_tri,1);
      for (i=0;i<num_val_per_tri;i++)
           OPTMS_MALLOC(gradient1[i],(double *),sizeof(double)*dimension,1);

      ind1=0; ind2=0; grad_ind = 0;
      for (i=0;i<num_tri;i++) {
            ind1=local_mesh->vtx_connectivity[i][0];
            ind2=local_mesh->vtx_connectivity[i][1];
            ierr = ComputeGradientValues2D(
                                  local_mesh->free_vtx,
                                  local_mesh->incident_vtx[ind1],
                                  local_mesh->incident_vtx[ind2],
                                  gradient1,&num_values);
                   OPTMS_CHKERR(ierr);
            for (j=0;j<num_values;j++) {
                OPTMS_COPY_VECTOR(gradient[grad_ind],gradient1[j],2);
                grad_ind++;
            }
        }
        for (i=0;i<num_val_per_tri;i++) OPTMS_FREE(gradient1[i]);
        OPTMS_FREE(gradient1);
    } else if (dimension == 3) {
        num_val_per_tri = 6;
        ComputeGradientValues3D = smooth_param->ComputeGradientValues3D;
        OPTMS_MALLOC(gradient1,(double **),sizeof(double *)*num_val_per_tri,1);
        for (i=0;i<num_val_per_tri;i++) {
            OPTMS_MALLOC(gradient1[i],(double *),sizeof(double)*dimension,1);
        }

        ind1=0; ind2=0; ind3 = 0; grad_ind=0;
        num_tri = local_mesh->num_tri;
        for (i=0;i<num_tri;i++) {
            ind1 = local_mesh->vtx_connectivity[i][0];
            ind2 = local_mesh->vtx_connectivity[i][1];
            ind3 = local_mesh->vtx_connectivity[i][2];
            ComputeGradientValues3D(local_mesh->free_vtx,
                       local_mesh->incident_vtx[ind1],
                       local_mesh->incident_vtx[ind2],
                       local_mesh->incident_vtx[ind3],
                       gradient1,&num_values);
            for (j=0;j<num_values;j++) {
                OPTMS_COPY_VECTOR(gradient[grad_ind],gradient1[j],dimension);
                grad_ind++;
            }
        }
        for (i=0;i<num_val_per_tri ;i++) OPTMS_FREE(gradient1[i]); /* was <3 */
        OPTMS_FREE(gradient1);
    } else {
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension must be 2 or 3");
    }
    SM_LOG_EVENT_END(__SM_GRADIENT__);
    return(ierr=0);
}
