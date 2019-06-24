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
#define __FUNC__ "SMminmaxOpt"
int SMminmaxOpt(SMlocal_mesh *local_mesh, SMparam *smooth_param, SMprocinfo *procinfo)
{
      int ierr;
#ifdef OPTMS_LOCALTEST
      FILE *fp;
#endif

      local_mesh->opt_done = OPTMS_TRUE;
      SM_LOG_EVENT_BEGIN(__SM_SMOOTH_OPT__);

      /* initialize the optimization structure */
      ierr = SMinitOpt(local_mesh->num_values, local_mesh->opt_info);
             OPTMS_CHKERR(ierr);
      ierr = SMinitMaxStepLength(local_mesh); OPTMS_CHKERR(ierr);

      local_mesh->opt_info->prev_active_values[0]=local_mesh->original_value;
      if ((local_mesh->lap_done == 0) ||
          (local_mesh->lap_info->lap_accepted == OPTMS_FALSE)) {
          OPTMS_COPY_VECTOR(local_mesh->opt_info->function,
                   local_mesh->original_function,local_mesh->num_values);
          local_mesh->opt_info->iter_count = 0;
      } else {
          OPTMS_COPY_VECTOR(local_mesh->opt_info->function,
                   local_mesh->lap_info->laplacian_function,
		   local_mesh->num_values);
          local_mesh->opt_info->prev_active_values[1]=
                                local_mesh->lap_info->laplacian_value;
          local_mesh->opt_info->iter_count = 1;
      }
      ierr = SMfindActiveSet(local_mesh->num_values,local_mesh->opt_info->function,
		      smooth_param->active_eps,local_mesh->opt_info->active);
             OPTMS_CHKERR(ierr);
      OPTMS_DEBUG_ACTION(3,{
 	 /* Print the active set */
 	 ierr = SMprintActiveSet(local_mesh->opt_info->active,
			 local_mesh->opt_info->function);
                OPTMS_CHKERR(ierr);
      });

     /* check for equilibrium point */
     /* compute the gradient */
     ierr = SMcomputeGradient(local_mesh,smooth_param,
                      local_mesh->opt_info->gradient); OPTMS_CHKERR(ierr);
     
     if (local_mesh->opt_info->active->num_active >= 2) {
	OPTMS_DEBUG_PRINT(3,"Testing for an equilibrium point \n");
	ierr = SMcheckEquilibrium(local_mesh->opt_info, 
               &local_mesh->opt_info->equilibrium_pt); OPTMS_CHKERR(ierr);

	OPTMS_DEBUG_ACTION(2,{
	    if (local_mesh->opt_info->equilibrium_pt) 
		fprintf(stdout,"Optimization Exiting: An equilibrium point \n");
        });
     }

    /* terminate if we have found an equilibrium point or if the step is
       too small to be worthwhile continuing */
    while ((local_mesh->opt_info->status != OPTMS_EQUILIBRIUM) && 
	   (local_mesh->opt_info->status != OPTMS_STEP_TOO_SMALL) &&
	   (local_mesh->opt_info->status != OPTMS_IMP_TOO_SMALL) &&
	   (local_mesh->opt_info->status != OPTMS_FLAT_NO_IMP) &&
           (local_mesh->opt_info->status != OPTMS_ZERO_SEARCH) &&
	   (local_mesh->opt_info->status != OPTMS_MAX_ITER_EXCEEDED)) {

	/* increase the iteration count by one */
        /* smooth_param->iter_count += 1; */
        local_mesh->opt_info->iter_count += 1;
        local_mesh->opt_info->opt_iter_count += 1;
        if (local_mesh->opt_info->iter_count > OPTMS_MAX_OPT_ITER)
	   local_mesh->opt_info->status = OPTMS_MAX_ITER_EXCEEDED;

	OPTMS_DEBUG_PRINT(3,"\n");
	OPTMS_DEBUG_ACTION(3,{ 
            fprintf(stdout,"ITERATION %d \n",local_mesh->opt_info->iter_count);
        });
	    
	/* compute the gradient */
	ierr = SMcomputeGradient(local_mesh,smooth_param,
                          local_mesh->opt_info->gradient); OPTMS_CHKERR(ierr);
        
	OPTMS_DEBUG_PRINT(3,"computing the search direction \n");
	ierr = SMsearchDirection(local_mesh); OPTMS_CHKERR(ierr);

	OPTMS_MATLAB_ON({
	    if (OPTMS_ISROOT(procinfo) && (dimension==2)) {
		/* plot out the local mesh */
	        OPTMS_DEBUG_PRINT(2,"Plotting the search \n");
		if ((fp = fopen("search.m","w")) == NULL) {
                    OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open search.m for writing\n");
                }
		ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
		ierr = SMwriteActiveSet(fp,local_mesh); OPTMS_CHKERR(ierr);  
		ierr = SMwriteSearch(fp,local_mesh); OPTMS_CHKERR(ierr);
		fclose(fp);
	    }
	OPTMS_MATLAB_OFF});

	/* if there are viable directions to search */
	if ((local_mesh->opt_info->status != OPTMS_ZERO_SEARCH) &&
            (local_mesh->opt_info->status != OPTMS_MAX_ITER_EXCEEDED)) {

	    OPTMS_DEBUG_PRINT(3,"Computing the projections of the gradients \n");
	    ierr = SMgetGradientProjections(local_mesh->opt_info);
                   OPTMS_CHKERR(ierr);

	    OPTMS_DEBUG_PRINT(3,"Computing the initial step size \n");
	    ierr = SMcomputeAlpha(local_mesh->opt_info); OPTMS_CHKERR(ierr);

	    OPTMS_DEBUG_PRINT(3,"Testing whether to accept this step \n");
	    ierr = SMstepAcceptance(local_mesh,smooth_param); OPTMS_CHKERR(ierr);
            OPTMS_DEBUG_ACTION(3,
              {printf("The new free vertex position is %f %f %f\n",
              local_mesh->free_vtx[0], local_mesh->free_vtx[1], 
		       local_mesh->free_vtx[2]);});

	    OPTMS_DEBUG_ACTION(3,{
     		/* Print the active set */
	     	ierr = SMprintActiveSet(local_mesh->opt_info->active,
				 local_mesh->opt_info->function);
                       OPTMS_CHKERR(ierr);
	    });

	    /* check for equilibrium point */
	    if (local_mesh->opt_info->active->num_active >= 2) {
		OPTMS_DEBUG_PRINT(3,"Testing for an equilibrium point \n");
                ierr = SMcheckEquilibrium(local_mesh->opt_info,
   		        &local_mesh->opt_info->equilibrium_pt); 
                       OPTMS_CHKERR(ierr);

		OPTMS_DEBUG_ACTION(2,{
		    if (local_mesh->opt_info->equilibrium_pt) 
			fprintf(stdout,"Optimization Exiting: An equilibrium point \n");
                });
	    }

	    /* record the values */
            local_mesh->current_active_value = 
                 local_mesh->opt_info->active->true_active_value;
	    OPTMS_RECORD_ITER_VALUE(local_mesh->opt_info);
	} else {
	    /* decrease the iteration count by one */
	    /* smooth_param->iter_count -= 1; */
	    local_mesh->opt_info->iter_count -= 1;
	    OPTMS_DEBUG_ACTION(2,{
		fprintf(stdout,"Optimization Exiting: No viable directions; equilibrium point \n");
		/* Print the old active set */
		ierr = SMprintActiveSet(local_mesh->opt_info->active,
			     local_mesh->opt_info->function);
                OPTMS_CHKERR(ierr);
	    });
	}
      }
      OPTMS_ASSERT_ON({
      OPTMS_DEBUG_PRINT(3,"Checking the validity of the mesh\n");
            ierr = SMvalidityCheck(local_mesh,smooth_param,&valid); OPTMS_CHKERR(ierr);
	    if (!valid) fprintf(stdout,"The final mesh is not valid\n");
      OPTMS_ASSERT_OFF});
      SM_LOG_EVENT_END(__SM_SMOOTH_OPT__);
      OPTMS_DEBUG_ACTION(2,{fprintf(stdout,"Number of optimization iterations %d\n",
                            local_mesh->opt_info->iter_count);});

/* added by DPS--assume that if we have reached this point without throwing an
   error that we are OK. */
      return(ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMstepAcceptance"
int SMstepAcceptance(SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
  int        ierr;
  int        step_accepted, i;
  int        num_values, num_steps;
  int        valid, step_status;
  int        accept_alpha;
  int        dimension;
  double     alpha;
  double     estimated_improvement;
  double     current_improvement = OPTMS_BIG_NEG_NMBR;
  double     previous_improvement = OPTMS_BIG_NEG_NMBR;
  double     current_percent_diff = OPTMS_BIG_POS_NMBR;
  double     original_point[OPTMS_MAX_DIM];
  SMactive   *active;
  SMoptimal  *opt_info;

  OPTMS_CHECK_NULL(local_mesh);
  OPTMS_CHECK_NULL(smooth_param);
    
  opt_info = local_mesh->opt_info;
  num_values = opt_info->num_values;
  active = opt_info->active;
  alpha = opt_info->alpha;

  dimension = local_mesh->dimension;

  step_accepted = 0;
  step_status = OPTMS_STEP_NOT_DONE;
  opt_info->status = 0;
  num_steps = 0;

  if (alpha < smooth_param->min_step_size) {
      opt_info->status = OPTMS_IMP_TOO_SMALL;
      step_status = OPTMS_STEP_DONE;
      OPTMS_DEBUG_PRINT(3,"Alpha starts too small, no improvement\n");
  }

  /* save the original function and active set */
  OPTMS_COPY_VECTOR(original_point,local_mesh->free_vtx,dimension);
  OPTMS_COPY_VECTOR(opt_info->original_function, opt_info->function,num_values);
  ierr = SMcopyActive(opt_info->active, opt_info->original_active); OPTMS_CHKERR(ierr);

  while (step_status == OPTMS_STEP_NOT_DONE) {

    num_steps++;  if (num_steps >= 100) step_status = OPTMS_STEP_DONE;

    accept_alpha = OPTMS_FALSE;
    while (!accept_alpha && alpha>smooth_param->min_step_size) {

      /* make the step */
      for (i=0;i<dimension;i++) {
         local_mesh->free_vtx[i] += alpha*opt_info->search[i];
      }

      /* assume alpha is acceptable */
      accept_alpha=OPTMS_TRUE;

      /* check for valid step */
      valid = 1;
      ierr = SMvalidMesh(local_mesh, &valid); OPTMS_CHKERR(ierr);

      /* never take a step that makes a valid mesh invalid */
      if (!valid && (local_mesh->validity==OPTMS_VALID_MESH)) {
            accept_alpha=OPTMS_FALSE;
            for (i=0;i<dimension;i++) {
               local_mesh->free_vtx[i] -= alpha*opt_info->search[i];
            }
            alpha = alpha/2;
            opt_info->alpha = alpha;
             OPTMS_DEBUG_ACTION(2,{
                 fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); 
             });

          if (alpha < smooth_param->min_step_size) {
 	        opt_info->status = OPTMS_STEP_TOO_SMALL;
                step_status = OPTMS_STEP_DONE;
                OPTMS_DEBUG_PRINT(2,"Step too small\n");
 	        /* get back the original point, function, and active set */
                OPTMS_COPY_VECTOR(local_mesh->free_vtx,original_point,dimension);
	        OPTMS_COPY_VECTOR(opt_info->function,opt_info->original_function,num_values);
	        ierr = SMcopyActive(opt_info->original_active, opt_info->active); 
                       OPTMS_CHKERR(ierr);
	  }
       }
    }     

    if ((valid || local_mesh->validity==OPTMS_INVALID_MESH) && 
        (alpha > smooth_param->min_step_size)) {
      /* compute the new function and active set */
      ierr = SMcomputeFunction(local_mesh,smooth_param,opt_info->function); 
             OPTMS_CHKERR(ierr);
      ierr = SMfindActiveSet(local_mesh->opt_info->num_values,opt_info->function,
		      smooth_param->active_eps,opt_info->active);
             OPTMS_CHKERR(ierr);
	
      /* estimate the minimum improvement by taking this step */
      ierr = SMgetMinEstimate(local_mesh->opt_info, &estimated_improvement);
             OPTMS_CHKERR(ierr);
      OPTMS_DEBUG_ACTION(3,{
           fprintf(stdout,"The estimated improvement for this step: %f\n",
		   estimated_improvement); 
      });
	
      /* calculate the actual increase */
      current_improvement = 
	opt_info->active->true_active_value -
	opt_info->prev_active_values[opt_info->iter_count-1];

      OPTMS_DEBUG_ACTION(3,{fprintf(stdout,"Actual improvement %f\n",current_improvement);});

      /* calculate the percent difference from estimated increase */
      current_percent_diff = fabs(current_improvement-estimated_improvement)/
	fabs(estimated_improvement);

      /* determine whether to accept a step */
      if ((previous_improvement > current_improvement) && 
	  (previous_improvement > 0)) {
	/* accept the previous step - it was better */
	     OPTMS_DEBUG_PRINT(2,"Accepting the previous step\n");
 
	/* add alpha in again (previous step) */
	for (i=0;i<dimension;i++) {
	  local_mesh->free_vtx[i] += alpha*opt_info->search[i];
	}

	/* does this make an invalid mesh valid? */
              valid = 1;
              ierr = SMvalidMesh(local_mesh, &valid); OPTMS_CHKERR(ierr);
              if (valid && local_mesh->validity==OPTMS_INVALID_MESH) {
                    local_mesh->validity=OPTMS_VALID_MESH;
              }

	/* copy test function and active set */
	OPTMS_COPY_VECTOR(opt_info->function,opt_info->test_function,num_values);
	ierr = SMcopyActive(opt_info->test_active, opt_info->active); OPTMS_CHKERR(ierr);

	opt_info->status = OPTMS_STEP_ACCEPTED;  step_status = OPTMS_STEP_DONE;
            
	/* check to see that we're still making good improvements */
	if (previous_improvement < smooth_param->min_acceptable_imp) {
	  opt_info->status = OPTMS_IMP_TOO_SMALL; step_status = OPTMS_STEP_DONE;
	  OPTMS_DEBUG_PRINT(2,"Optimization Exiting: Improvement too small\n");
	}

      } else if (((current_improvement > estimated_improvement) ||
		  (current_percent_diff < .1)) && (current_improvement>0)) {
	/* accept this step, exceeded estimated increase or was close */
	opt_info->status = OPTMS_STEP_ACCEPTED;  step_status = OPTMS_STEP_DONE;

	/* does this make an invalid mesh valid? */
              valid = 1;
              ierr = SMvalidMesh(local_mesh, &valid); OPTMS_CHKERR(ierr);
              if (valid && local_mesh->validity==OPTMS_INVALID_MESH) {
                    local_mesh->validity=OPTMS_VALID_MESH;
              }

            
	/* check to see that we're still making good improvements */
	if (current_improvement < smooth_param->min_acceptable_imp) {
	  OPTMS_DEBUG_PRINT(2,"Optimization Exiting: Improvement too small\n");
	  opt_info->status = OPTMS_IMP_TOO_SMALL; step_status = OPTMS_STEP_DONE;
	}
      } else if ((current_improvement < 0) && (previous_improvement < 0) &&
		 (fabs(current_improvement) < smooth_param->min_acceptable_imp) &&
		 (fabs(previous_improvement) < smooth_param->min_acceptable_imp)) {

	/* we are making no progress, quit */
	opt_info->status = OPTMS_FLAT_NO_IMP; step_status = OPTMS_STEP_DONE;
	OPTMS_DEBUG_PRINT(2,"Opimization Exiting: Flat no improvement\n");
           
	/* get back the original point, function, and active set */
	OPTMS_COPY_VECTOR(local_mesh->free_vtx,original_point,dimension);
	OPTMS_COPY_VECTOR(opt_info->function,opt_info->original_function,num_values);
	ierr = SMcopyActive(opt_info->original_active, opt_info->active); OPTMS_CHKERR(ierr);

      } else {
	/* halve alpha and try again */
	/* subtract out the old step */
	for (i=0;i<dimension;i++) 
                  local_mesh->free_vtx[i] -= alpha*opt_info->search[i];

	/* halve step size */
	alpha = alpha/2; 
	opt_info->alpha = alpha;
	OPTMS_DEBUG_ACTION(3,{fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); });

	if (alpha < smooth_param->min_step_size) {
	  /* get back the original point, function, and active set */
	  OPTMS_DEBUG_PRINT(2,"Optimization Exiting: Step too small\n");
	  OPTMS_COPY_VECTOR(local_mesh->free_vtx,original_point,dimension);
	  OPTMS_COPY_VECTOR(opt_info->function,opt_info->original_function,num_values);
	  ierr = SMcopyActive(opt_info->original_active, opt_info->active); OPTMS_CHKERR(ierr);
	  opt_info->status = OPTMS_STEP_TOO_SMALL;  step_status = OPTMS_STEP_DONE;
	} else {
	  OPTMS_COPY_VECTOR(opt_info->test_function, opt_info->function,num_values);
	  ierr = SMcopyActive(opt_info->active,opt_info->test_active); OPTMS_CHKERR(ierr);
	  previous_improvement = current_improvement;
	}
      }
    }
  }
  if (current_improvement<0 && opt_info->status==OPTMS_STEP_ACCEPTED) {
    OPTMS_DEBUG_ACTION(2,{printf("Accepted a negative step %f \n",current_improvement);
                       ierr = SMwrite_ordered_points(local_mesh); 
                       OPTMS_CHKERR(ierr);});
  }
  return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcheckEquilibrium"
int SMcheckEquilibrium(SMoptimal *opt_info, int *equil)
{
    int  ierr;
    int  i,j;
    int  num_active;
    int  ind1, ind2;  
    int dimension;
    double min;
    double **dir, **G;
    double mid_vec[OPTMS_MAX_DIM], mid_cos, test_cos;

    dimension = opt_info->dimension;

    num_active = opt_info->active->num_active;
    *equil = OPTMS_FALSE;
    ind1 = ind2 = -1;

    if (num_active==opt_info->num_values) {
         *equil = 1; opt_info->status = OPTMS_EQUILIBRIUM;
         OPTMS_DEBUG_PRINT(3,"All the function values are in the active set\n"); 
         return (ierr=0);
    }

    /* set up the active gradient directions */
    ierr = SMgetActiveDirections(num_active,opt_info->gradient,
                   opt_info->active->active_ind,dimension,&dir);
           OPTMS_CHKERR(ierr);

    /* normalize the active directions */
    for (j=0;j<num_active;j++) OPTMS_NORMALIZE(dir[j],dimension);

    if (dimension == 2) {
    /* form the grammian */
    ierr = SMformGrammian(num_active,dir,opt_info->G,dimension);  
           OPTMS_CHKERR(ierr);
    G = opt_info->G;

    /* find the minimum element in the upper triangular portion of G*/
    min = 1;
    for (i=0;i<num_active;i++) {
      for (j=i+1;j<num_active;j++) {
        if ( fabs(-1 - G[i][j]) < 1E-08 ) {
           *equil = 1; opt_info->status = OPTMS_EQUILIBRIUM;
           OPTMS_DEBUG_PRINT(3,"The gradients are antiparallel, eq. pt\n"); 
           
           /*DPS added */
           for (i=0;i<num_active;i++) OPTMS_FREE(dir[i]);
	   OPTMS_FREE(dir);

	   /* end DPS added */
           return (ierr=0);
         }
         if (G[i][j]  < min) {
           ind1 = i; ind2 = j;
           min = G[i][j];
        }
      }
    }

    if ((ind1 != -1) && (ind2 != -1)) {
      /* find the diagonal of the parallelepiped */
      for (j=0;j<dimension;j++) {
       mid_vec[j]=.5*(dir[ind1][j]+dir[ind2][j]);
      }
      OPTMS_NORMALIZE(mid_vec,dimension);
      OPTMS_DOT(mid_cos,dir[ind1],mid_vec,dimension);

      /* test the other vectors to be sure they lie in the cone */
      for (i=0;i<num_active;i++) {
         if ((i != ind1) && (i != ind2)) {
            OPTMS_DOT(test_cos,dir[i],mid_vec,dimension);
            if ((test_cos < mid_cos)  &&  fabs(test_cos-mid_cos) > OPTMS_MACHINE_EPS) {
              OPTMS_DEBUG_PRINT(3,"An equilibrium point \n");
              *equil = 1; opt_info->status = OPTMS_EQUILIBRIUM;
	      /*DPS added */
	      for (i=0;i<num_active;i++) OPTMS_FREE(dir[i]);
	      OPTMS_FREE(dir);

	      /* end DPS added */
  
              return (ierr=0);
            }
         }
       }
     }
    }
    if (dimension == 3) {
       ierr = SMconvexHullTest(dir,num_active,equil); OPTMS_CHKERR(ierr);
       if (*equil == 1) opt_info->status = OPTMS_EQUILIBRIUM;
    }
    /*DPS added */
    for (i=0;i<num_active;i++) OPTMS_FREE(dir[i]);
    OPTMS_FREE(dir);
    
    /* end DPS added */
  
    return (ierr=0);
}


#undef __FUNC__
#define __FUNC__ "SMgetMinEstimate"
int SMgetMinEstimate(SMoptimal *opt_info, double *final_est)
{
    int ierr;
    int    i;
    int    num_active;
    double alpha;
    double est_imp;

    num_active = opt_info->active->num_active;
    alpha = opt_info->alpha;
    
    *final_est = OPTMS_BIG_POS_NMBR;
    for (i=0;i<num_active;i++) {
	est_imp = alpha*opt_info->gs[opt_info->active->active_ind[i]];
        if (est_imp<*final_est) *final_est = est_imp;
    }
    if (*final_est == 0) {
	*final_est = OPTMS_BIG_POS_NMBR;
	for (i=0;i<opt_info->num_values;i++) {
	    est_imp = alpha*opt_info->gs[i];
	    if ((est_imp<*final_est) && (fabs(est_imp) > OPTMS_MACHINE_EPS)) {
		*final_est = est_imp;
	    }
	}
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeVerticalStep"
int SMcomputeVerticalStep(SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
    int       ierr;
    int       i;
    int       *PDG_ind;
    int       num_LI;   
    int       dimension;           
    double    **N, *R, *y, y_0;
    double    v[OPTMS_MAX_DIM];
    double    **gradient;
    SMoptimal *opt_info;

    opt_info = local_mesh->opt_info;
    num_LI = opt_info->num_LI;
    gradient = opt_info->gradient;
    PDG_ind = opt_info->PDG_ind;

    dimension = local_mesh->dimension;

    /* form the projection matrix */
    ierr = SMformVerticalMatrix(num_LI, opt_info->PDG, &N);
           OPTMS_CHKERR(ierr);

    /* form the RHS */
    ierr = SMformVerticalRHS(num_LI, opt_info->function, 
                          opt_info->PDG_ind, 
                          opt_info->active->true_active_value, &R);
           OPTMS_CHKERR(ierr);

    /* determine the step */
    switch(num_LI) {
    case 1:
        y_0 = R[0]/N[0][0];
        for (i=0;i<dimension;i++) {
            v[i] = gradient[PDG_ind[0]][i]*y_0;
        }
        break;
    case 2:
        ierr = SMsolve2x2(N[0][0],N[0][1],N[1][0],N[1][1],R[0],R[1],&y);
               OPTMS_CHKERR(ierr);
        for (i=0;i<dimension;i++) {
            v[i] = gradient[PDG_ind[0]][i]*y[0] +
                   gradient[PDG_ind[1]][i]*y[1];
	    }
        OPTMS_FREE(y);
        break;
    case 3:
        ierr =  SMsolveSPD3x3(N,R,&y); OPTMS_CHKERR(ierr);
        for (i=0;i<dimension;i++) {
            v[i] = gradient[PDG_ind[0]][i]*y[0] +
                   gradient[PDG_ind[1]][i]*y[1] +
                   gradient[PDG_ind[2]][i]*y[2];
	    }
        OPTMS_FREE(y);
    }
    
    /* add in the correction */
    for (i=0;i<dimension;i++) {
       local_mesh->free_vtx[i] += v[i]; 
    }

    /* free the variables */
    for(i=0;i<num_LI;i++) {OPTMS_FREE(N[i]);} OPTMS_FREE(N);
    OPTMS_FREE(R); 
     
    /* update the function and gradient */
    OPTMS_DEBUG_ACTION(2,{
    ierr = SMcomputeFunction(local_mesh,smooth_param,local_mesh->opt_info->test_function);
           OPTMS_CHKERR(ierr);
    ierr = SMfindActiveSet(local_mesh->opt_info->num_values,
		  local_mesh->opt_info->test_function,
		  smooth_param->active_eps,
		  local_mesh->opt_info->test_active); 
           OPTMS_CHKERR(ierr);
    ierr = SMprintActiveSet(local_mesh->opt_info->test_active,
	             local_mesh->opt_info->test_function);
           OPTMS_CHKERR(ierr);
    });
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMformVerticalRHS"
int SMformVerticalRHS(int n,double *function, int ind[OPTMS_MAX_DIM], 
                          double max_value, double **R)
{
    int ierr;
    int i;
  
    OPTMS_MALLOC((*R),(double *),sizeof(double)*n,1);
    for (i=0;i<n;i++) {
       (*R)[i] = max_value - function[ind[i]];
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMstepToCusp"
int SMstepToCusp(SMlocal_mesh *local_mesh, SMparam *smooth_param,
                  SMprocinfo *procinfo)
{
    int     ierr;
    int     ind0, ind1, i;
    double  f0,f1,proj0,proj1,denom;
    double  alpha_cusp, s_perp[OPTMS_MAX_DIM];
    SMoptimal  *opt_info;
#ifdef OPTMS_LOCALTEST
    FILE *fp;
#endif

    opt_info = local_mesh->opt_info;
    s_perp[OPTMS_XDIR] = -1.*opt_info->search[OPTMS_YDIR];
    s_perp[OPTMS_YDIR] = opt_info->search[OPTMS_XDIR];

    ind0 = opt_info->active->active_ind[0];
    ind1 = opt_info->active->active_ind[1];

    f0 = opt_info->function[ind0];
    f1 = opt_info->function[ind1];
    OPTMS_DOT(proj0,opt_info->gradient[ind0],s_perp,local_mesh->dimension);
    OPTMS_DOT(proj1,opt_info->gradient[ind1],s_perp,local_mesh->dimension);
    denom = proj0 - proj1;
    if (!OPTMS_LESS_THAN_MACHINE_EPS(denom)) {
        alpha_cusp = (f1 - f0)/denom;
        for (i=0;i<local_mesh->dimension;i++) {
           local_mesh->free_vtx[i] += alpha_cusp*s_perp[i];
        }
    }
	

    /*
    if (opt_info->function[ind0] == opt_info->active->true_active_value) {
	f_min      = opt_info->function[ind0];
	f_other    = opt_info->function[ind1];
	OPTMS_DOT(proj_min,opt_info->gradient[ind0],s_perp,opt_info->dimension);
	OPTMS_DOT(proj_other,opt_info->gradient[ind1],s_perp,opt_info->dimension);
    } else {
	f_min      = opt_info->function[ind1];
	f_other    = opt_info->function[ind0];
	OPTMS_DOT(proj_min,opt_info->gradient[ind1],s_perp,opt_info->dimension);
	OPTMS_DOT(proj_other,opt_info->gradient[ind0],s_perp,opt_info->dimension);
    }
    
    if ( fabs(proj_min - proj_other) > OPTMS_MACHINE_EPS) {
	alpha_cusp = (f_min - f_other)/(proj_other - proj_min);

        for (i=0;i<opt_info->dimension;i++) {
    	    local_mesh->free_vtx[i] += alpha_cusp*s_perp[i];
        }
    }
    */
    OPTMS_DEBUG_ACTION(2,{
        ierr = SMcomputeFunction(local_mesh,smooth_param,local_mesh->opt_info->test_function);
               OPTMS_CHKERR(ierr);
        ierr = SMfindActiveSet(local_mesh->opt_info->num_values,
         		local_mesh->opt_info->test_function,
		        smooth_param->active_eps,
		        local_mesh->opt_info->test_active);      
               OPTMS_CHKERR(ierr);
        ierr = SMprintActiveSet(local_mesh->opt_info->test_active,
	                 local_mesh->opt_info->test_function);
               OPTMS_CHKERR(ierr);
    });

    OPTMS_MATLAB_ON({
	if (OPTMS_ISROOT(procinfo) && (local_mesh->dimension==2)) {
	    fprintf(stdout,"Plotting the cusp step \n");
	    if ((fp = fopen("cusp.m","w")) == NULL) {
                OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open cusp.m for writing\n");
            }
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
    OPTMS_MATLAB_OFF});
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMgetGradientProjections"
int SMgetGradientProjections(SMoptimal *opt_info)
{
    int ierr;
    int i;

    OPTMS_CHECK_NULL(opt_info);
    for (i=0;i<opt_info->num_values;i++) {
	OPTMS_DOT(opt_info->gs[i],opt_info->gradient[i],opt_info->search,opt_info->dimension);
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeAlpha"
int SMcomputeAlpha(SMoptimal *opt_info)
{
    int       ierr;
    int       i, j, steepest;
    int       ind;
    int       num_values;
    int       *active;
    double    *function;
    double    steepest_function;
    double    steepest_grad;
    double    alpha_i;
    double    alpha;

    OPTMS_CHECK_NULL(opt_info);

    num_values = opt_info->num_values;
    alpha = OPTMS_BIG_POS_NMBR;

    OPTMS_MALLOC(active,(int *),sizeof(int)*num_values,1);
    
    for (i=0;i<num_values;i++) {
        active[i] = 0;
    }
    for (j=0;j<opt_info->active->num_active;j++) {
        ind = opt_info->active->active_ind[j];
        active[ind] = 1;
    }
    
    function = opt_info->function;
    steepest = opt_info->steepest;
    steepest_function = function[steepest];
    steepest_grad = opt_info->gs[steepest];
    for (i=0;i<num_values;i++) {
        /* if it's not active */
      if (i!=steepest) {
	    alpha_i = function[i]-steepest_function;
	   
	    if (fabs(opt_info->gs[steepest] - opt_info->gs[i])>1E-13) {
	       /* compute line intersection */
	       alpha_i = alpha_i/(steepest_grad - opt_info->gs[i]);
	    } else {
	       /* the lines don't intersect - it's not under consideration*/
	       alpha_i = 0;
	    }
	    if ((alpha_i > 0 ) && (fabs(alpha_i) < fabs(alpha))) {
	      alpha = fabs(alpha_i); 
	    }
	    
      }
    }

    /* if it never gets set, set it to the default */
    if (alpha == OPTMS_BIG_POS_NMBR) {
      alpha = opt_info->max_alpha;
      OPTMS_DEBUG_ACTION(3,{ fprintf(stdout,"Setting alpha to the maximum step length\n"); });
    }

    opt_info->alpha = alpha;

    OPTMS_FREE(active);
    OPTMS_DEBUG_ACTION(3,{ fprintf(stdout,"The initial step size: %f\n",alpha); });
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcopyActive"
int SMcopyActive(SMactive *active1, SMactive *active2)
{
    int ierr;
    int i;

    OPTMS_CHECK_NULL(active1);
    OPTMS_CHECK_NULL(active2);

    active2->num_active = active1->num_active;
    active2->num_equal  = active1->num_equal;
    active2->true_active_value = active1->true_active_value;
    for (i=0;i<active1->num_active;i++) {
	active2->active_ind[i] = active1->active_ind[i];
    }
    return(ierr=0);
}


