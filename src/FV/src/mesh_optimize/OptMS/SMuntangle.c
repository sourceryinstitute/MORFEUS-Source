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

/*@
   SMuntangle - This is the main routine that attempts to untangle a local submesh by 
      maximizing the minimum area of a triangle or tetrahedra in a local submesh.  The
      value of the function assumes a right-hand ordering of the vertices.  If the user is not
      sure if the mesh is invalid, a quality assessment should be done to compute the
      Jacobians of each of the mesh elements.  Once this information is accumulated,
      the routine SMinvalidMesh returns OPTMS_TRUE (1) if the mesh contains inverted elements.

   Input Parameters:
+  num_incident_vtx - the number of incident vertices in the local mesh 
.  num_tri - the number of incident triangles or tetrahedra
.  free_vtx - the coordinates of the free vertex in a 
           vector of length problem dimensions in x, y, z order
.  vtx_list - a matrix of the coordinates of the incident vtx;
           matrix dimensions are num_incident_vtx by problem dimension
.  vtx_connectivity - a matrix that gives the connectivity info for the
           incident vertices. matrix dimensions are num_incident_vtx by 
           the problem dimension.  Note: this assumes that the connectivity given 
           is for a right handed triangle or tetrahedra with the free vertex 
           ordered first
-  ext_smooth_data - data structure for mesh smoothing; created in 
          SMinitSmoothing and cast to the local data structure

   Output Parameter:
.  free_vtx - contains the new coordinates of the free vertex
 
   Note:
   This function can only be called after SMinitSmoothing has been called to create
    the smooth_data data structure.  Once the mesh has been untangled, SMsmooth
    should be called to improve the element quality as it is usually very poor at the conclusion
    of this step.

.seealso SMinitSmoothing(), SMinitQualityTable, SMinvalidMesh(), SMfinializeSmoothing()
@*/
#undef __FUNC__
#define __FUNC__ "SMuntangle"
int SMuntangle(int num_incident_vtx, int num_tri, double *free_vtx, 
              double **vtx_list, int **vtx_connectivity, 
              void *ext_smooth_data)
{
    int ierr;
    int do_linear_program;
    int degenerate=0;
    int point_removed;
    int valid;
    double drand48(void);
    SMsmooth_data *smooth_data;
    SMparam *smooth_param;
    SMlocal_mesh *local_mesh;
    SMprocinfo *procinfo;
    SMuntangle_param *untangle_param;
    double drand48(void);
#ifdef OPTMS_LOCALTEST
    FILE *fp;
#endif

    SM_LOG_EVENT_BEGIN(__SM_TOTAL__);
    SM_LOG_EVENT_BEGIN(__SM_UNTANGLE__);

    OPTMS_CHECK_NULL(ext_smooth_data);

    smooth_data = (SMsmooth_data *) ext_smooth_data;
    local_mesh = smooth_data->local_mesh;
    smooth_param = smooth_data->smooth_param;
    procinfo = smooth_data->smooth_procinfo;
    untangle_param = smooth_data->untangle_param;

    /* initialize the local mesh triangles and the smoothing parameters */
    ierr = SMinitLocalMesh(num_incident_vtx,num_tri,free_vtx,vtx_list,
           vtx_connectivity,local_mesh,smooth_data->smooth_param); OPTMS_CHKERR(ierr);

    OPTMS_MATLAB_ON({
	if (local_mesh->dimension==2) {
	    /* plot out the initial local mesh */
	    OPTMS_DEBUG_PRINT(1,"Plotting the initial tangled local mesh \n");
	    if ((fp = fopen("local_tangle.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open local_tangle.m for writing\n");
            } 
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
    OPTMS_MATLAB_OFF});

    /* the choices for technique are OPTMS_LAPLACE_ONLY, OPTMS_LINEAR_PROGRAM_ONLY, 
       OPTMS_COMBINED_UNTANGLING*/
    do_linear_program = OPTMS_FALSE;

    switch(untangle_param->untangle_technique) {
    case OPTMS_LAPLACIAN_ONLY:
        SM_LOG_EVENT_BEGIN(__SM_LAP_SMOOTH__);
        ierr = SMinitLap(local_mesh->num_values,local_mesh->lap_info); OPTMS_CHKERR(ierr);
        ierr = SMcentroidSmoothMesh(local_mesh->num_incident_vtx, 
                             local_mesh->incident_vtx, local_mesh->free_vtx,
                             local_mesh->dimension); OPTMS_CHKERR(ierr);
        SM_LOG_EVENT_END(__SM_LAP_SMOOTH__);
        break;
    case OPTMS_LINEAR_PROGRAM_ONLY:
        do_linear_program = OPTMS_TRUE;
        break;
    case OPTMS_COMBINED_UNTANGLING:
        /* compute the original function values */
        SM_LOG_EVENT_BEGIN(__SM_LAP_SMOOTH__);
        ierr = SMinitLap(local_mesh->num_values,local_mesh->lap_info); OPTMS_CHKERR(ierr);
        ierr = SMcentroidSmoothMesh(local_mesh->num_incident_vtx,local_mesh->incident_vtx, 
                             local_mesh->free_vtx,local_mesh->dimension); OPTMS_CHKERR(ierr);
        SM_LOG_EVENT_END(__SM_LAP_SMOOTH__);
        ierr = SMvalidMesh(local_mesh, &valid); OPTMS_CHKERR(ierr);
        if (!valid) do_linear_program = OPTMS_TRUE;
        break;
    default:
        OPTMS_SETERR(OPTMS_INPUT_ERR,0,
            "You have entered an invalid option for the mesh untangling technique\n");
    }

    /* untangle the mesh, or if it's a degenerate problem do laplace smoothing */
    if (do_linear_program) {  
       ierr = SMuntangle_mesh(local_mesh, &degenerate); OPTMS_CHKERR(ierr);
       if (degenerate) {

            ierr = SMremoveIdenticalVtx(local_mesh->dimension,&num_incident_vtx,
                                        &num_tri,&vtx_list,&vtx_connectivity, &point_removed);
                   OPTMS_CHKERR(ierr);
            if (point_removed) {
	      /* try again */
              ierr = SMinitLocalMesh(num_incident_vtx,num_tri,free_vtx,vtx_list,
                            vtx_connectivity,local_mesh,
                            smooth_data->smooth_param); OPTMS_CHKERR(ierr);
              ierr = SMuntangle_mesh(local_mesh, &degenerate); OPTMS_CHKERR(ierr);
            }
        }
    }

    OPTMS_MATLAB_ON({
	if (local_mesh->dimension==2) {
	    /* plot out the initial local mesh */
	    OPTMS_DEBUG_PRINT(1,"Plotting the untangled local mesh \n");
	    if ((fp = fopen("local_tangle.m","w")) == NULL) {
               OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open local_tangle.m for writing\n");
            } 
	    ierr = SMwriteLocalMesh(fp,local_mesh); OPTMS_CHKERR(ierr);
	    fclose(fp);
	}
    OPTMS_MATLAB_OFF});

    OPTMS_COPY_VECTOR(free_vtx,local_mesh->free_vtx,local_mesh->dimension);

    SM_LOG_EVENT_END(__SM_UNTANGLE__);
    SM_LOG_EVENT_END(__SM_TOTAL__);
    SM_LOG_GLOBAL_TIME(__SM_TOTAL__);
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMuntangle_mesh"
int SMuntangle_mesh(SMlocal_mesh *local_mesh, int *degenerate)
{
   int ierr;
   int i, j, num_constraints, num_active, num_free;
   int feasible, solved;
   double **Amat, **Amat_T, **Amat_T_O;
   double *pi, *b, *feasible_x, x_sum;
   SMlp *lp_info;

   /* In this linear program we will maximize the minumum area of the local
       submesh using the simplex measure.  This measure is always convex, regardless 
       of the initial submesh, and we therefore expect to converge, unless the problem 
       is degenerate.  If the problem is degenerate, we exit and perform laplacian 
       smoothing */

   /* set the number of active points and the number of constraints */
   num_constraints = local_mesh->num_tri;
   num_active = local_mesh->dimension;
   num_free = num_constraints - num_active;   
   *degenerate=0;

   ierr = SMinitLP(local_mesh);  OPTMS_CHKERR(ierr);
   lp_info = local_mesh->lp_info;

    /* initialize the constraint matrix, feasible point, rhs, and solution pi */
   Amat          = lp_info->Amat;
   Amat_T        = lp_info->Amat_T;
   Amat_T_O      = lp_info->Amat_T_O;
   feasible_x    = lp_info->feasible_x;
   b             = lp_info->b;
   pi            = lp_info->pi;

   /* compute the matrix for the constaints; dimension num_constraints x num_active*/
   SMcomputeConstraintMatrix(local_mesh,num_constraints,Amat,b);
   /* create the tanspose matrix for the phase one solution to find a feasible point for the LP */
   ierr = SMtransposeMatrix2(Amat,num_constraints,num_active,Amat_T); OPTMS_CHKERR(ierr);
  
   ierr = SMdegenerate(num_active, num_constraints, Amat_T, b, degenerate); OPTMS_CHKERR(ierr);
   if (!(*degenerate)) {
        
        /* entering the phase 1 solution to find a feasible point for the linear program */
        /* send in the matrix of constraints and it will return a feasible X */
        SM_LOG_EVENT_BEGIN(__SM_PHASE1__);
        ierr = SMphaseOneLP(num_constraints,num_active,Amat_T,feasible_x,lp_info,&feasible); 
               OPTMS_CHKERR(ierr);
        SM_LOG_EVENT_END(__SM_PHASE1__);
        
        if (feasible) {

            SM_LOG_EVENT_BEGIN(__SM_LINEAR_PROG__);
           /* set up for solving the LP */
           /* make the last x part of the active set; normalize x so that the components sum to 1;
           reverse the sign of b, and add a row of ones to the matrix for the slack variables*/
           feasible_x[num_constraints-1]=1;
           x_sum = 0;
           for (i=0;i<num_constraints;i++)   x_sum += feasible_x[i];
           for (i=0;i<num_constraints;i++) {
               feasible_x[i] = feasible_x[i]/x_sum;
               b[i] = -b[i];
           }   
           for (i=0;i<num_active+1;i++) {
               if (i<num_active) {
                   for (j=0;j<num_constraints;j++)  Amat_T_O[i][j] = -Amat_T[i][j];
               }
               if (i==num_active) {
                   for (j=0;j<num_constraints;j++)  Amat_T_O[i][j] = 1;
               }
           }

           /* now we can solve the linear program */
           num_active=num_active+1;
           ierr = SMsolveLP(num_constraints, num_active, Amat_T_O, feasible_x,b,pi,lp_info,&solved);
                  OPTMS_CHKERR(ierr);
           SM_LOG_EVENT_END(__SM_LINEAR_PROG__);

           if (solved) {
               /* Report on the results */
               if (pi[local_mesh->dimension] >= 0) {
	          OPTMS_DEBUG_ACTION(2,
                      {printf("SUCCESS! min area %e ",pi[local_mesh->dimension]);}
                  );
	          OPTMS_DEBUG_ACTION(3,{
                      printf("x %f y %f ",pi[0],pi[1]);
                      if (local_mesh->dimension==2)  printf("\n");
                      else printf("z %f \n",pi[2]);
	           });
                  local_mesh->validity=OPTMS_VALID_MESH;
               } else {
	          OPTMS_DEBUG_ACTION(2,
                       {printf("STILL TANGLED; min area %f ",pi[local_mesh->dimension]);}
                  );
	          OPTMS_DEBUG_ACTION(3,{
                       printf("x %f y %f ",pi[0],pi[1]);
                       if (local_mesh->dimension==2) printf("\n");
                       else printf("z %f \n",pi[2]);
                  });
              }
              for (i=0;i<local_mesh->dimension;i++) local_mesh->free_vtx[i] = pi[i];
           } else { /* LP not solved */
              OPTMS_DEBUG_PRINT(2,"Didn't solve the LP problem successfully\n");
           }
        } else {  /* not feasible */
          OPTMS_DEBUG_PRINT(2,"Didn't get a feasible phase one solution \n");
       }
       /* regardless, return the corresponding x,y,z coordinate in space why? */
       return(ierr=0);
    } else { 
       OPTMS_DEBUG_PRINT(2,"Degenerate problem\n");
       return(ierr=0); 
    }  
}

#undef __FUNC__
#define __FUNC__ "SMsolveLP"
int SMsolveLP(int num_constraints, int num_active, double **Amat, double *feasible_x,
                           double *b, double *pi, SMlp *lp_info, int *success)
{
   int ierr;
   int i,j,iter,done;
   int one = 1;
   int num_free, free_count, active_count;
   int *active_ind, *free_ind;
   int info;
   int min_s_ind, min_alpha_ind, in_ind, out_ind;
   int *ipivot;
   double *Bmat, *Bmat_T;
   double  mins, min_alpha, *alpha;
   double *c, *s, *step;

   *success=1;

   /* initialize the work variables contained in lp info */
   num_free = num_constraints - num_active;
   ipivot          = lp_info->ipivot;
   free_ind      = lp_info->free_ind;
   active_ind   = lp_info->active_ind;
   c                 = lp_info->c;
   s                 = lp_info->s;
   alpha          = lp_info->alpha;
   step            = lp_info->step;
   Bmat          = lp_info->Bmat;
   Bmat_T      = lp_info->Bmat_T;

   num_active=5;
   for (i=0;i<num_active;i++) {
     active_ind[i] = 0;     ipivot[i] = 0;     alpha[i] = 0;
   }
   for (i=0;i<num_constraints;i++) {
     free_ind[i] = 0;     c[i] = 0;     s[i] = 0;     step[i] = 0;
   }
   for (i=0;i<num_active*num_active;i++) {
     Bmat[i] = 0;     Bmat_T[i] = 0;
   }
     
   /* determine which are the active and which are the free unknowns */
   free_count = 0; active_count=0;
   for (i=0;i<num_constraints;i++) {
     if ((fabs(feasible_x[i]))<1E-15) {
          free_ind[free_count++]=i; 
     } else {
          active_ind[active_count++]=i;
     }
   }
   num_active = active_count;   num_free = free_count;

   /* initialize c=(0 0 0 ... 0 1) */
  for (i=0;i<num_active-1;i++)  c[i] = 0;
  c[num_active-1]=1;
  
  /* check that we in fact have a feasible point, using the array step as temp work space*/
  for (i=0;i<num_active-1;i++) {
    step[i]=0;
    for (j=0;j<num_constraints;j++) {
        step[i]+=Amat[i][j]*feasible_x[j];
     }
     if (fabs(step[i]) > 10*OPTMS_MACHINE_EPS) {
       OPTMS_DEBUG_ACTION(2,{fprintf(stderr,"Not a feasible point %f\n",step[i]);}); 
       *success=0;
       return(ierr=0);
     }
   }
  for (j=0;j<num_constraints;j++) {
     if ((feasible_x[j] < 0) && (fabs(feasible_x[j]) > 10*OPTMS_MACHINE_EPS)) {
          OPTMS_DEBUG_ACTION(2,{fprintf(stderr,"Not a feasible point \n");}); 
         *success=0;
         return(ierr=0);
     }
  }

  /* move from feasible point to feasible point, decreasing the objective function at
     each step.  the number of feasible points is finite so we're gauranteed to stop */
   done =0; iter=-1;
   while (!done) {
      SM_LOG_EVENT_BEGIN(__SM_LP_ITER__);
      iter++;
      if (iter == 100) {
        done=1;
        OPTMS_DEBUG_PRINT(2,"Exceeded max number of iterations in LP solve\n"); 
        *success=0;
        return(ierr=0);
      }

    /* get the current constraint matrix */
      ierr = SMgetActiveMatrix(Amat,num_active,active_ind,Bmat); OPTMS_CHKERR(ierr);
      ierr = SMgetActiveRHS(b,num_active,active_ind,pi); OPTMS_CHKERR(ierr);

      /* solve Bmat_T * pi = c */
      ierr = SMtransposeMatrix(Bmat,num_active,num_active,Bmat_T); OPTMS_CHKERR(ierr);

      /* solve the system using lapack */
      DGESV(&num_active, &one, Bmat_T, &num_active, ipivot, pi, &num_active, &info);
      if (info!=0) {
         *success=0;
         return(ierr=0);
      }

      /* compute the slack variables */
      mins = 1E300;
      for (i=0;i<num_free;i++) {
         s[i]=b[free_ind[i]];
         for (j=0;j<num_active;j++) {
           s[i] -= Amat[j][free_ind[i]]*pi[j];
         }
         if (s[i]<mins) {
            mins=s[i];  min_s_ind = i;
         }  
      }
 
      /* if all of the slack variables are greater than 0, done 
          the complemintarity condition is satisfied */
      if (mins>-10*OPTMS_MACHINE_EPS) {
         done = 1;
         OPTMS_DEBUG_ACTION(2,
              {fprintf(stdout,"Equilibrium point found in LP in %d iterations\n",iter);}
          );
      } else {
          /*  swap a point in and out of the basis */
          /* we'll be bringing in a new index to the active set.. in particular, the index
               associated with the minimum slack variable */
          in_ind = free_ind[min_s_ind];

          /* find the steps that we can take in each direction before we get to
          an infeasible region; the right hand side is the in_ind column of the 
          A matrix */
          for (i=0;i<num_active;i++)  step[i] = Amat[i][in_ind];
    
          DGESV(&num_active, &one, Bmat,  &num_active, ipivot, step, &num_active, &info);
          if (info!=0) {
            *success=0;
            return(ierr=0);
          }

         /* find the minimum step size */
         min_alpha = 1E300;
         for (i=0;i<num_active;i++){
             alpha[i] = feasible_x[active_ind[i]] / step[i];
             if ((alpha[i] < min_alpha) && alpha[i]>0) {
                 min_alpha = alpha[i];
                 out_ind = active_ind[i];
                 min_alpha_ind = i;
             }
         }

        /* update x, active_ind, and free_ind for the new basis */
        for (i=0;i<num_active;i++) {
          feasible_x[active_ind[i]] -= min_alpha*step[i];
        }
        feasible_x[in_ind] = min_alpha;
        free_ind[min_s_ind] = out_ind;
        active_ind[min_alpha_ind] = in_ind;
         
      } /* end else */
      SM_LOG_EVENT_END(__SM_LP_ITER__);
   } /* end while */

   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMphaseOneLP"
int SMphaseOneLP(int num_constraints, int num_active, double **Amat, double *feasible_x,
                  SMlp *lp_info, int *feasible)
{
   int ierr;
   int i, j, iter;
   int info, one = 1;
   int *ipivot;
   int done, free_count;
   int active;
   int num_free, *active_ind, *free_ind;
   int min_s_ind, min_alpha_ind,in_ind, out_ind;
   double  mins, min_alpha, *alpha;
   double temp;
   double b[OPTMS_MAX_NUM_TRI], *c, *pi, *s, *step;
   double *Bmat, *Bmat_T;
   double **AAmat;
   int degenerate;
   double vec1[3], vec2[3];
   double result[3];

   num_free = num_constraints-num_active;
   *feasible = 1;

   /* initialize the work variables contained in lp_info */
   ipivot       = lp_info->ipivot;
   free_ind     = lp_info->free_ind;
   active_ind   = lp_info->active_ind;
   c            = lp_info->c;
   pi           = lp_info->pi;
   s            = lp_info->s;
   alpha        = lp_info->alpha;
   step         = lp_info->step;
   Bmat         = lp_info->Bmat;
   Bmat_T       = lp_info->Bmat_T;
   AAmat        = lp_info->AAmat;

   for (i=0;i<num_active;i++) {
     active_ind[i] = 0;     ipivot[i] = 0;     alpha[i] = 0;     pi[i]=0;
     for (j=0;j<num_constraints;j++) AAmat[i][j]=0;
   }
   for (i=0;i<num_active*num_active;i++) {
     Bmat[i] = 0;     Bmat_T[i] = 0;
   }
   for (i=0;i<num_constraints;i++) {
     free_ind[i] = 0;     c[i] = 0;     s[i] = 0;     step[i] = 0;
   }

   /* initialize AAmat */
   for (i=0;i<num_active;i++) {
      for (j=0;j<num_constraints;j++) AAmat[i][j] = Amat[i][j];
   }

   /* initialize an active set */
   for (i=0;i<num_active-1;i++) active_ind[i]=i;   
   active_ind[num_active-1]=num_constraints-1;

   if (num_active == 3) {
      /* we need to be careful of singular systems */
      /* we know the last one isn't dependent on one or two other columns, other
          wise it would have been degenerate... make sure the first two aren't dependent
          on each other */
       degenerate = 0;
       vec1[0]=AAmat[0][active_ind[0]];    vec2[0]=AAmat[0][active_ind[1]];
       vec1[1]=AAmat[1][active_ind[0]];    vec2[1]=AAmat[1][active_ind[1]];
       vec1[2]=AAmat[2][active_ind[0]];    vec2[2]=AAmat[2][active_ind[1]];
       vCross(vec1,vec2,result); if (dMagnitude(result) < 1E-13) degenerate = 1;
       while (degenerate && active_ind[1]<num_constraints-1) {
          degenerate=0;
          active_ind[1]++; 
          vec1[0]=AAmat[0][active_ind[0]];    vec2[0]=AAmat[0][active_ind[1]];
          vec1[1]=AAmat[1][active_ind[0]];    vec2[1]=AAmat[1][active_ind[1]];
          vec1[2]=AAmat[2][active_ind[0]];    vec2[2]=AAmat[2][active_ind[1]];
          vCross(vec1,vec2,result); if (dMagnitude(result) < 1E-13) degenerate = 1;
       }
   }

   /* initialize the free indices */
   free_count = 0;
   for (i=0;i<num_constraints;i++) {
     active=0;
     for (j=0;j<num_active;j++) {
         if (i == active_ind[j]) active=1;
     }
     if (!active)  free_ind[free_count++] = i;
   }

   /* the initial guess for a feasible x is 1 at the active points, 0 elsewhere */
   /* c = (0 0 0... 1) */
   for (i=0;i<num_constraints;i++) {
      feasible_x[i] = 0;      c[i] = 0;
   }

   for (i=0;i<num_active;i++) {
        b[i] = -Amat[i][num_constraints-1];
        feasible_x[active_ind[i]]=1;
   }
   c[num_constraints-1] = 1;

   /* set up the matrix with the slack variable incorporated */
   for (i=0;i<num_active;i++) {
     AAmat[i][num_constraints-1] = b[i];
     for (j=0;j<num_constraints-1;j++) {
        AAmat[i][num_constraints-1] += -AAmat[i][j]*feasible_x[j];
       }
     }

   /* test that the first point is feasible */
   for (i=0;i<num_active;i++) {
      step[i]=0;
      for (j=0;j<num_constraints;j++) {
        step[i] += AAmat[i][j]*feasible_x[j];
      }
      if (fabs(step[i] +Amat[i][num_constraints-1])>10*OPTMS_MACHINE_EPS) {
            OPTMS_DEBUG_PRINT(2,"error in first point of phase 1\n");
            *feasible=0;
            return(ierr=0);
      }
   }

   /* Get ready to enter the loop to find a feasible point */
   done = 0;  iter = 0; 
   while (!done) {
      iter++;
      if (iter == 100) {
        done=1;
        OPTMS_DEBUG_PRINT(2,"Exceeded max number of iterations in phase 1\n");
        *feasible=0;
        return(ierr=0);
      }
      
      /* get the current constraint matrix */
      ierr = SMgetActiveMatrix(AAmat,num_active,active_ind,Bmat); OPTMS_CHKERR(ierr);
      ierr = SMgetActiveRHS(c,num_active,active_ind,pi); OPTMS_CHKERR(ierr);

      /* solve Bmat_T * pi = c */
      ierr = SMtransposeMatrix(Bmat,num_active,num_active,Bmat_T); OPTMS_CHKERR(ierr);

      /* note; on input pi contains rhs, on output pi contains the solution */
      DGESV(&num_active, &one, Bmat_T, &num_active, ipivot, pi, &num_active, &info);
      if (info != 0) {
         OPTMS_DEBUG_PRINT(2,"Problem in linear solve, phase 1\n"); 
         *feasible=0;
         return(ierr=0);
      }

      /* the slack variables are the dot products of the free columns of aa_T with pi */
      /* these are the distance from the current point to the constraint */
      mins = 1E300;
      for (i=0;i<num_free;i++) {
         s[i]=c[free_ind[i]];
         for (j=0;j<num_active;j++)  s[i] -= AAmat[j][free_ind[i]]*pi[j];
         if (s[i]<mins) {
            mins=s[i];  min_s_ind = i;
         }  
      }

       /* if all of the slack variables are greater than 0, done */
      if (mins>-10*OPTMS_MACHINE_EPS) {
         done = 1;
         OPTMS_DEBUG_PRINT(2,"Not Equilibrium in phase 1, Amat_T*(-pi) > 0\n");
         *feasible=1;
         return(ierr=0);
      } else {
      /* we'll be bringing in a new index to the active set.. in particular, the index
          associated with the minimum slack variable */
          in_ind = free_ind[min_s_ind];

          /* find the steps that we can take in each direction before we get to
          an infeasible region; the right hand side is the in_ind column of the 
          A matrix */
          for (i=0;i<num_active;i++)  step[i] = Amat[i][in_ind]; 

         /* solve the system to get the step sizes before each constraint goes infeasible */
          DGESV(&num_active, &one, Bmat, &num_active, ipivot, step, &num_active, &info);
          if (info != 0) {
              OPTMS_DEBUG_PRINT(2,"Problem in linear solve, phase 1\n"); 
              *feasible=0;
              return(ierr=0);
          }

         /* find the minimum step size */
         min_alpha = 1E300;
         for (i=0;i<num_active;i++){
             alpha[i] = feasible_x[active_ind[i]] / step[i];
             if ((alpha[i] < min_alpha) && alpha[i]>0) {
                 min_alpha = alpha[i];   min_alpha_ind = i;
                 out_ind = active_ind[i];
             }
         }

         /* stop if we've found a feasible point*/
          if (min_alpha*step[num_active-1] >= feasible_x[num_constraints-1]-10*OPTMS_MACHINE_EPS) {
               done=0;
               for (i=0;i<num_active;i++) {
               temp=0;
               for (j=0;j<num_constraints;j++) {
                  temp+=AAmat[i][j]*feasible_x[j];
               }
              if (fabs(temp+Amat[i][num_constraints-1]) < 10*OPTMS_MACHINE_EPS) done++;
             }
             if (done==num_active) {
                 done = 1;
                 min_alpha = feasible_x[num_constraints-1]/step[num_active-1];
                 out_ind = num_active-1;
              } 
        }

        /* update x, active_ind, and free_ind for the new basis */
        for (i=0;i<num_active;i++) {
          feasible_x[active_ind[i]] -= min_alpha*step[i];
        }
        feasible_x[in_ind] = min_alpha;
        free_ind[min_s_ind] = out_ind;
        active_ind[min_alpha_ind] = in_ind;

      } /* end else */
   } /* end while */

   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMdegnerate"
int  SMdegenerate(int num_active, int num_constraints, double **A, double *b, int *degenerate)
{
   int ierr;
   int i;
   int dependent_col;
   int column1, column2, column3;
   double **Atest;
   double temp;
   int icol1, icol2;
   int singular;
   double vec1[3], vec2[3], vec3[3];
   double result[3];

   /* determine which column to make the dependent one.. want to be sure that it
       isn't a multiple of any other column which leads to a degenerate case... start with
       num_constraints -1 and work your way up */

   OPTMS_MALLOC(Atest,(double **),sizeof(double *)*num_active,1);
   for (i=0;i<num_active;i++) OPTMS_MALLOC(Atest[i],(double *),sizeof(double)*num_active,1);
  
   if (num_active==2) {
       dependent_col = num_constraints -1;
       *degenerate = OPTMS_FALSE;
       for (i=0;i<num_constraints;i++) {
         if (i != dependent_col) {
            Atest[0][0] = A[0][i];  Atest[0][1] = A[0][dependent_col];
            Atest[1][0] = A[1][i];  Atest[1][1] = A[1][dependent_col];
            ierr = SMsingularTest(num_active,Atest,&singular); OPTMS_CHKERR(ierr);
            if (singular)  *degenerate = OPTMS_TRUE;
         }
       }
       while (*degenerate && dependent_col>0) {
          dependent_col = dependent_col-1;
          *degenerate = OPTMS_FALSE;
          for (i=0;i<num_constraints;i++) {
             if (i != dependent_col) {
               Atest[0][0] = A[0][i];  Atest[0][1] = A[0][dependent_col];
               Atest[1][0] = A[1][i];  Atest[1][1] = A[1][dependent_col];
               ierr = SMsingularTest(num_active,Atest,&singular); OPTMS_CHKERR(ierr);
               if (singular)  *degenerate = OPTMS_TRUE;
             }
         }
      }
    } else if (num_active == 3) {
       dependent_col = num_constraints -1;
       *degenerate = OPTMS_FALSE;
       for (icol1=0;icol1<num_constraints;icol1++) {
         for (icol2=0;icol2<num_constraints;icol2++) {
           if (!(*degenerate)) {
             if ((icol1 != dependent_col) && (icol2 != dependent_col) && (icol1 != icol2)) {
                column1 = icol1; column2 = icol2; column3 = dependent_col;
                vec1[0] = A[0][icol1];         vec2[0] = A[0][icol2];       vec3[0] = A[0][dependent_col];
                vec1[1] = A[1][icol1];         vec2[1] = A[1][icol2];       vec3[1] = A[1][dependent_col];
                vec1[2] = A[2][icol1];         vec2[2] = A[2][icol2];       vec3[2] = A[2][dependent_col];
                vCross(vec1,vec3,result);  if (dMagnitude(result) < 1E-13) *degenerate = OPTMS_TRUE;
                vCross(vec2,vec3,result);  if (dMagnitude(result) < 1E-13) *degenerate = OPTMS_TRUE;
                vCross(vec1,vec2,result);  
                if (dMagnitude(result) > 1E-13) {
                   Atest[0][0] = A[0][icol1];  Atest[0][1] = A[0][icol2]; Atest[0][2] = A[0][dependent_col];
                   Atest[1][0] = A[1][icol1];  Atest[1][1] = A[1][icol2]; Atest[1][2] = A[1][dependent_col];
                   Atest[2][0] = A[2][icol1];  Atest[2][1] = A[2][icol2]; Atest[2][2] = A[2][dependent_col];
                   ierr = SMsingularTest(num_active,Atest,&singular); OPTMS_CHKERR(ierr);
                   if (singular)  *degenerate = OPTMS_TRUE;
                }
              }
            }
          }
        }

       while (*degenerate && dependent_col>0) {
          dependent_col = dependent_col-1;
          *degenerate = OPTMS_FALSE;
          for (icol1=0;icol1<num_constraints;icol1++) {
             for (icol2=0;icol2<num_constraints;icol2++) {
                if (!(*degenerate)) {
                   if ((icol1 != dependent_col) && (icol2 != dependent_col) && (icol1 != icol2)) {
                      vec1[0] = A[0][icol1];     vec2[0] = A[0][icol2];     vec3[0] = A[0][dependent_col];
                      vec1[1] = A[1][icol1];     vec2[1] = A[1][icol2];     vec3[1] = A[1][dependent_col];
                      vec1[2] = A[2][icol1];     vec2[2] = A[2][icol2];     vec3[2] = A[2][dependent_col];
                      vCross(vec1,vec3,result);  
                      if (dMagnitude(result) < 1E-13) *degenerate = OPTMS_TRUE;
                      vCross(vec2,vec3,result);  
                      if (dMagnitude(result) < 1E-13) *degenerate = OPTMS_TRUE;
                      vCross(vec1,vec2,result);  
                      if (dMagnitude(result) > 1E-13) {
                        Atest[0][0]=A[0][icol1];  Atest[0][1]=A[0][icol2]; Atest[0][2]=A[0][dependent_col];
                        Atest[1][0]=A[1][icol1];  Atest[1][1]=A[1][icol2]; Atest[1][2]=A[1][dependent_col];
                        Atest[2][0]=A[2][icol1];  Atest[2][1]=A[2][icol2]; Atest[2][2]=A[2][dependent_col];
                        ierr = SMsingularTest(num_active,Atest,&singular); OPTMS_CHKERR(ierr);
                        if (singular)  *degenerate = OPTMS_TRUE;
                     }
                  }
                }
             }
          }
       }
    }

   /* free the test matrix */
   for (i=0;i<num_active;i++) OPTMS_FREE(Atest[i]);
   OPTMS_FREE(Atest);

   if (*degenerate) {
     OPTMS_DEBUG_ACTION(2,{ 
          printf("A degenerate problem, no idea what to do\n");
          for (i=0;i<num_constraints;i++) {
             printf("A(1,%d)=%e;\n",i+1,A[0][i]);
             printf("A(2,%d)=%e;\n",i+1,A[1][i]);
             printf("A(3,%d)=%e;\n",i+1,A[2][i]);
          }
    });
       return (ierr=0);
   }

   /* create the matrix with the non degenerate column last */
   if (dependent_col != num_constraints-1) {
      temp = b[dependent_col];
      b[dependent_col]=b[num_constraints-1];
      b[num_constraints-1] = temp;
      for (i=0;i<num_active;i++) {
             temp = A[i][dependent_col];
             A[i][dependent_col] = A[i][num_constraints-1];
             A[i][num_constraints-1] = temp;
        }
   }

   return (ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMgetActiveMatrix"
int SMgetActiveMatrix(double **Amat, int num_active,int *active_ind, double *Bmat)
{
    int ierr;
    int i,j,k;
    k=0;
    for (j=0;j<num_active;j++) {
      for (i=0;i<num_active;i++) {
          Bmat[k++] = Amat[i][active_ind[j]];
       }
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMgetActiveRHS"
int SMgetActiveRHS(double *c,int num_active,int *active_ind, double *c_active)
{
    int ierr;
    int i;
    for (i=0;i<num_active;i++) {
          c_active[i] = c[active_ind[i]];
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMcomputeConstraintMatrix"
int SMcomputeConstraintMatrix(SMlocal_mesh *local_mesh, int num_constriants,
                              double **Amat, double *b)
{
   int ierr;
   int i, ind1, ind2;
   double x1, x2, y1, y2;
   int **vtx_connectivity;
   double **incident_vtx;
   int ind3;
   double x3, y3;
   double z1, z2, z3;
   double a1, b1, c1, d1;

   vtx_connectivity = local_mesh->vtx_connectivity;
   incident_vtx = local_mesh->incident_vtx;

   if (local_mesh->dimension == 2 ) {
    /* compute the matrix and right hand side such that A*X-B=Tri Area for any X */
       for (i=0;i<local_mesh->num_tri;i++) {
           ind1 = vtx_connectivity[i][0];      
           ind2 = vtx_connectivity[i][1];
           x1 = incident_vtx[ind1][0];       y1 = incident_vtx[ind1][1];
           x2 = incident_vtx[ind2][0];       y2 = incident_vtx[ind2][1];

           Amat[i][0] = .5*(y1-y2);
           Amat[i][1] = .5*(x2-x1);
           b[i] = .5*(x2*y1-x1*y2);
       }
   } else if (local_mesh->dimension == 3) {
       for (i=0;i<local_mesh->num_tri;i++) {
           ind1 = vtx_connectivity[i][0];      
           ind2 = vtx_connectivity[i][1];
           ind3 = vtx_connectivity[i][2];
           x1 = incident_vtx[ind1][0];       y1 = incident_vtx[ind1][1];    z1 = incident_vtx[ind1][2];
           x2 = incident_vtx[ind2][0];       y2 = incident_vtx[ind2][1];    z2 = incident_vtx[ind2][2];
           x3 = incident_vtx[ind3][0];       y3 = incident_vtx[ind3][1];    z3 = incident_vtx[ind3][2];

           a1 = x1*(y2*z3-z2*y3) - x2*(y1*z3-y3*z1) + x3*(y1*z2 - z1*y2);
           b1 = (y2*z3-z2*y3) - (y1*z3-y3*z1) + (y1*z2-y2*z1);
           c1 = x1*(z3-z2) - x2*(z3-z1)+x3*(z2-z1);
           d1 = x1*(y2-y3) - x2*(y1-y3)+x3*(y1-y2);

           Amat[i][0] = -b1/6;
           Amat[i][1] = -c1/6;
           Amat[i][2] = -d1/6;
           b[i] = -a1/6;
       }
   }
   return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMremoveIdenticalVtx"
int SMremoveIdenticalVtx(int dimension, int *num_incident_vtx,int *num_tri, 
                         double ***vtx_list, int ***vtx_connectivity, int *point_removed)
{
    int ierr;
    int i,j,k;
    int id1, id2;
    int *ind;
    int count;
    int number_identical;
    int matches_id1, matches_id2;
    int num_original_vtx, num_original_tet;
    double *coord1;
    double *coord2;

    *point_removed=OPTMS_FALSE;
    
    OPTMS_MALLOC(coord1,(double *),sizeof(double)*dimension,1);
    OPTMS_MALLOC(coord2,(double *),sizeof(double)*dimension,1);
    OPTMS_MALLOC(ind,(int *),sizeof(int)*dimension,1);

    id1=-1;
    id2=-1;

    /* identify the identical vertices */
    for (i=0;i<(*num_incident_vtx);i++) {
      for (k=0;k<dimension;k++) coord1[k] = (*vtx_list)[i][k]; 
      for (j=0;j<i;j++) {
         number_identical = 0;
         for (k=0;k<dimension;k++) {
           if (fabs(coord1[k] - (*vtx_list)[j][k])<1E-8) number_identical++;
         }
         if (number_identical == dimension) {
	   /*            printf("Identical vertices %d %d : %f %f",i,j,coord1[0],coord1[1]);
            if (dimension==2) printf("\n");
            if (dimension==3) printf(" %f\n",coord1[2]); */
            id1 = j; id2 = i;
         }
       }
    }
      
    if ((id1 != -1) && (id2 != -1)) {
    /* remove one of the vertices from the vertex list */
      count = 0;
      num_original_vtx = (*num_incident_vtx);
      for (i=0;i<num_original_vtx;i++) {
        if (i!=id2) {
          for (k=0;k<dimension;k++) (*vtx_list)[count][k] = (*vtx_list)[i][k];
          count++;
        }
      }
      (*num_incident_vtx)--;
      *point_removed=OPTMS_TRUE;
    } else {
      *point_removed=OPTMS_FALSE;
      /*      printf("The problem is degenerate but I couldn't find identical vertices \n"); */
    }
     
    /* clean up the connectivity list
         - remove any tet that contains both id1 and id2
         - replace any occurance of id2 with id1
    */
    if ((id1 != -1) && (id2 != -1)) {
      count = 0;
      num_original_tet = (*num_tri);
      for (i=0;i<num_original_tet;i++) {
        matches_id1 = OPTMS_FALSE;
        matches_id2 = OPTMS_FALSE;
        for (k=0;k<dimension;k++) {
          if ((*vtx_connectivity)[i][k]==id1) matches_id1=OPTMS_TRUE;
          if ((*vtx_connectivity)[i][k]==id2) matches_id2=OPTMS_TRUE;
        }
        if (matches_id1 && matches_id2) {
	    /* remove the tetrahedra */
            (*num_tri)-- ;          
        } else if (matches_id2 && !matches_id1) {
            /* replace id2 with id2 */
            for (k=0;k<dimension;k++) {
              if ((*vtx_connectivity)[i][k]==id2) {
                (*vtx_connectivity)[count][k] = id1;
              } else {
                (*vtx_connectivity)[count][k] = (*vtx_connectivity)[i][k];
                if ((*vtx_connectivity)[count][k] > id2) (*vtx_connectivity)[count][k]--;
              }
            }
            count++;
        } else {
	    /* keep the tet */
            for (k=0;k<dimension;k++) {
                (*vtx_connectivity)[count][k] = (*vtx_connectivity)[i][k];
                if ((*vtx_connectivity)[count][k] > id2) (*vtx_connectivity)[count][k]--;
             }
            count++;
        }
      }
    }

    if (dimension == 3) {
      if (2*(*num_incident_vtx)-4 != (*num_tri)) {
        printf("Problems eliminating identical vertex; num_incident_vtx %d and num_tri %d\n",
              (*num_incident_vtx),(*num_tri));
      }
    } else if (dimension == 2) {
      if ((*num_incident_vtx) != (*num_tri)) {
        printf("Problems eliminating identical vertex; num_incident_vtx %d and num_tri %d\n",
              (*num_incident_vtx),(*num_tri));
      }
    }
  OPTMS_FREE(coord1);
  OPTMS_FREE(coord2);
  OPTMS_FREE(ind);
  return(ierr=0);
}

