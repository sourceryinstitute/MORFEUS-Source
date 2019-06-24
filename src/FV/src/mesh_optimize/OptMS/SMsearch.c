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
#define __FUNC__ "SMsearchDirection" 
int SMsearchDirection(SMlocal_mesh *local_mesh)
{
    int        ierr;
    int        i, num_active;
    int        *active_ind;
    int        viable;
    int        dimension;
    int        singular;
    double     a, b, c, denom;
    double     **gradient, **dir;
    double     R0, R1;
    double     **G, **P, *x;
    double     search_mag;
    SMoptimal  *opt_info;

    OPTMS_CHECK_NULL(local_mesh);

    SM_LOG_EVENT_BEGIN(__SM_SEARCH__); 
    opt_info = local_mesh->opt_info;
    gradient = opt_info->gradient;
    num_active = opt_info->active->num_active;
    active_ind = opt_info->active->active_ind;
 
    if (num_active==0) 
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"No active values in search");

    dimension = local_mesh->dimension;

    switch(num_active) {
    case 1: 
        OPTMS_COPY_VECTOR(opt_info->search,gradient[active_ind[0]],dimension);
        opt_info->steepest = active_ind[0];
        break;
    case 2:
        /* if there are two active points, move in the direction of the
	   intersection of the planes.  This is the steepest descent
           direction found by analytically solving the QP */
        
        /* set up the active gradient directions */
        ierr = SMgetActiveDirections(num_active,gradient,active_ind,dimension,&dir);
               OPTMS_CHKERR(ierr);

        /* form the grammian */
        ierr = SMformGrammian(num_active,dir,opt_info->G,dimension); OPTMS_CHKERR(ierr);
        ierr = SMformPDGrammian(opt_info); OPTMS_CHKERR(ierr);
        G = opt_info->G;

        denom = (G[0][0] + G[1][1] - 2*G[0][1]);
        viable = 1;
        if (fabs(denom) > OPTMS_MACHINE_EPS) {
	  /* gradients are LI, move along their intersection */
           b = (G[0][0] - G[0][1])/denom;  
           a = 1 - b;
           if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
           if (viable) {
             for (i=0;i<dimension;i++) {
               opt_info->search[i] = a*dir[0][i] + b*dir[1][i];
             }
           } else {
             /* the gradients are dependent, move along one face */
             OPTMS_COPY_VECTOR(opt_info->search,dir[0],dimension);
           }
        } else {
	   /* the gradients are dependent, move along one face */
           OPTMS_COPY_VECTOR(opt_info->search,dir[0],dimension);
        }
        opt_info->steepest = active_ind[0];

        for (i=0;i<num_active;i++) OPTMS_FREE(dir[i]);
	OPTMS_FREE(dir);

        break;
    default:
        /* as in case 2: solve the QP problem to find the steepest
           descent direction.  This can be done analytically - as
           is done in Gill, Murray and Wright 
             for 3 active points in 3 directions - test PD of G
             otherwise we know it's SP SD so search edges and faces */

        /* get the active gradient directions */
        ierr = SMgetActiveDirections(num_active,gradient,active_ind,dimension,&dir);
               OPTMS_CHKERR(ierr);

        /* form the entries of the grammian matrix */
        ierr = SMformGrammian(num_active,dir,opt_info->G,dimension); OPTMS_CHKERR(ierr);
        ierr = SMformPDGrammian(opt_info); OPTMS_CHKERR(ierr);
        G = opt_info->G;

        switch(dimension) {
        case 2:
            ierr = SMsearchEdgesFaces(num_active, G, dir, opt_info); OPTMS_CHKERR(ierr);
            break;
        case 3:
	  if (num_active == 3) {
              ierr = SMsingularTest(num_active,opt_info->G,&singular); OPTMS_CHKERR(ierr);
              if (!singular) {
	        /* form the entries of P=Z^T G Z where Z = [-1...-1; I ] */
                ierr = SMformReducedMatrix(num_active,G,&P); OPTMS_CHKERR(ierr);
                /* form  the RHS and solve the system for the coeffs */
                R0 = G[0][0] - G[1][0];  R1 = G[0][0] - G[2][0];
                ierr = SMsolve2x2(P[0][0],P[0][1],P[1][0],P[1][1],R0,R1,&x);
                       OPTMS_CHKERR(ierr);
                if (x!=NULL) {
                	a = 1 - x[0] - x[1];  b = x[0];  c = x[1];
                	for (i=0;i<dimension;i++) {
                    	opt_info->search[i] = a*dir[0][i] + b*dir[1][i] + 
                       	                      c*dir[2][i];
                	}
                	opt_info->steepest = active_ind[0];
                	for (i=0;i<num_active-1;i++)  OPTMS_FREE(P[i]);  
                	OPTMS_FREE(P);  OPTMS_FREE(x);
                } else { 
                  	ierr = SMsearchEdgesFaces(num_active, G, dir, opt_info);
                        OPTMS_CHKERR(ierr);
                }
	      } else {
                 ierr = SMsearchEdgesFaces(num_active, G, dir, opt_info);
                         OPTMS_CHKERR(ierr);
	      }
            } else {
                ierr = SMsearchEdgesFaces(num_active, G, dir, opt_info);
                       OPTMS_CHKERR(ierr);
            }
            break;
        }
        for (i=0;i<num_active;i++) OPTMS_FREE(dir[i]);
	OPTMS_FREE(dir);
    }
    /* if the search direction is essentially zero, equilibrium pt */
    OPTMS_DOT(search_mag,opt_info->search,opt_info->search,dimension);
    if (fabs(search_mag)<1E-13) opt_info->status = OPTMS_ZERO_SEARCH;
    else OPTMS_NORMALIZE(opt_info->search,dimension);

    SM_LOG_EVENT_END(__SM_SEARCH__); 
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMsearchEdgesFaces" 
int SMsearchEdgesFaces(int num_active, double **G, double **dir, 
                        SMoptimal *opt_info) 
{
    int ierr;
    int i,j,k;
    int viable;
    int dimension;
    double a,b,denom;
    double g_bar[OPTMS_MAX_DIM];
    double search[OPTMS_MAX_DIM];
    double projection, min_projection;

    dimension = opt_info->dimension;

    OPTMS_CHECK_NULL(opt_info);
    if ( (dimension != 2) && (dimension != 3)) {
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Dimension must be 2 or 3");
    }

    SM_LOG_EVENT_BEGIN(__SM_EDGE_FACE__); 
    /* initialize the search direction to 0,0 */
    for (i=0;i<dimension;i++) search[i] = 0;

    /* Check for viable faces */
    min_projection = OPTMS_BIG_POS_NMBR;
    for (i=0; i<num_active; i++) {
        /* FACE I */
        viable = 1;

        /* test the viability */
        for (j=0;j<num_active;j++) {       /* lagrange multipliers>0 */
             if (G[j][i] < 0) viable = 0;
        }
       
        /* find the minimum of viable directions */
        if ((viable) && (G[i][i] < min_projection)) {
            min_projection = G[i][i];
            OPTMS_COPY_VECTOR(search,dir[i],dimension);
            opt_info->steepest = opt_info->active->active_ind[i];
        }
    
       /* INTERSECTION IJ */
       for (j=i+1; j<num_active; j++) {
          viable = 1;

          /* find the coefficients of the intersection 
             and test the viability */
          denom = 2*G[i][j] - G[i][i] - G[j][j];
          a = b = 0;
          if (fabs(denom) > OPTMS_MACHINE_EPS) {
             b = (G[i][j] - G[i][i])/denom;
             a = 1 - b;
             if ((b < 0) || (b > 1)) viable=0;  /* 0 < b < 1 */
	     for (k=0;k<num_active;k++) {       /* lagrange multipliers>0 */
                 if ((a*G[k][i] + b*G[k][j]) <= 0) viable=0;
             }
          } else {
             viable = 0;                        /* Linearly dependent */
          }

          /* find the minimum of viable directions */
          if (viable) {
             for (k=0;k<dimension;k++) {
                g_bar[k] = a * dir[i][k] + b * dir[j][k];
             }
             OPTMS_DOT(projection,g_bar,g_bar,dimension);
             if (projection < min_projection) {
	        min_projection = projection;
                OPTMS_COPY_VECTOR(search,g_bar,dimension);
                opt_info->steepest = opt_info->active->active_ind[i];
             }
          }
       }
    }
    if (opt_info->status != OPTMS_EQUILIBRIUM) {
        OPTMS_COPY_VECTOR(opt_info->search,search,dimension);
    }
    SM_LOG_EVENT_END(__SM_EDGE_FACE__); 
    return(ierr=0);
}         

#undef __FUNC__
#define __FUNC__ "SMgetActiveDirections" 
int SMgetActiveDirections(int num_active, double **gradient,          
                          int *active_ind, int dimension, double ***dir)
{
    int ierr;
    int i;

    OPTMS_MALLOC((*dir),(double **),sizeof(double *)*num_active,1);
    for (i=0;i<num_active;i++) {
        OPTMS_MALLOC((*dir)[i],(double *),sizeof(double)*dimension,1);
        OPTMS_COPY_VECTOR((*dir)[i],gradient[active_ind[i]],dimension);
    }
    return(ierr=0);
}

