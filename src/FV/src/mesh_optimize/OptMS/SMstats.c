/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <stdlib.h>
#if defined(WIN32) || defined(_WIN32)
#define _USE_MATH_DEFINES
#include "wininclude/unistd.h"
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

extern int fclose (FILE *);

#undef __FUNC__
#define __FUNC__ "SMaccumulateStats"
int SMaccumulateStats(SMlocal_mesh *local_mesh, SMparam *smooth_param, 
		   SMstats *smooth_stats)
{
    int ierr;
    double imp;

    /* this cell was smoothed */
    smooth_stats->total_cells_smoothed += 1;

    /* Laplacian statistics */
    if (local_mesh->lap_done) {
       smooth_stats->num_cells_laplaced += 1;
       if (local_mesh->lap_info->lap_accepted == OPTMS_FALSE) {
           if (local_mesh->lap_info->lap_invalid == OPTMS_TRUE) {
              smooth_stats->num_lap_invalid += 1;
           } else {
              smooth_stats->num_lap_worse += 1;
           }
       }
    }

    /* if combined - how many was laplacian enough for */
    if ((smooth_param->smooth_technique != OPTMS_OPTIMIZATION_ONLY) &&
        (smooth_param->smooth_technique != OPTMS_SMART_LAPLACIAN_ONLY) &&
        (local_mesh->lap_done) &&
        (local_mesh->lap_info->lap_accepted) &&
        (local_mesh->lap_info->laplacian_value > 
         smooth_param->lap_accept_value) ) {
        smooth_stats->num_lap_enough += 1;
    }

    /* if optimization used */
    if (local_mesh->opt_done) {
        smooth_stats->num_cells_opted += 1;
        smooth_stats->opt_iter_count += local_mesh->opt_info->opt_iter_count;

        /* optimization termination status */
        switch(local_mesh->opt_info->status) {
        case OPTMS_EQUILIBRIUM:
           smooth_stats->num_equil += 1;
           if (local_mesh->opt_info->iter_count == 0)
  	      smooth_stats->num_started_equil += 1;
           break;
        case OPTMS_ZERO_SEARCH:
           smooth_stats->num_zero_search += 1;
           break;
        case OPTMS_IMP_TOO_SMALL:
           smooth_stats->num_imp_too_small += 1;
           break;
        case OPTMS_FLAT_NO_IMP:
           smooth_stats->num_flat_no_imp += 1;
           break;
        case OPTMS_STEP_TOO_SMALL:
           smooth_stats->num_step_too_small += 1;
           break;
        case OPTMS_MAX_ITER_EXCEEDED:
           smooth_stats->num_max_iter_exceeded += 1;
           break;
        default:
           fprintf(stderr,"An invalid termination status\n");
           OPTMS_WRITE_BINARY_ORDERED_PTS(local_mesh);
        }
    }

    /* active value in local submesh */
    smooth_stats->avg_active_val += local_mesh->current_active_value;

    /* global active value in the mesh */
    if (smooth_stats->global_minimum_val > local_mesh->current_active_value) {
      smooth_stats->global_minimum_val = local_mesh->current_active_value;
    }

    /* improvement in local submesh */
    imp = local_mesh->current_active_value - local_mesh->original_value;
    smooth_stats->avg_improvement += imp;
        
    /* was there actually any improvement */
    if (imp == 0) {
	smooth_stats->no_improvement += 1;
    }
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMprintStats"
int  SMprintStats(SMsmooth_data *smooth_data)
{
    int ierr = 0;
    SMstats *smooth_stats;
    int num_cells,num_cells_laplaced,num_cells_opted;
    double min_value;

    smooth_stats = smooth_data->smooth_stats;

    fprintf(stdout,"******************************************************************\n");
    fprintf(stdout,"SMOOTHING STATISTICS\n");
    fprintf(stdout,"******************************************************************\n");
    num_cells = smooth_stats->total_cells_smoothed;
    num_cells_laplaced = smooth_stats->num_cells_laplaced;
    num_cells_opted = smooth_stats->num_cells_opted;
    min_value = smooth_stats->global_minimum_val;
    ierr = SMconvertToDegrees(smooth_data->smooth_param->function_id,&min_value);
           OPTMS_CHKERR(ierr);

    printf("The approximate global minimum value                        %f \n",
            min_value);
    printf("The total number of nodes smoothed                          %d \n",
            num_cells);
    printf("The number of calls to Laplacian smoothing                  %d \n",
           num_cells_laplaced);
    if (num_cells_laplaced != 0) {
    printf("     Invalid Laplacian steps (percent)                      %d  (%f) \n",
             smooth_stats->num_lap_invalid,
             (100*smooth_stats->num_lap_invalid)/(float)num_cells_laplaced);
    printf("     Non-improvement Laplacian steps (percent)              %d  (%f) \n",
             smooth_stats->num_lap_worse,
             (100*smooth_stats->num_lap_worse)/(float)num_cells_laplaced);
    }
    printf("The number of calls to optimization smoothing               %d \n",
            num_cells_opted);
    if (num_cells_opted != 0) {
        printf("     Average number of iterations/optimization call     %f \n",
                smooth_stats->opt_iter_count/((float)num_cells_opted));
    }
    printf("The number of cells with no improvement                     %d \n",
                              smooth_stats->no_improvement);
    if (num_cells > 0) {
        printf("The average final active value                              %f \n",
                          (smooth_stats->avg_active_val/(float)num_cells));
	printf("The average improvement (over all cells)                    %f \n",
                          (smooth_stats->avg_improvement/(float)num_cells));
    }
    printf("The termination status: \n");
    printf("  Laplacian smoothing enough (percentage)     %d (%f) \n",
            smooth_stats->num_lap_enough,
            100*smooth_stats->num_lap_enough/(float)num_cells);
    printf("  Equilibrium pt found (percentage)           %d (%f) \n",
            smooth_stats->num_equil,100*smooth_stats->num_equil/(float)num_cells);
    printf("     Started at equilibrium (percentage)      %d (%f) \n",
            smooth_stats->num_started_equil,100*smooth_stats->num_started_equil/(float)num_cells);
    printf("  Zero search direction (percentage)          %d (%f) \n",
            smooth_stats->num_zero_search,100*smooth_stats->num_zero_search/(float)num_cells);
    printf("  Improvement too small (percentage)          %d (%f) \n",
            smooth_stats->num_imp_too_small,100*smooth_stats->num_imp_too_small/(float)num_cells);
    printf("  Flat no improvement (percentage)            %d (%f) \n",
            smooth_stats->num_flat_no_imp,100*smooth_stats->num_flat_no_imp/(float)num_cells);
    printf("  Step too small (percentage)                 %d (%f) \n",
            smooth_stats->num_step_too_small,100*smooth_stats->num_step_too_small/(float)num_cells);
    printf("  Max Iter Exceeded (percentage)              %d (%f) \n",
            smooth_stats->num_max_iter_exceeded,100*smooth_stats->num_max_iter_exceeded/(float)num_cells);
    fprintf(stdout,"******************************************************************\n");
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMwriteStatsFile"
int SMwriteStatsFile(SMstats *smooth_stats, int smooth_count)
{
    int ierr;
    FILE *fp;
    int  num_cells;

    num_cells = smooth_stats->total_cells_smoothed;

    if ((fp = fopen("OptMS.stats","a")) == NULL) {
         OPTMS_SETERR(OPTMS_FILE_OPEN_ERR,0,"Can't open OptMS.stats for appending\n");
    } 
    fprintf(fp,"\n");
    fprintf(fp,"******************************************************************\n");
    fprintf(fp,"                            ITERATION  %d                        \n",smooth_count);
    fprintf(fp,"******************************************************************\n");
    fprintf(fp,"  The number of cells smoothed                            %d \n",num_cells);
    fprintf(fp,"  The number of equilibrium pts found                     %d \n",smooth_stats->num_equil);
    fprintf(fp,"  The number of cells that started with an equil point    %d \n",smooth_stats->num_started_equil);
    if (smooth_stats->num_cells_opted == 0) {
        fprintf(fp,"  The average optimization count is 0\n");
    } else {
        fprintf(fp,"  The average optimization iteration count                %f \n",
                  smooth_stats->opt_count/((float)smooth_stats->num_cells_opted));
    }
    fprintf(fp,"  The number of cells with no improvement                 %d \n",smooth_stats->no_improvement);
    fprintf(fp,"  The average active value                                %f \n",
                          (smooth_stats->avg_active_val/(float)num_cells));
    fprintf(fp,"  The average improvement                                 %f \n",
                          (smooth_stats->avg_improvement/(float)num_cells));
    fprintf(fp,"******************************************************************\n");
    fclose(fp);
    return(ierr=0);
}













