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

#define SMinitActiveInfo(max_pts, active) \
{  \
    int i99; \
    active->num_active = 0; \
    active->num_equal = 0; \
    active->true_active_value = 0; \
    for (i99=0;i99<max_pts;i99++) { \
        active->active_ind[i99] = 0; \
    } \
}

#undef __FUNC__
#define __FUNC__ "SMfindActiveSet"
int SMfindActiveSet(int num_values, double *function, double active_eps,
		     SMactive *active_info)
{
    int         ierr;
    int         i, ind;
    int         num_active, num_equal;
    double      function_val;
    double      active_value0;
    double      true_active;
    double      temp;

    OPTMS_CHECK_NULL(active_info);
    if (num_values > OPTMS_MAX_NUM_TRI*OPTMS_DEFAULT_FUNC_PER_TRI)
       OPTMS_SETERR(OPTMS_INPUT_ERR,0,"Number of function values exceeds maximum");

    active_info->num_active = 0;
    active_info->num_equal = 0;
    active_info->true_active_value = 0;
    for (i=0;i<num_values;i++) {
        active_info->active_ind[i] = 0;
    }

    /* the first function value is our initial active value */
    num_active = 1;
    num_equal = 0;
    active_info->active_ind[0] = 0;
    true_active = function[0];

    /* first sort out the active set...
       all vals within active_eps of smallest val */

    for (i=1;i<num_values;i++) {
	function_val = function[i];
        true_active = OPTMS_MIN(function_val,true_active);
	active_value0 = function[active_info->active_ind[0]];
	temp = fabs(function_val - active_value0);
	if ( function_val < active_value0 ) {
	    if ( temp > active_eps) {
		num_active = 1;
		num_equal = 0;
		active_info->active_ind[0] = i;
	    } else if ( temp < active_eps) {
		num_active += 1;
		ind = num_active - 1;
		active_info->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < OPTMS_MACHINE_EPS) {
		    num_equal += 1;
		}
	    }
	} else {
	    if (temp < active_eps) {
		num_active += 1;
		ind = num_active - 1;
		active_info->active_ind[ind] = i;
		if (fabs(function_val - active_value0) < OPTMS_MACHINE_EPS) {
		    num_equal += 1;
		}
	    }
	}
    }
    active_info->true_active_value = true_active;
    active_info->num_active = num_active;
    active_info->num_equal  = num_equal;
    return(ierr=0);
}

#undef __FUNC__
#define __FUNC__ "SMprintActiveSet"
int SMprintActiveSet(SMactive *active, double *function)
{
    int ierr;
    int  i;

    OPTMS_CHECK_NULL(active);

    /* print the active set */
    for (i=0;i<active->num_active;i++) {
       OPTMS_DEBUG_ACTION(2,{
           fprintf(stdout,"Active value %d:   %f \n",
	               i+1,function[active->active_ind[i]]);
       });
    }
    return(ierr=0);
}
