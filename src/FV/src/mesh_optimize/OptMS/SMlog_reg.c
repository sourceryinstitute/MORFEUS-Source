/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#include <stdio.h>
#include "SMsmooth.h"

#undef __FUNC__
#define __FUNC__ "SMregisterEvents"
int SMregisterEvents()
{
    int ierr;
#ifdef SUMAA_LOG
    /* register the smoothing events */
    SUMAA_LOG_EVENT_REGISTER(__SM_TOTAL__,"Total Time");
    SUMAA_LOG_EVENT_REGISTER(__SM_SMOOTH__,"Smooth");
    SUMAA_LOG_EVENT_REGISTER(__SM_SMOOTH_INIT__,"Init Smooth");
    SUMAA_LOG_EVENT_REGISTER(__SM_INIT__,"Initialize");
    SUMAA_LOG_EVENT_REGISTER(__SM_SMOOTH_OPT__,"Optimize");
    SUMAA_LOG_EVENT_REGISTER(__SM_FUNCTION__,"Function");
    SUMAA_LOG_EVENT_REGISTER(__SM_GRADIENT__,"Gradient");
    SUMAA_LOG_EVENT_REGISTER(__SM_FIND_ACTIVE__,"Find Active Set");
    SUMAA_LOG_EVENT_REGISTER(__SM_COMP_ALPHA__,"Comp Alpha");
    SUMAA_LOG_EVENT_REGISTER(__SM_SEARCH__,"Search Dir");
    SUMAA_LOG_EVENT_REGISTER(__SM_STEP_ACCEPT__,"Step Accept");
    SUMAA_LOG_EVENT_REGISTER(__SM_MIN_EST__,"Est Min Imp");
    SUMAA_LOG_EVENT_REGISTER(__SM_EDGE_FACE__,"Edge/Face Srch");
    SUMAA_LOG_EVENT_REGISTER(__SM_CHK_EQUIL__,"Chk Equil");
    SUMAA_LOG_EVENT_REGISTER(__SM_CUSP__,"Step Cusp");
    SUMAA_LOG_EVENT_REGISTER(__SM_SMOOTH_FINAL__,"End Smooth");
    SUMAA_LOG_EVENT_REGISTER(__SM_VALID__,"Valid");
    SUMAA_LOG_EVENT_REGISTER(__SM_INIT_STATS__,"Init Stats");
    SUMAA_LOG_EVENT_REGISTER(__SM_GRAD_PROJ__,"Grad Proj");
    SUMAA_LOG_EVENT_REGISTER(__SM_GET_ACTIVE__,"Get Active Dir");
    SUMAA_LOG_EVENT_REGISTER(__SM_VERT_STEP__,"Vert Step");
    SUMAA_LOG_EVENT_REGISTER(__SM_FORM_GRAM__,"Form Grammian");
    SUMAA_LOG_EVENT_REGISTER(__SM_FORM_PDG__,"Form PDG");
    SUMAA_LOG_EVENT_REGISTER(__SM_FORM_REDUCED__,"Form Reduced");
    SUMAA_LOG_EVENT_REGISTER(__SM_SOLVE_2__,"Solve 2x2");
    SUMAA_LOG_EVENT_REGISTER(__SM_SOLVE_3__,"Solve 3x3");
    SUMAA_LOG_EVENT_REGISTER(__SM_NONSING_TEST__,"NonSing Test");
    SUMAA_LOG_EVENT_REGISTER(__SM_COPY_FCN__,"Copy Fcn");
    SUMAA_LOG_EVENT_REGISTER(__SM_COPY_ACT__,"Copy Act");
    SUMAA_LOG_EVENT_REGISTER(__SM_LAP_SMOOTH__,"Lap Smooth");
    SUMAA_LOG_EVENT_REGISTER(__SM_UNTANGLE__,"Untangle");
    SUMAA_LOG_EVENT_REGISTER(__SM_PHASE1__,"Phase 1 Soln");
    SUMAA_LOG_EVENT_REGISTER(__SM_LINEAR_PROG__,"LP Solve");
    SUMAA_LOG_EVENT_REGISTER(__SM_LP_ITER__,"LP Iter");
#endif
    return(ierr=0);
}
