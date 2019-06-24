/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_LOG_H
#define SM_LOG_H 1

#if defined(SUMAA_LOG)

#define SM_LOG_current                SUMAA_LOG_current
#define SM_LOG_flops                  SUMAA_LOG_total_flops
#define SM_LOG_time                   SUMAA_LOG_time
#define SM_LOG_initialized            SUMAA_LOG_initialized
#define SM_LOG_file                   SUMAA_LOG_file
#define SM_LOG_event                  SUMAA_LOG_event

#define SM_LOG_INIT(a,b)              SUMAAlogInit(a,b)
#define SM_REGISTER_LOG_EVENTS()      SMregisterEvents()
#define SM_LOG_EVENT_REGISTER(a,str)  SUMAAlogEventRegister(a,str)
#define SM_LOG_EVENT_BEGIN(a)         SUMAAlogEventBegin(a)
#define SM_LOG_EVENT_END(a)           SUMAAlogEventEnd(a)
#define SM_LOG_FLOPS(a,b)             SUMAA_LOG_FLOPS(a,b)
#define SM_LOG_GLOBAL_TIME(a)         SUMAA_LOG_GLOBAL_TIME(a)
#ifdef PARALLEL_LOG
   #define SM_LOG_PRINT(a,b,c)        SUMAAlogPrintParallel(a,b,c)
#else
   #define SM_LOG_PRINT()             SUMAAlogPrintSingle()
#endif
#define SM_LOG_FREE                   SUMAAlogFree()

#else

#define SM_LOG_current
#define SM_LOG_flops
#define SM_LOG_time
#define SM_LOG_initialized
#define SM_LOG_file
#define SM_LOG_event

#define SM_LOG_INIT(a,b)
#define SM_REGISTER_LOG_EVENTS()
#define SM_LOG_EVENT_REGISTER(a,str)
#define SM_LOG_EVENT_BEGIN(a)
#define SM_LOG_EVENT_END(a)
#define SM_LOG_FLOPS(a,b)
#define SM_LOG_GLOBAL_TIME(a)
#define SM_LOG_PRINT()
#define SM_LOG_FREE

#endif

/* log floating point operations */
#ifdef SM_LOG
#define SM_LOG_flop(num) SM_LOG_flops += num
#else
#define SM_LOG_flop(num)
#endif

/* the set of logging events for Smoothing code */
#ifdef SUMAA_LOG
#define __SM_SMOOTH__        0
#define __SM_SMOOTH_INIT__   1
#define __SM_INIT__          2
#define __SM_SMOOTH_OPT__    3
#define __SM_FUNCTION__      4
#define __SM_GRADIENT__      5
#define __SM_FIND_ACTIVE__   6
#define __SM_COMP_ALPHA__    7
#define __SM_SEARCH__        8
#define __SM_STEP_ACCEPT__   9
#define __SM_MIN_EST__       10
#define __SM_EDGE_FACE__     11
#define __SM_CHK_EQUIL__     12
#define __SM_CUSP__          13
#define __SM_SMOOTH_FINAL__  14
#define __SM_VALID__         15
#define __SM_INIT_STATS__    16
#define __SM_GRAD_PROJ__     17
#define __SM_GET_ACTIVE__    18
#define __SM_VERT_STEP__     19
#define __SM_FORM_GRAM__     20
#define __SM_FORM_PDG__      21
#define __SM_FORM_REDUCED__  22
#define __SM_SOLVE_2__       23
#define __SM_SOLVE_3__       24
#define __SM_NONSING_TEST__  25
#define __SM_COPY_FCN__      26
#define __SM_COPY_ACT__      27
#define __SM_LAP_SMOOTH__    28
#define __SM_UNTANGLE__    29
#define __SM_PHASE1__    30
#define __SM_LINEAR_PROG__    31
#define __SM_LP_ITER__    32
#define __SM_TOTAL__    33
#endif

#endif
