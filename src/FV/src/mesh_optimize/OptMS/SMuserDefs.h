/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_USER_DEFS_H
#define SM_USER_DEFS_H

#define OPTMS_DEFAULT -1

/* 2D function/gradient options included with Opt-MS
   The default is MAX_MIN_SINE */
#define OPTMS_MAX_MIN_ANGLE                     1
#define OPTMS_MIN_MAX_COSINE                    2
#define OPTMS_MAX_MIN_COSINE                    3
#define OPTMS_MAX_MIN_SINE                      4
#define OPTMS_MIN_MAX_ANGLE                     5
#define OPTMS_MIN_MAX_JACOBIAN_DIFF             6
#define OPTMS_MAX_MIN_SCALED_JACOBIAN           7
#define OPTMS_MAX_MIN_AREA_LENGTH_RATIO         8
#define OPTMS_MIN_MAX_LENGTH_AREA_RATIO         9
#define OPTMS_MAX_MIN_INTERIOR_ANGLE           10
#define OPTMS_MAX_MIN_INTERIOR_SINE            11
#define OPTMS_MIN_MAX_INTERIOR_COSINE          12
#define OPTMS_MAX_MIN_INTERIOR_SCALED_JACOBIAN 13
#define OPTMS_MIN_MAX_NORM_JAC_SQUARED_2D      14
#define OPTMS_MIN_MAX_CONDITION_2D             15
#define OPTMS_FUNCTION2D_DEFAULT                4

/* 3D function/gradient options included with Opt-MS
   The default is MAX_SINE_DIHEDRAL */
#define OPTMS_MAX_MIN_DIHEDRAL                 21
#define OPTMS_MIN_MAX_DIHEDRAL                 22
#define OPTMS_MAX_MIN_COSINE_DIHEDRAL          23
#define OPTMS_MIN_MAX_COSINE_DIHEDRAL          24
#define OPTMS_MAX_SINE_DIHEDRAL                25
#define OPTMS_MAX_MIN_SCALED_JACOBIAN_3D       26
#define OPTMS_MIN_MAX_SRMS_VOLUME_RATIO        27
#define OPTMS_MIN_MAX_CONDITION_3D             28
#define OPTMS_MIN_MAX_NORM_JAC_SQUARED_3D      29
#define OPTMS_FUNCTION3D_DEFAULT               25

/* user interface to the smoothing techniques */
#define OPTMS_LAPLACIAN_ONLY          1
#define OPTMS_SMART_LAPLACIAN_ONLY    2
#define OPTMS_OPTIMIZATION_ONLY       3
#define OPTMS_COMBINED                4
#define OPTMS_COMBINED1               5
#define OPTMS_COMBINED2               6
#define OPTMS_COMBINED3               7
#define OPTMS_FLOATING_THRESHOLD      8
#define OPTMS_TECHNIQUE_DEFAULT       4

/* untangling techniques */
#define OPTMS_LAPLACIAN_ONLY          1
#define OPTMS_LINEAR_PROGRAM_ONLY     2
#define OPTMS_COMBINED_UNTANGLING     3
#define OPTMS_UNTANGLING_DEFAULT      2

#endif
