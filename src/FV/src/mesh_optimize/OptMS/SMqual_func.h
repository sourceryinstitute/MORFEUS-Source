/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SM_QUAL_FUNC_H
#define SM_QUAL_FUNC_H 1

typedef int (*SMfunction_ptr2D)(double *, double *, double *, double *, int *);
typedef int (*SMgradfunc_ptr2D)(double *, double *, double *, double **, int *);

typedef int (*SMfunction_ptr3D)(const double adCoord0[3],
                        const double adCoord1[3], const double adCoord2[3],
                        const double adCoord3[3], double *adResult,
                        int* const piNGradient);

typedef int (*SMgradfunc_ptr3D)(const double adCoord0[3],
                        const double adCoord1[3],const double adCoord2[3],
                        const double adCoord3[3],double **adGradient,
                        int* const piNGradient);

#endif
