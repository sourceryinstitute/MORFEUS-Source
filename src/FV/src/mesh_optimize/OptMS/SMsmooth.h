/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
#ifndef SMOOTH_H
#define SMOOTH_H 1

/* define the structures associated with both geometric and octree
   representations of geometric entities... the void data type allows
   us to switch contexts
*/

#ifdef SUMAA_LOG
#include "SUMAA_log.h"
#endif
#if defined (IRIX)||(IRIX32)||(IRIX64)
#include <stdlib.h>
#include <unistd.h>
#endif
#include "SMqual_func.h"
#include "SMerror.h"
#include "SMdefs.h"
#include "SMinternalFunction.h"
#include "SMlog.h"
#include "OptMS.h"

#endif
