/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/
/*
  !
  !	NEMO - Numerical Engine (for) Multiphysics Operators
  ! Copyright (c) 2007, Stefano Toninel
  !                     Gian Marco Bianchi  University of Bologna
  !		      David P. Schmidt    University of Massachusetts - Amherst
  !		      Salvatore Filippone University of Rome Tor Vergata
  ! All rights reserved.
  !
  ! Redistribution and use in source and binary forms, with or without modification,
  ! are permitted provided that the following conditions are met:
  !
  !     1. Redistributions of source code must retain the above copyright notice,
  !        this list of conditions and the following disclaimer.
  !     2. Redistributions in binary form must reproduce the above copyright notice,
  !        this list of conditions and the following disclaimer in the documentation
  !        and/or other materials provided with the distribution.
  !     3. Neither the name of the NEMO project nor the names of its contributors
  !        may be used to endorse or promote products derived from this software
  !        without specific prior written permission.
  !
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  ! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  ! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
  ! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  ! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  ! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
  ! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  ! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !
  !---------------------------------------------------------------------------------
  !
  ! $Id: right_handed.c 333 2007-03-14 10:39:08Z stoninel $
  !
  ! Description:  checks that a triangle is right-handed
  !
*/


#include <stdlib.h>
#include <stdio.h>
#include "OptMS.h"

#ifdef LowerCase
#define right_handed2d    right_handed2d
#endif
#ifdef LowerUnderscore
#define right_handed2d    right_handed2d_
#endif
#ifdef LowerDoubleUnderscore
#define right_handed2d    right_handed2d__
#endif
#ifdef UpperCase
#define right_handed2d    RIGHT_HANDED2D
#endif
#ifdef UpperUndescore
#define right_handed2d    RIGHT_HANDED2D_
#endif
#ifdef UpperDoubleUnderscore
#define right_handed2d    RIGHT_HANDED2D__
#endif

/* prototypes */

int SMorient2D(double *vtx1, double *vtx2, double *free_vtx, int *valid);

int right_handed2d(double *vtx1, double *vtx2, double *vtx3)
{
  int ierr;
  int handedness;

  ierr=SMorient2D(vtx1,vtx2,vtx3,&handedness);

  return (handedness);

}
