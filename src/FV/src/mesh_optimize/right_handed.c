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
  ! $Id: right_handed.c 3064 2008-04-11 16:29:54Z sfilippo $
  !
  ! Description:  checks that a tet is right-handed
  ! 
*/


#include <stdlib.h>
#include <stdio.h>
#include "OptMS.h"

#ifdef LowerCase
#define right_handed    right_handed
#endif
#ifdef LowerUnderscore
#define right_handed    right_handed_
#endif
#ifdef LowerDoubleUnderscore
#define right_handed    right_handed__
#endif
#ifdef UpperCase
#define right_handed    RIGHT_HANDED
#endif
#ifdef UpperUndescore
#define right_handed    RIGHT_HANDED_
#endif
#ifdef UpperDoubleUnderscore
#define right_handed    RIGHT_HANDED__
#endif

/* prototypes */

int SMorient3D(double *vtx1, double *vtx2, double *vtx3, double *free_vtx, int *valid);

int right_handed(double *vtx1, double *vtx2, double *vtx3, double *vtx4)
{
  int ierr;
  int handedness;
  
  ierr=SMorient3D(vtx1,vtx2,vtx3,vtx4,&handedness);
   
  return (handedness);
  
}
