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
  ! $Id: call_smooth.c 3064 2008-04-11 16:29:54Z sfilippo $
  !
  ! Description: calls the smoother in OptMS
  !
*/


#include <stdlib.h>
#include <stdio.h>
#include "OptMS.h"

/* static variable for data structure used in smoothing*/
void *smooth_data;

#ifdef LowerCase
#define call_smooth    call_smooth
#endif
#ifdef LowerUnderscore
#define call_smooth    call_smooth_
#endif
#ifdef LowerDoubleUnderscore
#define call_smooth    call_smooth__
#endif
#ifdef UpperCase
#define call_smooth    CALL_SMOOTH
#endif
#ifdef UpperUndescore
#define call_smooth    CALL_SMOOTH_
#endif
#ifdef UpperDoubleUnderscore
#define call_smooth    CALL_SMOOTH__
#endif


int call_smooth(int *num_incident_vtx, int *num_incident_tet,
		double *free_pos,double incident_vtx[][3],
		int vtx_connectivity[][3],int *tangled)
{
    int ierr=0;
    double **pos_ptr;  /* sized [*num_incident_vtx][3] */
    int **conn_ptr;    /* sized [*num_incident_tet][3] */


    int i;

    pos_ptr=malloc(sizeof(double *)* (*num_incident_vtx));
    conn_ptr=malloc(sizeof(int *)* (*num_incident_tet));

    /* point each row of pos_ptr to point to the corresponding row of the input array */
    for (i=0;i<*num_incident_vtx;i++)
	pos_ptr[i]=&incident_vtx[i][0];

    /* point each row of conn_ptr to point to the corresponding row of the input array */
    for (i=0;i<*num_incident_tet;i++)
	conn_ptr[i]=&vtx_connectivity[i][0];

    if (*tangled > 0) /* untangling should always be followed by smoothing, since untangling
			 leaves the mesh in a valid, but low-quality state */
    {
	ierr = SMuntangle(*num_incident_vtx,*num_incident_tet,
		     free_pos, pos_ptr,
		     conn_ptr, smooth_data);
    }
    else
    {
	ierr = SMsmooth(*num_incident_vtx,*num_incident_tet,
		     free_pos, pos_ptr,
		     conn_ptr, smooth_data);

	if (ierr != 0 )   /* then we probably have a tangled mesh */
	{
	    ierr = 0;  /* try again, with untangling */
	    ierr = SMuntangle(*num_incident_vtx,*num_incident_tet,
		     free_pos, pos_ptr,
		     conn_ptr, smooth_data);


	}
    }

    OPTMS_CHKERR(ierr);

    free(pos_ptr);
    free(conn_ptr);

    return (ierr);

}
