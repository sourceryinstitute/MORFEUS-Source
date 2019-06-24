!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
!
!    NEMO - Numerical Engine (for) Multiphysics Operators
! Copyright (c) 2007, Stefano Toninel
!                     Gian Marco Bianchi  University of Bologna
!              David P. Schmidt    University of Massachusetts - Amherst
!              Salvatore Filippone University of Rome Tor Vergata
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
! $Id: wr_mtx_vector.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Writes a global dense vector in Matrix Market format
!
SUBMODULE (tools_output_basics) wr_mtx_vector_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE wr_mtx_vector
            USE class_psblas
            IMPLICIT NONE
            !
            INTEGER, PARAMETER :: vector = 10
            !
            INTEGER :: err_act, info
            INTEGER :: i, ncells_glob
            REAL(psb_dpk_), ALLOCATABLE :: glob_vect(:)

            CALL sw_out%tic()

            ! Sets error handling for PSBLAS-2 routines
            info = 0
            CALL psb_erractionsave(err_act)

            IF(mypnum_() == 0) THEN
                WRITE(*,*) 'Dumping dense  vector on file: ', TRIM(name)
                OPEN(unit=vector,file=name)
            END IF

            ncells_glob = psb_cd_get_global_cols(desc)

            ALLOCATE(glob_vect(ncells_glob),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            CALL psb_gather(glob_vect,loc_vect,desc,info,root=0)
            CALL psb_check_error(info,'wr_mtx_vector','psb_gather',icontxt_())


            IF(mypnum_() == 0) THEN
                ! Writes GLOB_VECT
                WRITE(vector,*) ncells_glob, 1, ncells_glob
                DO i = 1, ncells_glob
                    WRITE(vector,*) i, 1, glob_vect(i)
                END DO

                ! Closes file on P0
                CLOSE(vector)
            END IF

            ! Cleanups storage
            DEALLOCATE(glob_vect)

            ! Normal termination
            CALL psb_erractionrestore(err_act)

            CALL sw_out%toc()

100         FORMAT(' ERROR! Memory allocation failure in WR_MTX_VECTOR')

      END PROCEDURE wr_mtx_vector

END SUBMODULE wr_mtx_vector_implementation
