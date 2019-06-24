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
! $Id: rd_inp_bc_math.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(tools_bc) rd_inp_bc_math_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_inp_bc_math
        USE class_psblas
        USE tools_bc, ONLY: bc_dirichlet_, bc_neumann_, bc_robin_, &
            & bc_dirichlet_map_, bc_neumann_map_, bc_robin_map_
        USE tools_input

        IMPLICIT NONE
        !
        LOGICAL, PARAMETER :: debug = .FALSE.
        !
        INTEGER :: mypnum, icontxt
        INTEGER :: info, inp, na, nb, nc
        CHARACTER(len=15) :: par
        REAL(psb_dpk_) :: abc(3)

        icontxt = icontxt_()
        mypnum  = mypnum_()

        IF(mypnum == 0) THEN
            CALL open_file(input_file,inp)
            CALL find_section(sec,inp)

            WRITE(*,*) '- Reading ',TRIM(sec),' section: type MATH'

            READ(inp,'()')

            ! Reads BC section and the corresponding ID
            READ(inp,*) par, id

            ! Broadcasts the ID
            CALL psb_bcast(icontxt,id)

            ! Reads ABC parameters
            SELECT CASE(id)
            CASE(bc_dirichlet_)
                abc(1) = 1.d0
                abc(2) = 0.d0
                READ(inp,*) par, abc(3)
            CASE(bc_neumann_)
                abc(1) = 0.d0
                abc(2) = 1.d0
                READ(inp,*) par, abc(3)
            CASE(bc_robin_)
                READ(inp,*) par, abc(1), abc(2), abc(3)
            CASE default
                abc(:) = 0.d0
            END SELECT

            ! Broadcasts ABC parameters (only uniform BCs)
            CALL psb_bcast(icontxt,abc)

            CLOSE(inp)
        ELSE
            CALL psb_bcast(icontxt,id)
            CALL psb_bcast(icontxt,abc)
        END IF

        ! Allocates A, B, C pointers depending on the BC type
        SELECT CASE(id)
        CASE(  bc_dirichlet_, &
            & bc_neumann_,   &
            & bc_robin_ )
            na = 1; nb = 1; nc = 1
        CASE(  bc_dirichlet_map_, &
            & bc_neumann_map_)
            na = 1; nb = 1; nc = nbf
        CASE(  bc_robin_map_)
            na = nbf; nb = nbf; nc = nbf
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT
        ALLOCATE(a(na),b(nb),c(nc),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ! Sets BC values
        SELECT CASE(id)
        CASE(bc_dirichlet_, bc_neumann_, bc_robin_)
            a(1) = abc(1)
            b(1) = abc(2)
            c(1) = abc(3)
        CASE(bc_dirichlet_map_)
            a(1) = 1.d0
            b(1) = 0.d0
            c(:) = 0.d0 ! To be actually set at the first mapping
        CASE(bc_neumann_map_)
            a(1) = 0.d0
            b(1) = 1.d0
            c(:) = 0.d0 ! To be actually set at the first mapping
        CASE(bc_robin_map_)
            a(:) = 0.d0 ! To be actually set at the first mapping
            b(:) = 0.d0 ! To be actually set at the first mapping
            c(:) = 0.d0 ! To be actually set at the first mapping
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

        ! Debug
        IF(debug) THEN
            WRITE(*,*)
            WRITE(*,300) mypnum
            WRITE(*,400) TRIM(sec),' - Type: Math'
            WRITE(*,500) '  BC%id     = ', id
            WRITE(*,600) '  BC%a = ', a
            WRITE(*,600) '  BC%b = ', b
            WRITE(*,600) '  BC%c = ', c
            WRITE(*,*)
        END IF

100     FORMAT(' ERROR! Unsupported BC in RD_INP_BC_MATH')
200     FORMAT(' ERROR! Memory allocation failure in RD_INP_BC_MATH')
300     FORMAT(' ----- Process ID = ',i2,' -----')
400     FORMAT(1x,a,a)
500     FORMAT(1x,a,i2)
600     FORMAT(1x,a,es10.3)

        END PROCEDURE rd_inp_bc_math

END SUBMODULE rd_inp_bc_math_implementation