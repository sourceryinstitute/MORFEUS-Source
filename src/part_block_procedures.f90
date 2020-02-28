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
! $Id: part_block.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Block partitioning
!
SUBMODULE(part_block) part_block_procedures

    ! USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE bld_part_block
        !
        INTEGER :: dim_block, i, info
        INTEGER :: icontxt, mypnum
        REAL(psb_dpk_) :: load(nprocs)

        ! No need INTENT(OUT)
!!$    if(associated(part) .and. size(part) /= ncells) then
!!$       deallocate(part)
!!$       part => null()
!!$    end if
!!$
!!$    if(.not.associated(part)) then
        ALLOCATE(part(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
!!$    end if

        icontxt = icontxt_()
        mypnum = mypnum_()

        IF(mypnum == 0) THEN
            WRITE(*,200)
            WRITE(*,300,advance='no') 'BLOCK'
        END IF

        ! BLOCK partitioning
        dim_block = (ncells + nprocs - 1) / nprocs
        DO i = 1, ncells
            part(i) = (i - 1) / dim_block
        END DO

        ! Computes load of single processes
        load = 0.d0
        DO i=1,nprocs
            load(i) = COUNT(part==i-1)
            load(i) = REAL(load(i)) / REAL(ncells) * 100
        ENDDO
        WRITE(*,400) load



100     FORMAT(' ERROR! Memory allocation failure in BLD_PART_BLOCK')
200     FORMAT(' Domain Partitioning',5x,'Load')
300     FORMAT(1x,a)
400     FORMAT(16x,16(3x,f5.1,'%'))

    END PROCEDURE bld_part_block

END SUBMODULE part_block_procedures
