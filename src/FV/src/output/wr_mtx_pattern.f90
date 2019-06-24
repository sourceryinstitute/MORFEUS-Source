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
! $Id: wr_mtx_pattern.f90 3106 2008-05-05 15:07:27Z sfilippo $
!
! Description:
!    Dumps the pattern of a sparse matrix in Matrix Market format, starting
!    from the global C2C adjacency graph.
!
MODULE PROCEDURE wr_mtx_pattern
    USE class_psblas
    USE class_connectivity

    IMPLICIT NONE
    !
    INTEGER, PARAMETER :: pattern=10
    INTEGER :: i, j, n, ncells, nnz
    INTEGER, POINTER :: iconn(:) => NULL()
    CHARACTER(len=*), PARAMETER :: fmt='(3(i8,1x))'

    CALL tic(sw_out)

    IF(mypnum_() == 0) THEN
        WRITE(*,'(a)',advance='no') ' Now dumping sparsity pattern...'

        OPEN(unit=pattern,file=name)

        ncells = nel_(c2c)

        ! REMARK! The adjacency graph C2C contains only the non-diagonal elements.
        ! One has to add also the elements of the main diagonal in order to get
        ! the right number of non-zero coefficients.
        nnz = nconn_(c2c) + ncells

        WRITE(pattern,fmt) ncells, ncells, nnz

        DO i = 1, ncells

            CALL get_ith_conn(iconn,c2c,i)
            n = SIZE(iconn)

            ! Off-diagonal elements -> 1
            DO j = 1, n
                WRITE(pattern,fmt) i, iconn(j), 1
            END DO

            ! Main diagonal elements -> -1
            WRITE(pattern,fmt) i, i, -1
        END DO

        NULLIFY(iconn)

        WRITE(*,*) 'done.'
        WRITE(*,*)

        CLOSE(pattern)
    END IF

    CALL toc(sw_out)

END PROCEDURE wr_mtx_pattern
