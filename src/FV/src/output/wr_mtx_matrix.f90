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
! $Id: wr_mtx_matrix.f90 9092 2015-04-24 08:50:00Z sfilippo $
!
! Description:
!    Writes a global sparse matrix A in Matrix Market format.
!    WARNING! Only CSR is currently supported.
!
SUBMODULE (tools_output_basics) wr_mtx_matrix_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE wr_mtx_matrix
            USE class_psblas
            IMPLICIT NONE
            !
            INTEGER, PARAMETER :: matrix = 10
            !
            INTEGER :: info, err_act
            INTEGER :: icontxt, mypnum, nprocs
            INTEGER :: dim, i, j, ncells_glob, ncells_loc, nnz,nrloc
            INTEGER, ALLOCATABLE  :: ia(:), ja(:)
            INTEGER, ALLOCATABLE :: iloc_to_glob(:)
            REAL(psb_dpk_), ALLOCATABLE  :: acoo(:)

            ! Sets error handling for PSBLAS-2 routines
            info=0
            CALL psb_erractionsave(err_act)

            icontxt = icontxt_()
            mypnum = mypnum_()
            nprocs = nprocs_()

            CALL sw_out%tic()

            IF(mypnum == 0) THEN
                WRITE(*,*) 'Dumping sparse matrix on file: ', TRIM(name)
                OPEN(unit=matrix,file=name)
            END IF

            ncells_glob = desc%get_global_cols()

            ! Gets local to global list for cell indices
            CALL psb_get_loc_to_glob(desc,iloc_to_glob)

            ! Computes # of non-zero elements of local/global matrix
            ncells_loc = desc%get_local_rows()
            nnz        = a%get_nzeros()
            CALL psb_sum(icontxt,nnz)
            ! Opens file and writes header
            IF (mypnum == 0) THEN
                OPEN(unit=matrix,file=name)
                WRITE(matrix,*) ncells_glob, ncells_glob, nnz
            END IF

            ! NNZ = local also on P0
            nnz = A%get_nzeros()
            ! Allocates IA, JA, ACOO for storing non-zero elements in COO format
            dim = nnz
            CALL psb_amx(icontxt,dim)
            ALLOCATE(ia(dim),ja(dim),acoo(dim),stat=info)
            IF(info /= 0) THEN
                WRITE (*,100)
                CALL abort_psblas
            END IF

            nrloc = ncells_loc
            CALL a%csget(1,nrloc,nnz,ia,ja,acoo,info)
            CALL psb_loc_to_glob(ia(1:nnz),desc,info,'E')
            CALL psb_check_error(info,'wr_mtx_matrix','psloc_to_glob',icontxt)
            CALL psb_loc_to_glob(ja(1:nnz),desc,info,'E')
            CALL psb_check_error(info,'wr_mtx_matrix','psloc_to_glob',icontxt)

            DO i = 1, nprocs - 1
                IF(mypnum == i) THEN
                    CALL psb_snd(icontxt,nnz,0)
                    CALL psb_snd(icontxt,ia(1:nnz),0)
                    CALL psb_snd(icontxt,ja(1:nnz),0)
                    CALL psb_snd(icontxt,acoo(1:nnz),0)
                ELSE IF(mypnum == 0) THEN
                    CALL psb_rcv(icontxt,nnz,i)
                    CALL psb_rcv(icontxt,ia(1:nnz),i)
                    CALL psb_rcv(icontxt,ja(1:nnz),i)
                    CALL psb_rcv(icontxt,acoo(1:nnz),i)
                    DO j = 1, nnz
                        WRITE(matrix,*) ia(j), ja(j), acoo(j)
                    END DO
                END IF
            END DO

            ! Closes file on P0
            IF(mypnum == 0) CLOSE(matrix)

            ! Cleanups storage
            DEALLOCATE(ia,ja,acoo)
            DEALLOCATE(iloc_to_glob)

            ! Normal termination
            CALL psb_erractionrestore(err_act)

            CALL sw_out%toc()

100         FORMAT(' ERROR! Memory allocation failure in WR_MTX_MATRIX')

    END PROCEDURE wr_mtx_matrix

END SUBMODULE wr_mtx_matrix_implementation
