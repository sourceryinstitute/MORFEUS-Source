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
! $Id: part_graph.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    METIS partitioning
!
SUBMODULE(part_graph) part_graph_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE get_part_graph
        !
        INTEGER :: info, ncells

        ! Check on internal partitioning array PART_GLOB
        IF (.NOT.ALLOCATED(part_glob)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ncells = SIZE(part_glob)

        ALLOCATE(part(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        part(:) = part_glob(:)

        ! Deallocates memory
        CALL free_part_graph

100     FORMAT(' ERROR! Partitioning vector has not yet been built')
200     FORMAT(' ERROR! Memory allocation failure in GET_PART_GRAPH')

    END PROCEDURE get_part_graph


    MODULE PROCEDURE bld_part_graph
        !
        ! PSBLAS
        INTEGER :: icontxt, mypnum, nprocs
        !
        ! METIS
        INTEGER, ALLOCATABLE  :: vwgt(:), adjwgt(:)
        INTEGER, PARAMETER ::  wgtflag=0, numflag=1
        INTEGER :: nparts
        INTEGER :: options(0:4)
        INTEGER :: edgecut
        INTEGER, ALLOCATABLE :: part_loc(:)
        CHARACTER(len=32) :: methd
        !
        INTEGER :: i, info, j, n, ncells
        INTEGER, ALLOCATABLE :: work(:)
        REAL(psb_dpk_), ALLOCATABLE :: load(:)

        !-------------------------------------------------------------------

        ncells = SIZE(xadj_glob) - 1

        nprocs = nprocs_()

        ! Parallel job: preliminary checks
        IF    (ALLOCATED(part_glob) .OR. &
            & ALLOCATED(xadj)      .OR. &
            & ALLOCATED(adjncy)    .OR. &
            & ALLOCATED(perm)      .OR. &
            & ALLOCATED(pinv)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Allocates pointers PART_GLOB, PERM, PINV
        ALLOCATE(part_glob(ncells),perm(ncells),pinv(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,500)
            CALL abort_psblas
        END IF

        ! Trivial case: SERIAL job
        IF(nprocs == 1) THEN
            DO i = 1, ncells
                part_glob(i) = 0
            END DO
            RETURN
        END IF

        n = SIZE(xadj_glob) - 1
        nparts = nprocs
        options(0) = 0

        ALLOCATE(part_loc(n),work(ncells),load(nprocs),stat=info)
        IF(info /= 0) THEN
            WRITE(*,500)
            CALL abort_psblas
        END IF


        ! METIS partitioning
        methd = 'METIS_PartGraphKway'
        WRITE(*,600)
        WRITE(*,700,advance='no') TRIM(methd)

        CALL METIS_PartGraphKway(n,xadj_glob,adjncy_glob,vwgt,adjwgt,&
            & wgtflag,numflag,nparts,options,edgecut,part_glob)
        !
        load = 0.d0
        DO i=1,nprocs
            load(i) = COUNT(part_glob==i)
            load(i) = REAL(load(i)) / REAL(ncells) * 100
        ENDDO
        WRITE(*,800) edgecut, load


        ! Processes numbering starts from 0
        DO i = 1, ncells
            part_glob(i) = part_glob(i) - 1
        END DO

        DEALLOCATE(part_loc,work,load)

100     FORMAT(' ERROR! Pointers status in BLD_PART_GRAPH incompatible with IPART = 1')
200     FORMAT(' ERROR! Pointers status in BLD_PART_GRAPH incompatible with IPART > 1')
300     FORMAT(' ERROR! Partitioning of adaptive remeshing not yet supported')
400     FORMAT(' ERROR! Partitioning method not supported')
500     FORMAT(' ERROR! Memory allocation failure in BLD_PART_GRAPH')
600     FORMAT(' Domain Partitioning',6x,' Edgecut',2x,' Load')
700     FORMAT(1x,a)
800     FORMAT(i8,16(3x,f5.1,'%'))


    END PROCEDURE bld_part_graph


    SUBROUTINE free_part_graph

        DEALLOCATE(part_glob)
        IF(ALLOCATED(xadj))    DEALLOCATE(xadj)
        IF(ALLOCATED(adjncy))  DEALLOCATE(adjncy)
        IF(ALLOCATED(perm))    DEALLOCATE(perm)
        IF(ALLOCATED(pinv))    DEALLOCATE(pinv)

    END SUBROUTINE free_part_graph


END SUBMODULE part_graph_procedures
