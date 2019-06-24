!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
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
! $Id: laplacian_smooth.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description: moves all interior points to the average position of their neighbors.
! This is done by forming three linear systems, one for each of the x,y,z components
! of position and using PSBLAS linear solvers
!
SUBMODULE(tools_mesh_optimize) laplacian_smooth_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE laplacian_smooth
        USE class_psblas
        USE class_connectivity
        USE class_vector
        USE class_vertex
        USE tools_mesh_basics
        USE tools_math, ONLY : build_prec, solve_sys
        IMPLICIT NONE
        ! Local variables
        REAL(psb_dpk_), ALLOCATABLE :: bx(:),by(:),bz(:)! Right-hand side of the eqns
        REAL(psb_dpk_), ALLOCATABLE :: newx(:),newy(:),newz(:)  ! new locations
        REAL(psb_dpk_), ALLOCATABLE :: acoeffs(:)     ! non-zero coefficients of A
        INTEGER,        ALLOCATABLE :: ia(:),ja(:)    ! indexing arrays for A's elements
        TYPE(psb_dspmat_type)          :: Amatrix        ! the matrix for averaging for x,y,z
        INTEGER                        :: n_entries      ! number of non-zero elements of A
        INTEGER                        :: err_act, info, icontxt
        INTEGER                        :: i,j,iglob
        INTEGER                        :: nverts         ! number of local verts
        INTEGER,POINTER                :: iv2v(:)=>null()! given a vert, list neighbors
        REAL(psb_dpk_)                 :: rhs            ! an entry to be inserted in the RHS
        TYPE(psb_dprec_type)           :: Aprec          ! the preconditioner
        INTEGER                        :: iter
        REAL(psb_dpk_)                 :: err
        TYPE(vector)                   :: new_pos        ! new vertex position
        LOGICAL, ALLOCATABLE           :: fixed(:)       ! identifies verts to be fixed in pos.
        CHARACTER(len=32), PARAMETER   :: name_err = 'laplacian_smooth'

        ! Set error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Get PSBLAS context
        icontxt = icontxt_()

        nverts = psb_cd_get_local_rows(desc_v) ! local and shared vertices

        ALLOCATE(fixed(nverts), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        fixed = mixed .OR. on_boundary_(verts(1:nverts))
        !initialize so that all mixed-cell and boundary vertices are fixed

        IF (n_unconstrained > 0) THEN
            DO i = 1, n_unconstrained
                fixed( unconstrained(i) ) = .FALSE.
            ENDDO
        ENDIF

        !NEL is the number of diag. elements and NCONN the off-diag. elements
        n_entries = NEL_(v2v) + NCONN_(v2v) ! estimated non-zero values

        ! Allocate matrices
        CALL psb_spall(Amatrix, desc_v, info, n_entries)
        CALL psb_check_error(info,name_err,'psb_spall',icontxt)

        ! Allocate all three right-hand sides
        CALL psb_geall(bx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)

        CALL psb_geall(by, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)

        CALL psb_geall(bz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)

        ! Allocate all three solution vectors
        CALL psb_geall(newx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)

        CALL psb_geall(newy, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)

        CALL psb_geall(newz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geall',icontxt)


        ! Loop over all vertices and add equations to the linear systems
        ! Note the Amatrix is identical for all three systems
        DO i = 1, nverts

            iglob = loc_to_glob_(desc_v,i)

            IF ( fixed(i) ) THEN
                ! we simply apply a Dirichlet BC

                n_entries = 1

                ! insert 1.0 on diagonal
                CALL psb_spins(n_entries, (/iglob/), (/iglob/), (/1.0d0/), Amatrix, desc_v, info)
                CALL psb_check_error(info,name_err,'psb_spins',icontxt)

                ! insert vertex location on the right hand sides
                rhs = x_(verts(i))
                CALL psb_geins(1, (/iglob/), (/rhs/), bx, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

                rhs = y_(verts(i))
                CALL psb_geins(1, (/iglob/), (/rhs/), by, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

                rhs = z_(verts(i))
                CALL psb_geins(1, (/iglob/), (/rhs/), bz, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

            ELSE
                ! set this vertex's position to be the average of his neighbors
                ! thus the matrix has 1's on off-diagonals corresponding to neighbors
                ! and -1/N on the diagonal.  The RHS is zero
                ! -1/N*position_p + sum(positions_neighbors) = 0

                CALL get_ith_conn(iv2v, v2v,i) ! get neighbors
                n_entries = SIZE(iv2v) + 1 ! add 1 for the diagonal

                ALLOCATE(acoeffs(n_entries), ia(n_entries), ja(n_entries), stat=info)

                ! fill arrays with off-diag elements first
                DO j = 1, SIZE(iv2v)
                    acoeffs(j) = 1.0d0
                    ia(j)      = iglob
                    ja(j)      = loc_to_glob_(desc_v,iv2v(j))
                ENDDO

                acoeffs(n_entries) = -1.0d0*REAL(SIZE(iv2v))
                ia     (n_entries) = iglob
                ja     (n_entries) = iglob

                ! insert matrix coefficients
                CALL psb_spins(n_entries, ia, ja, acoeffs, Amatrix, desc_v, info)
                CALL psb_check_error(info,name_err,'psb_spins',icontxt)

                ! insert RHS
                CALL psb_geins(1, (/iglob/), (/0.0d0/), bx, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

                CALL psb_geins(1, (/iglob/), (/0.0d0/), by, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

                CALL psb_geins(1, (/iglob/), (/0.0d0/), bz, desc_v, info, psb_dupl_ovwrt_)
                CALL psb_check_error(info,name_err,'psb_geins',icontxt)

                DEALLOCATE(acoeffs, ia, ja)

            ENDIF ! end of if-test about interior vs. bndry point

        ENDDO

        ! Assemble matrices and right-hand sides

        CALL psb_spasb(Amatrix, desc_v, info, dupl= psb_dupl_ovwrt_)
        CALL psb_check_error(info,name_err,'psb_spasb',icontxt)

        CALL psb_geasb(bx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        CALL psb_geasb(by, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        CALL psb_geasb(bz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        ! Set initial guesses: Loop over all vertices
        DO i = 1, nverts

            iglob = loc_to_glob_(desc_v,i)

            CALL psb_geins(1, (/iglob/), (/x_(verts(i))/), newx, desc_v, info, psb_dupl_ovwrt_)
            CALL psb_check_error(info,name_err,'psb_geins',icontxt)

            CALL psb_geins(1, (/iglob/), (/y_(verts(i))/), newy, desc_v, info, psb_dupl_ovwrt_)
            CALL psb_check_error(info,name_err,'psb_geins',icontxt)

            CALL psb_geins(1, (/iglob/), (/z_(verts(i))/), newz, desc_v, info, psb_dupl_ovwrt_)
            CALL psb_check_error(info,name_err,'psb_geins',icontxt)

        ENDDO

        ! assemble guesses
        CALL psb_geasb(newx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        CALL psb_geasb(newy, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        CALL psb_geasb(newz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_geasb',icontxt)

        ! note the preconditioner can be re-used after each solve

        ! select and build preconditioner and solve
        CALL build_prec("BJAC",1,"BICGSTAB",Amatrix,desc_v,Aprec)
        CALL solve_sys(Amatrix, Aprec, bx, newx, desc_v, "BICGSTAB", 1.0D-6, 50, iter, err)
        CALL solve_sys(Amatrix, Aprec, by, newy, desc_v, "BICGSTAB", 1.0D-6, 50, iter, err)
        CALL solve_sys(Amatrix, Aprec, bz, newz, desc_v, "BICGSTAB", 1.0D-6, 50, iter, err)

        CALL psb_ovrl(newx,desc_v,info,update=psb_avg_)
        CALL psb_ovrl(newy,desc_v,info,update=psb_avg_)
        CALL psb_ovrl(newz,desc_v,info,update=psb_avg_)

        ! extract answers

        ! Loop over all vertices
        DO i = 1, nverts

            ! set position of vertices
            new_pos  = vector_( newx(i), newy(i), newz(i) )
            verts(i) = new_pos

        ENDDO

        ! free preconditioner
        CALL psb_precfree(Aprec,info)
        CALL psb_check_error(info,name_err,'psb_precfree',icontxt)

        ! free matrices and b vectors
        CALL psb_spfree(Amatrix, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_spfree',icontxt)

        CALL psb_gefree(bx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        CALL psb_gefree(by, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        CALL psb_gefree(bz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        CALL psb_gefree(newx, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        CALL psb_gefree(newy, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        CALL psb_gefree(newz, desc_v, info)
        CALL psb_check_error(info,name_err,'psb_gefree',icontxt)

        NULLIFY(iv2v)
        DEALLOCATE(fixed)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in LAPLACIAN_SMOOTH')

        END PROCEDURE laplacian_smooth

END SUBMODULE laplacian_smooth_implementation
