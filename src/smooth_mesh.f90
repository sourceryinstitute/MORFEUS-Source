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
! $Id: smooth_mesh.f90 3323 2008-08-28 15:44:18Z sfilippo $
!
! Description: sweeps through the mesh to find the optimum location of all the vertices
!
SUBMODULE (tools_mesh_optimize) smooth_mesh_implementation
    USE class_iterating, ONLY: iterating
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE smooth_mesh

        USE class_bc
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_keytable
        USE class_least_squares, ONLY : free_least_squares, set_least_squares
        USE class_mesh
        USE class_vector
        USE class_vertex
        USE tools_mesh_basics, ONLY : geom_cell, geom_diff, geom_face
        USE tools_mesh_check
        USE tools_mesh_optimize, nemo_protect_name => smooth_mesh

        IMPLICIT NONE

        ! Local variables

        !  Variables for list of unconstrained vertices
        !  note that unconstrained(:) is a list of integers, sized to ncells
        INTEGER, ALLOCATABLE :: unconstrained(:)  ! list of vertices that move
        ! (and connected only to tets)

        ! constrained lists a mobile, but constrained vertex
        INTEGER, ALLOCATABLE :: constrained(:)

        ! shared flags a vertex that is shared by another processor, sized to ncells
        INTEGER :: n_shared, tot_n_shared       ! number of shared vertices
        INTEGER, ALLOCATABLE :: shared(:)       !list of shared (overlap) vertices
        LOGICAL, ALLOCATABLE :: shared_flag(:)  !flags if a vertex is shared
        LOGICAL, ALLOCATABLE :: all_tets(:)     !flags if a vertex is used only by tets
        LOGICAL, ALLOCATABLE :: mixed(:)        !flags if a vertex is used only by tets

        REAL(psb_dpk_),ALLOCATABLE :: vpos(:,:) !  list of vertex positions

        INTEGER :: n_unconstrained     ! number of unconstrained vertices
        INTEGER :: n_constrained       ! number of constrained vertices

        TYPE(connectivity) ::   c2v ! given all vertices, finds connected cells
        TYPE(connectivity) ::   f2v ! given all vertices, finds connected faces
        TYPE(connectivity) ::   b2v ! given all verts, finds the boundaries to which they belong
        INTEGER,POINTER :: ib2v(:)=>null() ! given a vert, the bndry to which it belongs

        INTEGER :: iv,i             ! vertex number, cell id, and loop index
        INTEGER :: nverts           !number of local vertices in mesh
        INTEGER :: info
        INTEGER :: tangled     ! number of tangled cells in the initial mesh

        ! parameters for calling Opt-MS library
        INTEGER :: dims
        INTEGER :: technique, functionID
        TYPE(vector) :: new_pos  ! the new vertex position, after smoothing
        INTEGER,ALLOCATABLE :: idloc(:),idglob(:)

        ! for error handling
        INTEGER :: err_act,icontxt,mypnum
        CHARACTER(len=32), PARAMETER :: WHERE = 'smooth_mesh'

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()
        mypnum  = mypnum_()

        IF (msh%ncd == 2) RETURN   ! can't yet do 2d meshes

        nverts = psb_cd_get_local_rows(msh%desc_v) !local and shared vertices

        ! allocate storage for shared vertex positions
        CALL psb_geall(vpos,msh%desc_v,info,3)

        IF (info/=0) THEN
            WRITE(*,100)
            CALL abort_psblas
        ENDIF

        ! note shared (overlap) vertices
        CALL psb_get_overlap(shared,msh%desc_v,info)

        IF (info/=0) THEN
            WRITE(*,200)
            CALL abort_psblas
        ENDIF

        ALLOCATE(shared_flag(nverts), all_tets(nverts), mixed(nverts),stat=info)

        IF (info/=0) THEN
            WRITE(*,600)
            CALL abort_psblas
        ENDIF

        shared_flag = .FALSE.

        IF (ALLOCATED(shared)) THEN
            n_shared=SIZE(shared)
            DO i = 1, n_shared
                shared_flag( shared(i) ) = .TRUE.
            END DO
        ELSE
            n_shared=0  ! if there are no shared vertices, shared will still have size 1
        ENDIF
        tot_n_shared = n_shared
        CALL psb_sum(icontxt,tot_n_shared)

        IF (info/=0) THEN
            WRITE(*,200)
            CALL abort_psblas
        ENDIF

        !initialize OptMS parameters
        dims       = 3
        technique  = -1  ! -1 is default, 4 is OPTMS_COMBINED:  see SMuserDefs.h
        functionID = 27 ! Minimize Max SMRS Volume Ratio
        !(SMRS Vol. Ratio is Knupp's metric to the -3/2 power)

        ! pass parameters and initialize OptMS data structures for 3D smoothing
        info = initoptms(dims,technique,functionID)

        IF (info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        ENDIF

        !initialize OptMS parameters for 2D smoothing
        dims       = 2
        technique  = 3  ! -1 is default, 4 is OPTMS_COMBINED:  see SMuserDefs.h
        functionID = 7 ! was 7 ! Minimize Max Jacobian

        ! pass parameters and initialize OptMS data structures
        info = initoptms2d(dims,technique,functionID)

        ! prepare the c2v connectivity, required for calculating unconstrained vertices
        CALL msh%v2c%get_dual_conn(c2v)

        ! prepare the f2v connectivity, required for calculating surf. vtx movement
        CALL msh%v2f%get_dual_conn(f2v)

        ! prepare the b2v connectivity, required for calculating surf. vtx movement
        CALL msh%v2b%get_dual_conn(b2v)

        ALLOCATE(unconstrained(nverts),constrained(nverts),stat=info)
        IF(info /= 0) THEN
            WRITE(*,600)
            CALL abort_psblas
        END IF

        ! decide which vertices we will smooth and which are special boundary vertices
        CALL mobile_verts(msh,bc,c2v,shared_flag,unconstrained,n_unconstrained,constrained,  &
            &    n_constrained, all_tets, mixed)

        CALL check_right_handed(msh,shared,shared_flag,tangled,all_tets)

        IF (tangled > 0) THEN
            WRITE(6,*)
            WRITE(6,'(a,i4,a,i3,a)')"Warning: ",tangled, &
                & " inverted cells detected on processor ",mypnum_(),"."
            WRITE(6,'(a)')"Attempting to untangle mesh."
        ENDIF

        !           ============= repeatedly do surface sweeps ======================
        CALL surface_iter%reset()

        surface_iteration: DO

            CALL surface_iter%increment()

            IF (mypnum_() == 0) WRITE(6,'(a,i4,2a)')" - Iteration: ", &
                & surface_iter%current_iteration(), CHAR(27),CHAR(77)

            !  Loop over interior vertices
            DO i=1,n_constrained

                iv = constrained(i) ! get index of unconstrained vertex out of the list

                CALL b2v%get_ith_conn(ib2v,iv)

                IF ( SIZE(ib2v) == 1 ) THEN
                    ! a normal sliding vertex can only belong to one boundary,
                    ! so assume ib2v has only 1 elem.
                    IF (.NOT. mixed(iv) ) &
                        & CALL smooth_surf_vtx(iv,ib2v(1),msh,f2v,shared_flag, tangled, all_tets)

                ENDIF
                ! TBD: we could do some sort of averaging with points adjoining 2 sliding BC's

            END DO

            !copy all shared vertices to a dense vector and average among overlapping processors
            IF (tot_n_shared > 0) THEN ! there are shared vertices
                DO i = 1,n_shared

                    iv = shared(i)
                    vpos(iv,1) = msh%verts(iv)%x_()
                    vpos(iv,2) = msh%verts(iv)%y_()
                    vpos(iv,3) = msh%verts(iv)%z_()

                END DO ! end of loop over shared vertices

                CALL psb_ovrl(vpos,msh%desc_v,info,update=psb_avg_)

                IF(info /= 0) THEN
                    WRITE(*,700)
                    CALL abort_psblas
                END IF

                !now copy the averaged positions from vpos back to the vertices
                DO i = 1,n_shared
                    iv = shared(i)
                    new_pos = vector_(vpos(iv,1),vpos(iv,2),vpos(iv,3))
                    msh%verts(iv) = new_pos
                END DO

                CALL update_vertex_halo(msh%verts,msh%desc_v)

            ENDIF ! end of if-test for allocation of shared vertex storage

            tangled = 0 ! assume that after one surface sweep, things are untangled

            IF ( surface_iter%stop_iterating() ) EXIT surface_iteration
        END DO surface_iteration ! end of surface sweep loops

        ! Now that the surface is hopefully in good shape, do interior smoothing
        IF (mypnum_() == 0) WRITE(6,*)

        ! Free memory, because these will be recalculated later.  For now, since
        ! they are immediately out of date, they are worse than useless
        DEALLOCATE(msh%face_cntr, msh%af, msh%area,msh%cell_cntr,msh%vol,msh%df, &
            & msh%dist,msh%interp)

        CALL free_least_squares(msh%lsr)

        ! set all interior points (except at region boundaries) to the avg. position
        CALL laplacian_smooth (msh%desc_v, msh%v2v, n_unconstrained, unconstrained, msh%verts, mixed)

        CALL interior_iter%reset()

        IF (mypnum_() == 0) WRITE(6,'(a,i4,2a)')" - Iteration: ",  &
            &  interior_iter%current_iteration(),CHAR(27),CHAR(77)

        interior_iteration: DO
            !  Loop over interior vertices
            DO i=1,n_unconstrained

                iv = unconstrained(i) ! get index of unconstrained vertex out of the list

                ! Optimize the location of the interior vertex
                IF (.NOT. mixed(iv) ) &
                    & CALL smooth_interior_vtx(iv,msh,c2v,shared_flag, all_tets)

            END DO ! end of loop over unconstrained vertices
            !copy all shared vertices to a dense vector and average among overlapping processors
            IF (tot_n_shared > 0) THEN ! there are shared vertices
                DO i = 1,n_shared
                    iv = shared(i)
                    vpos(iv,1) = msh%verts(iv)%x_()
                    vpos(iv,2) = msh%verts(iv)%y_()
                    vpos(iv,3) = msh%verts(iv)%z_()
                END DO
                CALL psb_ovrl(vpos,msh%desc_v,info,update=psb_avg_)

                IF(info /= 0) THEN
                    WRITE(*,700)
                    CALL abort_psblas
                END IF

                !now copy the averaged positions from vpos back to the vertices
                DO i = 1,n_shared
                    iv = shared(i)
                    new_pos = vector_(vpos(iv,1),vpos(iv,2),vpos(iv,3))
                    msh%verts(iv) = new_pos
                END DO

                CALL update_vertex_halo(msh%verts,msh%desc_v)

            ENDIF ! end of if-test for allocation of shared vertex storage

            CALL check_right_handed(msh,shared,shared_flag,tangled,all_tets)

            !share largest value of tangled among processors

            CALL psb_amx(icontxt,tangled)
            CALL psb_check_error(info,TRIM(WHERE),'psb_amx',icontxt)

            IF (tangled > 0) THEN
                WRITE(6,*)
                WRITE(6,'(a,i4,a,i3,a)')"Warning: ",tangled, &
                    & " inverted cells detected on processor ",mypnum_(),"."
                WRITE(6,'(a)')"Attempting to untangle mesh."
                WRITE(6,*)

            ENDIF

            CALL interior_iter%increment()

            IF (mypnum_() == 0) WRITE(6,'(a,i4,2a)')" - Iteration: ",  &
                &  interior_iter%current_iteration(),CHAR(27),CHAR(77)


            IF ( interior_iter%stop_iterating() ) EXIT interior_iteration
        END DO interior_iteration

        IF (mypnum_() == 0) THEN
            WRITE(6,*)
            WRITE(6,*)
        ENDIF  !End of if-test for processor zero

        !now recalculate all mesh geometry

        ! Computes face-related metrics members MSH%FACE_CNTR, MSH%AF, MSH%AREA
        CALL geom_face(msh%verts,msh%v2f,msh%ncd, &
            & msh%face_cntr,msh%af,msh%area)

        ! Computes cell-related metrics members MSH%CELL_CNTR, MSH%VOL
        CALL geom_cell(msh%verts,msh%faces,msh%cells,msh%v2f,msh%v2c, &
            & msh%f2c,msh%ncd,msh%cell_cntr,msh%vol)

        ! Computes face-related metrics members MSH%DF, MSH%DIST, MSH%INTERP
        CALL geom_diff(msh%faces,msh%f2b,msh%face_cntr,msh%af,msh%cell_cntr, &
            & msh%df,msh%dist,msh%interp)

        ! Computes metrics for cell-centered least squares regression
        CALL set_least_squares(msh%lsr,msh%ncd,msh%desc_c,msh%c2c,msh%f2b, &
            & msh%faces,msh%cell_cntr,msh%face_cntr)

        ! free allocated memory

        CALL psb_gefree(vpos,msh%desc_v,info)

        IF ( ALLOCATED(shared) ) THEN
            DEALLOCATE (shared)
        ENDIF

        DEALLOCATE (shared_flag, all_tets)

        info=freeoptms()
        info=freeoptms2d()

        CALL free_conn(c2v)
        CALL free_conn(f2v)
        CALL free_conn(b2v)

        DEALLOCATE(unconstrained, constrained)

        IF ( ALLOCATED(idloc) )  DEALLOCATE(idloc)

        IF ( ALLOCATED(idglob) ) DEALLOCATE(idglob)

        IF(info /= 0) THEN
            WRITE(*,600)
            CALL abort_psblas
        END IF

        IF(info /= 0) THEN
            WRITE(*,600)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Failure to allocate dense matrix in SMOOTH_MESH')
200     FORMAT(' ERROR! Failed to get overlap points in SMOOTH_MESH')
600     FORMAT(' ERROR! Memory allocation failure in SMOOTH_MESH')
700     FORMAT(' ERROR! Failure to average shared vertex positions in SMOOTH_MESH')

        END PROCEDURE smooth_mesh

END SUBMODULE smooth_mesh_implementation
