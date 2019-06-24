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
SUBMODULE (tools_mesh) supplement_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE supplement_v2f 
        USE class_psblas
        USE class_connectivity
        USE class_face
        USE class_keytable
        IMPLICIT NONE
        !! $Id: supplement_v2f.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:  Creates supplementary connectivity info about cells connected to
        !!               overlap (shared) vertices.  The data are saved as keytables using
        !!               global indexing.  Since this routine is called prior to re-
        !!               allocating the mesh data, all indices (except the call to get
        !!               overlap vertex id's) use global ID numbers.
        !!
        ! Local variables
        INTEGER, PARAMETER :: max_tri = 15
        INTEGER :: i, j, info, err_act
        INTEGER :: n_shared, n_faces
        INTEGER :: lowest, highest
        INTEGER :: iface, iv, num_incident_tri
        INTEGER, ALLOCATABLE :: glob_v_ovrl(:)
        INTEGER, ALLOCATABLE :: v_ovrl(:)        ! overlap vertices
        INTEGER, POINTER :: if2v(:)  => NULL()   ! given a vert, finds connected faces
        INTEGER, POINTER :: iv2f(:)  => NULL()   ! given a face, finds connected verts
        INTEGER          :: ibf2v(1:max_tri)     ! given a face, finds connected bndry verts
        TYPE(connectivity) :: f2v
        TYPE(keytable) :: f2vtable, v2ftable


        ! sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        CALL bcast_conn(v2f)

        ! get overlap vertices--returns local numbering
        CALL psb_get_overlap(v_ovrl, desc_v, info)
        CALL psb_check_error(info,'supplement_v2f','psb_get_overlap',icontxt_())

        ! if there are no overlap vertices, then there is nothing for us to do!
        IF (.NOT. ALLOCATED(v_ovrl)) THEN
            IF(nprocs_() > 1 .AND. mypnum_() == 0) THEN
                WRITE(*,300)
                WRITE(*,200)
                WRITE(*,300)
                WRITE(*,*)
            END IF
            RETURN
        END IF

        IF ( mypnum_() == 0 ) THEN
            WRITE(*,*) 'Creating supplemental overlap vertex information (cont.)'
            WRITE(*,*)
        END IF

        CALL get_dual_conn(v2f,f2v)

        n_shared = SIZE (v_ovrl)

        ALLOCATE( glob_v_ovrl(n_shared), stat = info)
        IF (info/=0) THEN
            WRITE(*,100)
            CALL abort_psblas
        ENDIF

        ! make a non-pointer copy of local numbers to send to psblas routine
        glob_v_ovrl = v_ovrl

        ! we need global numbering: overwrites local numbers in glob_v_ovrl
        CALL psb_loc_to_glob(glob_v_ovrl,desc_v,info)
        CALL psb_check_error(info,'supplement_v2f','psb_loc_to_glob',icontxt_())

        ! note range for v2f_suppl
        lowest = MINVAL(glob_v_ovrl)
        highest = MAXVAL(glob_v_ovrl)

        ! instantiate keytable for f2v conversion using global ID's
        CALL alloc_keytable(f2vtable,lowest,highest)

        lowest  =  HUGE(lowest)
        highest = -HUGE(highest)

        ! loop over vertices and tally memory requirements for v2f_suppl
        DO i = 1, n_shared

            ! get global vertex number
            iv = glob_v_ovrl(i)

            ! get the faces connected to this vertex
            CALL get_ith_conn(if2v,f2v,iv)

            ! cull faces that do not lie on boundaries (interior faces)

            num_incident_tri = 0
            DO j = 1, SIZE(if2v)
                iface = if2v(j)

                IF ( flag_(faces(iface)) > 0 ) THEN
                    num_incident_tri = num_incident_tri + 1
                    ibf2v(num_incident_tri) = iface
                ENDIF

            ENDDO

            ! copy to pointer
            ALLOCATE(if2v(num_incident_tri))
            if2v = ibf2v(1:num_incident_tri)

            ! store in a keytable
            CALL set_kt_row(f2vtable, iv, if2v)

            lowest = MIN(lowest,MINVAL(if2v))
            highest = MAX(highest,MAXVAL(if2v))
            DEALLOCATE(if2v)
        ENDDO

        ! now make v2f keytable, again, based on global IDs
        CALL alloc_keytable(v2ftable,lowest,highest)

        ! loop over vertices and tally memory requirements for v2f_suppl
        DO i = 1, n_shared

            ! find what faces are connected to these vertices
            iv = glob_v_ovrl(i)

            ! get the faces from the keytable, since we need a list that excludes
            ! interior faces
            CALL get_kt_row(f2vtable, iv, if2v)
            n_faces = SIZE(if2v)

            DO j = 1 , n_faces
                iface = if2v(j)
                CALL get_ith_conn(iv2f, v2f, iface)
                CALL set_kt_row(v2ftable, iface, iv2f)
            ENDDO
        ENDDO

        NULLIFY(iv2f,if2v)
        DEALLOCATE(v_ovrl,glob_v_ovrl)
        CALL free_conn(f2v)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in SUPPLEMENT_V2F')
200     FORMAT(' * DEBUG: empty overlap set with NPROCS > 1 *')
300     FORMAT(' ********************************************')

        END PROCEDURE supplement_v2f


        MODULE PROCEDURE supplement_v2c 
        USE class_psblas
        USE class_connectivity
        USE class_keytable
        IMPLICIT NONE
        !! $Id: supplement_v2c.f90 3093 2008-04-22 14:51:09Z sfilippo $
        !!
        !! Description:  Creates supplementary connectivity info about cells connected to
        !!               overlap (shared) vertices.  The data are saved as keytables using
        !!               global indexing
        !!

        ! Local variables
        INTEGER :: i, j, info, err_act
        INTEGER :: n_shared, n_cells
        INTEGER :: lowest, highest
        INTEGER :: ic, iv
        INTEGER, ALLOCATABLE :: glob_v_ovrl(:)
        INTEGER, ALLOCATABLE :: v_ovrl(:) ! overlap vertices
        INTEGER, POINTER :: ic2v(:) => NULL()   ! given a vert, finds connected cells
        INTEGER, POINTER :: iv2c(:) => NULL()   ! given a cell, finds connected verts
        TYPE(connectivity) :: c2v
        TYPE(keytable) :: c2vtable, v2ctable


        ! sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        ! get overlap vertices--returns local numbering
        CALL psb_get_overlap(v_ovrl, desc_v, info)
        CALL psb_check_error(info,'supplement_v2c','psb_get_overlap',icontxt_())

        ! if there are no overlap vertices, then there is nothing for us to do!
        IF (.NOT. ALLOCATED(v_ovrl)) THEN
            IF(nprocs_() > 1 .AND. mypnum_() == 0) THEN
                WRITE(*,300)
                WRITE(*,200)
                WRITE(*,300)
                WRITE(*,*)
            END IF
            RETURN
        END IF

        IF ( mypnum_() == 0 ) THEN
            WRITE(*,*) 'Creating supplemental overlap vertex information'
        END IF

        CALL get_dual_conn(v2c,c2v)

        n_shared = SIZE (v_ovrl)

        ALLOCATE( glob_v_ovrl(n_shared), stat = info)
        IF (info/=0) THEN
            WRITE(*,100)
            CALL abort_psblas
        ENDIF

        ! make a non-pointer copy to send to psblas routine
        glob_v_ovrl = v_ovrl

        ! we need global numbering: overwrites local numbers in glob_v_ovrl
        CALL psb_loc_to_glob(glob_v_ovrl,desc_v,info)
        CALL psb_check_error(info,'supplement_v2c','psb_loc_to_glob',icontxt_())

        ! note range for v2c_suppl
        lowest = MINVAL(glob_v_ovrl)
        highest = MAXVAL(glob_v_ovrl)

        ! instantiate keytable
        CALL alloc_keytable(c2vtable,lowest,highest)

        lowest  =  HUGE(lowest)
        highest = -HUGE(highest)

        ! loop over vertices and tally memory requirements for v2c_suppl
        DO i = 1, n_shared

            ! get global vertex number
            iv = glob_v_ovrl(i)

            ! get the cells connected to this vertex
            CALL get_ith_conn(ic2v,c2v,iv)

            ! store in a key table
            CALL set_kt_row(c2vtable, iv, ic2v)

            lowest = MIN(lowest,MINVAL(ic2v))
            highest = MAX(highest,MAXVAL(ic2v))
        ENDDO

        ! now make v2c keytable
        CALL alloc_keytable(v2ctable,lowest,highest)

        ! loop over vertices and tally memory requirements for v2c_suppl
        DO i = 1, n_shared

            ! find what cells are connected to these vertices
            iv = glob_v_ovrl(i)

            CALL get_ith_conn(ic2v,c2v,iv)
            n_cells = SIZE(ic2v)

            DO j = 1 , n_cells
                ic = ic2v(j)
                CALL get_ith_conn(iv2c, v2c, ic)
                CALL set_kt_row(v2ctable, ic, iv2c)
            ENDDO
        ENDDO

        NULLIFY(iv2c,ic2v)
        DEALLOCATE(v_ovrl,glob_v_ovrl)
        CALL free_conn(c2v)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in SUPPLEMENT_V2C')
200     FORMAT(' * DEBUG: empty overlap set with NPROCS > 1 *')
300     FORMAT(' ********************************************')

        END PROCEDURE supplement_v2c

END SUBMODULE supplement_implementation
