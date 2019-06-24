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
! $Id: class_vertex.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    to be added...
!
SUBMODULE(class_vertex) class_vertex_procedures

    USE class_psblas
    USE class_vector

    IMPLICIT NONE

CONTAINS


    MODULE PROCEDURE nemo_vertex_sizeof
    IMPLICIT NONE

        nemo_vertex_sizeof = nemo_sizeof_int + vtx%position%nemo_sizeof()

    END PROCEDURE nemo_vertex_sizeof
    ! ----- Constructors -----

    ! Public default constructor
    MODULE PROCEDURE vertex_1_

        res%position = vector_(x,y,z)

        IF(PRESENT(on_boundary)) THEN
            res%on_boundary = on_boundary
        ELSE
            res%on_boundary = .FALSE.
        END IF

    END PROCEDURE vertex_1_

    MODULE PROCEDURE vertex_2_

        res%position = position

        IF(PRESENT(on_boundary)) THEN
            res%on_boundary = on_boundary
        ELSE
            res%on_boundary = .FALSE.
        END IF

    END PROCEDURE vertex_2_


    ! Array constructor
    MODULE PROCEDURE alloc_vertex
        !
        INTEGER :: i, info

        IF(ALLOCATED(verts))THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ALLOCATE(verts(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ! Initizialization
        DO i = 1, n
            verts(i)%position = vector_(0.d0,0.d0,0.d0)
            verts(i)%on_boundary = .FALSE.
        END DO

100     FORMAT(' ERROR! Illegal allocation: VERTEX array already associated')
200     FORMAT(' ERROR! Memory allocation failure in ALLOC_VERTEX')

    END PROCEDURE alloc_vertex


    ! ----- Destructor -----

    MODULE PROCEDURE free_vertex
        !
        INTEGER :: info

        IF (ALLOCATED(verts)) THEN
            DEALLOCATE(verts,stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF
100     FORMAT(' ERROR! Memory deallocation failure in FREE_VERTEX')

    END PROCEDURE free_vertex


    ! ----- Parallel Operations -----

    MODULE PROCEDURE bcast_vertex
        !
        INTEGER :: icontxt, mypnum
        INTEGER :: i, info, n
        LOGICAL, ALLOCATABLE :: lbuf(:)
        TYPE(vector), ALLOCATABLE :: vbuf(:)

        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! VERTS(:) is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            n = SIZE(verts)
            CALL psb_bcast(icontxt,n)

            IF(.NOT.ALLOCATED(verts)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        ELSE
            CALL psb_bcast(icontxt,n)

            IF(ALLOCATED(verts) .AND. SIZE(verts) /= n) THEN
                CALL free_vertex(verts)
            END IF

            IF(.NOT.ALLOCATED(verts)) THEN
                CALL alloc_vertex(verts,n)
            END IF
        END IF

        ALLOCATE(lbuf(n),vbuf(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        IF(mypnum == 0) THEN
            DO i = 1, n
                vbuf(i) = verts(i)%position
                lbuf(i) = verts(i)%on_boundary
            END DO

            CALL bcast_vector(vbuf)
            CALL psb_bcast(icontxt,lbuf)
        ELSE
            CALL bcast_vector(vbuf)
            CALL psb_bcast(icontxt,lbuf)

            DO i = 1, n
                verts(i)%position = vbuf(i)
                verts(i)%on_boundary = lbuf(i)
            END DO
        END IF

        DEALLOCATE(vbuf,lbuf)

100     FORMAT(' ERROR! Illegal Send: VERTEX array not associated')
200     FORMAT(' ERROR! Memory allocation failure in BCAST_VERTEX')

    END PROCEDURE bcast_vertex


    MODULE PROCEDURE g2l_vertex
        USE psb_base_mod
        !
        INTEGER :: info
        INTEGER :: iv_glob, iv_loc, n_loc
        INTEGER, ALLOCATABLE :: iloc_to_glob(:)
        LOGICAL, ALLOCATABLE :: lbuf(:)
        TYPE(vector), ALLOCATABLE :: vbuf(:)


        n_loc = psb_cd_get_local_cols(desc_v)

        ALLOCATE(lbuf(n_loc),vbuf(n_loc),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL psb_get_loc_to_glob(desc_v,iloc_to_glob)

        DO iv_loc = 1, n_loc
            iv_glob = iloc_to_glob(iv_loc)
            vbuf(iv_loc) = verts(iv_glob)%position
            lbuf(iv_loc) = verts(iv_glob)%on_boundary
        END DO

        ! Reallocates VERTS
        CALL free_vertex(verts)
        CALL alloc_vertex(verts,n_loc)

        verts(1:n_loc)%position = vbuf(1:n_loc)
        verts(1:n_loc)%on_boundary = lbuf(1:n_loc)

        DEALLOCATE(iloc_to_glob)
        DEALLOCATE(vbuf,lbuf)

100     FORMAT(' ERROR! Memory allocation failure in G2L_VERTEX')

    END PROCEDURE g2l_vertex


    MODULE PROCEDURE l2g_vertex
        USE psb_base_mod
        ! WARNING! The global results is allocated only on P0. After its usage
        ! it must be deallocated in the calling unit by means of the statement:
        ! "if(associated(glob_res)) deallocate(glob_res)"

        INTEGER :: err_act, info, icontxt
        INTEGER :: iv, n_glob, n_loc
        TYPE(vector), ALLOCATABLE :: vbuf_glob(:)
        TYPE(vector), ALLOCATABLE :: vbuf_loc(:)


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        n_loc = SIZE(verts_loc)             ! # local vertices
        n_glob = psb_cd_get_global_cols(desc_v) ! # global vertices

        ! Allocation of local objects
        CALL alloc_vector(vbuf_loc,n_loc)

        vbuf_loc(:) = verts_loc(:)%position
        CALL l2g_vector(vbuf_loc,vbuf_glob,desc_v)


        IF(mypnum_() == 0) THEN
            ! Allocation of global objects
            IF(ALLOCATED(verts_glob)) CALL free_vertex(verts_glob)
            CALL alloc_vertex(verts_glob,n_glob)

            DO iv = 1, n_glob
                verts_glob(iv)%position = vbuf_glob(iv)
            END DO

        END IF

        ! Frees memory
        IF(ALLOCATED(vbuf_glob)) CALL free_vector(vbuf_glob)
        CALL free_vector(vbuf_loc)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE l2g_vertex


    MODULE PROCEDURE update_vertex_halo
        ! synchronizes vertex positions
        ! (other data could be synchronized, too, if needed)

        CALL update_vector_halo(verts%position,desc)

    END PROCEDURE update_vertex_halo


    ! ----- Getters -----

    MODULE PROCEDURE position_

        position_ = vert%position

    END PROCEDURE position_


    MODULE PROCEDURE on_boundary_

        on_boundary_ = vert%on_boundary

    END PROCEDURE on_boundary_


    MODULE PROCEDURE get_vertex_x

        get_vertex_x = vert%position%x_()

    END PROCEDURE get_vertex_x


    MODULE PROCEDURE get_vertex_y

        get_vertex_y = vert%position%y_()

    END PROCEDURE get_vertex_y


    MODULE PROCEDURE get_vertex_z

        get_vertex_z = vert%position%z_()

    END PROCEDURE get_vertex_z


    ! ----- Setters -----

    MODULE PROCEDURE set_vertex_position

        vert%position = position

    END PROCEDURE set_vertex_position


    ! ----- Vector Algebra Operations -----

    MODULE PROCEDURE vert_sum_1

        vert_sum_1 = a%position + b%position

    END PROCEDURE vert_sum_1

    MODULE PROCEDURE vert_sum_2

        vert_sum_2 = a + b%position

    END PROCEDURE vert_sum_2


    MODULE PROCEDURE vert_diff

        vert_diff = a%position - b%position

    END PROCEDURE vert_diff


    MODULE PROCEDURE scalar_vertex_prod

        scalar_vertex_prod%position = alpha * v%position
        scalar_vertex_prod%on_boundary = v%on_boundary

    END PROCEDURE scalar_vertex_prod


    MODULE PROCEDURE dot_prod

        dot_prod = a%position .dot. b%position

    END PROCEDURE dot_prod


    MODULE PROCEDURE cross_prod

        cross_prod = a%position .cross. b%position

    END PROCEDURE cross_prod


    MODULE PROCEDURE vert_mag

        vert_mag = v%position%mag()

    END PROCEDURE vert_mag


END SUBMODULE class_vertex_procedures
