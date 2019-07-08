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
! $Id: class_face.F90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    face class and face management
!
! Provides:
!    FACE            class      (no public data)
!    FACE_           constructor (elemental)
!    FREE_FACE       destructor
!    NV_             get number of cell faces
!    MASTER_         get the master cell for a face
!    SLAVE_          get the slave cell for a face
!    FLAG_           indicates of a face is interior, on a processor/physical bndry
!    BCAST_FACE      broadcast face structures from process 0 to other processes
!    G2L_FACE        global to local renumbering and reallocation of faces
!    SET_FACE        resets number of vertices, master, slave, or flag for a face
!
SUBMODULE(class_face) class_face_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_face_sizeof
        USE psb_base_mod
        IMPLICIT NONE

        nemo_face_sizeof = 4 * nemo_sizeof_int

    END PROCEDURE nemo_face_sizeof


    ! ----- Constructors -----

    MODULE PROCEDURE face_
        IMPLICIT NONE
        !! Public default constructor

        !face_ = face(nv,master,slave,flag)
        !! Workaround for Intel 18 error #6053: Structure constructor may not have fields with a PRIVATE attribute
        face_%nv     = nv
        face_%master = master
        face_%slave  = slave
        face_%flag   = flag
    END PROCEDURE face_


    MODULE PROCEDURE alloc_face
        IMPLICIT NONE
        !! Array constructor

        INTEGER :: info

        IF(ALLOCATED(faces))THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ALLOCATE(faces(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Illegal allocation: FACE array already associated')
200     FORMAT(' ERROR! Memory allocation failure in ALLOC_FACE')

    END PROCEDURE alloc_face


    ! ----- Destructor -----

    MODULE PROCEDURE free_face
        IMPLICIT NONE
        !
        INTEGER :: info

        IF (ALLOCATED(faces)) THEN
            DEALLOCATE(faces,stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_FACE')

    END PROCEDURE free_face


    ! ----- Parallel Ops. -----

    MODULE PROCEDURE bcast_face
        IMPLICIT NONE
        !
        INTEGER, ALLOCATABLE :: intbuf(:,:)
        INTEGER :: i, info, n
        INTEGER :: icontxt, mypnum

        mypnum = mypnum_()
        icontxt = icontxt_()

        ! FACES(:) is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            IF(.NOT.ALLOCATED(faces)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            n = SIZE(faces)

            CALL psb_bcast(icontxt,n)
        ELSE
            IF(ALLOCATED(faces)) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            CALL psb_bcast(icontxt,n)

            CALL alloc_face(faces,n)
        END IF

        ALLOCATE(intbuf(4,n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ! Broadcasting
        IF(mypnum == 0) THEN
            DO i = 1, n
                intbuf(1,i) = faces(i)%nv
                intbuf(2,i) = faces(i)%master
                intbuf(3,i) = faces(i)%slave
                intbuf(4,i) = faces(i)%flag
            END DO

            CALL psb_bcast(icontxt,intbuf)
        ELSE
            CALL psb_bcast(icontxt,intbuf)

            DO i = 1, n
                faces(i)%nv = intbuf(1,i)
                faces(i)%master = intbuf(2,i)
                faces(i)%slave  = intbuf(3,i)
                faces(i)%flag   = intbuf(4,i)
            END DO
        END IF

        DEALLOCATE(intbuf)

100     FORMAT(' ERROR! Illegal Send: FACE array not associated')
200     FORMAT(' ERROR! Illegal Recv: FACE array already associated')
300     FORMAT(' ERROR! Memory allocation failure in BCAST_FACE')

    END PROCEDURE bcast_face


    MODULE PROCEDURE g2l_face
        USE psb_base_mod
        IMPLICIT NONE
        !!
        INTEGER :: flag, if_glob, if_loc
        INTEGER :: im_glob, im_loc, is_glob, is_loc
        INTEGER :: nfaces_loc
        INTEGER, ALLOCATABLE :: ic_glob_to_loc(:)
        INTEGER, ALLOCATABLE :: if_loc_to_glob(:)
        TYPE(face), ALLOCATABLE :: work(:)

        nfaces_loc  = psb_cd_get_local_cols(desc_f)

        CALL alloc_face(work,nfaces_loc)

        CALL psb_get_glob_to_loc(desc_c,ic_glob_to_loc)

        CALL psb_get_loc_to_glob(desc_f,if_loc_to_glob)

        DO if_loc = 1, nfaces_loc
            if_glob = if_loc_to_glob(if_loc)

            im_glob = faces(if_glob)%master
            is_glob = faces(if_glob)%slave
            flag    = faces(if_glob)%flag

            work(if_loc) = faces(if_glob)

            IF(flag == 0) THEN
                im_loc = ic_glob_to_loc(im_glob)
                is_loc = ic_glob_to_loc(is_glob)
                work(if_loc)%master = im_loc
                work(if_loc)%slave  = is_loc

                ! Faces having master/slave lying on another process, further
                ! the halo layer (ic_glob_to_loc < 0)
                IF(im_loc < 0 .OR. is_loc < 0) THEN
                    work(if_loc)%flag = -1
                END IF
            ELSE
                im_loc = ic_glob_to_loc(im_glob)
                work(if_loc)%master = im_loc
                work(if_loc)%slave  = 0
            END IF
        END DO
#if defined(HAVE_MOVE_ALLOC)
        CALL move_ALLOC(work,faces)
#else
        IF (SIZE(faces) /= SIZE(work)) THEN
            DEALLOCATE(faces)
            ALLOCATE(faces(SIZE(work)))
        ENDIF
        faces = work
        DEALLOCATE(work)
#endif

        ! Frees memory storage
        DEALLOCATE(if_loc_to_glob)
        DEALLOCATE(ic_glob_to_loc)

    END PROCEDURE g2l_face


    MODULE PROCEDURE l2g_face
        USE psb_base_mod
        IMPLICIT NONE
        !! WARNING! The global results is allocated only on P0. After its usage
        !! it must be deallocated in the calling unit by means of the statement:
        !! "if(associated(glob_res)) deallocate(glob_res)"
        INTEGER :: icontxt, info, err_act
        INTEGER :: IF, im_loc, is_loc, n_glob, n_loc
        INTEGER, ALLOCATABLE :: ibuf_glob(:,:)
        INTEGER, ALLOCATABLE :: ibuf_loc(:,:)
        INTEGER, ALLOCATABLE :: iloc_to_glob(:)

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        n_loc  = SIZE(faces_loc)            ! # local faces
        n_glob = psb_cd_get_global_cols(desc_f) ! # global faces

        ! Allocation of local vectors
        CALL psb_geall(ibuf_loc,desc_f,info,4)
        CALL psb_check_error(info,'l2g_face','psb_geall',icontxt)

        ALLOCATE(ibuf_glob(n_glob,4),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL psb_get_loc_to_glob(desc_c,iloc_to_glob)

        DO IF = 1, n_loc
            im_loc = faces_loc(IF)%master
            is_loc = faces_loc(IF)%slave

            IF(faces_loc(IF)%flag == 0) THEN
                ibuf_loc(IF,1) = iloc_to_glob(im_loc)
                ibuf_loc(IF,2) = iloc_to_glob(is_loc)
            ELSEIF(faces_loc(IF)%flag > 0) THEN
                ibuf_loc(IF,1) = iloc_to_glob(im_loc)
                ibuf_loc(IF,2)  = 0
            END IF

            ibuf_loc(IF,3) = faces_loc(IF)%nv
            ibuf_loc(IF,4) = faces_loc(IF)%flag
        END DO


        CALL psb_gather(ibuf_glob,ibuf_loc,desc_f,info,root=0)
        CALL psb_check_error(info,'l2g_face','psb_igatherm',icontxt)


        IF(mypnum_() == 0) THEN
            IF(ALLOCATED(faces_glob)) CALL free_face(faces_glob)
            CALL alloc_face(faces_glob,n_glob)
            DO IF = 1, n_glob
                faces_glob(IF)%master = ibuf_glob(IF,1)
                faces_glob(IF)%slave  = ibuf_glob(IF,2)
                faces_glob(IF)%nv     = ibuf_glob(IF,3)
                faces_glob(IF)%flag   = ibuf_glob(IF,4)
            END DO
        END IF

        ! Frees memory
        DEALLOCATE(iloc_to_glob)
        DEALLOCATE(ibuf_glob)

        CALL psb_gefree(ibuf_loc,desc_f,info)
        CALL psb_check_error(info,'l2g_face','psb_gefree',icontxt)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in L2G_FACE')

    END PROCEDURE l2g_face


    ! ----- Getters -----

    MODULE PROCEDURE get_face_nv
        IMPLICIT NONE

        get_face_nv = f%nv

    END PROCEDURE get_face_nv


    MODULE PROCEDURE get_face_master
        IMPLICIT NONE

        get_face_master = f%master

    END PROCEDURE get_face_master


    MODULE PROCEDURE get_face_slave
        IMPLICIT NONE

        get_face_slave = f%slave

    END PROCEDURE get_face_slave


    MODULE PROCEDURE get_face_flag
        IMPLICIT NONE

        get_face_flag = f%flag

    END PROCEDURE get_face_flag


    ! ----- Setter -----

    MODULE PROCEDURE set_face
        IMPLICIT NONE

        IF(PRESENT(nv))     f%nv     = nv
        IF(PRESENT(master)) f%master = master
        IF(PRESENT(slave))  f%slave  = slave
        IF(PRESENT(flag))   f%flag   = flag

    END PROCEDURE set_face


END SUBMODULE class_face_procedures
