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
! $Id: class_vector.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_vector) class_vector_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_vector_sizeof
        USE psb_base_mod

        nemo_vector_sizeof = 3*nemo_sizeof_dp

    END PROCEDURE nemo_vector_sizeof
    ! ----- Constructors -----

    ! Public default constructor
    MODULE PROCEDURE vector_

        vector_%x = x
        vector_%y = y
        vector_%z = z

    END PROCEDURE vector_


    ! Array constructor
    MODULE PROCEDURE alloc_vector
        !
        INTEGER :: i, info

        IF(ALLOCATED(vect))THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ALLOCATE(vect(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ! Initizializes to zero
        DO i = 1, n
            !vect(i) = vector(0.d0,0.d0,0.d0)
            !! Workaround for Intel 18 error #6053: Structure constructor may not have fields with a PRIVATE attribute
            vect(i)%x = 0.d0
            vect(i)%y = 0.d0
            vect(i)%z = 0.d0
        END DO

100     FORMAT(' ERROR! Illegal allocation: VECTOR array already associated')
200     FORMAT(' ERROR! Memory allocation failure in ALLOC_VECTOR')

    END PROCEDURE alloc_vector


    ! ----- Destructor -----

    MODULE PROCEDURE free_vector
        !
        INTEGER :: info

        IF (ALLOCATED(vect)) THEN
            DEALLOCATE(vect,stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_VECTOR')

    END PROCEDURE free_vector


    ! ----- Parallel Operations -----

    MODULE PROCEDURE bcast_vector
        !
        INTEGER :: icontxt, mypnum
        INTEGER :: i, info, n
        REAL(psb_dpk_), ALLOCATABLE :: realbuf(:,:)

        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! VECT(:) is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            n = SIZE(vect)
            CALL psb_bcast(icontxt,n)

            IF(.NOT.ALLOCATED(vect)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        ELSE
            CALL psb_bcast(icontxt,n)

            IF(ALLOCATED(vect) .AND. SIZE(vect) /= n) THEN
                CALL free_vector(vect)
            END IF

            IF(.NOT.ALLOCATED(vect)) THEN
                CALL alloc_vector(vect,n)
            END IF
        END IF

        ALLOCATE(realbuf(3,n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        IF(mypnum == 0) THEN
            DO i = 1, n
                realbuf(1,i) = vect(i)%x
                realbuf(2,i) = vect(i)%y
                realbuf(3,i) = vect(i)%z
            END DO

            CALL psb_bcast(icontxt,realbuf)
        ELSE
            CALL psb_bcast(icontxt,realbuf)

            DO i= 1, n
                vect(i)%x = realbuf(1,i)
                vect(i)%y = realbuf(2,i)
                vect(i)%z = realbuf(3,i)
            END DO
        END IF

        DEALLOCATE(realbuf)

100     FORMAT(' ERROR! Illegal Send: VECTOR array not associated')
200     FORMAT(' ERROR! Memory allocation failure BCAST_VECTOR')

    END PROCEDURE bcast_vector


    MODULE PROCEDURE g2l_vector
        USE psb_base_mod
        !
        INTEGER :: iv_glob, iv_loc
        INTEGER :: n_loc
        INTEGER, ALLOCATABLE  :: iloc_to_glob(:)
        TYPE(vector), ALLOCATABLE  :: work(:)


        n_loc  = psb_cd_get_local_cols(desc_v)

        CALL alloc_vector(work,n_loc)

        CALL psb_get_loc_to_glob(desc_v,iloc_to_glob)

        DO iv_loc = 1, n_loc
            iv_glob = iloc_to_glob(iv_loc)
            work(iv_loc) = verts(iv_glob)
        END DO

        ! Reallocates VERTS
        CALL free_vector(verts)
        CALL alloc_vector(verts,n_loc)
        verts(1:n_loc) = work(1:n_loc)

        DEALLOCATE(iloc_to_glob)
        CALL free_vector(work)

    END PROCEDURE g2l_vector


    MODULE PROCEDURE l2g_vector
        USE psb_base_mod
        ! WARNING! The global results is allocated only on P0. After its usage
        ! it must be deallocated in the calling unit by means of the statement:
        ! "if(associated(glob_res)) deallocate(glob_res)"

        !
        INTEGER :: err_act, info, icontxt
        INTEGER :: iv, n_glob, n_loc
        REAL(psb_dpk_), ALLOCATABLE  :: dbuf_glob(:,:)
        REAL(psb_dpk_), ALLOCATABLE  :: dbuf_loc(:,:)

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        n_loc = SIZE(verts_loc)             ! # local vertices
        n_glob = psb_cd_get_global_cols(desc_v) ! # global vertices

        ! Allocation of local objects
        CALL psb_geall(dbuf_loc,desc_v,info,3)
        CALL psb_check_error(info,'l2g_vector','psb_geall',icontxt)

        ALLOCATE(dbuf_glob(n_glob,3),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF


        DO iv = 1, n_loc
            dbuf_loc(iv,1) = verts_loc(iv)%x
            dbuf_loc(iv,2) = verts_loc(iv)%y
            dbuf_loc(iv,3) = verts_loc(iv)%z
        END DO


        CALL psb_ovrl(dbuf_loc,desc_v,info,update=psb_avg_)
        CALL psb_check_error(info,'l2g_vector','psb_ovrl',icontxt)

        CALL psb_gather(dbuf_glob,dbuf_loc,desc_v,info,root=0)
        CALL psb_check_error(info,'l2g_vector','psb_dgatherm',icontxt)


        IF(mypnum_() == 0) THEN
            ! Allocation of global object on P0
            IF(ALLOCATED(verts_glob)) CALL free_vector(verts_glob)
            CALL alloc_vector(verts_glob,n_glob)

            DO iv = 1, n_glob
                verts_glob(iv)%x = dbuf_glob(iv,1)
                verts_glob(iv)%y = dbuf_glob(iv,2)
                verts_glob(iv)%z = dbuf_glob(iv,3)
            END DO
        END IF


        ! Frees memory
        DEALLOCATE(dbuf_glob)

        CALL psb_gefree(dbuf_loc,desc_v,info)
        CALL psb_check_error(info,'l2g_vector','psb_gefree',icontxt)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in L2G_VECTOR')

    END PROCEDURE l2g_vector


    ! ----- Getters -----

    MODULE PROCEDURE get_vector_x

        get_vector_x = vert%x

    END PROCEDURE get_vector_x


    MODULE PROCEDURE get_vector_y

        get_vector_y = vert%y

    END PROCEDURE get_vector_y


    MODULE PROCEDURE get_vector_z

        get_vector_z = vert%z

    END PROCEDURE get_vector_z


    ! ----- Setters -----
    MODULE PROCEDURE set_vector_x

        vect%x = r

    END PROCEDURE set_vector_x

    MODULE PROCEDURE set_vector_y

        vect%y = r

    END PROCEDURE set_vector_y

    MODULE PROCEDURE set_vector_z

        vect%z = r

    END PROCEDURE set_vector_z



    MODULE PROCEDURE update_vector_halo
        USE psb_base_mod
        !
        CHARACTER(len=20), PARAMETER :: name_err = 'update_vector_halo'
        !
        INTEGER :: icontxt, info, err_act, ncells

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()
        ncells = psb_cd_get_local_cols(desc)

        IF(SIZE(v) /= ncells) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Updates halo elements
        CALL psb_halo(v%x,desc,info)
        CALL psb_check_error(info,name_err,'psb_halo',icontxt)

        CALL psb_halo(v%y,desc,info)
        CALL psb_check_error(info,name_err,'psb_halo',icontxt)

        CALL psb_halo(v%z,desc,info)
        CALL psb_check_error(info,name_err,'psb_halo',icontxt)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Illegal size of VECTOR array in UPDATE_VECTOR_HALO')

    END PROCEDURE update_vector_halo


    ! ----- Vector Algebra Operations -----

    MODULE PROCEDURE vec_sum

        vec_sum%x = a%x + b%x
        vec_sum%y = a%y + b%y
        vec_sum%z = a%z + b%z

    END PROCEDURE vec_sum


    MODULE PROCEDURE vec_diff

        vec_diff%x = a%x - b%x
        vec_diff%y = a%y - b%y
        vec_diff%z = a%z - b%z

    END PROCEDURE vec_diff

    MODULE PROCEDURE vec_minus

        vec_minus%x = -v%x
        vec_minus%y = -v%y
        vec_minus%z = -v%z
    END PROCEDURE vec_minus

    MODULE PROCEDURE scalar_vector_prod

        scalar_vector_prod%x = alpha * v%x
        scalar_vector_prod%y = alpha * v%y
        scalar_vector_prod%z = alpha * v%z

    END PROCEDURE scalar_vector_prod


    MODULE PROCEDURE dot_prod_v

        dot_prod_v = a%x * b%x + a%y * b%y + a%z * b%z

    END PROCEDURE dot_prod_v


    MODULE PROCEDURE dot_prod_t
        ! Used for GRAD .dot. DELTA in VECTOR_PDE_LAPLACIAN.
        ! GRAD is a tensor. It would need a proper class...

        dot_prod_t%x = a(1) .dot. b
        dot_prod_t%y = a(2) .dot. b
        dot_prod_t%z = a(3) .dot. b

    END PROCEDURE dot_prod_t


    MODULE PROCEDURE cross_prod

        cross_prod%x = a%y * b%z - a%z * b%y
        cross_prod%y = a%z * b%x - a%x * b%z
        cross_prod%z = a%x * b%y - a%y * b%x

    END PROCEDURE cross_prod


    MODULE PROCEDURE vec_mag
        !
        REAL(psb_dpk_) :: x, y

        x = MAX(ABS(v%x),ABS(v%y),ABS(v%z))
        IF (x == 0.0d0) THEN
            vec_mag = 0.0d0
        ELSE
            y = ((v%x/x)**2 + (v%y/x)**2 + (v%z/x)**2 )
            vec_mag = x*SQRT(y)
        ENDIF
    END PROCEDURE vec_mag


    MODULE PROCEDURE vec_unit
        ! Returns a unit vector in the direction of V
        !
        !
        REAL(psb_dpk_) :: vmag

        vmag = vec_mag(v)
        IF (vmag == 0.0d0) THEN

            vec_unit = v
        ELSE
            vec_unit = (1.0d0 / vmag) * v
        END IF

    END PROCEDURE vec_unit


    MODULE PROCEDURE vec_eq

        IF(    a%x == b%x .AND. &
            & a%y == b%y .AND. &
            & a%z == b%z) THEN
            vec_eq = .TRUE.
        ELSE
            vec_eq = .FALSE.
        END IF

    END PROCEDURE vec_eq


    MODULE PROCEDURE vector_to_array

        ! assume that the array is size 3... user beware

        a(1) = b%x
        a(2) = b%y
        a(3) = b%z

    END PROCEDURE vector_to_array

END SUBMODULE class_vector_procedures
